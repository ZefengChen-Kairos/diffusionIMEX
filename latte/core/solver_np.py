"""Numpy fp64 reference solver: 2D conservative active surface layer + 0D bulk ledger.

Coarse grid (N^2): conserved U = [M, Px, Py, Phi, Phi_1..Phi_K] (Phi_k = batch tags).
One face mass flux F drives every conserved field (Phi, tags and momentum are carried as
F * upwind(value)).  Fluxes: MUSCL (MC limiter) upwind advection + centred buoyant pressure
+ interface dissipation kappa grad(h) (kappa = 0.25 c' dx), SSP-RK2 in time.  Sources
(D, R, E, tau_j) are applied after the flux update with exact ledger bookkeeping; the drag
on the bulk is implicit.

Fine grid ((ratio N)^2): mass fractions Y = Phi/M and Y_k = Phi_k/M, semi-Lagrangian BFECC
with a local min/max clamp, sources applied in exact mixing form.  The fine field is the
material picture used for rendering and pattern metrics; the coarse Phi is the
conservative ledger (their consistency is a reported diagnostic).
"""
import math
import time
import numpy as np

from .grid import Grid, FineMap, Bilin, bilinear
from .params import PhysParams, NumParams
from closures.jet_stub import impact_event, footprint
from . import kernels_nb as kn


class ConservationError(RuntimeError):
    pass


def _mc(dl, dr):
    s = np.sign(dl)
    m = np.minimum(np.minimum(2 * np.abs(dl), 2 * np.abs(dr)), 0.5 * np.abs(dl + dr))
    return np.where(dl * dr > 0, s * m, 0.0)


def _minmod(dl, dr):
    return np.where(dl * dr > 0, np.sign(dl) * np.minimum(np.abs(dl), np.abs(dr)), 0.0)


class Bulk:
    """0D bulk ledger: milk-material / base-liquid / tag masses and horizontal momentum.

    bulk_return == "white" (default): Mm is the whipped (white) milk material in the bulk, i.e.
    what E entrained from the surface; liquid milk that drains or penetrates joins the base
    liquid Mc (its amount is counted in Ml).  R returns the white fraction Mm / M.
    bulk_return == "all_milk" (task-book form): every sunk milk mass goes to Mm and R returns it.
    """

    def __init__(self, P: PhysParams, K):
        self.Mm = 0.0
        self.Mc = P.rho_c * P.V_bulk0
        self.Ml = 0.0               # liquid milk merged into the base (part of Mc), white mode only
        self.Mk = np.zeros(K)
        self.Pm = np.zeros(2)       # bulk horizontal momentum
        self.P = P

    @property
    def M(self):
        return self.Mm + self.Mc

    @property
    def Y(self):
        return self.Mm / self.M

    @property
    def Yk(self):
        return self.Mk / self.M

    @property
    def rho(self):
        return self.M / (self.Mc / self.P.rho_c + self.Mm / self.P.rho_milk)


class Sim:
    def __init__(self, action=None, P: PhysParams = None, NP: NumParams = None, n_tags=None,
                 record_every=0.0, v_j_override=None, flags=None, backend="numba"):
        self.P = P or PhysParams()
        self.NP = NP or NumParams()
        P, NP = self.P, self.NP
        self.action = action
        self.backend = backend
        self.flags = dict(no_R=False, no_E=False, no_tau=False, no_D_momentum=False,
                          no_D_mass=False, kappa_off=False)
        if flags:
            self.flags.update(flags)
        self.v_j_override = v_j_override
        K = n_tags if n_tags is not None else (action.n_tags if action is not None else 1)
        self.K = max(K, 1)
        self.g = Grid(NP.L, NP.N, P.R_cup)
        self.gf = Grid(NP.L, NP.N * NP.ratio, P.R_cup)
        self.fm = FineMap(self.g, self.gf)
        N, Nf = self.g.N, self.gf.N
        wet = self.g.wet
        self.U = np.zeros((4 + self.K, N, N))
        self.U[0] = np.where(wet, P.rho_c * P.h0, 0.0)
        self.Yf = np.zeros((1 + self.K, Nf, Nf))
        self.bulk = Bulk(P, self.K)
        self.M_floor = 0.3 * P.rho_c * P.h0
        self.t = 0.0
        self.nstep = 0
        self.ubar_acc = np.zeros((2, N, N))
        self.dt_acc = 0.0
        self.k_acc = 0
        # ledgers (fp64)
        self.inj = dict(mass=0.0, milk=0.0, milk_total=0.0, tags=np.zeros(self.K), mom=np.zeros(2))
        self.cup_impulse = np.zeros(2)
        self.ret_mass = 0.0          # integral of D (retained new milk)
        self.E_milk = 0.0            # milk removed from the surface by E
        self.R_milk = 0.0            # milk returned by R
        self.E_tags = np.zeros(self.K)
        self.R_tags = np.zeros(self.K)
        self.X_sub = 0.0             # mass subducted by the thickness exchange
        self.X_up = 0.0              # bulk mass upwelled by the thickness exchange
        self.log = []                # per-step conservation errors
        self.frames = []             # (t, Yf copy, M copy)
        self.record_every = record_every
        self._next_rec = 0.0
        self.diag = []               # (t, r, v_j, D_j, ...)
        self.wall_time = 0.0
        self.E0 = self.totals()

    # ------------------------------------------------------------------ primitives
    def prims(self, U):
        P = self.P
        wet = self.g.wet
        M = np.where(wet, U[0], 1.0)
        u = np.where(wet, U[1] / M, 0.0)
        v = np.where(wet, U[2] / M, 0.0)
        Y = np.where(wet, U[3] / M, 0.0)
        Yk = np.where(wet, U[4:] / M, 0.0)
        h = (M - U[3]) / P.rho_c + U[3] / P.rho_m
        rb = M / h
        gp = np.maximum(P.g * (P.rho_c - rb) / P.rho_c, 0.0)
        if P.pressure == "buoyant":
            Ps = 0.5 * gp * rb * h * h
            c = np.sqrt(gp * h)
        else:
            Ps = 0.5 * P.K_s * (h * h - P.h0 * P.h0)
            c = np.sqrt(np.maximum(P.K_s * h * h / M, 0.0))
        Ps = np.where(wet, Ps, 0.0)
        c = np.where(wet, c, 0.0)
        return M, u, v, Y, Yk, h, rb, Ps, c

    # ------------------------------------------------------------------ fluxes
    def rhs(self, U):
        """dU/dt from the face fluxes; also returns the wall force (2,) on the layer."""
        if self.backend == "numba":
            return self._rhs_nb(U)
        g, P, NP = self.g, self.P, self.NP
        M, u, v, Y, Yk, h, rb, Ps, c = self.prims(U)
        W = np.concatenate([np.stack([M, u, v, Y]), Yk])        # primitive stack
        lim = _mc if NP.limiter == "mc" else _minmod
        dU = np.zeros_like(U)
        kap = 0.0 if self.flags["kappa_off"] else NP.kappa_fac * g.dx
        nu = P.nu_s
        dx = g.dx
        for ax in (2, 1):          # axis of W arrays: 2 = x (i), 1 = y (j)
            if ax == 2:
                op = g.ox
                Wl, Wr = W[:, :, :-1], W[:, :, 1:]
                d = Wr - Wl                                         # (nv, N, N-1)
                d = d * op
                z = np.zeros_like(W[:, :, :1])
                dL = np.concatenate([z, d], axis=2)
                dR = np.concatenate([d, z], axis=2)
                s = lim(dL, dR)
                WL = Wl + 0.5 * s[:, :, :-1]
                WR = Wr - 0.5 * s[:, :, 1:]
                hn = (h[:, :-1], h[:, 1:])
                rbf = 0.5 * (rb[:, :-1] + rb[:, 1:])
                cf = np.maximum(c[:, :-1], c[:, 1:])
                Pf = 0.5 * (Ps[:, :-1] + Ps[:, 1:])
                Mf = 0.5 * (M[:, :-1] + M[:, 1:])
                du = (u[:, 1:] - u[:, :-1], v[:, 1:] - v[:, :-1])
                ni, ti = 1, 2
            else:
                op = g.oy
                Wl, Wr = W[:, :-1, :], W[:, 1:, :]
                d = (Wr - Wl) * op
                z = np.zeros_like(W[:, :1, :])
                dL = np.concatenate([z, d], axis=1)
                dR = np.concatenate([d, z], axis=1)
                s = lim(dL, dR)
                WL = Wl + 0.5 * s[:, :-1, :]
                WR = Wr - 0.5 * s[:, 1:, :]
                hn = (h[:-1, :], h[1:, :])
                rbf = 0.5 * (rb[:-1, :] + rb[1:, :])
                cf = np.maximum(c[:-1, :], c[1:, :])
                Pf = 0.5 * (Ps[:-1, :] + Ps[1:, :])
                Mf = 0.5 * (M[:-1, :] + M[1:, :])
                du = (u[1:, :] - u[:-1, :], v[1:, :] - v[:-1, :])
                ni, ti = 2, 1
            un = 0.5 * (WL[ni] + WR[ni])
            F = un * np.where(un > 0, WL[0], WR[0])
            F = F - kap * cf * rbf * (hn[1] - hn[0]) / dx
            F = F * op
            up = F > 0
            Wup = np.where(up[None], WL, WR)
            flux = np.empty((U.shape[0],) + F.shape)
            flux[0] = F
            flux[1] = F * Wup[1]
            flux[2] = F * Wup[2]
            flux[ni] += Pf * op
            flux[1] -= nu * Mf * du[0] / dx * op
            flux[2] -= nu * Mf * du[1] / dx * op
            flux[3] = F * Wup[3]
            flux[4:] = F[None] * Wup[4:]
            if ax == 2:
                dU[:, :, :-1] -= flux / dx
                dU[:, :, 1:] += flux / dx
            else:
                dU[:, :-1, :] -= flux / dx
                dU[:, 1:, :] += flux / dx
        # wall pressure (reaction of the cup wall on closed faces)
        fx = -Ps * (g.wall_xr.astype(float) - g.wall_xl) / dx
        fy = -Ps * (g.wall_yt.astype(float) - g.wall_yb) / dx
        dU[1] += fx
        dU[2] += fy
        wall = np.array([fx.sum(), fy.sum()]) * g.dA
        dU[:, ~g.wet] = 0.0
        return dU, wall, (np.abs(u).max(), np.abs(v).max(), c.max())

    def _rhs_nb(self, U):
        g, P, NP = self.g, self.P, self.NP
        M, u, v, Y, Yk, h, rb, Ps, c = self.prims(U)
        W = np.ascontiguousarray(np.concatenate([np.stack([M, u, v, Y]), Yk]))
        kap = 0.0 if self.flags["kappa_off"] else NP.kappa_fac * g.dx
        dU = np.empty_like(U)
        wall = np.zeros(2)
        kn.rhs_kernel(U, W, h, rb, Ps, c, g.ox, g.oy, g.wall_xr, g.wall_xl, g.wall_yt, g.wall_yb,
                      g.wet, g.dx, kap, P.nu_s, NP.limiter == "mc", dU, wall)
        return dU, wall * g.dA, (np.abs(u).max(), np.abs(v).max(), c.max())

    # ------------------------------------------------------------------ totals
    def totals(self, U=None):
        U = self.U if U is None else U
        dA = self.g.dA
        return dict(M=U[0].sum() * dA, Px=U[1].sum() * dA, Py=U[2].sum() * dA,
                    Phi=U[3].sum() * dA, tags=U[4:].sum(axis=(1, 2)) * dA)

    # ------------------------------------------------------------------ sources
    def sources(self, dt, ev):
        """Apply D, R, E, tau_j to the coarse state; update the bulk ledger.

        Returns the fine-grid source description for the material step."""
        P, g, B, f = self.P, self.g, self.bulk, self.flags
        U = self.U
        K = self.K
        mdot = ev.mdot
        r = ev.r
        r0, sd, sk, sr = footprint(ev, P)
        x0 = ev.x
        # kernels
        jd, idd, Kd = g.kernel(x0, "ring", sd, r0)
        je, ie, Ke = g.kernel(x0, "gauss", sk)
        jr, ir, Kr = g.kernel(x0, "gauss", sr)
        Kt = Ke
        M = U[0]
        Yloc = U[3] / np.where(g.wet, M, 1.0)
        tag = ev.tag
        # --- retained deposit D
        Drate = 0.0 if f["no_D_mass"] else ev.phi_f * r * mdot
        mom_in = mdot * ev.v_h                              # incoming horizontal momentum rate
        dmom_surf = np.zeros(2)
        if Drate > 0:
            D = Drate * Kd
            U[0][jd, idd] += dt * D
            U[3][jd, idd] += dt * D
            if tag >= 0:
                U[4 + tag][jd, idd] += dt * D
            if not f["no_D_momentum"]:
                U[1][jd, idd] += dt * D * ev.u_d[0]
                U[2][jd, idd] += dt * D * ev.u_d[1]
                dmom_surf += dt * Drate * ev.u_d
        # --- penetrating part (1 - r) and drained liquid r (1 - phi_f) go to the bulk
        sunk = mdot * dt - Drate * dt
        white = P.bulk_return == "white"
        if white:
            B.Mc += sunk            # liquid milk joins the base liquid; it is not white material
            B.Ml += sunk
        else:
            B.Mm += sunk
            if tag >= 0:
                B.Mk[tag] += sunk
        # --- screened sink pair E / R (equal volumes)
        Qe = P.C_E * (1.0 - r) * ev.Q
        Erate_mass = Rrate_mass = 0.0
        if Qe > 0 and not f["no_E"]:
            Ms = M[je, ie]
            Mw = np.where(g.wet[je, ie], Ms, 1.0)
            hloc = (Mw - U[3][je, ie]) / P.rho_c + U[3][je, ie] / P.rho_m
            rho_s = np.where(g.wet[je, ie], Mw / hloc, 0.0)
            E = rho_s * Qe * Ke                            # kg/m^2/s
            # positivity: removal per step at most half the excess over the floor
            room = np.maximum(Ms - self.M_floor, 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                lim = np.where(E * dt > 0, 0.5 * room / (E * dt), np.inf)
            s = min(1.0, float(lim.min())) if lim.size else 1.0
            E *= s
            Qe_eff = Qe * s
            YE = Yloc[je, ie]
            Yk_E = U[4:, je, ie] / Mw
            uE = U[1][je, ie] / Mw
            vE = U[2][je, ie] / Mw
            Em = (E * dt)
            dA = g.dA
            e_mass = Em.sum() * dA
            e_milk = (Em * YE).sum() * dA
            e_tags = (Em[None] * Yk_E).sum(axis=(1, 2)) * dA
            e_mom = np.array([(Em * uE).sum(), (Em * vE).sum()]) * dA
            U[0][je, ie] -= Em
            U[3][je, ie] -= Em * YE
            U[4:, je, ie] -= Em[None] * Yk_E
            U[1][je, ie] -= Em * uE
            U[2][je, ie] -= Em * vE
            dmom_surf -= e_mom
            B.Mm += e_milk
            B.Mc += e_mass - e_milk
            B.Mk += e_tags
            self.E_milk += e_milk
            self.E_tags += e_tags
            Erate_mass = e_mass / dt
            if not f["no_R"]:
                Yb, Ybk, rhob = B.Y, B.Yk.copy(), B.rho
                R = rhob * Qe_eff * Kr
                Rm = R * dt
                r_mass = Rm.sum() * dA
                U[0][jr, ir] += Rm
                U[3][jr, ir] += Rm * Yb
                U[4:, jr, ir] += Rm[None] * Ybk[:, None, None]
                B.Mm -= r_mass * Yb
                B.Mc -= r_mass * (1.0 - Yb)
                B.Mk -= r_mass * Ybk
                self.R_milk += r_mass * Yb
                self.R_tags += r_mass * Ybk
                Rrate_mass = r_mass / dt
        # --- penetrating jet traction tau_j (implicit relaxation toward v_h)
        if not f["no_tau"]:
            beta = P.C_t * (1.0 - r) * mdot * Kt               # kg/m^2/s
            Mt = U[0][je, ie]
            Mtw = np.where(g.wet[je, ie], Mt, 1.0)
            for c_, vh in ((1, ev.v_h[0]), (2, ev.v_h[1])):
                Pold = U[c_][je, ie]
                Pnew = np.where(g.wet[je, ie], (Pold + dt * beta * vh) / (1.0 + dt * beta / Mtw), 0.0)
                imp = (Pnew - Pold).sum() * g.dA
                U[c_][je, ie] = Pnew
                dmom_surf[c_ - 1] += imp
        # the bulk receives the incoming horizontal momentum the surface did not keep
        # (mom_in - D u_d - tau_j) plus what E carried down (dmom_surf includes -E u)
        B.Pm += dt * mom_in - dmom_surf
        self.inj["mass"] += mdot * dt
        self.inj["milk_total"] += mdot * dt
        # material ledger input: white material only (white mode) or all milk (task-book form)
        dm = Drate * dt if white else mdot * dt
        self.inj["milk"] += dm
        if tag >= 0:
            self.inj["tags"][tag] += dm
        self.inj["mom"] += dt * mom_in
        self.ret_mass += Drate * dt
        return dict(x0=x0, r0=r0, sd=sd, sk=sk, sr=sr, Drate=Drate, tag=tag,
                    Qe=(Qe_eff if Erate_mass > 0 else 0.0), Rrate=Rrate_mass,
                    Yb=B.Y, Ybk=B.Yk.copy(), rhob=B.rho)

    def exchange(self, dt):
        """Vertical exchange relaxing the layer thickness to h0 (subduction / upwelling).

        Subduction (h > h0) removes the local mixture (local Y, tags, velocity) like E;
        upwelling (h < h0) brings bulk fluid (Y_b, u_b = 0) like R."""
        P, g, B = self.P, self.g, self.bulk
        U = self.U
        wet = g.wet
        M = np.where(wet, U[0], 1.0)
        h = (M - U[3]) / P.rho_c + U[3] / P.rho_m
        a = 1.0 - math.exp(-dt / P.tau_ex)
        dh = (h - P.h0) * a * wet                        # thickness to remove (>0) or add (<0)
        sub = np.maximum(dh, 0.0)
        frac = np.minimum(sub / h, 1.0 - self.M_floor / M)  # mass fraction subducted (local mixture)
        frac = np.maximum(frac, 0.0)
        dA = g.dA
        rem = U * frac[None]
        rem_tot = rem.sum(axis=(1, 2)) * dA
        U -= rem
        B.Mm += rem_tot[3]
        B.Mc += rem_tot[0] - rem_tot[3]
        B.Mk += rem_tot[4:]
        B.Pm += rem_tot[1:3]
        self.X_sub += rem_tot[0]
        # upwelling of bulk fluid
        up = np.maximum(-dh, 0.0)
        if up.any():
            rhob, Yb, Ybk = B.rho, B.Y, B.Yk.copy()
            add = rhob * up
            tot = add.sum() * dA
            U[0] += add
            U[3] += add * Yb
            U[4:] += add[None] * Ybk[:, None, None]
            B.Mm -= tot * Yb
            B.Mc -= tot * (1.0 - Yb)
            B.Mk -= tot * Ybk
            self.X_up += tot
            # fine grid: the upwelled bulk mixes into the material fractions
            self._up_acc = getattr(self, "_up_acc", 0.0) + add
            self._up_Yb = (Yb, Ybk)

    # ------------------------------------------------------------------ fine material
    def fine_sources(self, dt, src):
        """Exact mixing update of the fine mass fractions inside the source windows."""
        if src is None:
            return
        gf = self.gf
        Mc = self.g.extrapolate(self.U[0].copy())
        if src["Drate"] > 0:
            sj, si, Kd = gf.kernel(src["x0"], "ring", src["sd"], src["r0"])
            Mf = self.fm.to_fine_window(Mc, sj, si)
            D = src["Drate"] * Kd * dt
            Yw = self.Yf[:, sj, si]
            den = Mf + D
            Yw[0] = (Yw[0] * Mf + D) / den
            Yw[1:] = Yw[1:] * (Mf / den)[None]
            if src["tag"] >= 0:
                Yw[1 + src["tag"]] += D / den
            self.Yf[:, sj, si] = Yw
        if src["Rrate"] > 0:
            sj, si, Kr = gf.kernel(src["x0"], "gauss", src["sr"])
            Mf = self.fm.to_fine_window(Mc, sj, si)
            Rm = src["rhob"] * src["Qe"] * Kr * dt
            Yw = self.Yf[:, sj, si]
            den = Mf + Rm
            Yw[0] = (Yw[0] * Mf + Rm * src["Yb"]) / den
            Yw[1:] = (Yw[1:] * Mf[None] + Rm[None] * src["Ybk"][:, None, None]) / den[None]
            self.Yf[:, sj, si] = Yw

    def _fine_advect_nb(self, dt, ubar):
        gf, fm = self.gf, self.fm
        Nf = gf.N
        cvel = np.stack([fm.to_fine(ubar[0]), fm.to_fine(ubar[1])]) * (dt / gf.dx)
        cvel[:, gf.dry] = 0.0
        if not hasattr(self, "_dry_idx"):
            dj, di = np.nonzero(gf.dry)
            self._dry_idx = (dj, di, gf.near_j[dj, di], gf.near_i[dj, di])
        dj, di, nj, ni = self._dry_idx
        bj, bi = np.empty((Nf, Nf)), np.empty((Nf, Nf))
        fj, fi = np.empty((Nf, Nf)), np.empty((Nf, Nf))
        kn.trace_kernel(cvel, 1.0, gf.L, gf.dx, gf.R, bj, bi)
        kn.trace_kernel(cvel, -1.0, gf.L, gf.dx, gf.R, fj, fi)
        F = self.Yf.copy()
        kn.extrapolate_kernel(F, dj, di, nj, ni)
        f1 = np.empty_like(F)
        kn.sample_kernel(F, bj, bi, f1)
        kn.extrapolate_kernel(f1, dj, di, nj, ni)
        f0 = np.empty_like(F)
        kn.sample_kernel(f1, fj, fi, f0)
        corr = F + 0.5 * (F - f0)
        kn.extrapolate_kernel(corr, dj, di, nj, ni)
        out = np.empty_like(F)
        kn.bfecc_finish_kernel(F, corr, bj, bi, gf.wet, out)
        self.Yf = out

    def fine_advect(self, dt, ubar):
        if self.backend == "numba":
            return self._fine_advect_nb(dt, ubar)
        """Semi-Lagrangian BFECC with local min/max clamp on the fine grid.

        One RK2 back-trace (and one forward trace for the error estimate) per material
        step; the interpolation weights are shared by all material fields."""
        gf, fm = self.gf, self.fm
        Nf = gf.N
        dxf = gf.dx
        cu = fm.to_fine(ubar[0]) * (dt / dxf)
        cv = fm.to_fine(ubar[1]) * (dt / dxf)
        cu[gf.dry] = 0.0
        cv[gf.dry] = 0.0
        if not hasattr(self, "_jj"):
            self._jj, self._ii = np.meshgrid(np.arange(Nf, dtype=float), np.arange(Nf, dtype=float),
                                             indexing="ij")
        jj, ii = self._jj, self._ii
        rmax = gf.R - 0.5 * dxf
        cvel = np.stack([cu, cv])

        def trace(sign):
            m = Bilin(jj - 0.5 * sign * cv, ii - 0.5 * sign * cu, Nf)
            um, vm = m(cvel)
            jd = jj - sign * vm
            idp = ii - sign * um
            xc = (idp + 0.5) * dxf - 0.5 * gf.L
            yc = (jd + 0.5) * dxf - 0.5 * gf.L
            rr = np.sqrt(xc * xc + yc * yc)
            sc = np.where(rr > rmax, rmax / np.maximum(rr, 1e-30), 1.0)
            return Bilin((yc * sc + 0.5 * gf.L) / dxf - 0.5, (xc * sc + 0.5 * gf.L) / dxf - 0.5, Nf)

        bw = trace(+1.0)
        fw = trace(-1.0)
        F = gf.extrapolate(self.Yf.copy())
        f1 = gf.extrapolate(bw(F))
        f0 = fw(f1)
        corr = gf.extrapolate(F + 0.5 * (F - f0))
        out = bw(corr)
        lo, hi = bw.minmax(F)
        out = np.clip(out, lo, hi)
        out[1:] = np.minimum(out[1:], out[0][None])
        out[:, gf.dry] = 0.0
        self.Yf = out

    # ------------------------------------------------------------------ step
    def dt_cfl(self):
        M, u, v, Y, Yk, h, rb, Ps, c = self.prims(self.U)
        a = np.sqrt(u * u + v * v).max() + c.max()
        return min(self.NP.dt_max, self.NP.cfl * self.g.dx / max(a, 1e-9))

    def step(self, dt=None):
        tw = time.perf_counter()
        P, NP, g = self.P, self.NP, self.g
        if dt is None:
            dt = self.dt_cfl()
        smp = self.action.sample(self.t + 0.5 * dt) if self.action is not None else None
        vjo = self.v_j_override(self.t + 0.5 * dt) if callable(self.v_j_override) else self.v_j_override
        ev = impact_event(smp, P, vjo) if smp is not None else None
        T0 = self.totals()
        B = self.bulk
        bulk0 = (B.Mm, B.Mc, B.Mk.copy(), B.Pm.copy(), self.inj["milk"], self.inj["tags"].copy())
        cup0 = self.cup_impulse.copy()
        quiet = ev is None and not np.any(self.U[1:3]) and not np.any(self.U[3])
        if not quiet:
            # SSP-RK2 for the flux part
            U0 = self.U
            k1, w1, sp = self.rhs(U0)
            U1 = U0 + dt * k1
            k2, w2, _ = self.rhs(U1)
            self.U = 0.5 * (U0 + U1 + dt * k2)
            self.cup_impulse -= 0.5 * dt * (w1 + w2)       # reaction on the cup
            src = None
            if ev is not None:
                src = self.sources(dt, ev)
                self.diag.append((self.t, ev.r, ev.v_j, ev.D_j, ev.x[0], ev.x[1], ev.mdot))
            if P.tau_ex > 0:
                self.exchange(dt)
            # implicit drag (and optional yield) against the bulk
            Pv = self.U[1:3].copy()
            Pn = Pv / (1.0 + P.gamma * dt)
            if P.tau_y > 0:
                M = np.where(g.wet, self.U[0], 1.0)
                Yl = self.U[3] / M
                mag = np.sqrt(Pn[0] ** 2 + Pn[1] ** 2)
                red = P.tau_y * Yl * dt
                fac = np.where(mag > red, (mag - red) / np.maximum(mag, 1e-300), 0.0)
                Pn = Pn * fac[None]
            self.U[1:3] = Pn
            self.bulk.Pm += (Pv - Pn).sum(axis=(1, 2)) * g.dA
            self.fine_sources(dt, src)
            # accumulate the time-averaged velocity for the material step
            Mw = np.where(g.wet, self.U[0], 1.0)
            self.ubar_acc[0] += dt * np.where(g.wet, self.U[1] / Mw, 0.0)
            self.ubar_acc[1] += dt * np.where(g.wet, self.U[2] / Mw, 0.0)
        self.dt_acc += dt
        self.k_acc += 1
        self.t += dt
        self.nstep += 1
        if self.k_acc >= NP.mat_every:
            self.flush_material()
        if NP.check_every and self.nstep % NP.check_every == 0:
            self.check(T0, bulk0, cup0, dt, ev)
        if self.record_every and self.t >= self._next_rec - 1e-12:
            self.record()
            self._next_rec += self.record_every
        self.wall_time += time.perf_counter() - tw
        return dt

    def flush_material(self):
        up = getattr(self, "_up_acc", None)
        if up is not None and np.any(up):
            # mixing of upwelled bulk fluid into the fine fractions (exact mixing form)
            Yb, Ybk = self._up_Yb
            Mc = self.g.extrapolate(self.U[0].copy())
            Mf = self.fm.to_fine(Mc)
            A = self.fm.to_fine(up) * self.gf.wet
            Mold = np.maximum(Mf - A, 1e-12)
            self.Yf[0] = (self.Yf[0] * Mold + A * Yb) / (Mold + A)
            self.Yf[1:] = (self.Yf[1:] * Mold[None] + A[None] * Ybk[:, None, None]) / (Mold + A)[None]
            self._up_acc = 0.0
        if self.dt_acc > 0 and np.any(self.ubar_acc):
            self.fine_advect(self.dt_acc, self.ubar_acc / self.dt_acc)
        self.ubar_acc[:] = 0.0
        self.dt_acc = 0.0
        self.k_acc = 0

    # ------------------------------------------------------------------ checks
    def check(self, T0, bulk0, cup0, dt, ev):
        T1 = self.totals()
        B = self.bulk
        tol = self.NP.tol
        mdt = ev.mdot * dt if ev is not None else 0.0
        # (1) total mass: surface + bulk changes by mdot dt
        dtot = (T1["M"] - T0["M"]) + (B.Mm + B.Mc - bulk0[0] - bulk0[1]) - mdt
        e_mass = abs(dtot) / max(T0["M"] + bulk0[0] + bulk0[1], 1e-30)
        # (2) milk material: surface Phi + bulk material changes by the material input
        dm_in = self.inj["milk"] - bulk0[4]
        dmilk = (T1["Phi"] - T0["Phi"]) + (B.Mm - bulk0[0]) - dm_in
        e_milk = abs(dmilk) / max(T0["Phi"] + bulk0[0] + dm_in, 1e-12)
        # (3) tags
        dtag = (T1["tags"] - T0["tags"]) + (B.Mk - bulk0[2]) - (self.inj["tags"] - bulk0[5])
        e_tag = np.abs(dtag).max() / max(T0["Phi"] + bulk0[0] + dm_in, 1e-12)
        # (4) momentum: surface + bulk + cup changes by incoming momentum
        mom_in = ev.mdot * ev.v_h * dt if ev is not None else np.zeros(2)
        dP = (np.array([T1["Px"] - T0["Px"], T1["Py"] - T0["Py"]]) + (B.Pm - bulk0[3])
              + (self.cup_impulse - cup0) - mom_in)
        scale = max(T0["M"] * 1e-2, np.abs(mom_in).sum(), 1e-30)
        e_mom = np.abs(dP).max() / scale
        M = self.U[0][self.g.wet]
        pos = (M.min() >= self.M_floor * (1 - 1e-9),
               (self.U[3][self.g.wet] >= -1e-12 * M).all(),
               (self.U[3][self.g.wet] <= M * (1 + 1e-12)).all())
        self.log.append((self.t, e_mass, e_milk, e_tag, e_mom, float(M.min())))
        if max(e_mass, e_milk, e_tag, e_mom) > tol:
            raise ConservationError(f"t={self.t:.4f} mass {e_mass:.2e} milk {e_milk:.2e} "
                                    f"tag {e_tag:.2e} mom {e_mom:.2e}")
        if not pos[0]:
            self.pos_violation = getattr(self, "pos_violation", 0) + 1

    # ------------------------------------------------------------------ io
    def record(self):
        self.frames.append((self.t, self.Yf.copy(), self.U[0].copy()))

    def phi_fine(self, Yf=None, M=None):
        Yf = self.Yf if Yf is None else Yf
        M = self.U[0] if M is None else M
        Mf = self.fm.to_fine(np.where(self.g.wet, M, M[self.g.wet].mean()))
        return Yf[0] * Mf * self.gf.wet, Mf

    def run(self, T=None, progress=False):
        T = self.action.T_end if T is None else T
        if self.record_every and not self.frames:
            self.record()
            self._next_rec = self.record_every
        while self.t < T - 1e-12:
            dt = min(self.dt_cfl(), T - self.t)
            self.step(dt)
            if progress and self.nstep % 1000 == 0:
                print(f"  t={self.t:.3f} steps={self.nstep} wall={self.wall_time:.1f}s", flush=True)
        self.flush_material()
        return self
