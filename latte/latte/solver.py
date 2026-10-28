"""Time step of the layered non-hydrostatic cup solver (V3 Sec. 4, 6, 7).

One reference step (first-order splitting, V3 Sec. 7.4):
  1. sample the inlet at the step midpoint; build the surface arrival (same quadrature for
     volume, components and momentum)
  2. move the free surface with the stored divergence-free face fluxes U^n + arrival:
     D^{n+1} = D^n + dt (W_s + j) / A    (kinematic condition, column form)
  3. transport volume, components, mass and momentum with the *relative* face volumes
     (lateral: dt U; interfaces: dt W - A S dD; surface: the arrival only)
  4. implicit viscosity on the new geometry (backward Euler)
  5. one common pressure solve on the new geometry with p = p_air + sigma kappa on the
     free surface: face fluxes are projected exactly, cell velocities with the face-averaged
     net acceleration (gravity + pressure), as in Basilisk centered.h
  6. checks (bounds, liquid fraction, divergence, finiteness); any failure rolls the whole
     step back; numerical failures are retried with dt/2, input/geometry failures are not.
"""
from dataclasses import dataclass, field
import numpy as np

from .mesh import Geometry
from .operators import MetricOperator
from .transport import transport
from .linsolve import SPDSolver
from .inlet import surface_source, surface_source_subgrid, NoInlet, NO_FLOW
from .failures import SimFailure


@dataclass
class Model:
    rho: tuple = (1.0, 1.0, 1.0e-3)       # densities of B (base), M (milk liquid), G (effective gas)
    mu: float = 1.0e-3                     # constant dynamic viscosity (stage 1)
    g: float = 1.0
    sigma: float = 0.0                     # liquid-air surface tension
    liquid_min: float = 0.05               # required lower bound of C_B + C_M in every cell
    wall_bc: str = "no_slip"               # lateral walls: 'no_slip' | 'free_slip'
    bottom_bc: str = "no_slip"


@dataclass
class Numerics:
    cfl: float = 0.25
    wave_safety: float = 0.5
    cap_safety: float = 0.5
    scan_safety: float = 0.25
    max_dt: float = np.inf
    retry_limit: int = 12
    transport_scheme: str = "fct"          # 'fct' (selective joint FCT) | 'donor'
    ho_limiter: str = "superbee"           # high-order face composition: 'superbee' | 'vanleer' | 'lw'
    momentum_recon: str = "minmod"
    pressure_solver: str = "auto"
    pressure_rtol: float = 1e-12
    viscosity_solver: str = "auto"
    viscosity_rtol: float = 1e-12
    quad_nr: int = 8
    quad_ntheta: int = 32
    quad_max_spacing_cells: float = 0.25   # polar mode: inlet samples at most this many cells apart
    inlet_mapping: str = "subgrid"         # 'subgrid' (mesh-aligned sampling, exact totals) | 'polar'
    inlet_subgrid: int = 8
    connect_tol: float = None
    inlet_bc: str = "velocity"             # 'velocity': prescribed normal velocity on the patch, pressure from the
                                           # projection (V3 Sec. 5.2); 'surface_flux': p = p_air + flux source (comparison)
    inlet_patch_cover: float = 0.5         # columns covered at least this much form the velocity patch
    div_tol: float = 1e-8
    bound_tol: float = 1e-12
    max_surface_slope: float = 1.0


@dataclass
class State:
    t: float
    D: np.ndarray
    Vc: np.ndarray            # (2, nc) volumes of B and M
    Vtot: np.ndarray          # (nc,) material volume
    M: np.ndarray             # (nc,) mass
    P: np.ndarray             # (3, nc) momentum
    U: np.ndarray             # (nf,) projected absolute face volume fluxes [x, y, interfaces, surface]
    pi: np.ndarray            # (nc,) pressure (gauge p_air = 0)
    acc: dict = field(default_factory=dict)
    step: int = 0
    rejected: int = 0

    def copy(self):
        return State(self.t, self.D.copy(), self.Vc.copy(), self.Vtot.copy(), self.M.copy(),
                     self.P.copy(), self.U.copy(), self.pi.copy(),
                     {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.acc.items()},
                     self.step, self.rejected)

    def composition(self):
        VG = self.Vtot - self.Vc[0] - self.Vc[1]
        return np.stack([self.Vc[0], self.Vc[1], VG]) / self.Vtot[None, :]

    def velocity(self):
        return self.P / self.M[None, :]


def _new_acc():
    return {"in_vol": 0.0, "in_comp": np.zeros(3), "in_mass": 0.0, "in_mom": np.zeros(3),
            "out_vol": 0.0, "out_comp": np.zeros(3), "out_mass": 0.0, "out_mom": np.zeros(3),
            "imp_pressure": np.zeros(3), "imp_viscous": np.zeros(3), "adv_mom_residual": 0.0}


class Solver:
    def __init__(self, mesh, model=None, numerics=None, inlet=None, bc_flux=None):
        self.mesh = mesh
        self.model = model or Model()
        self.num = numerics or Numerics()
        self.inlet = inlet or NoInlet()
        self.bc_flux = bc_flux                 # optional prescribed boundary outflow per cell (tests)
        self.psolver = SPDSolver(self.num.pressure_solver, self.num.pressure_rtol, fail_code="PRESSURE_FAILED")
        self.vsolver = SPDSolver(self.num.viscosity_solver, self.num.viscosity_rtol, fail_code="VISCOSITY_FAILED")
        self.rho = np.asarray(self.model.rho, float)
        m = mesh
        self.n_tf = m.n_tf
        self.nf = m.n_tf + m.ncols
        self.log = []

    # ------------------------------------------------------------------ setup
    def initial_state(self, D0, C, v0=None, t0=0.0):
        m = self.mesh
        D0 = np.broadcast_to(np.asarray(D0, float), (m.ncols,)).copy()
        g0 = Geometry(m, D0)
        C = np.asarray(C, float).reshape(3, m.nc)
        if np.any(C < 0) or np.max(np.abs(C.sum(0) - 1)) > 1e-12:
            raise SimFailure("INPUT_INVALID", "initial composition must lie on the simplex")
        Vtot = g0.V.ravel().copy()
        Vc = C[:2] * Vtot[None, :]
        M = (self.rho @ C) * Vtot
        v0 = np.zeros((3, m.nc)) if v0 is None else np.asarray(v0, float).reshape(3, m.nc)
        st = State(t0, D0, Vc, Vtot, M, M[None, :] * v0, np.zeros(self.nf), np.zeros(m.nc), _new_acc())
        if np.any(v0 != 0):
            # divergence-free starting fluxes: project the interpolated initial velocity (no gravity)
            U, pi, _ = self._project(g0, M / Vtot, v0, dt=1.0, gravity=False)[:3]
            st.U = U
        return st

    # ------------------------------------------------------------------ pieces
    def _face_rho(self, g, rho_c):
        m = self.mesh
        r = rho_c.reshape(m.N, m.ncols)
        rx = 0.5 * (r[:, m.xm] + r[:, m.xp])
        ry = 0.5 * (r[:, m.ym] + r[:, m.yp])
        rz = np.empty_like(r)
        rz[:-1] = (g.h[:-1] * r[:-1] + g.h[1:] * r[1:]) / (g.h[:-1] + g.h[1:])
        rz[-1] = r[-1]
        return rx, ry, rz

    def _interp(self, g, v):
        """Face volume fluxes interpolated from cell-centred velocity (no gravity)."""
        m = self.mesh
        N, nc = m.N, m.ncols
        u, w_, vv = (v[d].reshape(N, nc) for d in (0, 2, 1))
        Ux = g.afx * 0.5 * (u[:, m.xm] + u[:, m.xp])
        Uy = g.afy * 0.5 * (vv[:, m.ym] + vv[:, m.yp])
        Uz = np.empty((N, nc))
        Uz[:-1] = m.A * (0.5 * (w_[:-1] + w_[1:]) - g.zx_int[:-1] * 0.5 * (u[:-1] + u[1:])
                         - g.zy_int[:-1] * 0.5 * (vv[:-1] + vv[1:]))
        Uz[-1] = m.A * (w_[-1] - g.zx_int[-1] * u[-1] - g.zy_int[-1] * vv[-1])
        return np.concatenate([Ux.ravel(), Uy.ravel(), Uz[:-1].ravel(), Uz[-1]])

    def _kappa(self, g):
        m = self.mesh
        ex, ey = g.etax_face, g.etay_face
        eyx = 0.5 * (g.etay_col[m.xm] + g.etay_col[m.xp])
        exy = 0.5 * (g.etax_col[m.ym] + g.etax_col[m.yp])
        fx = ex / np.sqrt(1 + ex ** 2 + eyx ** 2) * m.dy
        fy = ey / np.sqrt(1 + ey ** 2 + exy ** 2) * m.dx
        div = np.zeros(m.ncols)
        np.add.at(div, m.xm, fx); np.add.at(div, m.xp, -fx)
        np.add.at(div, m.ym, fy); np.add.at(div, m.yp, -fy)
        # wall columns: zero flux of the surface gradient (90 degree contact angle placeholder)
        return -div / m.A

    def _project(self, g, rho_c, v, dt, gravity=True, x0=None, inlet_cols=None, W_inlet=None):
        """Common pressure projection.  Free-surface faces are Dirichlet (p = sigma kappa) except on
        inlet-patch columns, whose surface face carries the prescribed absolute flux of the arriving
        material (V3 Sec. 5.2: prescribed normal velocity, pressure from the projection)."""
        m = self.mesh
        N, nc, ntf = m.N, m.ncols, m.n_tf
        smask = np.ones(nc, dtype=bool)
        inlet_cols = np.zeros(0, dtype=np.int64) if inlet_cols is None else np.asarray(inlet_cols, dtype=np.int64)
        smask[inlet_cols] = False
        rx, ry, rz = self._face_rho(g, rho_c)
        op = MetricOperator(m, g, 1.0 / rx, 1.0 / ry, 1.0 / rz, surface="dirichlet", surface_mask=smask)
        Ut = self._interp(g, v)
        Us = Ut.copy()
        if gravity and self.model.g != 0.0:
            Us[N * m.npx + N * m.npy:] -= self.model.g * m.A * dt
        Us[ntf + inlet_cols] = W_inlet if inlet_cols.size else Us[ntf + inlet_cols]
        sc = op.surface_cols
        Uop = np.concatenate([Us[:ntf], Us[ntf + sc]])
        fixed = np.zeros(m.nc)                                   # prescribed boundary fluxes (outward +)
        if inlet_cols.size:
            np.add.at(fixed, m.cell(N - 1, inlet_cols), W_inlet)
        if self.bc_flux is not None:
            fixed = fixed + self.bc_flux
        ps = self.model.sigma * self._kappa(g) if self.model.sigma > 0 else np.zeros(nc)
        div = -(op.G.T @ Uop) + fixed
        rhs = -div / dt - op.ghost_rhs(ps)
        pi, info = self.psolver.solve(op.L, rhs, x0=x0)
        F = op.flux(pi, ps)
        Uop_n = Uop - dt * F
        Un = Us.copy()
        Un[:ntf] = Uop_n[:ntf]
        Un[ntf + sc] = Uop_n[ntf:]
        divn = -(op.G.T @ Uop_n) + fixed
        info["e_div"] = float(dt * np.max(np.abs(divn) / g.V.ravel()))
        info["div_abs"] = float(np.max(np.abs(divn)))
        return Un, pi, info, Ut, op

    def _cell_acceleration(self, g, Un, Ut, dt):
        """Cell-centred net acceleration (gravity + pressure) from face accelerations."""
        m = self.mesh
        N, nc = m.N, m.ncols
        a = (Un - Ut) / dt
        nx_, ny_ = N * m.npx, N * m.npy
        ax_f = a[:nx_].reshape(N, m.npx) / g.afx
        ay_f = a[nx_:nx_ + ny_].reshape(N, m.npy) / g.afy
        az = a[nx_ + ny_:]
        az_f = np.concatenate([az[:(N - 1) * nc].reshape(N - 1, nc), az[(N - 1) * nc:][None, :]]) / m.A
        ax = np.zeros((N, nc)); ay = np.zeros((N, nc))
        np.add.at(ax.T, m.xm, 0.5 * ax_f.T); np.add.at(ax.T, m.xp, 0.5 * ax_f.T)
        np.add.at(ay.T, m.ym, 0.5 * ay_f.T); np.add.at(ay.T, m.yp, 0.5 * ay_f.T)
        ax_i = np.empty((N, nc)); ay_i = np.empty((N, nc))
        ax_i[:-1] = 0.5 * (ax[:-1] + ax[1:]); ax_i[-1] = ax[-1]
        ay_i[:-1] = 0.5 * (ay[:-1] + ay[1:]); ay_i[-1] = ay[-1]
        aw_f = az_f + g.zx_int * ax_i + g.zy_int * ay_i
        aw = 0.5 * aw_f
        aw[1:] += 0.5 * aw_f[:-1]                     # bottom wall face contributes zero
        return np.stack([ax.ravel(), ay.ravel(), aw.ravel()])

    def _viscosity(self, g, M, v, dt):
        m = self.mesh
        mu = self.model.mu
        N = m.N
        op = MetricOperator(m, g, np.full((N, m.npx), mu), np.full((N, m.npy), mu), np.full((N, m.ncols), mu),
                            surface="neumann", walls=(self.model.wall_bc == "no_slip"),
                            bottom=(self.model.bottom_bc == "no_slip"), c_wall=mu)
        from scipy import sparse
        A = (sparse.diags(M / dt) + op.L).tocsr()
        xs, info = self.vsolver.solve_many(A, [M * v[d] / dt for d in range(3)], [v[d] for d in range(3)])
        return np.stack(xs), {"iters": info["iters"], "rel_res": info["rel_res"]}

    def _source(self, inl, D, connect_tol):
        num, m = self.num, self.mesh
        if num.inlet_mapping == "subgrid":
            return surface_source_subgrid(inl, m, D, self.rho, num.inlet_subgrid, connect_tol, num.quad_nr, num.quad_ntheta)
        return surface_source(inl, m, D, self.rho, num.quad_nr, num.quad_ntheta, connect_tol,
                              num.quad_max_spacing_cells * min(m.dx, m.dy))

    # ------------------------------------------------------------------ time step control
    def stable_dt(self, st):
        m, num, mod = self.mesh, self.num, self.model
        hmin = min(m.dx, m.dy)
        nx_, ny_ = m.N * m.npx, m.N * m.npy
        rate = np.abs(st.U[:m.n_tf])
        Ws = st.U[m.n_tf:]
        dDdt = np.abs(Ws) / m.A
        try:
            inl = self.inlet.sample(st.t)
        except SimFailure:
            inl = None
        if inl is not None and inl.status != NO_FLOW:
            src = self._source(inl, st.D, None)
            dDdt = dDdt + src.vol_rate / m.A
        rz = rate[nx_ + ny_:].reshape(m.N - 1, m.ncols) + m.A * m.S[1:-1, None] * dDdt[None, :]
        tot = (np.bincount(m.tf_cm[:nx_ + ny_], weights=rate[:nx_ + ny_], minlength=m.nc)
               + np.bincount(m.tf_cp[:nx_ + ny_], weights=rate[:nx_ + ny_], minlength=m.nc))
        tot[:-m.ncols] += rz.ravel()
        tot[m.ncols:] += rz.ravel()
        tot[-m.ncols:] += np.abs(Ws) + (src.vol_rate if (inl is not None and inl.status != NO_FLOW) else 0.0)
        cands = {"max_dt": num.max_dt}
        if np.any(tot > 0):
            cands["cfl"] = num.cfl * float(np.min(st.Vtot[tot > 0] / tot[tot > 0]))
        if mod.g > 0:
            cands["wave"] = num.wave_safety * np.sqrt(hmin / (np.pi * mod.g))
        if mod.sigma > 0:
            rho_min = float(np.min(st.M / st.Vtot))
            cands["capillary"] = num.cap_safety * np.sqrt(rho_min * hmin ** 3 / (np.pi * mod.sigma))
        if inl is not None and inl.status != NO_FLOW and inl.scan_speed > 0:
            cands["scan"] = num.scan_safety * min(hmin, inl.section.a, inl.section.b) / inl.scan_speed
        ev = self.inlet.next_event(st.t)
        if np.isfinite(ev):
            cands["event"] = ev - st.t
        key = min(cands, key=cands.get)
        return float(cands[key]), key

    # ------------------------------------------------------------------ the step
    def _step(self, st, dt):
        m, mod, num = self.mesh, self.model, self.num
        N, A = m.N, m.A
        inl = self.inlet.sample(st.t + 0.5 * dt)
        src = self._source(inl, st.D, num.connect_tol)
        Ws = st.U[m.n_tf:].copy()
        if num.inlet_bc == "velocity":
            inlet_cols = np.nonzero((src.vol_rate > 0) & (src.cover >= num.inlet_patch_cover))[0]
        elif num.inlet_bc == "surface_flux":
            inlet_cols = np.zeros(0, dtype=np.int64)
        else:
            raise ValueError(num.inlet_bc)
        W_in = src.abs_rate[inlet_cols]
        # the free surface moves with the SAME stored, projected fluxes that transport uses (exact GCL);
        # a patch column therefore follows the section motion with a one-step lag
        dVJ = src.vol_rate * dt
        dD = (dt * Ws + dVJ) / A
        D1 = st.D + dD
        if np.any(D1 <= 0):
            raise SimFailure("BOUNDS_FAILED", "non-positive liquid depth")
        dV = dt * st.U[:m.n_tf].copy()
        nz0 = N * m.npx + N * m.npy
        dV[nz0:] -= (A * m.S[1:N, None] * dD[None, :]).ravel()
        hit = np.nonzero(dVJ > 0)[0]
        inflow = {"cell": m.cell(N - 1, hit), "vol": dVJ[hit], "comp": src.comp_rate[:, hit] * dt,
                  "mass": src.mass_rate[hit] * dt, "mom": src.mom_rate[:, hit] * dt}
        outflow = None
        if self.bc_flux is not None:
            o = np.nonzero(self.bc_flux > 0)[0]
            outflow = {"cell": o, "vol": self.bc_flux[o] * dt}
        tr = transport(m, st.Vc, st.Vtot, st.P, st.M, dV, self.rho, inflow, outflow,
                       num.transport_scheme, num.momentum_recon, num.bound_tol, num.ho_limiter)
        g1 = Geometry(m, D1)
        if g1.max_slope > num.max_surface_slope:
            # retried with smaller dt first; if it persists it is reported as OUT_OF_SCOPE
            raise SimFailure("BOUNDS_FAILED", f"free-surface slope {g1.max_slope:.3f} exceeds the layered representation",
                             geometry=True, slope=g1.max_slope)
        adv_res = float(np.max(np.abs(tr.P.sum(1) - (st.P.sum(1) + tr.in_mom.sum(1) - tr.out_mom.sum(1)))))
        v = tr.P / tr.M[None, :]
        rho_c = tr.M / tr.Vtot
        vinfo = {"iters": 0, "rel_res": 0.0}
        v_adv = v
        if mod.mu > 0:
            v, vinfo = self._viscosity(g1, tr.M, v, dt)
        imp_visc = (tr.M[None, :] * (v - v_adv)).sum(1)
        Un, pi, pinfo, Ut, op = self._project(g1, rho_c, v, dt, gravity=True, x0=st.pi,
                                              inlet_cols=inlet_cols, W_inlet=W_in)
        if pinfo["e_div"] > num.div_tol:
            raise SimFailure("PRESSURE_FAILED", f"divergence error {pinfo['e_div']:.2e} > {num.div_tol}")
        acc = self._cell_acceleration(g1, Un, Ut, dt)
        vn = v + dt * acc
        Pn = tr.M[None, :] * vn
        if not (np.all(np.isfinite(Pn)) and np.all(np.isfinite(Un))):
            raise SimFailure("BOUNDS_FAILED", "non-finite state")
        C = np.stack([tr.Vc[0], tr.Vc[1], tr.Vtot - tr.Vc[0] - tr.Vc[1]]) / tr.Vtot[None, :]
        liq = C[0] + C[1]
        if liq.min() < mod.liquid_min:
            raise SimFailure("OUT_OF_SCOPE", f"liquid fraction {liq.min():.3f} below the model bound {mod.liquid_min}")

        new = State(st.t + dt, D1, tr.Vc, tr.Vtot, tr.M, Pn, Un, pi,
                    {k: (v_.copy() if isinstance(v_, np.ndarray) else v_) for k, v_ in st.acc.items()},
                    st.step + 1, st.rejected)
        a = new.acc
        a["in_vol"] += float(tr.in_vol.sum()); a["in_comp"] += tr.in_comp.sum(1)
        a["in_mass"] += float(tr.in_mass.sum()); a["in_mom"] += tr.in_mom.sum(1)
        a["out_vol"] += float(tr.out_vol.sum()); a["out_comp"] += tr.out_comp.sum(1)
        a["out_mass"] += float(tr.out_mass.sum()); a["out_mom"] += tr.out_mom.sum(1)
        a["imp_pressure"] += (tr.M[None, :] * dt * acc).sum(1)
        a["imp_viscous"] += imp_visc
        a["adv_mom_residual"] = max(a["adv_mom_residual"], adv_res)
        geo_vol = g1.V.ravel()
        diag = {"t": new.t, "dt": dt, "e_div": pinfo["e_div"], "p_iters": pinfo["iters"], "p_res": pinfo["rel_res"],
                "v_iters": vinfo["iters"], "v_res": vinfo["rel_res"], "cfl": tr.cfl,
                "vol_mismatch": float(np.max(np.abs(tr.Vtot - geo_vol) / geo_vol)),
                "mass_consistency": tr.mass_consistency, "fct_limited": tr.fct_limited_fraction, "Q": src.Q,
                "inlet_patch_cols": int(inlet_cols.size), "max_slope": g1.max_slope,
                "cmin": float(C.min()), "cmax": float(C.max()), "liq_min": float(liq.min())}
        return new, diag

    def advance(self, st, dt_try=None):
        dt, limiter = self.stable_dt(st)
        if dt_try is not None:
            dt = min(dt, dt_try)
        rejected = 0
        last = None
        for _ in range(self.num.retry_limit):
            try:
                new, diag = self._step(st, dt)
                new.rejected = st.rejected + rejected
                diag["limiter"] = limiter
                diag["rejected"] = rejected
                self.log.append(diag)
                return new, diag
            except SimFailure as e:
                if not e.retryable:
                    raise
                last = e
                rejected += 1
                dt *= 0.5
        code = "OUT_OF_SCOPE" if last.context.get("geometry") else last.code
        raise SimFailure(code, f"retries exhausted: {last.reason}", **last.context)

    def run(self, st, t_end, callback=None, max_steps=10 ** 7):
        steps = 0
        while st.t < t_end - 1e-12 and steps < max_steps:
            st, diag = self.advance(st, dt_try=t_end - st.t)
            steps += 1
            if callback is not None:
                callback(st, diag)
        return st


# ---------------------------------------------------------------------- diagnostics
def totals(solver, st):
    m = solver.mesh
    C = st.composition()
    Vj = C * st.Vtot[None, :]
    v = st.velocity()
    g = Geometry(m, st.D)
    x_liq = C[0] + C[1]
    xM = np.where(x_liq > 0, C[1] / x_liq, 0.0)
    return {
        "V": float(st.Vtot.sum()), "V_geo": float(g.V.sum()), "Vj": Vj.sum(1), "M": float(st.M.sum()),
        "Mj": solver.rho * Vj.sum(1), "P": st.P.sum(1),
        "KE": float(0.5 * (st.M * (v ** 2).sum(0)).sum()),
        "PE": float((st.M * solver.model.g * g.zc.ravel()).sum()),
        "mixing": float((x_liq * st.Vtot * xM * (1 - xM)).sum()),
        "vmax": float(np.abs(v).max()),
    }
