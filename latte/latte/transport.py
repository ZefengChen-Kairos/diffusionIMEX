"""Conservative transport of volume, components and momentum with one set of face fluxes
(V3 Sec. 4.3, 4.4, 7.2).

Inputs are time-integrated *relative* volume fluxes dV through every internal face
(x faces, y faces, interior interfaces; positive from the minus cell to the plus cell),
plus boundary inflows (surface inlet, or prescribed boundary faces in tests) and
boundary outflows.  Component face volumes are built from the same dV:

  low order   : donor-cell composition (bounded, simplex preserving under CFL <= 1)
  high order  : Lax-Wendroff face composition  C_d + (1 - nu)/2 (C_a - C_d)
  limiter     : Zalesak FCT on each of B, M, G; ONE factor per face = min over the three
                components (joint limiter: antidiffusive component fluxes sum to zero, so
                the face volume, the simplex and every bound are preserved; no clamping).

The mass flux is F_rho = sum_j rho_j F_j and the momentum flux is F_rho * v_f with a
minmod-limited upwind face velocity, so a uniform velocity is preserved exactly across
density jumps.
"""
import numpy as np

from .failures import SimFailure


def _net_out(mesh, f):
    """Net outflow per cell of a face quantity f (positive from tf_cm to tf_cp)."""
    nc = mesh.nc
    return (np.bincount(mesh.tf_cm, weights=f, minlength=nc)
            - np.bincount(mesh.tf_cp, weights=f, minlength=nc))


def _minmod(a, b):
    return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


class TransportResult:
    pass


def _phi(r, kind):
    if kind == "superbee":
        return np.maximum(0.0, np.maximum(np.minimum(2 * r, 1.0), np.minimum(r, 2.0)))
    if kind == "vanleer":
        return (r + np.abs(r)) / (1 + np.abs(r))
    raise ValueError(kind)


def transport(mesh, Vc, Vtot, P, M, dV, rho_j, inflow=None, outflow=None,
              scheme="fct", momentum_recon="minmod", bound_tol=1e-12, ho_limiter="superbee", dc_floor=1e-12):
    """One explicit transport step.

    Vc     : (2, nc) component volumes of B and M (G is derived: Vtot - V_B - V_M)
    Vtot   : (nc,) material volume
    P, M   : (3, nc) momentum and (nc,) mass
    dV     : (n_tf,) relative face volumes for this step
    inflow : dict(cell, vol, comp (3,n), mass, mom (3,n)) volumes entering through boundaries
    outflow: dict(cell, vol) volumes leaving through open boundaries (donor = cell)
    """
    nc = mesh.nc
    cm, cp = mesh.tf_cm, mesh.tf_cp
    rho_j = np.asarray(rho_j, float)
    VG = Vtot - Vc[0] - Vc[1]
    V3 = np.stack([Vc[0], Vc[1], VG])
    C = V3 / Vtot[None, :]

    pos = dV > 0
    donor = np.where(pos, cm, cp)
    accept = np.where(pos, cp, cm)
    adV = np.abs(dV)

    in_vol = np.zeros(nc)
    in_comp = np.zeros((3, nc))
    in_mass = np.zeros(nc)
    in_mom = np.zeros((3, nc))
    if inflow is not None and len(inflow["cell"]):
        c = inflow["cell"]
        np.add.at(in_vol, c, inflow["vol"])
        for j in range(3):
            np.add.at(in_comp[j], c, inflow["comp"][j])
        np.add.at(in_mass, c, inflow["mass"])
        for d in range(3):
            np.add.at(in_mom[d], c, inflow["mom"][d])
    out_vol = np.zeros(nc)
    if outflow is not None and len(outflow["cell"]):
        np.add.at(out_vol, outflow["cell"], outflow["vol"])

    # positivity of the low-order update (donor CFL): total outgoing volume <= V^n
    out_total = (np.bincount(donor, weights=adV, minlength=nc) + out_vol)
    if np.any(out_total > Vtot * (1 + 1e-12)):
        i = int(np.argmax(out_total / Vtot))
        raise SimFailure("BOUNDS_FAILED", "transport CFL exceeded (outflow > cell volume)",
                         cell=i, ratio=float(out_total[i] / Vtot[i]))

    Vnew = Vtot - _net_out(mesh, dV) + in_vol - out_vol
    if np.any(Vnew <= 0):
        raise SimFailure("BOUNDS_FAILED", "non-positive cell volume after transport")

    FL = dV[None, :] * C[:, donor]                                       # (3, nf)
    VL = (V3 - np.stack([_net_out(mesh, FL[j]) for j in range(3)])
          + in_comp - out_vol[None, :] * C)
    CL = VL / Vnew[None, :]

    if scheme == "donor":
        F = FL
    elif scheme in ("fct", "ho_unlimited"):
        nu = np.minimum(adV / Vtot[donor], 1.0)
        dC = C[:, accept] - C[:, donor]
        if ho_limiter == "lw":
            phi = np.ones(dV.size)
        else:
            # flux-limited Lax-Wendroff; ONE phi per face (min over components) keeps sum_j A_j = 0
            upup = np.where(pos, mesh.tf_um, mesh.tf_up)
            Cu = np.where((upup >= 0)[None, :], C[:, np.maximum(upup, 0)], C[:, donor])
            num_ = C[:, donor] - Cu
            # components with a negligible jump (incl. round-off of the derived component) do not
            # constrain the joint limiter; FCT below still guarantees their bounds
            sig = np.abs(dC) > dc_floor
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(sig, num_ / np.where(sig, dC, 1.0), 0.0)
            ph = _phi(r, ho_limiter)
            ph = np.where(sig, ph, np.inf)
            phi = np.min(ph, axis=0)
            phi = np.where(np.isfinite(phi), phi, 0.0)
        Aj = dV[None, :] * 0.5 * ((1.0 - nu) * phi)[None, :] * dC                      # (3, nf)
        # admissible bounds: old and low-order values over the face neighbourhood (+ inflow state)
        cmax = np.maximum(C, CL)
        cmin = np.minimum(C, CL)
        nbmax, nbmin = cmax.copy(), cmin.copy()
        for j in range(3):
            np.maximum.at(nbmax[j], cm, cmax[j, cp])
            np.maximum.at(nbmax[j], cp, cmax[j, cm])
            np.minimum.at(nbmin[j], cm, cmin[j, cp])
            np.minimum.at(nbmin[j], cp, cmin[j, cm])
        if inflow is not None and len(inflow["cell"]):
            Cin = inflow["comp"] / np.maximum(inflow["vol"], 1e-300)[None, :]
            for j in range(3):
                np.maximum.at(nbmax[j], inflow["cell"], Cin[j])
                np.minimum.at(nbmin[j], inflow["cell"], Cin[j])
        np.clip(nbmin, 0.0, 1.0, out=nbmin)
        np.clip(nbmax, 0.0, 1.0, out=nbmax)
        if scheme == "ho_unlimited":          # diagnostic only: high-order flux without any limiting
            F = FL + Aj
            lam_tot = np.ones(dV.size)
        else:
            # standard Zalesak factor (joint over components): bounded on its own
            lam_z = np.ones(dV.size)
            for j in range(3):
                a = Aj[j]
                ap, an = np.maximum(a, 0.0), np.maximum(-a, 0.0)
                Pp = np.bincount(cp, weights=ap, minlength=nc) + np.bincount(cm, weights=an, minlength=nc)
                Pm = np.bincount(cm, weights=ap, minlength=nc) + np.bincount(cp, weights=an, minlength=nc)
                Qp = np.maximum(nbmax[j] - CL[j], 0.0) * Vnew
                Qm = np.maximum(CL[j] - nbmin[j], 0.0) * Vnew
                Rp = np.where(Pp > 0, np.minimum(1.0, Qp / np.where(Pp > 0, Pp, 1.0)), 1.0)
                Rm = np.where(Pm > 0, np.minimum(1.0, Qm / np.where(Pm > 0, Pm, 1.0)), 1.0)
                lj = np.where(a > 0, np.minimum(Rp[cp], Rm[cm]), np.where(a < 0, np.minimum(Rp[cm], Rm[cp]), 1.0))
                lam_z = np.minimum(lam_z, lj)
            # selective FCT: keep the (TVD-limited) high-order flux except on faces of cells that
            # would leave the admissible bounds; extend the limited set until no cell violates.
            # Worst case all faces carry lam_z, i.e. standard Zalesak, which is bounded.
            limited = np.zeros(dV.size, dtype=bool)
            tol = 0.1 * bound_tol
            for _ in range(64):
                lam_tot = np.where(limited, lam_z, 1.0)
                F = FL + lam_tot[None, :] * Aj
                Vt = VL - np.stack([_net_out(mesh, lam_tot * Aj[j]) for j in range(3)])
                Ct = Vt / Vnew[None, :]
                bad = np.any((Ct < nbmin - tol) | (Ct > nbmax + tol), axis=0)
                if not bad.any():
                    break
                newl = bad[cm] | bad[cp]
                if np.all(limited | ~newl):
                    limited[:] = True
                else:
                    limited |= newl
            else:
                lam_tot = lam_z                  # fallback: full Zalesak (bounded)
                F = FL + lam_tot[None, :] * Aj
    else:
        raise ValueError(scheme)

    Vn3 = (V3 - np.stack([_net_out(mesh, F[j]) for j in range(3)])
           + in_comp - out_vol[None, :] * C)
    VB, VM = Vn3[0], Vn3[1]
    VGd = Vnew - VB - VM                     # derived third component (V3 Sec. 4.2)
    Cn = np.stack([VB, VM, VGd]) / Vnew[None, :]
    lo, hi = Cn.min(), Cn.max()
    if lo < -bound_tol or hi > 1 + bound_tol:
        raise SimFailure("BOUNDS_FAILED", "component fraction out of [0, 1] after transport",
                         min=float(lo), max=float(hi))

    # mass and momentum with the same component fluxes
    Frho = rho_j @ F
    rho_c = M / Vtot
    Mnew = M - _net_out(mesh, Frho) + in_mass - out_vol * rho_c
    v = P / M[None, :]
    vd = v[:, donor]
    if momentum_recon == "minmod":
        upup = np.where(pos, mesh.tf_um, mesh.tf_up)
        has = upup >= 0
        vu = np.where(has[None, :], v[:, np.maximum(upup, 0)], vd)
        vf = vd + 0.5 * _minmod(v[:, accept] - vd, vd - vu)
    elif momentum_recon == "upwind":
        vf = vd
    else:
        raise ValueError(momentum_recon)
    Pnew = (P - np.stack([_net_out(mesh, Frho * vf[d]) for d in range(3)])
            + in_mom - (out_vol * rho_c)[None, :] * v)

    r = TransportResult()
    r.Vc = np.stack([VB, VM])
    r.Vtot = Vnew
    r.M = Mnew
    r.P = Pnew
    r.face_comp = F
    r.face_mass = Frho
    r.in_vol, r.in_comp, r.in_mass, r.in_mom = in_vol, in_comp, in_mass, in_mom
    r.out_vol = out_vol
    r.out_comp = out_vol[None, :] * C
    r.out_mass = out_vol * rho_c
    r.out_mom = (out_vol * rho_c)[None, :] * v
    r.mass_consistency = float(np.max(np.abs(Mnew - rho_j @ np.stack([VB, VM, VGd]))) / max(np.max(np.abs(Mnew)), 1e-300))
    r.cfl = float(np.max(out_total / Vtot))
    r.fct_lambda = lam_tot if scheme != 'donor' else None
    r.fct_limited_fraction = float(np.mean(lam_tot < 1.0)) if scheme != 'donor' else None
    r.antidiffusive = Aj if scheme != 'donor' else None
    return r
