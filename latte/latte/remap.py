"""Conservative, bounded vertical remap between layer distributions (V3 Sec. 7.4, T10).

Not used inside the stage-1 time step (layers are fixed sigma fractions, interlayer
exchange is carried by the relative interface fluxes), but provided for changing the layer
distribution and tested in isolation.

Reconstruction per old layer: piecewise linear with central slopes and ONE limiter factor
theta per cell shared by all components (Barth-Jespersen type: edge values stay within the
neighbouring cell means).  Central slopes of compositions that sum to one in every cell sum
to zero, so the reconstruction stays on the simplex pointwise; the remapped cell values are
exact integrals of that reconstruction and are therefore conservative and bounded.
"""
import numpy as np


def _limited_slopes(zc, q):
    """zc (N,), q (m, N) cell means -> limited slopes (m, N)."""
    sl = np.zeros_like(q)
    sl[:, 1:-1] = (q[:, 2:] - q[:, :-2]) / (zc[2:] - zc[:-2])[None, :]
    qmax = q.copy(); qmin = q.copy()
    qmax[:, 1:] = np.maximum(qmax[:, 1:], q[:, :-1]); qmax[:, :-1] = np.maximum(qmax[:, :-1], q[:, 1:])
    qmin[:, 1:] = np.minimum(qmin[:, 1:], q[:, :-1]); qmin[:, :-1] = np.minimum(qmin[:, :-1], q[:, 1:])
    return sl, qmax, qmin


def remap_column(z_old, z_new, comp_vol, extra=None, order=2):
    """Remap one column.

    z_old, z_new : interface heights (N_old+1,), (N_new+1,) with equal end points
    comp_vol     : (3, N_old) component volumes per unit area (thickness * fraction)
    extra        : (k, N_old) other conserved densities-times-thickness (e.g. momentum), remapped
                   with their own limiter
    Returns component volumes (3, N_new) and extra (k, N_new).
    """
    z_old = np.asarray(z_old, float); z_new = np.asarray(z_new, float)
    if abs(z_old[0] - z_new[0]) > 1e-14 * max(1, abs(z_old[-1])) or abs(z_old[-1] - z_new[-1]) > 1e-14 * max(1, abs(z_old[-1])):
        raise ValueError("remap requires identical column end points")
    h = np.diff(z_old)
    zc = 0.5 * (z_old[1:] + z_old[:-1])
    C = comp_vol / h[None, :]
    groups = [(C, True)]
    if extra is not None:
        groups.append((np.asarray(extra, float) / h[None, :], False))
    out = []
    for q, joint in groups:
        if order == 1 or q.shape[1] < 3:
            sl = np.zeros_like(q)
        else:
            sl, qmax, qmin = _limited_slopes(zc, q)
            dev = sl * (0.5 * h)[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                th_up = np.where(dev > 0, (qmax - q) / dev, np.where(dev < 0, (qmin - q) / dev, 1.0))
                th_dn = np.where(dev > 0, (q - qmin) / dev, np.where(dev < 0, (q - qmax) / dev, 1.0))
            th = np.clip(np.minimum(th_up, th_dn), 0.0, 1.0)
            if joint:
                th = np.broadcast_to(th.min(0, keepdims=True), th.shape)
            sl = sl * th
        # exact overlap integrals of the linear pieces
        lo = np.maximum(z_old[:-1][:, None], z_new[:-1][None, :])          # (N_old, N_new)
        hi = np.minimum(z_old[1:][:, None], z_new[1:][None, :])
        L = np.maximum(hi - lo, 0.0)
        mid = 0.5 * (hi + lo)
        w_const = L
        w_lin = L * (mid - zc[:, None])
        res = q @ w_const + sl @ w_lin
        out.append(res)
    return out[0], (out[1] if extra is not None else None)
