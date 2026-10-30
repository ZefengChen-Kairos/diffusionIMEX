"""Pattern metrics on the fine material field."""
import numpy as np
from scipy import ndimage

from render.render import alpha


def alpha_field(sim, Yf=None, M=None):
    phi, _ = sim.phi_fine(Yf, M)
    return alpha(phi, sim.P.Phi_star) * sim.gf.wet


def white(sim, Yf=None, M=None, thr=0.5):
    return alpha_field(sim, Yf, M) > thr


def area(mask, grid):
    return float(mask.sum() * grid.dA)


def centroid(w, grid):
    s = w.sum()
    if s <= 0:
        return np.array([np.nan, np.nan])
    return np.array([(w * grid.X).sum() / s, (w * grid.Y).sum() / s])


def iou(a, b):
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else 1.0


def mirror_iou(mask):
    return iou(mask, mask[:, ::-1])


def boundary(mask):
    return mask & ~ndimage.binary_erosion(mask)


def boundary_displacement(m0, m1, grid):
    """Hausdorff-type max distance between the white boundaries of two masks [m]."""
    b0, b1 = boundary(m0), boundary(m1)
    if not b0.any() or not b1.any():
        return np.inf
    d1 = ndimage.distance_transform_edt(~b1) * grid.dx
    d0 = ndimage.distance_transform_edt(~b0) * grid.dx
    return float(max(d1[b0].max(), d0[b1].max()))


def mean_boundary_displacement(m0, m1, grid):
    b0, b1 = boundary(m0), boundary(m1)
    if not b0.any() or not b1.any():
        return np.inf
    d1 = ndimage.distance_transform_edt(~b1) * grid.dx
    return float(d1[b0].mean())


def profile(a, grid, x0=0.0, axis="y", width=0.0015):
    """Mean alpha along a line (x = x0 for axis 'y'), averaged over +-width."""
    if axis == "y":
        sel = np.abs(grid.xc - x0) <= width
        return grid.xc, a[:, sel].mean(axis=1)
    sel = np.abs(grid.xc - x0) <= width
    return grid.xc, a[sel, :].mean(axis=0)


def stripes(prof, amp=0.15):
    """Indices of white peaks that rise > amp above the trough on their left and fall > amp
    on their right (hysteresis peak finder)."""
    p = np.asarray(prof, float)
    peaks = []
    if p.size == 0:
        return peaks
    looking = True
    run_min = p[0]
    cand_v, cand_i, left = -np.inf, None, p[0]
    cur_min = p[0]
    for k, v in enumerate(p):
        if looking:
            if v > cand_v:
                cand_v, cand_i, left = v, k, run_min
            run_min = min(run_min, v)
            if cand_i is not None and cand_v - left > amp and cand_v - v > amp:
                peaks.append(cand_i)
                looking = False
                cur_min = v
        else:
            cur_min = min(cur_min, v)
            if v - cur_min > amp:
                looking = True
                run_min = cur_min
                cand_v, cand_i, left = v, k, cur_min
    return peaks


def n_stripes(peaks):
    """A single white body is not a stripe pattern: stripes need at least two peaks."""
    return len(peaks) if len(peaks) >= 2 else 0


def components(mask, min_area_px=20):
    lab, n = ndimage.label(mask)
    if n == 0:
        return lab, 0
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    return lab, int((sizes >= min_area_px).sum())


def contrast(prof, peaks):
    """Mean peak-to-neighbouring-trough contrast of the stripes."""
    p = np.asarray(prof)
    if len(peaks) < 2:
        return 0.0
    c = []
    for a, b in zip(peaks[:-1], peaks[1:]):
        c.append(0.5 * (p[a] + p[b]) - p[a:b + 1].min())
    return float(np.mean(c))
