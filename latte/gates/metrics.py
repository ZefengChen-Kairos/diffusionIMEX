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
    """Indices of white peaks separated by dips deeper than amp (hysteresis)."""
    p = np.asarray(prof, float)
    peaks = []
    i_hi, v_hi = None, -np.inf
    v_lo = np.inf
    rising = True
    for k, v in enumerate(p):
        if rising:
            if v > v_hi:
                v_hi, i_hi = v, k
            if v_hi - v > amp and i_hi is not None and v_hi - v_lo > amp:
                peaks.append(i_hi)
                rising = False
                v_lo = v
        else:
            if v < v_lo:
                v_lo = v
            if v - v_lo > amp:
                rising = True
                v_hi, i_hi = v, k
    # a final peak that ends at the pattern edge
    if rising and i_hi is not None and v_hi - v_lo > amp and (len(p) - 1 - i_hi) > 0 and v_hi - p[-1] > amp:
        peaks.append(i_hi)
    return peaks


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
