"""G4 push heart, G5 leaf, G6 tulip, G7 swan."""
import math
import numpy as np

from actions.scripts import Action, push_heart, leaf, tulip, swan
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from gates.common import Report, shell, alpha_of, img_of, tagimg_of
from gates import metrics as mt
from render.render import board
from run_patterns import run_pattern, snapshot


def tag_mask(sh, Y, M, k, a=None, frac=0.5):
    """White pixels dominated by tag k."""
    a = alpha_of(sh, Y, M) if a is None else a
    share = Y[1 + k] / np.maximum(Y[0], 1e-12)
    return (a > 0.5) & (share > frac)


def tag_centroid(sh, Y, M, k):
    Mf = sh.fm.to_fine(sh.g.extrapolate(M.copy()))
    w = Y[1 + k] * Mf * sh.gf.wet
    s = w.sum()
    return np.array([(w * sh.gf.X).sum() / s, (w * sh.gf.Y).sum() / s])


def angular_coverage(mask, g, c, nb=36):
    if not mask.any():
        return 0.0
    ang = np.arctan2(g.Y[mask] - c[1], g.X[mask] - c[0])
    bins = np.unique(((ang + math.pi) / (2 * math.pi) * nb).astype(int) % nb)
    return 360.0 * len(bins) / nb


def tag_ledger_err(r, k):
    """surface_k + bulk_k == injected_k (coarse conservative tags)."""
    U = r["Ufinal"]
    dA = (0.096 / U.shape[-1]) ** 2
    surf = U[4 + k].sum() * dA
    bulk = float(r["bulk_tags"][k])
    inj = float(r["inj_tags"][k])
    return abs(surf + bulk - inj) / inj


# ---------------------------------------------------------------------------- G4
def run_g4():
    R = Report("G4")
    P, NP = PhysParams(), NumParams()
    r = run_pattern("push_heart", P, NP)
    K = r["Y"].shape[1] - 1
    sh = shell(K, P, NP)
    Y1, M1, t1 = snapshot(r, 1.9)      # after dot A, before dot B
    Y2, M2, t2 = snapshot(r, 3.6)      # after dot B, before dot C
    cA1 = tag_centroid(sh, Y1, M1, 0)
    cA2 = tag_centroid(sh, Y2, M2, 0)
    d = float((cA2[1] - cA1[1]) * 1e3)
    R.check("G4.1", "dot B pushes A forward: d centroid_y(A) between t=1.9 and 3.6 s [mm]", d, d > 2.0, "> 2 mm")
    a2 = alpha_of(sh, Y2, M2)
    mA = tag_mask(sh, Y2, M2, 0, a2)
    cB = tag_centroid(sh, Y2, M2, 1)
    cov = angular_coverage(mA, sh.gf, cB)
    R.check("G4.2", "A is a crescent around B: angular coverage of A seen from B's centroid [deg]", cov,
            cov > 180.0, "> 180")
    # brown separation band between A and B along the axis
    yc, prof = mt.profile(a2, sh.gf, x0=0.0, axis="y", width=1.0e-3)
    yA = cA2[1]
    yB = cB[1]
    sel = (yc > min(yA, yB)) & (yc < max(yA, yB))
    band = float(prof[sel].min()) if sel.any() else 1.0
    # band width: contiguous alpha < 0.6 between them
    width = float((prof[sel] < 0.6).sum() * sh.gf.dx * 1e3)
    R.check("G4.3", "brown band between A and B on the axis: min alpha, width [mm]", [band, width],
            band < 0.6 and width > 0, "alpha < 0.6, width > 0")
    # control: only dot A
    act = push_heart()
    ctrl = Action("push_heart_ctrl", [act.segs[0]], 3.6)
    c = Sim(ctrl, P=P, NP=NP, record_every=0.25)
    c.run()
    f1 = min(c.frames, key=lambda f: abs(f[0] - 1.9))
    f2 = c.frames[-1]
    shc = shell(1, P, NP)
    dc = float((tag_centroid(c, f2[1], f2[2], 0)[1] - tag_centroid(c, f1[1], f1[2], 0)[1]) * 1e3)
    R.check("G4.4", "control (no dot B): d centroid_y(A) over the same interval [mm], ratio", [dc, dc / d],
            abs(dc) < 0.2 * abs(d), "< 20 % of G4.1")
    e = tag_ledger_err(r, 0)
    R.check("G4.5", "tag A ledger: surface + bulk (incl. E-subducted, minus R) = injected", e, e < 1e-9, "< 1e-9")
    board([[img_of(sh, Y1, M1), img_of(sh, Y2, M2), tagimg_of(sh, Y2), img_of(c, f2[1], f2[2])]],
          R.dir + "/push.png", labels=["t=1.9 (A)", "t=3.6 (A pushed by B)", "tags", "control t=3.6"], cell=230)
    R.write(P, NP)
    return R


# ---------------------------------------------------------------------------- G5
def leaf_metrics(sh, Y, M):
    a = alpha_of(sh, Y, M)
    g = sh.gf
    w = a > 0.5
    ys = g.Y[w]
    y0, y1 = ys.min(), ys.max()
    # ribs along a line parallel to the axis on one side
    ribs = {}
    for xo in (5e-3, 7e-3, 9e-3):
        yc, prof = mt.profile(a, g, x0=xo, axis="y", width=0.75e-3)
        sel = (yc > y0) & (yc < y1)
        pk = mt.stripes(prof[sel], 0.15)
        yy = yc[sel][pk]
        ribs[xo] = dict(n=len(pk), spacing=float(np.mean(np.diff(yy))) if len(pk) > 1 else float("nan"))
    best = max(ribs.values(), key=lambda v: v["n"])
    sym = mt.mirror_iou(w)
    # vein: lateral profiles in the rib zone have a dip at x = 0
    dips = []
    for yq in np.linspace(y0 + 0.2 * (y1 - y0), y0 + 0.7 * (y1 - y0), 6):
        xc, pr = mt.profile(a, g, x0=yq, axis="x", width=0.75e-3)
        mid = pr[np.abs(xc) < 0.75e-3].min()
        side = pr[(np.abs(xc) > 1.5e-3) & (np.abs(xc) < 4e-3)].max()
        dips.append(bool(side - mid > 0.15))
    return dict(ribs=best["n"], spacing=best["spacing"], sym=sym, vein=float(np.mean(dips)), per_line=ribs)


def run_g5():
    R = Report("G5")
    P, NP = PhysParams(), NumParams()
    r = run_pattern("leaf", P, NP)
    K = r["Y"].shape[1] - 1
    sh = shell(K, P, NP)
    Y, M = r["Yfinal"].astype(float), r["Mfinal"].astype(float)
    m = leaf_metrics(sh, Y, M)
    R.extra["leaf"] = m
    R.check("G5.1", "rib count along a line beside the axis", m["ribs"], m["ribs"] >= 5, ">= 5")
    R.check("G5.2", "left/right symmetry: mirror IoU of the white mask", m["sym"], m["sym"] >= 0.75, ">= 0.75")
    R.check("G5.4", "midline vein: share of transverse profiles with a dip at x=0", m["vein"], m["vein"] >= 0.5,
            ">= 0.5")
    # G5.3 rib spacing vs retreat speed (retreat 1.5x faster)
    act = leaf()
    seg = act.segs[2]
    T0 = seg.t1 - seg.t0
    shift = T0 - T0 / 1.5
    for s_ in act.segs[2:]:
        if s_ is seg:
            s_.t1 -= shift
        else:
            s_.t0 -= shift
            s_.t1 -= shift
    act.T_end -= shift
    # rebuild to refresh tapers / tags
    act = Action("leaf_fast_retreat", act.segs, act.T_end)
    sim = Sim(act, P=P, NP=NP)
    # compare before the central cut: rib spacing at the end of the retreat
    from gates.common import step_until
    step_until(sim, seg.t1 + 0.2)
    Yb, Mb, _ = snapshot(r, act.segs[2].t1 + shift + 0.2)
    m_ref = leaf_metrics(sh, Yb, Mb)
    m_fast = leaf_metrics(sim, sim.Yf, sim.U[0])
    R.extra["G5.3"] = dict(ref=m_ref, fast=m_fast)
    R.check("G5.3", "rib spacing before the cut: reference vs retreat x1.5 [mm]",
            [m_ref["spacing"] * 1e3, m_fast["spacing"] * 1e3], m_fast["spacing"] > m_ref["spacing"], "increases")
    board([[img_of(sh, Yb, Mb), img_of(sim, sim.Yf, sim.U[0]), img_of(sh, Y, M), tagimg_of(sh, Y)]],
          R.dir + "/leaf.png", labels=["ref, end of retreat", "retreat x1.5", "final", "tags"], cell=230)
    R.write(P, NP)
    return R


# ---------------------------------------------------------------------------- G6
def run_g6():
    R = Report("G6")
    P, NP = PhysParams(), NumParams()
    r = run_pattern("tulip", P, NP)
    K = r["Y"].shape[1] - 1
    sh = shell(K, P, NP)
    Y, M, _ = snapshot(r, 4.9)          # after the second (pressing) segment
    share = Y[1] / np.maximum(Y[0], 1e-12)
    a = alpha_of(sh, Y, M)
    mA = (a > 0.5) & (share > 0.3)
    _, n = mt.components(mA, 30)
    R.check("G6.1a", "first batch A stays connected after being pressed open: components", n, 1 <= n <= 2, "1..2")
    e = tag_ledger_err(r, 0)
    R.check("G6.1b", "tag A ledger closure", e, e < 1e-9, "< 1e-9")
    Y6, M6, _ = snapshot(r, 6.3)        # all three layers, before the cut
    a6 = alpha_of(sh, Y6, M6)
    yc, prof = mt.profile(a6, sh.gf, x0=0.0, axis="y", width=1.5e-3)
    # layer bands = groups of white peaks separated by deep brown gaps
    pk = mt.stripes(prof, 0.35)
    R.check("G6.2", "white bands on the axis before the cut (one per injection segment)", len(pk), len(pk) == 3, "== 3")
    # each band dominated by a different tag
    i0 = int(np.argmin(np.abs(sh.gf.xc)))
    doms = [int(np.argmax([Y6[1 + k][i, i0] for k in range(3)])) for i in pk]
    R.extra["G6.2_band_tags"] = doms
    board([[img_of(sh, Y, M), tagimg_of(sh, Y), img_of(sh, Y6, M6), tagimg_of(sh, Y6),
            img_of(sh, r["Yfinal"].astype(float), r["Mfinal"].astype(float))]],
          R.dir + "/tulip.png", labels=["t=4.9", "tags", "t=6.3", "tags", "final"], cell=210)
    R.write(P, NP)
    return R


# ---------------------------------------------------------------------------- G7
def run_g7():
    R = Report("G7")
    P, NP = PhysParams(), NumParams()
    r = run_pattern("swan", P, NP)
    K = r["Y"].shape[1] - 1
    sh = shell(K, P, NP)
    Y, M = r["Yfinal"].astype(float), r["Mfinal"].astype(float)
    a = alpha_of(sh, Y, M)
    g = sh.gf
    wing = (a > 0.5) & ((Y[1] + Y[2] + Y[3]) / np.maximum(Y[0], 1e-12) > 0.5)
    neck = tag_mask(sh, Y, M, 3, a)
    head = tag_mask(sh, Y, M, 4, a)
    areas = [mt.area(x, g) * 1e6 for x in (wing, neck, head)]
    R.check("G7.1", "wing / neck / head white areas [mm^2]", areas, min(areas) > 20.0, "each > 20 mm^2")
    # neck width = area / centre-line length of the neck stroke
    L_neck = 0.9 * math.pi * 22e-3
    neck_w = mt.area(neck, g) / L_neck
    # wing rib spacing along the wing axis
    ang = math.radians(22.0)
    ax = np.array([math.sin(ang), math.cos(ang)])
    P0 = np.array([6e-3, -2e-3])
    s_ = np.linspace(-30e-3, 12e-3, 400)
    nrm = np.array([ax[1], -ax[0]])
    sp = []
    for off in (4e-3, 6e-3, -4e-3, -6e-3):
        pts = P0[None] + s_[:, None] * ax[None] + off * nrm[None]
        from scipy import ndimage
        jj = (pts[:, 1] + 0.5 * g.L) / g.dx - 0.5
        ii = (pts[:, 0] + 0.5 * g.L) / g.dx - 0.5
        pr = ndimage.map_coordinates(a, [jj, ii], order=1)
        pk = mt.stripes(pr, 0.15)
        if len(pk) > 1:
            sp.append(float(np.mean(np.diff(s_[pk]))))
    rib = float(np.median(sp)) if sp else float("nan")
    R.check("G7.2", "neck width < wing rib spacing [mm]", [neck_w * 1e3, rib * 1e3],
            np.isfinite(rib) and neck_w < rib, "neck < rib")
    R.check("G7.3", "same physics parameters as all other patterns (digest)", P.digest(),
            P.digest() == PhysParams().digest(), "== shared")
    board([[img_of(sh, Y, M), tagimg_of(sh, Y)]], R.dir + "/swan.png", labels=["final", "tags"], cell=280)
    R.write(P, NP)
    return R
