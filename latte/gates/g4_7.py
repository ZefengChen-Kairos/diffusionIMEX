"""G4 push heart, G5 leaf, G6 tulip, G7 swan."""
import math
import numpy as np

from actions.scripts import ALL
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


def front_edge(sh, Y, M, k, q=0.98):
    """y of the front of tag k's white region (q-quantile of y over its pixels)."""
    m = tag_mask(sh, Y, M, k)
    if not m.any():
        return float("nan")
    return float(np.quantile(sh.gf.Y[m], q))


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
    pa = ALL["push_heart"]()
    kA, kB = pa.tag("A"), pa.tag("B")
    T1, T2 = round(pa.end("A") + 0.25, 2), round(pa.end("B") + 0.25, 2)   # after part A / after part B
    Y1, M1, t1 = snapshot(r, T1)
    Y2, M2, t2 = snapshot(r, T2)
    cA2 = tag_centroid(sh, Y2, M2, kA)
    # the crescent's arms wrap backwards around B, so the push is measured on A's front edge
    fA1 = front_edge(sh, Y1, M1, kA)
    fA2 = front_edge(sh, Y2, M2, kA)
    d = float((fA2 - fA1) * 1e3)
    R.check("G4.1", f"dot B pushes A forward: front edge of A, t={T1} -> {T2} s [mm]", d, d > 2.0, "> 2 mm")
    a2 = alpha_of(sh, Y2, M2)
    mA = tag_mask(sh, Y2, M2, kA, a2)
    cB = tag_centroid(sh, Y2, M2, kB)
    cov = angular_coverage(mA, sh.gf, cB)
    R.check("G4.2", "A is a crescent around B: angular coverage of A seen from B's centroid [deg]", cov,
            cov > 180.0, "> 180")
    # brown separation band between A and B along the axis: A wraps around B as a crescent, so the band
    # is looked for on the axis from B's centroid forward to A's front edge
    yc, prof = mt.profile(a2, sh.gf, x0=0.0, axis="y", width=1.0e-3)
    yA = fA2
    yB = cB[1]
    sel = (yc > min(yA, yB)) & (yc < max(yA, yB))
    band = float(prof[sel].min()) if sel.any() else 1.0
    # band width: contiguous alpha < 0.6 between them
    width = float((prof[sel] < 0.6).sum() * sh.gf.dx * 1e3)
    R.check("G4.3", "brown band between B (centroid) and A's front edge on the axis: min alpha, width [mm]", [band, width],
            band < 0.6 and width > 0, "alpha < 0.6, width > 0")
    # control: the same pour, but the flow runs out at the end of part A (no part B)
    ctrl = ALL["push_heart"]().truncated(pa.end("A"), T_end=T2, name="push_heart_ctrl")
    c = Sim(ctrl, P=P, NP=NP, record_every=0.25)
    c.run()
    f1 = min(c.frames, key=lambda f: abs(f[0] - T1))
    f2 = c.frames[-1]
    shc = shell(1, P, NP)
    dc = float((front_edge(c, f2[1], f2[2], kA) - front_edge(c, f1[1], f1[2], kA)) * 1e3)
    R.check("G4.4", "control (flow runs out after part A): front edge of A over the same interval [mm], ratio", [dc, dc / d],
            abs(dc) < 0.2 * abs(d), "< 20 % of G4.1")
    e = tag_ledger_err(r, kA)
    R.check("G4.5", "tag A ledger: surface + bulk (incl. E-subducted, minus R) = injected", e, e < 1e-9, "< 1e-9")
    board([[img_of(sh, Y1, M1), img_of(sh, Y2, M2), tagimg_of(sh, Y2), img_of(c, f2[1], f2[2])]],
          R.dir + "/push.png", labels=[f"t={T1} (after A)", f"t={T2} (A pushed by B)", "tags", f"control t={T2}"],
          cell=230)
    R.write(P, NP)
    return R


# ---------------------------------------------------------------------------- G5
def leaf_metrics(sh, Y, M, yr=None):
    """yr: optional (y_lo, y_hi) window [m] for the rib count / spacing (e.g. the retreat zone)."""
    a = alpha_of(sh, Y, M)
    g = sh.gf
    w = a > 0.5
    ys = g.Y[w]
    y0, y1 = ys.min(), ys.max()
    if yr is not None:
        y0, y1 = max(y0, yr[0]), min(y1, yr[1])
    # ribs along a line parallel to the axis on one side
    ribs = {}
    for xo in (5e-3, 7e-3, 9e-3):
        yc, prof = mt.profile(a, g, x0=xo, axis="y", width=0.75e-3)
        sel = (yc > y0) & (yc < y1)
        pk = mt.stripes(prof[sel], 0.15)
        yy = yc[sel][pk]
        ribs[xo] = dict(n=mt.n_stripes(pk), spacing=float(np.mean(np.diff(yy))) if len(pk) > 1 else float("nan"))
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
    base = ALL["leaf"]()
    r0, r1 = base.marks["retreat0"], base.marks["retreat1"]
    shift = (r1 - r0) * (1 - 1 / 1.5)
    act = base.time_warped(r0, r1, 1.5, name="leaf_fast_retreat")
    sim = Sim(act, P=P, NP=NP)
    # compare before the pull-through: rib spacing at the end of the retreat
    from gates.common import step_until
    step_until(sim, r1 - shift + 0.2)
    Yb, Mb, _ = snapshot(r, r1 + 0.2)
    # ribs laid down during the retreat: between the spout positions at the end and the start of the retreat
    yr = (base.pos_mm(r1)[1] * 1e-3, base.pos_mm(r0)[1] * 1e-3)
    m_ref = leaf_metrics(sh, Yb, Mb, yr)
    m_fast = leaf_metrics(sim, sim.Yf, sim.U[0], yr)
    R.extra["G5.3"] = dict(ref=m_ref, fast=m_fast)
    R.check("G5.3", "rib spacing in the retreat zone before the pull-through: reference vs retreat x1.5 [mm]",
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
    ta = ALL["tulip"]()
    kA, kB, kC = ta.tag("A"), ta.tag("B"), ta.tag("C")
    tA, tB, tC = ta.end("A") + 0.2, ta.end("B") + 0.2, ta.end("C") + 0.1
    Y, M, _ = snapshot(r, tB)           # after the second (pressing) segment
    Y0, M0, _ = snapshot(r, tA)         # end of the first segment
    share = Y[1 + kA] / np.maximum(Y[0], 1e-12)
    a = alpha_of(sh, Y, M)
    mA_white = (a > 0.5) & (share > 0.3)
    _, n_white = mt.components(mA_white, 30)
    # material of batch A (not its whiteness): Y_A > 0.2
    _, n_mat = mt.components(Y[1 + kA] > 0.2, 30)

    def massA(Yx, Mx):
        Mf = sh.fm.to_fine(sh.g.extrapolate(Mx.copy()))
        return float((Yx[1 + kA] * Mf * sh.gf.wet).sum() * sh.gf.dA)
    keep = massA(Y, M) / massA(Y0, M0)
    R.extra["G6.1_white_components_of_A"] = n_white
    R.check("G6.1a", f"batch A is pressed open, not erased and redrawn: surface A mass kept {tA:.1f} -> {tB:.1f} s; "
            "components of A material (Y_A > 0.2)", [keep, n_mat], keep >= 0.8 and 1 <= n_mat <= 2,
            "kept >= 80 %, 1..2 pieces")
    e = tag_ledger_err(r, kA)
    R.check("G6.1b", "tag A ledger closure", e, e < 1e-9, "< 1e-9")
    Y6, M6, _ = snapshot(r, tC)         # all three layers, before the cut
    a6 = alpha_of(sh, Y6, M6)
    yc, prof = mt.profile(a6, sh.gf, x0=0.0, axis="y", width=1.5e-3)
    # layer bands = groups of white peaks separated by deep brown gaps
    pk = mt.stripes(prof, 0.35)
    R.check("G6.2", "white bands on the axis before the pull-through (one per part A, B, C)", len(pk), len(pk) == 3, "== 3")
    # each band dominated by a different tag
    i0 = int(np.argmin(np.abs(sh.gf.xc)))
    doms = [ta.labels()[[kA, kB, kC][int(np.argmax([Y6[1 + k][i, i0] for k in (kA, kB, kC)]))]] for i in pk]
    R.extra["G6.2_band_tags"] = doms
    board([[img_of(sh, Y, M), tagimg_of(sh, Y), img_of(sh, Y6, M6), tagimg_of(sh, Y6),
            img_of(sh, r["Yfinal"].astype(float), r["Mfinal"].astype(float))]],
          R.dir + "/tulip.png", labels=[f"t={tB:.1f}", "tags", f"t={tC:.1f}", "tags", "final"], cell=210)
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
    act = ALL["swan"]()
    wing = (a > 0.5) & (Y[1 + act.tag("wing")] / np.maximum(Y[0], 1e-12) > 0.5)
    neck = tag_mask(sh, Y, M, act.tag("neck"), a)
    head = tag_mask(sh, Y, M, act.tag("head"), a)
    areas = [mt.area(x, g) * 1e6 for x in (wing, neck, head)]
    R.check("G7.1", "wing / neck / head white areas [mm^2]", areas, min(areas) > 20.0, "each > 20 mm^2")
    # neck width = area / centre-line length of the neck stroke
    L_neck = act.meta["neck_len_mm"] * 1e-3
    neck_w = mt.area(neck, g) / L_neck
    # wing rib spacing along the wing axis
    ang = math.radians(act.meta["wing_ang_deg"])
    ax = np.array([math.sin(ang), math.cos(ang)])
    P0 = np.array(act.meta["wing_P0"]) * 1e-3
    s_ = np.linspace(act.meta["wing_s"][0] * 1e-3, act.meta["wing_s"][1] * 1e-3, 400)
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
