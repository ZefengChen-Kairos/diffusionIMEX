"""G2: lateral wobble -> layered heart stripes (plotter relation lambda ~ v_rel / f)."""
import numpy as np

from actions.scripts import Action, Segment
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from gates.common import Report, alpha_of, img_of, step_until
from gates import metrics as mt
from render.render import board


def wobble_action(a=9.0, f=1.5, T_w=3.8, y0=-4.0, y1=-20.0, T_end=None):
    """The layered-heart drawing phase without the cut: retreat while wobbling, then stop."""
    S = Segment
    return Action(f"wobble_a{a:g}_f{f:g}", [
        S(0.0, 0.6, (0, y0), None, 11.0, None, 0.014),
        S(0.6, 0.6 + T_w, (0, y0), (0, y1), 11.0, None, 0.012, a=a, f=f),
    ], T_end or 0.6 + T_w)


def axis_stripes(a, g, x0=0.0, amp=0.15):
    yc, prof = mt.profile(a, g, x0=x0, axis="y", width=1.0e-3)
    pk = mt.stripes(prof, amp)
    ys = yc[pk]
    lam = float(np.mean(np.diff(ys))) if len(ys) >= 2 else float("nan")
    return pk, ys, lam, prof


def run_g2():
    R = Report("G2")
    P, NP = PhysParams(), NumParams()
    imgs, labels = [], []
    # G2.1 control a = 0
    s0 = Sim(wobble_action(a=0.0), P=P, NP=NP)
    s0.run()
    a0 = alpha_of(s0, s0.Yf, s0.U[0])
    pk0, *_ = axis_stripes(a0, s0.gf)
    R.check("G2.1", "wobble amplitude 0: number of stripes on the axis", len(pk0), len(pk0) == 0, "== 0")
    imgs.append(img_of(s0, s0.Yf, s0.U[0]))
    labels.append("a=0")

    # G2.2 / G2.3 frequency scan at fixed amplitude and retreat speed
    res = {}
    for f in (1.0, 1.5, 2.0):
        s = Sim(wobble_action(a=9.0, f=f), P=P, NP=NP)
        s.run()
        a = alpha_of(s, s.Yf, s.U[0])
        pk, ys, lam, prof = axis_stripes(a, s.gf)
        res[f] = dict(n=len(pk), lam_mm=lam * 1e3, ys_mm=(ys * 1e3).round(2).tolist(),
                      contrast=mt.contrast(prof, pk))
        imgs.append(img_of(s, s.Yf, s.U[0]))
        labels.append(f"f={f:g} Hz")
        if f == 1.5:
            s15 = s
    R.extra["freq_scan"] = res
    R.check("G2.2", "a=9 mm, f=1.5 Hz: number of stripes on the axis", res[1.5]["n"], res[1.5]["n"] >= 4, ">= 4")
    lams = [res[f]["lam_mm"] for f in (2.0, 1.5, 1.0)]
    ok = all(np.isfinite(lams)) and lams[0] < lams[1] < lams[2]
    R.check("G2.3", "mean stripe spacing increases with v_rel/f (f = 2, 1.5, 1 Hz) [mm]", lams, ok, "increasing")

    # G2.5 curvature: stripes are forward-convex arcs (ridge at x=0 ahead of the ridge at x = +-4 mm)
    a = alpha_of(s15, s15.Yf, s15.U[0])
    pk, ys, lam, prof = axis_stripes(a, s15.gf)
    conv = []
    for y in ys:
        best = []
        for xo in (-4e-3, 4e-3):
            yc, pr = mt.profile(a, s15.gf, x0=xo, axis="y", width=0.75e-3)
            sel = np.abs(yc - y) < 0.6 * (lam if np.isfinite(lam) else 3e-3)
            if sel.any():
                best.append(yc[sel][np.argmax(pr[sel])])
        if len(best) == 2:
            conv.append(y - 0.5 * (best[0] + best[1]))
    frac = float(np.mean(np.array(conv) > 0)) if conv else 0.0
    R.check("G2.5", "share of stripes whose axis point leads the off-axis points (forward-convex)", frac,
            frac >= 0.6, ">= 0.6")

    # G2.4 contrast decay over 10 s after the stop
    s = Sim(wobble_action(a=9.0, f=1.5, T_end=4.4 + 10.5), P=P, NP=NP)
    step_until(s, 4.4 + 0.5)
    a1 = alpha_of(s, s.Yf, s.U[0])
    pk1, _, _, p1 = axis_stripes(a1, s.gf)
    c1 = mt.contrast(p1, pk1)
    step_until(s, 4.4 + 10.5)
    a2 = alpha_of(s, s.Yf, s.U[0])
    # same stripe positions (frozen pattern): contrast on the same peaks
    _, p2 = mt.profile(a2, s.gf, x0=0.0, axis="y", width=1.0e-3)
    c2 = mt.contrast(p2, pk1)
    dec = 1 - c2 / c1 if c1 > 0 else 1.0
    R.check("G2.4", "stripe contrast loss from stop+0.5 s to stop+10.5 s", [c1, c2, dec], dec < 0.3, "< 30 %")
    imgs.append(img_of(s, s.Yf, s.U[0]))
    labels.append("f=1.5, +10 s")
    board([imgs], R.dir + "/wobble_scan.png", labels=labels, cell=230)
    R.write(P, NP)
    return R
