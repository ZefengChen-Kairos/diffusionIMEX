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
    # median spacing: robust to a single gap where stripes are unresolved (the mean is reported too)
    lam = float(np.median(np.diff(ys))) if len(ys) >= 2 else float("nan")
    return pk, ys, lam, prof


def ridge_drop(a, g, j0, x_end=8e-3):
    """Follow a stripe ridge from the axis (row j0) sideways; return y(x_end) - y(0) for both sides
    (negative = the stripe bends backwards = forward-convex arc), or None if the ridge is lost."""
    i0 = int(np.argmin(np.abs(g.xc)))
    n = int(round(x_end / g.dx))
    out = []
    for sgn in (-1, 1):
        j = j0
        for k in range(1, n + 1):
            i = i0 + sgn * k
            lo, hi = max(j - 2, 0), min(j + 3, g.N)
            j = lo + int(np.argmax(a[lo:hi, i]))
            if a[j, i] < 0.5:
                return None
        out.append(g.xc[j] - g.xc[j0])
    return out


def run_g2():
    R = Report("G2")
    P, NP = PhysParams(), NumParams()
    imgs, labels = [], []
    # G2.1 control a = 0
    s0 = Sim(wobble_action(a=0.0), P=P, NP=NP)
    s0.run()
    a0 = alpha_of(s0, s0.Yf, s0.U[0])
    pk0, *_ = axis_stripes(a0, s0.gf)
    R.check("G2.1", "wobble amplitude 0: number of stripes on the axis", mt.n_stripes(pk0),
            mt.n_stripes(pk0) == 0, "== 0")
    imgs.append(img_of(s0, s0.Yf, s0.U[0]))
    labels.append("a=0")

    # G2.2 / G2.3 frequency scan at fixed amplitude and retreat speed
    res = {}
    for f in (1.0, 1.5, 2.0):
        s = Sim(wobble_action(a=9.0, f=f), P=P, NP=NP)
        s.run()
        a = alpha_of(s, s.Yf, s.U[0])
        pk, ys, lam, prof = axis_stripes(a, s.gf)
        res[f] = dict(n=mt.n_stripes(pk), lam_mm=lam * 1e3,
                      lam_mean_mm=float(np.mean(np.diff(ys)) * 1e3) if len(ys) > 1 else float("nan"),
                      ys_mm=(ys * 1e3).round(2).tolist(),
                      contrast=mt.contrast(prof, pk))
        imgs.append(img_of(s, s.Yf, s.U[0]))
        labels.append(f"f={f:g} Hz")
        if f == 1.5:
            s15 = s
    R.extra["freq_scan"] = res
    R.check("G2.2", "a=9 mm, f=1.5 Hz: number of stripes on the axis", res[1.5]["n"], res[1.5]["n"] >= 4, ">= 4")
    lams = [res[f]["lam_mm"] for f in (2.0, 1.5, 1.0)]
    ok = all(np.isfinite(lams)) and lams[0] < lams[1] < lams[2]
    R.check("G2.3", "median stripe spacing increases with v_rel/f (f = 2, 1.5, 1 Hz) [mm]", lams, ok, "increasing")

    # G2.5 curvature: follow each axis stripe's ridge to x = +-8 mm; forward-convex if both ends fall back
    a = alpha_of(s15, s15.Yf, s15.U[0])
    pk, ys, lam, prof = axis_stripes(a, s15.gf)
    drops = [ridge_drop(a, s15.gf, j) for j in pk]
    tracked = [d for d in drops if d is not None]
    frac = float(np.mean([max(d) < 0 for d in tracked])) if tracked else 0.0
    R.extra["G2.5_ridge_drops_mm"] = [None if d is None else [x * 1e3 for x in d] for d in drops]
    R.check("G2.5", "share of tracked stripes that bend backwards at x = +-8 mm (forward-convex arcs)",
            [frac, len(tracked)], frac >= 0.6 and len(tracked) >= 3, ">= 0.6 of >= 3 tracked")

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
