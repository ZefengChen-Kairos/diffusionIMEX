"""G1: fixed-point injection -> big white heart."""
import numpy as np

from actions.scripts import Action, Segment, heart
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from gates.common import Report, shell, alpha_of, img_of, step_until, save_run
from gates import metrics as mt
from render.render import board, save
from run_patterns import run_pattern


def dwell(Q, T_inj=3.0, T_end=None, x=(0, -10), h=0.012):
    return Action(f"dwell_Q{Q:g}", [Segment(0.0, T_inj, x, None, Q, None, h)], T_end or T_inj)


def heart_metrics(a, g):
    """Lobes / notch / tip descriptors of a white mask after the cut (axis +y)."""
    w = a > 0.5
    ys = g.Y[w]
    if ys.size == 0:
        return dict(lobes=0, notch=False, tip_len=0.0, ncomp=0)
    y0, y1 = ys.min(), ys.max()
    # widths per row
    rows = np.arange(g.N)
    width = w.sum(axis=1) * g.dx
    yc = g.xc
    wmax = width.max()
    j_wmax = int(np.argmax(width))
    # tip: forward length where the width is below 30% of the max width
    fwd = (yc > yc[j_wmax]) & (width > 0)
    narrow = fwd & (width < 0.3 * wmax)
    tip_len = float(narrow.sum() * g.dx)
    # lobes / notch: lateral profile in the rear third
    jr = (yc > y0) & (yc < y0 + 0.35 * (y1 - y0))
    prof = a[jr].mean(axis=0)
    mid = np.abs(g.xc) < 1.5e-3
    side = (np.abs(g.xc) > 3e-3) & (np.abs(g.xc) < 15e-3)
    notch = bool(prof[mid].mean() < 0.8 * np.max(prof[side]) and prof[mid].mean() < 0.6)
    pk = mt.stripes(prof, amp=0.15)
    _, ncomp = mt.components(w, 40)
    return dict(lobes=len(pk), notch=notch, tip_len=tip_len, ncomp=ncomp, width_max=float(wmax),
                y_extent=[float(y0), float(y1)])


def run_g1():
    R = Report("G1")
    P, NP = PhysParams(), NumParams()
    sh = None

    # ---- the heart itself (cached pattern run)
    r = run_pattern("heart", P, NP)
    K = r["Y"].shape[1] - 1
    sh = shell(K, P, NP)
    g = sh.gf
    t = r["t"]
    areas = np.array([mt.area(alpha_of(sh, r["Y"][i].astype(float), r["M"][i].astype(float)) > 0.5, g)
                      for i in range(len(t))])
    inj = (t >= 0.8) & (t <= 5.2)
    dA = np.diff(areas[inj])
    tol = 2 * g.dA
    mono = bool((dA >= -tol).all())
    R.check("G1.1", "A(t) non-decreasing during the low dwell [0.8, 5.2] s (min dA, mm^2)",
            float(dA.min() * 1e6), mono, f">= -{tol * 1e6:.2f}")
    R.extra["heart_area_mm2"] = dict(t=t.tolist(), A=(areas * 1e6).round(1).tolist())
    rate = float(np.polyfit(t[(t >= 2.0) & (t <= 5.0)], areas[(t >= 2.0) & (t <= 5.0)], 1)[0] * 1e6)
    R.extra["heart_dA_dt_mm2_per_s"] = rate
    R.extra["video_dA_dt_mm2_per_s"] = 350.0

    # ---- G1.2 dA/dt vs Q (dwell runs)
    rates = {}
    snaps = {}
    for Q in (6.0, 9.0, 12.0):
        s = Sim(dwell(Q), P=P, NP=NP, record_every=0.25)
        s.run()
        A = np.array([mt.area(alpha_of(s, f[1], f[2]) > 0.5, s.gf) for f in s.frames])
        tt = np.array([f[0] for f in s.frames])
        sel = tt >= 1.0
        rates[Q] = float(np.polyfit(tt[sel], A[sel], 1)[0] * 1e6)
        snaps[Q] = s
    vals = [rates[q] for q in (6.0, 9.0, 12.0)]
    R.check("G1.2", "mean dA/dt [mm^2/s] for Q = 6, 9, 12 ml/s strictly increasing", vals,
            vals[0] < vals[1] < vals[2], "increasing")
    from closures.jet_stub import jet_speed, retention, jet_diameter
    R.extra["G1.2_r_of_Q"] = {q: retention(jet_speed(q * 1e-6, 0.012, P), 0, P) for q in (6.0, 9.0, 12.0)}

    # ---- G1.3 / G1.4 hole and horseshoe (Q = 9 dwell, mid-injection snapshot)
    s9 = snaps[9.0]
    fr = min(s9.frames, key=lambda f: abs(f[0] - 2.0))
    a = alpha_of(s9, fr[1], fr[2])
    xp = np.array([0.0, -0.010])
    i0 = int(np.argmin(np.abs(s9.gf.xc - xp[0])))
    j0 = int(np.argmin(np.abs(s9.gf.xc - xp[1])))
    a_xp = float(a[j0, i0])
    row = a[j0]
    # lateral brown gap containing x_p
    left = i0
    while left > 0 and row[left] < 0.5:
        left -= 1
    right = i0
    while right < len(row) - 1 and row[right] < 0.5:
        right += 1
    hole_w = float((right - left - 1) * s9.gf.dx) if a_xp < 0.5 else 0.0
    R.check("G1.3", "impact hole: alpha(x_p) < 0.5 and lateral brown width 3-10 mm (alpha, mm)",
            [a_xp, hole_w * 1e3], a_xp < 0.5 and 3e-3 <= hole_w <= 10e-3, "alpha<0.5, 3..10 mm")
    c = mt.centroid(a > 0.5, s9.gf)
    R.check("G1.4", "horseshoe: white centroid ahead of x_p (dy mm)", float((c[1] - xp[1]) * 1e3),
            c[1] - xp[1] > 0, "> 0")
    save(img_of(s9, fr[1], fr[2]), R.dir + "/dwell_Q9_t2.png")

    # ---- G1.5 / G1.6 freeze and creep: dwell 3 s then 6 s free
    s = Sim(dwell(9.0, 3.0, 9.0), P=P, NP=NP)
    step_until(s, 3.0)
    ts, us = [], []
    while s.t < 4.5 - 1e-12:
        s.step(min(s.dt_cfl(), 4.5 - s.t))
        M, u, v, Y, *_ = s.prims(s.U)
        wmask = (Y > 0.2) & s.g.wet
        ts.append(s.t - 3.0)
        us.append(float(np.sqrt(u[wmask] ** 2 + v[wmask] ** 2).max()) if wmask.any() else 0.0)
    ts, us = np.array(ts), np.array(us)
    sel = (us > 1e-9) & (ts > 0.02) & (ts < 1.0)
    tau = float(-1.0 / np.polyfit(ts[sel], np.log(us[sel]), 1)[0])
    R.check("G1.5", "post-stop velocity e-folding time of the white region [s]", tau, 0.1 <= tau <= 0.5,
            "0.1..0.5")
    R.extra["G1.5_u_decay"] = dict(t=ts[::20].tolist(), umax=us[::20].tolist())
    step_until(s, 4.0)
    w1 = mt.white(s)
    step_until(s, 9.0)
    w2 = mt.white(s)
    dmax = mt.boundary_displacement(w1, w2, s.gf)
    dmean = mt.mean_boundary_displacement(w1, w2, s.gf)
    R.check("G1.6", "creep: boundary displacement 1 s -> 6 s after stop, max / mean [mm]",
            [dmax * 1e3, dmean * 1e3], dmax < 2e-3, "max < 2 mm")

    # ---- G1.7 heart topology after the cut; G1.8 low-cut control
    Yf, Mf = r["Yfinal"].astype(float), r["Mfinal"].astype(float)
    a = alpha_of(sh, Yf, Mf)
    hm = heart_metrics(a, sh.gf)
    R.extra["G1.7_metrics"] = hm
    ok7 = hm["notch"] and hm["lobes"] >= 2 and hm["tip_len"] >= 5e-3 and hm["ncomp"] <= 2
    R.check("G1.7", "after the cut: rear notch, two lobes, tip >= 5 mm, connected", hm, ok7,
            "notch & lobes>=2 & tip>=5mm & ncomp<=2")
    act = heart()
    for sg in act.segs:
        if sg.t0 >= 5.2:
            sg.h0 = 0.012
            sg.h1 = 0.012
    c8 = Sim(act, P=P, NP=NP)
    c8.run()
    a8 = alpha_of(c8, c8.Yf, c8.U[0])
    hm8 = heart_metrics(a8, c8.gf)
    R.extra["G1.8_metrics"] = hm8
    R.check("G1.8", "control (cut at h=1.2 cm): no rear notch / no narrow tip", hm8,
            (not hm8["notch"]) or hm8["tip_len"] < 0.5 * hm["tip_len"], "no notch or tip < half")
    board([[img_of(sh, Yf, Mf), img_of(c8, c8.Yf, c8.U[0])]], R.dir + "/heart_vs_lowcut.png",
          labels=["heart (cut h=10 cm)", "control (cut h=1.2 cm)"], cell=260)
    R.write(P, NP)
    return R
