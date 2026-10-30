"""Parameter sensitivity: each calibrated parameter at -30 % / +30 % -> G1 and G2 metrics.

    python -m gates.sensitivity [--N 96]

G1 metrics come from a 3 s low dwell (Q = 9 ml/s): white-area growth rate, impact-hole width,
forward offset of the white centroid; plus the post-stop velocity e-folding time.
G2 metrics come from the layered-heart drawing phase (a = 9 mm, f = 1.5 Hz): stripe count,
mean spacing, contrast.  Output: <out>/gates/sensitivity/sensitivity.md / .json
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.params import PhysParams, NumParams          # noqa: E402
from core.solver_np import Sim                         # noqa: E402
from gates.common import gate_dir, alpha_of, step_until  # noqa: E402
from gates.g1 import dwell                             # noqa: E402
from gates.g2 import wobble_action, axis_stripes       # noqa: E402
from gates import metrics as mt                        # noqa: E402

PARAMS = ["gamma", "chi", "k_h", "nu_s", "phi_f", "v_e", "dv", "C_E", "C_t", "v_c", "sig_k_min", "Phi_star"]


def g1_metrics(P, NP):
    s = Sim(dwell(9.0, 3.0, 3.6), P=P, NP=NP, record_every=0.25)
    step_until(s, 3.0)
    A = np.array([mt.area(alpha_of(s, f[1], f[2]) > 0.5, s.gf) for f in s.frames])
    tt = np.array([f[0] for f in s.frames])
    sel = tt >= 1.0
    rate = float(np.polyfit(tt[sel], A[sel], 1)[0] * 1e6)
    fr = min(s.frames, key=lambda f: abs(f[0] - 2.0))
    a = alpha_of(s, fr[1], fr[2])
    i0 = int(np.argmin(np.abs(s.gf.xc)))
    j0 = int(np.argmin(np.abs(s.gf.xc + 0.010)))
    row = a[j0]
    hole = 0.0
    if row[i0] < 0.5:
        l = i0
        while l > 0 and row[l] < 0.5:
            l -= 1
        r = i0
        while r < len(row) - 1 and row[r] < 0.5:
            r += 1
        hole = (r - l - 1) * s.gf.dx * 1e3
    c = mt.centroid(a > 0.5, s.gf)
    ts, us = [], []
    while s.t < 3.6 - 1e-12:
        s.step(min(s.dt_cfl(), 3.6 - s.t))
        M, u, v, Y, *_ = s.prims(s.U)
        w = (Y > 0.2) & s.g.wet
        ts.append(s.t - 3.0)
        us.append(float(np.sqrt(u[w] ** 2 + v[w] ** 2).max()) if w.any() else 0.0)
    ts, us = np.array(ts), np.array(us)
    sel = (us > 1e-9) & (ts > 0.02)
    tau = float(-1.0 / np.polyfit(ts[sel], np.log(us[sel]), 1)[0]) if sel.sum() > 3 else float("nan")
    return dict(dAdt=rate, hole_mm=hole, offset_mm=float((c[1] + 0.010) * 1e3), tau_s=tau)


def g2_metrics(P, NP):
    s = Sim(wobble_action(a=9.0, f=1.5), P=P, NP=NP)
    s.run()
    a = alpha_of(s, s.Yf, s.U[0])
    pk, ys, lam, prof = axis_stripes(a, s.gf)
    return dict(n=mt.n_stripes(pk), lam_mm=lam * 1e3, contrast=mt.contrast(prof, pk))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=96)
    ap.add_argument("--ratio", type=int, default=3)
    ap.add_argument("--params", nargs="*", default=PARAMS)
    a = ap.parse_args()
    NP = NumParams(N=a.N, ratio=a.ratio)
    P0 = PhysParams()
    d = gate_dir("sensitivity")
    res = {"base": dict(g1=g1_metrics(P0, NP), g2=g2_metrics(P0, NP))}
    print("base", res["base"], flush=True)
    for p in a.params:
        for f in (0.7, 1.3):
            P = P0.with_(**{p: getattr(P0, p) * f})
            res[f"{p}x{f}"] = dict(g1=g1_metrics(P, NP), g2=g2_metrics(P, NP))
            print(p, f, res[f"{p}x{f}"], flush=True)
        json.dump(res, open(os.path.join(d, "sensitivity.json"), "w"), indent=1)
    b = res["base"]
    keys = [("g1", "dAdt"), ("g1", "hole_mm"), ("g1", "offset_mm"), ("g1", "tau_s"),
            ("g2", "n"), ("g2", "lam_mm"), ("g2", "contrast")]
    lines = [f"# Sensitivity (+-30 %, N={a.N}, fine {a.N * a.ratio})", "",
             "Base: " + ", ".join(f"{k}={b[g][k]:.3g}" for g, k in keys), "",
             "Entries: value at -30 % / +30 % (relative change of the continuous metrics in brackets).", "",
             "| parameter | " + " | ".join(k for _, k in keys) + " |", "|---|" + "---|" * len(keys)]
    for p in a.params:
        cells = []
        for g, k in keys:
            lo, hi = res[f"{p}x0.7"][g][k], res[f"{p}x1.3"][g][k]
            base = b[g][k]
            if k == "n":
                cells.append(f"{lo} / {hi}")
            else:
                rel = lambda v: (v - base) / base * 100 if base else float("nan")
                cells.append(f"{lo:.3g} / {hi:.3g} ({rel(lo):+.0f}%/{rel(hi):+.0f}%)")
        lines.append(f"| {p} | " + " | ".join(cells) + " |")
    open(os.path.join(d, "sensitivity.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
