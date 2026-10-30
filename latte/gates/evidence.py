"""Evidence runs for the closure decisions (docs/MODEL.md section 3, docs/REPORT.md failure analysis).

    python -m gates.evidence [--N 96]

F1  spec main pressure closure (buoyant, coffee pressureless):
    (a) snowplow: max coffee areal mass (Y < 0.05) during a low dwell, vs the fallback closure
    (b) creep: boundary displacement 1 s -> 6 s after the stop, buoyant vs fallback
F2  spec retention calibration (all retained milk surfaces, phi_f = 1, v_e = 1.05): flooded cup
F3  v_e = 1.05 with phi_f calibrated: no impact hole, no stripes (needs a penetrating share at low pour)
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
from actions.scripts import ALL                        # noqa: E402
from gates.common import gate_dir, alpha_of, img_of, step_until  # noqa: E402
from gates.g1 import dwell                             # noqa: E402
from gates.g2 import wobble_action, axis_stripes       # noqa: E402
from gates import metrics as mt                        # noqa: E402
from render.render import board                        # noqa: E402


def snowplow(P, NP, T=3.0):
    s = Sim(dwell(9.0, T, T), P=P, NP=NP)
    mx = []
    while s.t < T - 1e-12:
        s.step(min(s.dt_cfl(), T - s.t))
        if s.nstep % 50 == 0:
            M = s.U[0]
            Y = s.U[3] / np.where(s.g.wet, M, 1)
            coffee = (Y < 0.05) & s.g.wet
            mx.append((s.t, float(M[coffee].max() / (P.rho_c * P.h0))))
    s.flush_material()
    return s, mx


def creep(P, NP):
    s = Sim(dwell(9.0, 3.0, 9.0), P=P, NP=NP)
    step_until(s, 4.0)
    w1 = mt.white(s)
    img1 = img_of(s, s.Yf, s.U[0])
    step_until(s, 9.0)
    w2 = mt.white(s)
    return s, mt.boundary_displacement(w1, w2, s.gf), mt.area(w2, s.gf) / max(mt.area(w1, s.gf), 1e-30), img1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=96)
    a = ap.parse_args()
    NP = NumParams(N=a.N, ratio=3)
    d = gate_dir("evidence")
    P0 = PhysParams()
    Pb = P0.with_(pressure="buoyant")
    res = {}
    # F1a snowplow
    s_b, mx_b = snowplow(Pb, NP)
    s_e, mx_e = snowplow(P0, NP)
    res["F1a_snowplow_max_coffee_M_over_rho_c_h0"] = dict(buoyant=mx_b[-1][1], elastic=mx_e[-1][1],
                                                            buoyant_series=mx_b[::4], elastic_series=mx_e[::4])
    # F1b creep
    cb = creep(Pb, NP)
    ce = creep(P0, NP)
    res["F1b_creep_boundary_disp_mm"] = dict(buoyant=cb[1] * 1e3, elastic=ce[1] * 1e3,
                                              area_ratio_buoyant=cb[2], area_ratio_elastic=ce[2])
    board([[img_of(s_b, s_b.Yf, s_b.U[0]), img_of(s_e, s_e.Yf, s_e.U[0]), cb[3], img_of(cb[0], cb[0].Yf, cb[0].U[0]),
            ce[3], img_of(ce[0], ce[0].Yf, ce[0].U[0])]], os.path.join(d, "F1_pressure_closure.png"),
          labels=["buoyant t=3", "fallback t=3", "buoyant stop+1s", "buoyant stop+6s", "fallback stop+1s",
                  "fallback stop+6s"], cell=200)
    # F2 flood with the spec retention calibration (all retained milk surfaces)
    rows = []
    s2 = Sim(ALL["heart"](), P=P0.with_(phi_f=1.0, v_e=1.05), NP=NP)
    step_until(s2, 5.2)
    res["F2_flood_white_area_fraction"] = mt.area(mt.white(s2), s2.gf) / (np.pi * P0.R_cup ** 2)
    s2b = Sim(ALL["heart"](), P=P0, NP=NP)
    step_until(s2b, 5.2)
    res["F2_calibrated_white_area_fraction"] = mt.area(mt.white(s2b), s2b.gf) / (np.pi * P0.R_cup ** 2)
    # F3 no penetration at low pour
    s3 = Sim(wobble_action(a=9.0, f=1.5), P=P0.with_(v_e=1.05, phi_f=0.10), NP=NP)
    s3.run()
    a3 = alpha_of(s3, s3.Yf, s3.U[0])
    pk3, *_ = axis_stripes(a3, s3.gf)
    s3b = Sim(wobble_action(a=9.0, f=1.5), P=P0, NP=NP)
    s3b.run()
    a3b = alpha_of(s3b, s3b.Yf, s3b.U[0])
    pk3b, *_ = axis_stripes(a3b, s3b.gf)
    res["F3_stripes"] = dict(v_e_105=len(pk3), calibrated=len(pk3b))
    board([[img_of(s2, s2.Yf, s2.U[0]), img_of(s2b, s2b.Yf, s2b.U[0]), img_of(s3, s3.Yf, s3.U[0]),
            img_of(s3b, s3b.Yf, s3b.U[0])]], os.path.join(d, "F2_F3_retention.png"),
          labels=["phi_f=1, v_e=1.05 (t=5.2)", "calibrated (t=5.2)", "wobble, v_e=1.05", "wobble, calibrated"],
          cell=220)
    json.dump(res, open(os.path.join(d, "evidence.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
