"""G9: similarity to the reference videos (task priority 3).

Per pattern: white-mask IoU between the simulated final pattern and the rectified final frame of the
reference video (gates/reference.py), maximised over rotation and mirror.  Reported per pattern with
an aligned overlay (white = both, blue = simulation only, orange = video only).  The reference
masks are rebuilt from the videos when available.
"""
import os
import numpy as np

from core.params import PhysParams, NumParams
from gates.common import Report, shell, out_root
from gates import reference as ref
from run_patterns import run_pattern

NAMES = ["heart", "push_heart", "layered_heart", "tulip", "leaf", "swan"]


def run_g9():
    R = Report("G9")
    P, NP = PhysParams(), NumParams()
    rdir = os.path.join(out_root(), "reference")
    if not os.path.exists(os.path.join(rdir, "heart_mask.npy")):
        ref.build(rdir)
    runs = {}
    for n in NAMES:
        r = run_pattern(n, P, NP)
        runs[n] = (r["Yfinal"].astype(float), r["Mfinal"].astype(float), shell(r["Y"].shape[1] - 1, P, NP))
    res = ref.score_runs(runs, R.dir, P, NP, tag="similarity")
    for n in NAMES:
        if n in res:
            v = res[n]
            R.check(f"G9-{n}", f"{n}: white IoU vs the rectified video frame (rotation {v['angle']} deg, "
                    f"mirror {v['mirror']}); white fraction sim / video",
                    [v["iou"], v["white_sim"], v["white_ref"]], True, "report")
    R.extra["mean_iou"] = float(np.mean([v["iou"] for v in res.values()]))
    R.extra["per_pattern"] = res
    R.write(P, NP)
    return R
