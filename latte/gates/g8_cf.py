"""G8 generalisation (action perturbations, dt halving, cost) and the counterfactual experiments."""
import copy
import json
import os
import time
import resource
import numpy as np

from actions.scripts import ALL, Sample
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from closures.jet_stub import jet_speed
from gates.common import Report, alpha_of, img_of, step_until, out_root
from gates import metrics as mt
from gates.g1 import heart_metrics
from gates.g4_7 import leaf_metrics
from render.render import board
from run_patterns import run_pattern


def heart():
    return ALL["heart"]()


def leaf():
    return ALL["leaf"]()


def perturb(act, rng, eps=0.10):
    """+-eps random scaling of the pour's timing, flow and positions (the hand motion stays smooth)."""
    return act.perturbed(rng, eps)


def classify(name, sim):
    a = alpha_of(sim, sim.Yf, sim.U[0])
    if name == "heart":
        m = heart_metrics(a, sim.gf)
        return m["notch"] and m["lobes"] >= 2 and m["tip_len"] >= 5e-3, m
    m = leaf_metrics(sim, sim.Yf, sim.U[0])
    return m["ribs"] >= 5 and m["sym"] >= 0.6, m


def run_g8():
    R = Report("G8")
    P, NP = PhysParams(), NumParams()
    rng = np.random.default_rng(1234)
    imgs, labels = [], []
    for name in ("heart", "leaf"):
        ref = run_pattern(name, P, NP)
        K = ref["Y"].shape[1] - 1
        sref = Sim(None, P=P, NP=NP, n_tags=K)
        wref = alpha_of(sref, ref["Yfinal"].astype(float), ref["Mfinal"].astype(float)) > 0.5
        for k in range(2):
            act = perturb(ALL[name](), rng)
            s = Sim(act, P=P, NP=NP)
            s.run()
            ok_cls, m = classify(name, s)
            io = mt.iou(wref, mt.white(s))
            R.check(f"G8.1-{name}-{k}", f"{name} +-10 % perturbation: class kept, white IoU vs nominal",
                    [ok_cls, io], ok_cls and io >= 0.7, "class & IoU >= 0.7")
            imgs.append(img_of(s, s.Yf, s.U[0]))
            labels.append(f"{name} pert {k}")
    # dt halving
    ref = run_pattern("heart", P, NP)
    t0 = time.time()
    s = Sim(heart(), P=P, NP=NP.with_(cfl=0.175))
    s.run()
    wall_half = time.time() - t0
    K = ref["Y"].shape[1] - 1
    sref = Sim(None, P=P, NP=NP, n_tags=K)
    wref = alpha_of(sref, ref["Yfinal"].astype(float), ref["Mfinal"].astype(float)) > 0.5
    io = mt.iou(wref, mt.white(s))
    R.check("G8.2", "dt halved (CFL 0.35 -> 0.175): heart white IoU", io, io >= 0.9, ">= 0.9")
    imgs.append(img_of(s, s.Yf, s.U[0]))
    labels.append("heart, dt/2")
    # cost report
    tj = os.path.join(out_root(), "patterns", "timing.json")
    if os.path.exists(tj):
        R.extra["timing"] = json.load(open(tj))["runs"]
    R.extra["heart_dt_half_wall_s"] = wall_half
    R.extra["max_rss_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    R.check("G8.3", "cost reported (s per physical s, peak RSS MB)", R.extra["max_rss_MB"], True, "report")
    board([imgs], R.dir + "/perturbations.png", labels=labels, cell=200)
    R.write(P, NP)
    return R


# ---------------------------------------------------------------------------- counterfactuals
class HeightRelabel:
    """Same jet (Q, v_j) but the action carries a different height label."""

    def __init__(self, base, factor):
        self.base, self.factor = base, factor
        self.name, self.T_end, self.n_tags = base.name + "_relabel", base.T_end, base.n_tags

    def sample(self, t, dt_fd=None):
        s = self.base.sample(t)
        return Sample(s.x, s.v, s.Q, s.h * self.factor, s.e, s.tag)

    def to_json(self):
        return self.base.to_json()


def run_cf():
    R = Report("CF")
    P, NP = PhysParams(), NumParams()
    # CF1: colour only -> no push.  Paint a tag disc mid-run; the dynamics must be bitwise unchanged.
    a = Sim(leaf(), P=P, NP=NP)
    b = Sim(leaf(), P=P, NP=NP)
    step_until(a, 2.0)
    step_until(b, 2.0)
    g = b.g
    disc = (g.X ** 2 + (g.Y - 0.01) ** 2) < 0.008 ** 2
    k = b.K - 1
    b.U[4 + k][disc] = np.minimum(b.U[3][disc], b.U[4 + k][disc] + 0.5 * b.U[3][disc])
    b.U[4 + k] = np.minimum(b.U[4 + k], b.U[3])
    # fine colour as well
    gf = b.gf
    discf = (gf.X ** 2 + (gf.Y - 0.01) ** 2) < 0.008 ** 2
    b.Yf[1 + k][discf] = b.Yf[0][discf]
    # (the paint is a pure label: its own tag ledger is not checked after painting)
    b.NP = b.NP.with_(check_every=0)
    step_until(a, 4.0)
    step_until(b, 4.0)
    same = bool(np.array_equal(a.U[:4], b.U[:4]) and np.array_equal(a.Yf[0], b.Yf[0]))
    R.check("CF1a", "colour-only paint (tag field) mid-run: M, P, Phi and Y bitwise unchanged", same, same, "bitwise")
    # colour-only injection: D with no mass and no momentum produces no flow from a static start
    c = Sim(heart(), P=P, NP=NP, flags=dict(no_D_mass=True, no_D_momentum=True, no_E=True, no_tau=True))
    step_until(c, 3.0)
    umax = float(np.abs(c.U[1:3]).max())
    R.check("CF1b", "injection with no volume / momentum input (D, E, tau off): max |P| after 3 s", umax,
            umax == 0.0, "== 0")
    # CF2: same jet data, height label changed -> identical output
    base = heart()
    vj = lambda t: jet_speed(base.sample(t).Q, base.sample(t).h, P) if base.sample(t).Q > 0 else 1.0
    s1 = Sim(base, P=P, NP=NP, v_j_override=vj)
    s2 = Sim(HeightRelabel(heart(), 3.0), P=P, NP=NP, v_j_override=vj)
    step_until(s1, 6.5)
    step_until(s2, 6.5)
    same2 = bool(np.array_equal(s1.U, s2.U) and np.array_equal(s1.Yf, s2.Yf))
    R.check("CF2", "height label x3 with the same (Q, v_j): output bitwise identical", same2, same2, "bitwise")
    # CF3: no screening (R = 0) -> uniform shrinkage
    ref = run_pattern("heart", P, NP)
    K = ref["Y"].shape[1] - 1
    sh = Sim(None, P=P, NP=NP, n_tags=K)
    w_ref = alpha_of(sh, ref["Yfinal"].astype(float), ref["Mfinal"].astype(float)) > 0.5
    n = Sim(heart(), P=P, NP=NP, flags=dict(no_R=True))
    n.run()
    w_n = mt.white(n)
    far = np.abs(sh.gf.X) > 10e-3
    A_ref, A_n = mt.area(w_ref & far, sh.gf), mt.area(w_n & far, sh.gf)
    Atot_ref, Atot_n = mt.area(w_ref, sh.gf), mt.area(w_n, sh.gf)
    shrink_far = 1 - A_n / A_ref
    R.extra["CF3"] = dict(A_far_ref_mm2=A_ref * 1e6, A_far_noR_mm2=A_n * 1e6, A_ref_mm2=Atot_ref * 1e6,
                          A_noR_mm2=Atot_n * 1e6, M_mean_ref=float(ref["Mfinal"][sh.g.wet].mean()),
                          M_mean_noR=float(n.U[0][n.g.wet].mean()))
    R.check("CF3", "R = 0 (bare sink): white area far from the cut line (|x|>10 mm) shrinks [fraction]",
            shrink_far, shrink_far > 0.1, "> 10 %")
    board([[img_of(sh, ref["Yfinal"].astype(float), ref["Mfinal"].astype(float)), img_of(n, n.Yf, n.U[0])]],
          R.dir + "/cf3_screening.png", labels=["screened sink (E+R)", "bare sink (R=0)"], cell=260)
    R.write(P, NP)
    return R
