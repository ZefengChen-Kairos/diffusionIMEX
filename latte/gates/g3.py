"""G3: low / high pour switch (draw vs cut), and separation of C_t (drag) and C_E (collapse)."""
import numpy as np

from actions.scripts import Action, Segment
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from gates.common import Report, alpha_of, img_of, step_until
from render.render import board

T_DISC = 2.5
T_PASS0, T_PASS1 = 3.0, 4.2


def pass_action(h, Q=5.0, x0=(0, -30), x1=(0, 30)):
    S = Segment
    return Action(f"pass_h{h:g}", [
        S(0.0, T_DISC, (0, -6), None, 9.0, None, 0.012),               # white disc, tag A
        S(T_PASS0, T_PASS1, x0, x1, Q, None, h),                        # straight pass, tag B
    ], 5.5)


def tag_centroid(sim, k):
    g = sim.gf
    Mf = sim.fm.to_fine(sim.g.extrapolate(sim.U[0].copy()))
    w = sim.Yf[1 + k] * Mf * g.wet
    s = w.sum()
    return np.array([(w * g.X).sum() / s, (w * g.Y).sum() / s]), s * g.dA


def pass_run(P, NP, h, flags=None):
    s = Sim(pass_action(h), P=P, NP=NP, flags=flags)
    step_until(s, T_PASS0)
    cA0, mA0 = tag_centroid(s, 0)
    phi0 = s.totals()["Phi"]
    ret0 = s.ret_mass
    E0, R0 = s.E_milk, s.R_milk
    # time-integrated surface fluxes within 10 mm of the pass line during the pass
    g = s.g
    near = (np.abs(g.X) < 10e-3) & g.wet
    inflow = 0.0      # inward (towards the line) lateral mass flux
    drag = 0.0        # along-path mass flux
    while s.t < 5.5 - 1e-12:
        dt = min(s.dt_cfl(), 5.5 - s.t)
        s.step(dt)
        if s.t <= T_PASS1 + 0.3:
            inflow += dt * float((-np.sign(g.X) * s.U[1])[near].sum() * g.dA)
            drag += dt * float(s.U[2][near].sum() * g.dA)
    s.flush_material()
    cA1, mA1 = tag_centroid(s, 0)
    phi1 = s.totals()["Phi"]
    return dict(sim=s, dcA=(cA1 - cA0), retained=s.ret_mass - ret0, dphi=phi1 - phi0,
                E=s.E_milk - E0, R=s.R_milk - R0, inflow=inflow, drag=drag)


def run_g3():
    R = Report("G3")
    P, NP = PhysParams(), NumParams()
    lo = pass_run(P, NP, 0.015)
    hi = pass_run(P, NP, 0.10)
    ratio = lo["retained"] / max(hi["retained"], 1e-30)
    R.check("G3.1", "retained milk mass ratio low (h=1.5 cm) / high (h=10 cm), same path", ratio, ratio >= 10, ">= 10")
    dA = float(np.linalg.norm(hi["dcA"]) * 1e3)
    R.check("G3.2a", "high pass drags the old white: |d centroid(A)| [mm]", dA, dA > 0.5, "> 0.5 mm")
    R.check("G3.2b", "high pass adds ~no white: retained mass relative to the low pass",
            hi["retained"] / lo["retained"], hi["retained"] / lo["retained"] < 0.1, "< 0.1")
    s = hi["sim"]
    closure = abs(hi["dphi"] - (hi["retained"] - hi["E"] + hi["R"])) / max(abs(hi["E"]), 1e-30)
    R.check("G3.3", "milk removed by E is exactly the Phi ledger change: |dPhi - (D - E + R)| / E",
            closure, closure < 1e-9, "< 1e-9")
    R.extra["low"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in lo.items() if k != "sim"}
    R.extra["high"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in hi.items() if k != "sim"}
    # G3.4 separation of C_t and C_E with identical (1 - r) mdot
    nt = pass_run(P.with_(C_t=0.0), NP, 0.10)
    ne = pass_run(P.with_(C_E=0.0), NP, 0.10)
    R.extra["C_t=0"] = dict(drag=nt["drag"], inflow=nt["inflow"], dcA=nt["dcA"].tolist())
    R.extra["C_E=0"] = dict(drag=ne["drag"], inflow=ne["inflow"], dcA=ne["dcA"].tolist())
    R.check("G3.4a", "C_t = 0: along-path drag flux < 50 % of the reference", [nt["drag"], hi["drag"]],
            abs(nt["drag"]) < 0.5 * abs(hi["drag"]), "< 50 %")
    R.check("G3.4b", "C_E = 0: inward (collapse) flux < 50 % of the reference", [ne["inflow"], hi["inflow"]],
            abs(ne["inflow"]) < 0.5 * abs(hi["inflow"]), "< 50 %")
    board([[img_of(x["sim"], x["sim"].Yf, x["sim"].U[0]) for x in (lo, hi, nt, ne)]], R.dir + "/passes.png",
          labels=["low h=1.5cm", "high h=10cm", "high, C_t=0", "high, C_E=0"], cell=230)
    R.write(P, NP)
    return R
