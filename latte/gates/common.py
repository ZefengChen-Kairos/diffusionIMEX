"""Shared gate utilities: output directory, runs, reports."""
import json
import os
import time
import numpy as np

from core.solver_np import Sim
from core.params import PhysParams, NumParams
from render.render import render, tag_image, crop_cup, save, board


def out_root():
    r = os.environ.get("LATTE_OUT")
    if r:
        return r
    if os.path.isdir("/root/autodl-tmp"):
        return "/root/autodl-tmp/latte_out"
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")


def gate_dir(name):
    d = os.path.join(out_root(), "gates", name)
    os.makedirs(d, exist_ok=True)
    return d


class Report:
    def __init__(self, gate):
        self.gate = gate
        self.checks = []
        self.extra = {}
        self.dir = gate_dir(gate)
        self.t0 = time.time()

    def check(self, cid, desc, value, ok, thr=""):
        v = value
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        if isinstance(v, np.ndarray):
            v = v.tolist()
        self.checks.append(dict(id=cid, desc=desc, value=v, threshold=thr, PASS=bool(ok)))
        print(f"  [{'PASS' if ok else 'FAIL'}] {cid} {desc}: {v} ({thr})", flush=True)
        return ok

    @property
    def passed(self):
        return all(c["PASS"] for c in self.checks)

    def write(self, P=None, NP=None):
        out = dict(gate=self.gate, passed=self.passed, checks=self.checks, extra=self.extra,
                   wall_s=round(time.time() - self.t0, 1))
        if P is not None:
            out["params"] = P.to_dict()
            out["params_digest"] = P.digest()
        if NP is not None:
            out["numerics"] = NP.to_dict()
        with open(os.path.join(self.dir, "metrics.json"), "w") as f:
            json.dump(out, f, indent=1, default=str)
        return out


def run(action, P=None, NP=None, T=None, record_every=0.0, flags=None, backend="numba", **kw):
    s = Sim(action, P=P or PhysParams(), NP=NP or NumParams(), record_every=record_every,
            flags=flags, backend=backend, **kw)
    s.run(T=T)
    return s


def save_ledger(sim, path):
    L = np.array(sim.log) if sim.log else np.zeros((0, 6))
    np.savetxt(path, L, delimiter=",", header="t,err_mass,err_milk,err_tag,err_mom,M_min", comments="")


def frame_images(sim, with_tags=True):
    imgs, tags, ts = [], [], []
    for (t, Yf, M) in sim.frames:
        phi, _ = sim.phi_fine(Yf, M)
        imgs.append(crop_cup(render(phi, sim.gf, sim.P.Phi_star), sim.gf))
        if with_tags:
            tags.append(crop_cup(tag_image(Yf, sim.gf), sim.gf))
        ts.append(t)
    return imgs, tags, ts


def save_run(sim, d, name, board_frames=True):
    os.makedirs(d, exist_ok=True)
    save_ledger(sim, os.path.join(d, f"{name}_ledger.csv"))
    phi, _ = sim.phi_fine()
    save(crop_cup(render(phi, sim.gf, sim.P.Phi_star), sim.gf), os.path.join(d, f"{name}_final.png"))
    if board_frames and sim.frames:
        imgs, tags, ts = frame_images(sim)
        board([imgs, tags], os.path.join(d, f"{name}_frames.png"), labels=[f"t={t:.2f}s" for t in ts],
              cell=170)
    if sim.action is not None and hasattr(sim.action, "to_json"):
        with open(os.path.join(d, f"{name}_action.json"), "w") as f:
            f.write(sim.action.to_json())


_SHELLS = {}


def shell(K=1, P=None, NP=None):
    """Grid-only Sim for analysing cached snapshots."""
    P = P or PhysParams()
    NP = NP or NumParams()
    k = (K, P.digest(), NP.N, NP.ratio)
    if k not in _SHELLS:
        _SHELLS[k] = Sim(None, P=P, NP=NP, n_tags=K)
    return _SHELLS[k]


def alpha_of(sh, Y, M):
    from render.render import alpha
    phi, _ = sh.phi_fine(Y, M)
    return alpha(phi, sh.P.Phi_star) * sh.gf.wet


def img_of(sh, Y, M):
    phi, _ = sh.phi_fine(Y, M)
    return crop_cup(render(phi, sh.gf, sh.P.Phi_star), sh.gf)


def tagimg_of(sh, Y):
    return crop_cup(tag_image(Y, sh.gf), sh.gf)


def step_until(sim, T):
    while sim.t < T - 1e-12:
        sim.step(min(sim.dt_cfl(), T - sim.t))
    sim.flush_material()
    return sim
