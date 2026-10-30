"""Run the six patterns with the single shared PhysParams; cache snapshots; build boards.

    python run_patterns.py [names...] [--N 128] [--force]

Outputs (under $LATTE_OUT or /root/autodl-tmp/latte_out or ./out):
  patterns/<name>.npz         snapshots every 0.25 s (float16 Y, tags; float32 M), ledger, diag
  patterns/<name>_final.png   final render
  patterns/<name>_stages.png  4 stage frames + tag map
  patterns/board.png          six patterns x (4 stages, final, reference)
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import json
import time
import hashlib
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.params import PhysParams, NumParams                      # noqa: E402
from core.solver_np import Sim                                       # noqa: E402
from actions.scripts import ALL                                      # noqa: E402
from render.render import render, tag_image, crop_cup, save, board   # noqa: E402
from gates.common import out_root                                    # noqa: E402

NAMES = ["heart", "push_heart", "layered_heart", "tulip", "leaf", "swan"]
LABEL_ZH = {"heart": "大白心", "push_heart": "推推乐", "layered_heart": "千层心", "tulip": "郁金香",
            "leaf": "树叶", "swan": "天鹅"}
SNAP_DT = 0.25


def stage_times(action, t_end):
    """Three stages before the end (the ends of the drawing parts, or fractions of a single drawing
    part, up to the start of the pull-through), and the end."""
    if action.name == "swan":            # wing, neck, head
        return [round(action.end(l) + 0.1, 2) for l in ("wing", "neck", "head")] + [t_end]
    t_cut = action.start("cut")
    parts = [l for l in action.labels() if l != "cut" and not l.startswith("lift")]
    if len(parts) >= 3:                  # push heart, tulip: after A, after B, before the pull-through
        return [round(action.end(l), 2) for l in parts[:2]] + [round(t_cut - 0.05, 2), t_end]
    return [round(0.4 * t_cut, 2), round(0.7 * t_cut, 2), round(t_cut - 0.05, 2), t_end]


def code_hash():
    """Hash of the solver / closure / grid sources: cached runs are invalidated by code changes."""
    root = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha1()
    for f in ("core/solver_np.py", "core/grid.py", "core/kernels_nb.py", "closures/jet_stub.py",
              "actions/scripts.py"):
        if not os.path.exists(os.path.join(root, f)):
            continue
        h.update(open(os.path.join(root, f), "rb").read())
    return h.hexdigest()[:12]


def key(name, P, NP, action):
    h = hashlib.sha1((P.digest() + json.dumps(NP.to_dict(), sort_keys=True) + action.to_json()
                      + code_hash()).encode())
    return h.hexdigest()[:12]


def pattern_dir():
    d = os.path.join(out_root(), "patterns")
    os.makedirs(d, exist_ok=True)
    return d


def run_pattern(name, P=None, NP=None, action=None, force=False, tag=None, verbose=True):
    """Run (or load from cache) one pattern; returns dict with snapshots and the Sim (if run)."""
    P = P or PhysParams()
    NP = NP or NumParams()
    action = action or ALL[name]()
    k = key(name, P, NP, action)
    tag = tag or name
    path = os.path.join(pattern_dir(), f"{tag}_{k}.npz")
    if os.path.exists(path) and not force:
        d = dict(np.load(path, allow_pickle=True))
        d["path"] = path
        return d
    sim = Sim(action, P=P, NP=NP, record_every=SNAP_DT)
    t0 = time.time()
    sim.run()
    wall = time.time() - t0
    ts = np.array([f[0] for f in sim.frames])
    Y = np.array([f[1] for f in sim.frames]).astype(np.float16)
    M = np.array([f[2] for f in sim.frames]).astype(np.float32)
    L = np.array(sim.log)
    out = dict(t=ts, Y=Y, M=M, ledger=L.astype(np.float32), diag=np.array(sim.diag, np.float32),
               Yfinal=sim.Yf.astype(np.float32), Mfinal=sim.U[0].astype(np.float32),
               Ufinal=sim.U.astype(np.float64), wall=wall, steps=sim.nstep, T=sim.t,
               params=json.dumps(P.to_dict()), numerics=json.dumps(NP.to_dict()),
               action=action.to_json(), bulk=np.array([sim.bulk.Mm, sim.bulk.Mc]),
               ret_mass=sim.ret_mass, inj_mass=sim.inj["mass"],
               E_milk=sim.E_milk, R_milk=sim.R_milk, E_tags=sim.E_tags, R_tags=sim.R_tags,
               inj_tags=sim.inj["tags"], bulk_tags=sim.bulk.Mk)
    np.savez_compressed(path, **out)
    if verbose:
        print(f"{tag}: {sim.nstep} steps, {sim.t:.2f} s physical, wall {wall:.1f} s "
              f"({wall / sim.t:.1f} s per physical s), max ledger err {L[:, 1:5].max():.2e}", flush=True)
    out["path"] = path
    out["sim"] = sim
    return out


def make_sim_shell(P=None, NP=None, K=1):
    """A Sim used only for grids / rendering of cached snapshots."""
    return Sim(None, P=P or PhysParams(), NP=NP or NumParams(), n_tags=K)


def snapshot(d, t):
    i = int(np.argmin(np.abs(d["t"] - t)))
    return d["Y"][i].astype(np.float64), d["M"][i].astype(np.float64), float(d["t"][i])


def render_snap(shell, Y, M, P):
    phi, _ = shell.phi_fine(Y, M)
    return crop_cup(render(phi, shell.gf, P.Phi_star), shell.gf)


def make_outputs(names, runs, P, NP):
    from PIL import Image
    d = pattern_dir()
    rows, labels = [], []
    ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
    for n in names:
        r = runs[n]
        K = r["Y"].shape[1] - 1
        shell = make_sim_shell(P, NP, K)
        st = stage_times(ALL[n](), float(r["t"][-1]))
        imgs, tags = [], []
        for t in st:
            Y, M, tt = snapshot(r, t)
            imgs.append(render_snap(shell, Y, M, P))
            tags.append(crop_cup(tag_image(Y, shell.gf), shell.gf))
        board([imgs, tags], os.path.join(d, f"{n}_stages.png"), labels=[f"t={t:.1f}s" for t in st], cell=220)
        save(imgs[-1], os.path.join(d, f"{n}_final.png"))
        ref = os.path.join(ref_dir, f"final_{n}.png")
        refim = None
        if os.path.exists(ref):
            im = Image.open(ref).convert("RGB")
            w, h = im.size
            refim = np.array(im.crop((0, (h - w) // 2, w, (h - w) // 2 + w)))
        rows.append(imgs + [refim])
        labels.append(f"{n}")
    board(rows, os.path.join(d, "board.png"), labels=["stage 1", "stage 2", "stage 3", "final", "reference video"],
          row_labels=labels, cell=220)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", default=NAMES)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--ratio", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    P = PhysParams()
    NP = NumParams(N=a.N, ratio=a.ratio)
    runs = {}
    for n in a.names:
        runs[n] = run_pattern(n, P, NP, force=a.force)
    make_outputs(a.names, runs, P, NP)
    rep = {n: dict(wall_s=float(runs[n]["wall"]), steps=int(runs[n]["steps"]),
                   s_per_phys_s=float(runs[n]["wall"]) / float(runs[n]["T"]),
                   max_ledger_err=float(np.max(runs[n]["ledger"][:, 1:5]))) for n in a.names}
    with open(os.path.join(pattern_dir(), "timing.json"), "w") as f:
        json.dump(dict(params_digest=P.digest(), params=P.to_dict(), numerics=NP.to_dict(), runs=rep), f, indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
