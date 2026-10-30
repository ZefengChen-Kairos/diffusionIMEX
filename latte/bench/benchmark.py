"""Hot-spot profile and backend speed table.

    python bench/benchmark.py [--steps 300] [--N 128]

Writes <out>/bench/benchmark.md and benchmark.json: per-operator share of the step time for the
numpy reference and the numba backend, ms/step, steps/s, s per physical second, speedup,
plus the environment self-check (GPU / Warp availability; the Warp backend is phase two).
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import json
import time
import cProfile
import pstats
import argparse
import platform
import resource

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                   # noqa: E402
from core.params import PhysParams, NumParams         # noqa: E402
from core.solver_np import Sim                        # noqa: E402
from actions.scripts import ALL                       # noqa: E402
from gates.common import out_root                     # noqa: E402

OPS = {
    "rhs": "coarse fluxes (MUSCL + pressure + kappa), x2 per step (RK2)",
    "_rhs_nb": "coarse fluxes (numba kernel), x2 per step",
    "fine_advect": "fine material BFECC (numpy)",
    "_fine_advect_nb": "fine material BFECC (numba kernels)",
    "fine_sources": "fine source mixing in windows",
    "sources": "coarse D/E/R/tau_j with windowed kernels",
    "kernel": "kernel evaluation (exp on windows)",
    "prims": "primitive variables / pressure",
    "check": "conservation assertion",
    "dt_cfl": "CFL time step",
}


def env_report():
    rep = dict(python=platform.python_version(), numpy=np.__version__, machine=platform.machine(),
               cpus=os.cpu_count())
    try:
        import subprocess
        out = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=10)
        rep["nvidia_smi"] = out.stdout.strip() or out.stderr.strip()
    except Exception as e:
        rep["nvidia_smi"] = f"not available ({e.__class__.__name__})"
    try:
        import warp as wp
        wp.init()
        rep["warp"] = wp.config.version
        rep["warp_devices"] = [str(d) for d in wp.get_devices()]
    except Exception as e:
        rep["warp"] = f"not installed ({e.__class__.__name__}); GPU backend is phase two"
    return rep


def profile(backend, N, steps, warm=200):
    P = PhysParams()
    NP = NumParams(N=N)
    s = Sim(ALL["leaf"](), P=P, NP=NP, backend=backend)
    for _ in range(warm):          # into the wobble phase; also JIT warm-up
        s.step()
    pr = cProfile.Profile()
    t0 = time.perf_counter()
    t_phys0 = s.t
    pr.enable()
    for _ in range(steps):
        s.step()
    pr.disable()
    wall = time.perf_counter() - t0
    st = pstats.Stats(pr)
    rows = []
    tot = sum(v[3] for k, v in st.stats.items() if k[2] == "step")
    for (fn, line, name), v in st.stats.items():
        if name in OPS and "latte" in fn:
            rows.append(dict(op=name, desc=OPS[name], calls=v[1], cum_s=v[3], share=v[3] / tot))
    rows.sort(key=lambda r: -r["share"])
    return dict(backend=backend, N=N, steps=steps, wall_s=wall, ms_per_step=1e3 * wall / steps,
                steps_per_s=steps / wall, s_per_phys_s=wall / (s.t - t_phys0), ops=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--N", type=int, default=128)
    a = ap.parse_args()
    d = os.path.join(out_root(), "bench")
    os.makedirs(d, exist_ok=True)
    env = env_report()
    res = [profile("numpy", a.N, max(a.steps // 3, 50)), profile("numba", a.N, a.steps)]
    base = res[0]["ms_per_step"]
    lines = ["# Benchmark (phase one: CPU reference)", "", "## Environment", "", "```",
             json.dumps(env, indent=1), "```", "",
             "## Speed", "", "| backend | N | ms/step | steps/s | s per physical s | speedup vs numpy |",
             "|---|---|---|---|---|---|"]
    for r in res:
        lines.append(f"| {r['backend']} | {r['N']} | {r['ms_per_step']:.2f} | {r['steps_per_s']:.1f} | "
                     f"{r['s_per_phys_s']:.1f} | {base / r['ms_per_step']:.2f}x |")
    lines.append("| warp (GPU) | - | phase two | | | |")
    for r in res:
        lines += ["", f"## Hot spots, {r['backend']} backend (cumulative share of step time)", "",
                  "| operator | what | calls | cum. s | share |", "|---|---|---|---|---|"]
        for o in r["ops"]:
            lines.append(f"| {o['op']} | {o['desc']} | {o['calls']} | {o['cum_s']:.3f} | {100 * o['share']:.1f} % |")
    lines += ["", f"peak RSS: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB"]
    txt = "\n".join(lines) + "\n"
    open(os.path.join(d, "benchmark.md"), "w").write(txt)
    json.dump(dict(env=env, results=res), open(os.path.join(d, "benchmark.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
