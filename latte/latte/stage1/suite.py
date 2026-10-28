"""Stage-1 verification suite runner.

    python -m latte.stage1.suite [--only T02,T05] [--out runs/stage1] [--report reports/stage1]
"""
import argparse
import glob
import json
import os
import shutil
import sys
import time
import traceback

import yaml

from ..runio import Run, _jsonable, model_id, git_info
from . import t_inlet, t_transport, t_dynamics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

TESTS = [
    ("T01", None, "V0.5 legacy volumetric-source analytic check; V0.5 baseline excluded by user decision (2026-10-28)"),
    ("T02", t_inlet.T02, "inlet flux integral, moving section"),
    ("T03", t_inlet.T03, "inlet input contract"),
    ("T04", None, "V0.5 rename-equivalence check; V0.5 baseline excluded by user decision (2026-10-28)"),
    ("T05", t_transport.T05, "static stripes (full solver)"),
    ("T06", t_transport.T06, "periodic translation"),
    ("T07", t_transport.T07, "time-reversed vortex"),
    ("T08", t_transport.T08, "piston displacement"),
    ("T09", t_transport.T09, "component simplex"),
    ("T10", t_transport.T10, "pure vertical remap"),
    ("T11", t_dynamics.T11, "hydrostatic balance"),
    ("T12", t_dynamics.T12, "common pressure projection"),
    ("T13", t_dynamics.T13, "two-layer viscous exchange"),
    ("T14", t_transport.T14, "interlayer material exchange"),
    ("T15", t_dynamics.T15, "injection and liquid level"),
    ("T16", t_dynamics.T16, "density / buoyancy"),
    ("T17", t_dynamics.T17, "surface / wall / inlet switch-off"),
    ("T17b", "blocked", "cup-edge meniscus / contact line: contact-angle and slip model not chosen (V3 Sec. 6, to seal before real-cup tests)"),
    ("W01", t_dynamics.W01, "extra: non-hydrostatic gravity-wave dispersion"),
    ("W02", t_dynamics.W02, "extra: capillary-gravity wave dispersion"),
]


def _short(v):
    if v is None:
        return "null"
    if isinstance(v, float):
        return f"{v:.2e}"
    if isinstance(v, list):
        return "[" + ", ".join(_short(x) for x in v) + "]"
    return str(v)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(ROOT, "configs", "stage1", "suite.yaml"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--out", default=os.path.join(ROOT, "runs", "stage1"))
    ap.add_argument("--report", default=os.path.join(ROOT, "reports", "stage1"))
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    only = set(a.only.split(",")) if a.only else None
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join(a.report, "figures"), exist_ok=True)
    summary = []
    for tid, fn, title in TESTS:
        if only and tid not in only:
            continue
        sec = {"test": cfg["tests"].get(tid.rstrip("b")), "model": cfg["model"], "numerics": cfg["numerics"],
               "thresholds": cfg["thresholds"], "seed": cfg["seed"], "suite_version": cfg["suite_version"]}
        run = Run(a.out, tid, sec, cfg["seed"])
        t0 = time.time()
        if fn is None or fn == "blocked":
            rec = run.finish({"status": "not_run" if fn is None else "blocked", "verdict": "not_applicable", "reason": title})
        else:
            try:
                res = fn(cfg, run)
            except Exception as e:                       # record the first failure with traceback
                with open(run.path("failure", "traceback.txt"), "w") as f:
                    traceback.print_exc(file=f)
                res = {"status": "failed", "verdict": "fail", "reason": f"{type(e).__name__}: {e}"}
            res["wall_seconds"] = round(time.time() - t0, 2)
            rec = run.finish(res)
        rec["title"] = title
        rec["run_dir"] = os.path.relpath(run.dir, ROOT)
        rec["figures"] = sorted(glob.glob(run.path("figures", "*.png")))
        for fpath in rec["figures"]:
            shutil.copy(fpath, os.path.join(a.report, "figures", os.path.basename(fpath)))
        summary.append(rec)
        print(f"{tid:4s} {rec['status']:9s} {rec['verdict']:14s} {rec.get('wall_seconds', 0):7.1f}s  {rec['reason']}", flush=True)
    keep = ("test_id", "title", "status", "verdict", "reason", "mass_residual", "volume_residual", "simplex_error",
            "bounds_violation", "divergence_error", "solver_residual", "accepted_steps", "rejected_steps", "wall_seconds", "run_dir",
            "figures")
    for r in summary:
        r["figures"] = [os.path.relpath(f, ROOT) for f in r["figures"]]
    sp = os.path.join(a.report, "summary.json")
    if only and os.path.exists(sp):                     # partial rerun: merge into the existing summary
        old = {r["test_id"]: r for r in json.load(open(sp))["results"]}
        new = {r["test_id"]: r for r in summary}
        order = [t for t, *_ in TESTS]
        summary = [new.get(t, old.get(t)) for t in order if (t in new or t in old)]
    with open(os.path.join(a.report, "summary.json"), "w") as f:
        json.dump(_jsonable({"model_id": model_id(), "git": git_info(), "config": os.path.relpath(a.config, ROOT),
                             "results": [{k: r.get(k) for k in keep} for r in summary]}), f, indent=2)
    lines = ["# Stage-1 verification results (auto-generated)", "",
             f"model_id `{model_id()}`, config `{os.path.relpath(a.config, ROOT)}`, git `{git_info().get('commit')}`", "",
             "| test | what | status | verdict | key result | volume res. | mass res. | div. err. | steps (rej.) | wall s |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in summary:
        steps = "" if r.get("accepted_steps") is None else f"{r['accepted_steps']} ({r.get('rejected_steps') or 0})"
        lines.append(f"| {r['test_id']} | {r['title']} | {r['status']} | **{r['verdict']}** | {r['reason']} | "
                     f"{_short(r.get('volume_residual'))} | {_short(r.get('mass_residual'))} | {_short(r.get('divergence_error'))} | "
                     f"{steps} | {r.get('wall_seconds', '')} |")
    lines += ["", "## Figures", ""]
    for r in summary:
        for fpath in r.get("figures") or []:
            lines.append(f"**{r['test_id']}** `{os.path.basename(fpath)}`  ")
            lines.append(f"![{r['test_id']}](figures/{os.path.basename(fpath)})")
            lines.append("")
    with open(os.path.join(a.report, "RESULTS.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    bad = [r["test_id"] for r in summary if r["verdict"] == "fail"]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
