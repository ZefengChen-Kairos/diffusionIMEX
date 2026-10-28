"""Run directories, manifests and metric records (V3 Sec. 15)."""
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys

import numpy as np

REVIEW_FIELDS = ("test_id", "model_id", "input_id", "status", "verdict", "reason", "mass_residual",
                 "volume_residual", "simplex_error", "bounds_violation", "divergence_error",
                 "mixture_indicator", "solver_residual", "accepted_steps", "rejected_steps")
STATUS = ("not_run", "blocked", "completed", "failed")
VERDICT = ("pass", "fail", "review", "not_applicable")

PKG_DIR = os.path.dirname(os.path.abspath(__file__))


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating,)):
        x = float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, float) and not np.isfinite(x):
        return str(x)
    return x


def source_hashes():
    out = {}
    for root, _, files in os.walk(PKG_DIR):
        for f in sorted(files):
            if f.endswith(".py"):
                p = os.path.join(root, f)
                out[os.path.relpath(p, PKG_DIR)] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return out


def model_id():
    h = hashlib.sha256(json.dumps(source_hashes(), sort_keys=True).encode()).hexdigest()[:12]
    return f"latte-fvm-stage1-{h}"


def config_hash(cfg):
    return hashlib.sha256(json.dumps(_jsonable(cfg), sort_keys=True).encode()).hexdigest()[:12]


def git_info():
    try:
        root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=PKG_DIR, text=True).strip()
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain", "--", "latte"], cwd=root, text=True).strip())
        return {"commit": commit, "latte_dirty": dirty}
    except Exception as e:                                  # pragma: no cover
        return {"commit": None, "reason": str(e)}


def environment():
    import scipy
    try:
        import pyamg
        pv = pyamg.__version__
    except Exception:                                       # pragma: no cover
        pv = None
    return {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
            "pyamg": pv, "platform": platform.platform(), "processor": platform.processor(),
            "cpu_count": os.cpu_count(), "precision": "float64"}


class Run:
    def __init__(self, root, test_id, cfg, seed=None):
        self.test_id = test_id
        self.cfg = cfg
        self.input_id = config_hash(cfg)
        self.run_id = f"{test_id}_{self.input_id}"
        self.dir = os.path.join(root, self.run_id)
        for sub in ("diagnostics", "figures", "failure"):
            os.makedirs(os.path.join(self.dir, sub), exist_ok=True)
        self.seed = seed
        man = {"run_id": self.run_id, "test_id": test_id, "model_id": model_id(), "input_id": self.input_id,
               "seed": seed, "git": git_info(), "environment": environment(), "source_sha256": source_hashes()}
        self._write("manifest.json", man)
        import yaml
        with open(os.path.join(self.dir, "config_resolved.yaml"), "w") as f:
            yaml.safe_dump(_jsonable(cfg), f, sort_keys=False)

    def _write(self, name, obj):
        with open(os.path.join(self.dir, name), "w") as f:
            json.dump(_jsonable(obj), f, indent=2)

    def path(self, *p):
        return os.path.join(self.dir, *p)

    def write_log(self, rows, name="diagnostics/solver.csv"):
        if not rows:
            return
        keys = sorted({k for r in rows for k in r})
        with open(self.path(name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: _jsonable(r.get(k)) for k in keys})

    def finish(self, result):
        rec = {k: None for k in REVIEW_FIELDS}
        rec.update({"test_id": self.test_id, "model_id": model_id(), "input_id": self.input_id})
        rec.update(result)
        nulls = rec.setdefault("null_reasons", {})
        for k in REVIEW_FIELDS:
            if rec.get(k) is None and k not in nulls:
                nulls[k] = "not applicable to this test"
        assert rec["status"] in STATUS and rec["verdict"] in VERDICT, (rec["status"], rec["verdict"])
        self._write("metrics.json", rec)
        return rec
