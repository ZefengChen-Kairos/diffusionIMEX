"""Run the gates in order G0 -> G8 + counterfactuals; write out/gates/summary.md.

    python -m gates.run_gates [G0 G1 ...] [--quick]
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gates.common import out_root                     # noqa: E402

ORDER = ["G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "CF", "G9"]


def run(name, quick=False):
    if name == "G0":
        from gates.g0 import run_g0
        return run_g0(full=not quick)
    if name == "G1":
        from gates.g1 import run_g1
        return run_g1()
    if name == "G2":
        from gates.g2 import run_g2
        return run_g2()
    if name == "G3":
        from gates.g3 import run_g3
        return run_g3()
    if name in ("G4", "G5", "G6", "G7"):
        from gates import g4_7
        return getattr(g4_7, "run_" + name.lower())()
    if name == "G8":
        from gates.g8_cf import run_g8
        return run_g8()
    if name == "G9":
        from gates.g9 import run_g9
        return run_g9()
    if name == "CF":
        from gates.g8_cf import run_cf
        return run_cf()
    raise KeyError(name)


def summary():
    root = os.path.join(out_root(), "gates")
    lines = ["# Gate summary", "", "| gate | check | value | threshold | result |", "|---|---|---|---|---|"]
    for g in ORDER:
        p = os.path.join(root, g, "metrics.json")
        if not os.path.exists(p):
            lines.append(f"| {g} | (not run) | | | |")
            continue
        d = json.load(open(p))
        for c in d["checks"]:
            v = c["value"]
            if isinstance(v, float):
                v = f"{v:.4g}"
            elif isinstance(v, list):
                v = "[" + ", ".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in v) + "]"
            elif isinstance(v, dict):
                v = "{" + ", ".join(f"{k}: {x:.3g}" if isinstance(x, float) else f"{k}: {x}" for k, x in v.items()
                                    if not isinstance(x, (list, dict))) + "}"
            lines.append(f"| {g} | {c['id']} {c['desc']} | {v} | {c['threshold']} | "
                         f"{'PASS' if c['PASS'] else '**FAIL**'} |")
    with open(os.path.join(root, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return "\n".join(lines)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    quick = "--quick" in sys.argv
    names = args or ORDER
    for n in names:
        t0 = time.time()
        print(f"=== {n}", flush=True)
        try:
            R = run(n, quick)
            print(f"=== {n}: {'PASS' if R.passed else 'FAIL'} ({time.time() - t0:.0f} s)", flush=True)
        except Exception as e:   # keep going; the failure is recorded
            import traceback
            traceback.print_exc()
            os.makedirs(os.path.join(out_root(), "gates", n), exist_ok=True)
            with open(os.path.join(out_root(), "gates", n, "metrics.json"), "w") as f:
                json.dump(dict(gate=n, passed=False, checks=[dict(id=n, desc=f"exception: {e!r}", value="",
                                                                    threshold="", PASS=False)]), f)
    print(summary())


if __name__ == "__main__":
    main()
