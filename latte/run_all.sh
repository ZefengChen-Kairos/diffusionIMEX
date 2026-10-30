#!/usr/bin/env bash
# One-shot: environment -> G0 gate -> six patterns + board -> remaining gates -> benchmark.
# Headless (matplotlib Agg / PIL only).  Outputs: /root/autodl-tmp/latte_out (or ./out).
set -e
cd "$(dirname "$0")"
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
if [ -d /root/autodl-tmp ]; then export LATTE_OUT=${LATTE_OUT:-/root/autodl-tmp/latte_out}; else export LATTE_OUT=${LATTE_OUT:-$(pwd)/out}; fi
mkdir -p "$LATTE_OUT"
[ "${SKIP_SETUP:-0}" = "1" ] || bash setup.sh
echo "== G0 (must pass before anything else)"
python3 -m gates.run_gates G0
python3 - <<'PY'
import json, os, sys
d = json.load(open(os.path.join(os.environ["LATTE_OUT"], "gates", "G0", "metrics.json")))
sys.exit(0 if d["passed"] else 1)
PY
echo "== pour inputs (continuous pitcher motion; feasibility checked)"
python3 -m actions.plot "$LATTE_OUT/pours.png"
echo "== six patterns (shared parameters) + comparison board"
python3 run_patterns.py
echo "== gates G1..G8 + counterfactuals"
python3 -m gates.run_gates G1 G2 G3 G4 G5 G6 G7 G8 CF G9
echo "== sensitivity table"
python3 -m gates.sensitivity
echo "== benchmark"
python3 bench/benchmark.py
echo "done: $LATTE_OUT"
