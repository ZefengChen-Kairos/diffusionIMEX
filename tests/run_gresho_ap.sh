#!/bin/bash
# Gresho vortex, low-Mach AP check (B7 Theorem 3, devlog V11 §11.8).
# IMEX: cases 50-55 (M = 1e-1 ... 1e-6), the same fixed dt = 0.0025 for every M (no start-up ramp), 80 x 80, one turn.
# Explicit (HLLC, acoustic CFL 0.4): cases 50-51 only (M = 0.01 already needs 5e4 steps).
# usage: bash tests/run_gresho_ap.sh OUTDIR [noex]   then   python3 tests/plot_gresho_ap.py OUTDIR docs/figures/gresho_ap.png
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/gresho_ap}; WHAT=${2:-all}
mkdir -p "$OUT"
EX="$H/src/spb6_ex"; IM="$H/src/spb6_imex"
g++ -O2 -std=c++17 "$H/src/SPB6_explicit.cpp" -o "$EX" || exit 1
g++ -O2 -std=c++17 -I/usr/include/eigen3 "$H/src/SPB6_TV_IMEX.cpp" -o "$IM" || exit 1
jobs_list=()
for c in 50 51 52 53 54 55; do jobs_list+=("im_$c|I|case=$c N=80 ny=80 dt=0.0025"); done
if [ "$WHAT" != noex ]; then for c in 50 51; do jobs_list+=("ex_$c|E|case=$c N=80 ny=80"); done; fi
run_one() { local name=${1%%|*} rest=${1#*|}; local kind=${rest%%|*} args=${rest#*|}; local bin=$IM; [ "$kind" = E ] && bin=$EX
  (cd "$OUT" && "$bin" $args out="$OUT/$name.dat" > "$OUT/$name.log" 2>&1); }
export -f run_one; export OUT EX IM
printf '%s\n' "${jobs_list[@]}" | xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
for j in "${jobs_list[@]}"; do n=${j%%|*}; echo "$n: $(tail -1 "$OUT/$n.log" | cut -c1-170)"; done
