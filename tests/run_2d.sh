#!/bin/bash
# 2D paper tests (V13): explicit (E) and IMEX (I) on the grids of the papers (the defaults set in src/spb6_physics.h test_setup),
# profiles written at the times of the paper figures (snap=) and, for the underwater explosion, the wall pressure history (probe=).
# usage: bash tests/run_2d.sh OUTDIR [REGEX] [extra args, e.g. N=300 ny=100]   (jobs run one after another; each uses all OpenMP threads)
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/2d}; SEL=${2:-.}; shift $(( $# > 2 ? 2 : $# )); EXTRA="$*"
mkdir -p "$OUT"
EX="$H/src/spb6_ex"; IM="$H/src/spb6_imex"
g++ -O2 -std=c++17 -fopenmp "$H/src/SPB6_explicit.cpp" -o "$EX" || exit 1
g++ -O2 -std=c++17 -fopenmp -I/usr/include/eigen3 "$H/src/SPB6_TV_IMEX.cpp" -o "$IM" || exit 1
jobs_list=()
both() { jobs_list+=("${1}_E|E|$2" "${1}_I|I|$2"); }
both he    "case=70 snap=134e-6,274e-6,344e-6,554e-6"                  # SPB 2009 Figs. 13-14: shock / helium bubble
both rm    "case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3"             # SPB 2009 Figs. 15-17: cavitating Richtmyer-Meshkov
both bub   "case=72 snap=0.4e-3,0.8e-3,1.2e-3"                         # P-S 2014 Fig. 5: dodecane vapour bubble compression (relax 2)
both undex "case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6"   # Pelanti 2022 Figs. 2-3: underwater explosion near a wall
for j in "${jobs_list[@]}"; do
  n=${j%%|*}; rest=${j#*|}; k=${rest%%|*}; a=${rest#*|}; b=$IM; [ "$k" = E ] && b=$EX
  echo "$n" | grep -qE "^($SEL)" || continue
  [ -s "$OUT/$n.dat" ] && continue
  (cd "$OUT" && "$b" $a prec=8 $EXTRA out="$OUT/$n.dat" > "$OUT/$n.log" 2>&1)
  echo "$n: $(tail -1 "$OUT/$n.log" | cut -c1-170)"
done
