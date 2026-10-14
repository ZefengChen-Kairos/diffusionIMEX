#!/bin/bash
# Archive runs (V14): IMEX (src/SPB6_TV_IMEX.cpp) on every compared paper test, on the paper grid, with the per-test settings of
# src/spb6_physics.h test_setup (material CFL, Rusanov speed, relaxation level, boundaries). Raw outputs go to OUTDIR;
# tests/make_archive.py then writes data/<test>/ (.dat files + note.md).
# usage: bash tests/run_archive.sh OUTDIR [1d|2d|all] [REGEX]   (existing outputs are kept; 1D jobs 4 in parallel, 2D jobs one at a time on 4 threads)
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/archive}; WHAT=${2:-all}; SEL=${3:-.}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
IM="$H/src/spb6_imex"
g++ -O2 -std=c++17 -fopenmp -I/usr/include/eigen3 "$H/src/SPB6_TV_IMEX.cpp" -o "$IM" || exit 1
J1=(
  "spb41|case=4 N=200"                                   # SPB 2009 4.1, Fig. 5
  "spb42|case=60 N=500"                                  # SPB 2009 4.2, Figs. 6-7
  "spb431|case=5 N=1000"                                 # SPB 2009 4.3.1, Fig. 8
  "spb431_n100|case=5 N=100"                             # SPB 2009 4.3.1, Fig. 9
  "spb432|case=6 N=1000"                                 # SPB 2009 4.3.2, Fig. 10
  "spb45|case=0 N=1000"                                  # SPB 2009 4.5, Fig. 12
  "spb51|case=12 N=500"                                  # SPB 2009 5.1, Fig. 19
  "spb52|case=13 N=500"                                  # SPB 2009 5.2, Fig. 20
  "pel811|case=61 N=500"                                 # Pelanti 2022 8.1.1, Fig. 1
  "ps1_r0|case=1 N=5000 relax=0" "ps1_r1|case=1 N=5000 relax=1" "ps1_r2|case=1 N=5000 relax=2" "ps1_r3|case=1 N=5000 relax=3"   # P-S 2014 Fig. 1
  "ps2_r0|case=2 N=5000 relax=0" "ps2_r1|case=2 N=5000 relax=1" "ps2_r2|case=2 N=5000 relax=2" "ps2_r3|case=2 N=5000 relax=3"   # P-S 2014 Fig. 2
  "ps3_r0|case=3 N=2000 relax=0" "ps3_r2|case=3 N=2000 relax=2" "ps3_r3|case=3 N=2000 relax=3"                                # P-S 2014 Fig. 3
  "pel821_p|case=30 N=5000 relax=0" "pel821_pTf|case=30 N=5000 relax=1 theta=2000" "pel821_pT|case=30 N=5000 relax=1"
  "pel821_pTgf|case=30 N=5000 relax=2 nu=1e-4" "pel821_pTg|case=30 N=5000 relax=2"                                           # Pelanti 2022 Figs. 4-5
  "pel822_r0|case=31 N=2000 relax=0" "pel822_r2|case=31 N=2000 relax=2"                                                    # Pelanti 2022 Fig. 6
)
J2=(
  "he|case=70 snap=134e-6,274e-6,344e-6,554e-6"                              # SPB 2009 4.6, Figs. 13-14
  "rm|case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3"                         # SPB 2009 4.7 / P-S 2014 5.3, Figs. 15-17
  "bub|case=72 snap=0.4e-3,0.8e-3,1.2e-3"                                    # P-S 2014 5.4, Fig. 5
  "undex|case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6"  # Pelanti 2022 8.1.2, Figs. 2-3
)
run_one() { local name=${1%%|*} args=${1#*|}
  [ -s "$OUT/$name.dat" ] && return 0
  (cd "$OUT" && "$IM" $args out="$OUT/$name.dat" > "$OUT/$name.log" 2>&1); }
export -f run_one; export OUT IM
if [ "$WHAT" != 2d ]; then
  printf '%s\n' "${J1[@]}" | grep -E "^($SEL)" | OMP_NUM_THREADS=1 xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
fi
if [ "$WHAT" != 1d ]; then
  for j in "${J2[@]}"; do echo "$j" | grep -qE "^($SEL)" || continue; OMP_NUM_THREADS=4 run_one "$j"; done
fi
for j in "${J1[@]}" "${J2[@]}"; do n=${j%%|*}; [ -f "$OUT/$n.log" ] && echo "$n: $(head -1 "$OUT/$n.log" | cut -c1-150)"; done
