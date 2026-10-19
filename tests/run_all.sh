#!/bin/bash
# Every result of the paper (V22), resumable: phase A = tests/jobs_all.txt (1D paper tests, Bartak, Taylor-Green, two-phase vortex,
# Schropff Riemann problems; 4 jobs x 1 thread), phase B = 2D on the grids of the original papers (SPB 4.7: 900 x 400; Pelanti 8.1.2: 721 x 400; one job at a time, 4 threads).
# Lines of jobs_all.txt: SUBDIR|NAME|ARGUMENTS.  Outputs: OUTDIR/SUBDIR/NAME.dat/.log; a job with a finished log is skipped.
# usage: bash tests/run_all.sh [OUTDIR=run/v22]
#   then: python3 tests/make_archive.py OUTDIR data; python3 tests/archive_extra.py OUTDIR data
H=$(cd "$(dirname "$0")/.." && pwd); R=${1:-$H/run/v22}; mkdir -p "$R"; R=$(cd "$R" && pwd); BIN=$H/src/spb6_imex
g++ -O2 -std=c++17 -fopenmp -I/usr/include/eigen3 "$H/src/SPB6_TV_IMEX.cpp" -o "$BIN" || exit 1
one() { IFS='|' read -r d n a <<< "$1"; mkdir -p "$R/$d"; cd "$R/$d" || return 1
  [ -s "$n.log" ] && grep -q "^case" "$n.log" && return 0
  OMP_NUM_THREADS=${OMPT:-1} "$BIN" $a out="$PWD/$n.dat" > "$n.log" 2>&1; }
export -f one; export R BIN
grep -v '^#' "$H/tests/jobs_all.txt" | xargs -P 4 -I{} bash -c 'one "$@"' _ {}
echo PHASE_A_DONE
for j in ".|he|case=70 snap=134e-6,274e-6,344e-6,554e-6" ".|rm|case=71 N=900 ny=400 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3" \
         ".|bub|case=72 snap=0.4e-3,0.8e-3,1.2e-3" ".|undex|case=73 N=721 ny=400 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6"; do
  OMPT=4 one "$j"; done
echo ALL_DONE
