#!/bin/bash
# Run the test suite of src/SPB6_TV_IMEX.cpp and print metrics against exact solutions.
# usage: bash tests/run_suite.sh OUTDIR [mech|relax|all]   (runs 4 jobs in parallel)
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/suite}; WHAT=${2:-all}
mkdir -p "$OUT"
BIN="$H/src/spb6"
g++ -O2 -std=c++17 "$H/src/SPB6_TV_IMEX.cpp" -o "$BIN" || exit 1
jobs_list=()
add() { jobs_list+=("$1|$2"); }   # name|args
if [ "$WHAT" = mech ] || [ "$WHAT" = all ]; then
  add c4_imex   "case=4 N=1000 mode=1 rk=2"
  add c5_ex     "case=5 N=1000 mode=0 cfl=0.4"
  add c5_rk1    "case=5 N=1000 mode=1 rk=1"
  add c5_rk2    "case=5 N=1000 mode=1 rk=2"
  add c5_rk3    "case=5 N=1000 mode=1 rk=3"
  add c5_rk2c05 "case=5 N=1000 mode=1 rk=2 cfl=0.05"
  add c6_ex     "case=6 N=1000 mode=0 cfl=0.4"
  add c6_rk2    "case=6 N=1000 mode=1 rk=2"
  add c0_ex     "case=0 N=1000 mode=0 cfl=0.4"
  add c0_rk1    "case=0 N=1000 mode=1 rk=1"
  add c0_rk2    "case=0 N=1000 mode=1 rk=2"
  add c3_ex     "case=3 N=1000 mode=0 cfl=0.4 relax=0"
  add c3_rk1    "case=3 N=1000 mode=1 rk=1 relax=0"
  add c3_rk2    "case=3 N=1000 mode=1 rk=2 relax=0"
  add c12_ex    "case=12 N=500 mode=0 cfl=0.4"
  add c12_rk2   "case=12 N=500 mode=1 rk=2"
  add c13_ex    "case=13 N=500 mode=0 cfl=0.4"
  add c13_rk2   "case=13 N=500 mode=1 rk=2"
fi
if [ "$WHAT" = relax ] || [ "$WHAT" = all ]; then
  for r in 0 1 2 3; do
    add c1_ex_r$r  "case=1 N=1000 mode=0 cfl=0.4 relax=$r"
    add c1_rk2_r$r "case=1 N=1000 mode=1 rk=2 relax=$r"
    add c2_ex_r$r  "case=2 N=1000 mode=0 cfl=0.4 relax=$r"
    add c2_rk2_r$r "case=2 N=1000 mode=1 rk=2 relax=$r"
    add c3_ex_r$r  "case=3 N=1000 mode=0 cfl=0.4 relax=$r"
    add c3_rk2_r$r "case=3 N=1000 mode=1 rk=2 relax=$r"
  done
fi
run_one() { local name=${1%%|*} args=${1#*|}; (cd "$OUT" && "$BIN" $args out="$OUT/$name.dat" > "$OUT/$name.log" 2>&1); }
export -f run_one; export OUT BIN
printf '%s\n' "${jobs_list[@]}" | xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
for j in "${jobs_list[@]}"; do n=${j%%|*}; echo "$n: $(tail -1 "$OUT/$n.log" | cut -c1-170)"; done
