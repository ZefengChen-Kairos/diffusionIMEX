#!/bin/bash
# Run the test suite of src/SPB6_explicit.cpp (E) and src/SPB6_TV_IMEX.cpp (I), and print the statistics line of each run.
# All runs use the same scheme settings: explicit = acoustic CFL 0.4, IMEX = material CFL 0.2 with LSDIRK2 (rk=2).
# usage: bash tests/run_suite.sh OUTDIR [mech|relax|all]   (4 jobs in parallel; metrics: python3 tests/metrics.py CASE files)
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/suite}; WHAT=${2:-all}
mkdir -p "$OUT"
EX="$H/src/spb6_ex"; IM="$H/src/spb6_imex"
g++ -O2 -std=c++17 "$H/src/SPB6_explicit.cpp" -o "$EX" || exit 1
g++ -O2 -std=c++17 "$H/src/SPB6_TV_IMEX.cpp" -o "$IM" || exit 1
jobs_list=()
add() { jobs_list+=("$1|$2|$3"); }   # name|E or I|args
both() { add "${1}_ex" E "$2"; add "${1}_im" I "$2"; }
if [ "$WHAT" = mech ] || [ "$WHAT" = all ]; then
  both c4  "case=4 N=1000"            # SPB 4.1 interface advection (WB)
  both c5  "case=5 N=1000"            # SPB 4.3.1 water-air 1 GPa
  both c6  "case=6 N=1000"            # SPB 4.3.2 water-air 1 TPa
  both c0  "case=0 N=1000"            # SPB 4.5 cavitation (vacuum)
  both c3  "case=3 N=1000 relax=0"    # P-S dodecane, p-relaxation only
  both c12 "case=12 N=500"            # SPB 5.1 epoxy-spinel
  both c13 "case=13 N=500"            # SPB 5.2 epoxy-spinel
fi
if [ "$WHAT" = relax ] || [ "$WHAT" = all ]; then
  for r in 0 1 2 3; do
    both c1_r$r "case=1 N=1000 relax=$r"   # P-S water cavitation tube +-2 m/s (p, pT, pTg, p-pTg)
    both c2_r$r "case=2 N=1000 relax=$r"   # P-S water cavitation tube +-500 m/s
    both c3_r$r "case=3 N=1000 relax=$r"   # P-S dodecane liquid-vapour shock tube
  done
  both p30_p    "case=30 N=1000 relax=0"              # Pelanti 2022 8.2.1, five relaxation levels
  both p30_pTf  "case=30 N=1000 relax=1 theta=2000"
  both p30_pT   "case=30 N=1000 relax=1"
  both p30_pTgf "case=30 N=1000 relax=2 nu=1e-4"
  both p30_pTg  "case=30 N=1000 relax=2"
fi
run_one() { local name=${1%%|*} rest=${1#*|}; local kind=${rest%%|*} args=${rest#*|}; local bin=$IM; [ "$kind" = E ] && bin=$EX
  (cd "$OUT" && "$bin" $args out="$OUT/$name.dat" > "$OUT/$name.log" 2>&1); }
export -f run_one; export OUT EX IM
printf '%s\n' "${jobs_list[@]}" | xargs -P 4 -I{} bash -c 'run_one "$@"' _ {}
for j in "${jobs_list[@]}"; do n=${j%%|*}; echo "$n: $(tail -1 "$OUT/$n.log" | cut -c1-170)"; done
