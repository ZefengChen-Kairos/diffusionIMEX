#!/bin/bash
# Runs for the reproductions of the paper figures (tests/plot_paper_figs.py): explicit (E, src/SPB6_explicit.cpp) and
# IMEX (I, src/SPB6_TV_IMEX.cpp) on the grid of each paper figure. The scheme settings are those of tests/run_suite.sh
# (explicit: acoustic CFL 0.4; IMEX: material CFL 0.2, LSDIRK2); only the grid and the relaxation level change.
# usage: bash tests/run_paper_figs.sh OUTDIR [NPAR] [REGEX]   (existing outputs are kept; REGEX selects job names)
set -u
H=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$H/run/figs}; NPAR=${2:-4}; SEL=${3:-.}
mkdir -p "$OUT"
EX="$H/src/spb6_ex"; IM="$H/src/spb6_imex"
g++ -O2 -std=c++17 "$H/src/SPB6_explicit.cpp" -o "$EX" || exit 1
g++ -O2 -std=c++17 -I/usr/include/eigen3 "$H/src/SPB6_TV_IMEX.cpp" -o "$IM" || exit 1
jobs_list=()
both() { jobs_list+=("${1}_E|E|$2" "${1}_I|I|$2"); }
both spb5  "case=4 N=200"                                          # SPB 2009 Fig. 5: interface advection
both spb8  "case=5 N=1000"                                         # Fig. 8: water-air shock tube, 1 GPa
both spb10 "case=6 N=1000"                                         # Fig. 10: water-air shock tube, 1 TPa
both spb12 "case=0 N=1000"                                         # Fig. 12: expansion tube (cavitation)
both spb19 "case=12 N=500"                                         # Fig. 19: epoxy-spinel, 1e10 Pa
both spb20 "case=13 N=500"                                         # Fig. 20: epoxy-spinel, 2e11 Pa
for r in 0 2 3; do both ps3_r$r "case=3 N=2000 relax=$r"; done     # P-S 2014 Fig. 3: dodecane (SG), t = 473 us
for r in 0 2; do both pel6_r$r "case=31 N=2000 relax=$r"; done     # Pelanti 2022 Fig. 6: dodecane (NASG), t = 473 us
for r in 0 1 2 3; do both ps2_r$r "case=2 N=5000 relax=$r"; done   # P-S 2014 Fig. 2: +-500 m/s, t = 0.58 ms
both pel4_p    "case=30 N=5000 relax=0"                             # Pelanti 2022 Figs. 4-5 (8.2.1), t = 3 ms (explicit: ~40 min each)
both pel4_pTf  "case=30 N=5000 relax=1 theta=2000"
both pel4_pT   "case=30 N=5000 relax=1"
both pel4_pTgf "case=30 N=5000 relax=2 nu=1e-4"
both pel4_pTg  "case=30 N=5000 relax=2"
for r in 0 1 2 3; do both ps1_r$r "case=1 N=5000 relax=$r"; done   # P-S 2014 Fig. 1: water cavitation tube +-2 m/s, t = 3.2 ms (explicit: ~45 min each)
run_one() { local name=${1%%|*} rest=${1#*|}; local kind=${rest%%|*} args=${rest#*|}; local bin=$IM; [ "$kind" = E ] && bin=$EX
  [ -s "$OUT/$name.dat" ] && return 0
  (cd "$OUT" && "$bin" $args out="$OUT/$name.dat" > "$OUT/$name.log" 2>&1); }
export -f run_one; export OUT EX IM
printf '%s\n' "${jobs_list[@]}" | grep -E "^($SEL)" | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
for j in "${jobs_list[@]}"; do n=${j%%|*}; [ -f "$OUT/$n.log" ] && echo "$n: $(tail -1 "$OUT/$n.log" | cut -c1-170)"; done
