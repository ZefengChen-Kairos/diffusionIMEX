#!/bin/bash
# usage: suite.sh TAG   (options passed through the environment)
S=${RUN:-$(pwd)/run}; mkdir -p $S/su; T=$1; B=${SPB:-$(pwd)/spb}
r(){ out=$1; shift; "$@" > $S/su/${out}_$T.out 2>&1; }
CASE=5 r wa1000 $B 1 1000 0.2 $S/su/wa1000_$T.dat 5 2; CASE=5 r wa100 $B 1 100 0.2 $S/su/wa100_$T.dat 5 2; CASE=6 r w2 $B 1 1000 0.2 $S/su/w2_$T.dat 5 2
CASE=4 r adv $B 1 200 0.2 $S/su/adv_$T.dat 5 2; CASE=3 RELAX=0 r d0 $B 1 1000 0.2 $S/su/d0_$T.dat 5 2; CASE=3 RELAX=2 r d2 $B 1 1000 0.2 $S/su/d2_$T.dat 5 2
r cav $B 1 1000 0.2 $S/su/cav_$T.dat 5 2; CASE=1 RELAX=2 r c1 $B 1 1000 0.2 $S/su/c1_$T.dat 5 2; CASE=2 RELAX=2 r c2 $B 1 1000 0.2 $S/su/c2_$T.dat 5 2
CASE=12 r e1 $B 1 500 0.2 $S/su/e1_$T.dat 5 2; CASE=13 r e2 $B 1 500 0.2 $S/su/e2_$T.dat 5 2
for c in 0.2 0.1 0.05 0.003125; do CASE=21 $B 1 200 $c $S/su/sm${c}_$T.dat 20 2 >/dev/null 2>&1; done
