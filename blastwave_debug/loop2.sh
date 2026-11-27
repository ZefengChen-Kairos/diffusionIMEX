#!/bin/bash
set -u
D=/tmp/claude-0/-home-user-diffusionIMEX/fa41aa46-3215-5a47-a5d7-66870c4d191b/scratchpad/outtest
P=/home/user/diffusionIMEX/blastwave_debug/patch
cd $D
say(){ echo "[$(date +%T)] $*" | tee -a loop2.log; }
if [ ! -x mhdloop2 ]; then
  say "构建 128x64 field loop (A0=1)"
  sed -e 's/^Test test = testblastwave;/Test test = testFieldLoop;/' \
      -e 's/^const int nxCells = 1000;/const int nxCells = 128;/' \
      -e 's/^const int nyCells = 1000;/const int nyCells = 64;/' \
      -e 's/^const double LOOP_A0 = 1e-3;/const double LOOP_A0 = 1.0;/' $P/Final_emb_omp_patched.cpp > tloop2.cpp
  g++ -O2 -march=native -fopenmp -std=c++17 -I/usr/include/eigen3 tloop2.cpp -o mhdloop2 2>&1|tail -2
fi
mkdir -p e2
for spec in "base:" "noad:hllemoff=1" "epsB0:epsBmax=0"; do
  t=${spec%%:*}; opt=${spec#*:}
  if [ -s e2/$t/DONE ]; then say "$t 已完成, 跳过"; continue; fi
  say "$t 开始 ($opt)"; rm -rf e2/$t; mkdir -p e2/$t; cd e2/$t
  OMP_NUM_THREADS=4 $D/mhdloop2 cfl=0.2 senv=2 tag=$t $opt >run.log 2>/dev/null && touch DONE
  cd $D; say "$t 完成 steps=$(grep -c '^step ' e2/$t/run.log)"
done
say "全部完成"
