#!/bin/bash
export OMP_NUM_THREADS=2
./ref_mhd 256 1 0.4 > ref256.log 2>&1
./ref_mhd 512 1 0.4 > ref512.log 2>&1
./ref_mhd 1024 1 0.4 > ref1024.log 2>&1
./ref_mhd 512 0 0.4 > ref512_sharp.log 2>&1
echo done > ref.done
