# SPB2009_4.7_cavitating_Richtmyer_Meshkov

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.7; grid of Pelanti & Shyue (2014), Sec. 5.3

**Figure(s):** SPB Figs. 15-17

**Problem:** Water column (x < 1.2 + sqrt(0.36 - (y-0.5)^2)) and gas, everything moving at u = -200 m/s onto the left wall, p = 101325 Pa (water 1000, gas 100 kg/m3).

## Setting

| item | value |
|---|---|
| test case in the code | `case=71` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=71`) |
| cells | 300 x 100 (P-S 2014 grid; SPB used 900 x 400: N=900 ny=400) |
| domain | [0, 3] x [0, 1] m |
| boundaries | left wall, right transmissive, bottom/top wall |
| EOS | water SG (4.4, 6e8); gas ideal (1.8, SPB; P-S use 1.4) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | global maximum speed, `rus=0` |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V16: energy-defect split from the shock relations (II.5), devlog V16), code commit `0c13222`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (global maximum speed, `rus=0`); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `cavitating_RM_t2ms.dat.gz` | t = 2 ms | `case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3` | FAILED: case 71 mode 1 order 2 rk 2 N 300 CFL 0.2 / steps 2202 (rk1 60, mood 23) rejects 0 floor 178 / stages 4366 newton it 15014 (3.44/stage) maxres 4.82e-13 Tmax 270.5 1.31e+07 / closure fails 0 / dM1/M -3.72e-02 dM2/M -7.87e-03 dE/E -4.86e-02 / wall 1568.30 s |
| `cavitating_RM_t3.1ms.dat.gz` | t = 3.1 ms | `case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3` | FAILED: case 71 mode 1 order 2 rk 2 N 300 CFL 0.2 / steps 2202 (rk1 60, mood 23) rejects 0 floor 178 / stages 4366 newton it 15014 (3.44/stage) maxres 4.82e-13 Tmax 270.5 1.31e+07 / closure fails 0 / dM1/M -3.72e-02 dM2/M -7.87e-03 dE/E -4.86e-02 / wall 1568.30 s |
| `cavitating_RM_t4.8ms.dat.gz` | t = 4.8 ms | `case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3` | FAILED: case 71 mode 1 order 2 rk 2 N 300 CFL 0.2 / steps 2202 (rk1 60, mood 23) rejects 0 floor 178 / stages 4366 newton it 15014 (3.44/stage) maxres 4.82e-13 Tmax 270.5 1.31e+07 / closure fails 0 / dM1/M -3.72e-02 dM2/M -7.87e-03 dE/E -4.86e-02 / wall 1568.30 s |
| `cavitating_RM_t6.4ms.dat.gz` | t = 6.4 ms | `case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3` | FAILED: case 71 mode 1 order 2 rk 2 N 300 CFL 0.2 / steps 2202 (rk1 60, mood 23) rejects 0 floor 178 / stages 4366 newton it 15014 (3.44/stage) maxres 4.82e-13 Tmax 270.5 1.31e+07 / closure fails 0 / dM1/M -3.72e-02 dM2/M -7.87e-03 dE/E -4.86e-02 / wall 1568.30 s |
| `cavitating_RM_t8.6ms.dat.gz` | t = 8.6 ms | `case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3` | FAILED: case 71 mode 1 order 2 rk 2 N 300 CFL 0.2 / steps 2202 (rk1 60, mood 23) rejects 0 floor 178 / stages 4366 newton it 15014 (3.44/stage) maxres 4.82e-13 Tmax 270.5 1.31e+07 / closure fails 0 / dM1/M -3.72e-02 dM2/M -7.87e-03 dE/E -4.86e-02 / wall 1568.30 s |

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
