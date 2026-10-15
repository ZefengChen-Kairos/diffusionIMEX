# Pelanti2022_8.1.2_underwater_explosion

**Paper:** Pelanti (2022), Sec. 8.1.2

**Figure(s):** Figs. 2-3

**Problem:** Gas bubble (JWL state of 8.1.1, radius 0.05 m, centre (0, -0.32) m) in water below a rigid wall at y = 0.

## Setting

| item | value |
|---|---|
| test case in the code | `case=73` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=73`) |
| cells | 361 x 200 (half the paper grid 721 x 400: the IMEX run on the paper grid takes about 14 h; `N=721 ny=400` reproduces it) |
| domain | [-1.35, 1.35] x [-1.5, 0] m |
| boundaries | top (y = 0) wall, others transmissive |
| EOS | as Pelanti2022_8.1.1_JWL_water |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V16: energy-defect split from the shock relations (II.5), devlog V16), code commit `0c13222`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `UNDEX_t0.03ms.dat.gz` | t = 0.03 ms | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |
| `UNDEX_t0.06ms.dat.gz` | t = 0.06 ms | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |
| `UNDEX_t0.08ms.dat.gz` | t = 0.08 ms | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |
| `UNDEX_t0.21ms.dat.gz` | t = 0.21 ms | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |
| `UNDEX_t0.3ms.dat.gz` | t = 0.3 ms | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |
| `wall_pressure_history.dat` | Fig. 3: pressure at the wall centre (0, 0) vs time (columns: t [s], p [Pa]) | `case=73 snap=0.03e-3,0.06e-3,0.08e-3,0.21e-3,0.3e-3 probe=0,-1e-6` | case 73 mode 1 order 2 rk 2 N 361 CFL 0.2 / steps 536 (rk1 46, mood 0) rejects 0 floor 0 / stages 1026 newton it 3104 (3.03/stage) maxres 9.23e-14 Tmax 0 2.989e+04 / closure fails 0 / dM1/M -6.23e-07 dM2/M -1.21e-03 dE/E -2.55e-03 / wall 4753.58 s |

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
