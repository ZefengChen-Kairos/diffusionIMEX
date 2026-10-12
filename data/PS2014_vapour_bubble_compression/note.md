# PS2014_vapour_bubble_compression

**Paper:** Pelanti & Shyue (2014), Sec. 5.4

**Figure(s):** Fig. 5

**Problem:** Dodecane vapour bubble (radius 0.2 m, centred) in liquid dodecane at 1e5 Pa and 600 K, compressed by a piston (left wall moving at 100 m/s); p-pT-pTg relaxation.

## Setting

| item | value |
|---|---|
| test case in the code | `case=72` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=72`) |
| cells | 200 x 200 |
| domain | [0, 1] x [0, 1] m |
| boundaries | left wall moving at +100 m/s (fixed grid: the piston face stays at x = 0), other sides walls |
| EOS | as PS2014_dodecane_shock_tube |
| relaxation | 2 = p-pT-pTg |
| material CFL | 0.2 (default) |
| Rusanov speed | global maximum speed, `rus=0` |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `86e7d7a`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (global maximum speed, `rus=0`); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

P-S use a moving boundary; here the grid is fixed and the liquid enters through the piston face (approximation, devlog V13).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `bubble_compression_t0.4ms.dat.gz` | t = 0.4 ms | `case=72 snap=0.4e-3,0.8e-3,1.2e-3` | case 72 mode 1 order 2 rk 2 N 200 CFL 0.2 / steps 760 (rk1 61, mood 17) rejects 4 / stages 1480 newton it 3968 (2.68/stage) maxres 9.96e-14 / closure fails 0 / dM1/M 1.29e-01 dM2/M 8.42e+00 dE/E 1.52e-01 / wall 425.15 s |
| `bubble_compression_t0.8ms.dat.gz` | t = 0.8 ms | `case=72 snap=0.4e-3,0.8e-3,1.2e-3` | case 72 mode 1 order 2 rk 2 N 200 CFL 0.2 / steps 760 (rk1 61, mood 17) rejects 4 / stages 1480 newton it 3968 (2.68/stage) maxres 9.96e-14 / closure fails 0 / dM1/M 1.29e-01 dM2/M 8.42e+00 dE/E 1.52e-01 / wall 425.15 s |
| `bubble_compression_t1.2ms.dat.gz` | t = 1.2 ms | `case=72 snap=0.4e-3,0.8e-3,1.2e-3` | case 72 mode 1 order 2 rk 2 N 200 CFL 0.2 / steps 760 (rk1 61, mood 17) rejects 4 / stages 1480 newton it 3968 (2.68/stage) maxres 9.96e-14 / closure fails 0 / dM1/M 1.29e-01 dM2/M 8.42e+00 dE/E 1.52e-01 / wall 425.15 s |

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
