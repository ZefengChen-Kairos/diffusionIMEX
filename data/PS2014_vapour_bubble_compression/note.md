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

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `c4d60d3`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (global maximum speed, `rus=0`); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

P-S use a moving boundary; here the grid is fixed and the liquid enters through the piston face (approximation, devlog V13).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
