# Pelanti2022_8.1.2_underwater_explosion

**Paper:** Pelanti (2022), Sec. 8.1.2

**Figure(s):** Figs. 2-3

**Problem:** Gas bubble (JWL state of 8.1.1, radius 0.05 m, centre (0, -0.32) m) in water below a rigid wall at y = 0.

## Setting

| item | value |
|---|---|
| test case in the code | `case=73` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=73`) |
| cells | 721 x 400 |
| domain | [-1.35, 1.35] x [-1.5, 0] m |
| boundaries | top (y = 0) wall, others transmissive |
| EOS | as Pelanti2022_8.1.1_JWL_water |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `c4d60d3`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
