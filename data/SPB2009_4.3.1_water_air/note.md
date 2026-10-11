# SPB2009_4.3.1_water_air

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.3.1

**Figure(s):** Figs. 8-9

**Problem:** Water-air shock tube: water (x < 0.75, alpha = 1-1e-6, 1000 kg/m3, 1e9 Pa) | air (1 kg/m3, 1e5 Pa); t = 240 us. Fig. 8 uses 1000 cells, Fig. 9 100 cells.

## Setting

| item | value |
|---|---|
| test case in the code | `case=5` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=5`) |
| cells | 1000 (Fig. 8) and 100 (Fig. 9) |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water SG (4.4, 6e8); air ideal gas (1.4) |
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
| `water_air_N1000_t240us.dat` | N = 1000 | `case=5 N=1000` | case 5 mode 1 order 2 rk 2 N 1000 CFL 0.2 / steps 595 (rk1 6, mood 0) rejects 0 / stages 1184 newton it 2979 (2.52/stage) maxres 9.98e-14 / closure fails 0 / dM1/M -6.36e-15 dM2/M -1.82e-15 dE/E -1.14e-14 / wall 10.34 s |
| `water_air_N100_t240us.dat` | N = 100 | `case=5 N=100` | case 5 mode 1 order 2 rk 2 N 100 CFL 0.2 / steps 62 (rk1 6, mood 0) rejects 0 / stages 118 newton it 313 (2.65/stage) maxres 9.11e-14 / closure fails 0 / dM1/M 1.30e-05 dM2/M -9.93e-05 dE/E 2.51e-05 / wall 0.10 s |

1D/2D profile columns (`x y p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2`, one line per cell, full precision):

| column | meaning | unit |
|---|---|---|
| x, y | cell centre | m |
| p | mixture pressure $\alpha_1p_1+\alpha_2p_2$ | Pa |
| u, v | velocity | m/s |
| rho | mixture density | kg/m$^3$ |
| alpha2 | volume fraction of phase 2 (phase 1 = $1-\alpha_2$) | – |
| Y2 | mass fraction of phase 2 | – |
| rho1, rho2 | phase densities | kg/m$^3$ |
| c_frozen, c_wood | frozen and Wood sound speeds | m/s |
| T1, T2 | phase temperatures | K |
| p1, p2 | phase pressures | Pa |
| g1-g2 | Gibbs free energy difference (SG/NASG only; NaN otherwise) | J/kg |
