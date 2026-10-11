# SPB2009_4.3.2_water_air_TPa

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.3.2

**Figure(s):** Fig. 10

**Problem:** Water-air shock tube at 1 TPa: water (x < 0.75, 1000 kg/m3, 1e12 Pa) | air (10 kg/m3, 1e5 Pa); t = 8.3 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=6` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=6`) |
| cells | 1000 |
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
| `water_air_TPa_t8.3us.dat` | t = 8.3 us | `case=6 N=1000` | case 6 mode 1 order 2 rk 2 N 1000 CFL 0.2 / steps 1391 (rk1 34, mood 0) rejects 7 / stages 2755 newton it 7600 (2.76/stage) maxres 7.53e-14 / closure fails 0 / dM1/M -7.36e-03 dM2/M -9.99e-01 dE/E -2.85e-02 / wall 27.16 s |

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
