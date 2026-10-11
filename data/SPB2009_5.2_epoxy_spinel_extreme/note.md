# SPB2009_5.2_epoxy_spinel_extreme

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 5.2

**Figure(s):** Fig. 20

**Problem:** Same mixture, p = 2e11 Pa (x < 0.6) | 1e5 Pa; t = 34 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=13` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=13`) |
| cells | 500 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | global maximum speed, `rus=0` |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `c4d60d3`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (global maximum speed, `rus=0`); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `epoxy_spinel_extreme_t34us.dat` | t = 34 us | `case=13 N=500` | case 13 mode 1 order 2 rk 2 N 500 CFL 0.2 / steps 364 (rk1 5, mood 0) rejects 0 / stages 723 newton it 2170 (3.00/stage) maxres 6.04e-14 / closure fails 0 / dM1/M 1.77e-08 dM2/M 1.77e-08 dE/E 3.56e-08 / wall 3.37 s |

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
