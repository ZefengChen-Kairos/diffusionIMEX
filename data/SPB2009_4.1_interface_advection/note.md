# SPB2009_4.1_interface_advection

**Paper:** Saurel, Petitpas & Berry, J. Comput. Phys. 228 (2009) 1678-1712, Sec. 4.1

**Figure(s):** Fig. 5

**Problem:** Advection of a water-air interface at u = 100 m/s (water x < 0.5, alpha_water = 1-1e-8; water 1000 kg/m3, air 10 kg/m3, p = 1e5 Pa); t = 2.79 ms.

## Setting

| item | value |
|---|---|
| test case in the code | `case=4` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=4`) |
| cells | 200 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water SG (gamma 4.4, pi 6e8 Pa); air ideal gas (gamma 1.4) |
| relaxation | p (mechanical) |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `86e7d7a`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `interface_advection_t2.79ms.dat` | t = 2.79 ms | `case=4 N=200` | case 4 mode 1 order 2 rk 2 N 200 CFL 0.2 / steps 287 (rk1 12, mood 0) rejects 0 / stages 562 newton it 0 (0.00/stage) maxres 2.28e-14 / closure fails 0 / dM1/M 5.58e-01 dM2/M -5.58e-01 dE/E 5.58e-01 / wall 0.23 s |

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
