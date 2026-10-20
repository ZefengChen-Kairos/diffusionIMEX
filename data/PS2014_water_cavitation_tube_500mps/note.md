# PS2014_water_cavitation_tube_500mps

**Paper:** Pelanti & Shyue (2014), Sec. 5.1

**Figure(s):** Fig. 2

**Problem:** As the 2 m/s case with u = -500 | +500 m/s; t = 0.58 ms.

## Setting

| item | value |
|---|---|
| test case in the code | `case=2` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=2`) |
| cells | 5000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | as PS2014_water_cavitation_tube_2mps |
| relaxation | 0 = p, 1 = p-pT, 2 = p-pT-pTg, 3 = p-pTg |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_PC_IMEX.cpp` (final scheme: Zha-Bilgen split, volume fraction advanced entirely in the explicit step by the path-conservative Rusanov scheme; implicit pressure, pressure work, incremental fluxes and the three relaxations in one local closure; energy defect split by the shock relations (II.5); implicit Laplace dissipation of p and u), code commit `99bd20e`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped start-up steps first order; an inadmissible step is redone with half the time step (LSDIRK2 kept).
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit Laplace dissipation of pressure and velocity weighted by (1 - c tau/dx)+.
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order start-up steps, rejects = halved steps, Newton, conservation, wall time) |
|---|---|---|---|
| `p_relax.dat` | p | `case=2 N=5000 relax=0` | case 2 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 7252 (rk1 4) rejects 0 floor 0 / stages 14500 newton it 14505 (1.00/stage) maxres 8.66e-14 Tmax 388.7 354.7 / closure fails 0 / dM1/M -5.80e-01 dM2/M -5.80e-01 dE/E -5.80e-01 / wall 485.66 s |
| `p-pT_relax.dat` | p-pT | `case=2 N=5000 relax=1` | case 2 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 7252 (rk1 4) rejects 0 floor 0 / stages 14500 newton it 14504 (1.00/stage) maxres 8.26e-14 Tmax 388.7 388.7 / closure fails 0 / dM1/M -5.80e-01 dM2/M -5.80e-01 dE/E -5.80e-01 / wall 328.89 s |
| `p-pT-pTg_relax.dat` | p-pT-pTg | `case=2 N=5000 relax=2` | case 2 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 7252 (rk1 4) rejects 0 floor 0 / stages 14500 newton it 15990 (1.10/stage) maxres 9.99e-14 Tmax 386.5 386.5 / closure fails 0 / dM1/M -5.80e-01 dM2/M 8.92e+00 dE/E -5.80e-01 / wall 508.68 s |
| `p-pTg_relax.dat` | p-pTg | `case=2 N=5000 relax=3` | case 2 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 7252 (rk1 4) rejects 0 floor 0 / stages 14500 newton it 17813 (1.23/stage) maxres 1.00e-13 Tmax 386.6 386.5 / closure fails 0 / dM1/M -5.80e-01 dM2/M 8.93e+00 dE/E -5.80e-01 / wall 714.87 s |

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
