# Pelanti2022_8.2.2_dodecane_NASG

**Paper:** Pelanti (2022), Sec. 8.2.2

**Figure(s):** Fig. 6

**Problem:** Dodecane shock tube with NASG parameters (Table 4): liquid (x < 0.75, 500 kg/m3, 1e8 Pa) | vapour (2 kg/m3, 1e5 Pa); t = 473 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=31` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=31`) |
| cells | 2000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | dodecane liquid NASG (1.206, 1681e5, q -996054, cv 2532, b 7.51e-4); vapour (1.021, 0, q -384592, q' -4301, cv 2274) |
| relaxation | p; pTg |
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
| `p_relax.dat` | p | `case=31 N=2000 relax=0` | case 31 mode 1 order 2 rk 2 N 2000 CFL 0.2 / steps 1144 (rk1 5) rejects 0 floor 0 / stages 2283 newton it 4575 (2.00/stage) maxres 7.71e-14 Tmax 642 1.047e+06 / closure fails 0 / dM1/M -7.76e-15 dM2/M 3.30e-15 dE/E -3.77e-14 / wall 49.09 s |
| `pTg_relax.dat` | pTg | `case=31 N=2000 relax=2` | case 31 mode 1 order 2 rk 2 N 2000 CFL 0.2 / steps 1326 (rk1 5) rejects 0 floor 0 / stages 2647 newton it 5305 (2.00/stage) maxres 8.83e-14 Tmax 642 1.047e+06 / closure fails 0 / dM1/M -5.40e-04 dM2/M 4.05e-01 dE/E -3.98e-14 / wall 58.14 s |

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
