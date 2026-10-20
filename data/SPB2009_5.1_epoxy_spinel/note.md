# SPB2009_5.1_epoxy_spinel

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 5.1

**Figure(s):** Fig. 19

**Problem:** Epoxy-spinel mixture shock tube: alpha_epoxy = 0.5954, rho 1185 | 3622 kg/m3 on both sides, p = 1e10 Pa (x < 0.6) | 1e5 Pa; t = 80 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=12` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=12`) |
| cells | 500 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9) |
| relaxation | p |
| material CFL | 0.05 |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_PC_IMEX.cpp` (final scheme: Zha-Bilgen split, volume fraction advanced entirely in the explicit step by the path-conservative Rusanov scheme; implicit pressure, pressure work, incremental fluxes and the three relaxations in one local closure; energy defect split by the shock relations (II.5); implicit Laplace dissipation of p and u), code commit `99bd20e`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped start-up steps first order; an inadmissible step is redone with half the time step (LSDIRK2 kept).
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit Laplace dissipation of pressure and velocity weighted by (1 - c tau/dx)+.
- Material CFL: **0.05**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order start-up steps, rejects = halved steps, Newton, conservation, wall time) |
|---|---|---|---|
| `epoxy_spinel_t80us.dat` | t = 80 us | `case=12 N=500` | case 12 mode 1 order 2 rk 2 N 500 CFL 0.05 / steps 443 (rk1 5) rejects 0 floor 0 / stages 881 newton it 2209 (2.51/stage) maxres 9.81e-14 Tmax 9029 6.724e+04 / closure fails 0 / dM1/M 9.90e-16 dM2/M -8.42e-15 dE/E 1.55e-15 / wall 4.83 s |

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
