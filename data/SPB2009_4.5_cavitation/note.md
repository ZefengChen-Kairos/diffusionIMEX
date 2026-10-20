# SPB2009_4.5_cavitation

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.5

**Figure(s):** Fig. 12

**Problem:** Expansion (cavitation) tube: liquid water with alpha_air = 1e-2, u = -100 | +100 m/s about x = 0.5, p = 1e5 Pa (air density 1 kg/m3, not given in the paper); t = 1.85 ms.

## Setting

| item | value |
|---|---|
| test case in the code | `case=0` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=0`) |
| cells | 1000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water SG (4.4, 6e8); air ideal gas (1.4) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_PC_IMEX.cpp` (final scheme: Zha-Bilgen split, volume fraction advanced entirely in the explicit step by the path-conservative Rusanov scheme; implicit pressure, pressure work, incremental fluxes and the three relaxations in one local closure; energy defect split by the shock relations (II.5); implicit Laplace dissipation of p and u), code commit `99bd20e`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped start-up steps first order; an inadmissible step is redone with half the time step (LSDIRK2 kept).
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit Laplace dissipation of pressure and velocity weighted by (1 - c tau/dx)+.
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

The 5-equation exact solution has a vacuum in [0.326, 0.674] m.

## Files

| file | content | arguments | run statistics (steps, first-order start-up steps, rejects = halved steps, Newton, conservation, wall time) |
|---|---|---|---|
| `cavitation_t1.85ms.dat` | t = 1.85 ms | `case=0 N=1000` | case 0 mode 1 order 2 rk 2 N 1000 CFL 0.2 / steps 933 (rk1 12) rejects 0 floor 0 / stages 1854 newton it 1889 (1.02/stage) maxres 9.56e-14 Tmax 177.6 250 / closure fails 0 / dM1/M -3.70e-01 dM2/M -3.70e-01 dE/E -3.70e-01 / wall 11.89 s |

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
