# SPB2009_5.2_epoxy_spinel_extreme

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 5.2

**Figure(s):** Fig. 20

**Problem:** Same mixture, p = 2e11 Pa (x < 0.6) | 1e5 Pa; t = 34 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=13` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=13`) |
| cells | 500 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | global maximum speed, `rus=0` |

## Scheme

- Program: `src/SPB6_PC_IMEX.cpp` (final scheme: Zha-Bilgen split, volume fraction advanced entirely in the explicit step by the path-conservative Rusanov scheme; implicit pressure, pressure work, incremental fluxes and the three relaxations in one local closure; energy defect split by the shock relations (II.5); implicit Laplace dissipation of p and u), code commit `99bd20e`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped start-up steps first order; an inadmissible step is redone with half the time step (LSDIRK2 kept).
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (global maximum speed, `rus=0`); implicit Laplace dissipation of pressure and velocity weighted by (1 - c tau/dx)+.
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order start-up steps, rejects = halved steps, Newton, conservation, wall time) |
|---|---|---|---|
| `epoxy_spinel_extreme_t34us.dat` | t = 34 us | `case=13 N=500` | case 13 mode 1 order 2 rk 2 N 500 CFL 0.2 / steps 356 (rk1 5) rejects 0 floor 0 / stages 707 newton it 2124 (3.00/stage) maxres 5.09e-14 Tmax 1.212e+05 1.518e+05 / closure fails 0 / dM1/M 1.72e-08 dM2/M 1.72e-08 dE/E 3.45e-08 / wall 4.88 s |

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

## Ablation of the energy-defect partition (`share_ablation/`)

Same run with the partition of the energy defect changed (command-line `share`, default 0), figure `docs/figures/spb20_share.png`
(`python3 tests/plot_share_ablation.py data/SPB2009_5.2_epoxy_spinel_extreme/share_ablation docs/figures/spb20_share.png`).

| file | partition | post-shock mixture density (plateau) | steps / rejects |
|---|---|---|---|
| `spb52_share0.dat` | shock relations: physical share Q_k/Q up to Q, error share lambda^E (the scheme; bitwise equal to the main run) | 3770.0 | 356 / 0 |
| `spb52_share1.dat` | no physical share: lambda^E for the whole defect | 4078.2 | 368 / 0 |
| `spb52_share2.dat` | mass fractions Y_k | 3960.0 | 364 / 0 |

Exact (5-equation, SPB (II.5)): 3766.7.
