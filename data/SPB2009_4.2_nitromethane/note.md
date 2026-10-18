# SPB2009_4.2_nitromethane

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.2

**Figure(s):** Figs. 6-7

**Problem:** Nitromethane shock tube with the Cochran-Chan EOS (same fluid as both phases): 20 GPa | 0.2 MPa, rho 1134 | 1200 kg/m3, interface at x = 0.5 m; t = 67 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=60` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=60`) |
| cells | 500 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | Cochran-Chan (Gamma 1.19, rho_ref 1134, A1 0.819181e9, A2 1.50835e9, E1 4.52969, E2 1.42144) for both phases |
| relaxation | p (mechanical) |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V21: energy-defect split from the shock relations (II.5), devlog V16; trace-phase-positive increment donors and compact shock-budget stencil, devlog V21), code commit `9aa22e2`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

Exact solution: `python3 tests/exact_mg.py 60 0 1 N` (p* = 9.827 GPa, u* = 1770.6 m/s).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `nitromethane_t67us.dat` | t = 67 us | `case=60 N=500` | case 60 mode 1 order 2 rk 2 N 500 CFL 0.2 / steps 300 (rk1 5, mood 0) rejects 0 floor 0 / stages 595 newton it 1788 (3.01/stage) maxres 5.50e-14 Tmax 0 0 / closure fails 0 / dM1/M 1.72e-08 dM2/M 1.16e-15 dE/E 5.35e-08 / wall 126.94 s |

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
