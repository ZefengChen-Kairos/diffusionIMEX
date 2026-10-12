# PS2014_water_cavitation_tube_2mps

**Paper:** Pelanti & Shyue, J. Comput. Phys. 259 (2014) 331-357, Sec. 5.1

**Figure(s):** Fig. 1

**Problem:** Liquid water with alpha_vapour = 1e-2 at T = 354.728 K, p = 1e5 Pa, u = -2 | +2 m/s about x = 0.5; t = 3.2 ms. Relaxation levels p, p-pT, p-pT-pTg, p-pTg.

## Setting

| item | value |
|---|---|
| test case in the code | `case=1` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=1`) |
| cells | 5000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water liquid SG (2.35, 1e9, q -1167e3, cv 1816); vapour (1.43, 0, q 2030e3, q' -23.4e3, cv 1040) |
| relaxation | 0 = p, 1 = p-pT, 2 = p-pT-pTg, 3 = p-pTg (instantaneous) |
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
| `p_relax.dat` | p | `case=1 N=5000 relax=0` | case 1 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 188 (rk1 33, mood 0) rejects 0 / stages 343 newton it 719 (2.10/stage) maxres 9.91e-14 / closure fails 0 / dM1/M -1.28e-02 dM2/M -1.28e-02 dE/E -1.28e-02 / wall 14.69 s |
| `p-pT_relax.dat` | p-pT | `case=1 N=5000 relax=1` | case 1 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 188 (rk1 33, mood 0) rejects 0 / stages 343 newton it 695 (2.03/stage) maxres 6.44e-14 / closure fails 0 / dM1/M -1.28e-02 dM2/M -1.28e-02 dE/E -1.28e-02 / wall 15.71 s |
| `p-pT-pTg_relax.dat` | p-pT-pTg | `case=1 N=5000 relax=2` | case 1 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 188 (rk1 33, mood 0) rejects 0 / stages 343 newton it 785 (2.29/stage) maxres 9.86e-14 / closure fails 0 / dM1/M -1.28e-02 dM2/M 4.30e-01 dE/E -1.28e-02 / wall 18.71 s |
| `p-pTg_relax.dat` | p-pTg | `case=1 N=5000 relax=3` | case 1 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 188 (rk1 33, mood 0) rejects 0 / stages 343 newton it 893 (2.60/stage) maxres 9.99e-14 / closure fails 0 / dM1/M -1.28e-02 dM2/M 4.00e-01 dE/E -1.28e-02 / wall 22.61 s |

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
