# Pelanti2022_8.2.1_water_cavitation

**Paper:** Pelanti (2022), Sec. 8.2.1

**Figure(s):** Figs. 4-5

**Problem:** Liquid water (NASG) with alpha_vapour = 1e-2 at 353 K, 1e5 Pa, u = -2 | +2 m/s about x = 0.5; t = 3 ms. Five relaxation levels.

## Setting

| item | value |
|---|---|
| test case in the code | `case=30` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=30`) |
| cells | 5000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water liquid NASG (1.187, 7028e5, q -1177788, cv 3610, b 6.61e-4); vapour (1.467, 0, q 2077616, q' 14317, cv 955) |
| relaxation | p; pT(f) theta = 2000; pT; pTg(f) nu = 1e-4; pTg |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `622bda6`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `p_relax.dat` | p | `case=30 N=5000 relax=0` | case 30 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 178 (rk1 33, mood 0) rejects 0 / stages 323 newton it 867 (2.68/stage) maxres 9.97e-14 / closure fails 0 / dM1/M -1.20e-02 dM2/M -1.20e-02 dE/E -1.20e-02 / wall 10.57 s |
| `pT_finite_theta2000.dat` | pT(f), theta = 2000 | `case=30 N=5000 relax=1 theta=2000` | case 30 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 178 (rk1 33, mood 0) rejects 0 / stages 323 newton it 865 (2.68/stage) maxres 9.89e-14 / closure fails 0 / dM1/M -1.20e-02 dM2/M -1.20e-02 dE/E -1.20e-02 / wall 14.69 s |
| `pT_relax.dat` | pT | `case=30 N=5000 relax=1` | case 30 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 178 (rk1 33, mood 0) rejects 0 / stages 323 newton it 666 (2.06/stage) maxres 9.95e-14 / closure fails 0 / dM1/M -1.20e-02 dM2/M -1.20e-02 dE/E -1.20e-02 / wall 10.93 s |
| `pTg_finite_nu1e-4.dat` | pTg(f), nu = 1e-4 | `case=30 N=5000 relax=2 nu=1e-4` | case 30 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 178 (rk1 33, mood 0) rejects 0 / stages 323 newton it 666 (2.06/stage) maxres 9.94e-14 / closure fails 0 / dM1/M -1.20e-02 dM2/M 2.36e-01 dE/E -1.20e-02 / wall 15.82 s |
| `pTg_relax.dat` | pTg | `case=30 N=5000 relax=2` | case 30 mode 1 order 2 rk 2 N 5000 CFL 0.2 / steps 178 (rk1 33, mood 0) rejects 0 / stages 323 newton it 782 (2.42/stage) maxres 9.73e-14 / closure fails 0 / dM1/M -1.20e-02 dM2/M 3.53e-01 dE/E -1.20e-02 / wall 15.63 s |

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
