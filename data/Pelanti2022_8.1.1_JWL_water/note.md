# Pelanti2022_8.1.1_JWL_water

**Paper:** Pelanti, Int. J. Multiphase Flow 153 (2022) 104097, Sec. 8.1.1

**Figure(s):** Fig. 1

**Problem:** Detonation gases (JWL, x < 0, 2000 kg/m3, 4.6406e10 Pa) | liquid water (NASG, 1044 kg/m3, 1e5 Pa); t = 0.9 ms.

## Setting

| item | value |
|---|---|
| test case in the code | `case=61` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=61`) |
| cells | 500 |
| domain | [-10, 10] m |
| boundaries | transmissive |
| EOS | JWL (Gamma 0.25, rho0 1840, A 854.5e9, B 20.5e9, R1 4.6, R2 1.35); water NASG (1.187, 7028e5, q -1177788, cv 3610, b 6.61e-4) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `86e7d7a`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

Exact solution: `python3 tests/exact_mg.py 61 -10 10 N` (p* = 2.191e10 Pa, u* = 2397 m/s).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `JWL_water_t0.9ms.dat` | t = 0.9 ms | `case=61 N=500` | case 61 mode 1 order 2 rk 2 N 500 CFL 0.2 / steps 278 (rk1 4, mood 0) rejects 0 / stages 552 newton it 1872 (3.39/stage) maxres 9.41e-14 / closure fails 0 / dM1/M 4.54e-15 dM2/M 1.90e-15 dE/E -4.42e-15 / wall 46.03 s |

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
