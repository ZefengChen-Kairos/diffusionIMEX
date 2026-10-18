# PS2014_dodecane_shock_tube

**Paper:** Pelanti & Shyue (2014), Sec. 5.2

**Figure(s):** Fig. 3

**Problem:** Dodecane liquid (x < 0.75, alpha_l = 1-1e-8, 500 kg/m3, 1e8 Pa) | vapour (2 kg/m3, 1e5 Pa); t = 473 us.

## Setting

| item | value |
|---|---|
| test case in the code | `case=3` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=3`) |
| cells | 2000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | dodecane liquid SG (2.35, 4e8, q -775.269e3, cv 1077.7); vapour (1.025, 0, q -237.547e3, q' -24.4e3, cv 1956.45) |
| relaxation | 0 = p, 2 = p-pT-pTg, 3 = p-pTg |
| material CFL | 0.05 |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V21: energy-defect split from the shock relations (II.5), devlog V16; trace-phase-positive increment donors and compact shock-budget stencil, devlog V21), code commit `9aa22e2`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.05**. Newton tolerance 1e-13.

## Remarks

Material CFL 0.05: the 1e8 Pa rarefaction in the liquid is the structure to be resolved; at CFL 0.2 its acoustic Courant number is about 0.9 (devlog V14).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `p_relax.dat` | p | `case=3 N=2000 relax=0` | case 3 mode 1 order 2 rk 2 N 2000 CFL 0.05 / steps 2675 (rk1 3, mood 0) rejects 0 floor 0 / stages 5347 newton it 10717 (2.00/stage) maxres 7.08e-14 Tmax 687.3 1.022e+06 / closure fails 0 / dM1/M 4.12e-08 dM2/M 2.61e-15 dE/E 5.18e-08 / wall 91.03 s |
| `p-pT-pTg_relax.dat` | p-pT-pTg | `case=3 N=2000 relax=2` | case 3 mode 1 order 2 rk 2 N 2000 CFL 0.05 / steps 5896 (rk1 3, mood 0) rejects 0 floor 0 / stages 11789 newton it 23587 (2.00/stage) maxres 9.79e-14 Tmax 687.3 1.022e+06 / closure fails 0 / dM1/M -4.22e-03 dM2/M 3.16e+00 dE/E 6.47e-08 / wall 212.42 s |
| `p-pTg_relax.dat` | p-pTg | `case=3 N=2000 relax=3` | case 3 mode 1 order 2 rk 2 N 2000 CFL 0.05 / steps 5896 (rk1 3, mood 0) rejects 0 floor 0 / stages 11789 newton it 23587 (2.00/stage) maxres 6.01e-14 Tmax 687.3 1.022e+06 / closure fails 0 / dM1/M -4.22e-03 dM2/M 3.16e+00 dE/E 6.47e-08 / wall 214.71 s |

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
