# SPB2009_4.5_cavitation

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.5

**Figure(s):** Fig. 12

**Problem:** Expansion (cavitation) tube: liquid water with alpha_air = 1e-2, u = -100 | +100 m/s about x = 0.5, p = 1e5 Pa (air density 1 kg/m3, not given in the paper); t = 1.85 ms.

## Setting

| item | value |
|---|---|
| test case in the code | `case=0` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=0`) |
| cells | 1000 |
| domain | [0, 1] m |
| boundaries | transmissive |
| EOS | water SG (4.4, 6e8); air ideal gas (1.4) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V21: energy-defect split from the shock relations (II.5), devlog V16; trace-phase-positive increment donors and compact shock-budget stencil, devlog V21), code commit `9aa22e2`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

The 5-equation exact solution has a vacuum in [0.326, 0.674] m.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `cavitation_t1.85ms.dat` | t = 1.85 ms | `case=0 N=1000` | case 0 mode 1 order 2 rk 2 N 1000 CFL 0.2 / steps 933 (rk1 12, mood 0) rejects 0 floor 0 / stages 1854 newton it 1889 (1.02/stage) maxres 9.54e-14 Tmax 177.5 250 / closure fails 0 / dM1/M -3.70e-01 dM2/M -3.70e-01 dE/E -3.70e-01 / wall 10.95 s |

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
