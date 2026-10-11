# Gresho_vortex_low_Mach

**Paper:** Gresho vortex as in the user's 2DEuler_TV.cpp (testGreshoVortex); asymptotic-preserving check (devlog V11.8)

**Figure(s):** —

**Problem:** Single-phase limit (both phases ideal gas 1.4, alpha = 0.5, rho = 1), maximum velocity 1 at r = 0.2, p0 = 1/(gamma M^2); one turn t = 0.4 pi; the same fixed dt = 0.0025 for every M.

## Setting

| item | value |
|---|---|
| test case in the code | `case=50` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=50`) |
| cells | 80 x 80 |
| domain | [-0.5, 0.5]^2 |
| boundaries | periodic |
| EOS | ideal gas 1.4 (both phases) |
| relaxation | p |
| material CFL | fixed dt = 0.0025 |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V14), code commit `c4d60d3`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **fixed dt = 0.0025**. Newton tolerance 1e-13.

## Remarks

Cases 50-55, fixed time step `dt=0.0025` (no CFL, no start-up ramp).

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `M1e-1.dat` | M = 1e-1 | `case=50 dt=0.0025` | case 50 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 2012 (2.00/stage) maxres 3.19e-16 / closure fails 0 / dM1/M -2.84e-15 dM2/M -2.84e-15 dE/E 2.96e-14 / wall 58.20 s |
| `M1e-2.dat` | M = 1e-2 | `case=51 dt=0.0025` | case 51 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 1039 (1.03/stage) maxres 9.85e-14 / closure fails 0 / dM1/M 4.12e-15 dM2/M 4.12e-15 dE/E 4.73e-14 / wall 64.11 s |
| `M1e-3.dat` | M = 1e-3 | `case=52 dt=0.0025` | case 52 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 1011 (1.00/stage) maxres 6.67e-14 / closure fails 0 / dM1/M -4.97e-15 dM2/M -4.97e-15 dE/E 3.39e-14 / wall 80.97 s |
| `M1e-4.dat` | M = 1e-4 | `case=53 dt=0.0025` | case 53 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 1011 (1.00/stage) maxres 5.67e-14 / closure fails 0 / dM1/M -7.11e-15 dM2/M -7.11e-15 dE/E -3.48e-14 / wall 81.28 s |
| `M1e-5.dat` | M = 1e-5 | `case=54 dt=0.0025` | case 54 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 1011 (1.00/stage) maxres 7.60e-14 / closure fails 0 / dM1/M 5.37e-14 dM2/M 5.37e-14 dE/E 1.16e-14 / wall 81.79 s |
| `M1e-6.dat` | M = 1e-6 | `case=55 dt=0.0025` | case 55 mode 1 order 2 rk 2 N 80 CFL 0.2 / steps 503 (rk1 0, mood 0) rejects 0 / stages 1006 newton it 1014 (1.01/stage) maxres 8.37e-14 / closure fails 0 / dM1/M 0.00e+00 dM2/M 0.00e+00 dE/E 5.25e-16 / wall 79.57 s |

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
