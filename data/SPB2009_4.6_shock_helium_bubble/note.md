# SPB2009_4.6_shock_helium_bubble

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.6 (experiment of Layes & Le Metayer)

**Figure(s):** Figs. 13-14

**Problem:** Mach 1.5 shock in air (1.29 kg/m3, 101325 Pa) moving left onto a helium bubble (0.167 kg/m3, diameter 0.04 m, centre (0.2, 0.04) m); shock initially at x = 0.221 m (not given in the paper).

## Setting

| item | value |
|---|---|
| test case in the code | `case=70` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_PC_IMEX.cpp` or `case=70`) |
| cells | 600 x 160 |
| domain | [0, 0.3] x [0, 0.08] m |
| boundaries | left/right transmissive, bottom/top wall |
| EOS | air ideal gas (1.4); helium ideal gas (1.67) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_PC_IMEX.cpp` (final scheme: Zha-Bilgen split, volume fraction advanced entirely in the explicit step by the path-conservative Rusanov scheme; implicit pressure, pressure work, incremental fluxes and the three relaxations in one local closure; energy defect split by the shock relations (II.5); implicit Laplace dissipation of p and u), code commit `99bd20e`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped start-up steps first order; an inadmissible step is redone with half the time step (LSDIRK2 kept).
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit Laplace dissipation of pressure and velocity weighted by (1 - c tau/dx)+.
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

Sub-cell initialisation 4 x 4 in cells cut by the bubble.

## Files

| file | content | arguments | run statistics (steps, first-order start-up steps, rejects = halved steps, Newton, conservation, wall time) |
|---|---|---|---|
| `shock_helium_t134us.dat.gz` | t = 134 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3392 (rk1 8) rejects 0 floor 0 / stages 6776 newton it 14287 (2.11/stage) maxres 9.96e-14 Tmax 384 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 5224.11 s |
| `shock_helium_t274us.dat.gz` | t = 274 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3392 (rk1 8) rejects 0 floor 0 / stages 6776 newton it 14287 (2.11/stage) maxres 9.96e-14 Tmax 384 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 5224.11 s |
| `shock_helium_t344us.dat.gz` | t = 344 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3392 (rk1 8) rejects 0 floor 0 / stages 6776 newton it 14287 (2.11/stage) maxres 9.96e-14 Tmax 384 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 5224.11 s |
| `shock_helium_t554us.dat.gz` | t = 554 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3392 (rk1 8) rejects 0 floor 0 / stages 6776 newton it 14287 (2.11/stage) maxres 9.96e-14 Tmax 384 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 5224.11 s |

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
