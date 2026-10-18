# SPB2009_4.6_shock_helium_bubble

**Paper:** Saurel, Petitpas & Berry (2009), Sec. 4.6 (experiment of Layes & Le Metayer)

**Figure(s):** Figs. 13-14

**Problem:** Mach 1.5 shock in air (1.29 kg/m3, 101325 Pa) moving left onto a helium bubble (0.167 kg/m3, diameter 0.04 m, centre (0.2, 0.04) m); shock initially at x = 0.221 m (not given in the paper).

## Setting

| item | value |
|---|---|
| test case in the code | `case=70` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=70`) |
| cells | 600 x 160 |
| domain | [0, 0.3] x [0, 0.08] m |
| boundaries | left/right transmissive, bottom/top wall |
| EOS | air ideal gas (1.4); helium ideal gas (1.67) |
| relaxation | p |
| material CFL | 0.2 (default) |
| Rusanov speed | local per face (default) |

## Scheme

- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V21: energy-defect split from the shock relations (II.5), devlog V16; trace-phase-positive increment donors and compact shock-budget stencil, devlog V21), code commit `9aa22e2`; run script `tests/run_archive.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed (local per face); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **0.2 (default)**. Newton tolerance 1e-13.

## Remarks

Sub-cell initialisation 4 x 4 in cells cut by the bubble.

## Files

| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |
|---|---|---|---|
| `shock_helium_t134us.dat.gz` | t = 134 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3402 (rk1 40, mood 0) rejects 0 floor 0 / stages 6764 newton it 14297 (2.11/stage) maxres 9.97e-14 Tmax 298.6 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 2716.04 s |
| `shock_helium_t274us.dat.gz` | t = 274 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3402 (rk1 40, mood 0) rejects 0 floor 0 / stages 6764 newton it 14297 (2.11/stage) maxres 9.97e-14 Tmax 298.6 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 2716.04 s |
| `shock_helium_t344us.dat.gz` | t = 344 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3402 (rk1 40, mood 0) rejects 0 floor 0 / stages 6764 newton it 14297 (2.11/stage) maxres 9.97e-14 Tmax 298.6 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 2716.04 s |
| `shock_helium_t554us.dat.gz` | t = 554 us | `case=70 snap=134e-6,274e-6,344e-6,554e-6` | case 70 mode 1 order 2 rk 2 N 600 CFL 0.2 / steps 3402 (rk1 40, mood 0) rejects 0 floor 0 / stages 6764 newton it 14297 (2.11/stage) maxres 9.97e-14 Tmax 298.6 2331 / closure fails 0 / dM1/M 5.27e-01 dM2/M 4.23e-06 dE/E 8.56e-01 / wall 2716.04 s |

2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `x y p u v rho alpha2 Y2 T1 T2`, 8 significant digits, one line per cell.
