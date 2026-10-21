# Schropff2026_Riemann_problems

**Paper:** Schropff, Petitpas & Daniel, J. Comput. Phys. 547 (2026) 114545, Secs. 5.1 (Table 2), 5.1.1 and 5.2 (Table 3).

**Problem:** nine 1D Riemann problems on [0, 1] m with transmissive boundaries (cases 80-88 in `src/spb6_physics.h`). Single-phase cases are two identical phases in the same state (alpha1 = 1 - 1e-8); two-phase cases use the p relaxation (5-equation limit). EOS: air ideal gas 1.4; water SG (4.4, 6e8); epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9).

| case | test | left (rho, u, p) or (alpha1, rho1, rho2, u, p) | right | x0 | t_end |
|---|---|---|---|---|---|
| 80 | low-Mach air advection | (1000, 1, 1e5) | (0.01, 1, 1e5) | 0.25 | 0.5 s |
| 81 | low-Mach water shock tube | (1000, 0, 1e8) | (1000, 15, 0.98e8) | 0.5 | 1e-4 s |
| 82 | Lax shock tube (air) | (0.445, 0.698, 3.528) | (0.5, 0, 0.571) | 0.5 | 0.14 |
| 83 | strong air shock tube | (1, 0, 1e10) | (0.125, 0, 0.1) | 0.5 | 2.5e-6 s |
| 84 | two rarefactions (air) | (1, -1, 0.4) | (1, 1, 0.4) | 0.5 | 0.15 |
| 85 | liquid-liquid shock tube (water, Sec. 5.1.1) | (1000, 0, 1e6) | (1000, 0, 1e5) | 0.5 | 0.15 ms |
| 86 | water(1)-air(2) advection | (1-1e-7, 1000, 0.01, 1, 1e5) | (1e-7, 1000, 0.01, 1, 1e5) | 0.25 | 0.25 s |
| 87 | water-air shock tube | (1-1e-7, 1000, 50, 0, 1e9) | (1e-7, 1000, 50, 0, 1e5) | 0.7 | 241 us |
| 88 | epoxy-spinel shock tube | (0.5954, 1185, 3622, 0, 1e10) | (0.5954, 1185, 3622, 0, 1e5) | 0.6 | 30 us |

## Setting

Program `src/SPB6_PC_IMEX.cpp` (final scheme: ZB split, path-conservative Rusanov), LSDIRK2, material CFL 0.2, minmod, local Rusanov speed, Newton tolerance 1e-13 (the defaults; no per-test setting). Command: `src/bin_SPB6_PC_IMEX case=C N=N out=cC_N.dat`. Figures `docs/figures/sch80.png` ... `sch88.png` by `python3 tests/plot_schropff.py run/pc/sch docs/figures` (also prints the L1 table, Eq. (51) of the paper).

## Run statistics (steps, first-order ramp steps, rejects, Newton iterations per stage, wall time in s with 5 runs in parallel)

| file | steps | first order | rejects | Newton/stage | wall (s) |
|---|---|---|---|---|---|
| `c80_1000.dat` | 2537 | 42 | 0 | 0.00 | 28.40 |
| `c80_200.dat` | 537 | 42 | 0 | 0.00 | 1.01 |
| `c81_1000.dat` | 25 | 23 | 0 | 2.81 | 0.50 |
| `c81_200.dat` | 17 | 17 | 0 | 2.88 | 0.09 |
| `c82_1000.dat` | 1082 | 4 | 0 | 2.51 | 40.21 |
| `c82_200.dat` | 218 | 4 | 0 | 2.55 | 1.58 |
| `c83_1000.dat` | 1503 | 37 | 8 | 4.09 | 72.57 |
| `c83_200.dat` | 317 | 37 | 8 | 4.53 | 4.11 |
| `c84_1000.dat` | 750 | 0 | 0 | 2.00 | 29.45 |
| `c84_200.dat` | 150 | 0 | 0 | 2.02 | 0.98 |
| `c85_1000.dat` | 27 | 27 | 0 | 2.00 | 0.28 |
| `c85_200.dat` | 18 | 18 | 0 | 2.00 | 0.05 |
| `c86_1000.dat` | 1287 | 42 | 0 | 0.00 | 12.97 |
| `c86_200.dat` | 287 | 42 | 0 | 0.00 | 0.36 |
| `c87_1000.dat` | 588 | 6 | 0 | 3.00 | 23.55 |
| `c87_200.dat` | 120 | 6 | 0 | 3.00 | 0.72 |
| `c87_5000.dat` | 2923 | 6 | 0 | 3.00 | 434.73 |
| `c88_1000.dat` | 98 | 11 | 0 | 3.02 | 3.39 |
| `c88_200.dat` | 24 | 11 | 0 | 3.05 | 0.17 |

The 8 rejects of case 83 all occur in the first 23 steps (t < 2e-8 s, pressure ratio 1e11), on both meshes.

## Pure-acoustic problems with the acoustic CFL condition (`acoustic_cfl/`)

In cases 81 (low-Mach water shock tube) and 85 (liquid-liquid shock tube) |u| is so small that the material CFL step exceeds
t_end, so every step is a start-up first-order step and the error hardly decreases with the grid. The same cases with the
acoustic CFL number 0.4 (`acdt=1 cfl=0.4`); figures `docs/figures/sch81_cfl.png`, `sch85_cfl.png`
(`python3 tests/plot_acoustic_cfl.py run/pc/sch data/Schropff2026_Riemann_problems/acoustic_cfl docs/figures`).
Density error of Eq. (51) (steps, first-order steps):

| case | cells | material CFL 0.2 | acoustic CFL 0.4 |
|---|---|---|---|
| 81 | 200 | 3.738e-4 (17, 17) | 1.083e-4 (89, 0) |
| 81 | 1000 | 2.686e-4 (25, 23) | 3.509e-5 (442, 0) |
| 81 | 5000 | - | 1.101e-5 (2210, 0) |
| 85 | 200 | 1.995e-5 (18, 18) | 4.723e-6 (122, 0) |
| 85 | 1000 | 1.903e-5 (27, 27) | 1.646e-6 (610, 0) |
| 85 | 5000 | - | 5.681e-7 (3050, 0) |

Orders with the acoustic CFL: case 81 0.70, 0.72; case 85 0.65, 0.66 (shock tubes, first-order-limited at the discontinuities).
