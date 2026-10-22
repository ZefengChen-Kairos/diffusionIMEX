# Data of the IMEX results (final scheme, `src/SPB6_PC_IMEX.cpp`)

One folder per test; each has the data files and `note.md` (paper, figure, cells, CFL and all settings, command, run statistics).
Produced by `bash tests/run_all.sh run/pc`, `python3 tests/make_archive.py run/pc data` and `python3 tests/archive_extra.py run/pc data`; The 2D folders SPB2009_4.7_cavitating_Richtmyer_Meshkov, PS2014_vapour_bubble_compression and Pelanti2022_8.1.2_underwater_explosion still hold the V22 data (SPB6_TV_IMEX) until the 2D reruns of SPB6_PC_IMEX on the grids of Sec. 4 finish. Figures by `python3 tests/plot_paper_figs.py run/pc docs/figures` and `python3 tests/plot_2d_paper.py run/pc docs/figures`.

| test | paper | figure | cells | material CFL |
|---|---|---|---|---|
| [SPB2009_4.1_interface_advection](SPB2009_4.1_interface_advection/note.md) | Saurel, Petitpas & Berry 2009 | Fig. 5 | 200 | 0.2 (default) |
| [SPB2009_4.2_nitromethane](SPB2009_4.2_nitromethane/note.md) | Saurel, Petitpas & Berry 2009 | Figs. 6-7 | 500 | 0.2 (default) |
| [SPB2009_4.3.1_water_air](SPB2009_4.3.1_water_air/note.md) | Saurel, Petitpas & Berry 2009 | Figs. 8-9 | 1000 (Fig. 8) and 100 (Fig. 9) | 0.2 (default) |
| [SPB2009_4.3.2_water_air_TPa](SPB2009_4.3.2_water_air_TPa/note.md) | Saurel, Petitpas & Berry 2009 | Fig. 10 | 1000 | 0.2 (default) |
| [SPB2009_4.5_cavitation](SPB2009_4.5_cavitation/note.md) | Saurel, Petitpas & Berry 2009 | Fig. 12 | 1000 | 0.2 (default) |
| [SPB2009_5.1_epoxy_spinel](SPB2009_5.1_epoxy_spinel/note.md) | Saurel, Petitpas & Berry 2009 | Fig. 19 | 500 | 0.05 |
| [SPB2009_5.2_epoxy_spinel_extreme](SPB2009_5.2_epoxy_spinel_extreme/note.md) | Saurel, Petitpas & Berry 2009 | Fig. 20 | 500 | 0.2 (default) |
| [SPB2009_4.6_shock_helium_bubble](SPB2009_4.6_shock_helium_bubble/note.md) | Saurel, Petitpas & Berry 2009 | Figs. 13-14 | 600 x 160 | 0.2 (default) |
| [SPB2009_4.7_cavitating_Richtmyer_Meshkov](SPB2009_4.7_cavitating_Richtmyer_Meshkov/note.md) | Saurel, Petitpas & Berry 2009 | SPB Figs. 15-17 | 900 x 400 (grid of SPB 2009) | 0.2 (default) |
| [PS2014_water_cavitation_tube_2mps](PS2014_water_cavitation_tube_2mps/note.md) | Pelanti & Shyue 2014 | Fig. 1 | 5000 | 0.2 (default) |
| [PS2014_water_cavitation_tube_500mps](PS2014_water_cavitation_tube_500mps/note.md) | Pelanti & Shyue 2014 | Fig. 2 | 5000 | 0.2 (default) |
| [PS2014_dodecane_shock_tube](PS2014_dodecane_shock_tube/note.md) | Pelanti & Shyue 2014 | Fig. 3 | 2000 | 0.05 |
| [PS2014_vapour_bubble_compression](PS2014_vapour_bubble_compression/note.md) | Pelanti & Shyue 2014 | Fig. 5 | 200 x 200 | 0.2 (default) |
| [Pelanti2022_8.1.1_JWL_water](Pelanti2022_8.1.1_JWL_water/note.md) | Pelanti 2022 | Fig. 1 | 500 | 0.2 (default) |
| [Pelanti2022_8.1.2_underwater_explosion](Pelanti2022_8.1.2_underwater_explosion/note.md) | Pelanti 2022 | Figs. 2-3 | 721 x 400 (grid of Pelanti 2022) | 0.2 (default) |
| [Pelanti2022_8.2.1_water_cavitation](Pelanti2022_8.2.1_water_cavitation/note.md) | Pelanti 2022 | Figs. 4-5 | 5000 | 0.2 (default) |
| [Pelanti2022_8.2.2_dodecane_NASG](Pelanti2022_8.2.2_dodecane_NASG/note.md) | Pelanti 2022 | Fig. 6 | 2000 | 0.2 (default) |
| [Pelanti2022_8.2.4_Bartak](Pelanti2022_8.2.4_Bartak) | Pelanti 2022 | Fig. 8 (probe at x = 48 mm) | 1000 | 0.2 (default) |
| [Schropff2026_Riemann_problems](Schropff2026_Riemann_problems/note.md) | Schropff, Petitpas & Daniel 2026 | Secs. 5.1-5.2 | 200, 1000 (and 5000) | 0.2 (default) |
| [allmach_convergence](allmach_convergence) | Battisti & Boscheri 2025 | Secs. 4.1, 4.4 (tables; fields in `fields/`) | 32-256 | fixed dt / 0.2 |
