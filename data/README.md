# Data of the IMEX reproductions of the paper figures (V14)

One folder per test; each has the data files and `note.md` (paper, figure, cells, CFL and all settings, command, run statistics).
Produced by `bash tests/run_archive.sh run/archive` and `python3 tests/make_archive.py run/archive data`; figures by `python3 tests/plot_paper_figs.py run/archive docs/figures` and `python3 tests/plot_2d_paper.py run/archive docs/figures`.

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
| [SPB2009_4.7_cavitating_Richtmyer_Meshkov](SPB2009_4.7_cavitating_Richtmyer_Meshkov/note.md) | Saurel, Petitpas & Berry 2009 | SPB Figs. 15-17 | 300 x 100 (P-S 2014 grid; SPB used 900 x 400: N=900 ny=400) | 0.2 (default) |
| [PS2014_water_cavitation_tube_2mps](PS2014_water_cavitation_tube_2mps/note.md) | Pelanti & Shyue 2014 | Fig. 1 | 5000 | 0.2 (default) |
| [PS2014_water_cavitation_tube_500mps](PS2014_water_cavitation_tube_500mps/note.md) | Pelanti & Shyue 2014 | Fig. 2 | 5000 | 0.2 (default) |
| [PS2014_dodecane_shock_tube](PS2014_dodecane_shock_tube/note.md) | Pelanti & Shyue 2014 | Fig. 3 | 2000 | 0.05 |
| [PS2014_vapour_bubble_compression](PS2014_vapour_bubble_compression/note.md) | Pelanti & Shyue 2014 | Fig. 5 | 200 x 200 | 0.2 (default) |
| [Pelanti2022_8.1.1_JWL_water](Pelanti2022_8.1.1_JWL_water/note.md) | Pelanti 2022 | Fig. 1 | 500 | 0.2 (default) |
| [Pelanti2022_8.1.2_underwater_explosion](Pelanti2022_8.1.2_underwater_explosion/note.md) | Pelanti 2022 | Figs. 2-3 | 361 x 200 (half the paper grid 721 x 400: the IMEX run on the paper grid takes about 14 h; `N=721 ny=400` reproduces it) | 0.2 (default) |
| [Pelanti2022_8.2.1_water_cavitation](Pelanti2022_8.2.1_water_cavitation/note.md) | Pelanti 2022 | Figs. 4-5 | 5000 | 0.2 (default) |
| [Pelanti2022_8.2.2_dodecane_NASG](Pelanti2022_8.2.2_dodecane_NASG/note.md) | Pelanti 2022 | Fig. 6 | 2000 | 0.2 (default) |
