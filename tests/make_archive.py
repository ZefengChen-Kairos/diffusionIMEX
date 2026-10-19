"""Write data/<test>/ from the raw IMEX runs of tests/run_all.sh: the profiles of every variable plotted in the paper figures
(.dat; 2D snapshots as .dat.gz with the plotted columns) and note.md with the paper, figure, grid, CFL and all other settings.
usage: python3 tests/make_archive.py RUNDIR DATADIR"""
import os, sys, re, gzip, shutil, subprocess
import numpy as np

RUN, DATA = sys.argv[1], sys.argv[2]
H = os.path.dirname(os.path.abspath(__file__))
COMMIT = subprocess.run(['git', '-C', H, 'log', '-1', '--format=%h', '--', '../src'], capture_output=True, text=True).stdout.strip()   # last commit of the source code
COLS1 = 'x y p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2'
COLS_DOC = """| column | meaning | unit |
|---|---|---|
| x, y | cell centre | m |
| p | mixture pressure $\\alpha_1p_1+\\alpha_2p_2$ | Pa |
| u, v | velocity | m/s |
| rho | mixture density | kg/m$^3$ |
| alpha2 | volume fraction of phase 2 (phase 1 = $1-\\alpha_2$) | – |
| Y2 | mass fraction of phase 2 | – |
| rho1, rho2 | phase densities | kg/m$^3$ |
| c_frozen, c_wood | frozen and Wood sound speeds | m/s |
| T1, T2 | phase temperatures | K |
| p1, p2 | phase pressures | Pa |
| g1-g2 | Gibbs free energy difference (SG/NASG only; NaN otherwise) | J/kg |"""
COLS2 = [0, 1, 2, 3, 4, 5, 6, 7, 12, 13]   # 2D: x y p u v rho alpha2 Y2 T1 T2
COLS2_NAMES = 'x y p u v rho alpha2 Y2 T1 T2'

SCHEME = """- Program: `src/SPB6_TV_IMEX.cpp` (TV-split IMEX, V22: energy-defect split from the shock relations (II.5) with the exact constant share of the shock profile, devlog V22; trace-phase-positive increment donors, compact shock-budget stencil and admissible face states, devlog V21-V22), code commit `{commit}`; run script `tests/run_all.sh`, this folder written by `tests/make_archive.py`.
- Time integration: LSDIRK2(2,2,2) (`rk=2`); start-up ramp from the acoustic CFL 0.4 (first step), at most x1.2 per step, ramped steps first order; an inadmissible step is redone once at the same step with backward Euler (MOOD), then halved.
- Space: second-order MUSCL, minmod (shared limiter for the volume-transported quantities, velocity separately), Rusanov with the material speed ({rus}); implicit pressure with the acoustic complement (B4.1).
- Material CFL: **{cfl}**. Newton tolerance 1e-13."""

# folder, paper, figure, description, case, settings dict, list of (raw name, archived file, label), extra text
T = [
 ('SPB2009_4.1_interface_advection', 'Saurel, Petitpas & Berry, J. Comput. Phys. 228 (2009) 1678-1712, Sec. 4.1', 'Fig. 5',
  'Advection of a water-air interface at u = 100 m/s (water x < 0.5, alpha_water = 1-1e-8; water 1000 kg/m3, air 10 kg/m3, p = 1e5 Pa); t = 2.79 ms.',
  4, dict(N='200', dom='[0, 1] m', bc='transmissive', eos='water SG (gamma 4.4, pi 6e8 Pa); air ideal gas (gamma 1.4)', relax='p (mechanical)'),
  [('spb41', 'interface_advection_t2.79ms.dat', 't = 2.79 ms')], ''),
 ('SPB2009_4.2_nitromethane', 'Saurel, Petitpas & Berry (2009), Sec. 4.2', 'Figs. 6-7',
  'Nitromethane shock tube with the Cochran-Chan EOS (same fluid as both phases): 20 GPa | 0.2 MPa, rho 1134 | 1200 kg/m3, interface at x = 0.5 m; t = 67 us.',
  60, dict(N='500', dom='[0, 1] m', bc='transmissive', eos='Cochran-Chan (Gamma 1.19, rho_ref 1134, A1 0.819181e9, A2 1.50835e9, E1 4.52969, E2 1.42144) for both phases', relax='p (mechanical)'),
  [('spb42', 'nitromethane_t67us.dat', 't = 67 us')], 'Exact solution: `python3 tests/exact_mg.py 60 0 1 N` (p* = 9.827 GPa, u* = 1770.6 m/s).'),
 ('SPB2009_4.3.1_water_air', 'Saurel, Petitpas & Berry (2009), Sec. 4.3.1', 'Figs. 8-9',
  'Water-air shock tube: water (x < 0.75, alpha = 1-1e-6, 1000 kg/m3, 1e9 Pa) | air (1 kg/m3, 1e5 Pa); t = 240 us. Fig. 8 uses 1000 cells, Fig. 9 100 cells.',
  5, dict(N='1000 (Fig. 8) and 100 (Fig. 9)', dom='[0, 1] m', bc='transmissive', eos='water SG (4.4, 6e8); air ideal gas (1.4)', relax='p'),
  [('spb431', 'water_air_N1000_t240us.dat', 'N = 1000'), ('spb431_n100', 'water_air_N100_t240us.dat', 'N = 100')], ''),
 ('SPB2009_4.3.2_water_air_TPa', 'Saurel, Petitpas & Berry (2009), Sec. 4.3.2', 'Fig. 10',
  'Water-air shock tube at 1 TPa: water (x < 0.75, 1000 kg/m3, 1e12 Pa) | air (10 kg/m3, 1e5 Pa); t = 8.3 us.',
  6, dict(N='1000', dom='[0, 1] m', bc='transmissive', eos='water SG (4.4, 6e8); air ideal gas (1.4)', relax='p'),
  [('spb432', 'water_air_TPa_t8.3us.dat', 't = 8.3 us')], ''),
 ('SPB2009_4.5_cavitation', 'Saurel, Petitpas & Berry (2009), Sec. 4.5', 'Fig. 12',
  'Expansion (cavitation) tube: liquid water with alpha_air = 1e-2, u = -100 | +100 m/s about x = 0.5, p = 1e5 Pa (air density 1 kg/m3, not given in the paper); t = 1.85 ms.',
  0, dict(N='1000', dom='[0, 1] m', bc='transmissive', eos='water SG (4.4, 6e8); air ideal gas (1.4)', relax='p'),
  [('spb45', 'cavitation_t1.85ms.dat', 't = 1.85 ms')], 'The 5-equation exact solution has a vacuum in [0.326, 0.674] m.'),
 ('SPB2009_5.1_epoxy_spinel', 'Saurel, Petitpas & Berry (2009), Sec. 5.1', 'Fig. 19',
  'Epoxy-spinel mixture shock tube: alpha_epoxy = 0.5954, rho 1185 | 3622 kg/m3 on both sides, p = 1e10 Pa (x < 0.6) | 1e5 Pa; t = 80 us.',
  12, dict(N='500', dom='[0, 1] m', bc='transmissive', eos='epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9)', relax='p'),
  [('spb51', 'epoxy_spinel_t80us.dat', 't = 80 us')], ''),
 ('SPB2009_5.2_epoxy_spinel_extreme', 'Saurel, Petitpas & Berry (2009), Sec. 5.2', 'Fig. 20',
  'Same mixture, p = 2e11 Pa (x < 0.6) | 1e5 Pa; t = 34 us.',
  13, dict(N='500', dom='[0, 1] m', bc='transmissive', eos='epoxy SG (2.43, 5.3e9); spinel SG (1.62, 141e9)', relax='p'),
  [('spb52', 'epoxy_spinel_extreme_t34us.dat', 't = 34 us')], ''),
 ('SPB2009_4.6_shock_helium_bubble', 'Saurel, Petitpas & Berry (2009), Sec. 4.6 (experiment of Layes & Le Metayer)', 'Figs. 13-14',
  'Mach 1.5 shock in air (1.29 kg/m3, 101325 Pa) moving left onto a helium bubble (0.167 kg/m3, diameter 0.04 m, centre (0.2, 0.04) m); shock initially at x = 0.221 m (not given in the paper).',
  70, dict(N='600 x 160', dom='[0, 0.3] x [0, 0.08] m', bc='left/right transmissive, bottom/top wall', eos='air ideal gas (1.4); helium ideal gas (1.67)', relax='p'),
  [('he_s0', 'shock_helium_t134us.dat.gz', 't = 134 us'), ('he_s1', 'shock_helium_t274us.dat.gz', 't = 274 us'),
   ('he_s2', 'shock_helium_t344us.dat.gz', 't = 344 us'), ('he_s3', 'shock_helium_t554us.dat.gz', 't = 554 us')], 'Sub-cell initialisation 4 x 4 in cells cut by the bubble.'),
 ('SPB2009_4.7_cavitating_Richtmyer_Meshkov', 'Saurel, Petitpas & Berry (2009), Sec. 4.7; grid of Pelanti & Shyue (2014), Sec. 5.3', 'SPB Figs. 15-17',
  'Water column (x < 1.2 + sqrt(0.36 - (y-0.5)^2)) and gas, everything moving at u = -200 m/s onto the left wall, p = 101325 Pa (water 1000, gas 100 kg/m3).',
  71, dict(N='900 x 400 (grid of SPB 2009)', dom='[0, 3] x [0, 1] m', bc='left wall, right transmissive, bottom/top wall', eos='water SG (4.4, 6e8); gas ideal (1.8, SPB; P-S use 1.4)', relax='p'),
  [('rm_s0', 'cavitating_RM_t2ms.dat.gz', 't = 2 ms'), ('rm_s1', 'cavitating_RM_t3.1ms.dat.gz', 't = 3.1 ms'), ('rm_s2', 'cavitating_RM_t4.8ms.dat.gz', 't = 4.8 ms'),
   ('rm_s3', 'cavitating_RM_t6.4ms.dat.gz', 't = 6.4 ms'), ('rm_s4', 'cavitating_RM_t8.6ms.dat.gz', 't = 8.6 ms')], ''),
 ('PS2014_water_cavitation_tube_2mps', 'Pelanti & Shyue, J. Comput. Phys. 259 (2014) 331-357, Sec. 5.1', 'Fig. 1',
  'Liquid water with alpha_vapour = 1e-2 at T = 354.728 K, p = 1e5 Pa, u = -2 | +2 m/s about x = 0.5; t = 3.2 ms. Relaxation levels p, p-pT, p-pT-pTg, p-pTg.',
  1, dict(N='5000', dom='[0, 1] m', bc='transmissive', eos='water liquid SG (2.35, 1e9, q -1167e3, cv 1816); vapour (1.43, 0, q 2030e3, q\' -23.4e3, cv 1040)', relax='0 = p, 1 = p-pT, 2 = p-pT-pTg, 3 = p-pTg (instantaneous)'),
  [('ps1_r0', 'p_relax.dat', 'p'), ('ps1_r1', 'p-pT_relax.dat', 'p-pT'), ('ps1_r2', 'p-pT-pTg_relax.dat', 'p-pT-pTg'), ('ps1_r3', 'p-pTg_relax.dat', 'p-pTg')], ''),
 ('PS2014_water_cavitation_tube_500mps', 'Pelanti & Shyue (2014), Sec. 5.1', 'Fig. 2',
  'As the 2 m/s case with u = -500 | +500 m/s; t = 0.58 ms.',
  2, dict(N='5000', dom='[0, 1] m', bc='transmissive', eos='as PS2014_water_cavitation_tube_2mps', relax='0 = p, 1 = p-pT, 2 = p-pT-pTg, 3 = p-pTg'),
  [('ps2_r0', 'p_relax.dat', 'p'), ('ps2_r1', 'p-pT_relax.dat', 'p-pT'), ('ps2_r2', 'p-pT-pTg_relax.dat', 'p-pT-pTg'), ('ps2_r3', 'p-pTg_relax.dat', 'p-pTg')], ''),
 ('PS2014_dodecane_shock_tube', 'Pelanti & Shyue (2014), Sec. 5.2', 'Fig. 3',
  'Dodecane liquid (x < 0.75, alpha_l = 1-1e-8, 500 kg/m3, 1e8 Pa) | vapour (2 kg/m3, 1e5 Pa); t = 473 us.',
  3, dict(N='2000', dom='[0, 1] m', bc='transmissive', eos='dodecane liquid SG (2.35, 4e8, q -775.269e3, cv 1077.7); vapour (1.025, 0, q -237.547e3, q\' -24.4e3, cv 1956.45)', relax='0 = p, 2 = p-pT-pTg, 3 = p-pTg'),
  [('ps3_r0', 'p_relax.dat', 'p'), ('ps3_r2', 'p-pT-pTg_relax.dat', 'p-pT-pTg'), ('ps3_r3', 'p-pTg_relax.dat', 'p-pTg')],
  'Material CFL 0.05: the 1e8 Pa rarefaction in the liquid is the structure to be resolved; at CFL 0.2 its acoustic Courant number is about 0.9 (devlog V14).'),
 ('PS2014_vapour_bubble_compression', 'Pelanti & Shyue (2014), Sec. 5.4', 'Fig. 5',
  'Dodecane vapour bubble (radius 0.2 m, centred) in liquid dodecane at 1e5 Pa and 600 K, compressed by a piston (left wall moving at 100 m/s); p-pT-pTg relaxation.',
  72, dict(N='200 x 200', dom='[0, 1] x [0, 1] m', bc='left wall moving at +100 m/s (fixed grid: the piston face stays at x = 0), other sides walls', eos='as PS2014_dodecane_shock_tube', relax='2 = p-pT-pTg'),
  [('bub_s0', 'bubble_compression_t0.4ms.dat.gz', 't = 0.4 ms'), ('bub_s1', 'bubble_compression_t0.8ms.dat.gz', 't = 0.8 ms'), ('bub_s2', 'bubble_compression_t1.2ms.dat.gz', 't = 1.2 ms')],
  'P-S use a moving boundary; here the grid is fixed and the liquid enters through the piston face (approximation, devlog V13).'),
 ('Pelanti2022_8.1.1_JWL_water', 'Pelanti, Int. J. Multiphase Flow 153 (2022) 104097, Sec. 8.1.1', 'Fig. 1',
  'Detonation gases (JWL, x < 0, 2000 kg/m3, 4.6406e10 Pa) | liquid water (NASG, 1044 kg/m3, 1e5 Pa); t = 0.9 ms.',
  61, dict(N='500', dom='[-10, 10] m', bc='transmissive', eos='JWL (Gamma 0.25, rho0 1840, A 854.5e9, B 20.5e9, R1 4.6, R2 1.35); water NASG (1.187, 7028e5, q -1177788, cv 3610, b 6.61e-4)', relax='p'),
  [('pel811', 'JWL_water_t0.9ms.dat', 't = 0.9 ms')], 'Exact solution: `python3 tests/exact_mg.py 61 -10 10 N` (p* = 2.191e10 Pa, u* = 2397 m/s).'),
 ('Pelanti2022_8.1.2_underwater_explosion', 'Pelanti (2022), Sec. 8.1.2', 'Figs. 2-3',
  'Gas bubble (JWL state of 8.1.1, radius 0.05 m, centre (0, -0.32) m) in water below a rigid wall at y = 0.',
  73, dict(N='721 x 400 (grid of Pelanti 2022)', dom='[-1.35, 1.35] x [-1.5, 0] m', bc='top (y = 0) wall, others transmissive', eos='as Pelanti2022_8.1.1_JWL_water', relax='p'),
  [('undex_s0', 'UNDEX_t0.03ms.dat.gz', 't = 0.03 ms'), ('undex_s1', 'UNDEX_t0.06ms.dat.gz', 't = 0.06 ms'), ('undex_s2', 'UNDEX_t0.08ms.dat.gz', 't = 0.08 ms'),
   ('undex_s3', 'UNDEX_t0.21ms.dat.gz', 't = 0.21 ms'), ('undex_s4', 'UNDEX_t0.3ms.dat.gz', 't = 0.3 ms'), ('undex.dat.probe', 'wall_pressure_history.dat', 'Fig. 3: pressure at the wall centre (0, 0) vs time (columns: t [s], p [Pa])')], ''),
 ('Pelanti2022_8.2.1_water_cavitation', 'Pelanti (2022), Sec. 8.2.1', 'Figs. 4-5',
  'Liquid water (NASG) with alpha_vapour = 1e-2 at 353 K, 1e5 Pa, u = -2 | +2 m/s about x = 0.5; t = 3 ms. Five relaxation levels.',
  30, dict(N='5000', dom='[0, 1] m', bc='transmissive', eos='water liquid NASG (1.187, 7028e5, q -1177788, cv 3610, b 6.61e-4); vapour (1.467, 0, q 2077616, q\' 14317, cv 955)',
           relax='p; pT(f) theta = 2000; pT; pTg(f) nu = 1e-4; pTg'),
  [('pel821_p', 'p_relax.dat', 'p'), ('pel821_pTf', 'pT_finite_theta2000.dat', 'pT(f), theta = 2000'), ('pel821_pT', 'pT_relax.dat', 'pT'),
   ('pel821_pTgf', 'pTg_finite_nu1e-4.dat', 'pTg(f), nu = 1e-4'), ('pel821_pTg', 'pTg_relax.dat', 'pTg')], ''),
 ('Pelanti2022_8.2.2_dodecane_NASG', 'Pelanti (2022), Sec. 8.2.2', 'Fig. 6',
  'Dodecane shock tube with NASG parameters (Table 4): liquid (x < 0.75, 500 kg/m3, 1e8 Pa) | vapour (2 kg/m3, 1e5 Pa); t = 473 us.',
  31, dict(N='2000', dom='[0, 1] m', bc='transmissive', eos='dodecane liquid NASG (1.206, 1681e5, q -996054, cv 2532, b 7.51e-4); vapour (1.021, 0, q -384592, q\' -4301, cv 2274)', relax='p; pTg'),
  [('pel822_r0', 'p_relax.dat', 'p'), ('pel822_r2', 'pTg_relax.dat', 'pTg')], ''),
]
CFL = {3: '0.05', 12: '0.05'}
RUS = {13: 'global maximum speed, `rus=0`', 71: 'global maximum speed, `rus=0`', 72: 'global maximum speed, `rus=0`'}
ARGS = {}   # command line per raw name, from tests/jobs_all.txt (1D) and tests/run_all.sh (2D)
for line in open(os.path.join(H, 'jobs_all.txt')):
    p = line.strip().split('|')
    if len(p) == 3 and p[0] == '.': ARGS[p[1]] = p[2]
for line in open(os.path.join(H, 'run_all.sh')):
    for m in re.finditer(r'"\.\|(\w+)\|([^"]+)"', line): ARGS[m.group(1)] = m.group(2)


def stat(name):
    f = os.path.join(RUN, name + '.log')
    if not os.path.exists(f): return None
    l = open(f).read().splitlines()
    return l[0] if l and l[0].startswith('case') else ('FAILED: ' + (l[-1] if l else 'empty log'))


os.makedirs(DATA, exist_ok=True)
index = []
for folder, paper, fig, desc, case, st, files, extra in T:
    d = os.path.join(DATA, folder); os.makedirs(d, exist_ok=True)
    rows, missing = [], []
    for raw, out, lab in files:
        src = os.path.join(RUN, raw if raw.endswith('.probe') else raw + '.dat')
        if not os.path.exists(src): missing.append(raw); continue
        dst = os.path.join(d, out)
        if out.endswith('.gz'):
            a = np.loadtxt(src)[:, COLS2]
            with gzip.open(dst, 'wt') as g: np.savetxt(g, a, fmt='%.8e', header=COLS2_NAMES + '   (snapshot, %s)' % lab)
        elif raw.endswith('.probe'):
            shutil.copyfile(src, dst)
        else:
            shutil.copyfile(src, dst)
        base = re.sub(r'_s\d+$', '', raw.replace('.dat.probe', ''))
        rows.append('| `%s` | %s | `%s` | %s |' % (out, lab, ARGS.get(base, ''), (stat(base) or '').replace('|', '/')))
    if missing: print(folder, 'missing', missing)
    cfl = CFL.get(case, '0.2 (default)') if case != 50 else 'fixed dt = 0.0025'
    note = ['# %s' % folder, '', '**Paper:** %s' % paper, '', '**Figure(s):** %s' % fig, '', '**Problem:** %s' % desc, '',
            '## Setting', '', '| item | value |', '|---|---|',
            '| test case in the code | `case=%d` (`Test` enum in `src/spb6_physics.h`; switch with `Test test = ...;` in `src/SPB6_TV_IMEX.cpp` or `case=%d`) |' % (case, case),
            '| cells | %s |' % st['N'], '| domain | %s |' % st['dom'], '| boundaries | %s |' % st['bc'], '| EOS | %s |' % st['eos'],
            '| relaxation | %s |' % st['relax'], '| material CFL | %s |' % cfl, '| Rusanov speed | %s |' % RUS.get(case, 'local per face (default)'), '',
            '## Scheme', '', SCHEME.format(commit=COMMIT, cfl=cfl, rus=RUS.get(case, 'local per face')), '']
    if extra: note += ['## Remarks', '', extra, '']
    note += ['## Files', '', '| file | content | arguments | run statistics (steps, first-order steps, MOOD steps, rejects, Newton, conservation, wall time) |', '|---|---|---|---|'] + rows + ['']
    if any(o.endswith('.gz') for _, o, _ in files):
        note += ['2D snapshots are gzipped text (`numpy.loadtxt` reads them directly), columns `%s`, 8 significant digits, one line per cell.' % COLS2_NAMES, '']
    if any(not o.endswith('.gz') and not r.endswith('.probe') for r, o, _ in files):
        note += ['1D/2D profile columns (`%s`, one line per cell, full precision):' % COLS1, '', COLS_DOC, '']
    open(os.path.join(d, 'note.md'), 'w').write('\n'.join(note))
    short = {'SPB2009': 'Saurel, Petitpas & Berry 2009', 'PS2014': 'Pelanti & Shyue 2014', 'Pelanti2022': 'Pelanti 2022'}[folder.split('_')[0]]
    index.append('| [%s](%s/note.md) | %s | %s | %s | %s |' % (folder, folder, short, fig, st['N'], cfl))
open(os.path.join(DATA, 'README.md'), 'w').write('\n'.join([
    '# Data of the IMEX results (V22)', '',
    'One folder per test; each has the data files and `note.md` (paper, figure, cells, CFL and all settings, command, run statistics).',
    'Produced by `bash tests/run_all.sh run/v22`, `python3 tests/make_archive.py run/v22 data` and `python3 tests/archive_extra.py run/v22 data`; figures by '
    '`python3 tests/plot_paper_figs.py run/v22 docs/figures` and `python3 tests/plot_2d_paper.py run/v22 docs/figures`.', '',
    '| test | paper | figure | cells | material CFL |', '|---|---|---|---|---|'] + index + [
    '| [Pelanti2022_8.2.4_Bartak](Pelanti2022_8.2.4_Bartak) | Pelanti 2022 | Fig. 8 (probe at x = 48 mm) | 1000 | 0.2 (default) |',
    '| [Schropff2026_Riemann_problems](Schropff2026_Riemann_problems/note.md) | Schropff, Petitpas & Daniel 2026 | Secs. 5.1-5.2 | 200, 1000 (and 5000) | 0.2 (default) |',
    '| [allmach_convergence](allmach_convergence) | Battisti & Boscheri 2025 | Secs. 4.1, 4.4 (tables; fields in `fields/`) | 32-256 | fixed dt / 0.2 |', '']))
print('wrote', DATA)
