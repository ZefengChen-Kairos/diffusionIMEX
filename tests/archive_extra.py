"""Archive the results that tests/make_archive.py does not cover, from the runs of tests/run_all.sh:
Bartak probes (Pelanti 2022 8.2.4), Schropff et al. 2026 Riemann problems, Taylor-Green (BB 4.4) and two-phase vortex (BB 4.1)
tables and fields.  Figures go to docs/figures.  2D fields are stored as .dat.gz with the columns x y p u v rho alpha2.
usage: python3 tests/archive_extra.py RUNDIR DATADIR"""
import os, sys, glob, gzip, shutil, subprocess
import numpy as np
RUN, DATA = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
H = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(H, '..', 'docs', 'figures')
py = lambda *a, **k: subprocess.run([sys.executable, os.path.join(H, a[0])] + list(a[1:]), check=True, **k)

# Bartak
py('plot_bartak.py', os.path.join(RUN, 'bartak'), os.path.join(DATA, 'Pelanti2022_8.2.4_Bartak'), os.path.join(FIG, 'bartak.png'))

# Schropff et al. 2026
d = os.path.join(DATA, 'Schropff2026_Riemann_problems'); os.makedirs(d, exist_ok=True)
for f in sorted(glob.glob(os.path.join(RUN, 'sch', 'c*.dat'))): shutil.copy(f, d)
out = py('plot_schropff.py', os.path.join(RUN, 'sch'), FIG, capture_output=True, text=True).stdout
open(os.path.join(d, 'L1_table.md'), 'w').write('Relative L1 density error, Eq. (51) of Schropff et al. (2026), from `python3 tests/plot_schropff.py`:\n\n```\n' + out + '```\n')

# all-Mach: tables
d = os.path.join(DATA, 'allmach_convergence'); os.makedirs(d, exist_ok=True)
open(os.path.join(d, 'taylor_green.md'), 'w').write(py('allmach_table.py', os.path.join(RUN, 'allmach'), capture_output=True, text=True).stdout)
vroot = os.path.join(RUN, 'vortex_table_root'); os.makedirs(vroot, exist_ok=True)
if not os.path.exists(os.path.join(vroot, 'vortex')): os.symlink(os.path.join(RUN, 'vortex'), os.path.join(vroot, 'vortex'))
open(os.path.join(d, 'vortex_p0_1_to_1e4.md'), 'w').write(py('allmach_table.py', vroot, capture_output=True, text=True).stdout)
py('vortex_bb_table.py', os.path.join(FIG, 'vortex_table'), os.path.join(RUN, 'vortex'))
for f in glob.glob(os.path.join(FIG, 'vortex_table_*.md')) + glob.glob(os.path.join(FIG, 'vortex_table_*.tex')): shutil.copy(f, d)

# all-Mach: fields
for sub in ['allmach', 'vortex']:
    od = os.path.join(d, 'fields', 'taylor_green' if sub == 'allmach' else 'vortex'); os.makedirs(od, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(RUN, sub, '*.dat'))):
        a = np.loadtxt(f)[:, :7]
        with gzip.open(os.path.join(od, os.path.basename(f) + '.gz'), 'wt') as g:
            np.savetxt(g, a, fmt='%.10e', header='x y p u v rho alpha2   (%s)' % next((l for l in open(f[:-4] + '.log').read().splitlines() if l.startswith('case')), '')[:300])
print('archive_extra done')
