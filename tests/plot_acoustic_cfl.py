"""Schropff et al. 2026 pure-acoustic Riemann problems (cases 81, 85): material CFL 0.2 (default) against acoustic CFL 0.4
(acdt=1 cfl=0.4). Density error of Eq. (51) for every grid and one figure per case with the 1000-cell solutions.
usage: python3 tests/plot_acoustic_cfl.py MATDIR ACDIR OUTDIR
MATDIR holds cC_N.dat (run/pc/sch), ACDIR holds cC_N_ac.dat."""
import os, sys, glob, runpy, tempfile
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
MAT, AC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
argv = sys.argv; sys.argv = ['x', tempfile.mkdtemp(), tempfile.mkdtemp()]
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()): g = runpy.run_path(os.path.join(H, 'plot_schropff.py'))
sys.argv = argv
CASES, exact = g['CASES'], g['exact']

def load(f):
    if not os.path.exists(f): return None
    d = np.loadtxt(f); return {'x': d[:, 0], 'p': d[:, 2], 'u': d[:, 3], 'rho': d[:, 5]}

def steps(f):
    s = open(f.replace('.dat', '.log')).read(); import re
    m = re.search(r'steps (\d+) \(rk1 (\d+)', s); return (int(m.group(1)), int(m.group(2))) if m else (0, 0)

print('| case | cells | material CFL 0.2: error (steps, first-order steps) | acoustic CFL 0.4: error (steps) |\n|---|---|---|---|')
plt.rcParams.update({'font.size': 9, 'legend.fontsize': 7.5})
for c in (81, 85):
    title, _, _, L, R, x0, t, _ = CASES[c]
    for n in (200, 1000, 5000):
        row = []
        for f in (os.path.join(MAT, 'c%d_%d.dat' % (c, n)), os.path.join(AC, 'c%d_%d_ac.dat' % (c, n))):
            s = load(f)
            if s is None: row.append('—'); continue
            dx = 1.0 / n; xc = (np.arange(n) + 0.5) * dx
            e = np.sum(np.abs(s['rho'] - exact(c, xc, t)['rho'])) / np.sum(np.abs(exact(c, xc, 0)['rho']))
            st = steps(f); row.append('%.3e (%d, %d)' % (e, *st))
        print('| %d | %d | %s | %s |' % (c, n, row[0], row[1]))
    xe = np.linspace(0, 1, 4001); ex = exact(c, xe, t)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    for k, (v, lab) in enumerate([('rho', r'$\rho$ (kg/m$^3$)'), ('u', '$u$ (m/s)'), ('p', '$p$ (Pa)')]):
        ax[k].plot(xe, ex[v], 'k-', lw=1.0, label='exact')
        m = load(os.path.join(MAT, 'c%d_1000.dat' % c)); a = load(os.path.join(AC, 'c%d_1000_ac.dat' % c))
        ax[k].plot(m['x'], m[v], color='tab:red', ls='--', lw=1.1, label='material CFL 0.2')
        ax[k].plot(a['x'], a[v], color='k', ls='none', marker='o', mfc='none', ms=3.5, mew=0.8, markevery=0.012, label='acoustic CFL 0.4')
        ax[k].set_xlabel('$x$ (m)'); ax[k].set_ylabel(lab); ax[k].grid(alpha=0.25)
        if v == 'rho' or v == 'p':
            d = np.abs(ex[v] - ex[v][0]); j = d > 0.02 * d.max()
            if j.any(): lo, hi = xe[j].min(), xe[j].max(); w = hi - lo; ax[k].set_xlim(max(0, lo - 0.6 * w), min(1, hi + 0.6 * w))
    ax[0].legend(loc='best')
    fig.suptitle('%s (Schropff et al. 2026), t = %g s, 1000 cells: material against acoustic CFL' % (title, t), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(OUT, 'sch%d_cfl.png' % c), dpi=130); plt.close(fig)
