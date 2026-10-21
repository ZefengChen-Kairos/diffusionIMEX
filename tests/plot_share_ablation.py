"""SPB 2009 Fig. 20 (epoxy-spinel, 2e11 Pa | 1e5 Pa, t = 34 us, 500 cells): partition of the energy defect, ablation.
share=0: shock relations (physical share Q_k/Q, error share lambda^E); share=1: no physical share (lambda^E only);
share=2: mass fractions Y_k.   usage: python3 tests/plot_share_ablation.py RUNDIR OUT.png
RUNDIR holds spb52_share0.dat, spb52_share1.dat, spb52_share2.dat (src/SPB6_PC_IMEX.cpp case=13 share=0|1|2)"""
import os, sys, runpy
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
RUN, OUT = sys.argv[1], sys.argv[2]
argv = sys.argv; sys.argv = ['x', RUN, os.path.dirname(os.path.abspath(OUT)), 'none']
g = runpy.run_path(os.path.join(H, 'plot_paper_figs.py')); sys.argv = argv
P = g['PH']['spb_es']
ex = g['exact5'](P, (0.5954, 1185., 3622., 0., 2e11), (0.5954, 1185., 3622., 0., 1e5), 0.6, 34e-6)
def load(m):
    d = np.loadtxt(os.path.join(RUN, 'spb52_share%d.dat' % m))
    s = {'x': d[:, 0], 'p': d[:, 2], 'u': d[:, 3], 'rho': d[:, 5], 'a2': d[:, 6], 'r1': d[:, 8], 'r2': d[:, 9]}
    s['a1'] = 1 - s['a2']; return s
runs = [(0, 'shock relations (this work)', dict(color='k', ls='none', marker='o', mfc='none', ms=3.5, mew=0.8, markevery=0.012)),
        (1, r'no physical share ($\lambda^{\mathcal{E}}_k$ only)', dict(color='tab:red', ls='--', lw=1.2)),
        (2, 'mass fractions $Y_k$', dict(color='tab:blue', ls='-.', lw=1.2))]
var = [('rho', 'mixture density (kg/m$^3$)'), ('a1', 'volume fraction of epoxy'),
       ('r1', 'density of epoxy (kg/m$^3$)'), ('r2', 'density of spinel (kg/m$^3$)')]
plt.rcParams.update({'font.size': 9, 'legend.fontsize': 7.5})
fig, axs = plt.subplots(2, 2, figsize=(9.2, 6.0))
for ax, (v, lab) in zip(axs.flat, var):
    ax.plot(ex['x'], ex[v], 'k-', lw=1.0, label='exact (5-eq.)')
    for m, l, kw in runs:
        s = load(m); ax.plot(s['x'], s[v], label=l, **kw)
    ax.set_xlim(0.55, 1.0); ax.set_xlabel('x (m)'); ax.set_ylabel(lab); ax.grid(alpha=0.25)
axs[0, 0].legend(loc='upper left')
fig.suptitle('SPB 2009 Fig. 20: epoxy-spinel shock tube, 2e11 Pa | 1e5 Pa, t = 34 us, 500 cells: partition of the energy defect', fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(OUT, dpi=130); print('wrote', OUT)
