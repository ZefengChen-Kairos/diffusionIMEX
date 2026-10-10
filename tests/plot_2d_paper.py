"""Figures of the 2D paper tests from tests/run_2d.sh: rows = snapshot times, columns = explicit | IMEX for each variable.
usage: python3 tests/plot_2d_paper.py RUNDIR OUTDIR [he rm bub undex]"""
import os, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from plot_2d import field, LAB, CMAP

RUN, OUT = sys.argv[1], sys.argv[2]; ONLY = sys.argv[3:]
os.makedirs(OUT, exist_ok=True)
SETUP = 'explicit = HLLC + relaxation, Heun, CFL 0.4;  IMEX = TV-split, LSDIRK2, material CFL 0.2'
FIGS = {   # name: (title, times, variables, figure file)
    'he': ('SPB 2009 Figs. 13-14: Mach 1.5 shock in air / helium bubble, 600 x 160', ['134', '274', '344', '554'], ['schlieren', 'a1'], 'spb13'),
    'rm': ('SPB 2009 Figs. 15-17: cavitating Richtmyer-Meshkov, 300 x 100 (P-S grid)', ['2', '3.1', '4.8', '6.4', '8.6'], ['rho', 'a1'], 'spb16'),
    'bub': ('P-S 2014 Fig. 5: dodecane vapour bubble compressed by a piston (100 m/s), relax 2, 200 x 200', ['0.4', '0.8', '1.2'], ['Y2', 'pbar'], 'ps5'),
    'undex': ('Pelanti 2022 Fig. 2: underwater explosion near a wall, 721 x 400', ['0.03', '0.06', '0.08', '0.21', '0.3'], ['logp', 'a1'], 'pel2'),
}
UNIT = {'he': 'us', 'rm': 'ms', 'bub': 'ms', 'undex': 'ms'}

for name, (title, times, vars_, fig_name) in FIGS.items():
    if ONLY and name not in ONLY: continue
    files = [(k, os.path.join(RUN, '%s_%s_s%d.dat' % (name, K, i))) for i in range(len(times)) for k, K in (('explicit', 'E'), ('IMEX', 'I'))]
    if not all(os.path.exists(f) for _, f in files): print('missing snapshots for', name); continue
    xs, ys, _ = field(files[0][1]); asp = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    nc = 2 * len(vars_); w = 3.3
    fig, axs = plt.subplots(len(times), nc, figsize=(w * nc, w * asp * len(times) + 1.0), squeeze=False)
    for r, tt in enumerate(times):
        for kk, (lab, K) in enumerate((('explicit', 'E'), ('IMEX', 'I'))):
            xs, ys, F = field(os.path.join(RUN, '%s_%s_s%d.dat' % (name, K, r)))
            for c, v in enumerate(vars_):
                ax = axs[r, 2 * c + kk]
                im = ax.pcolormesh(xs, ys, F[v], shading='auto', cmap=CMAP.get(v, 'turbo'))
                ax.set_aspect('equal'); ax.tick_params(labelsize=6)
                if v != 'schlieren': fig.colorbar(im, ax=ax, fraction=0.046 * max(asp, 0.3), pad=0.02).ax.tick_params(labelsize=6)
                ax.set_title('%s, t = %s %s: %s' % (lab, tt, UNIT[name], LAB.get(v, v)), fontsize=7)
    fig.suptitle(title + '\n' + SETUP, fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(OUT, fig_name + '.png'), dpi=100); plt.close(fig); print('wrote', fig_name)

if (not ONLY or 'undex' in ONLY):   # Pelanti Fig. 3: pressure history at the wall centre (0, 0)
    pr = [(lab, os.path.join(RUN, 'undex_%s.dat.probe' % K)) for lab, K in (('explicit', 'E'), ('IMEX', 'I'))]
    if all(os.path.exists(f) for _, f in pr):
        fig, ax = plt.subplots(figsize=(6, 3.4))
        for (lab, f), kw in zip(pr, (dict(color='tab:blue'), dict(color='tab:red', ls='--'))):
            d = np.loadtxt(f); ax.plot(d[:, 0] * 1e3, d[:, 1] / 1e9, label=lab, **kw)
        ax.set_xlabel('t (ms)'); ax.set_ylabel('p at (0, 0) (GPa)'); ax.grid(alpha=0.3); ax.legend()
        ax.set_title('Pelanti 2022 Fig. 3: wall pressure at the centre (0, 0)\n' + SETUP, fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, 'pel3.png'), dpi=110); plt.close(fig); print('wrote pel3')
