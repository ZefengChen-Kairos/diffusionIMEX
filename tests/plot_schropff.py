"""Schropff, Petitpas & Daniel (2026, JCP 547, 114545) Secs. 5.1, 5.1.1, 5.2 with the IMEX program (cases 80-88):
density, velocity, pressure and Mach number against the exact solution (tests/exact5.py: 5-equation model; single-phase
cases are two identical phases), and the relative L1 density error of their Eq. (51).
usage: python3 tests/plot_schropff.py RUNDIR OUTDIR   (RUNDIR holds c<case>_<N>.dat)"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import exact5 as e5

RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, '..', 'run', 'sch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H, '..', 'docs', 'figures')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 7, 'lines.linewidth': 1.1})
AIR, WATER, EPOXY, SPINEL = (1.4, 0.0), (4.4, 6e8), (2.43, 5.3e9), (1.62, 141e9)
E8, E7 = 1e-8, 1e-7
# case: (title, EOS phase 1, EOS phase 2, L (alpha1, rho1, rho2, u, p), R, x_inter, t_end, grids)
CASES = {
    80: ('Low-Mach air advection', AIR, AIR, (1 - E8, 1000, 1000, 1, 1e5), (1 - E8, 0.01, 0.01, 1, 1e5), 0.25, 0.5, (200, 1000)),
    81: ('Low-Mach water shock tube', WATER, WATER, (1 - E8, 1000, 1000, 0, 1e8), (1 - E8, 1000, 1000, 15, 0.98e8), 0.5, 1e-4, (200, 1000)),
    82: ('Lax shock tube (air)', AIR, AIR, (1 - E8, .445, .445, .698, 3.528), (1 - E8, .5, .5, 0, .571), 0.5, 0.14, (200, 1000)),
    83: ('Strong air shock tube', AIR, AIR, (1 - E8, 1, 1, 0, 1e10), (1 - E8, .125, .125, 0, .1), 0.5, 2.5e-6, (200, 1000)),
    84: ('Two rarefaction fans (air)', AIR, AIR, (1 - E8, 1, 1, -1, .4), (1 - E8, 1, 1, 1, .4), 0.5, 0.15, (200, 1000)),
    85: ('Liquid-liquid shock tube (water)', WATER, WATER, (1 - E8, 1000, 1000, 0, 1e6), (1 - E8, 1000, 1000, 0, 1e5), 0.5, 0.15e-3, (200, 1000)),
    86: ('Water-air advection', WATER, AIR, (1 - E7, 1000, .01, 1, 1e5), (E7, 1000, .01, 1, 1e5), 0.25, 0.25, (200, 1000)),
    87: ('Water-air shock tube', WATER, AIR, (1 - E7, 1000, 50, 0, 1e9), (E7, 1000, 50, 0, 1e5), 0.7, 241e-6, (200, 1000, 5000)),
    88: ('Epoxy-spinel shock tube', EPOXY, SPINEL, (.5954, 1185, 3622, 0, 1e10), (.5954, 1185, 3622, 0, 1e5), 0.6, 30e-6, (200, 1000)),
}


def exact(c, x, t):
    _, P1, P2, L, R, x0, _, _ = CASES[c]; G = [P1[0], P2[0]]; PI = [P1[1], P2[1]]
    side = lambda a, r1, r2, u, p: e5.Side([a, 1 - a], [r1, r2], u, p, G, PI)
    if t == 0: return {'rho': np.where(x < x0, side(*L).rho, side(*R).rho)}
    o, _, _ = e5.solve(side(*L), side(*R), x0, t, x)
    a = [o['a1'], 1 - o['a1']]; r = [o['r1'], o['r2']]
    with np.errstate(invalid='ignore', divide='ignore'):
        c2 = [G[k] * (o['p'] + PI[k]) / r[k] for k in (0, 1)]
        o['cw'] = np.sqrt(1 / (o['rho'] * (a[0] / (r[0] * c2[0]) + a[1] / (r[1] * c2[1]))))
    return o


def load(c, n, tag=''):
    f = os.path.join(RUN, 'c%d_%d%s.dat' % (c, n, tag))
    if not os.path.exists(f): return None
    d = np.loadtxt(f); return {'x': d[:, 0], 'p': d[:, 2], 'u': d[:, 3], 'rho': d[:, 5], 'cw': d[:, 11]}


rows = []
for c, (title, P1, P2, L, R, x0, t, grids) in CASES.items():
    xe = np.linspace(0, 1, 4001); ex = exact(c, xe, t)
    fig, ax = plt.subplots(2, 2, figsize=(9, 6)); ax = ax.ravel()
    q = [('rho', r'$\rho$ (kg/m$^3$)'), ('u', '$u$ (m/s)'), ('p', '$p$ (Pa)'), ('M', 'Mach number $|u|/c_W$')]
    ex['M'] = np.abs(ex['u']) / ex['cw']
    adv = c in (80, 86)   # uniform u, p: plot the deviations u - u0, p/p0 - 1 (the exact ones are zero)
    if adv:
        q[1] = ('du', '$u-u_0$ (m/s)'); q[2] = ('dp', '$p/p_0-1$'); ex['du'] = 0 * xe; ex['dp'] = 0 * xe
    for k, (v, lab) in enumerate(q): ax[k].plot(xe, ex[v], 'k-', lw=1.0, label='exact')
    runs = [(n, '', 'IMEX %d cells' % n) for n in grids]
    PLOT_N = 1000   # one relaxation level: the figure shows the 1000-cell solution (black hollow circles); the table has all grids
    for n, tag, name in runs:
        s = load(c, n, tag)
        if s is None: continue
        s['M'] = np.abs(s['u']) / s['cw']; s['du'] = s['u'] - L[3]; s['dp'] = s['p'] / L[4] - 1
        mk = 'o'
        if n == PLOT_N:
            for k, (v, lab) in enumerate(q): ax[k].plot(s['x'], s[v], ls='none', marker=mk, mfc='none', ms=3.5, mew=0.8, color='k', markevery=0.012, label=name)
        dx = 1.0 / n; xc = (np.arange(n) + 0.5) * dx
        num = np.sum(np.abs(s['rho'] - exact(c, xc, t)['rho'])) * dx; den = np.sum(np.abs(exact(c, xc, 0)['rho'])) * dx
        rows.append((c, title + (' (dt = 0.5 dx/c)' if tag else ''), n, num / den))
    for k, (v, lab) in enumerate(q):
        ax[k].set_xlabel('$x$ (m)'); ax[k].set_ylabel(lab); ax[k].grid(alpha=.3)
        if (v == 'rho' and c in (80, 83, 86)) or (v == 'p' and c == 83): ax[k].set_yscale('log')
    ax[0].legend(); fig.suptitle('%s (Schropff et al. 2026), $t$ = %g s; IMEX ZB split, LSDIRK2, material CFL 0.2' % (title, t))
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'sch%d.png' % c), dpi=130); plt.close(fig)
print('| case | test | cells | L1(rho) relative, Eq. (51) |\n|---|---|---|---|')
for c, title, n, e in rows: print('| %d | %s | %d | %.3e |' % (c, title, n, e))
