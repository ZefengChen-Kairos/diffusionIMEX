"""Reproductions of the paper figures: explicit (E, src/SPB6_explicit.cpp) and IMEX (I, src/SPB6_TV_IMEX.cpp) results with the
variables, units and relaxation levels of each paper figure, against the exact solutions of the relaxed limit models:
  5-eq. (p-relaxed, tests/exact5.py), 4-eq. (pT-relaxed) and activation-limited pT -> pTg (tests/exact_relaxed.py).
usage: python3 tests/plot_paper_figs.py RUNDIR OUTDIR [FIG ...]   (RUNDIR from tests/run_paper_figs.sh)"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import exact5 as e5, exact_relaxed as er

RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, '..', 'run', 'figs')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H, '..', 'docs', 'figures')
ONLY = sys.argv[3:]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 7, 'lines.linewidth': 1.1})
EXPL = 'explicit (HLLC + relaxation, Heun, CFL 0.4)'
IMEX = 'IMEX (TV-split, LSDIRK2, material CFL 0.2)'
BOTH = 'explicit = HLLC + relaxation, Heun, CFL 0.4;  IMEX = TV-split, LSDIRK2, material CFL 0.2'

# ---------------------------------------------------------------------------------------------------- EOS of the cases (src/spb6_common.h)
PH = {'spb_wa': [(4.4, 6e8, 0, 0, 1000., 0), (1.4, 0, 0, 0, 1000., 0)],
      'spb_es': [(2.43, 5.3e9, 0, 0, 1000., 0), (1.62, 141e9, 0, 0, 1000., 0)],
      'ps_w': [(2.35, 1e9, -1167e3, 0.0, 1816.0, 0.0), (1.43, 0.0, 2030e3, -23.4e3, 1040.0, 0.0)],
      'ps_d': [(2.35, 4e8, -775.269e3, 0.0, 1077.7, 0.0), (1.025, 0.0, -237.547e3, -24.4e3, 1956.45, 0.0)],
      'pel_w': [(1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4), (1.467, 0.0, 2077616.0, 14317.0, 955.0, 0.0)],
      'pel_d': [(1.206, 1681e5, -996054.0, 0.0, 2532.0, 7.51e-4), (1.021, 0.0, -384592.0, -4301.0, 2274.0, 0.0)]}


def gibbs(P, k, p, T):
    g, pi, q, qp, cv, b = P[k]
    with np.errstate(invalid='ignore', divide='ignore'):
        return g * cv * T + b * p + q - T * (cv * (g * np.log(T) - (g - 1) * np.log(p + pi)) + qp)


def load(name):
    f = os.path.join(RUN, name + '.dat')
    if not os.path.exists(f): print('  missing', f); return None
    d = np.loadtxt(f); return {'x': d[:, 0], 'p': d[:, 2], 'u': d[:, 3], 'rho': d[:, 5], 'a2': d[:, 6], 'Y2': d[:, 7], 'r1': d[:, 8],
                               'r2': d[:, 9], 'cf': d[:, 10], 'cw': d[:, 11], 'T1': d[:, 12], 'T2': d[:, 13], 'p1': d[:, 14], 'p2': d[:, 15]}


def derived(s, P):   # quantities plotted in the papers
    s['pbar'] = s['p'] / 1e5; s['a1'] = 1 - s['a2']
    with np.errstate(invalid='ignore', divide='ignore'):
        if 'cf' in s: s['Mf'] = np.abs(s['u']) / s['cf']
        if 'cw' in s: s['Mw'] = np.abs(s['u']) / s['cw']
        if 'T1' in s:
            s['dT'] = s['T2'] - s['T1']
            p1 = s.get('p1', s['p']); p2 = s.get('p2', s['p'])
            s['g1'] = gibbs(P, 0, p1, s['T1']); s['g2'] = gibbs(P, 1, p2, s['T2'])
    return s


def exact5(P, L, R, x0, t, n=4000):
    """L, R = (alpha_0, rho_0, rho_1, u, p): 5-eq. exact solution with the c_f, c_W and T_k of the exact state"""
    x = np.linspace(0, 1, n); G = [P[0][0], P[1][0]]; PI = [P[0][1], P[1][1]]; B = [P[0][5], P[1][5]]
    side = lambda a0, r0, r1, u, p: e5.Side([a0, 1 - a0], [r0, r1], u, p, G, PI, B)
    o, ps, us = e5.solve(side(*L), side(*R), x0, t, x)
    s = {'x': x, 'p': o['p'], 'u': o['u'], 'rho': o['rho'], 'a2': 1 - o['a1'], 'r1': o['r1'], 'r2': o['r2']}
    with np.errstate(invalid='ignore', divide='ignore'):
        a = [o['a1'], 1 - o['a1']]; r = [o['r1'], o['r2']]; c2 = [P[k][0] * (s['p'] + P[k][1]) / (r[k] * (1 - r[k] * P[k][5])) for k in (0, 1)]
        s['Y2'] = a[1] * r[1] / s['rho']
        s['cf'] = np.sqrt((1 - s['Y2']) * c2[0] + s['Y2'] * c2[1]); s['cw'] = np.sqrt(1 / (s['rho'] * (a[0] / (r[0] * c2[0]) + a[1] / (r[1] * c2[1]))))
        s['T1'], s['T2'] = [(1 / r[k] - P[k][5]) * (s['p'] + P[k][1]) / ((P[k][0] - 1) * P[k][4]) for k in (0, 1)]
    vac = s['rho'] <= 0   # vacuum of the 5-eq. solution: p = 0, rho = 0, alpha_2 = 1 (as plotted in SPB Fig. 12); the rest undefined
    for k in ('r1', 'r2', 'Y2', 'cf', 'cw', 'T1', 'T2'): s[k] = np.where(vac, np.nan, s[k])
    return derived(s, P)


def fan(P, p0, T0, av, U0, t, chem):
    F = er.Fan(P, p0, T0, av); (x, p, u, Y, a, T), info = F.full(U0, t, chem)
    s = {'x': x, 'p': p, 'u': u, 'Y2': Y, 'a2': a, 'T1': T, 'T2': T}; derived(s, P); return s, info


def ic_rho(P, p, T): return [1 / er.Fan(P, p, T, 0.5).v(k, p, T) for k in (0, 1)]


# ---------------------------------------------------------------------------------------------------- figure helpers
def panel(ax, v, curves, lab, log=False, ylim=None):
    for s, kw in curves:
        if s is not None and v in s: ax.plot(s['x'], s[v], **kw)
    ax.set_xlabel('x (m)'); ax.set_ylabel(lab); ax.set_xlim(0, 1)
    if log: ax.set_yscale('log')
    if ylim: ax.set_ylim(*ylim); ax.ticklabel_format(axis='y', useOffset=False)
    ax.grid(alpha=0.25)


def single(name, title, cases, var, exact, ncol=2):
    """one numerical setting: E, I and the exact solution overlaid in each panel (SPB 2009 figures)"""
    E, I = load(cases + '_E'), load(cases + '_I'); P = exact[1]
    E = derived(E, P) if E else None; I = derived(I, P) if I else None
    nr = (len(var) + ncol - 1) // ncol; fig, axs = plt.subplots(nr, ncol, figsize=(4.6 * ncol, 2.9 * nr), squeeze=False)
    curves = [(exact[0], dict(color='k', lw=1.0, label='exact (5-eq.)')), (E, dict(color='tab:blue', label='explicit')),
              (I, dict(color='tab:red', ls='--', label='IMEX'))]
    for ax, (v, lab, *opt) in zip(axs.flat, var): panel(ax, v, curves, lab, *opt)
    for ax in list(axs.flat)[len(var):]: ax.axis('off')
    axs.flat[0].legend(loc='best'); fig.suptitle(title + '\n' + BOTH, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)


def levels(name, title, prefix, lev, var, exacts, P):
    """several relaxation levels: column E | column I, levels in the colours of the paper, exact solutions in black"""
    fig, axs = plt.subplots(len(var), 2, figsize=(10, 2.75 * len(var)), squeeze=False)
    for j, (kind, head) in enumerate((('E', EXPL), ('I', IMEX))):
        curves = [(s, dict(kw, label=l)) for s, kw, l in exacts]
        for key, l, kw in lev:
            s = load('%s_%s_%s' % (prefix, key, kind)); curves.append((derived(s, P) if s else None, dict(kw, label=l)))
        for i, (v, lab, *opt) in enumerate(var):
            panel(axs[i, j], v, curves, lab, *opt)
            if v in ('T12', 'g12'): two_phase(axs[i, j], v, curves)
        axs[0, j].set_title(head); axs[0, j].legend(loc='best', ncol=2)
    fig.suptitle(title, fontsize=10); fig.tight_layout(rect=(0, 0, 1, 1 - 0.35 / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)


def two_phase(ax, v, curves):   # liquid solid, vapour dashed (Pelanti Figs. 5-6)
    k1, k2 = ('T1', 'T2') if v == 'T12' else ('g1', 'g2')
    for s, kw in curves:
        if s is None or k1 not in s: continue
        kw = dict(kw); l = kw.pop('label', ''); kw.pop('ls', None)
        ax.plot(s['x'], s[k1], ls='-', label=l + ' (l)', **kw); ax.plot(s['x'], s[k2], ls='--', label=l + ' (v)', **kw)


def want(n): return not ONLY or n in ONLY


# ---------------------------------------------------------------------------------------------------- SPB 2009 (Saurel, Petitpas, Berry)
P = PH['spb_wa']
if want('spb5'):   # Fig. 5: interface advection, 200 cells, t = 2.79 ms (exact: translated initial data)
    x = np.linspace(0, 1, 4000); L = x < 0.5 + 100 * 2.79e-3; a1 = np.where(L, 1 - 1e-8, 1e-8)
    s = {'x': x, 'p': np.full_like(x, 1e5), 'u': np.full_like(x, 100.), 'a2': 1 - a1, 'r1': np.full_like(x, 1000.), 'r2': np.full_like(x, 10.)}
    s['rho'] = a1 * 1000 + (1 - a1) * 10; c2 = [P[k][0] * (1e5 + P[k][1]) / r for k, r in ((0, 1000.), (1, 10.))]; s['Y2'] = (1 - a1) * 10 / s['rho']
    s['cf'] = np.sqrt((1 - s['Y2']) * c2[0] + s['Y2'] * c2[1]); s['cw'] = np.sqrt(1 / (s['rho'] * (a1 / (1000 * c2[0]) + (1 - a1) / (10 * c2[1]))))
    single('spb5', 'SPB 2009 Fig. 5: advection of a water-air interface, u = 100 m/s, t = 2.79 ms, 200 cells', 'spb5',
           [('pbar', 'pressure (bar)'), ('u', 'velocity (m/s)'), ('rho', 'mixture density (kg/m$^3$)'), ('a1', 'volume fraction of water'),
            ('Mf', 'Mach number, frozen sound speed'), ('Mw', 'Mach number, Wood sound speed')], (derived(s, P), P))
if want('spb8'):   # Fig. 8: water-air shock tube 1 GPa / 0.1 MPa, 1000 cells, t = 240 us
    ex = exact5(P, (1 - 1e-6, 1000., 1., 0., 1e9), (1e-6, 1000., 1., 0., 1e5), 0.75, 240e-6)
    single('spb8', 'SPB 2009 Fig. 8: water-air shock tube, 1 GPa | 0.1 MPa, t = 240 us, 1000 cells', 'spb8',
           [('pbar', 'pressure (bar)'), ('u', 'velocity (m/s)'), ('rho', 'mixture density (kg/m$^3$)'), ('a1', 'volume fraction of water')], (ex, P))
if want('spb10'):  # Fig. 10: water-air shock tube 1 TPa, air 10 kg/m3, 1000 cells, t = 8.3 us
    ex = exact5(P, (1 - 1e-6, 1000., 10., 0., 1e12), (1e-6, 1000., 10., 0., 1e5), 0.75, 8.3e-6)
    single('spb10', 'SPB 2009 Fig. 10: water-air shock tube, 1 TPa | 0.1 MPa, t = 8.3 us, 1000 cells', 'spb10',
           [('pbar', 'pressure (bar)'), ('u', 'velocity (m/s)'), ('Mw', 'Mach number (Wood sound speed)'), ('rho', 'mixture density (kg/m$^3$)'),
            ('a1', 'volume fraction of water')], (ex, P))
if want('spb12'):  # Fig. 12: expansion tube, alpha_air = 1e-2, -+100 m/s, 1000 cells, t = 1.85 ms (5-eq. solution has a vacuum)
    ex = exact5(P, (0.99, 1000., 1., -100., 1e5), (0.99, 1000., 1., 100., 1e5), 0.5, 1.85e-3)
    single('spb12', 'SPB 2009 Fig. 12: expansion tube with cavitation, u = -+100 m/s, t = 1.85 ms, 1000 cells', 'spb12',
           [('pbar', 'pressure (bar)'), ('u', 'velocity (m/s)'), ('rho', 'mixture density (kg/m$^3$)'), ('a2', 'volume fraction of air')], (ex, P))
P = PH['spb_es']
for nm, pl, t, fig in (('spb19', 1e10, 80e-6, 19), ('spb20', 2e11, 34e-6, 20)):
    if not want(nm): continue
    ex = exact5(P, (0.5954, 1185., 3622., 0., pl), (0.5954, 1185., 3622., 0., 1e5), 0.6, t)
    single(nm, 'SPB 2009 Fig. %d: epoxy-spinel shock tube, %.0e Pa | 1e5 Pa, t = %.0f us, 500 cells' % (fig, pl, t * 1e6), nm,
           [('pbar', 'pressure (bar)'), ('u', 'velocity (m/s)'), ('Mw', 'Mach number (Wood sound speed)'), ('rho', 'mixture density (kg/m$^3$)'),
            ('a1', 'volume fraction of epoxy'), ('r2', 'density of spinel (kg/m$^3$)')], (ex, P))

# ---------------------------------------------------------------------------------------------------- Pelanti & Shyue 2014
P = PH['ps_w']; T0 = 354.728; r0 = ic_rho(P, 1e5, T0)
PS_LEV = [('r0', 'p', dict(color='m')), ('r1', 'p-pT', dict(color='r', ls=':')), ('r2', 'p-pT-pTg', dict(color='b')), ('r3', 'p-pTg', dict(color='g', ls='--'))]
for nm, U0, t, fig in (('ps1', 2., 3.2e-3, 1), ('ps2', 500., 0.58e-3, 2)):
    if not want(nm): continue
    ex5 = exact5(P, (0.99, r0[0], r0[1], -U0, 1e5), (0.99, r0[0], r0[1], U0, 1e5), 0.5, t)
    ex4, i4 = fan(P, 1e5, T0, 1e-2, U0, t, False); exg, ig = fan(P, 1e5, T0, 1e-2, U0, t, True)
    print('  %s exact centre p: 5-eq %.4g, 4-eq %s, pT->pTg %.4g (Y %.4g)' % (nm, np.nanmin(ex5['p']), '%.4g' % i4['p_end1'] if i4['centre'] else 'vacuum', exg['p'][len(exg['p']) // 2], exg['Y2'][len(exg['p']) // 2]))
    levels(nm, 'Pelanti & Shyue 2014 Fig. %d: water cavitation tube, |u| = %g m/s, t = %g ms, 5000 cells (SG water, T = 354.728 K, alpha_v = 1e-2)' % (fig, U0, t * 1e3), nm,
           PS_LEV,
           [('p', 'p (Pa)'), ('u', 'u (m/s)'), ('a2', r'$\alpha_v$'), ('Y2', '$Y_v$')],
           [(ex5, dict(color='k', lw=0.9), 'exact 5-eq. (p)'), (ex4, dict(color='k', ls='--', lw=0.9), 'exact 4-eq. (pT)'),
            (exg, dict(color='k', ls=':', lw=1.2), 'exact pT->pTg')], P)
P = PH['ps_d']
if want('ps3'):    # Fig. 3: dodecane liquid-vapour shock tube, 2000 cells, t = 473 us
    ex = exact5(P, (1 - 1e-8, 500., 2., 0., 1e8), (1e-8, 500., 2., 0., 1e5), 0.75, 473e-6)
    levels('ps3', 'Pelanti & Shyue 2014 Fig. 3: dodecane liquid-vapour shock tube (SG), t = 473 us, 2000 cells', 'ps3',
           [('r0', 'p', dict(color='m')), ('r2', 'p-pT-pTg', dict(color='b')), ('r3', 'p-pTg', dict(color='g', ls='--'))],
           [('rho', r'$\rho$ (kg/m$^3$)', True), ('u', 'u (m/s)'), ('a2', r'$\alpha_v$'), ('Y2', '$Y_v$'), ('pbar', 'p (bar)', True), ('dT', '$T_v - T_l$ (K)')],
           [(ex, dict(color='k', lw=0.9), 'exact 5-eq. (p)')], P)

# ---------------------------------------------------------------------------------------------------- Pelanti 2022 (IJMF)
P = PH['pel_w']; r0 = ic_rho(P, 1e5, 353.0)
PEL_LEV = [('p', 'p-relax', dict(color='darkviolet')), ('pTf', 'pT(f)', dict(color='limegreen', ls='--')), ('pT', 'pT', dict(color='darkgreen')),
           ('pTgf', 'pTg(f)', dict(color='deepskyblue', ls='--')), ('pTg', 'pTg', dict(color='navy'))]
if want('pel4') or want('pel5'):
    ex5 = exact5(P, (0.99, r0[0], r0[1], -2., 1e5), (0.99, r0[0], r0[1], 2., 1e5), 0.5, 3e-3)
    ex4, i4 = fan(P, 1e5, 353.0, 1e-2, 2., 3e-3, False); exg, ig = fan(P, 1e5, 353.0, 1e-2, 2., 3e-3, True)
    print('  pel4 exact centre p: 5-eq %.5g, 4-eq %.5g, pT->pTg %.5g' % (np.nanmin(ex5['p']), i4['p_end1'], ig['p_end2']))
if want('pel4'):   # Fig. 4: p, u, alpha_2, Y_2 for five relaxation levels, 5000 cells, t = 3 ms
    levels('pel4', 'Pelanti 2022 Fig. 4: water cavitation tube (NASG), |u| = 2 m/s, t = 3 ms, 5000 cells; pT(f): theta = 2000, pTg(f): nu = 1e-4', 'pel4',
           PEL_LEV, [('p', 'p (Pa)'), ('u', 'u (m/s)'), ('a2', r'$\alpha_2$'), ('Y2', '$Y_2$')],
           [(ex5, dict(color='k', lw=0.9), 'exact 5-eq.'), (ex4, dict(color='k', ls='--', lw=0.9), 'exact 4-eq.'),
            (exg, dict(color='k', ls=':', lw=1.2), 'exact pT->pTg')], P)
if want('pel5'):   # Fig. 5: temperature and chemical potentials (cases with temperature equilibrium)
    levels('pel5', 'Pelanti 2022 Fig. 5: temperatures and chemical potentials g_k (l solid, v dashed), t = 3 ms, 5000 cells\n'
           'T axis as in the paper (352.92-353.01 K); the pTg dip goes off scale (exact pT->pTg: 352.63 K at x = 0.5)', 'pel4',
           [l for l in PEL_LEV if l[0] in ('pT', 'pTgf', 'pTg')], [('T12', 'T (K)', False, (352.92, 353.01)), ('g12', '$g_k$ (J/kg)')],
           [(ex4, dict(color='k', lw=0.9), 'exact 4-eq.'), (exg, dict(color='grey', lw=0.9), 'exact pT->pTg')], P)
P = PH['pel_d']
if want('pel6'):   # Fig. 6: dodecane (NASG), p-relax and pTg-relax, 2000 cells, t = 473 us
    ex = exact5(P, (1 - 1e-8, 500., 2., 0., 1e8), (1e-8, 500., 2., 0., 1e5), 0.75, 473e-6)
    levels('pel6', 'Pelanti 2022 Fig. 6: dodecane shock tube (NASG), t = 473 us, 2000 cells', 'pel6',
           [('r0', 'p-relax', dict(color='darkviolet')), ('r2', 'pTg-relax', dict(color='b'))],
           [('rho', r'$\rho$ (kg/m$^3$)', True), ('u', 'u (m/s)'), ('Y2', '$Y_2$'), ('p', 'p (Pa)', True), ('a2', r'$\alpha_2$'), ('T12', '$T_k$ (K)', False, (300, 1200))],
           [(ex, dict(color='k', lw=0.9), 'exact 5-eq.')], P)
