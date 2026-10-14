"""Reproductions of the paper figures with the IMEX program (src/SPB6_TV_IMEX.cpp, V16): the variables, units and relaxation levels
of each paper figure, against the exact solutions of the relaxed limit models:
  5-eq. (p-relaxed, tests/exact5.py), 4-eq. (pT-relaxed) and activation-limited pT -> pTg (tests/exact_relaxed.py).
usage: python3 tests/plot_paper_figs.py RUNDIR OUTDIR [FIG ...]   (RUNDIR from tests/run_archive.sh)"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import exact5 as e5, exact_relaxed as er, exact_mg as emg

RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, '..', 'run', 'figs')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H, '..', 'docs', 'figures')
ONLY = sys.argv[3:]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 7, 'lines.linewidth': 1.1})
# per-test settings of src/spb6_physics.h test_setup (material CFL; Rusanov speed local unless noted)
SET = {'spb19': 'material CFL 0.05', 'ps3': 'material CFL 0.05', 'spb20': 'material CFL 0.2, global Rusanov speed'}
def imex_lab(fig): return 'IMEX (TV-split, LSDIRK2, %s)' % SET.get(fig, 'material CFL 0.2')
# archive names (tests/run_archive.sh) of the figure runs
ALIAS = {'spb5': 'spb41', 'spb6': 'spb42', 'spb8': 'spb431', 'spb10': 'spb432', 'spb12': 'spb45', 'spb19': 'spb51', 'spb20': 'spb52',
         'pel1': 'pel811', 'pel4': 'pel821', 'pel6': 'pel822'}

# ---------------------------------------------------------------------------------------------------- EOS of the cases (src/spb6_physics.h)
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
    if name.endswith('_I'): name = name[:-2]
    head = name.split('_', 1); name = ALIAS.get(head[0], head[0]) + ('_' + head[1] if len(head) > 1 else '')
    f = os.path.join(RUN, name + '.dat')
    if not os.path.exists(f): print('  missing', f); return None
    d = np.loadtxt(f); return {'x': d[:, 0], 'p': d[:, 2], 'u': d[:, 3], 'rho': d[:, 5], 'a2': d[:, 6], 'Y2': d[:, 7], 'r1': d[:, 8],
                               'r2': d[:, 9], 'cf': d[:, 10], 'cw': d[:, 11], 'T1': d[:, 12], 'T2': d[:, 13], 'p1': d[:, 14], 'p2': d[:, 15]}


def derived(s, P):   # quantities plotted in the papers
    s['pbar'] = s['p'] / 1e5; s['a1'] = 1 - s['a2']; s['pGPa'] = s['p'] / 1e9
    with np.errstate(invalid='ignore', divide='ignore'):
        if 'cf' in s: s['Mf'] = np.abs(s['u']) / s['cf']
        if 'cw' in s: s['Mw'] = np.abs(s['u']) / s['cw']
        if 'T1' in s and P is not None:
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
    xs = [s['x'] for s, kw in curves if s is not None and v in s]
    ax.set_xlabel('x (m)'); ax.set_ylabel(lab); ax.set_xlim(min(x[0] for x in xs), max(x[-1] for x in xs)) if xs else ax.set_xlim(0, 1)
    if log: ax.set_yscale('log')
    if ylim: ax.set_ylim(*ylim); ax.ticklabel_format(axis='y', useOffset=False)
    ax.grid(alpha=0.25)


def single(name, title, cases, var, exact, ncol=2):
    """one numerical setting: IMEX and the exact solution overlaid in each panel (SPB 2009 figures)"""
    I = load(cases + '_I'); P = exact[1]
    I = derived(I, P) if I else None
    nr = (len(var) + ncol - 1) // ncol; fig, axs = plt.subplots(nr, ncol, figsize=(4.6 * ncol, 2.9 * nr), squeeze=False)
    curves = [(exact[0], dict(color='k', lw=1.0, label='exact (5-eq.)' if P is not None else 'exact (Euler)')),
              (I, dict(color='tab:red', marker='.', ms=2.5, lw=0.9, label='IMEX'))]
    for ax, (v, lab, *opt) in zip(axs.flat, var): panel(ax, v, curves, lab, *opt)
    for ax in list(axs.flat)[len(var):]: ax.axis('off')
    axs.flat[0].legend(loc='best'); fig.suptitle(title + '\n' + imex_lab(name), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)


def levels(name, title, prefix, lev, var, exacts, P):
    """several relaxation levels (IMEX), levels in the colours of the paper, exact solutions in black"""
    fig, axs = plt.subplots(len(var), 1, figsize=(6.6, 3.1 * len(var)), squeeze=False)
    for j, (kind, head) in enumerate((('I', imex_lab(name)),)):
        curves = [(s, dict(kw, label=l)) for s, kw, l in exacts]
        for key, l, kw in lev:
            s = load('%s_%s' % (prefix, key)); curves.append((derived(s, P) if s else None, dict(kw, label=l)))
        for i, (v, lab, *opt) in enumerate(var):
            panel(axs[i, j], v, curves, lab, *opt)
            if v in ('T12', 'g12'): two_phase(axs[i, j], v, curves)
        axs[0, j].set_title(head); axs[0, j].legend(loc='best', ncol=2)
    import textwrap; title = textwrap.fill(title, 95)
    fig.suptitle(title, fontsize=9); fig.tight_layout(rect=(0, 0, 1, 1 - (0.2 + 0.15 * title.count('\n')) / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)



def paper_style(name, title, prefix, lev, var, P, ncol=2):
    """the paper's own figure: only the numerical curves of the relaxation levels the paper shows (no exact solutions)"""
    nr = (len(var) + ncol - 1) // ncol; fig, axs = plt.subplots(nr, ncol, figsize=(4.8 * ncol, 3.0 * nr), squeeze=False)
    curves = []
    for key, l, kw in lev:
        s_ = load('%s_%s' % (prefix, key)); curves.append((derived(s_, P) if s_ else None, dict(kw, label=l)))
    for ax, (v, lab, *opt) in zip(axs.flat, var): panel(ax, v, curves, lab, *opt)
    for ax in list(axs.flat)[len(var):]: ax.axis('off')
    axs.flat[0].legend(loc='best'); fig.suptitle(title + '\n' + imex_lab(name), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.45 / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)


def vs_exact(name, title, s, P, ex, exlab, var, centre):
    """one relaxation level against its exact solution; pressure also on a log axis; centre values written in each panel"""
    s = derived(s, P) if s else None
    fig, axs = plt.subplots(1, len(var), figsize=(3.6 * len(var), 3.2), squeeze=False)
    for ax, (v, lab, *opt) in zip(axs.flat, var):
        panel(ax, v, [(ex, dict(color='k', ls=':', lw=1.3, label=exlab)), (s, dict(color='b', lw=1.0, label='IMEX'))], lab, *opt)
        if centre and s is not None and v in s:
            i = len(s['x']) // 2; j = len(ex['x']) // 2
            ax.text(0.5, 0.04, 'x = 0.5: IMEX %.4g, exact %.4g' % (s[v][i], ex[v][j]), transform=ax.transAxes, ha='center', fontsize=7,
                    bbox=dict(fc='w', ec='0.7', lw=0.5))
    axs.flat[0].legend(loc='upper right'); fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.35 / fig.get_figheight())); fig.savefig(os.path.join(OUT, name + '.png'), dpi=110); plt.close(fig); print('wrote', name)

def two_phase(ax, v, curves):   # liquid solid, vapour dashed (Pelanti Figs. 5-6)
    k1, k2 = ('T1', 'T2') if v == 'T12' else ('g1', 'g2')
    for s, kw in curves:
        if s is None or k1 not in s: continue
        kw = dict(kw); l = kw.pop('label', ''); kw.pop('ls', None)
        ax.plot(s['x'], s[k1], ls='-', label=l + ' (l)', **kw); ax.plot(s['x'], s[k2], ls='--', label=l + ' (v)', **kw)


def want(n): return not ONLY or n in ONLY


def exact_mg(case, x0, x1, n=4000):
    """Euler exact solution with the MG EOS of the case (tests/exact_mg.py)"""
    eL, sL, eR, sR, xd, t = emg.CASES[case]
    L = emg.Side(eL(), *sL, -1); R = emg.Side(eR(), *sR, +1); ps, us, rL, rR = emg.solve(L, R)
    x = np.linspace(x0, x1, n); o = emg.sample(L, R, ps, us, rL, rR, (x - xd) / t)
    a1 = np.where(x < xd + us * t, 1.0, 0.0)
    return {'x': x, 'rho': o[:, 0], 'u': o[:, 1], 'p': o[:, 2], 'a2': 1 - a1, 'pbar': o[:, 2] / 1e5, 'a1': a1, 'pGPa': o[:, 2] / 1e9}


def gpa(s):
    if s: s['pGPa'] = s['p'] / 1e9
    return s


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
    ex = exact5(P, (1 - 1e-6, 1000., 10., 0., 1e12), (1e-6, 1000., 10., 0., 1e5), 0.6, 8.3e-6)   # interface 0.6 m (inferred from the paper's Fig. 10)
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
    # as in the paper: only the four relaxation levels
    paper_style(nm + '_paper', 'Pelanti & Shyue 2014 Fig. %d as in the paper: |u| = %g m/s, t = %g ms, 5000 cells, four relaxation levels' % (fig, U0, t * 1e3), nm,
                PS_LEV, [('p', 'p (Pa)'), ('u', 'u (m/s)'), ('a2', r'$\alpha_v$'), ('Y2', '$Y_v$')], P)
    # p-pT-pTg against its exact solution (pT -> pTg fan), pressure also on a log axis
    vs_exact(nm + '_exact', 'Pelanti & Shyue 2014 Fig. %d, p-pT-pTg level against the exact pT->pTg fan (IMEX, 5000 cells)' % fig, load(nm + '_r2'), P, exg,
             'exact pT->pTg', [('p', 'p (Pa)'), ('p', 'p (Pa), log axis', True), ('Y2', '$Y_v$')], True)
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
if want('pel4'):
    paper_style('pel4_paper', 'Pelanti 2022 Fig. 4 as in the paper: NASG water, |u| = 2 m/s, t = 3 ms, 5000 cells, five relaxation levels', 'pel4',
                PEL_LEV, [('p', 'p (Pa)'), ('u', 'u (m/s)'), ('a2', r'$\alpha_2$'), ('Y2', '$Y_2$')], P)
    vs_exact('pel4_exact', 'Pelanti 2022 Fig. 4, pTg level against the exact pT->pTg fan (IMEX, 5000 cells)', load('pel4_pTg'), P, exg,
             'exact pT->pTg', [('p', 'p (Pa)'), ('a2', r'$\alpha_2$'), ('Y2', '$Y_2$')], True)
if want('pel5'):   # Fig. 5: temperature and chemical potentials; T_l = T_v in all three levels, so one T line per level
    L5 = [('pT', 'pT', dict(color='darkgreen')), ('pTgf', 'pTg(f)', dict(color='deepskyblue')), ('pTg', 'pTg', dict(color='navy'))]
    S5 = [(k, l, kw, derived(load('pel4_' + k), P)) for k, l, kw in L5]
    fig = plt.figure(figsize=(10.5, 6.6)); gs = fig.add_gridspec(2, 2, width_ratios=(1.6, 1))
    a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[0, 1]); a2 = fig.add_subplot(gs[1, :])
    for ax, xl, yl in ((a0, (0, 1), (352.92, 353.01)), (a1, (0.48, 0.52), (352.55, 353.01))):
        ax.plot(ex4['x'], ex4['T1'], color='k', ls='--', lw=0.9, label='exact 4-eq. (pT)')
        ax.plot(exg['x'], exg['T1'], color='0.5', ls=':', lw=1.4, label='exact pT->pTg')
        for k, l, kw, s_ in S5:
            if s_: ax.plot(s_['x'], s_['T1'], label=l, lw=1.2, **kw)
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.ticklabel_format(axis='y', useOffset=False); ax.grid(alpha=0.25); ax.set_xlabel('x (m)'); ax.set_ylabel('T (K)')
    a0.set_title('T (liquid = vapour), axis of the paper', fontsize=9); a1.set_title('centre, full range', fontsize=9); a0.legend(loc='lower left', fontsize=7)
    for k, l, kw, s_ in S5:
        if s_: a2.plot(s_['x'], s_['g1'], ls='-', label=l + ' liquid', **kw); a2.plot(s_['x'], s_['g2'], ls='--', label=l + ' vapour', **kw)
    a2.set_xlim(0, 1); a2.set_xlabel('x (m)'); a2.set_ylabel('$g_k$ (J/kg)'); a2.grid(alpha=0.25); a2.legend(loc='lower left', ncol=3, fontsize=7)
    fig.suptitle('Pelanti 2022 Fig. 5: temperature and chemical potentials, NASG water, t = 3 ms, 5000 cells\n' + imex_lab('pel5'), fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(OUT, 'pel5.png'), dpi=110); plt.close(fig); print('wrote pel5')
P = PH['pel_d']
if want('pel6'):   # Fig. 6: dodecane (NASG), p-relax and pTg-relax, 2000 cells, t = 473 us
    ex = exact5(P, (1 - 1e-8, 500., 2., 0., 1e8), (1e-8, 500., 2., 0., 1e5), 0.75, 473e-6)
    levels('pel6', 'Pelanti 2022 Fig. 6: dodecane shock tube (NASG), t = 473 us, 2000 cells', 'pel6',
           [('r0', 'p-relax', dict(color='darkviolet')), ('r2', 'pTg-relax', dict(color='b'))],
           [('rho', r'$\rho$ (kg/m$^3$)', True), ('u', 'u (m/s)'), ('Y2', '$Y_2$'), ('p', 'p (Pa)', True), ('a2', r'$\alpha_2$'), ('T12', '$T_k$ (K)', False, (300, 1200))],
           [(ex, dict(color='k', lw=0.9), 'exact 5-eq.')], P)


# ---------------------------------------------------------------------------------------------------- Mie-Grueneisen (V12)
if want('spb6'):   # SPB 2009 Figs. 6-7: nitromethane (Cochran-Chan), 20 GPa | 0.2 MPa, 500 cells, t = 67 us
    ex = exact_mg(60, 0, 1)
    single('spb6', 'SPB 2009 Fig. 6: nitromethane shock tube (Cochran-Chan), 20 GPa | 0.2 MPa, t = 67 us, 500 cells', 'spb6',
           [('pGPa', 'pressure (GPa)'), ('u', 'velocity (m/s)'), ('rho', 'density (kg/m$^3$)'), ('a1', 'volume fraction of "phase" 1')],
           (ex, None))
    I = gpa(load('spb6_I'))
    if I:
        fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.0))
        for ax, (v, lab, ylim) in zip(axs, [('pGPa', 'pressure (GPa)', (9.6, 10.0)), ('u', 'velocity (m/s)', (1740, 1790))]):
            for s_, kw in ((ex, dict(color='k', lw=1.0, label='exact')), (I, dict(color='tab:red', ls='--', marker='x', ms=3, label='IMEX'))):
                ax.plot(s_['x'], s_[v], **kw)
            ax.set_xlim(0.55, 0.75); ax.set_ylim(*ylim); ax.set_xlabel('x (m)'); ax.set_ylabel(lab); ax.grid(alpha=0.25)
        axs[0].legend(); fig.suptitle('SPB 2009 Fig. 7: magnified around the contact (x = 0.5 + u* t = 0.619 m)\n' + imex_lab('spb7'), fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, 'spb7.png'), dpi=110); plt.close(fig); print('wrote spb7')
if want('pel1'):   # Pelanti 2022 Fig. 1: detonation gases (JWL) | water (NASG), 500 cells, t = 0.9 ms
    ex = exact_mg(61, -10, 10)
    single('pel1', 'Pelanti 2022 Fig. 1: JWL detonation gases | NASG water, 4.64e10 Pa | 1e5 Pa, t = 0.9 ms, 500 cells', 'pel1',
           [('rho', 'density (kg/m$^3$)'), ('u', 'velocity (m/s)'), ('pGPa', 'pressure (GPa)'), ('a1', 'volume fraction of gas')],
           (ex, None))
