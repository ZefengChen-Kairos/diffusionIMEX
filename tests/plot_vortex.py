#!/usr/bin/env python3
"""Battisti-Boscheri 4.1 two-phase isentropic vortex (case 56), IMEX: Fig. 1 analogue + convergence.
usage: python3 tests/plot_vortex.py RUNDIR OUT.png"""
import sys, os, re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d, out = sys.argv[1], sys.argv[2]
def exact_p(x, y, p0, g=1.4, eps=5.0):   # 精确解 = 初值（点值）
    dT = -(g - 1) * eps**2 / (8 * g * np.pi**2) * np.exp(1 - x * x - y * y); T = p0 + dT
    return (T / p0)**(1 / (g - 1)) * T
fig, ax = plt.subplots(2, 3, figsize=(15, 9))
for k, p0 in enumerate(['1', '100', '10000']):
    A = np.loadtxt(os.path.join(d, f'lim0_p{p0}_N160.dat')); n = 160
    X = A[:, 0].reshape(n, n); Y = A[:, 1].reshape(n, n); P = A[:, 2].reshape(n, n)
    a = ax[0, k]; c = a.contourf(X, Y, P - float(p0), 21, cmap='viridis'); a.contour(X, Y, P - float(p0), 21, colors='k', linewidths=0.3)
    fig.colorbar(c, ax=a); a.set_aspect('equal'); a.set_title(f'$p-p_0$, $p_0={p0}$, $160^2$, t = 5 (IMEX, limiter=0)')
for k, p0 in enumerate(['1', '100', '10000']):
    A = np.loadtxt(os.path.join(d, f'lim0_p{p0}_N160.dat')); n = 160
    j = n // 2; X = A[:, 0].reshape(n, n); P = A[:, 2].reshape(n, n); y = A[:, 1].reshape(n, n)[j, 0]
    xe = np.linspace(-5, 5, 801)
    ax[1, 0].plot(xe, exact_p(xe, y, float(p0)) - float(p0), 'k-', lw=0.8, label='exact' if k == 0 else None)
    ax[1, 0].plot(X[j], P[j] - float(p0), 'o', ms=2.5, label=f'$p_0={p0}$')
ax[1, 0].set_title(f'cut $y={y:.3f}$: $p-p_0$'); ax[1, 0].legend(); ax[1, 0].set_xlabel('x')
paper = {'a2r2': [4.2250e-01, 1.6787e-01, 4.6463e-02, 1.1957e-02], 'a2r2u': [1.0690, 3.4438e-01, 8.7884e-02, 2.5092e-02], 'a2': [1.6308e-01, 6.3713e-02, 1.6970e-02, 4.2119e-03]}
lab = {'a2r2': r'$\alpha_2\rho_2$', 'a2r2u': r'$\alpha_2\rho_2u$', 'a2': r'$\alpha_2$'}
Ns = np.array([20, 40, 80, 160])
def errs(tag, p0, c):
    e = []
    for N in Ns:
        s = open(os.path.join(d, f'{tag}_p{p0}_N{N}.log')).read(); e.append(float(re.search(rf' {c} ([0-9.eE+-]+)', re.search(r'^L2 .*$', s, re.M).group()).group(1)))
    return np.array(e)
for k, c in enumerate(['a2r2', 'a2r2u', 'a2']):
    a = ax[1, 1] if k < 2 else ax[1, 2]
    a.loglog(Ns, errs('lim0', '1', c), 'o-', label=f'{lab[c]} IMEX limiter=0')
    a.loglog(Ns, errs('mm', '1', c), 's--', label=f'{lab[c]} IMEX minmod')
    a.loglog(Ns, paper[c], 'x:', label=f'{lab[c]} Battisti-Boscheri SIFV-BN')
for a in ax[1, 1:]:
    a.loglog(Ns, 2.0 * (Ns / 20.0)**-2, 'k-', lw=0.6, label='slope 2'); a.set_xlabel('N'); a.set_ylabel('L2 error ($p_0=1$)'); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(out, dpi=110)
