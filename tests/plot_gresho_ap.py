"""Gresho vortex low-Mach AP check (B7 Theorem 3, devlog V11 §11.8): table + figure.
usage: python3 tests/plot_gresho_ap.py OUTDIR [FIG.png]
OUTDIR holds im_50..im_55 (IMEX, M = 1e-1..1e-6, same fixed dt) and optionally ex_50, ex_51 (explicit HLLC)."""
import os, re, sys, numpy as np

MS = {50: 1e-1, 51: 1e-2, 52: 1e-3, 53: 1e-4, 54: 1e-5, 55: 1e-6}


def vel(r): return np.where(r < 0.2, 5 * r, np.where(r < 0.4, 2 - 5 * r, 0.0))


def pex(r):   # incompressible p - p_infinity (rho = 1)
    p = np.where(r < 0.2, 12.5 * r * r, np.where(r < 0.4, 12.5 * r * r + 4 * (1 - 5 * r - np.log(0.2) + np.log(np.maximum(r, 1e-300))), -2 + 4 * np.log(2)))
    return p - (-2 + 4 * np.log(2))


def load(f):
    d = np.loadtxt(f); o = np.lexsort((d[:, 0], d[:, 1])); return d[o]


def stats(log):
    if not os.path.exists(log): return ''
    s = open(log).read().strip().splitlines()[-1]
    m = re.search(r'steps (\d+).*?rejects (\d+)', s); n = re.search(r'\(([\d.]+)/stage\)', s)
    return '%s/%s' % (m.group(1), m.group(2)) + (' %s' % n.group(1) if n else '') if m else s[:40]


def diag(d):
    x, y, p, u, v, rho = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    r = np.hypot(x, y); ut = (-u * y + v * x) / np.maximum(r, 1e-300)
    ke = np.sum(0.5 * rho * (u * u + v * v)) / np.sum(0.5 * vel(r) ** 2)
    p2 = p - p[r > 0.45].mean()
    return r, ut, p2, ke, np.sqrt(np.mean((ut - vel(r)) ** 2)), np.abs(rho - 1).max(), np.sqrt(np.mean((p2 - pex(r)) ** 2))


out = sys.argv[1]
im = {c: load('%s/im_%d.dat' % (out, c)) for c in MS if os.path.exists('%s/im_%d.dat' % (out, c))}
ex = {c: load('%s/ex_%d.dat' % (out, c)) for c in (50, 51) if os.path.exists('%s/ex_%d.dat' % (out, c))}
print('%-7s %-15s %-8s %-10s %-11s %-10s %-11s %-11s' % ('M', 'steps/rej N/st', 'KE/KE0', 'L2(u_th)', 'max|rho-1|/M^2', 'L2(p2)', '|du| Cauchy', '|dp2| Cauchy'))
cs = sorted(im); D = {c: diag(im[c]) for c in cs}
cau = {}
for c in cs:
    M = MS[c]; r, ut, p2, ke, eu, er, ep = D[c]
    if c + 1 in im:
        a, b = im[c], im[c + 1]
        du = np.sqrt(np.mean((a[:, 3] - b[:, 3]) ** 2 + (a[:, 4] - b[:, 4]) ** 2)); dp = np.sqrt(np.mean((D[c][2] - D[c + 1][2]) ** 2))
        cau[c] = (du, dp); cs_txt = '%-11.2e %-11.2e' % (du, dp)
    else: cs_txt = ''
    print('%-7.0e %-15s %-8.4f %-10.3e %-14.3f %-10.3e %s' % (M, stats('%s/im_%d.log' % (out, c)), ke, eu, er / M ** 2, ep, cs_txt))
for c in sorted(ex):
    r, ut, p2, ke, eu, er, ep = diag(ex[c])
    print('ex %-4.0e %-15s %-8.4f %-10.3e %-14.3f %-10.3e' % (MS[c], stats('%s/ex_%d.log' % (out, c)), ke, eu, er / MS[c] ** 2, ep))
print('(Cauchy: RMS difference between M and M/10 at the same dt and grid; p2 = p - mean(p, r > 0.45))')

if len(sys.argv) > 2:
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10})
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.6)); rr = np.linspace(0, 0.7, 700)
    a1, a2, a3, a4 = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]
    a1.plot(rr, vel(rr), 'k-', lw=1.2, label='incompressible exact'); a2.plot(rr, pex(rr), 'k-', lw=1.2, label='incompressible $p_2$')
    mk = ['o', 's', '^', 'v', 'D', 'x']
    for i, c in enumerate(cs):
        r, ut, p2 = D[c][0], D[c][1], D[c][2]; m = r < 0.7
        a1.plot(r[m][::7], ut[m][::7], mk[i], ms=2.2, mfc='none', alpha=0.6, label='IMEX $M=10^{%d}$' % round(np.log10(MS[c])))
        a2.plot(r[m][::7], p2[m][::7], mk[i], ms=2.2, mfc='none', alpha=0.6, label='IMEX $M=10^{%d}$' % round(np.log10(MS[c])))
    for c in sorted(ex):
        r, ut = diag(ex[c])[:2]; m = r < 0.7
        a1.plot(r[m][::7], ut[m][::7], '+', ms=3, alpha=0.6, label='explicit HLLC $M=10^{%d}$' % round(np.log10(MS[c])))
    a1.set_xlabel('$r$'); a1.set_ylabel(r'$u_\theta$'); a1.set_xlim(0, 0.7); a1.legend(fontsize=7, markerscale=2)
    a2.set_xlabel('$r$'); a2.set_ylabel(r'$p-p_\infty$ (units of $\rho U^2$)'); a2.set_xlim(0, 0.7); a2.legend(fontsize=7, markerscale=2)
    Ms = [MS[c] for c in cs]
    a3.semilogx(Ms, [D[c][3] for c in cs], 'o-', label='TV-IMEX (fixed $\\Delta t=0.0025$)')
    if ex: a3.semilogx([MS[c] for c in sorted(ex)], [diag(ex[c])[3] for c in sorted(ex)], 's--', label='explicit HLLC (acoustic CFL 0.4)')
    a3.set_xlabel('$M$'); a3.set_ylabel('kinetic energy KE/KE$_0$'); a3.set_ylim(0, 1); a3.invert_xaxis(); a3.legend(fontsize=8); a3.grid(alpha=0.3)
    cc = sorted(cau); Mc = [MS[c] for c in cc]
    a4.loglog(Mc, [cau[c][0] for c in cc], 'o-', label=r'$\|u_M-u_{M/10}\|_2$')
    a4.loglog(Mc, [cau[c][1] for c in cc], 's-', label=r'$\|p_{2,M}-p_{2,M/10}\|_2$')
    a4.loglog(Ms, [D[c][5] / MS[c] ** 2 for c in cs], '^-', label=r'$\max|\rho-1|/M^2$')
    a4.set_xlabel('$M$'); a4.invert_xaxis(); a4.legend(fontsize=8); a4.grid(alpha=0.3, which='both')
    fig.suptitle('Gresho vortex, $80^2$, one turn: the same $\\Delta t$ for $M=10^{-1}\\ldots10^{-6}$ (low-Mach AP)')
    fig.tight_layout(); fig.savefig(sys.argv[2], dpi=130)
