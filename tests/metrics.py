"""Metrics against exact solutions (exact_sg: two-SG Euler; exact5: 5-eq./Kapila model with SPB (II.5) shock relations).
Output columns of spb6:
0 x 1 y 2 p 3 u 4 v 5 rho 6 alpha2 7 Y2 8 rho1 9 rho2 10 cf 11 cw 12 T1 13 T2 14 p1 15 p2 16 g1-g2
usage: python3 metrics.py CASE file [file ...]"""
import numpy as np, sys, os
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import exact_sg as esg, exact5 as e5
def L1(a, b, m=None):
    if m is None: return np.mean(np.abs(a - b))
    return np.mean(np.abs(a - b)[m])
case = int(sys.argv[1])
for f in sys.argv[2:]:
    d = np.loadtxt(f); x = d[:, 0]; tag = os.path.basename(f)
    if case in (5, 6):   # SPB 4.3.1 / 4.3.2: nearly pure phases -> two-SG Euler exact solution
        pl, ra, t = (1e9, 1.0, 240e-6) if case == 5 else (1e12, 10.0, 8.3e-6)
        r, u, p, ps, us = esg.solve(1000, 0, pl, 4.4, 6e8, ra, 0, 1e5, 1.4, 0, 0.75, t, x)
        g = (x > 0.75 + us * t + 2e-3) & (r < 100)            # shocked gas (between contact and shock)
        xs = x[np.where(d[:, 3] > 0.5 * us)[0].max()]          # numerical shock: last cell with u > u*/2
        xse = x[np.where(u > 0.5 * us)[0].max()] if (u > 0.5 * us).any() else np.nan
        print('%-24s L1 u %8.2f rho %7.3f p/p* %.4f | max u/u* %.4f | gas L1 p/p* %s | shock x %.4f (exact %.4f) | min p %.3e' % (
            tag, L1(d[:, 3], u), L1(d[:, 5], r), L1(d[:, 2], p) / ps, np.nanmax(d[:, 3]) / us,
            ('%.3f' % (L1(d[:, 2], p, g) / ps)) if g.any() else ' n/a ', xs, xse, np.nanmin(d[:, 2])))
    elif case in (12, 13):   # SPB 5.1 / 5.2 epoxy-spinel: 5-eq. exact solution
        a = [0.5954, 0.4046]; rr = [1185, 3622]; ES = ([2.43, 1.62], [5.3e9, 141e9]); pl, t = (1e10, 80e-6) if case == 12 else (2e11, 34e-6)
        o, ps, us = e5.solve(e5.Side(a, rr, 0, pl, *ES), e5.Side(a, rr, 0, 1e5, *ES), 0.6, t, x)
        m = x > 0.65
        print('%-24s L1 u %7.2f rho %7.2f p/p* %.4f | max p/p* %.3f | rho* exact %.1f num(max) %.1f | alpha1 L1 %.2e' % (
            tag, L1(d[:, 3], o['u']), L1(d[:, 5], o['rho']), L1(d[:, 2], o['p']) / ps, d[m, 2].max() / ps, o['rho'][m].max(), d[m, 5].max(), L1(1 - d[:, 6], o['a1'])))
    elif case == 0:          # SPB 4.5 cavitation: 5-eq. exact solution has a vacuum (escape velocity < 100 m/s)
        Ls = e5.Side([0.99, 0.01], [1000., 1.], -100., 1e5, [4.4, 1.4], [6e8, 0.]); Rs = e5.Side([0.99, 0.01], [1000., 1.], 100., 1e5, [4.4, 1.4], [6e8, 0.])
        o, ps, us = e5.solve(Ls, Rs, 0.5, 1.85e-3, x)
        xv = 0.5 + (Ls.u - Ls.f(1e-9)) * 1.85e-3            # left vacuum edge
        m = (x < xv - 0.005) | (x > 1 - xv + 0.005)             # outside the vacuum
        print('%-24s fan: L1 u %.4f p %.3e rho %.3f | vacuum edge %.4f: centre p %.3e rho %.3e alpha_g %.6f' % (
            tag, L1(d[:, 3], o['u'], m), L1(d[:, 2], o['p'], m), L1(d[:, 5], o['rho'], m), xv, d[len(x)//2, 2], d[len(x)//2, 5], d[len(x)//2, 6]))
    elif case in (3, 31):    # dodecane liquid-vapour tube; RELAX 0 (p only) vs 5-eq. exact solution
        G = [2.35, 1.025]; PI = [4e8, 0.]; B = [0., 0.]
        if case == 31: G = [1.206, 1.021]; PI = [1681e5, 0.]; B = [7.51e-4, 0.]   # NASG: liquid covolume
        o, ps, us = e5.solve(e5.Side([1 - 1e-8, 1e-8], [500., 2.], 0., 1e8, G, PI, B), e5.Side([1e-8, 1 - 1e-8], [500., 2.], 0., 1e5, G, PI, B), 0.75, 473e-6, x)
        w = (0.2, 0.8) if case == 3 else (0.5, 0.85)          # plateau between the rarefaction tail and the contact
        m = (x > w[0]) & (x < w[1])
        print('%-24s (p-relax exact) L1 p/1e5 %.3f u %.3f rho %.3f | plateau[%.2f,%.2f] L1 p/p* %.4f max|p-p*|/p* %.3f | min p %.3e' % (
            tag, L1(d[:, 2], o['p']) / 1e5, L1(d[:, 3], o['u']), L1(d[:, 5], o['rho']), w[0], w[1], L1(d[:, 2], ps, m) / ps, np.max(np.abs(d[m, 2] - ps)) / ps, d[:, 2].min()))
