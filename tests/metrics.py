"""Metrics against exact / reference solutions.  Output columns of spb6:
0 x 1 y 2 p 3 u 4 v 5 rho 6 alpha2 7 Y2 8 rho1 9 rho2 10 cf 11 cw 12 T1 13 T2 14 p1 15 p2 16 g1-g2
usage: python3 metrics.py CASE file [file ...]"""
import numpy as np, sys, os
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import exact_sg as esg, exact5 as e5
REF = H + '/ref/'
def L(f): return np.loadtxt(f)
def ref_oracle(name):   # oracle explicit N=4000 averaged to 1000 cells: x p u rho alpha_g Y_g
    return np.loadtxt(REF + name).reshape(-1, 4, 6).mean(1)
case = int(sys.argv[1])
for f in sys.argv[2:]:
    d = L(f); x = d[:, 0]; tag = os.path.basename(f)
    if case in (5, 6):
        pl, ra, t = (1e9, 1.0, 240e-6) if case == 5 else (1e12, 10.0, 8.3e-6)
        r, u, p, ps, us = esg.solve(1000, 0, pl, 4.4, 6e8, ra, 0, 1e5, 1.4, 0, 0.75, t, x)
        g = x > 0.75 + us * t + 0.01
        print('%-26s L1 u %7.2f rho %7.3f | gas max u %6.0f (exact %.0f) | gas L1 p/p* %.3f | min p %.3e (p* %.3e)' % (tag, np.mean(abs(d[:, 3] - u)), np.mean(abs(d[:, 5] - r)), d[g, 3].max(), us, np.mean(abs(d[:, 2] - p)[g]) / ps, d[:, 2].min(), ps))
    elif case in (12, 13):
        a = [0.5954, 0.4046]; rr = [1185, 3622]; ES = ([2.43, 1.62], [5.3e9, 141e9]); pl, t = (1e10, 80e-6) if case == 12 else (2e11, 34e-6)
        o, ps, us = e5.solve(e5.Side(a, rr, 0, pl, *ES), e5.Side(a, rr, 0, 1e5, *ES), 0.6, t, x)
        print('%-26s L1 u %7.2f rho %7.2f p/p* %.4f | max p/p* %.3f | rho* exact %.1f num(max) %.1f' % (tag, np.mean(abs(d[:, 3] - o['u'])), np.mean(abs(d[:, 5] - o['rho'])), np.mean(abs(d[:, 2] - o['p'])) / ps, d[:, 2][x > 0.65].max() / ps, o['rho'][x > 0.65].max(), d[:, 5][x > 0.65].max()))
    elif case == 0:
        R = ref_oracle('ref0_0.dat')
        print('%-26s L1 u %.3f alpha %.3e p %.3e | centre p %.4e (ref %.4e) alpha %.4f (ref %.4f) | min p %.3e' % (tag, np.mean(abs(d[:, 3] - R[:, 2])), np.mean(abs(d[:, 6] - R[:, 4])), np.mean(abs(d[:, 2] - R[:, 1])), d[len(x)//2, 2], R[len(x)//2, 1], d[len(x)//2, 6], R[len(x)//2, 4], d[:, 2].min()))
    elif case == 3:
        for k in (0, 2, 3):
            pass
        R0 = ref_oracle('ref3_0.dat')
        print('%-26s vs ref3_0: L1 u %.3f p/1e5 %.3f Y %.2e | min p[0.1,0.4] %.3e (ref %.3e) | max p in gas %.3e' % (tag, np.mean(abs(d[:, 3] - R0[:, 2])), np.mean(abs(d[:, 2] - R0[:, 1])) / 1e5, np.mean(abs(d[:, 7] - R0[:, 5])), d[100:400, 2].min(), R0[100:400, 1].min(), d[x > 0.76, 2].max()))
