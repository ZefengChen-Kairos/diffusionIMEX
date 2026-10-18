"""All-Mach tests, L2 errors of mixture density rho, velocity vector |u - u_ex| and mixture pressure p (L2 = sqrt(sum e^2 dx dy)).
Reference = exact cell averages written by the solver itself (init=1; for both tests the exact solution at t_end equals the initial state:
Taylor-Green is steady, the vortex returns after one period).  Order = log2(e_{N/2} / e_N).
  Taylor-Green (cases 57/58, <root>/tg<case>_p<p0>_lim<l>_N<N>.log/.dat): Mach sweep on 50^2 (+ divergences) and convergence N = 32..256
  two-phase vortex (case 56, <root>/vortex/lim0_p<p0>_N<N>.log/.dat): BB Table 1, p0 = 1 ... 1e4, convergence N = 32..256
usage: python3 tests/allmach_table.py [root=run/v21e/allmach]"""
import re, glob, os, sys, math, subprocess
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
root = sys.argv[1] if len(sys.argv) > 1 else 'run/v21e/allmach'
BIN = os.path.join(H, 'src', 'spb6_imex')   # tests/run_archive.sh builds it; only used with init=1 (exact initial cell averages)
if not os.path.exists(BIN): subprocess.run(['g++', '-O2', '-std=c++17', '-fopenmp', '-I/usr/include/eigen3', os.path.join(H, 'src', 'SPB6_TV_IMEX.cpp'), '-o', BIN], check=True)
num = r'([-+0-9.eE]+)'

def load(f): d = np.loadtxt(f); return d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]   # x y p u v rho
def ref(case, N, p0):
    os.makedirs(os.path.join(root, 'init'), exist_ok=True); f = os.path.join(root, 'init', 'init_c%d_N%d_p%s.dat' % (case, N, p0))
    if not os.path.exists(f): subprocess.run([BIN, 'case=%d' % case, 'N=%d' % N, 'ny=%d' % N, 'p0=%s' % p0, 'init=1', 'out=' + os.path.abspath(f)], check=True, stdout=subprocess.DEVNULL)
    return load(f)
def errors(dat, case, N, p0):
    x, y, p, u, v, r = load(dat); _, _, p0_, u0, v0, r0 = ref(case, N, p0)
    A = (x.max() - x.min()) / (N - 1) * (y.max() - y.min()) / (N - 1)
    L2 = lambda e: math.sqrt(float(np.sum(e * e)) * A)
    return L2(r - r0), L2(np.hypot(u - u0, v - v0)), L2(p - p0_)
def order(e, e0): return '%.2f' % math.log2(e0 / e) if e0 and e > 0 else '—'
def runinfo(s):
    c = re.search(r'steps (\d+).*?rejects (\d+).*?\((\S+)/stage\)', s); return (c[1], c[2], c[3]) if c else None

TG = {}
for f in glob.glob(os.path.join(root, 'tg*_N*.log')):
    m = re.search(r'tg(\d+)_p(\S+)_lim(\d)_N(\d+)\.log', f); s = open(f).read(); ri = runinfo(s)
    if not ri or not os.path.exists(f[:-4] + '.dat'): continue
    a = re.search(r'TG face-div L1 \S+ L2 ' + num, s); b = re.search(r'TG p0 \S+ M ' + num + r' N \d+ \| div L1 \S+ L2 ' + num, s)
    case, N = int(m[1]), int(m[4])
    TG[(case, int(m[3]), m[2], N)] = dict(M=float(b[1]), fd=float(a[1]), wd=float(b[2]), e=errors(f[:-4] + '.dat', case, N, m[2]), ri=ri)
if TG:
    print('### Taylor-Green (BB 4.4), t = 0.2, Mach sweep on 50^2 (dt = 2.5e-3), L2')
    print('| case | lim | p0 | M | rho | velocity | p | face-div | wide-central div | rejects | Newton/stage |\n|---|---|---|---|---|---|---|---|---|---|---|')
    for k in sorted((k for k in TG if k[3] == 50), key=lambda k: (k[0], k[1], float(k[2]))):
        r = TG[k]; e = r['e']
        print('| %d | %d | %s | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %s | %s |' % (k[0], k[1], k[2], r['M'], e[0], e[1], e[2], r['fd'], r['wd'], r['ri'][1], r['ri'][2]))
    print('\n### Taylor-Green convergence (dt = 0.08/N), L2')
    print('| case | lim | p0 | N | rho | (order) | velocity | (order) | p | (order) | face-div | rejects |\n|---|---|---|---|---|---|---|---|---|---|---|---|')
    for k in sorted((k for k in TG if k[3] != 50), key=lambda k: (k[0], k[1], float(k[2]), k[3])):
        r = TG[k]; r0 = TG.get((k[0], k[1], k[2], k[3] // 2)); e = r['e']
        cols = ' | '.join('%.2e | %s' % (e[i], order(e[i], r0['e'][i]) if r0 else '—') for i in range(3))
        print('| %d | %d | %s | %d | %s | %.2e | %s |' % (k[0], k[1], k[2], k[3], cols, r['fd'], r['ri'][1]))

V = {}
for f in glob.glob(os.path.join(root, 'vortex', '*.log')):
    m = re.search(r'(\w+?)_p(\S+)_N(\d+)\.log', os.path.basename(f)); ri = runinfo(open(f).read())
    if m and ri and os.path.exists(f[:-4] + '.dat'): V[(m[1], m[2], int(m[3]))] = (errors(f[:-4] + '.dat', 56, int(m[3]), m[2]), ri)
if V:
    print('\n### two-phase isentropic vortex (case 56), t = 5 (one period), L2 vs exact cell averages')
    print('| init | p0 | N | rho | (order) | velocity | (order) | p | (order) | steps | rejects | Newton/stage |\n|---|---|---|---|---|---|---|---|---|---|---|---|')
    for k in sorted(V, key=lambda k: (k[0], float(k[1]), k[2])):
        (e, ri), prev = V[k], V.get((k[0], k[1], k[2] // 2))
        cols = ' | '.join('%.3e | %s' % (e[i], order(e[i], prev[0][i]) if prev else '—') for i in range(3))
        print('| %s | %s | %d | %s | %s | %s | %s |' % (k[0], k[1], k[2], cols, *ri))
