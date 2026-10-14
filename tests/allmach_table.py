"""All-Mach tests, L2 norms only.  Taylor-Green (cases 57/58, run/allmach): Mach sweep on 50^2 and convergence N = 32..256;
low-Mach two-phase vortex (case 56, run/allmach/vortex): convergence N = 32..256.  Order = log2(e_{N/2} / e_N)."""
import re, glob, os, sys, math
root = sys.argv[1] if len(sys.argv) > 1 else 'run/allmach'
num = r'([-+0-9.eE]+)'
def tg(f):
    s = open(f).read()
    a = re.search(r'TG face-div L1 \S+ L2 ' + num, s)
    b = re.search(r'TG p0 \S+ M ' + num + r' N (\d+) \| div L1 \S+ L2 ' + num + r'.*?u-err L1 \S+ L2 ' + num + r'.*?p2 L2 ' + num, s)
    c = re.search(r'rejects (\d+).*?\((\S+)/stage\)', s)
    if not (a and b and c): return None
    return dict(M=float(b[1]), fd=float(a[1]), wd=float(b[3]), ue=float(b[4]), p2=float(b[5]), rej=int(c[1]), nw=c[2])
def order(e, e0): return '%.2f' % math.log2(e0 / e) if e0 and e > 0 else '—'
TG = {}
for f in glob.glob(os.path.join(root, 'tg*_N*.log')):
    m = re.search(r'tg(\d+)_p(\S+)_lim(\d)_N(\d+)\.log', f); r = tg(f)
    if r: TG[(int(m[1]), int(m[3]), m[2], int(m[4]))] = r
print('### Taylor-Green (BB 4.4), t = 0.2, Mach sweep on 50^2 (dt = 2.5e-3)')
print('| case | lim | p0 | M | face-div | wide-central div | u-err | p - p0 | rejects | Newton/stage |\n|---|---|---|---|---|---|---|---|---|---|')
for k in sorted(k for k in TG if k[3] == 50):
    r = TG[k]; c, lim, ps, N = k
    print('| %d | %d | %s | %.2e | %.2e | %.2e | %.2e | %.2e | %d | %s |' % (c, lim, ps, r['M'], r['fd'], r['wd'], r['ue'], r['p2'], r['rej'], r['nw']))
print('\n### Taylor-Green convergence (dt = 0.08/N)')
print('| case | lim | p0 | N | face-div | wide-central div | (order) | u-err | (order) | p - p0 | (order) | rejects |\n|---|---|---|---|---|---|---|---|---|---|---|---|')
for k in sorted((k for k in TG if k[3] != 50), key=lambda k: (k[0], k[1], float(k[2]), k[3])):
    r = TG[k]; c, lim, ps, N = k; r0 = TG.get((c, lim, ps, N // 2)) if N > 32 else None
    o = lambda q: order(r[q], r0[q]) if r0 else '—'
    print('| %d | %d | %s | %d | %.2e | %.2e | %s | %.2e | %s | %.2e | %s | %d |' % (c, lim, ps, N, r['fd'], r['wd'], o('wd'), r['ue'], o('ue'), r['p2'], o('p2'), r['rej']))
print('\n### two-phase isentropic vortex (case 56), L2 errors vs exact (initial) solution, t = 5')
print('| init | p0 | N | a2r2 | (order) | a2r2u | (order) | p | (order) | u | (order) | steps | rejects | Newton/stage |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
V = {}
for f in glob.glob(os.path.join(root, 'vortex', '*.log')):
    s = open(f).read(); m = re.search(r'(\w+?)_p(\S+)_N(\d+)\.log', os.path.basename(f))
    e = re.search(r'a2r2 ' + num + ' a2r2u ' + num + ' a2 ' + num + ' p ' + num + ' u ' + num, s); c = re.search(r'steps (\d+).*?rejects (\d+).*?\((\S+)/stage\)', s)
    if e and c: V[(m[1], m[2], int(m[3]))] = ([float(e[i]) for i in (1, 2, 4, 5)], c)
for k in sorted(V, key=lambda k: (k[0], float(k[1]), k[2])):
    (err, c), prev = V[k], V.get((k[0], k[1], k[2] // 2))
    cols = ' | '.join('%.3e | %s' % (err[i], order(err[i], prev[0][i]) if prev else '—') for i in range(4))
    print('| %s | %s | %d | %s | %s | %s | %s |' % (k[0], k[1], k[2], cols, c[1], c[2], c[3]))
