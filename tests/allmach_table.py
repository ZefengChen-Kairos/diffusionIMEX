"""L2 norms only. Tabulate the all-Mach tests: Taylor-Green (cases 57/58, run/allmach) and low-Mach two-phase vortex (case 56, run/allmach/vortex)."""
import re, glob, os, sys
root = sys.argv[1] if len(sys.argv) > 1 else 'run/allmach'
num = r'([-+0-9.eE]+)'
def tg(f):
    s = open(f).read()
    a = re.search(r'TG face-div L1 \S+ L2 ' + num, s); b = re.search(r'TG p0 \S+ M ' + num + r' N (\d+) \| div L1 \S+ L2 ' + num + r'.*?u-err L1 \S+ L2 ' + num + r'.*?p2 L2 ' + num, s)
    c = re.search(r'rejects (\d+).*?\((\S+)/stage\)', s)
    if not (a and b and c): return None
    return dict(fd=float(a[1]), M=float(b[1]), N=int(b[2]), wd=float(b[3]), ue=float(b[4]), p2=float(b[5]), rej=int(c[1]), nw=c[2])
print('### Taylor-Green (BB 4.4), t = 0.2')
print('| case | lim | N | p0 | M | face-div L2 | wide-central div L2 | u-err L2 | p2 L2 | rejects | Newton/stage |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
rows = []
for f in glob.glob(os.path.join(root, 'tg*_N*.log')):
    m = re.search(r'tg(\d+)_p(\S+)_lim(\d)_N(\d+)\.log', f); r = tg(f)
    if r: rows.append((int(m[1]), int(m[3]), int(m[4]), float(m[2]), m[2], r))
for c, lim, N, p0, ps, r in sorted(rows):
    print('| %d | %d | %d | %s | %.2e | %.2e | %.2e | %.2e | %.2e | %d | %s |' % (c, lim, N, ps, r['M'], r['fd'], r['wd'], r['ue'], r['p2'], r['rej'], r['nw']))
print('\n### two-phase isentropic vortex (case 56), L2 errors')
print('| init | N | p0 | a2r2 | a2r2u | p | u | steps | rejects | Newton/stage |')
print('|---|---|---|---|---|---|---|---|---|---|')
rows = []
for f in glob.glob(os.path.join(root, 'vortex', '*.log')):
    s = open(f).read(); m = re.search(r'(\w+?)_p(\S+)_N(\d+)\.log', os.path.basename(f))
    e = re.search(r'a2r2 ' + num + ' a2r2u ' + num + ' a2 ' + num + ' p ' + num + ' u ' + num, s); c = re.search(r'steps (\d+).*?rejects (\d+).*?\((\S+)/stage\)', s)
    if e and c: rows.append((m[1], int(m[3]), float(m[2]), m[2], e, c))
for t, N, p0, ps, e, c in sorted(rows):
    print('| %s | %d | %s | %.3e | %.3e | %.3e | %.3e | %s | %s | %s |' % (t, N, ps, float(e[1]), float(e[2]), float(e[4]), float(e[5]), c[1], c[2], c[3]))
