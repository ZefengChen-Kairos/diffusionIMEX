"""Battisti-Boscheri Table-1 style convergence tables for the two-phase isentropic vortex (case 56), t = 5, L2 = sqrt(sum e^2 dx dy)
against the exact cell averages (written by the solver with init=1).  Variables: alpha2 rho2, alpha2 rho2 u, alpha2 (as BB Table 1),
mixture density rho, velocity vector |u - u_ex|, mixture pressure p.  One table per initial data (lim0 = cell averages as BB,
wp = discretely well-prepared velocity), blocks p0 = 1 ... 1e10, N = 32 ... 256, h = 10/N.
usage: python3 tests/vortex_bb_table.py OUT_PREFIX     (writes OUT_PREFIX_{lim0,wp}.png/.tex/.md)"""
import os, sys, math, glob, subprocess
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(H, 'src', 'spb6_imex')   # tests/run_archive.sh builds it; only used with init=1 (exact initial cell averages)
if not os.path.exists(BIN): subprocess.run(['g++', '-O2', '-std=c++17', '-fopenmp', '-I/usr/include/eigen3', os.path.join(H, 'src', 'SPB6_TV_IMEX.cpp'), '-o', BIN], check=True); REF = os.path.join(H, 'run', 'vortex_ref'); os.makedirs(REF, exist_ok=True)
DIRS = [os.path.join(H, 'run', 'v21e', 'vortex'), os.path.join(H, 'run', 'v21e', 'allmach', 'vortex')]
P0 = ['1', '10', '100', '1000', '10000', '1e6', '1e8', '1e10']; NS = [32, 64, 128, 256]
def load(f): d = np.loadtxt(f); return d
def ref(N, p0):
    f = os.path.join(REF, 'init_N%d_p%s.dat' % (N, p0))
    if not os.path.exists(f): subprocess.run([BIN, 'case=56', 'N=%d' % N, 'ny=%d' % N, 'p0=%s' % p0, 'init=1', 'out=' + f], check=True, stdout=subprocess.DEVNULL)
    return load(f)
def find(tag, p0, N):
    for d in DIRS:
        for name in ['%s_p%s_N%d.dat' % (tag, p0, N), '%s_p%g_N%d.dat' % (tag, float(p0), N)]:
            f = os.path.join(d, name)
            if os.path.exists(f) and os.path.exists(f[:-4] + '.log') and '\ncase' in '\n' + open(f[:-4] + '.log').read(): return f
    return None
def errors(f, N, p0):
    d, e = load(f), ref(N, p0); h = 10.0 / N
    # columns: x y p u v rho alpha2 Y2 rho1 rho2 ...
    a2r2 = lambda q: q[:, 6] * q[:, 9]; L2 = lambda v: math.sqrt(float(np.sum(v * v)) * h * h)
    return [L2(a2r2(d) - a2r2(e)), L2(a2r2(d) * d[:, 3] - a2r2(e) * e[:, 3]), L2(d[:, 6] - e[:, 6]),
            L2(d[:, 5] - e[:, 5]), L2(np.hypot(d[:, 3] - e[:, 3], d[:, 4] - e[:, 4])), L2(d[:, 2] - e[:, 2])]
HEAD = [r'$N$', r'$h$', r'$L_2(\alpha_2\rho_2)$', r'$\mathcal{O}$', r'$L_2(\alpha_2\rho_2u)$', r'$\mathcal{O}$', r'$L_2(\alpha_2)$', r'$\mathcal{O}$',
        r'$L_2(\rho)$', r'$\mathcal{O}$', r'$L_2(\mathbf{u})$', r'$\mathcal{O}$', r'$L_2(p)$', r'$\mathcal{O}$']
TITLE = {'lim0': 'cell-average initial data (as Battisti-Boscheri)', 'wp': 'discretely well-prepared initial velocity'}
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, 'docs', 'figures', 'vortex_table')
for tag in ['lim0', 'wp']:
    blocks = []
    for p0 in P0:
        rows, prev = [], None
        for N in NS:
            f = find(tag, p0, N)
            if not f: continue
            e = errors(f, N, p0); r = [str(N), '%.4g' % (10.0 / N)]
            for k in range(6): r += ['%.4E' % e[k], '%.2f' % math.log2(prev[k] / e[k]) if prev else '–']
            rows.append(r); prev = e
        if rows: blocks.append((p0, rows))
    # markdown + latex
    md = ['| ' + ' | '.join(h.replace('$', '').replace('\\mathcal{O}', 'O').replace('\\', '') for h in HEAD) + ' |', '|' + '---|' * len(HEAD)]
    tex = [r'\begin{tabular}{rr' + 'rr' * 6 + '}', r'\hline', ' & '.join(HEAD) + r' \\', r'\hline']
    for p0, rows in blocks:
        lab = '$p_0=%s$' % ('10^{%d}' % round(math.log10(float(p0))) if float(p0) > 1 else '1')
        md.append('| ' + lab.replace('$', '') + ' |' + ' |' * (len(HEAD) - 1)); tex.append(r'\multicolumn{%d}{c}{%s} \\' % (len(HEAD), lab))
        for r in rows: md.append('| ' + ' | '.join(r) + ' |'); tex.append(' & '.join(r).replace('–', '--') + r' \\')
    tex += [r'\hline', r'\end{tabular}']
    open(out + '_' + tag + '.md', 'w').write('\n'.join(md) + '\n'); open(out + '_' + tag + '.tex', 'w').write('\n'.join(tex) + '\n')
    # figure
    nrow = sum(len(r) + 1 for _, r in blocks) + 1
    fig, ax = plt.subplots(figsize=(15, 0.27 * nrow + 1.0)); ax.axis('off')
    y = 1.0; dy = 1.0 / (nrow + 1); xs = np.cumsum([0, 0.045, 0.06, 0.095, 0.045, 0.095, 0.045, 0.095, 0.045, 0.095, 0.045, 0.095, 0.045, 0.095, 0.045])
    xs = xs / xs[-1] * 0.98 + 0.01
    for k, hd in enumerate(HEAD): ax.text(xs[k + 1] - 0.005, y, hd, ha='right', va='center', fontsize=9)
    ax.plot([0.01, 0.99], [y - dy / 2] * 2, 'k-', lw=0.8); y -= dy
    for p0, rows in blocks:
        lab = r'$p_0=%s$' % ('10^{%d}' % round(math.log10(float(p0))) if float(p0) > 1 else '1')
        ax.text(0.5, y, lab, ha='center', va='center', fontsize=9, style='italic'); y -= dy
        for r in rows:
            for k, c in enumerate(r): ax.text(xs[k + 1] - 0.005, y, c, ha='right', va='center', fontsize=8.5, family='monospace')
            y -= dy
    ax.plot([0.01, 0.99], [y + dy / 2] * 2, 'k-', lw=0.8)
    ax.set_ylim(y, 1.0 + dy); ax.set_xlim(0, 1)
    ax.set_title('Two-phase isentropic vortex (Battisti-Boscheri 4.1), t = 5, L2 errors vs exact cell averages, IMEX (V16), %s' % TITLE[tag], fontsize=10)
    fig.tight_layout(); fig.savefig(out + '_' + tag + '.png', dpi=130); plt.close(fig)
    print('wrote', out + '_' + tag)
