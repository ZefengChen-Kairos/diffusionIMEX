#!/usr/bin/env python3
"""x=0 and y=0 line-outs of several runs (RK_*.dat)."""
import sys, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_compare import load_dir
args = sys.argv[1:]; png = args[-1]; dirs = args[:-1][0::2]; labels = args[:-1][1::2]
data = [load_dir(d) for d in dirs]
fields = [('rho', r'$\rho$', False), ('p', r'$p$', True), ('u2', r'$|u|^2$', False), ('m', r'$m$', False), ('vy', r'$v_y$', False), ('Bx', r'$B_x$', False)]
fig, axs = plt.subplots(2, len(fields), figsize=(4.2*len(fields), 8))
for c, (f, lab, logy) in enumerate(fields):
    for d, l in zip(data, labels):
        jmid = len(d['y'])//2; imid = len(d['x'])//2
        axs[0, c].plot(d['y'], 0.5*(d[f][imid, :] + d[f][imid-1, :]), '-', lw=1, label=l)
        axs[1, c].plot(d['x'], 0.5*(d[f][:, jmid] + d[f][:, jmid-1]), '-', lw=1, label=l)
    axs[0, c].set_title('%s along x=0 (vs y)' % lab); axs[1, c].set_title('%s along y=0 (vs x)' % lab)
    for r in range(2):
        axs[r, c].grid(alpha=.3); axs[r, c].legend(fontsize=7)
        if logy: axs[r, c].set_yscale('log')
plt.tight_layout(); plt.savefig(png, dpi=100); print('wrote', png)
