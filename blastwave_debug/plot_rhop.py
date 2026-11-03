#!/usr/bin/env python3
"""Density / pressure oscillations: zoomed maps with TIGHT colour scales (interior and wings) and roughness metrics.
   rough(f, region) = ||Lap_h f||_2 / ||grad_h f||_2 * h ;  usage: plot_rhop.py run1 label1 run2 label2 ... out.png"""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk
from plot_lines2 import load_ref
def load(p):
    if p.endswith('.dat'):
        d = load_ref(p); return d['x'], d['rho'].T, d['p'].T
    vt = sorted(glob.glob(os.path.join(p, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))[-1]
    v = read_vtk(vt); nx = v['nx']; x = -0.5 + (np.arange(nx)+0.5)/nx
    return x, v['rho'], v['p']
def rough(f, mask):
    lap = f[2:,1:-1] + f[:-2,1:-1] + f[1:-1,2:] + f[1:-1,:-2] - 4*f[1:-1,1:-1]
    gx = 0.5*(f[2:,1:-1] - f[:-2,1:-1]); gy = 0.5*(f[1:-1,2:] - f[1:-1,:-2]); m = mask[1:-1,1:-1]
    return np.sqrt(np.sum((lap*m)**2)) / np.sqrt(np.sum(((gx*gx+gy*gy)*m)) + 1e-300)
def tv1d(line):   # total variation / range of a line-out: 1 for monotone / single-peaked smooth profiles
    return np.sum(np.abs(np.diff(line))) / (line.max() - line.min() + 1e-300)
args = sys.argv[1:]; png = args[-1]; paths = args[:-1][0::2]; labels = args[:-1][1::2]
regions = {'interior': lambda X, Y: (np.abs(X) < 0.1) & (np.abs(Y) < 0.13), 'wings': lambda X, Y: (np.abs(Y) > 0.17) & (np.abs(Y) < 0.33) & (np.abs(X) < 0.2), 'lobes(x)': lambda X, Y: (np.abs(X) > 0.15) & (np.abs(X) < 0.28) & (np.abs(Y) < 0.08)}
print('%-24s %8s %8s %8s %8s %8s %8s   %8s %8s' % ('run', 'rho_int', 'p_int', 'rho_wing', 'p_wing', 'rho_lobe', 'p_lobe', 'TVrho_x0', 'TVp_x0'))
fig, axs = plt.subplots(4, len(paths), figsize=(4.3*len(paths), 16.5), squeeze=False)
for c, (p, l) in enumerate(zip(paths, labels)):
    x, rho, pr = load(p); X, Y = np.meshgrid(x, x); n = len(x); imid = n//2
    r = {k: (rough(rho, m(X, Y)), rough(pr, m(X, Y))) for k, m in regions.items()}
    iy = np.abs(x) < 0.45; lr = 0.5*(rho[:, imid] + rho[:, imid-1])[iy]; lp = 0.5*(pr[:, imid] + pr[:, imid-1])[iy]
    print('%-24s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f   %8.2f %8.2f' % (l, r['interior'][0], r['interior'][1], r['wings'][0], r['wings'][1], r['lobes(x)'][0], r['lobes(x)'][1], tv1d(lr), tv1d(lp)))
    panels = [('rho interior', rho, 0.17, 0.30, 0.16), ('p interior', pr, 95, 165, 0.16), ('rho wings', rho, 0.9, 1.5, 0.45), ('p wings', pr, 0.0, 12.0, 0.45)]
    for rr, (name, f, vmin, vmax, zoom) in enumerate(panels):
        ax = axs[rr, c]; im = ax.pcolormesh(x, x, f, shading='auto', cmap='jet', vmin=vmin, vmax=vmax); plt.colorbar(im, ax=ax)
        ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom); ax.set_aspect('equal'); ax.set_title('%s %s' % (l, name), fontsize=9)
plt.tight_layout(); plt.savefig(png, dpi=90); print('wrote', png)
