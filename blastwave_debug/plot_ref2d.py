#!/usr/bin/env python3
"""2D zoomed maps: reference solver .dat files and code-2 VTK files side by side (p, |u|^2, m, rho)."""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk
from plot_lines2 import load_ref
def load_any2d(p):
    if p.endswith('.dat'):
        d = load_ref(p); return {'x': d['x'], 'y': d['y'], 'p': d['p'].T, 'u2': d['u2'].T, 'm': d['m'].T, 'rho': d['rho'].T}
    vtks = sorted(glob.glob(os.path.join(p, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))
    v = read_vtk(vtks[-1]); nx, ny = v['nx'], v['ny']
    x = -0.5 + (np.arange(nx)+0.5)/nx; y = -0.5 + (np.arange(ny)+0.5)/ny
    return {'x': x, 'y': y, 'p': v['p'], 'u2': v['velX']**2 + v['velY']**2, 'm': 0.5*(v['B_X']**2+v['B_Y']**2+v['B_Z']**2), 'rho': v['rho']}
args = sys.argv[1:]; png = args[-1]; paths = args[:-1][0::2]; labels = args[:-1][1::2]
data = [load_any2d(p) for p in paths]
rows = [('p', 0, 220), ('u2', 0, 250), ('m', 80, 400), ('rho', 0.1, 3.6)]
zoom = float(os.environ.get('ZOOM', '0.3'))
fig, axs = plt.subplots(len(rows), len(data), figsize=(4.3*len(data), 4.0*len(rows)), squeeze=False)
for c, (d, lab) in enumerate(zip(data, labels)):
    for r, (name, vmin, vmax) in enumerate(rows):
        ax = axs[r, c]; f = d[name]
        im = ax.pcolormesh(d['x'], d['y'], f, shading='auto', cmap='jet', vmin=vmin, vmax=vmax); plt.colorbar(im, ax=ax)
        ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom); ax.set_aspect('equal'); ax.set_title('%s %s [%.3g,%.3g]' % (lab, name, f.min(), f.max()), fontsize=8)
plt.tight_layout(); plt.savefig(png, dpi=95); print('wrote', png)
