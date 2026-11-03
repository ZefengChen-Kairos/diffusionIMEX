#!/usr/bin/env python3
"""Zoomed interior comparison of several runs from their final VTK files: p, div u, |u|^2, eps_M."""
import sys, glob, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from read_vtk import read_vtk
args = sys.argv[1:]; png = args[-1]; runs = args[:-1]
dirs = runs[0::2]; labels = runs[1::2]
rows = [('p', None), ('div u', None), ('|u|^2', None), ('eps_M', None), ('rho', None)]
fig, axs = plt.subplots(len(rows), len(dirs), figsize=(4.3*len(dirs), 4.0*len(rows)), squeeze=False)
zoom = 0.3
for c, (d, lab) in enumerate(zip(dirs, labels)):
    vtks = sorted(glob.glob(os.path.join(d, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))
    v = read_vtk(vtks[-1]); nx, ny = v['nx'], v['ny']
    x = -0.5 + (np.arange(nx)+0.5)/nx; y = -0.5 + (np.arange(ny)+0.5)/ny; dx = 1.0/nx
    u, w = v['velX'], v['velY']
    divu = np.zeros_like(u); divu[1:-1,1:-1] = (u[1:-1,2:]-u[1:-1,:-2])/(2*dx) + (w[2:,1:-1]-w[:-2,1:-1])/(2*dx)
    F = {'p': v['p'], 'div u': divu, '|u|^2': u*u+w*w, 'eps_M': v['eps_M'], 'rho': v['rho']}
    for r, (name, _) in enumerate(rows):
        ax = axs[r, c]; f = F[name]
        if name == 'div u': vmin, vmax = -400, 400
        elif name == 'p': vmin, vmax = 0, 220
        elif name == '|u|^2': vmin, vmax = 0, 250
        elif name == 'eps_M': vmin, vmax = 0, 1
        else: vmin, vmax = 0.1, 3.6
        im = ax.pcolormesh(x, y, f, shading='auto', cmap='jet', vmin=vmin, vmax=vmax); plt.colorbar(im, ax=ax)
        ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom); ax.set_aspect('equal')
        ax.set_title('%s %s t=%.3g [%.3g,%.3g]' % (lab, name, v['t'], f.min(), f.max()), fontsize=8)
plt.tight_layout(); plt.savefig(png, dpi=95); print('wrote', png)
