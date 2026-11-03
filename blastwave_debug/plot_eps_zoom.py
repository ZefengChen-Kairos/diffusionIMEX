#!/usr/bin/env python3
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from read_vtk import read_vtk
fn = sys.argv[1]; png = sys.argv[2]
d = read_vtk(fn); nx, ny = d['nx'], d['ny']
x = -0.5 + (np.arange(nx) + 0.5)/nx; y = -0.5 + (np.arange(ny) + 0.5)/ny
rho = d['rho']; p = d['p']; u = d['velX']; v = d['velY']; Bx = d['B_X']; By = d['B_Y']
divu = np.zeros_like(u); divu[1:-1,1:-1] = (u[1:-1,2:]-u[1:-1,:-2])/(2*x[1]-2*x[0]) + (v[2:,1:-1]-v[:-2,1:-1])/(2*y[1]-2*y[0])
fields = [('eps_M', d['eps_M'], None), ('eps_B', d['eps_B'], None), ('p', p, None), ('velY', v, None), ('velX', u, None), ('div u', divu, None),
          ('rho', rho, None), ('B_Y', By, None), ('|u|^2 (zoom |x|,|y|<0.25)', u**2+v**2, 0.25)]
fig, axs = plt.subplots(3, 3, figsize=(15, 13))
for ax, (name, f, zoom) in zip(axs.flat, fields):
    im = ax.pcolormesh(x, y, f, shading='auto', cmap='jet'); plt.colorbar(im, ax=ax); ax.set_aspect('equal')
    ax.set_title('%s t=%.4g [%.3g, %.3g]' % (name, d['t'], f.min(), f.max()), fontsize=9)
    if zoom: ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom)
plt.tight_layout(); plt.savefig(png, dpi=100); print('wrote', png)
