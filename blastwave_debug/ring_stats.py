#!/usr/bin/env python3
"""Where are the density/pressure streaks?  Quadrant zoom (x<0, y>0.1) of rho, p, eps_M, eps_B, |u|, M_A and ring statistics."""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk
from plot_lines2 import load_ref
args = sys.argv[1:]; png = args[-1]; paths = args[:-1][0::2]; labels = args[:-1][1::2]
rows = [('rho', 0.9, 1.5), ('p', 0, 12), ('eps_M', 0, 1), ('eps_B', 0, 1), ('|u|', 0, 16), ('M_A', 0, 1.5)]
fig, axs = plt.subplots(len(rows), len(paths), figsize=(4.2*len(paths), 3.9*len(rows)), squeeze=False)
for c, (p, l) in enumerate(zip(paths, labels)):
    if p.endswith('.dat'):
        d = load_ref(p); x = d['x']; F = {'rho': d['rho'].T, 'p': d['p'].T, 'u': d['vx'].T, 'v': d['vy'].T, 'Bx': d['Bx'].T, 'By': d['By'].T, 'eps_M': None, 'eps_B': None}
    else:
        vt = sorted(glob.glob(os.path.join(p, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))[-1]
        v = read_vtk(vt); nx = v['nx']; x = -0.5 + (np.arange(nx)+0.5)/nx
        F = {'rho': v['rho'], 'p': v['p'], 'u': v['velX'], 'v': v['velY'], 'Bx': v['B_X'], 'By': v['B_Y'], 'eps_M': v['eps_M'], 'eps_B': v['eps_B']}
    X, Y = np.meshgrid(x, x); r = np.sqrt(X*X+Y*Y); th = np.degrees(np.arctan2(Y, np.abs(X)))
    U = np.sqrt(F['u']**2 + F['v']**2); cB = np.sqrt((F['Bx']**2+F['By']**2)/F['rho']); a = np.sqrt(1.4*F['p']/F['rho'])
    F['|u|'] = U; F['M_A'] = U/np.maximum(cB, 1e-12); M = U/a
    ring = (r > 0.30) & (r < 0.42) & (th > 35) & (th < 85)
    txt = '%-22s ring(0.30<r<0.42, 35<th<85): M med=%.2f  M_A med=%.2f  beta med=%.3g  rho[%.2f,%.2f]  p[%.2f,%.2f]' % (l, np.median(M[ring]), np.median(F['M_A'][ring]), np.median(F['p'][ring]/(0.5*(F['Bx']**2+F['By']**2)))[ring] if False else np.median((F['p']/(0.5*(F['Bx']**2+F['By']**2)))[ring]), F['rho'][ring].min(), F['rho'][ring].max(), F['p'][ring].min(), F['p'][ring].max())
    if F['eps_M'] is not None: txt += '  eps_M med=%.2f frac>0.5=%.2f  eps_B med=%.2f frac>0.5=%.2f' % (np.median(F['eps_M'][ring]), np.mean(F['eps_M'][ring] > 0.5), np.median(F['eps_B'][ring]), np.mean(F['eps_B'][ring] > 0.5))
    print(txt)
    for rr, (name, vmin, vmax) in enumerate(rows):
        ax = axs[rr, c]; f = F[name]
        if f is None: ax.axis('off'); continue
        im = ax.pcolormesh(x, x, f, shading='auto', cmap='jet', vmin=vmin, vmax=vmax); plt.colorbar(im, ax=ax)
        ax.set_xlim(-0.45, 0.0); ax.set_ylim(0.1, 0.45); ax.set_aspect('equal'); ax.set_title('%s %s' % (l, name), fontsize=9)
plt.tight_layout(); plt.savefig(png, dpi=90); print('wrote', png)
