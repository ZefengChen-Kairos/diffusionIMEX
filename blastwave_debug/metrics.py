#!/usr/bin/env python3
"""Roughness metrics of the blast-wave interior (|x|,|y| < 0.13 at t = 0.01), one line per run.
   wiggle(v_y | x=0) = TV(v_y along x=0, |y|<0.13) / (max-min)   (1 for a monotone profile; >1 = oscillations)
   rough(v_y)        = ||Lap_h v_y||_2 / ||grad_h v_y||_2 * h     (grid-scale roughness, 2D interior)"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk
from plot_lines2 import load_ref
from plot_compare import load_dir
def load(p):
    if p.endswith('.dat'):
        d = load_ref(p); return d['x'], d['y'], d['vy'], d['vx'], d['p']
    vt = sorted(glob.glob(os.path.join(p, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))
    if vt:
        v = read_vtk(vt[-1]); nx, ny = v['nx'], v['ny']; x = -0.5 + (np.arange(nx)+0.5)/nx; y = -0.5 + (np.arange(ny)+0.5)/ny
        return x, y, v['velY'].T, v['velX'].T, v['p'].T
    d = load_dir(p); return d['x'], d['y'], d['vy'], d['vx'], d['p']
def metrics(p):
    x, y, vy, vx, pr = load(p); nx = len(x); h = x[1]-x[0]
    ix = np.abs(x) < 0.13; iy = np.abs(y) < 0.13
    imid = nx//2; line = 0.5*(vy[imid, :] + vy[imid-1, :])[iy]
    tv = np.sum(np.abs(np.diff(line))); rng = line.max() - line.min()
    f = vy[np.ix_(ix, iy)]
    lap = f[2:,1:-1] + f[:-2,1:-1] + f[1:-1,2:] + f[1:-1,:-2] - 4*f[1:-1,1:-1]
    gx = 0.5*(f[2:,1:-1] - f[:-2,1:-1]); gy = 0.5*(f[1:-1,2:] - f[1:-1,:-2])
    rough = np.sqrt(np.sum(lap**2)) / np.sqrt(np.sum(gx**2 + gy**2) + 1e-300)
    u2 = (vx**2 + vy**2); l2 = 0.5*(u2[imid, :] + u2[imid-1, :])[iy]
    return tv/rng, rough, l2.max()
if __name__ == '__main__':
    args = sys.argv[1:]; paths = args[0::2]; labels = args[1::2]
    print('%-32s %10s %10s %12s' % ('run', 'wiggle_vy', 'rough_vy', 'max|u|2(x=0)'))
    for p, l in zip(paths, labels):
        w, r, m = metrics(p); print('%-32s %10.3f %10.3f %12.1f' % (l, w, r, m))
