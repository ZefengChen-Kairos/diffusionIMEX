#!/usr/bin/env python3
"""Spectrum of div u along x inside the blast (|x|<0.1, |y|<0.05) at t=0.01: is the striping a fixed number of cells?"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk
from plot_lines2 import load_ref
def divu_of_dir(p):
    if p.endswith('.dat'):
        d = load_ref(p); x = d['x']; u = d['vx'].T; w = d['vy'].T   # rows = y
    else:
        vt = sorted(glob.glob(os.path.join(p, 'plt_*.vtk')), key=lambda f: int(f.rsplit('_',1)[1].split('.')[0]))[-1]
        v = read_vtk(vt); nx = v['nx']; x = -0.5 + (np.arange(nx)+0.5)/nx; u = v['velX']; w = v['velY']
    h = x[1]-x[0]; d = np.zeros_like(u); d[1:-1,1:-1] = (u[1:-1,2:]-u[1:-1,:-2])/(2*h) + (w[2:,1:-1]-w[:-2,1:-1])/(2*h)
    return x, d
args = sys.argv[1:]; paths = args[0::2]; labels = args[1::2]
print('%-26s %8s %14s %12s' % ('run', 'N', 'peak period', 'rms(div u)'))
for p, l in zip(paths, labels):
    x, d = divu_of_dir(p); N = len(x)
    ix = np.abs(x) < 0.1; iy = np.abs(x) < 0.05
    block = d[np.ix_(iy, ix)]
    block = block - block.mean(axis=1, keepdims=True)
    spec = np.mean(np.abs(np.fft.rfft(block, axis=1))**2, axis=0)
    k = np.fft.rfftfreq(block.shape[1], d=1.0)   # cycles per cell
    spec[0] = 0; kk = k[1:]; ss = spec[1:]
    # keep only the sub-physical range (period <= 16 cells) for the peak
    sel = kk >= 1.0/16.0
    kp = kk[sel][np.argmax(ss[sel])]
    print('%-26s %8d %10.2f cells  %10.1f   (power in period<=8 cells: %.0f%%)' % (l, N, 1.0/kp, np.sqrt(np.mean(block**2)), 100*ss[kk >= 1.0/8.0].sum()/ss.sum()))
