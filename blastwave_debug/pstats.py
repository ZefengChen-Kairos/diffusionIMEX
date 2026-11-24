#!/usr/bin/env python3
"""Pressure statistics of a blast wave output, for checking a figure against the data.

    python3 pstats.py path/to/plt_blast_wave_xxx.vtk [ambient_p] [blast_radius]

Prints what the log10 p colour scale should be, and the far-field distribution
(percentiles, not p_min: the minimum is set by a single cell and is the worst
possible summary of this defect -- see REPORT.md 6.4).
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_vtk import read_vtk

f = sys.argv[1]
p_amb = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
r_far = float(sys.argv[3]) if len(sys.argv) > 3 else 0.33

v = read_vtk(f)
p, nx, ny = v['p'], v['nx'], v['ny']
print('%s\n  t = %.6g   mesh = %d x %d' % (os.path.basename(f), v.get('t', float('nan')), nx, ny))

lo, hi = p.min(), p.max()
print('  p in [%.6g, %.6g]' % (lo, hi))
print('  log10 p in [%s, %.4f]   <- what an auto-scaled colour bar shows'
      % ('undefined (p < 0!)' if lo <= 0 else '%.4f' % np.log10(lo), np.log10(hi)))
if lo <= 0:
    print('  WARNING: %d cells (%.2f%%) have p <= 0; np.log10 gives NaN and matplotlib'
          % ((p <= 0).sum(), 100*(p <= 0).mean()))
    print('           silently masks them, so they leave NO trace on the figure.')

# consistency: is the stored p really (g-1)(U - KE - |B|^2/2) ?
if all(k in v for k in ('U', 'rho', 'mom_X', 'B_X')):
    ke = 0.5*(v['mom_X']**2 + v['mom_Y']**2 + v['mom_Z']**2)/np.maximum(v['rho'], 1e-300)
    mb = 0.5*(v['B_X']**2 + v['B_Y']**2 + v['B_Z']**2)
    for g in (1.4, 5.0/3.0):
        d = np.abs(p - (g - 1.0)*(v['U'] - ke - mb)).max()
        print('  gamma = %.4f -> max |p_stored - (g-1)(U - KE - MB)| = %.3g%s'
              % (g, d, '   <- this is the run gamma' if d < 1e-8 else ''))

x = np.linspace(-0.5 + 0.5/nx, 0.5 - 0.5/nx, nx)
y = np.linspace(-0.5 + 0.5/ny, 0.5 - 0.5/ny, ny)
r = np.hypot(*np.meshgrid(x, y))
q = p[r > r_far]
print('  far field (r > %.2f), exact answer p = %g everywhere, %d cells:' % (r_far, p_amb, q.size))
print('     percentiles  1%%: %-11.5g 10%%: %-11.5g 50%%: %-11.5g 90%%: %.5g' %
      tuple(np.percentile(q, k) for k in (1, 10, 50, 90)))
print('     frac(p < %g) = %.4f   frac(p < %g) = %.4f   frac(p < 0) = %.4f' %
      (p_amb, (q < p_amb).mean(), 0.5*p_amb, (q < 0.5*p_amb).mean(), (q < 0).mean()))
print('  reference (independent explicit HLL, 1024^2): p in [0.0982, 225.5], log10 [-1.01, 2.353],')
print('            far-field median 0.1, frac(p < 0.1) = 0.005, frac(p < 0.05) = 0')
