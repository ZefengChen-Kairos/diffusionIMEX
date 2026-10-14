"""Evaporation-front velocity (Pelanti 2022 8.2.3 / Simoes-Moreira & Shepherd 1999), case 33.
Front x_f(t): where the vapour volume fraction alpha2 first crosses 0.5 (linear interpolation), at t = 400, 425, 450, 475, 500 us.
Front speed into the liquid in the laboratory frame U_lab = -dx_f/dt (least squares over the snapshots; the experiments measure this),
and relative to the incoming liquid U_rel = u_liq - dx_f/dt, u_liq = mean u over [x_f - 2 mm, x_f - 1 mm].
Also the Saurel et al. (2008) jump formula (laboratory frame, sign flipped so that positive = into the liquid) U_S = ((rho u)_a - (rho u)_b) / (rho_a - rho_b), b = liquid ahead, a = first cell with alpha2 > A (A = 0.9, 0.99).
usage: python3 tests/evap_front.py RUNDIR [T ...]   (files t<T>.dat, t<T>_s0..3.dat)"""
import numpy as np, sys, os
run = sys.argv[1]; Ts = sys.argv[2:] or sorted({f[1:].split('.')[0].split('_')[0] for f in os.listdir(run) if f.startswith('t') and f.endswith('.dat')})
times = [400e-6, 425e-6, 450e-6, 475e-6, 500e-6]
print('| T (K) | U_lab (tracking, m/s) | U_rel | U_S (a: alpha2>0.9) | U_S (a: alpha2>0.99) | u_liq (m/s) | rho_liq |')
print('|---|---|---|---|---|---|---|')
for T in Ts:
    fs = [os.path.join(run, 't%s_s%d.dat' % (T, k)) for k in range(4)] + [os.path.join(run, 't%s.dat' % T)]
    if not all(os.path.exists(f) for f in fs): continue
    xf, ul, US = [], [], {0.9: [], 0.99: []}
    for f in fs:
        d = np.loadtxt(f); x, u, r, a = d[:, 0], d[:, 3], d[:, 5], d[:, 6]
        i = np.argmax(a > 0.5); xf.append(x[i - 1] + (0.5 - a[i - 1]) / (a[i] - a[i - 1]) * (x[i] - x[i - 1]))
        m = (x > xf[-1] - 2e-3) & (x < xf[-1] - 1e-3); ub, rb = u[m].mean(), r[m].mean(); ul.append(ub)
        for A in US: j = np.argmax(a > A); US[A].append((r[j] * u[j] - rb * ub) / (r[j] - rb))
    v = np.polyfit(times, xf, 1)[0]
    print('| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.1f |' % (T, -v, np.mean(ul) - v, -np.mean(US[0.9]), -np.mean(US[0.99]), np.mean(ul), rb))
