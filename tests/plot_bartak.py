"""Bartak blowdown (Pelanti 2022 8.2.4, case 32): probe histories at x = 48 mm against the exact HEM simple wave (tests/exact_hem_bartak.py).
Copies the probe files and the exact curve to DATADIR.
usage: python3 tests/plot_bartak.py RUNDIR/bartak DATADIR OUT.png"""
import os, sys, shutil
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
RUN, DATA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
argv, sys.argv = sys.argv, ['x']; ns = {}
exec(open(os.path.join(H, 'exact_hem_bartak.py')).read().split("print('psat(T0)")[0], ns); sys.argv = argv
tab, ps, x = ns['tab'], ns['ps'], 0.048
ts = np.linspace(1e-5, 17.5e-3, 400); pe, Ye = [], []
for t in ts:
    xi = x / t
    if xi >= tab[0, 3]: pe.append(ps); Ye.append(0.0)
    else: j = np.argmax(tab[:, 3] <= xi); pe.append(tab[j, 0]); Ye.append(tab[j, 4])
os.makedirs(DATA, exist_ok=True)
np.savetxt(os.path.join(DATA, 'exact_hem_probe.txt'), np.c_[ts, pe, Ye], header='t p Y  exact HEM simple wave at x = 48 mm (tests/exact_hem_bartak.py)')
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for n, lab, c in [('nu', r'finite $\nu$ (Pelanti eq. 84)', 'tab:red'), ('inf', r'instantaneous ($pTg$)', 'tab:blue'), ('noMT', r'no mass transfer ($pT$)', 'tab:green')]:
    f = os.path.join(RUN, '%s.dat.probe' % n); shutil.copy(f, DATA)
    d = np.loadtxt(f); ax[0].plot(d[:, 0] * 1e3, d[:, 1] / 1e6, c, label=lab); ax[1].plot(d[:, 0] * 1e3, d[:, 2], c, label=lab)
ax[0].plot(ts * 1e3, np.array(pe) / 1e6, 'k--', lw=1, label='exact HEM (instantaneous)'); ax[1].plot(ts * 1e3, Ye, 'k--', lw=1, label='exact HEM')
ax[0].plot([0.3, 17.5], [4.3, 4.3], 'k:', lw=1, label='Pelanti Fig. 8 finite rate ~4.3, experiment ~4.5 (read)')
ax[0].set_xlabel('t (ms)'); ax[0].set_ylabel('p at x = 48 mm (MPa)'); ax[1].set_xlabel('t (ms)'); ax[1].set_ylabel('vapour mass fraction Y at x = 48 mm')
for a in ax: a.grid(alpha=.3)
ax[0].legend(fontsize=7); ax[1].legend(fontsize=7)
fig.suptitle('Barták blowdown (Pelanti 2022 8.2.4): NASG water (Table 5), 1.7 m tube, p = 1e5 Pa at the open end, 1000 cells, IMEX', fontsize=9); fig.tight_layout()
fig.savefig(OUT, dpi=110)
