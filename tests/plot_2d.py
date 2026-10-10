"""2D fields from profile files (columns: x y p u v rho alpha2 Y2 ...): numerical schlieren, volume fraction, pressure, ...
usage: python3 tests/plot_2d.py OUT.png TITLE VARS 'label:file' ['label:file' ...]
  VARS: comma list from schlieren, rho, a1, a2, Y2, p, pbar, logp, u, v, speed (one column per variable, one row per file)"""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

COL = {'p': 2, 'u': 3, 'v': 4, 'rho': 5, 'a2': 6, 'Y2': 7}


def field(f):
    d = np.loadtxt(f); o = np.lexsort((d[:, 0], d[:, 1])); d = d[o]
    xs, ys = np.unique(d[:, 0]), np.unique(d[:, 1]); nx, ny = len(xs), len(ys)
    F = {k: d[:, c].reshape(ny, nx) for k, c in COL.items()}
    F['a1'] = 1 - F['a2']; F['pbar'] = F['p'] / 1e5; F['speed'] = np.hypot(F['u'], F['v'])
    with np.errstate(divide='ignore', invalid='ignore'): F['logp'] = np.log10(np.maximum(F['p'], 1.0))
    gy, gx = np.gradient(F['rho'], ys, xs); g = np.hypot(gx, gy)
    F['schlieren'] = np.exp(-15 * g / (g.max() + 1e-300))   # Quirk-type numerical schlieren
    return xs, ys, F


LAB = {'schlieren': 'numerical schlieren', 'rho': 'density (kg/m$^3$)', 'a1': r'$\alpha_1$', 'a2': r'$\alpha_2$', 'Y2': '$Y_2$',
       'p': 'p (Pa)', 'pbar': 'p (bar)', 'logp': r'$\log_{10}p$ (Pa)', 'u': 'u (m/s)', 'v': 'v (m/s)', 'speed': '|u| (m/s)'}
CMAP = {'schlieren': 'gray', 'a1': 'viridis', 'a2': 'viridis', 'Y2': 'viridis'}

if __name__ == '__main__':
    out, title, vars_ = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
    items = [a.split(':', 1) for a in sys.argv[4:]]
    xs, ys, _ = field(items[0][1]); asp = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    w = 4.2; fig, axs = plt.subplots(len(items), len(vars_), figsize=(w * len(vars_), w * asp * len(items) + 0.9), squeeze=False)
    for r, (lab, f) in enumerate(items):
        xs, ys, F = field(f)
        for c, v in enumerate(vars_):
            ax = axs[r, c]
            im = ax.pcolormesh(xs, ys, F[v], shading='auto', cmap=CMAP.get(v, 'turbo'))
            ax.set_aspect('equal'); ax.tick_params(labelsize=7)
            if v != 'schlieren': fig.colorbar(im, ax=ax, fraction=0.046 * max(asp, 0.3), pad=0.02).ax.tick_params(labelsize=7)
            ax.set_title('%s  %s' % (lab, LAB.get(v, v)), fontsize=8)
    fig.suptitle(title, fontsize=9); fig.tight_layout(); fig.savefig(out, dpi=110)
