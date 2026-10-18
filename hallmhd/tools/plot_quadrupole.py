"""Hall quadrupole plot (B_z with in-plane field lines) from the final RK_*.dat of a run.

usage: python3 plot_quadrupole.py <run out= dir> <title> <out.png>
"""
import sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(d, name):
    f = glob.glob(d + '/**/RK_' + name + '.dat', recursive=True)[0]
    a = np.loadtxt(f)
    nx, ny = len(np.unique(a[:, 0].round(9))), len(np.unique(a[:, 1].round(9)))
    # rebuild exactly equidistant coordinates (the text output is not, at round-off level)
    x = np.linspace(a[:, 0].min(), a[:, 0].max(), nx)
    y = np.linspace(a[:, 1].min(), a[:, 1].max(), ny)
    return x, y, a[:, 2].reshape(nx, ny).T                   # [j, i]


def main():
    d, title, out = sys.argv[1], sys.argv[2], sys.argv[3]
    x, y, Bz = load(d, 'Bz')
    _, _, Bx = load(d, 'Bx')
    _, _, By = load(d, 'By')
    dx, dy = x[1] - x[0], y[1] - y[0]
    # flux function A_z (Bx = dAz/dy, By = -dAz/dx): integrate up the left column, then along x
    Az = np.zeros_like(Bx)
    Az[1:, 0] = np.cumsum(0.5 * (Bx[1:, 0] + Bx[:-1, 0]) * dy)
    Az[:, 1:] = Az[:, [0]] - np.cumsum(0.5 * (By[:, 1:] + By[:, :-1]) * dx, axis=1)
    # X point: saddle of A_z on the sheet row = maximum of A_z along x there
    j0 = np.argmin(np.abs(y))
    iX = np.argmax(Az[j0, :])

    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    lim = np.abs(Bz).max()
    pm = ax.pcolormesh(x, y, Bz, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='nearest')
    ax.streamplot(x, y, Bx, By, color='k', linewidth=0.8, density=1.6, arrowsize=0.6)
    ax.plot(x[iX], y[j0], marker='*', ms=14, mfc='white', mec='k')
    ax.set_xlim(x[0] - dx / 2, x[-1] + dx / 2)
    ax.set_ylim(y[0] - dy / 2, y[-1] + dy / 2)
    ax.set_xlabel(r'$x\ [d_i]$')
    ax.set_ylabel(r'$y\ [d_i]$')
    ax.set_title(title)
    fig.colorbar(pm, ax=ax, label=r'$B_z$')
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print('X point at x=%.3f y=%.3f, |Bz|max=%.4f' % (x[iX], y[j0], lim))


if __name__ == '__main__':
    main()
