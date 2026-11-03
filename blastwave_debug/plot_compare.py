#!/usr/bin/env python3
"""Compare blast-wave results of the two codes (RK_*.dat: x y value, written at the final time)."""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_dir(d):
    out = {}
    for name in ['rho','p','vx','vy','Bx','By','Bz','E']:
        f = os.path.join(d, 'RK_%s.dat' % name)
        if not os.path.exists(f): continue
        a = np.loadtxt(f)
        xs = np.unique(a[:,0]); ys = np.unique(a[:,1])
        nx, ny = len(xs), len(ys)
        # file order: i outer, j inner
        v = a[:,2].reshape(nx, ny)
        out[name] = v
        out['x'] = xs; out['y'] = ys
    if 'vx' in out:
        out['u2'] = out['vx']**2 + out['vy']**2
    if 'Bx' in out:
        out['m'] = 0.5*(out['Bx']**2 + out['By']**2 + out.get('Bz', 0*out['Bx'])**2)
    return out

def main(dirs, labels, png):
    data = [load_dir(d) for d in dirs]
    fields = [('rho', r'$\rho$'), ('p', r'$p$'), ('u2', r'$|u|^2$'), ('m', r'$m=|B|^2/2$')]
    n = len(data)
    fig, axs = plt.subplots(len(fields), n + 1, figsize=(4.2*(n+1), 3.6*len(fields)))
    for r, (f, lab) in enumerate(fields):
        vmin = min(np.nanmin(d[f]) for d in data); vmax = max(np.nanmax(d[f]) for d in data)
        for c, d in enumerate(data):
            ax = axs[r, c]
            im = ax.pcolormesh(d['x'], d['y'], d[f].T, shading='auto', vmin=vmin, vmax=vmax, cmap='jet')
            ax.set_aspect('equal'); ax.set_title('%s %s  [%.3g, %.3g]' % (labels[c], lab, np.nanmin(d[f]), np.nanmax(d[f])), fontsize=9)
            plt.colorbar(im, ax=ax)
        ax = axs[r, n]
        for c, d in enumerate(data):
            jmid = len(d['y'])//2; imid = len(d['x'])//2
            ax.plot(d['x'], d[f][:, jmid], '-', lw=1, label='%s y=0' % labels[c])
            ax.plot(d['y'], d[f][imid, :], '--', lw=1, label='%s x=0' % labels[c])
        ax.set_title('line-outs ' + lab, fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=.3)
        if f in ('p',): ax.set_yscale('log')
    plt.tight_layout(); plt.savefig(png, dpi=110)
    print('wrote', png)
    for c, d in enumerate(data):
        print(labels[c], 'rho[min,max]=%.4g %.4g  p[min,max]=%.4g %.4g  |u|max=%.4g  m[min,max]=%.4g %.4g' % (
            d['rho'].min(), d['rho'].max(), d['p'].min(), d['p'].max(), np.sqrt(d['u2'].max()), d['m'].min(), d['m'].max()))

if __name__ == '__main__':
    args = sys.argv[1:]
    png = args[-1]; rest = args[:-1]
    dirs = rest[0::2]; labels = rest[1::2]
    main(dirs, labels, png)
