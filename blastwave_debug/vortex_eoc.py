#!/usr/bin/env python3
"""EOC table of the travelling (and stationary) vortex from the VORTEX lines of the run logs."""
import re, math, sys, os
S = sys.argv[1]
def table(prefix, cfgs, Ns, title):
    res = {}
    for cfg in cfgs:
        for N in Ns:
            fn = f'{S}/{prefix}_{cfg}_{N}.log'
            if not os.path.exists(fn): continue
            m = re.search(r'VORTEX N=\d+ rho0=\S+ t=\S+ L2: rho=(\S+) u=(\S+) p=(\S+) Bx=(\S+) Az=(\S+)', open(fn).read())
            if m: res[(cfg, N)] = [float(v) for v in m.groups()]
    names = ['rho', 'u', 'p', 'Bx', 'Az']
    print(title); print('%-5s %4s ' % ('cfg', 'N') + ' '.join('%11s %5s' % (n, 'EOC') for n in names))
    for cfg in cfgs:
        prev = None
        for N in Ns:
            if (cfg, N) not in res: continue
            e = res[(cfg, N)]; row = '%-5s %4d ' % (cfg, N)
            for k in range(5):
                row += ' %11.4e %5.2f' % (e[k], (math.log(prev[k]/e[k])/math.log(2) if prev else float('nan')))
            print(row); prev = e
    print()
table('v', ['orig', 'A', 'AC'], [32, 64, 128, 256], 'translating vortex (paper 2 sec 4.1: rho0=1, Bt=1, vt=sqrt(2pi), u0=1, L=5, t=1, CFL 0.8), L2 errors vs exact cell averages')
table('sv', ['orig', 'AC'], [32, 64, 128], 'stationary vortex (u0=0, L=10, t=1, CFL 0.8)')
