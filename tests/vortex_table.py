#!/usr/bin/env python3
"""Battisti-Boscheri 4.1 two-phase isentropic vortex (case 56): L2 errors vs exact solution and orders.
usage: python3 tests/vortex_table.py RUNDIR TAG   (reads RUNDIR/TAG_p{p0}_N{N}.log written by the IMEX program)"""
import sys, re, math, os, glob
d, tag = sys.argv[1], sys.argv[2]
cols = ['a2r2', 'a2r2u', 'a2', 'p', 'u']
print('| p0 | N | ' + ' | '.join(f'L2({c}) | ord' for c in cols) + ' | steps |')
print('|' + '---|' * (3 + 2 * len(cols)))
for p0 in ['1', '10', '100', '1000', '10000']:
    prev = None
    for N in sorted(int(re.search(r'_N(\d+)\.log$', g).group(1)) for g in glob.glob(os.path.join(d, f'{tag}_p{p0}_N*.log'))):
        f = os.path.join(d, f'{tag}_p{p0}_N{N}.log')
        if not os.path.exists(f): continue
        s = open(f).read(); m = re.search(r'^L2 .*$', s, re.M)
        if not m: print(f'| {p0} | {N} | (failed / running) |'); continue
        e = {c: float(re.search(rf' {c} ([0-9.eE+-]+)', m.group()).group(1)) for c in cols}
        st = re.search(r'steps (\d+)', s).group(1)
        row = []
        for c in cols:
            o = '–' if prev is None else f'{math.log2(prev[c] / e[c]):.2f}'
            row += [f'{e[c]:.4e}', o]
        print(f'| {p0} | {N} | ' + ' | '.join(row) + f' | {st} |'); prev = e
