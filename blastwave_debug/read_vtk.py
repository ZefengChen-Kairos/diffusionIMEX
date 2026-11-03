#!/usr/bin/env python3
"""Minimal reader for the ASCII RECTILINEAR_GRID VTK files written by the codes (FIELD data, cell data)."""
import numpy as np
def read_vtk(fn):
    with open(fn) as f: txt = f.read().split('\n')
    i = 0; out = {}; nx = ny = None
    while i < len(txt):
        line = txt[i].strip()
        if line.startswith('DIMENSIONS'):
            p = line.split(); nx = int(p[1]) - 1; ny = int(p[2]) - 1
        if line.startswith('TIME 1 1 double'):
            out['t'] = float(txt[i+1])
        if nx is not None and len(line.split()) == 4 and line.split()[1] == '1' and line.split()[3] == 'double' and not line.startswith('TIME'):
            name = line.split()[0]; n = int(line.split()[2])
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(float(v) for v in txt[j].split()); j += 1
            out[name] = np.array(vals).reshape(ny, nx)   # rows = j (y), cols = i (x)
            i = j; continue
        i += 1
    out['nx'] = nx; out['ny'] = ny
    return out
