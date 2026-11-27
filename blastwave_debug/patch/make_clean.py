#!/usr/bin/env python3
"""Build the clean single-source version from Final_emb_omp_original.cpp.

    python3 make_clean.py Final_emb_omp_original.cpp Final_emb_omp_clean.cpp

STEP 1 (this file, so far): remove dead code only.  Deleting code that is never
reached must not change a single bit of output, which is how each step is checked.
"""
import sys, re

def cut(s, old, tag, count=1):
    assert s.count(old) >= count, f"{tag}: pattern not found\n{old[:160]}"
    return s.replace(old, '', count)

def dead_code(s):
    n0 = len(s)

    # --- 1. six RusanovK_* forward declarations: no definition, no call anywhere ---
    for f in ('void RusanovK_x(double *QL, double *QR, double *flux);\n',
              'void RusanovK_y(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);\n',
              'void RusanovK_Pressure_x(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Pressure_y(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Update_Pressure(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);\n'):
        s = cut(s, f, 'RusanovK decl')

    # --- 2. epsilon_min / epsilon_max: declared, referenced nowhere else ---
    # item 3 of the porting notes describes exactly these two variables; items
    # 1, 2 and 4 document real code and stay.
    item3 = '''//   3. epsilon_min/epsilon_max didn't exist in this file -- added as
//      globals, defaulted to [0,1] (no artificial floor). Note: the WENO3
//      characteristic reconstruction's condition number was measured to
//      blow up for eps_min below ~0.05-0.1 (see earlier conversation) --
//      if that turns out to matter in practice, raising epsilon_min is a
//      one-line fix here, not a structural change.
'''
    s = cut(s, item3, 'porting note 3')
    m = re.search(r'// EPSILON BOUNDS \(per author instruction\).*?double epsilon_max = 1\.0;\n', s, re.S)
    assert m, 'epsilon bounds block not found'
    s = s[:m.start()] + s[m.end():]
    assert s.count('epsilon_min') == 0 and s.count('epsilon_max') == 0, \
        f"epsilon_min/max still referenced after removing the block: {s.count('epsilon_min')}/{s.count('epsilon_max')}"   # the block held the declarations AND every mention

    # --- 3. Bx_amb: stale B0 = 70 constant, never read ---
    assert s.count('Bx_amb') == 1, f'Bx_amb referenced {s.count("Bx_amb")} times'
    s = s.replace('double rho_amb = 1.0, p_amb = 0.1, Bx_amb = 70.0/sqrt(4.0*M_PI);',
                  'double rho_amb = 1.0, p_amb = 0.1;')

    # --- 4. [removed: ...] comment stubs left by earlier edits ---
    s = re.sub(r'^// \[removed: [^\]]*\]\n', '', s, flags=re.M)

    print(f'  dead code: {n0 - len(s)} chars removed')
    return s

if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    s = open(src).read()
    s = dead_code(s)
    open(dst, 'w').write(s)
    print('wrote', dst)
