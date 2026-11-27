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

def bake_diagnostics(s):
    """STEP 2: remove the diagnostic switches, writing in the chosen value.

    The decision and the evidence for each is in blastwave_debug/PATCHES.md.
      g_chisplit  FALSIFIED: NaN on Orszag-Tang with the density point-value
            guard fired 1.25e8 times, against zero without it.  The division
            guard set the indicator to one where the eps-weighted wave speed
            vanishes, which at a magnetic null is a large region.  Not here.
      g_jguard    the switch goes; the two J-keyed safeguards stay, now
            unconditional as they were before.
      g_rho_scale written in as ON.  The scaling form perturbs only by the
            excess beyond the floor, and the guard it governs fires zero times
            on both 2D tests, so it is risk-free.
      g_rho_thr   written in as the named constant RHO_PT_FLOOR.
      g_idp, g_hllemoff  probes only.
    """
    n0 = len(s)

    # ---- fix I: gone entirely ----
    s = cut(s, 'double switchS(double M);   // [fix I] defined below; chi needs the split weights\n', 'fix I decl')
    s = s.replace('double switchS(double M) {', 'inline double switchS(double M) {')
    m = re.search(r'[\t ]*double cf;\n[\t ]*if\(g_chisplit\)\{.*?\}else\{\n([\t ]*cf = std::sqrt\(a2 \+ cA2\);\n)[\t ]*\}\n', s, re.S)
    assert m, 'fix I cf block not found'
    s = s[:m.start()] + m.group(1) + s[m.end():]
    m = re.search(r'[\t ]*double denominator = shock_kappa1 \* cf_min_nbr;\n[\t ]*if\(!\(denominator > 0\.0\)\)\{[^\n]*\n[\t ]*shock_eta_tilde\(i, j\) = \(div_u < 0\.0\) \? 1\.0 : 0\.0;\n[\t ]*\}else\{\n(.*?)\n[\t ]*\}\n', s, re.S)
    assert m, 'fix I chi block not found'
    body = '\n'.join(l[1:] if l.startswith('\t\t\t\t') else l for l in m.group(1).split('\n'))
    s = (s[:m.start()]
         + '\t\t\tdouble numerator   = h * div_u + shock_kappa1 * cf_min_nbr;\n'
         + '\t\t\tdouble denominator = shock_kappa1 * cf_min_nbr;\n\t\t\t\n'
         + body + '\n' + s[m.end():])

    # ---- fix K: the switch goes, the safeguards stay ----
    s = s.replace('if(g_jguard){ const double J = ', '{ const double J = ')
    s = s.replace('if(g_jguard && Jw > 0.5 && g_env_cap > 0.0){   // [fix K]',
                  'if(Jw > 0.5 && g_env_cap > 0.0){')

    # ---- fix J: written in as ON, with a named floor ----
    s = s.replace('const double thr = g_rho_thr*in(i,j,0);   // [fix J]',
                  'const double thr = RHO_PT_FLOOR*in(i,j,0);')
    s = s.replace('out(i,j,0) = g_rho_scale ? thr : in(i,j,0); ng++; } }',
                  'out(i,j,0) = thr; ng++; } }   // scale the increment to the floor')
    s = s.replace('if(!(Qpt(i,j,0) >= g_rho_thr*Qbar(i,j,0))){ Qpt(i,j,0) = g_rho_scale ? g_rho_thr*Qbar(i,j,0) : Qbar(i,j,0); n++; }   // [fix J]',
                  'if(!(Qpt(i,j,0) >= RHO_PT_FLOOR*Qbar(i,j,0))){ Qpt(i,j,0) = RHO_PT_FLOOR*Qbar(i,j,0); n++; }')

    # ---- the globals of fixes I/J/K ----
    s = cut(s, 'int g_jguard = 1;      // [fix K] 1: keep the two J-keyed safeguards of HLLEM_x; 0: both off\n'
               'int g_chisplit = 0;    // [fix I] 1: gate chi on sqrt(eps_M a^2 + eps_B c_A^2) instead of sqrt(a^2 + c_A^2)\n'
               'int g_rho_scale = 0;   // [fix J] 1: scale the deconvolution increment to the threshold instead of discarding it\n'
               'double g_rho_thr = 0.5;// [fix J] threshold of the density point-value guard (original: 0.5)\n',
            'I/J/K globals')
    s = s.replace('double shock_kappa1 = 0.4;',
                  '// Floor of the density point value relative to its cell average, enforced by\n'
                  '// scaling the deconvolution increment down to it.  It bounds the derived point\n'
                  '// velocity m_pt/rho_pt by 1/RHO_PT_FLOOR times the cell-average one.\n'
                  'const double RHO_PT_FLOOR = 0.5;\n'
                  'double shock_kappa1 = 0.4;', 1)

    # ---- runtime keys of the removed switches ----
    s = s.replace(' else if(key=="idp") g_idp = (int)v;\n'
                  '\t\telse if(key=="chisplit") g_chisplit = (int)v; else if(key=="rhoscale") g_rho_scale = (int)v; else if(key=="rhothr") g_rho_thr = v;\n'
                  '\t\telse if(key=="jguard") g_jguard = (int)v;', '')
    s = s.replace(' else if(key=="hllemoff") g_hllemoff = (int)v;', '')

    print('  diagnostics baked: %d chars removed' % (n0 - len(s)))
    return s


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    step = sys.argv[3] if len(sys.argv) > 3 else 'all'
    s = open(src).read()
    if step in ('all', '1'):
        s = dead_code(s)
    if step in ('all', '2'):
        s = bake_diagnostics(s)
    open(dst, 'w').write(s)
    print('wrote', dst)
