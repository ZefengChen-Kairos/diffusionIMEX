#!/usr/bin/env python3
"""Build the clean single-source version.

    python3 make_clean.py Final_emb_omp_patched.cpp Final_emb_omp_clean.cpp [step]

Input is the PATCHED source, so every mechanism is present and working; this
script removes what is dead and writes in the chosen value of what is kept.
Each step is checked the same way: build both, run a 128^2 blast wave, compare
the final VTK byte for byte against the patched build at the chosen settings.

  step 1  remove dead code                      (must be bit-identical)
  step 2  remove the falsified / probe switches  (must be bit-identical at their
                                                  default values)
  step 3  write in the flattener (rung L5)       (must be bit-identical to
                                                  flat1=4 hllemcap=0)

Decisions and the measurements behind them: blastwave_debug/PATCHES.md.
"""
import sys, re

def cut(s, old, tag, count=1):
    assert s.count(old) >= count, tag + ': pattern not found\n' + old[:200]
    return s.replace(old, '', count)

def sub1(s, old, new, tag):
    assert s.count(old) == 1, tag + ': expected 1 occurrence, found %d\n%s' % (s.count(old), old[:200])
    return s.replace(old, new)


def step1_dead_code(s):
    """Declared but never defined, never called, or unreachable."""
    n0 = len(s)
    for f in ('void RusanovK_x(double *QL, double *QR, double *flux);\n',
              'void RusanovK_y(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);\n',
              'void RusanovK_Pressure_x(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Pressure_y(double *QL, double *QR, double *flux);\n',
              'void RusanovK_Update_Pressure(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);\n'):
        s = cut(s, f, 'RusanovK decl')
    item3 = ("//   3. epsilon_min/epsilon_max didn't exist in this file -- added as\n"
             "//      globals, defaulted to [0,1] (no artificial floor). Note: the WENO3\n"
             "//      characteristic reconstruction's condition number was measured to\n"
             "//      blow up for eps_min below ~0.05-0.1 (see earlier conversation) --\n"
             "//      if that turns out to matter in practice, raising epsilon_min is a\n"
             "//      one-line fix here, not a structural change.\n")
    s = cut(s, item3, 'porting note 3')
    m = re.search(r'// EPSILON BOUNDS \(per author instruction\).*?double epsilon_max = 1\.0;\n', s, re.S)
    assert m, 'epsilon bounds block not found'
    s = s[:m.start()] + s[m.end():]
    assert s.count('epsilon_min') == 0 and s.count('epsilon_max') == 0, 'epsilon_min/max still referenced'
    assert s.count('Bx_amb') == 1, 'Bx_amb referenced %d times' % s.count('Bx_amb')
    s = sub1(s, 'double rho_amb = 1.0, p_amb = 0.1, Bx_amb = 70.0/sqrt(4.0*M_PI);',
                'double rho_amb = 1.0, p_amb = 0.1;', 'Bx_amb')
    s = re.sub(r'^// \[removed: [^\]]*\]\n', '', s, flags=re.M)
    print('  step 1 dead code: %d chars removed' % (n0 - len(s)))
    return s


def step2_probes(s):
    """fix I (falsified: NaN on Orszag-Tang) and fix J (written in as ON)."""
    n0 = len(s)
    # --- fix I out ---
    s = cut(s, 'double switchS(double M);   // [fix I] defined below; chi needs the split weights\n', 'fix I decl')
    s = sub1(s, 'double switchS(double M) {', 'inline double switchS(double M) {', 'switchS inline')
    m = re.search(r'[\t ]*double cf;\n[\t ]*if\(g_chisplit\)\{.*?\}else\{\n([\t ]*cf = std::sqrt\(a2 \+ cA2\);\n)[\t ]*\}\n', s, re.S)
    assert m, 'fix I cf block not found'
    s = s[:m.start()] + '\t\t\t\t\tdouble cf  = std::sqrt(a2 + cA2);\n' + s[m.end():]
    m = re.search(r'[\t ]*double denominator = shock_kappa1 \* cf_min_nbr;\n[\t ]*if\(!\(denominator > 0\.0\)\)\{[^\n]*\n[\t ]*shock_eta_tilde\(i, j\) = \(div_u < 0\.0\) \? 1\.0 : 0\.0;\n[\t ]*\}else\{\n(.*?)\n[\t ]*\}\n', s, re.S)
    assert m, 'fix I chi block not found'
    # restore exactly the pre-fix-I form; the kept body declared numerator itself
    s = (s[:m.start()]
         + '\t\t\tdouble numerator   = h * div_u + shock_kappa1 * cf_min_nbr;\n'
         + '\t\t\tdouble denominator = shock_kappa1 * cf_min_nbr;\n\t\t\t\n'
         + '\t\t\tshock_eta_tilde(i, j) = std::min(1.0,\n'
         + '\t\t\t\tstd::max(0.0, -numerator / denominator));\n'
         + s[m.end():])
    # --- fix J written in as ON, with named floors ---
    assert s.count('const double thr = g_rho_thr*in(i,j,0);   // [fix J]') == 2, 'fix J thr sites'
    s = s.replace('const double thr = g_rho_thr*in(i,j,0);   // [fix J]',
                  'const double thr = RHO_PT_FLOOR*in(i,j,0);')
    assert s.count('out(i,j,0) = g_rho_scale ? thr : in(i,j,0); ng++; } }') == 2, 'fix J out sites'
    s = s.replace('out(i,j,0) = g_rho_scale ? thr : in(i,j,0); ng++; } }',
                  'out(i,j,0) = thr; ng++; } }   // scale the increment down to the floor')
    s = sub1(s, 'if(!(Qpt(i,j,0) >= g_rho_thr*Qbar(i,j,0))){ Qpt(i,j,0) = g_rho_scale ? g_rho_thr*Qbar(i,j,0) : Qbar(i,j,0); n++; }   // [fix J]',
                'if(!(Qpt(i,j,0) >= RHO_PT_FLOOR*Qbar(i,j,0))){ Qpt(i,j,0) = RHO_PT_FLOOR*Qbar(i,j,0); n++; }', 'fix J pt guard')
    # --- the globals and the runtime keys of the removed switches ---
    s = cut(s, 'int g_flat1 = 0;       // [fix L] 1: one unified face flattener sigma for delta_k AND the envelope\n'
               'int g_jguard = 1;      // [fix K] 1: keep the two J-keyed safeguards of HLLEM_x; 0: both off\n'
               'int g_chisplit = 0;    // [fix I] 1: gate chi on sqrt(eps_M a^2 + eps_B c_A^2) instead of sqrt(a^2 + c_A^2)\n'
               'int g_rho_scale = 0;   // [fix J] 1: scale the deconvolution increment to the threshold instead of discarding it\n'
               'double g_rho_thr = 0.5;// [fix J] threshold of the density point-value guard (original: 0.5)\n',
            'I/J/K/L globals')
    s = sub1(s, 'double shock_kappa1 = 0.4;',
             '// Floor of the density point value relative to its cell average, reached by\n'
             '// scaling the deconvolution increment down to it rather than discarding the\n'
             '// increment.  It bounds the derived point velocity m_pt/rho_pt by\n'
             '// 1/RHO_PT_FLOOR times the cell-average one.\n'
             'const double RHO_PT_FLOOR = 0.5;\n'
             '// Ends of the flattener ramp in the relative density jump across a face.\n'
             'const double J0_FLAT = 0.5, J1_FLAT = 1.0;\n'
             'double shock_kappa1 = 0.4;', 'named constants')
    s = sub1(s, ' else if(key=="idp") g_idp = (int)v;\n'
                '\t\telse if(key=="chisplit") g_chisplit = (int)v; else if(key=="rhoscale") g_rho_scale = (int)v; else if(key=="rhothr") g_rho_thr = v;\n'
                '\t\telse if(key=="jguard") g_jguard = (int)v; else if(key=="flat1") g_flat1 = (int)v;', '', 'I/J/K/L keys')
    s = sub1(s, ' else if(key=="hllemoff") g_hllemoff = (int)v;', '', 'hllemoff key')
    print('  step 2 probes: %d chars removed' % (n0 - len(s)))
    return s


def step3_flattener(s):
    """Write in rung L5: one face flattener, no cap, no pressure sensor."""
    n0 = len(s)
    old_sigma = s[s.index('\tdouble g_sigma_face = 0.0, g_sigma_env = 0.0;'):
                  s.index('\t{\tconst double Jw = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);')]
    s = cut(s, old_sigma, 'fix L sigma block')
    s = sub1(s, '\t{\tconst double Jw = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);\n'
                '\t\tconst double th = g_flat1 ? g_sigma_env : fmin(1.0, (Jw - 0.5)/0.5);\n'
                '\t\tif((g_flat1 ? (th > 0.0) : (g_jguard && Jw > 0.5)) && g_env_cap > 0.0){   // [fix L]',
             '\t// ---- Flattener: revert the AP-HLLEM flux to HLLE across an under-resolved jump ----\n'
             '\t// J_flat is the relative density jump across the face: dimensionless and\n'
             '\t// scale-free, so J = O(h) in smooth flow and the whole mechanism deactivates\n'
             '\t// under refinement, while J = O(1) only across a contact the mesh does not\n'
             '\t// resolve.  theta_flat ramps over [J0_FLAT, J1_FLAT].  Two actions, both of\n'
             '\t// them the same statement -- revert to HLLE:\n'
             '\t//   delta_k  <- (1 - max(chi, theta)) delta_k   anti-diffusion off, which is\n'
             '\t//               the carbuncle cure of Dumbser & Balsara;\n'
             '\t//   (sL, sR) <- widened by theta towards the CLASSICAL Einfeldt envelope of\n'
             '\t//               eq. (52), capped by the speed that sets the time step, so the\n'
             '\t//               time step is unchanged and the asymptotic scaling is kept.\n'
             '\t// chi is Balsara\'s shock flattener and enters the anti-diffusion only: gating\n'
             '\t// the envelope on chi as well costs 8% of the shock sharpness for nothing.\n'
             '\tconst double J_flat = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);\n'
             '\tconst double theta_flat = fmin(1.0, fmax(0.0, (J_flat - J0_FLAT)/(J1_FLAT - J0_FLAT)));\n'
             '\t{\tconst double th = theta_flat;\n'
             '\t\tif(th > 0.0 && g_env_cap > 0.0){', 'envelope block')
    s = sub1(s, '\t\t\tif(g_jguard && !g_flat1){ const double J = fabs(rhoR - rhoL)/fmin(rhoL, rhoR); if(J > 0.5){ delta_C *= fmax(0.0, (1.0 - J)/0.5);   // [fix K]\n'
                '#pragma omp atomic\n'
                '\t\t\t\tg_ent_lim++; } }\n'
                '\t\t\tif(g_hllemoff || g_hllem_lo){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }\n'
                '\t\t\telse if(g_flat1){ const double fl = fmax(0.0, 1.0 - g_sigma_face); delta_C *= fl; delta_Ap *= fl; delta_Am *= fl;   // [fix L] one flattener\n'
                '#pragma omp atomic\n'
                '\t\t\t\tg_ent_lim += (g_sigma_face > 0.0); }\n'
                '\t\t\telse if(g_hllemflat){ const double fl = fmax(0.0, 1.0 - g_face_chi); delta_C *= fl; delta_Ap *= fl; delta_Am *= fl; }   // [fix C] flattener\n'
                '\t\t\tif(g_hllemcap > 0.0){ const double cp = 1.0 - g_hllemcap; delta_C = fmin(delta_C, cp); delta_Ap = fmin(delta_Ap, cp); delta_Am = fmin(delta_Am, cp); }   // [fix C] cap\n',
             '\t\t\tif(g_hllem_lo){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }\n'
             '\t\t\telse { const double fl = fmax(0.0, 1.0 - fmax(g_face_chi, theta_flat));   // revert to HLLE\n'
             '\t\t\t\tdelta_C *= fl; delta_Ap *= fl; delta_Am *= fl;\n'
             '#pragma omp atomic\n'
             '\t\t\t\tg_ent_lim += (fl < 1.0); }\n', 'delta block')
    print('  step 3 flattener: %d chars removed' % (n0 - len(s)))
    return s


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    upto = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    s = open(src).read()
    if upto >= 1: s = step1_dead_code(s)
    if upto >= 2: s = step2_probes(s)
    if upto >= 3: s = step3_flattener(s)
    open(dst, 'w').write(s)
    print('wrote', dst)
