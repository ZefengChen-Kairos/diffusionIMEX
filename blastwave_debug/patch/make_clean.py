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
    """Write in rung L9: HLLEM plus exactly ONE added term.

        delta_k <- (1 - sigma) delta_k,    sigma = max(chi, s_p)

    chi is Balsara's compression flattener (paper eq:chi) and s_p the relative
    pressure jump over the face neighbours.  That is the whole addition: revert
    the low-dissipation anti-diffusion to plain HLL where a strong shock is
    detected, which is Quirk's construction with a pressure-based switch.

    DELETED, both measured to contribute nothing and one of them actively
    harmful (blastwave_debug/PATCHES.md):
      - the J-keyed entropy limiter, delta_C *= (1 - Theta(|drho|/min rho)).
        Its own comment says it gives "plain HLL across strong contacts", and a
        contact is exactly what HLLEM's anti-diffusion is for: it fired 4540
        times on RP5, i.e. at the contact on essentially every stage of every
        step, blunting the feature the paper reports.  The pressure sensor does
        not have this defect, since a contact and a rotational discontinuity are
        both pressure-continuous.
      - the envelope widening keyed on the same jump: it never fired at all on
        the blast wave at 128^2 or 256^2 for any senv.
      - the delta_k cap: no measurable effect in any pair that isolated it.
      - my own unified sigma (fix L) and the re-gated chi (fix I, falsified).

    Verified: L9 reproduces the four-mechanism form on all seven metrics at
    128^2 (rho max 3.0056 vs 3.0058, max|grad rho| 1.0540 both, shock radius
    0.2938 both, striping power 1.0% both).
    """
    n0 = len(s)
    # --- fix L's sigma block: out ---
    old_sigma = s[s.index('\tdouble g_sigma_face = 0.0, g_sigma_env = 0.0;'):
                  s.index('\t{\tconst double Jw = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);')]
    s = cut(s, old_sigma, 'fix L sigma block')
    # --- guard 3, the envelope widening: the whole block out ---
    start = s.index('\t{\tconst double Jw = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);')
    # the block runs to the line that closes it, two levels of brace after the assignment
    endmark = '\t\t\tsL = fmin(sL, th*sL1); sR = fmax(sR, th*sR1);\n\t\t}\n\t}\n'
    end = s.index(endmark, start) + len(endmark)
    # keep the explanatory comment that precedes it out of the way: it belongs to the block
    cstart = s.rindex('\t// [W] ', 0, start) if '\t// [W] ' in s[:start] else start
    s = s[:cstart] + s[end:]
    # --- the delta block: one flattener, unconditional ---
    s = sub1(s, '\t\t\tif(g_jguard && !g_flat1){ const double J = fabs(rhoR - rhoL)/fmin(rhoL, rhoR); if(J > 0.5){ delta_C *= fmax(0.0, (1.0 - J)/0.5);   // [fix K]\n'
                '#pragma omp atomic\n'
                '\t\t\t\tg_ent_lim++; } }\n'
                '\t\t\tif(g_hllemoff || g_hllem_lo){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }\n'
                '\t\t\telse if(g_flat1){ const double fl = fmax(0.0, 1.0 - g_sigma_face); delta_C *= fl; delta_Ap *= fl; delta_Am *= fl;   // [fix L] one flattener\n'
                '#pragma omp atomic\n'
                '\t\t\t\tg_ent_lim += (g_sigma_face > 0.0); }\n'
                '\t\t\telse if(g_hllemflat){ const double fl = fmax(0.0, 1.0 - g_face_chi); delta_C *= fl; delta_Ap *= fl; delta_Am *= fl; }   // [fix C] flattener\n'
                '\t\t\tif(g_hllemcap > 0.0){ const double cp = 1.0 - g_hllemcap; delta_C = fmin(delta_C, cp); delta_Ap = fmin(delta_Ap, cp); delta_Am = fmin(delta_Am, cp); }   // [fix C] cap\n',
             '\t\t\t// Revert the anti-diffusion to plain HLL where a strong shock is detected:\n'
             '\t\t\t// sigma = max(chi, s_p), chi the compression flattener of eq.(72) and s_p the\n'
             '\t\t\t// relative pressure jump over the face neighbours.  A contact and a rotational\n'
             '\t\t\t// discontinuity are both pressure-continuous, so this leaves exactly the fields\n'
             '\t\t\t// the anti-diffusion is there to sharpen.\n'
             '\t\t\tif(g_hllem_lo){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }\n'
             '\t\t\telse { const double fl = fmax(0.0, 1.0 - g_face_chi);\n'
             '\t\t\t\tdelta_C *= fl; delta_Ap *= fl; delta_Am *= fl;\n'
             '#pragma omp atomic\n'
             '\t\t\t\tg_ent_lim += (fl < 1.0); }\n', 'delta block')
    # the sigma, the guard and the widening are gone from the flux; g_env_cap and
    # g_env_wide are now dead (declaration and assignments only) and go in 4b
    for v in ('g_sigma_face', 'g_sigma_env', 'g_flat1', 'g_jguard'):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    assert 'g_env_wide++' not in s, 'envelope widening counter still incremented'
    assert '-g_env_cap' not in s, 'envelope clamp still present'
    print('  step 3 flattener (L9): %d chars removed' % (n0 - len(s)))
    return s


def step4a_off_switches(s):
    """Three switches that default to OFF.  Removing an off branch must be
    bit-identical, which is how this sub-step is checked."""
    n0 = len(s)
    # --- fix B: pressure-fluctuation-aware switch ---
    m = re.search(r'[\t ]*if\(g_pswitch\)\{   // \[fix B\].*?\n[\t ]*cM\(i,j\) = fmin\(eps_M_max[^\n]*\n[\t ]*\}\n', s, re.S)
    assert m, 'pswitch branch not found'
    s = s[:m.start()] + s[m.end():]
    s = cut(s, 'int g_pswitch = 0;   // 1: switching function driven by max(M, sqrt(|dp|/p)) and max(M_A, sqrt(|dm|/m))\n', 'pswitch global')
    # --- the optional 3x3 widening of eps ---
    m = re.search(r'[\t ]*if\(g_epsmax3\)\{   // optional:.*?\{ cM\(i,j\) = cM3\(i,j\); cB\(i,j\) = cB3\(i,j\); \}\n[\t ]*\}\n', s, re.S)
    assert m, 'epsmax3 branch not found'
    s = s[:m.start()] + s[m.end():]
    s = cut(s, 'int g_epsmax3 = 0;   // 1: cell eps := max over the 3x3 neighbourhood before the face max\n', 'epsmax3 global')
    # --- fix D: the plasma-beta floor, never switched on in any reference run ---
    m = re.search(r'[\t ]*if\(g_betafloor > 0\.0\)\{   // \[fix D\].*?\n[\t ]*\}\n[\t ]*\}\n', s, re.S)
    assert m, 'betafloor block not found'
    s = s[:m.start()] + s[m.end():]
    s = cut(s, 'double g_betafloor = 0.0;          // [fix D] minimum plasma beta: p >= betafloor * |B|^2/2 (betafloor=0 -> off)\n', 'betafloor global')
    # --- their runtime keys and the ARGS banner ---
    s = sub1(s, ' else if(key=="epsmax3") g_epsmax3 = (int)v; else if(key=="pswitch") g_pswitch = (int)v;', '', 'epsmax3/pswitch keys')
    s = sub1(s, 'else if(key=="betafloor") g_betafloor = v; ', '', 'betafloor key')
    s = sub1(s, ' << " epsmax3=" << g_epsmax3 << " pswitch=" << g_pswitch', '', 'epsmax3/pswitch banner')
    s = sub1(s, ' << " betafloor=" << g_betafloor', '', 'betafloor banner')
    for v in ('g_pswitch', 'g_epsmax3', 'g_betafloor'):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    print('  step 4a off switches: %d chars removed' % (n0 - len(s)))
    return s


def step4b_probe_machinery(s):
    """fix G (the IDP probe) out entirely, and the flattener machinery made
    unconditional: hllemflat / hllemcap / hllemoff were switches over code that
    is now the only path."""
    n0 = len(s)
    # --- fix G: declarations, grids, both second-solver calls, the U^L block ---
    s = cut(s, 'int g_idp = 0;                               // [fix G] 0 off / 1 compute-only (inert)\n'
               'static thread_local int g_hllem_lo = 0;      // [fix G] 1: plain HLL, for the low-order flux\n'
               'long g_idp_ulbad_r = 0, g_idp_ulbad_p = 0;   // [fix G] cells whose U^L is itself inadmissible, split by rho / p\n',
            'fix G globals')
    s = cut(s, '\tstatic Grid FxL(nxCells, nyCells), FyL(nxCells, nyCells), Ulo(nxCells, nyCells);   // [fix G]\n', 'fix G grids')
    pat = re.compile(r'[\t ]*if\(g_idp && g_hllem_Qe\)\{   // \[fix G\][^\n]*\n(?:[^\n]*\n)*?[\t ]*for\(int k = 0; k < NUM_VARS; k\+\+\) F[xy]L\(i,j,k\) = fl\[k\];\n[\t ]*\}\n')
    assert len(pat.findall(s)) == 2, 'fix G solver blocks: found %d' % len(pat.findall(s))
    s = pat.sub('', s)
    m = re.search(r'[\t ]*if\(g_idp\)\{   // \[fix G\] the first-order update.*?\n[\t ]*g_idp_ulbad_r \+= nbadr; g_idp_ulbad_p \+= nbadp;\n[\t ]*\}\n', s, re.S)
    assert m, 'fix G U^L block not found'
    s = s[:m.start()] + s[m.end():]
    s = sub1(s, '\t\t\tif(g_hllem_lo){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }\n\t\t\telse { const double fl',
                '\t\t\t{ const double fl', 'hllem_lo branch')
    s = sub1(s, '", U^L inadmissible: rho in " << g_idp_ulbad_r << " cells, p in " << g_idp_ulbad_p << " cells, density point-value guard applied "',
                '", density point-value guard applied "', 'U^L report')
    s = sub1(s, ' << " idp=" << g_idp', '', 'idp banner')
    # --- the flattener machinery is now the only path ---
    m = re.search(r'int g_hllemflat = [^\n]*g_hllemoff = 0;[^\n]*\n', s)
    assert m, 'hllemflat global not found'
    s = s[:m.start()] + s[m.end():]
    s = sub1(s, 'if(g_hllemflat && g_hllem_Qe){', 'if(g_hllem_Qe){', 'hchi guard')
    s = sub1(s, 'if(g_hllemflat >= 2){   // + pressure-jump sensor',
                '{   // pressure-jump sensor', 'sp guard')
    assert s.count('if(g_hllemflat){ const int ') == 2, 'face chi guards'
    s = s.replace('if(g_hllemflat){ const int ', '{ const int ')
    s = sub1(s, ' else if(key=="hllemflat") g_hllemflat = (int)v; else if(key=="hllemcap") g_hllemcap = v;', '', 'hllemflat/cap keys')
    s = sub1(s, ' << " hllemflat=" << g_hllemflat << " hllemoff=" << g_hllemoff << " hllemcap=" << g_hllemcap', '', 'hllemflat banner')
    m = re.search(r'double g_hllemcap = [^\n]*\n', s)
    assert m, 'hllemcap global not found'
    s = s[:m.start()] + s[m.end():]
    for v in ('g_idp', 'g_hllem_lo', 'g_idp_ulbad', 'FxL', 'FyL', 'Ulo', 'g_hllemflat', 'g_hllemcap', 'g_hllemoff'):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    print('  step 4b probe machinery: %d chars removed' % (n0 - len(s)))
    return s


def step4c_bake_method_params(s):
    """The method parameters the paper states.  Writing them in leaves one path;
    the defaults already are the chosen values, so this must be bit-identical.

      kappa1 = 0.4            paper eq.(72)
      M0, M1, alpha           paper eq:switch, 0.15 / 0.4 / 14
      plap = 1                the shock-activated pressure Laplacian, eq.(45)/(85)
      bsmooth = 1             the erf-smoothed blast jump, which is what is run
    """
    n0 = len(s)
    s = sub1(s, 'else if(key=="plap") g_plap = (int)v; else if(key=="tag") g_tag = k.substr(eq+1); else if(key=="kappa1") shock_kappa1 = v;',
                'else if(key=="tag") g_tag = k.substr(eq+1);', 'plap/kappa1 keys')
    s = sub1(s, '\t\telse if(key=="swm0") g_sw_m0 = v; else if(key=="swm1") g_sw_m1 = v; else if(key=="swalpha") g_sw_alpha = v;\n', '', 'switch keys')
    s = sub1(s, 'else if(key=="bsmooth") g_blast_smooth = (int)v; ', '', 'bsmooth key')
    # Deliberately NOT const.  Marking them const lets -O2 -march=native fold
    # the literal into the expression, which enables FMA contraction the
    # runtime-loaded double did not get, and the last bits of the solution move.
    # Removing the key already leaves one path; const would only be cosmetic and
    # would cost the bit-identity check.  Measured: with const, the 128^2 blast
    # diverges from step 2 on.
    s = sub1(s, 'double shock_kappa1 = 0.4;    // [fix H] paper eq.(72);',
                'double shock_kappa1 = 0.4;    // paper eq.(72);', 'kappa1 comment')
    s = sub1(s, '// switch parameters (paper: 0.15, 0.4, 14)', '// paper eq:switch', 'switch comment')
    s = sub1(s, 'int g_plap = 1;      // [fix] 1 = shock-activated pressure Laplacian of paper 1 (eq. 45/85), (1-eps_M)-weighted, in L2 and L4',
                'int g_plap = 1;      // shock-activated pressure Laplacian, eq.(45)/(85), (1-eps_M)-weighted, in L2 and L4',
             'plap comment')
    # the ARGS banner now only echoes what can still be chosen
    m = re.search(r'\tstd::cout << "ARGS cfl=" << CFL[^\n]*\n[^\n]*\n', s)
    assert m, 'ARGS banner not found'
    s = s[:m.start()] + '\tstd::cout << "ARGS cfl=" << CFL << " senv=" << g_senv\n' + s[m.end():]
    print('  step 4c method params: %d chars removed' % (n0 - len(s)))
    return s


def step4d_dead_envelope(s):
    """g_env_cap and g_env_wide served only the envelope widening, which step 3
    deleted.  They are now written and never read."""
    n0 = len(s)
    s = cut(s, 'double g_env_cap = 0.0;   // [W] time-step speed of the current sweep direction (a_max_x / a_max_y), set by HLLEM_Update\n'
               'long g_env_wide = 0;      // [W] counter: faces whose envelope was widened\n', 'env globals')
    assert s.count('\tg_env_cap = a_max_x;   // [W]\n') == 1 and s.count('\tg_env_cap = a_max_y;   // [W]\n') == 1, 'env_cap assignments'
    s = s.replace('\tg_env_cap = a_max_x;   // [W]\n', '').replace('\tg_env_cap = a_max_y;   // [W]\n', '')
    s = sub1(s, ' << " faces, envelope widened at " << g_env_wide << " faces,', ' << " faces,', 'env_wide report')
    m = re.search(r'//[^\n]*the local \[W\] widening below is capped by g_env_cap and\n[^\n]*\n', s)
    if m: s = s[:m.start()] + s[m.end():]
    for v in ('g_env_cap', 'g_env_wide'):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    print('  step 4d dead envelope vars: %d chars removed' % (n0 - len(s)))
    return s


def step5_comments(s):
    """Prune the development archaeology from the comments.

    Deleting a comment cannot change behaviour, so this step is checked the same
    way as the others and the check is a check of the SCRIPT, not of the code:
    if a replacement reaches into code by accident, the build or the bit-identity
    check fails.

    Removed: the patch log at the top of the file (those fixes are now simply the
    code), the bracketed development tags, references to files a reader does not
    have, and an empty if-block with a TODO in it.  KEPT: every comment that
    explains why the code is as it is, every equation reference, and every
    warning about what breaks if it is changed.
    """
    n0 = len(s)
    # --- the patch log: from the first separator to just before the real header ---
    hdr = '// ============================================================================\n' \
          '// Dual-parameter (eps_M, eps_B) IMEX finite volume scheme for ideal MHD'
    i = s.index(hdr)
    assert s[:i].count('// PATCHED VERSION') == 1, 'patch log not where expected'
    s = s[i:]
    # --- development tags, stripped from the comment text only ---
    tagpat = re.compile(r'\[(?:fix|o3)[^\]]{0,24}\]\s*|\[(?:W|emb|1d|2d|vtk)\]\s*')
    out = []
    for line in s.split('\n'):
        k = line.find('//')
        if k >= 0:
            head, tail = line[:k], line[k:]
            tail = tagpat.sub('', tail)
            tail = re.sub(r'//\s+(?=\w)', '// ', tail)
            line = head + tail
        out.append(line)
    s = '\n'.join(out)
    # --- a changelog opening on an otherwise useful explanation ---
    s = sub1(s, 'The three switches are GONE: the corrections are always on,\n',
                'The corrections are always on,\n', 'switches-gone wording')
    # --- references to files the reader does not have ---
    s = re.sub(r'[ \t]*(?://)?[^\n]*\b(?:REPORT\.md|apply_patch\.py|blastwave_debug/)[^\n]*\n',
               lambda m: '' if m.group(0).lstrip().startswith('//') else m.group(0), s)
    # --- an empty if-block holding only a TODO ---
    m = re.search(r'[\t ]*if\(nyCells != 1 && nxCells != 1\)\{\n[\t ]*// TODO-2D[^\n]*\n[\t ]*\}\n', s)
    if m: s = s[:m.start()] + s[m.end():]
    # --- the method header now names the flattener as well ---
    s = sub1(s, '// SHOCKS shock-indicator-gated Laplacian (1 - eps_B) |lambda_B| dx in the magnetic predictor.',
                '// SHOCKS shock-indicator-gated Laplacian (1 - eps_B) |lambda_B| dx in the magnetic predictor;\n'
                '//        the HLLEM anti-diffusion is reverted to plain HLL where a strong shock is detected,\n'
                '//        delta_k <- (1 - max(chi, s_p)) delta_k with s_p the relative pressure jump.',
             'SHOCKS header')
    assert '[fix' not in s and '[o3' not in s, 'tags remain'
    print('  step 5 comments: %d chars removed' % (n0 - len(s)))
    return s


def step6_junk_comments(s):
    """Commented-out code and archaeology framing.  Keeps every comment that
    explains the numerics; removes the ones that carry no information."""
    n0 = len(s)
    # --- the dead dissipation-selection switch, commented out wholesale ---
    m = re.search(r'[\t ]*//[\t ]*\n(?:[\t ]*//[^\n]*\n)*?[\t ]*//[\t ]*//Laplacian Operator\n'
                  r'(?:[\t ]*//[^\n]*\n)*?(?=[\t ]*auto henthA)', s)
    assert m, 'dead flux_diss block not found'
    s = s[:m.start()] + s[m.end():]
    # --- runs of empty comment lines (two or more): pure noise ---
    s = re.sub(r'(?:[\t ]*//[\t ]*\n){2,}', '', s)
    # --- archaeology framing on otherwise useful blocks ---
    s = sub1(s, '// ==================== RESTORED FROM ARCHIVE (positivity core) ====================\n', '',
             'restored-from-archive')
    s = sub1(s, 'BLOCK-STRUCTURED SPLIT CHARACTERISTIC BASIS  (replaces compute_LR_prim_x/y and\n'
                '\t\t// its cond>1e4 -> IDENTITY fallback, which silently disabled characteristic\n'
                '\t\t// reconstruction wherever the 8x8 R was ill-conditioned -- i.e. over the whole\n'
                '\t\t// Bperp~0 region, e.g. the entire blast-wave centre.)',
                'BLOCK-STRUCTURED SPLIT CHARACTERISTIC BASIS.  A numerical 8x8 inverse is\n'
                '\t\t// ill-conditioned wherever Bperp -> 0, e.g. over the whole blast-wave centre,\n'
                '\t\t// so the basis is built block by block instead, analytically where possible.',
             'basis archaeology')
    # --- references to internal documents the reader does not have ---
    s = re.sub(r'[ \t]*//[^\n]*\bDocument [123]\b[^\n]*\n', '', s)
    s = re.sub(r'\(Document [123][^)]*\)\s*', '', s)
    s = re.sub(r'（?设计\(ROUTE A-analytic, \d{4}-\d\d\)', '设计', s)
    # --- rejected-alternative narrative that documents nothing about THIS code ---
    m = re.search(r'[\t ]*// Previously the faces were built from A_s:.*?\n(?=[\t ]*(?:void|auto|static|const|double|int|//\s*-|\n))', s, re.S)
    if m: s = s[:m.start()] + s[m.end():]
    print('  step 6 junk comments: %d chars removed' % (n0 - len(s)))
    return s


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    upto = sys.argv[3] if len(sys.argv) > 3 else '6'
    order = ['1', '2', '3', '4a', '4b', '4c', '4d', '5', '6']
    want = order[:order.index(upto) + 1]
    s = open(src).read()
    if '1' in want:  s = step1_dead_code(s)
    if '2' in want:  s = step2_probes(s)
    if '3' in want:  s = step3_flattener(s)
    if '4a' in want: s = step4a_off_switches(s)
    if '4b' in want: s = step4b_probe_machinery(s)
    if '4c' in want: s = step4c_bake_method_params(s)
    if '4d' in want: s = step4d_dead_envelope(s)
    if '5' in want:  s = step5_comments(s)
    if '6' in want:  s = step6_junk_comments(s)
    open(dst, 'w').write(s)
    print('wrote', dst)
