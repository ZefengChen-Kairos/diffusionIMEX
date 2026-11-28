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
    # fix J REVERTED: the original and the patched default both DISCARD the increment
    # (rho_pt <- rho_bar).  Baking in "scale to the floor" was a deviation no test could
    # detect, since the guard fires zero times in every test run, and it disagreed with the
    # second guard in fv_cell_lim_to, which keeps the cell average.
    s = s.replace('out(i,j,0) = g_rho_scale ? thr : in(i,j,0); ng++; } }',
                  'out(i,j,0) = in(i,j,0); ng++; } }   // discard the deconvolution increment')
    s = sub1(s, 'if(!(Qpt(i,j,0) >= g_rho_thr*Qbar(i,j,0))){ Qpt(i,j,0) = g_rho_scale ? g_rho_thr*Qbar(i,j,0) : Qbar(i,j,0); n++; }   // [fix J]',
                'if(!(Qpt(i,j,0) >= RHO_PT_FLOOR*Qbar(i,j,0))){ Qpt(i,j,0) = Qbar(i,j,0); n++; }', 'fix J pt guard')
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
    # the hllemflat=3 isolation hook (s_p alone, to measure chi's contribution) collapses to L9
    s = sub1(s,
        '\t\t\tfor(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++)\n'
        '\t\t\t\thchi(i,j) = (g_hllemflat >= 3) ? sp(i,j) : fmax(hchi(i,j), sp(i,j));   // 3: s_p ALONE, to isolate chi\'s contribution\n',
        '\t\t\tfor(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++) hchi(i,j) = fmax(hchi(i,j), sp(i,j));\n',
        'sp merge (drop the hllemflat=3 isolation hook)')
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
    s = s[:m.start()] + '\tstd::cout << "ARGS cfl=" << CFL\n' + s[m.end():]
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


def step4e_senv(s):
    """Fix the envelope lower bound to the material speed of the interface and
    drop the switch.

    This is the one step that is NOT bit-identical to the previous rung by
    construction: the old default was a three-way option.  The reference is the
    patched source run with senv=2.

    Why the bound is needed: as eps_M, eps_B -> 0 the magnetosonic quartic tends
    to lambda^2 (lambda - u)^2 (SplitMagnetosonicCoeffs at eps = 0 gives
    c3 = -2u, c2 = u^2, c1 = c0 = 0), so its roots tend to {0, 0, u, u} and the
    envelope to [min(0,u), max(0,u)].  Where u keeps one sign across the
    interface one endpoint is EXACTLY zero, s_L s_R = 0, and the HLL part carries
    no dissipation at all -- not only at isolated stagnation points.  Without
    the bound the scheme converges to a WRONG WEAK SOLUTION -- measured on RP6:
    post-shock rho 3.884 against the exact 3.982, shock at -0.383 against -0.371,
    unchanged from 250 to 2000 cells.

    Why local and not domain-global: both close the degeneracy, but a global bound
    drives delta_ent = 1 - |u|/u_max -> 1 wherever |u| << max_Omega|u|, which
    REMOVES the dissipation from the three linearly degenerate fields in the quiet
    ambient of a large-dynamic-range flow.  Measured on the 128^2 blast wave,
    cells that are a local extremum of relative amplitude > 5%: rho 0.10% local,
    1.03% global; > 1%: 1.32% local, 6.30% global, 2.20% with no bound.  The local
    bound is also the Rusanov speed of paper eq:rusanov in the stiff limit, so it
    needs no citation of its own.

    Paper: eq:einfeldt carries the bound, eq:umaxface defines it, and the AP
    paragraph derives the stiff limit from it (eq:aplimit).

    BIT-IDENTITY: this step is semantically a no-op at senv=2 -- the branch that
    ran there becomes unconditional with the same arithmetic, and the deleted
    g_senv_u writes were never read -- but it is NOT bit-identical under
    -O3 -march=native.  Removing the branch lets the compiler contract the
    fmax/fabs chain into FMAs it did not use before, so the last bits move
    (measured after 4 steps of the 128^2 blast: rho 5e-12 relative, eps_M still
    exact; after 618 steps that has grown to 5e-4, which is ordinary chaotic
    amplification, not a semantic change).  Verified instead with
    -ffp-contract=off, under which the clean source and the patched source run
    with senv=2 agree bit for bit over the whole run.  Same family as the const
    note in step 4c: do not try to "fix" this by reintroducing a branch.
    """
    n0 = len(s)
    old = ('\tif(g_senv == 1 && g_senv_u > 0.0){ sL = fmin(sL, -g_senv_u); sR = fmax(sR, g_senv_u); }   // [fix E] global |u| floor\n'
           '\telse if(g_senv == 2){   // [fix E2] LOCAL material-speed floor: max|u_n| over {L, Roe, R} of THIS face\n'
           '\t\tconst double um = fmax(fabs(vxL), fmax(fabs(vx_t), fabs(vxR)));\n'
           '\t\tif(um > 0.0){ sL = fmin(sL, -um); sR = fmax(sR, um); }\n'
           '\t}\n')
    new = ('\t// Lower bound on the envelope, paper eq:einfeldt with eq:umaxface: the largest material\n'
           '\t// speed over the same three states {L, Roe, R}; since the Roe state is a convex combination\n'
           '\t// this equals max(|u_L|, |u_R|), the Rusanov speed of paper 1 eq.(25).  As eps_M, eps_B -> 0\n'
           '\t// the magnetosonic quartic tends to lambda^2 (lambda - u)^2, its roots to {0, 0, u, u} and\n'
           '\t// the envelope to [min(0,u), max(0,u)]: where u keeps one sign across the interface one\n'
           '\t// endpoint is EXACTLY zero, s_L s_R = 0, and the HLL part carries no dissipation at all.\n'
           '\t// The scheme then converges to a wrong weak solution (RP6 without the bound:\n'
           '\t// post-shock rho 3.884 against the exact 3.982, shock -0.383 against -0.371, unchanged\n'
           '\t// from 250 to 2000 cells).\n'
           '\t// It is taken over the states of THIS interface: a domain-global bound closes the same\n'
           '\t// degeneracy but drives delta_ent = 1 - |u|/u_max -> 1 wherever |u| << max|u|, removing\n'
           '\t// the dissipation from the linearly degenerate fields in the quiet ambient of a flow with\n'
           '\t// a large dynamic range (128^2 blast, cells that are a local extremum of amplitude > 5%:\n'
           '\t// rho 0.10% with this bound, 1.03% with the global one).\n'
           '\t{\tconst double um = fmax(fabs(vxL), fmax(fabs(vx_t), fabs(vxR)));\n'
           '\t\tif(um > 0.0){ sL = fmin(sL, -um); sR = fmax(sR, um); }\n'
           '\t}\n')
    s = sub1(s, old, new, 'senv floor -> unconditional interface-local')
    m2 = re.search(r'// \[fix E\] HLL wave-speed envelope:.*?\nint g_senv = 2;\n', s, re.S)
    assert m2, 'senv comment block + declaration not found'
    s = s[:m2.start()] + s[m2.end():]
    s = sub1(s, '\t\telse if(key=="senv") g_senv = (int)v;\n', '', 'senv key')
    # g_senv_u served only the global branch and is now written and never read
    s = cut(s, 'double g_senv_u = 0.0;             // global max|u| of the current sweep direction\n', 'g_senv_u declaration')
    s = sub1(s, '\tg_senv_u  = vx_max;    // [fix E]\n', '', 'g_senv_u x assignment')
    s = sub1(s, '\tg_senv_u  = vy_max;    // [fix E]\n', '', 'g_senv_u y assignment')
    assert s.count('g_senv') == 0, 'g_senv still referenced %d times' % s.count('g_senv')
    print('  step 4e envelope bound: %d chars removed' % (n0 - len(s)))
    return s


def step4f_drop_conversion_limiter(s):
    """Delete the minmod limiter from the average<->point / transverse conversions
    and from the CT slopes of grad A, and with it the whole eps_M blending.

    What goes: d2lim and mm3, the two kappa = 0.5 constants, the two global
    grid_absmax / reduction scales that fed their thresholds, the eps_M blend
    weights w and wq, the globals g_epsM_cell / g_epsM_fx / g_epsM_fy that exist
    only to carry those weights, and the CT lambdas mm2g / sl with S_gx, S_gy.
    What stays: the conversions themselves, which become paper eq:convdeconv
    verbatim, and the positivity guard of fv_cell_lim_to (its "redo fully
    limited" retry is dropped, since there is no limited branch to redo).

    Measured, limiter on vs off:
      low-Mach vortex rho_0 = 1e-2 (eps_M ~ 0): bit-identical at 32, 64, 128 --
        the eps_M blend had already switched the limiter off there.
      vortex rho_0 = 1 (eps_M ~ 1): third order either way; the 64->128 rates are
        3.13/2.94/2.94/3.01/3.00 with and 3.03/2.95/2.94/3.02/3.02 without, the
        32->64 rate of rho is better without (2.67 vs 2.48), and at N = 128 the
        two are BIT-IDENTICAL -- the threshold carries one power of h less than
        the quantity it gates, so the gate opens everywhere under refinement.
      blast 128^2: cells that are a local extremum above 5% 0.101 either way,
        above 1% 1.316 with and 1.373 without.
      RP5 800x1: 9 cells above 1% with (median 1.92%), 8 without (median 2.18%).
    """
    n0 = len(s)
    # --- fv_cell_lim_to: drop the scale, kappa, w, the two-pass retry and the blend ---
    pat = ('\tdouble S[NUM_VARS]; grid_absmax(tmp, ghost_cells, nxCells+ghost_cells, ghost_cells, nyCells+ghost_cells, S);\n'
           '\tconst double kappa = 0.5;\n')
    assert s.count(pat) == 2, 'scale + kappa: found %d' % s.count(pat)
    s = s.replace(pat, '')   # fv_cell_lim_to and fv_transverse_lim
    s = sub1(s,
        '\t\t\tdouble w = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(std::max(ghost_cells,std::min(nxCells+ghost_cells-1,i)), std::max(ghost_cells,std::min(nyCells+ghost_cells-1,j))))) : 1.0;\n'
        '\t\t\tfor(int pass = 0; pass < 2; pass++){\n'
        '\t\t\t\tfor(int k = 0; k < NUM_VARS; k++){\n'
        '\t\t\t\t\tdouble d2 = w * d2lim(tmp(i-2,j,k),tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k),tmp(i+2,j,k), kappa*dx*S[k]) + (1.0 - w) * d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k));\n'
        '\t\t\t\t\td2 += w * d2lim(tmp(i,j-2,k),tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k),tmp(i,j+2,k), kappa*dy*S[k]) + (1.0 - w) * d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k));\n'
        '\t\t\t\t\tout(i,j,k) = tmp(i,j,k) + sgn * d2 / 24.0;\n'
        '\t\t\t\t}\n',
        '\t\t\t{\n'
        '\t\t\t\tfor(int k = 0; k < NUM_VARS; k++){\n'
        '\t\t\t\t\tconst double d2 = d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k)) + d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k));\n'
        '\t\t\t\t\tout(i,j,k) = tmp(i,j,k) + sgn * d2 / 24.0;\n'
        '\t\t\t\t}\n',
        'fv_cell_lim_to blend + pass loop')
    s = sub1(s, '\t\t\t\t\t\tif(pass == 0 && w < 1.0){ w = 1.0; continue; }\n', '', 'fv_cell_lim_to retry')
    # the pass loop's trailing `break` would now bind to the enclosing i-loop and abandon the
    # rest of the row (measured: NaN from step 1, density guard firing 1.4e8 times)
    s = sub1(s,
        '\t\t\t\t\t\tfor(int k = 0; k < NUM_VARS; k++) out(i,j,k) = tmp(i,j,k);\n'
        '\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\tbreak;\n\t\t\t}\n',
        '\t\t\t\t\t\tfor(int k = 0; k < NUM_VARS; k++) out(i,j,k) = tmp(i,j,k);\n'
        '\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n',
        'fv_cell_lim_to pass-loop break')
    # --- fv_cell_lim_s ---
    s = sub1(s, '\tdouble S = 0.0;\n'
                '#pragma omp parallel for collapse(2) schedule(static) reduction(max:S)\n'
                '\tfor(int j = ghost_cells; j < nyCells+ghost_cells; j++) for(int i = ghost_cells; i < nxCells+ghost_cells; i++) S = fmax(S, fabs(in(i,j)));\n'
                '\tconst double kappa = 0.5;\n', '', 'fv_cell_lim_s scale + kappa')
    s = sub1(s,
        '\t\t\tdouble w = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(std::max(ghost_cells,std::min(nxCells+ghost_cells-1,i)), std::max(ghost_cells,std::min(nyCells+ghost_cells-1,j))))) : 1.0;\n'
        '\t\t\tdouble d2 = w * d2lim(in(i-2,j),in(i-1,j),in(i,j),in(i+1,j),in(i+2,j), kappa*dx*S) + (1.0 - w) * d2lin(in(i-1,j),in(i,j),in(i+1,j));\n'
        '\t\t\td2 += w * d2lim(in(i,j-2),in(i,j-1),in(i,j),in(i,j+1),in(i,j+2), kappa*dy*S) + (1.0 - w) * d2lin(in(i,j-1),in(i,j),in(i,j+1));\n',
        '\t\t\tdouble d2 = d2lin(in(i-1,j),in(i,j),in(i+1,j)) + d2lin(in(i,j-1),in(i,j),in(i,j+1));\n',
        'fv_cell_lim_s blend')
    # --- fv_transverse_lim ---
    s = sub1(s,
        "\t\t\tArray* ef = (dir == 'x') ? g_epsM_fx : g_epsM_fy;\n"
        "\t\t\tconst int fi = (dir == 'x') ? i - off : i, fj = (dir == 'y') ? j - off : j;   // [o3q] eps of the face this state belongs to\n"
        '\t\t\tconst double w = ef ? fmin(1.0, fmax(0.0, (*ef)(std::max(ghost_cells-1,std::min(nxCells+ghost_cells-1,fi)), std::max(ghost_cells-1,std::min(nyCells+ghost_cells-1,fj))))) : 1.0;\n'
        '\t\t\tfor(int k = 0; k < NUM_VARS; k++){\n'
        "\t\t\t\tconst double thr = kappa * ((dir == 'x') ? dy : dx) * S[k];\n"
        "\t\t\t\tdouble d2 = (dir == 'x')\n"
        '\t\t\t\t          ? w * d2lim(tmp(i,j-2,k),tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k),tmp(i,j+2,k), thr) + (1.0 - w) * d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k))\n'
        '\t\t\t\t          : w * d2lim(tmp(i-2,j,k),tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k),tmp(i+2,j,k), thr) + (1.0 - w) * d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k));\n',
        '\t\t\tfor(int k = 0; k < NUM_VARS; k++){\n'
        "\t\t\t\tconst double d2 = (dir == 'x') ? d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k))\n"
        '\t\t\t\t                               : d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k));\n',
        'fv_transverse_lim blend')
    # --- CT: the same limiter on the slopes of G = grad A ---
    s = sub1(s, '\tauto mm2g = [&](double a, double b, double thr)->double{ if(fabs(a) < thr && fabs(b) < thr) return 0.5*(a+b); return (a*b > 0.0) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0; };\n',
                '', 'CT mm2g')
    s = sub1(s, '\t\tdouble S_gx = 0.0, S_gy = 0.0;\n'
                '#pragma omp parallel for collapse(2) schedule(static) reduction(max:S_gx,S_gy)\n'
                '\t\tfor(int j = ghost_cells; j < nyCells+ghost_cells; j++) for(int i = ghost_cells; i < nxCells+ghost_cells; i++){ S_gx = fmax(S_gx, fabs(Gx(i,j))); S_gy = fmax(S_gy, fabs(Gy(i,j))); }\n',
                '', 'CT global slope scales')
    s = sub1(s, '\t\t\t\tdouble wq = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(i,j))) : 1.0;\n'
                '\t\t\t\tauto sl = [&](double a, double b, double thr)->double{ return wq * mm2g(a, b, thr) + (1.0 - wq) * 0.5*(a+b); };\n'
                '\t\t\t\t// [o3t] the slopes of G = grad A use the SAME eps_M-blended, threshold-gated limiter as the slopes of u:\n'
                '\t\t\t\t// the unconditional minmod zeroed them at every smooth extremum of G (O(h^2) in h-wide bands, B residual\n'
                '\t\t\t\t// ~2.5th order on the stationary vortex); the gate still limits across a jump (|dG| > S_G/2).\n'
                '\t\t\t\tSxx(i,j) = sl(Gx(i+1,j)-Gx(i,j), Gx(i,j)-Gx(i-1,j), 0.5*S_gx)/dx; Sxy(i,j) = sl(Gx(i,j+1)-Gx(i,j), Gx(i,j)-Gx(i,j-1), 0.5*S_gx)/dy;\n'
                '\t\t\t\tSyx(i,j) = sl(Gy(i+1,j)-Gy(i,j), Gy(i,j)-Gy(i-1,j), 0.5*S_gy)/dx; Syy(i,j) = sl(Gy(i,j+1)-Gy(i,j), Gy(i,j)-Gy(i,j-1), 0.5*S_gy)/dy;\n',
                '\t\t\t\tauto sl = [](double a, double b)->double{ return 0.5*(a+b); };   // central slope of G, unlimited (paper eq:Aupdate)\n'
                '\t\t\t\tSxx(i,j) = sl(Gx(i+1,j)-Gx(i,j), Gx(i,j)-Gx(i-1,j))/dx; Sxy(i,j) = sl(Gx(i,j+1)-Gx(i,j), Gx(i,j)-Gx(i,j-1))/dy;\n'
                '\t\t\t\tSyx(i,j) = sl(Gy(i+1,j)-Gy(i,j), Gy(i,j)-Gy(i-1,j))/dx; Syy(i,j) = sl(Gy(i,j+1)-Gy(i,j), Gy(i,j)-Gy(i,j-1))/dy;\n',
                'CT slope limiter -> central')
    # --- the limiter helpers and the eps_M plumbing that only fed the blend weights ---
    m = re.search(r'static inline double mm3\(double a, double b, double c\)\{.*?\n\}\n', s, re.S)
    assert m, 'mm3 not found'
    s = s[:m.start()] + s[m.end():]
    m = re.search(r'static inline double d2lim\(double qm2[^\n]*\{\n(?:[^\n]*\n)*?\}\n', s)
    assert m, 'd2lim not found'
    s = s[:m.start()] + s[m.end():]
    # the rationale comments for the gate and the blend go with them
    m = re.search(r'// limited undivided second difference at the centre point of a 5-point line\n'
                  r'(?:// \[o3b\][^\n]*\n|//[^\n]*\n)*?'
                  r'Array \*g_epsM_cell = nullptr, \*g_epsM_fx = nullptr, \*g_epsM_fy = nullptr;\n', s)
    assert m, 'gate/blend rationale block not found'
    s = s[:m.start()] + s[m.end():]
    s = sub1(s, '\tg_epsM_cell = &eps_M;\n', '', 'g_epsM_cell assignment')
    s = sub1(s, '\t\t\tg_epsM_fx = &eps_M_grid_x; g_epsM_fy = &eps_M_grid_y;\n', '', 'face eps assignment')
    s = sub1(s, '\t\t\tg_epsM_fx = nullptr; g_epsM_fy = nullptr;   // the stage\'s face-eps arrays go out of scope below\n', '', 'face eps reset')
    s = sub1(s, '\tconst void* key[3] = {&A_e, g_emb_Qe, g_epsM_cell};', '\tconst void* key[3] = {&A_e, g_emb_Qe, nullptr};', 'CT cache key')
    for v in ('d2lim', 'mm3', 'g_epsM_cell', 'g_epsM_fx', 'g_epsM_fy', 'mm2g', 'S_gx', 'S_gy',
              'const double kappa', 'kappa*dx', 'kappa*dy', 'kappa * (('):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    print('  step 4f conversion limiter: %d chars removed' % (n0 - len(s)))
    return s

def step4g_drop_chi_and_rho_guard(s):
    """Two deletions the measurements licensed.

    (a) chi out of the FLATTENER.  chi stays in the code and in the paper for the
    shock-activated Laplacians (eq:chi, eq:mag_stab, eq:p_stab, four ShockDetector
    call sites); only the fifth one, which fed the HLLEM flattener, goes, so the
    flattener is sigma = s_p, exactly what eq:flattener states.  Measured: on the
    128^2 blast chi is identically zero (kappa_1 c_f = 11.3 against a compression
    of 0.8, because the indicator uses the magnetosonic rather than the sound
    speed that Balsara calibrated 0.4 against); on Orszag-Tang chi exceeds s_p on
    0.317% / 0.255% / 0.189% of cells at N = 128 / 192 / 256, a share that FALLS
    under refinement, and every quality metric is unchanged to four digits; on
    RP5 with the envelope bound it exceeds s_p nowhere at all.

    (b) the density point-value guard.  Measured on Orszag-Tang to t = 0.35, past
    the t = 0.3234 at which the comment recorded a 1024^2 failure:
    min rho_pt/rho_bar = 0.9589 / 0.9599 / 0.9595 at N = 128 / 192 / 256 against
    a threshold of 0.5, with no trend, so the extrapolation to the paper's 1000^2
    is 0.96 and the guard is dead code there.  The recorded failure predates the
    kappa_1 = 0.4 indicator, the split envelope and its bound.  Its action also
    duplicated the positivity guard inside fv_cell_lim_to, which keeps the cell
    average on the same 0.5 test, so two mechanisms were doing one job.
    """
    n0 = len(s)
    # (a) chi out of the flattener
    s = sub1(s, '\t\tShockDetector(hchi, *g_hllem_Qe);   // [fix C] chi of the explicit stage state (ghosts filled by bound_cond_array)\n', '', 'flattener ShockDetector call')
    s = sub1(s, '\t\t\tfor(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++) hchi(i,j) = fmax(hchi(i,j), sp(i,j));\n',
                '\t\t\tfor(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++) hchi(i,j) = sp(i,j);\n',
                'flattener sensor = s_p')
    s = sub1(s, '\tstatic Array hchi(nxCells, nyCells);\n', '\tstatic Array hchi(nxCells, nyCells);   // the flattener sensor s_p of eq:sp\n', 'hchi comment')
    # (b) the density point-value guard
    s = sub1(s, '// The deconvolution rho - (d2x + d2y) rho / 24', '// REMOVED (measured dead at the paper resolution): the deconvolution rho - (d2x + d2y) rho / 24', 'guard comment head') if '// The deconvolution rho - (d2x + d2y) rho / 24' in s else s
    m = re.search(r'[\t ]*// [^\n]*deconvolution[^\n]*\n(?:[\t ]*//[^\n]*\n)*const double RHO_PT_FLOOR = 0\.5;\n', s)
    assert m, 'RHO_PT_FLOOR block not found'
    s = s[:m.start()] + s[m.end():]
    m = re.search(r'[\t ]*// [^\n]*rho_pt_guard[^\n]*\n(?:[\t ]*//[^\n]*\n)*static void rho_pt_guard\(Grid& Qpt, Grid& Qbar\)\{.*?\n\}\n', s, re.S)
    if not m:
        m = re.search(r'static void rho_pt_guard\(Grid& Qpt, Grid& Qbar\)\{.*?\n\}\n', s, re.S)
    assert m, 'rho_pt_guard not found'
    s = s[:m.start()] + s[m.end():]
    for sig in ('static void fv_cell_to(Grid &out, const Grid &in, double sgn, bool rho_guard = false){',
                'template<int N, typename In> static void fvN_to(FieldN<N> &out, const In &in, double sgn, bool rho_guard = false){'):
        assert s.count(sig) == 1, 'signature: ' + sig[:60]
        s = s.replace(sig, sig.replace(', bool rho_guard = false', ''))
    body = ('\t\tif(rho_guard) for(int i = 0; i < NX; i++){ const double thr = RHO_PT_FLOOR*in(i,j,0);\n'
            '\t\t\tif(!(out(i,j,0) >= thr)){ out(i,j,0) = in(i,j,0); ng++; } }   // discard the deconvolution increment\n')
    assert s.count(body) == 2, 'guard bodies: found %d' % s.count(body)
    s = s.replace(body, '')
    for old, new in (('\tlong ng = 0;\n#pragma omp parallel for schedule(static) reduction(+:ng)\n',
                      '#pragma omp parallel for schedule(static)\n'),
                     ('\tg_rho_guard += ng;\n', ''),
                     ('\tif(rho_guard) g_rho_guard += ng;\n', '')):
        if old in s: s = s.replace(old, new)
    for call, fixed in ((', -1.0, true);   // [emb] linear in the momentum; density point-value guard', ', -1.0);   // [emb] linear in the momentum'),
                        ('fv_cell_to(Qp, Q, -1.0, true);', 'fv_cell_to(Qp, Q, -1.0);'),
                        ('fv_cell_to(*g_mfv_R,  u_con, -1.0, true);', 'fv_cell_to(*g_mfv_R,  u_con, -1.0);'),
                        ('fv_cell_to(Qp, u_test, -1.0, true);', 'fv_cell_to(Qp, u_test, -1.0);'),
                        (', -1.0, true);   // with the density point-value guard', ', -1.0);')):
        assert s.count(call) == 1, 'call site: ' + call[:50]
        s = s.replace(call, fixed)
    s = sub1(s, 'extern long g_rho_guard; long g_ent_lim = 0;   // counters: density point-value guard, entropy anti-diffusion limiter',
                'long g_ent_lim = 0;   // counter: faces where the flattener reduced the anti-diffusion', 'counter decl')
    s = sub1(s, 'long g_rho_guard = 0;\n', '', 'g_rho_guard definition')
    s = sub1(s, '", density point-value guard applied " << g_rho_guard << " times, entropy anti-diffusion limited at "',
                '", entropy anti-diffusion limited at "', 'guard report')
    s = sub1(s, '// rho_guard: in the same pass, the density point-value guard of rho_pt_guard (out, in) on every cell.\n', '', 'fv_cell_to doc comment')
    s = sub1(s, '// fv_cell_to for N components; In(i,j,c) is any accessor (a FieldN, or N slots of a Grid); rho_guard as in fv_cell_to (component 0)',
                '// fv_cell_to for N components; In(i,j,c) is any accessor (a FieldN, or N slots of a Grid)', 'fvN_to doc comment')
    for v in ('RHO_PT_FLOOR', 'rho_guard', 'g_rho_guard', 'rho_pt_guard'):
        assert s.count(v) == 0, '%s still referenced %d times' % (v, s.count(v))
    print('  step 4g chi + density guard: %d chars removed' % (n0 - len(s)))
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
    order = ['1', '2', '3', '4a', '4b', '4c', '4d', '4e', '4f', '4g', '5', '6']
    want = order[:order.index(upto) + 1]
    s = open(src).read()
    if '1' in want:  s = step1_dead_code(s)
    if '2' in want:  s = step2_probes(s)
    if '3' in want:  s = step3_flattener(s)
    if '4a' in want: s = step4a_off_switches(s)
    if '4b' in want: s = step4b_probe_machinery(s)
    if '4c' in want: s = step4c_bake_method_params(s)
    if '4d' in want: s = step4d_dead_envelope(s)
    if '4e' in want: s = step4e_senv(s)
    if '4f' in want: s = step4f_drop_conversion_limiter(s)
    if '4g' in want: s = step4g_drop_chi_and_rho_guard(s)
    if '5' in want:  s = step5_comments(s)
    if '6' in want:  s = step6_junk_comments(s)
    open(dst, 'w').write(s)
    print('wrote', dst)
