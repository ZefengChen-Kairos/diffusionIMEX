#!/usr/bin/env python3
"""Apply the blast-wave fixes to Final_emb_omp.cpp (paper-2 code).

    python3 apply_patch.py Final_emb_omp_original.cpp Final_emb_omp_patched.cpp

Changes (all other lines are untouched):
 1. Shock-activated pressure Laplacian of paper 1 (arXiv:2608.15837, eq. (45)/(85)),
        V'(p) * dt * div( chi * |u|_max * dx * (1-eps_M)_face * grad p ),
    added to the implicit pressure operators: the matrix-free second-order operator
    (Pressure_Subsystem_AX_ref), its assembled stencil (pst_assemble) and the fourth-order
    operator (Pressure_Subsystem_AX_HF), so that the defect correction converges to the
    fourth-order operator WITH the dissipation.  The original code computed `lap_p` in
    Pressure_Subsystem_AX_ref but never added it to Ax (dead code).  Switch: g_plap (default 1).
 2. Command-line options (no recompilation needed): cfl=, epsMmin=, epsMmax=, epsBmin=, epsBmax=,
    plap=, kappa1=, swm0=, swm1=, swalpha=, epsmax3=, tag=  (tag is appended to the output name).
 3. The switching-function constants (M0, M1, alpha) and an optional 3x3 widening of the cell
    splitting parameters are runtime variables (defaults = paper values; epsmax3 default 0).
 4. Pressure-fluctuation-aware switch (pswitch=1, default 0): the switching function is driven by
    M_eff = max(M, sqrt(max_nbr |dp|/p)) and M_A,eff = max(M_A, sqrt(max_nbr |dm|/m)), m = |B|^2/2,
    over the four face neighbours.  In a smooth low-Mach flow dp/p = O(M^2), so the limit is unchanged;
    a strong pressure jump (dp/p >~ 0.16) is handed to the explicit HLLEM flux.
"""
import sys

def patch(s):
    n0 = len(s)
    def rep(old, new, count=1):
        nonlocal s
        assert old in s, "pattern not found:\n" + old[:200]
        s = s.replace(old, new, count)

    # ---- CFL as a runtime variable ----
    rep('const double CFL = 0.1;', 'double CFL = 0.1;')

    # ---- globals (defined before test_name(), which uses g_tag) ----
    rep('std::string g_outdir = "";   // output directory, set in main to plot/<test>/',
        'std::string g_outdir = "";   // output directory, set in main to plot/<test>/\n'
        'int g_plap = 1;      // [fix] 1 = shock-activated pressure Laplacian of paper 1 (eq. 45/85), (1-eps_M)-weighted, in L2 and L4\n'
        'std::string g_tag = "";')

    # ---- switching function parameters as runtime variables ----
    rep("""inline double switchS(double M) {
	double SW_EPS0  = 0.15;
	double SW_EPS1  = 0.40;
	double SW_ALPHA = 14.0;
""", """double g_sw_m0 = 0.15, g_sw_m1 = 0.40, g_sw_alpha = 14.0;   // switch parameters (paper: 0.15, 0.4, 14)
int g_epsmax3 = 0;   // 1: cell eps := max over the 3x3 neighbourhood before the face max
int g_pswitch = 0;   // 1: switching function driven by max(M, sqrt(|dp|/p)) and max(M_A, sqrt(|dm|/m))
inline double switchS(double M) {
	double SW_EPS0  = g_sw_m0;
	double SW_EPS1  = g_sw_m1;
	double SW_ALPHA = g_sw_alpha;
""")
    rep("""#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			eps_M_grid_x(i,j) = fmax(cM(i,j), cM(i+1,j));""",
        """	if(g_epsmax3){   // optional: widen eps by one cell in every direction (3x3 max), as the shock indicator chi
		static Array cM3(nxCells, nyCells), cB3(nxCells, nyCells);
		bound_cond_epsilon(cM); bound_cond_epsilon(cB);
#pragma omp parallel for schedule(static)
		for(int j = ghost_cells - 1; j <= ghost_cells + nyCells; j++)
			for(int i = ghost_cells - 1; i <= ghost_cells + nxCells; i++){
				double m = 0.0, b = 0.0;
				for(int dj = -1; dj <= 1; dj++) for(int di = -1; di <= 1; di++){ m = fmax(m, cM(i+di,j+dj)); b = fmax(b, cB(i+di,j+dj)); }
				cM3(i,j) = m; cB3(i,j) = b;
			}
#pragma omp parallel for schedule(static)
		for(int j = ghost_cells - 1; j <= ghost_cells + nyCells; j++)
			for(int i = ghost_cells - 1; i <= ghost_cells + nxCells; i++){ cM(i,j) = cM3(i,j); cB(i,j) = cB3(i,j); }
	}
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			eps_M_grid_x(i,j) = fmax(cM(i,j), cM(i+1,j));""")

    rep("""			double con[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++) con[k] = u_con(i,j,k);
			cM(i,j) = f_epsM(con); cB(i,j) = f_epsB(con);
		}""", """			double con[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++) con[k] = u_con(i,j,k);
			cM(i,j) = f_epsM(con); cB(i,j) = f_epsB(con);
			if(g_pswitch){   // [fix B] pressure-fluctuation-aware switch: M_eff = max(M, sqrt(max_nbr |dp|/p)), M_A,eff = max(M_A, sqrt(max_nbr |dm|/m))
				double pri[NUM_VARS]; ini_con2pri(pri, con);
				const double pc = fmax(pri[4], 1e-300), mc = fmax(0.5*(pri[5]*pri[5]+pri[6]*pri[6]+pri[7]*pri[7]), 1e-300);
				double dp = 0.0, dm = 0.0;
				const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
				for(int q = 0; q < 4; q++){
					const int ii = std::max(0, std::min(nxCells + 2*ghost_cells - 1, i + di[q])), jj = std::max(0, std::min(nyCells + 2*ghost_cells - 1, j + dj[q]));
					double cn[NUM_VARS], pn[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) cn[k] = u_con(ii,jj,k); ini_con2pri(pn, cn);
					const double pnb = fmax(pn[4], 1e-300), mnb = fmax(0.5*(pn[5]*pn[5]+pn[6]*pn[6]+pn[7]*pn[7]), 1e-300);
					dp = fmax(dp, fabs(pnb - pc)/fmin(pnb, pc)); dm = fmax(dm, fabs(mnb - mc)/fmin(mnb, mc));
				}
				const double Mp = sqrt(dp), Mm = sqrt(dm);
				const double eM = 1.0 - switchS(fmax(calMach(con), Mp)), eB = 1.0 - switchS(fmax(calAlfvenMach(con), Mm));
				cM(i,j) = fmin(eps_M_max, fmax(eps_M_min, eM)); cB(i,j) = fmin(eps_B_max, fmax(eps_B_min, eB));
			}
		}""")

    # ---- (1a) matrix-free second-order operator: apply the Laplacian ----
    rep("""				lap_p =  nu_x/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * (getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * (pii - getp(i - 1,j))) + nu_y/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * (pii - getp(i,j-1)));
			}
""", """				lap_p =  nu_x/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * (getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * (pii - getp(i - 1,j))) + nu_y/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * (pii - getp(i,j-1)));
			}
			if(g_plap){   // [fix] paper-1 eq.(45)/(85): V'(p) * dt * div(chi |u|_max dx grad p), (1-eps_M) per face, compact (M-matrix)
				const double a_ = 1.0/(gamma_plasma - 1.0);
				const double gxp = 0.5*(shock_eta(i+1,j)+shock_eta(i,j))*(1.0-eps_M_grid_x(i,j)),   gxm = 0.5*(shock_eta(i,j)+shock_eta(i-1,j))*(1.0-eps_M_grid_x(i-1,j));
				const double gyp = 0.5*(shock_eta(i,j+1)+shock_eta(i,j))*(1.0-eps_M_grid_y(i,j)),   gym = 0.5*(shock_eta(i,j)+shock_eta(i,j-1))*(1.0-eps_M_grid_y(i,j-1));
				const double lp = nu_x/(dx*dx)*(gxp*(getp(i+1,j)-pii) - gxm*(pii-getp(i-1,j))) + nu_y/(dy*dy)*(gyp*(getp(i,j+1)-pii) - gym*(pii-getp(i,j-1)));
				Ax(M_P(i,j)) -= a_*h*lp;
			}
""")
    # ---- (1b) assembled stencil: same term ----
    rep("""	const bool emb = g_emb_L2c && g_emb_Qe;
	const double a = 1.0/(gamma_plasma - 1.0), hh = h*h;""",
        """	const bool emb = g_emb_L2c && g_emb_Qe;
	const double a = 1.0/(gamma_plasma - 1.0), hh = h*h;
	static Array pst_eta(nxCells, nyCells); static long pst_eta_stamp = -1;
	if(g_plap && pst_eta_stamp != g_stage_stamp){ ShockDetector(pst_eta, Qe); pst_eta_stamp = g_stage_stamp; }""")
    rep("""			int nb[4]; double sg[4];
			pst_nb(i+1, j, nb[0], sg[0]); pst_nb(i-1, j, nb[1], sg[1]); pst_nb(i, j+1, nb[2], sg[2]); pst_nb(i, j-1, nb[3], sg[3]);""",
        """			if(g_plap){   // [fix] shock-activated pressure Laplacian (same as the matrix-free form)
				const double nu_x = vx_max*dx, nu_y = vy_max*dy;
				const double gxp = 0.5*(pst_eta(i+1,j)+pst_eta(i,j))*(1.0-eMx(i,j)),   gxm = 0.5*(pst_eta(i,j)+pst_eta(i-1,j))*(1.0-eMx(i-1,j));
				const double gyp = 0.5*(pst_eta(i,j+1)+pst_eta(i,j))*(1.0-eMy(i,j)),   gym = 0.5*(pst_eta(i,j)+pst_eta(i,j-1))*(1.0-eMy(i,j-1));
				const double kxp2 = a*h*nu_x/(dx*dx)*gxp, kxm2 = a*h*nu_x/(dx*dx)*gxm, kyp2 = a*h*nu_y/(dy*dy)*gyp, kym2 = a*h*nu_y/(dy*dy)*gym;
				c0 += kxp2 + kxm2 + kyp2 + kym2; cc[0] -= kxp2; cc[1] -= kxm2; cc[2] -= kyp2; cc[3] -= kym2;
			}
			int nb[4]; double sg[4];
			pst_nb(i+1, j, nb[0], sg[0]); pst_nb(i-1, j, nb[1], sg[1]); pst_nb(i, j+1, nb[2], sg[2]); pst_nb(i, j-1, nb[3], sg[3]);""")
    # ---- (1c) fourth-order operator: same compact term, so the defect correction converges to L4 + Laplacian ----
    rep("""void Pressure_Subsystem_AX_HF(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,Array& eps_M_grid_x, Array& eps_M_grid_y,Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	build_Hbar(Qe, u_con);   // [o3f] (no-op if already built for this stage)""",
        """void Pressure_Subsystem_AX_HF(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,Array& eps_M_grid_x, Array& eps_M_grid_y,Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	build_Hbar(Qe, u_con);   // [o3f] (no-op if already built for this stage)
	static Array hf_eta(nxCells, nyCells); static long hf_eta_stamp = -1;
	if(g_plap && hf_eta_stamp != g_stage_stamp){ ShockDetector(hf_eta, Qe); hf_eta_stamp = g_stage_stamp; }""")
    rep("""			double xdivA = (FXf(i,j) - FXf(i-1,j)) / dx;
			double ydivA = (FYf(i,j) - FYf(i,j-1)) / dy;
			Ax(M_P(i,j)) = 1/(gamma_plasma - 1) * x(M_P(i,j)) - h * h * (xdivA + ydivA);
		}
	}
}""", """			double xdivA = (FXf(i,j) - FXf(i-1,j)) / dx;
			double ydivA = (FYf(i,j) - FYf(i,j-1)) / dy;
			Ax(M_P(i,j)) = 1/(gamma_plasma - 1) * x(M_P(i,j)) - h * h * (xdivA + ydivA);
			if(g_plap){   // [fix] shock-activated pressure Laplacian, identical to the inner operator's
				const double a_ = 1.0/(gamma_plasma - 1.0), nu_x = vx_max*dx, nu_y = vy_max*dy, pii = x(M_P(i,j));
				const double gxp = 0.5*(hf_eta(i+1,j)+hf_eta(i,j))*(1.0-eps_M_grid_x(i,j)),   gxm = 0.5*(hf_eta(i,j)+hf_eta(i-1,j))*(1.0-eps_M_grid_x(i-1,j));
				const double gyp = 0.5*(hf_eta(i,j+1)+hf_eta(i,j))*(1.0-eps_M_grid_y(i,j)),   gym = 0.5*(hf_eta(i,j)+hf_eta(i,j-1))*(1.0-eps_M_grid_y(i,j-1));
				const double lp = nu_x/(dx*dx)*(gxp*(x(M_P(i+1,j))-pii) - gxm*(pii-x(M_P(i-1,j)))) + nu_y/(dy*dy)*(gyp*(x(M_P(i,j+1))-pii) - gym*(pii-x(M_P(i,j-1))));
				Ax(M_P(i,j)) -= a_*h*lp;
			}
		}
	}
}""")

    # ---- (2) command-line options ----
    rep('int main(int argc, char *argv[]){\n\tswitch (test) {',
        'int main(int argc, char *argv[]){\n'
        '\tfor(int a = 1; a < argc; a++){ std::string k(argv[a]); auto eq = k.find(\'=\'); if(eq == std::string::npos) continue; std::string key = k.substr(0,eq); double v = atof(k.c_str()+eq+1);\n'
        '\t\tif(key=="cfl") CFL = v; else if(key=="epsMmin") eps_M_min = v; else if(key=="epsMmax") eps_M_max = v; else if(key=="epsBmin") eps_B_min = v; else if(key=="epsBmax") eps_B_max = v;\n'
        '\t\telse if(key=="plap") g_plap = (int)v; else if(key=="tag") g_tag = k.substr(eq+1); else if(key=="kappa1") shock_kappa1 = v;\n'
        '\t\telse if(key=="swm0") g_sw_m0 = v; else if(key=="swm1") g_sw_m1 = v; else if(key=="swalpha") g_sw_alpha = v; else if(key=="epsmax3") g_epsmax3 = (int)v; else if(key=="pswitch") g_pswitch = (int)v;\n'
        '\t\telse { std::cout << "unknown arg " << k << std::endl; return 1; } }\n'
        '\tstd::cout << "ARGS cfl=" << CFL << " epsM=[" << eps_M_min << "," << eps_M_max << "] epsB=[" << eps_B_min << "," << eps_B_max << "] plap=" << g_plap << " kappa1=" << shock_kappa1 << " switch(m0,m1,alpha)=" << g_sw_m0 << "," << g_sw_m1 << "," << g_sw_alpha << " epsmax3=" << g_epsmax3 << " pswitch=" << g_pswitch << std::endl;\n'
        '\tswitch (test) {')
    rep('\tos << "_CFL" << CFL;\n\ttestname = os.str();', '\tos << "_CFL" << CFL; if(!g_tag.empty()) os << "_" << g_tag;\n\ttestname = os.str();')
    return s

if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    s = open(src).read()
    open(dst, 'w').write(patch(s))
    print('wrote', dst)
