// =====================================================================================================================
//  SPB6_TV_IMEX.cpp
//  SPB 6+1 方程两相模型的 TV 分裂 IMEX（docs/derivation/B_TV_IMEX.md），代码结构仿照 2DEuler_TV.cpp
//  （Grid、Qe/Qs/Hs 的半隐式 IMEX-RK、TVDReconstruction、caldt、bound_cond_u）。
//    显式 = 体积输运形式的 MUSCL-Rusanov（与线性路径、面速度 ubar 的 path-conservative Rusanov 等价），
//           耗散速度取物质速度、且不小于面态速度（B3）；BB (36)-(38) 都包含在对 alpha_k, m_k, eps_k, rho u, E 的 Rusanov 耗散里；
//    面   = 面速度/面压力（B4.1）：半格动量平衡 + 声学 Riemann 解的迎风耗散，权重 (1 - c tau/dx)^+；
//           增量通量对总面速度 u_f 迎风（B4.2）；
//    隐式 = 每格一个未知量（混合压力 P）的非线性方程 G(P) = pbar_cl(P) - P = 0，pbar_cl 由统一的局部闭合给出：
//           离散 Gibbs 关系 + 后向 Euler 的松弛（速率 mu, theta, nu 为参数）+ 混合能量一致性（SPB 的 "+1" 在闭合内）；
//    时间 = 半隐式 IMEX-RK（2DEuler_TV 的 Ae/Ai/Hs 结构），RKs = 1, 2 (LSDIRK2), 3 (SSP2(3,3,2))；
//           瞬时热/化学松弛的激活每步判断一次并投影到平衡流形（B6.7）。
//
//  编译: g++ -O2 -std=c++17 SPB6_TV_IMEX.cpp -o spb6_imex
//  运行: ./spb6_imex key=value ...   共用键见 spb6_common.h；另有 rk, cflac, ramp, tol, dt
// =====================================================================================================================
#include "spb6_common.h"

// ------------------------------------------------------------------ IMEX 参数
int RKs = 2;                // IMEX-RK 级数: 1 (后向/前向 Euler), 2 (LSDIRK2(2,2,2)), 3 (SSP2(3,3,2))
double CFL_ac = 0.4;        // 起步用的声学 CFL
int RAMP = 1;               // 从静止/间断起算：首步声学步长，之后每步至多放大 1.2 倍；被 1.2 倍限制的步用一阶（后向 Euler）
double NEWTON_TOL = 1e-13;  // 相对 |P| + rho c_W^2 + P_E（B8）
int NEWTON_MAX = 60;
double DTFIX = 0.0;         // > 0：固定时间步（时间收敛阶试验用；不做起步斜坡）
long st_newton_it = 0, st_newton_calls = 0, st_stage = 0;
double st_maxres = 0;

// B6.7：激活与投影（每步一次，P-S 的算子顺序）。瞬时松弛（速率 = inf）在激活时是状态的跳跃；若在级内激活，
// 第一级把状态投影到平衡流形上（与 dt 无关的 O(1) 变化，K_1 = O(1/dt)），IMEX-RK 的级组合（LSDIRK2 的 Q^s_2 = U^n + ((1-g)/g)(U_1 - U^n)，
// 系数 2.41）会把这个跳跃外推，微量相的含量变负。刚性精确 IMEX-RK 的 AP 性质要求初值在流形上（well-prepared）。
// 所以：每步开始时按 U^n 判断激活（等价于 P-S 在双曲步之后判断），把 U^n 投影到该闭合的流形上（s = 0 的闭合，守恒 m 或 rho、rho u、E），
// 本步各级冻结这些类型。有限速率的类型不投影（速率有限，状态以有限速度松弛）。
std::vector<int> STEP_TYPE;   // 本步各内部格的闭合类型（空 = 各级自行判断）
static int cell_type_U(const Grid& U, int i, int j) {
	const double rho = U(i, j, IM1) + U(i, j, IM2);
	const double rhoe = U(i, j, IEN) - 0.5 * (U(i, j, IMU) * U(i, j, IMU) + U(i, j, IMV) * U(i, j, IMV)) / rho;
	return closure_type_of(U(i, j, IA1), U(i, j, IM1), U(i, j, IM2), U(i, j, IE1), U(i, j, IE2), rhoe);
}
static bool project_cell(Grid& U, int i, int j, int type) {
	if (!(type == CL_PT || type == CL_PTG)) return true;
	const double rho = U(i, j, IM1) + U(i, j, IM2);
	ClIn in{}; in.atr[0] = U(i, j, IA1); in.atr[1] = 1.0 - U(i, j, IA1); in.mtr[0] = U(i, j, IM1); in.mtr[1] = U(i, j, IM2);
	in.etr[0] = U(i, j, IE1); in.etr[1] = U(i, j, IE2); in.s = 0.0;
	in.rhoeC = U(i, j, IEN) - 0.5 * (U(i, j, IMU) * U(i, j, IMU) + U(i, j, IMV) * U(i, j, IMV)) / rho;
	in.tau = 0.0; in.type = type; in.Z[0] = in.Z[1] = 1.0; in.nodelta = false;
	in.pguess = std::max(p_mix(in.atr[0], in.mtr[0], in.mtr[1], in.rhoeC), 1.0);
	ClOut o; if (!closure(in, o)) return false;
	U(i, j, IA1) = o.a[0]; U(i, j, IM1) = o.m[0]; U(i, j, IM2) = o.m[1]; U(i, j, IE1) = o.e[0]; U(i, j, IE2) = o.e[1];
	return true;
}
// 对整个场：判断类型、投影；投影失败的格退回力学类型（不投影）
static void activate_and_project(Grid& U) {
	const int nx = U.nx, ny = U.ny;
	STEP_TYPE.assign((size_t)nx * ny, CL_P);
	for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
		int t = cell_type_U(U, IG(i), IG(j));
		if (!project_cell(U, IG(i), IG(j), t)) t = std::isinf(mu_rate) ? CL_P : CL_MU;
		STEP_TYPE[(size_t)i + (size_t)j * nx] = t;
	}
}

// ==================================================================================================== TV 分裂 IMEX
// ---- B3 显式算子：共用限制因子的 minmod（在守恒量上）；MUSCL-Rusanov（中心部分与耗散都用重构面值），耗散速度取物质速度
static void TVDReconstruction(const Grid& Q, int dir, std::vector<double>& sig, std::vector<double>* phis = nullptr) {
	const int g = ghost_cells, NXT = Q.nx + 2 * g, NYT = Q.ny + 2 * g;
	sig.assign((size_t)NXT * NYT * NUM_VARS, 0.0);
	if (phis) phis->assign((size_t)NXT * NYT, 0.0);
	if (SpaceOrder < 2) return;
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		int ia = i, ja = j, ib = i, jb = j;
		if (dir == 0) { ia = i - 1; ib = i + 1; if (ia < 0 || ib >= NXT) continue; } else { ja = j - 1; jb = j + 1; if (ja < 0 || jb >= NYT) continue; }
		double phi = 1.0;
		for (int k = 0; k < NUM_VARS; k++) {
			const double xa = Q(ia, ja, k), xc = Q(i, j, k), xb = Q(ib, jb, k);
			const double a = xc - xa, b = xb - xc, tiny = 1e-12 * (std::fabs(xa) + std::fabs(xc) + std::fabs(xb)) + 1e-300;
			if (std::fabs(a) <= tiny && std::fabs(b) <= tiny) continue;          // 该量在模板内（到舍入）为常数：不约束
			if (a * b <= 0) { phi = 0.0; break; }                                // 极值点：minmod 斜率为 0
			const double c = 0.5 * (a + b), m = std::fabs(a) < std::fabs(b) ? a : b;
			phi = std::min(phi, m / c);
		}
		for (int k = 0; k < NUM_VARS; k++) sig[((size_t)i + (size_t)j * NXT) * NUM_VARS + k] = phi * 0.5 * (Q(ib, jb, k) - Q(ia, ja, k));
		if (phis) (*phis)[(size_t)i + (size_t)j * NXT] = phi;
	}
}

struct StageData {
	int nx, ny, NXT, NYT; double rx, ry, tau;
	Grid Xs;                       // 显式更新后的 X*
	Array a2s;                     // alpha2*（不等于 1 - alpha1*：两者之和 = 1 - r delta[ubar^e]）
	std::vector<double> cont;      // 每格单位体积内容 (alpha1, alpha2, m1, m2, e1, e2) / sum alpha*
	Array ust, vst;                // u*, v*（供体速度）
	Array ubsx, kpx, ubex, ubsy, kpy, ubey;   // 面量（x 面 (i,j) = i+1/2）：u_f = ubs - kp (P_R - P_L)，ube = 显式中心通量的面速度
	Array wLx, wLy, pvx, pvy;      // 面压力：p_f = wL P_L + (1 - wL) P_R - pv（B4.1 的半格/声学 Riemann 解）
	std::vector<int> type; std::vector<double> Z0, Z1;
	std::vector<double> Kw;        // 残量尺度 K + P_E（X* 处，B8）
	StageData(int x, int y) : nx(x), ny(y), NXT(x + 2 * ghost_cells), NYT(y + 2 * ghost_cells), Xs(x, y), a2s(x, y), ust(x, y), vst(x, y),
		ubsx(x, y), kpx(x, y), ubex(x, y), ubsy(x, y), kpy(x, y), ubey(x, y), wLx(x, y, 0.5), wLy(x, y, 0.5), pvx(x, y), pvy(x, y) {
		cont.assign((size_t)NXT * NYT * 6, 0.0); Kw.assign((size_t)NXT * NYT, 0.0); type.assign((size_t)NXT * NYT, CL_P); Z0.assign((size_t)NXT * NYT, 1.0); Z1.assign((size_t)NXT * NYT, 1.0);
	}
	inline size_t id(int i, int j) const { return (size_t)i + (size_t)j * NXT; }
};

// 面上的隐式增量（B4.2'）：总体积通量对 u_f 迎风。显式部分已按 ubar^e 从其上游格扫过 ubar^e tau，增量补足到 u_f tau：
//   I_f[X] = Xhat_L wL + Xhat_R wR,  wL = (u_f)^+ - (ubar^e)^+,  wR = (u_f)^- - (ubar^e)^-
// 同号时就是 Xhat_up du（du < 0 退回的是显式扫过来的同一种物质）；反向时退回显式扫过的全部、再从另一侧扫 u_f。
// wL + wR = du，故 sum_k alphahat 的增量 = du（B4.4 不变）；I_f 对 P 连续（分段线性）。
struct FaceVals { double uf, du, pf, wL, wR, utL, utR; const double* cL; const double* cR; };
inline void face_weights(FaceVals& f, double ube) { f.wL = std::max(f.uf, 0.0) - std::max(ube, 0.0); f.wR = std::min(f.uf, 0.0) - std::min(ube, 0.0); }
inline bool is_wall_x(const StageData& S, int i) { return x_bc == x_Reflective && (i == ghost_cells - 1 || i == S.nx + ghost_cells - 1); }
inline bool is_wall_y(const StageData& S, int j) { return y_bc == y_Reflective && (j == ghost_cells - 1 || j == S.ny + ghost_cells - 1); }
// x 面 (i+1/2, j)，左格压力 PL，右格压力 PR（B4.1, B4.2）
inline void face_x(const StageData& S, int i, int j, double PL, double PR, FaceVals& f) {
	f.pf = S.wLx(i, j) * PL + (1.0 - S.wLx(i, j)) * PR - S.pvx(i, j);
	f.cL = &S.cont[S.id(i, j) * 6]; f.cR = &S.cont[S.id(i + 1, j) * 6]; f.utL = S.vst(i, j); f.utR = S.vst(i + 1, j);
	if (is_wall_x(S, i)) { f.uf = 0.0; f.du = 0.0; f.wL = f.wR = 0.0; return; }
	f.uf = S.ubsx(i, j) - S.kpx(i, j) * (PR - PL);
	f.du = f.uf - S.ubex(i, j);
	face_weights(f, S.ubex(i, j));
}
inline void face_y(const StageData& S, int i, int j, double PL, double PR, FaceVals& f) {
	f.pf = S.wLy(i, j) * PL + (1.0 - S.wLy(i, j)) * PR - S.pvy(i, j);
	f.cL = &S.cont[S.id(i, j) * 6]; f.cR = &S.cont[S.id(i, j + 1) * 6]; f.utL = S.ust(i, j); f.utR = S.ust(i, j + 1);
	if (is_wall_y(S, j)) { f.uf = 0.0; f.du = 0.0; f.wL = f.wR = 0.0; return; }
	f.uf = S.ubsy(i, j) - S.kpy(i, j) * (PR - PL);
	f.du = f.uf - S.ubey(i, j);
	face_weights(f, S.ubey(i, j));
}
inline double inc(const FaceVals& f, int k) { return f.cL[k] * f.wL + f.cR[k] * f.wR; }
// 增量质量 M = ML + MR；法向动量通量 u_f M；切向 utL ML + utR MR；总能 sum eps + 1/2 u_f^2 M + 1/2 (utL^2 ML + utR^2 MR) + p_f u_f
struct FaceFlux { double M, Mt, E; };
inline FaceFlux face_flux(const FaceVals& f) {
	const double ML = (f.cL[2] + f.cL[3]) * f.wL, MR = (f.cR[2] + f.cR[3]) * f.wR, M = ML + MR;
	FaceFlux F; F.M = M; F.Mt = f.utL * ML + f.utR * MR;
	F.E = inc(f, 4) + inc(f, 5) + 0.5 * f.uf * f.uf * M + 0.5 * (f.utL * f.utL * ML + f.utR * f.utR * MR) + f.pf * f.uf;
	return F;
}

struct CellOut { ClOut cl; double mx, my, E; int type; };
// 单格残量 G = pbar_cl - Pc（B4.3-B4.6, B5.4）。邻格压力由调用者给（边界虚格按边界条件取值）
static bool cell_residual(const StageData& S, int i, int j, double Pc, double PW, double PE, double PS, double PN, double& G, CellOut* out) {
	FaceVals fe, fw, fn, fs;
	face_x(S, i, j, Pc, PE, fe); face_x(S, i - 1, j, PW, Pc, fw);
	const bool twoD = S.ny > 1;
	if (twoD) { face_y(S, i, j, Pc, PN, fn); face_y(S, i, j - 1, PS, Pc, fs); }
	double X[6] = {S.Xs(i, j, IA1), S.a2s(i, j), S.Xs(i, j, IM1), S.Xs(i, j, IM2), S.Xs(i, j, IE1), S.Xs(i, j, IE2)};
	for (int k = 0; k < 6; k++) {
		X[k] -= S.rx * (inc(fe, k) - inc(fw, k));
		if (twoD) X[k] -= S.ry * (inc(fn, k) - inc(fs, k));
	}
	const double s = 1.0 - X[0] - X[1];
	const FaceFlux FEe = face_flux(fe), FEw = face_flux(fw);
	double mx = S.Xs(i, j, IMU) - S.rx * (fe.pf - fw.pf) - S.rx * (fe.uf * FEe.M - fw.uf * FEw.M);
	double my = S.Xs(i, j, IMV) - S.rx * (FEe.Mt - FEw.Mt);
	double E = S.Xs(i, j, IEN) - S.rx * (FEe.E - FEw.E);
	if (twoD) {
		const FaceFlux FEn = face_flux(fn), FEs = face_flux(fs);
		my -= S.ry * (fn.pf - fs.pf) + S.ry * (fn.uf * FEn.M - fs.uf * FEs.M);
		mx -= S.ry * (FEn.Mt - FEs.Mt);
		E -= S.ry * (FEn.E - FEs.E);
	}
	const double rho = X[2] + X[3];
	if (!(rho > 0)) return false;
	ClIn in{}; in.atr[0] = X[0]; in.atr[1] = X[1]; in.mtr[0] = X[2]; in.mtr[1] = X[3]; in.etr[0] = X[4]; in.etr[1] = X[5];
	in.s = s; in.rhoeC = E - 0.5 * (mx * mx + my * my) / rho; in.tau = S.tau; in.type = S.type[S.id(i, j)];
	in.Z[0] = S.Z0[S.id(i, j)]; in.Z[1] = S.Z1[S.id(i, j)]; in.pguess = Pc; in.nodelta = false;
	ClOut o; if (!closure(in, o)) {
		if (VERBOSE >= 3) fprintf(stderr, "    closure fail cell (%d,%d) type %d atr %.6e %.6e mtr %.6e %.6e etr %.6e %.6e s %.3e rhoeC %.6e Pc %.6e\n", i - ghost_cells, j - ghost_cells, in.type, in.atr[0], in.atr[1], in.mtr[0], in.mtr[1], in.etr[0], in.etr[1], s, in.rhoeC, Pc);
		return false; }
	G = o.pbar - Pc;
	if (out) { out->cl = o; out->mx = mx; out->my = my; out->E = E; out->type = in.type; }
	return std::isfinite(G);
}

// 邻格压力（边界：透射/固壁 = 本格，周期 = 绕回）
struct PField { int nx, ny; std::vector<double> P; inline double& at(int i, int j) { return P[(size_t)i + (size_t)j * nx]; } inline double at(int i, int j) const { return P[(size_t)i + (size_t)j * nx]; } };
inline double Pnb(const PField& P, int i, int j, int di, int dj, double Pc) {
	int ii = i + di, jj = j + dj;
	if (ii < 0 || ii >= P.nx) { if (x_bc == x_Periodic) ii = (ii + P.nx) % P.nx; else return Pc; }
	if (jj < 0 || jj >= P.ny) { if (y_bc == y_Periodic && P.ny > 1) jj = (jj + P.ny) % P.ny; else return Pc; }
	return P.at(ii, jj);
}
static bool residual_at(const StageData& S, const PField& P, int i, int j, double dc, int which, double dn, double& G, CellOut* out = nullptr) {
	// which: -1 不扰动邻格；0 W, 1 E, 2 S, 3 N 扰动 dn
	const double Pc = P.at(i, j) + dc;
	double PW = Pnb(P, i, j, -1, 0, Pc), PE = Pnb(P, i, j, 1, 0, Pc), PS = Pnb(P, i, j, 0, -1, Pc), PN = Pnb(P, i, j, 0, 1, Pc);
	const bool perW = (i > 0 || x_bc == x_Periodic), perE = (i < P.nx - 1 || x_bc == x_Periodic), perS = (j > 0 || (y_bc == y_Periodic && P.ny > 1)), perN = (j < P.ny - 1 || (y_bc == y_Periodic && P.ny > 1));
	if (which == 0 && perW) PW += dn;
	if (which == 1 && perE) PE += dn;
	if (which == 2 && perS) PS += dn;
	if (which == 3 && perN) PN += dn;
	return cell_residual(S, IG(i), IG(j), Pc, PW, PE, PS, PN, G, out);
}

// 稀疏矩阵（每行至多 5 个非零元）+ BiCGStab（Jacobi 预条件）；一维非周期用 Thomas
struct SpRow { int n = 0; int col[5]; double val[5]; };
static bool bicgstab(const std::vector<SpRow>& A, const std::vector<double>& b, std::vector<double>& x, double rtol, int maxit) {
	const int n = (int)b.size(); std::vector<double> dinv(n, 1.0);
	for (int r = 0; r < n; r++) for (int q = 0; q < A[r].n; q++) if (A[r].col[q] == r && A[r].val[q] != 0) dinv[r] = 1.0 / A[r].val[q];
	auto mv = [&](const std::vector<double>& v, std::vector<double>& y) { for (int r = 0; r < n; r++) { double s = 0; for (int q = 0; q < A[r].n; q++) s += A[r].val[q] * v[A[r].col[q]]; y[r] = s; } };
	auto dot = [&](const std::vector<double>& a, const std::vector<double>& c) { double s = 0; for (int r = 0; r < n; r++) s += a[r] * c[r]; return s; };
	x.assign(n, 0.0); for (int r = 0; r < n; r++) x[r] = b[r] * dinv[r];
	std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n), Ax(n);
	mv(x, Ax); for (int q = 0; q < n; q++) r[q] = b[q] - Ax[q]; r0 = r;
	const double bn = std::sqrt(dot(b, b)) + 1e-300; double rho = 1, alpha = 1, omega = 1;
	for (int it = 0; it < maxit; it++) {
		if (std::sqrt(dot(r, r)) <= rtol * bn) return true;
		const double rho1 = dot(r0, r); if (rho1 == 0) return false;
		const double beta = (rho1 / rho) * (alpha / omega); rho = rho1;
		for (int q = 0; q < n; q++) p[q] = r[q] + beta * (p[q] - omega * v[q]);
		for (int q = 0; q < n; q++) ph[q] = dinv[q] * p[q];
		mv(ph, v); alpha = rho / dot(r0, v);
		for (int q = 0; q < n; q++) s[q] = r[q] - alpha * v[q];
		if (std::sqrt(dot(s, s)) <= rtol * bn) { for (int q = 0; q < n; q++) x[q] += alpha * ph[q]; return true; }
		for (int q = 0; q < n; q++) sh[q] = dinv[q] * s[q];
		mv(sh, t); omega = dot(t, s) / dot(t, t);
		for (int q = 0; q < n; q++) { x[q] += alpha * ph[q] + omega * sh[q]; r[q] = s[q] - omega * t[q]; }
	}
	return std::sqrt(dot(r, r)) <= 1e3 * rtol * bn;
}

// B8：全局 Newton，G_i(P) = pbar_cl,i - P_i
static bool newton_solve(const StageData& S, PField& P, int& iters, double& resfinal) {
	const int nx = P.nx, ny = P.ny, n = nx * ny;
	std::vector<double> G(n), Gn(n), dP(n); PField Pn = P;
	auto evalAll = [&](const PField& Q, std::vector<double>& R) -> bool {
		for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) { double g; if (!residual_at(S, Q, i, j, 0, -1, 0, g)) return false; R[i + j * nx] = g; }
		return true; };
	// 残量范数：|G| / (|P| + K + P_E)（见 imex_stage 中 Kw 的注释）；按 |P| 归一在 p << pi 时达不到
	auto norm = [&](const std::vector<double>& R, const PField& Q) { double m = 0; for (int q = 0; q < n; q++) { const int i = q % nx, j = q / nx; m = std::max(m, std::fabs(R[q]) / (std::fabs(Q.P[q]) + S.Kw[S.id(IG(i), IG(j))] + 1.0)); } return m; };
	if (!evalAll(P, G)) { iters = 0; resfinal = INF; return false; }
	double res = norm(G, P);
	const bool tri = (ny == 1 && x_bc != x_Periodic);
	std::vector<SpRow> A(n); std::vector<double> lo(n), di(n), up(n), rhs(n);
	for (int it = 0; it < NEWTON_MAX; it++) {
		iters = it; resfinal = res;
		if (VERBOSE >= 3) { int qm = 0; for (int q = 0; q < n; q++) if (std::fabs(G[q]) > std::fabs(G[qm])) qm = q; fprintf(stderr, "    newton it %d res %.3e  worst cell %d G %.4e P %.8e\n", it, res, qm, G[qm], P.P[qm]); }
		if (res < NEWTON_TOL) return true;
		// Jacobian（逐格差分）
		for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
			const int r = i + j * nx; const double Pc = P.at(i, j), h = 1e-7 * (std::fabs(Pc) + 1.0);
			double g1; if (!residual_at(S, P, i, j, h, -1, 0, g1)) return false;
			SpRow& row = A[r]; row.n = 0; row.col[row.n] = r; row.val[row.n++] = (g1 - G[r]) / h;
			const int dd[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
			for (int w = 0; w < 4; w++) {
				if (ny == 1 && w >= 2) break;
				int ii = i + dd[w][0], jj = j + dd[w][1];
				if (ii < 0 || ii >= nx) { if (x_bc == x_Periodic) ii = (ii + nx) % nx; else continue; }
				if (jj < 0 || jj >= ny) { if (y_bc == y_Periodic && ny > 1) jj = (jj + ny) % ny; else continue; }
				const double hn = 1e-7 * (std::fabs(P.at(ii, jj)) + 1.0);
				double g2; if (!residual_at(S, P, i, j, 0, w, hn, g2)) return false;
				const int c = ii + jj * nx; bool merged = false;
				for (int q = 0; q < row.n; q++) if (row.col[q] == c) { row.val[q] += (g2 - G[r]) / hn; merged = true; }
				if (!merged) { row.col[row.n] = c; row.val[row.n++] = (g2 - G[r]) / hn; }
			}
		}
		for (int q = 0; q < n; q++) rhs[q] = -G[q];
		if (tri) {
			for (int r = 0; r < n; r++) { lo[r] = up[r] = 0; for (int q = 0; q < A[r].n; q++) { const int c = A[r].col[q]; if (c == r) di[r] = A[r].val[q]; else if (c == r - 1) lo[r] = A[r].val[q]; else if (c == r + 1) up[r] = A[r].val[q]; } }
			std::vector<double> cp(n), dp(n);
			cp[0] = up[0] / di[0]; dp[0] = rhs[0] / di[0];
			for (int r = 1; r < n; r++) { const double m = di[r] - lo[r] * cp[r - 1]; cp[r] = up[r] / m; dp[r] = (rhs[r] - lo[r] * dp[r - 1]) / m; }
			dP[n - 1] = dp[n - 1]; for (int r = n - 2; r >= 0; r--) dP[r] = dp[r] - cp[r] * dP[r + 1];
		} else {
			if (!bicgstab(A, rhs, dP, 1e-13, 4000)) return false;
		}
		// 线搜索：可容许且残量下降
		double lam = 1.0; bool acc = false;
		for (int ls = 0; ls < 40; ls++) {
			for (int q = 0; q < n; q++) Pn.P[q] = P.P[q] + lam * dP[q];
			if (evalAll(Pn, Gn)) { const double rn = norm(Gn, Pn); if (rn < res * (1.0 - 1e-4 * lam) || rn < NEWTON_TOL) { acc = true; res = rn; break; } }
			lam *= 0.5;
		}
		if (!acc) { resfinal = res; return false; }
		std::swap(P.P, Pn.P); std::swap(G, Gn);
		st_newton_it++;
	}
	iters = NEWTON_MAX; resfinal = res; return res < NEWTON_TOL;
}

// 一级 IMEX（B10）：U = Phi(Qe, Qs; tau)；P 进入时为初猜，返回时为解
static bool imex_stage(const Grid& Qe_in, const Grid& Qs_in, double tau, Grid& U, PField& P) {
	const int g = ghost_cells, nx = U.nx, ny = U.ny;
	Grid Qe = Qe_in, Qs = Qs_in; bound_cond_u(Qe); bound_cond_u(Qs);
	StageData S(nx, ny); S.tau = tau; S.rx = tau / dx; S.ry = tau / dy;
	const int NXT = S.NXT, NYT = S.NYT;
	// 显式级态的速度（只需 rho > 0）
	Array ue(nx, ny), ve(nx, ny);
	double ax = 0, ay = 0;
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		const double rho = Qe(i, j, IM1) + Qe(i, j, IM2); if (!(rho > 0)) return false;
		ue(i, j) = Qe(i, j, IMU) / rho; ve(i, j) = Qe(i, j, IMV) / rho;
		if (i >= g && i < nx + g && j >= g && j < ny + g) { ax = std::max(ax, std::fabs(ue(i, j))); ay = std::max(ay, std::fabs(ve(i, j))); }
	}
	// 显式 Rusanov 通量 (B3.1)
	std::vector<double> sx, sy, phx, phy; TVDReconstruction(Qe, 0, sx, &phx); if (ny > 1) TVDReconstruction(Qe, 1, sy, &phy);
	// MUSCL-Rusanov：中心部分与耗散都用重构面值（对线性输运就是重构值的迎风通量，TVD，微量相保正）。
	// BB 的“中心部分用格均值”形式对线性输运给出 u X_i + 1/4 u (sig_i + sig_{i+1})，含下游斜率，不是 TVD（B3）。
	auto face_states = [&](int i0, int j0, int i1, int j1, const std::vector<double>& sig, int dir, double* XL, double* XR, double& uL, double& uR) {
		const size_t k0 = ((size_t)i0 + (size_t)j0 * NXT) * NUM_VARS, k1 = ((size_t)i1 + (size_t)j1 * NXT) * NUM_VARS;
		for (int k = 0; k < NUM_VARS; k++) { XL[k] = Qe(i0, j0, k) + 0.5 * sig[k0 + k]; XR[k] = Qe(i1, j1, k) - 0.5 * sig[k1 + k]; }
		const int kn = dir == 0 ? IMU : IMV;
		uL = XL[kn] / (XL[IM1] + XL[IM2]); uR = XR[kn] / (XR[IM1] + XR[IM2]);
	};
	// 耗散速度 a_Omega 必须不小于中心部分所用面态的速度（B3 的单调性条件）：把格更新写成面值的组合，下游面值 X^R_{i+1/2} 的系数是
	// r (a - u^R_{i+1/2}) / 2，a < u^R 时下游格的含量以负系数进入，微量相含量变负。MUSCL 的面速度 (rho u)^{L,R} / rho^{L,R} 不是格速度的
	// 凸组合，可能超出 max|u_cell|（水-气界面：超出 0.34 m/s 就使界面格的微量气体质量变负）。
	for (int j = g; j < ny + g; j++) for (int i = g - 1; i < nx + g; i++) { double XL[NUM_VARS], XR[NUM_VARS], uL, uR; face_states(i, j, i + 1, j, sx, 0, XL, XR, uL, uR); ax = std::max(ax, std::max(std::fabs(uL), std::fabs(uR))); }
	if (ny > 1) for (int j = g - 1; j < ny + g; j++) for (int i = g; i < nx + g; i++) { double XL[NUM_VARS], XR[NUM_VARS], uL, uR; face_states(i, j, i, j + 1, sy, 1, XL, XR, uL, uR); ay = std::max(ay, std::max(std::fabs(uL), std::fabs(uR))); }
	auto flux_face = [&](int i0, int j0, int i1, int j1, const std::vector<double>& sig, int dir, double a, double* F, double& ubar) {
		double XL[NUM_VARS], XR[NUM_VARS], uL, uR; face_states(i0, j0, i1, j1, sig, dir, XL, XR, uL, uR);
		for (int k = 0; k < NUM_VARS; k++) F[k] = 0.5 * (XL[k] * uL + XR[k] * uR) - 0.5 * a * (XR[k] - XL[k]);
		ubar = 0.5 * (uL + uR);   // = sum_k F[alpha_k]（alpha2 面值 = 1 - alpha1 面值）
	};
	std::vector<double> FXv((size_t)NXT * NYT * NUM_VARS, 0.0), FYv; 
	for (int j = g; j < ny + g; j++) for (int i = g - 1; i < nx + g; i++) { double ub; flux_face(i, j, i + 1, j, sx, 0, ax, &FXv[((size_t)i + (size_t)j * NXT) * NUM_VARS], ub); S.ubex(i, j) = ub; }
	if (ny > 1) { FYv.assign((size_t)NXT * NYT * NUM_VARS, 0.0); for (int j = g - 1; j < ny + g; j++) for (int i = g; i < nx + g; i++) { double ub; flux_face(i, j, i, j + 1, sy, 1, ay, &FYv[((size_t)i + (size_t)j * NXT) * NUM_VARS], ub); S.ubey(i, j) = ub; } }
	// X* = Qs - r delta[F]；alpha2* 用 alpha2 的通量 ubar^e - F[alpha1]
	S.Xs = Qs;
	for (int j = g; j < ny + g; j++) for (int i = g; i < nx + g; i++) {
		const double* fp = &FXv[((size_t)i + (size_t)j * NXT) * NUM_VARS]; const double* fm = &FXv[((size_t)(i - 1) + (size_t)j * NXT) * NUM_VARS];
		for (int k = 0; k < NUM_VARS; k++) S.Xs(i, j, k) = Qs(i, j, k) - S.rx * (fp[k] - fm[k]);
		double a2 = (1.0 - Qs(i, j, IA1)) - S.rx * ((S.ubex(i, j) - fp[IA1]) - (S.ubex(i - 1, j) - fm[IA1]));
		if (ny > 1) {
			const double* gp = &FYv[((size_t)i + (size_t)j * NXT) * NUM_VARS]; const double* gm = &FYv[((size_t)i + (size_t)(j - 1) * NXT) * NUM_VARS];
			for (int k = 0; k < NUM_VARS; k++) S.Xs(i, j, k) -= S.ry * (gp[k] - gm[k]);
			a2 -= S.ry * ((S.ubey(i, j) - gp[IA1]) - (S.ubey(i, j - 1) - gm[IA1]));
		}
		S.a2s(i, j) = a2;
	}
	bound_cond_u(S.Xs); bound_cond_array(S.a2s);
	// 单位体积内容、供体速度、面量、闭合类型
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		const double sa = S.Xs(i, j, IA1) + S.a2s(i, j), rho = S.Xs(i, j, IM1) + S.Xs(i, j, IM2);
		if (!(sa > 0 && rho > 0)) return false;
		double* c = &S.cont[S.id(i, j) * 6];
		c[0] = S.Xs(i, j, IA1) / sa; c[1] = S.a2s(i, j) / sa; c[2] = S.Xs(i, j, IM1) / sa; c[3] = S.Xs(i, j, IM2) / sa; c[4] = S.Xs(i, j, IE1) / sa; c[5] = S.Xs(i, j, IE2) / sa;
		S.ust(i, j) = S.Xs(i, j, IMU) / rho; S.vst(i, j) = S.Xs(i, j, IMV) / rho;
	}
	auto rhoS = [&](int i, int j) { return S.Xs(i, j, IM1) + S.Xs(i, j, IM2); };
	// B4.1 声学项用的声速：冻结声速 c_f（各松弛闭合声速的上界，次特征条件，所以补足后的耗散不低于任何闭合的 Godunov 耗散）。
	// 在显式预测态 X* 上计算；X* 不可容许时取 c = inf，于是 phi = 0，该面不加声学项（只剩半格动量平衡）。
	// （早先试过的“半格质量取 rho min(c tau, dx/2)”的面平衡已弃用：大步长时耗散随 a_ii 变化，见 docs/devlog.md V3。）
	auto cfz = [&](int i, int j) -> double {
		const double a1 = S.Xs(i, j, IA1), a2 = S.a2s(i, j), m1 = S.Xs(i, j, IM1), m2 = S.Xs(i, j, IM2);
		if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return INF;
		const double p1 = p_ame(0, a1, m1, S.Xs(i, j, IE1)), p2 = p_ame(1, a2, m2, S.Xs(i, j, IE2));
		if (!(p1 + PH[0].pi > 0 && p2 + PH[1].pi > 0)) return INF;
		const double c21 = c2_rp(0, m1 / a1, p1), c22 = c2_rp(1, m2 / a2, p2);
		const double c2 = (m1 * c21 + m2 * c22) / (m1 + m2);
		return c2 > 0 ? std::sqrt(c2) : INF; };
	std::vector<double> cz((size_t)NXT * NYT, INF);
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) cz[S.id(i, j)] = cfz(i, j);
	// 残量尺度（B8）：Wood 体积模量 K = rho c_W^2（G/K 即体积约束残量）+ 能量的压力尺度 P_E = sum(|eps_k| + |m_k q_k|) / sum de_k/dp
	// （守恒总能的舍入误差 eps_mach |E| 经 lambda_k delta 进入各相压力，故 G 的可达精度是 eps_mach P_E）
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		const double a1 = S.Xs(i, j, IA1), a2 = S.a2s(i, j), m1 = S.Xs(i, j, IM1), m2 = S.Xs(i, j, IM2), sa = a1 + a2;
		double K = std::max(PH[0].pi, PH[1].pi), PE = 0.0;
		if (a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0) {
			const double p1 = p_ame(0, a1, m1, S.Xs(i, j, IE1)), p2 = p_ame(1, a2, m2, S.Xs(i, j, IE2));
			if (p1 + PH[0].pi > 0 && p2 + PH[1].pi > 0) K = 1.0 / ((a1 / sa) / (m1 / a1 * c2_rp(0, m1 / a1, p1)) + (a2 / sa) / (m2 / a2 * c2_rp(1, m2 / a2, p2)));
			const double dd = dedp_am(0, a1, m1) + dedp_am(1, a2, m2);
			if (dd > 0) PE = (std::fabs(S.Xs(i, j, IE1)) + std::fabs(m1 * PH[0].q) + std::fabs(S.Xs(i, j, IE2)) + std::fabs(m2 * PH[1].q)) / dd;
		}
		S.Kw[S.id(i, j)] = K + PE;
	}
	// 上一次的压力场（Newton 初猜）与 u*：只用来算二阶限制因子 psi（显式量，残量对 P 仍是线性的）
	Array Pold(nx, ny);
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		int ii = i - g, jj = j - g;
		if (x_bc == x_Periodic) ii = (ii % nx + nx) % nx; else ii = std::min(std::max(ii, 0), nx - 1);
		if (y_bc == y_Periodic && ny > 1) jj = (jj % ny + ny) % ny; else jj = std::min(std::max(jj, 0), ny - 1);
		Pold(i, j) = P.at(ii, jj);
	}
	auto minmod = [](double a, double b) { return a * b <= 0 ? 0.0 : (std::fabs(a) < std::fabs(b) ? a : b); };
	// psi = 1 - (sig_L + sig_R)/(2 D)，sig 为 minmod 斜率：与 D 同号且 |sig| <= |D|，故 psi in [0,1]；光滑区 psi = O(dx)（耗散为二阶），间断处 psi -> 1（BB (35) 的做法：耗散作用在重构跳跃上）
	auto psi_face = [&](double xm, double xL, double xR, double xp) {
		if (SpaceOrder < 2) return 1.0;
		const double D = xR - xL; if (D == 0.0) return 1.0;
		return 1.0 - 0.5 * (minmod(xL - xm, D) + minmod(D, xp - xR)) / D; };
	// B4.1：u_f = ubar* - tau dP/(rhobar h) - phi psi_P dP/(Z_L + Z_R)，p_f = (rho_R P_L + rho_L P_R)/(rho_L + rho_R) - phi psi_u Z_L Z_R/(Z_L + Z_R) du*，
	// phi = (1 - c tau/h)^+：后向 Euler 的声学耗散 c^2 tau/2 不足声学 Godunov 的 c h/2 时补足差额（c tau/h >= 1 时为 0，低马赫/大步长不受影响）
	auto set_face = [&](int iL, int jL, int iR, int jR, double h, bool ydir, double& ubs, double& kp, double& wL, double& pv) {
		const double rL = rhoS(iL, jL), rR = rhoS(iR, jR);
		const double uL = ydir ? S.vst(iL, jL) : S.ust(iL, jL), uR = ydir ? S.vst(iR, jR) : S.ust(iR, jR);
		ubs = 0.5 * (uL + uR); kp = tau / h / (0.5 * (rL + rR)); wL = rR / (rL + rR); pv = 0.0;
		const double cL = cz[S.id(iL, jL)], cR = cz[S.id(iR, jR)];
		if (!(std::isfinite(cL) && std::isfinite(cR))) return;
		const double phi = std::max(0.0, 1.0 - tau * std::min(cL, cR) / h);
		if (phi == 0.0) return;
		const int di = iR - iL, dj = jR - jL, im = iL - di, jm = jL - dj, ip = iR + di, jp = jR + dj;
		double psiP = 1.0, psiU = 1.0;
		if (im >= 0 && jm >= 0 && ip < NXT && jp < NYT) {
			// 只在所有被输运量都光滑处（共用限制因子 phi > 0，B3）才把耗散降到二阶：接触间断处 p、u 连续但阻抗跳跃，phi = 0 保持一阶耗散
			const std::vector<double>& ph = ydir ? phy : phx;
			const double phis = std::min(ph[S.id(iL, jL)], ph[S.id(iR, jR)]);
			psiP = 1.0 - (1.0 - psi_face(Pold(im, jm), Pold(iL, jL), Pold(iR, jR), Pold(ip, jp))) * phis;
			const Array& V = ydir ? S.vst : S.ust;
			psiU = 1.0 - (1.0 - psi_face(V(im, jm), uL, uR, V(ip, jp))) * phis;
		}
		const double ZL = rL * cL, ZR = rR * cR;
		kp += phi * psiP / (ZL + ZR); pv = phi * psiU * ZL * ZR / (ZL + ZR) * (uR - uL);
	};
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT - 1; i++) set_face(i, j, i + 1, j, dx, false, S.ubsx(i, j), S.kpx(i, j), S.wLx(i, j), S.pvx(i, j));
	if (ny > 1) for (int j = 0; j < NYT - 1; j++) for (int i = 0; i < NXT; i++) set_face(i, j, i, j + 1, dy, true, S.ubsy(i, j), S.kpy(i, j), S.wLy(i, j), S.pvy(i, j));
	for (int j = g; j < ny + g; j++) for (int i = g; i < nx + g; i++) {
		const double* c = &S.cont[S.id(i, j) * 6];
		const double rho = S.Xs(i, j, IM1) + S.Xs(i, j, IM2), sa = S.Xs(i, j, IA1) + S.a2s(i, j);
		const double rhoe = (S.Xs(i, j, IEN) - 0.5 * (S.Xs(i, j, IMU) * S.Xs(i, j, IMU) + S.Xs(i, j, IMV) * S.Xs(i, j, IMV)) / rho) / sa;
		S.type[S.id(i, j)] = STEP_TYPE.empty() ? closure_type_of(c[0], c[2], c[3], c[4], c[5], rhoe) : STEP_TYPE[(size_t)(i - g) + (size_t)(j - g) * nx];
		if (S.type[S.id(i, j)] == CL_MU) {
			for (int k = 0; k < 2; k++) { const double a = c[k], m = c[2 + k], pk = p_ame(k, a, m, c[4 + k]), r = m / a; const double Z = r * std::sqrt(std::max(c2_rp(k, r, pk), 1e-30)); (k == 0 ? S.Z0 : S.Z1)[S.id(i, j)] = Z; }
		}
	}
	// 隐式：全局 Newton
	int iters; double res;
	st_newton_calls++;
	const bool ok = newton_solve(S, P, iters, res);
	if (ok) st_maxres = std::max(st_maxres, res);
	if (!ok) { if (VERBOSE) fprintf(stderr, "  Newton failed: it %d res %.3e tau %.3e\n", iters, res, tau); return false; }
	// 终态
	U = Qs;
	for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
		double Gv; CellOut o; if (!residual_at(S, P, i, j, 0, -1, 0, Gv, &o)) return false;
		const int I = IG(i), J = IG(j);
		U(I, J, IA1) = o.cl.a[0]; U(I, J, IM1) = o.cl.m[0]; U(I, J, IM2) = o.cl.m[1]; U(I, J, IMU) = o.mx; U(I, J, IMV) = o.my;
		U(I, J, IE1) = o.cl.e[0]; U(I, J, IE2) = o.cl.e[1]; U(I, J, IEN) = o.E;
	}
	bound_cond_u(U);
	st_stage++;
	return true;
}

// ==================================================================================================== main（2DEuler_TV 的时间循环）
int main(int argc, char** argv) {
	CommonArgs A;
	for (int a = 1; a < argc; a++) {
		char key[64]; char val[256]; if (sscanf(argv[a], "%63[^=]=%255s", key, val) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); return 1; }
		std::string k(key); const double d = atof(val);
		if (parse_common(k, val, d, A)) continue;
		if (k == "rk") RKs = (int)d; else if (k == "cflac") CFL_ac = d; else if (k == "ramp") RAMP = (int)d; else if (k == "tol") NEWTON_TOL = d; else if (k == "dt") DTFIX = d;
		else { fprintf(stderr, "unknown key %s\n", key); return 1; }
	}
	setup_common(A);
	Grid U(nxCells, nyCells); init_case(U);
	double M0[3]; conserved_totals(U, M0);

	// Butcher 表（2DEuler_TV）
	double Ae[3][3] = {{0}}, Ai[3][3] = {{0}};
	auto set_tab = [&](int rk) {
		for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) Ae[a][b] = Ai[a][b] = 0;
		if (rk == 1) { Ai[0][0] = 1.0; }
		else if (rk == 2) { const double bb = 1.0 - 1.0 / std::sqrt(2.0), cc = 1.0 / (2.0 * bb); Ae[1][0] = cc; Ai[0][0] = bb; Ai[1][0] = 1 - bb; Ai[1][1] = bb; }
		else { Ae[1][0] = 0.5; Ae[2][0] = 0.5; Ae[2][1] = 0.5; Ai[0][0] = 0.25; Ai[1][1] = 0.25; Ai[2][0] = Ai[2][1] = Ai[2][2] = 1.0 / 3.0; }
	};
	auto t0 = std::chrono::steady_clock::now();
	t = 0; double dtprev = -1; long nrk1 = 0;
	PField P{nxCells, nyCells, std::vector<double>((size_t)nxCells * nyCells, 0.0)};
	for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) P.at(i, j) = pbar_cell(U, IG(i), IG(j));
	while (t < time_stop * (1 - 1e-12)) {
		bool ramping = false;
		if (DTFIX > 0) dt = DTFIX;
		else {
			dt = CFL * caldt(U, false);
			if (RAMP) { if (dtprev < 0) dtprev = CFL_ac * caldt(U, true) / 1.2; if (dt > 1.2 * dtprev) { dt = 1.2 * dtprev; ramping = true; } }
			if (!std::isfinite(dt)) dt = CFL_ac * caldt(U, true);
		}
		dt = std::min(dt, time_stop - t);
		if (!(dt > 1e-12 * time_stop)) { fprintf(stderr, "FAILED: time step underflow (dt = %.3e) at t = %.6e (step %ld)\n", dt, t, st_steps); output(U, A.out.c_str()); return 2; }
		bool done = false;
		for (int attempt = 0; attempt < 30 && !done; attempt++) {
			const int rk = (RAMP && ramping) ? 1 : RKs; if (rk == 1) nrk1++;
			set_tab(rk);
			if (RELAX_LEVEL > 0) { activate_and_project(U); bound_cond_u(U); }
			Grid Un = U; std::vector<Grid> Hs(rk, Grid(nxCells, nyCells)); Grid Qe = U, Qs = U, Ust = U; PField Pst = P;
			bool ok = true;
			for (int rks = 0; rks < rk && ok; rks++) {
				Qe = Un; Qs = Un;
				for (int r = 0; r < rks; r++) for (size_t q = 0; q < Qe.data.size(); q++) { Qe.data[q] += Ae[rks][r] * dt * Hs[r].data[q]; Qs.data[q] += Ai[rks][r] * dt * Hs[r].data[q]; }
				bound_cond_u(Qe); bound_cond_u(Qs);
				const double adt = Ai[rks][rks] * dt;
				ok = imex_stage(Qe, Qs, adt, Ust, Pst);
				if (ok) for (size_t q = 0; q < Ust.data.size(); q++) Hs[rks].data[q] = (Ust.data[q] - Qs.data[q]) / adt;
			}
			if (ok) { U = Ust; P = Pst; done = true; }
			if (!done) { st_reject++; dt *= 0.5; if (VERBOSE) fprintf(stderr, "reject: t %.4e new dt %.3e\n", t, dt); }
		}
		if (!done) { fprintf(stderr, "FAILED at t = %.6e (step %ld)\n", t, st_steps); output(U, A.out.c_str()); return 2; }
		t += dt; dtprev = dt; st_steps++;
		if (VERBOSE > 1 && st_steps % 50 == 0) fprintf(stderr, "step %ld t %.4e dt %.3e newton it %ld\n", st_steps, t, dt, st_newton_it);
	}
	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double M1[3]; conserved_totals(U, M1);
	if (RELAX_LEVEL > 0) { activate_and_project(U); bound_cond_u(U); }   // 输出前同样做激活投影（P-S：输出的是松弛后的态）
	output(U, A.out.c_str());
	printf("case %d mode 1 order %d rk %d N %d CFL %g | steps %ld (rk1 %ld) rejects %ld | stages %ld newton it %ld (%.2f/stage) maxres %.2e | closure fails %ld | dM1/M %.2e dM2/M %.2e dE/E %.2e | wall %.2f s\n",
		TESTCASE, SpaceOrder, RKs, nxCells, CFL, st_steps, nrk1, st_reject, st_stage, st_newton_it, st_stage ? (double)st_newton_it / st_stage : 0.0, st_maxres, st_cl_fail,
		(M1[0] - M0[0]) / M0[0], (M1[1] - M0[1]) / M0[1], (M1[2] - M0[2]) / M0[2], wall);
	return 0;
}
