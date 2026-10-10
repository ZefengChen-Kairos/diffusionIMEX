// =====================================================================================================================
//  SPB6_TV_IMEX.cpp
//  Saurel-Petitpas-Berry (JCP 228 (2009) 1678) 6+1 方程两相模型（单速度、双压力 + 冗余总能），
//  加 Pelanti-Shyue (2014) / Pelanti (2022) 的热、化学松弛。代码结构仿照 2DEuler_TV.cpp。
//
//  MODE 0  全显式参考格式（docs/derivation/A_model_EOS_explicit.md）：
//          HLLC + Godunov（SPB 3.1-3.2），力学松弛 + 用总能重初始化（SPB 3.3-3.4），P-S/Pelanti 的热/化学松弛。
//  MODE 1  TV 分裂 IMEX（docs/derivation/B_TV_IMEX.md）：
//          显式 = 体积输运形式的 MUSCL-Rusanov（与线性路径、面速度 ubar 的 path-conservative Rusanov 等价），
//                 耗散速度取物质速度（BB 2025 的做法；BB (36)-(39) 都包含在对 alpha_k, m_k, eps_k, rho u, E 的 Rusanov 耗散里）；
//          面   = 面速度/面压力（B4.1）：半格动量平衡 + 声学 Riemann 解的迎风耗散，权重 (1 - c tau/dx)^+；
//                 增量通量对总面速度 u_f 迎风（B4.2）；
//          隐式 = 每格一个未知量（混合压力 P）的非线性方程 G(P) = pbar_cl(P) - P = 0，
//                 pbar_cl 由统一的局部闭合给出：离散 Gibbs 关系 + 后向 Euler 的松弛（速率 mu, theta, nu 为参数）
//                 + 混合能量一致性（SPB 的 "+1" 在闭合内）。
//          时间 = 半隐式 IMEX-RK（2DEuler_TV 的 Ae/Ai/Hs 结构），RKs = 1, 2 (LSDIRK2), 3 (SSP2(3,3,2))；
//                 瞬时热/化学松弛的激活每步判断一次并投影到平衡流形（B6.7，保证 IMEX-RK 的初值在流形上）。
//
//  守恒变量（NUM_VARS = 8）: 0 alpha1, 1 m1=(alpha rho)_1, 2 m2, 3 rho u, 4 rho v, 5 eps1=(alpha rho e)_1, 6 eps2, 7 E
//  eps_k 含生成焓 q_k；E = eps1 + eps2 + 1/2 rho |u|^2。
//
//  编译: g++ -O2 -std=c++17 SPB6_TV_IMEX.cpp -o spb6
//  运行: ./spb6 key=value ...   (见 main 里的参数表)
// =====================================================================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <limits>
#include <chrono>

#define NUM_VARS 8
enum { IA1 = 0, IM1 = 1, IM2 = 2, IMU = 3, IMV = 4, IE1 = 5, IE2 = 6, IEN = 7 };
static const double INF = std::numeric_limits<double>::infinity();

// ------------------------------------------------------------------ 全局参数（2DEuler_TV 风格：全局变量，main 里可由命令行覆盖）
enum x_boundary_condition { x_Transmissive, x_Reflective, x_Periodic };
enum y_boundary_condition { y_Transmissive, y_Reflective, y_Periodic };
x_boundary_condition x_bc = x_Transmissive;
y_boundary_condition y_bc = y_Transmissive;

int nxCells = 1000, nyCells = 1;
const int ghost_cells = 2;
int MODE = 1;               // 0 显式参考, 1 IMEX
int RKs = 2;                // IMEX-RK 级数: 1 (后向/前向 Euler), 2 (LSDIRK2(2,2,2)), 3 (SSP2(3,3,2))
int SpaceOrder = 2;         // 1 或 2（minmod）
double CFL = 0.2;           // MODE 1: 物质 CFL; MODE 0: 声学 CFL
double CFL_ac = 0.4;        // 起步用的声学 CFL
int RAMP = 1;               // 从静止/间断起算：首步声学步长，之后每步至多放大 1.2 倍；斜坡期间用一阶（后向 Euler）
double time_stop = 1e-3, t = 0, dt = 0;
double x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1, dx = 1, dy = 1;
int TESTCASE = 5;

// 松弛设置
double mu_rate = INF, theta_rate = 0.0, nu_rate = 0.0;   // 力学 / 热 / 化学 速率（INF = 瞬时）
int RELAX_LEVEL = 0;        // 0: p;  1: p-pT;  2: p-pT-pTg;  3: p-pTg（P-S 的第四档）
double EPS_I = 0.0;         // 热/化学松弛只在 eps_I <= alpha1 <= 1 - eps_I 的格子打开
int SPB_EXACT = 1;          // MODE 0 的力学松弛: 1 = SPB 原文（松弛后用总能重初始化，分两步）; 0 = 统一闭合（s = 0）
// 以下开关只用于对照试验（默认值即推导 B 中的格式）
int CENTRAL_CELLAVG = 0;     // 1 = BB 2025 的中心部分用格均值（非 TVD，B3）
int ISEN = 0;                // 1 = 体积应变按各相精确等熵线 (p+pi)(alpha-mb)^gamma = const（对照用）
int ACOUSTIC = 1;            // 1 = 声学 Riemann 解的迎风耗散，权重 (1 - c tau/dx)^+（B4.1）；0 = 只有 BB 的物质速度耗散
int ACORD = 2;               // 声学耗散的空间阶：2 = 乘限制因子 psi（B4.1）；1 = 用格心跳跃
int CSPEED = 0;              // 声学耗散用的声速：0 冻结（满足所有闭合的次特征条件）；1 Wood（对照用）
int INCDONOR = 1;            // 增量通量：1 = 总通量对 u_f 迎风（B4.2）；0 = 按 du 的符号取供体（旧，有相间逆流）
int STEPACT = 1;             // 1 = 每步开始按 U^n 激活并投影到平衡流形、级内冻结类型（B6.7）；0 = 每级按预测态判断（旧）
double NEWTON_TOL = 1e-13;   // 相对 |P| + rho c_W^2 + P_E
int NEWTON_MAX = 60;
int VERBOSE = 0;
int DUMP = 0;                // 调试：前 DUMP 步每步输出 out.stepN

// 统计
long st_steps = 0, st_reject = 0, st_newton_it = 0, st_newton_calls = 0, st_cl_fail = 0, st_stage = 0;
double st_maxres = 0;

// ------------------------------------------------------------------ 数据结构（同 2DEuler_TV：AoS，i 为快指标；用 vector 避免手写拷贝）
struct Grid {
	int nx, ny; std::vector<double> data;
	Grid(int x = 1, int y = 1) : nx(x), ny(y), data((size_t)(x + 2 * ghost_cells) * (y + 2 * ghost_cells) * NUM_VARS, 0.0) {}
	inline double  operator()(int i, int j, int k) const { return data[NUM_VARS * ((size_t)i + (size_t)j * (nx + 2 * ghost_cells)) + k]; }
	inline double& operator()(int i, int j, int k)       { return data[NUM_VARS * ((size_t)i + (size_t)j * (nx + 2 * ghost_cells)) + k]; }
};
struct Array {
	int nx, ny; std::vector<double> data;
	Array(int x = 1, int y = 1, double v = 0.0) : nx(x), ny(y), data((size_t)(x + 2 * ghost_cells) * (y + 2 * ghost_cells), v) {}
	inline double  operator()(int i, int j) const { return data[(size_t)i + (size_t)j * (nx + 2 * ghost_cells)]; }
	inline double& operator()(int i, int j)       { return data[(size_t)i + (size_t)j * (nx + 2 * ghost_cells)]; }
};
struct Tensor {   // 各级的级导数 Hs(r,i,j,k)
	std::vector<Grid> g;
	Tensor(int r, int x, int y) : g(r, Grid(x, y)) {}
	inline double& operator()(int r, int i, int j, int k) { return g[r](i, j, k); }
};
inline void copy(Grid& a, const Grid& b) { a.data = b.data; }
inline int IG(int i) { return i + ghost_cells; }   // 内部格 i (0-based) -> 带虚格下标

// ==================================================================================================== 状态方程（NASG；b = 0 为 SG）
struct Phase { double g, pi, q, qp, cv, b; };
Phase PH[2];

inline double eps_amp(int k, double a, double m, double p) { const Phase& P = PH[k]; return (p + P.g * P.pi) * (a - m * P.b) / (P.g - 1.0) + m * P.q; }   // (A2.2)
inline double p_ame(int k, double a, double m, double e) { const Phase& P = PH[k]; return (e - m * P.q) * (P.g - 1.0) / (a - m * P.b) - P.g * P.pi; }
inline double T_amp(int k, double a, double m, double p) { const Phase& P = PH[k]; return (a / m - P.b) * (p + P.pi) / ((P.g - 1.0) * P.cv); }
inline double v_pT(int k, double p, double T) { const Phase& P = PH[k]; return (P.g - 1.0) * P.cv * T / (p + P.pi) + P.b; }
inline double e_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.cv * T * (p + P.g * P.pi) / (p + P.pi) + P.q; }
inline double s_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.cv * (P.g * std::log(T) - (P.g - 1.0) * std::log(p + P.pi)) + P.qp; }
inline double g_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.g * P.cv * T + P.b * p + P.q - T * s_pT(k, p, T); }
inline double c2_rp(int k, double r, double p) { const Phase& P = PH[k]; return P.g * (p + P.pi) / (r * (1.0 - r * P.b)); }
inline double dedp_am(int k, double a, double m) { const Phase& P = PH[k]; return (a - m * P.b) / (P.g - 1.0); }    // d eps_k / d p |_{alpha,m}

// 混合压力 (A2.3)：给定 alpha_k, m_k 与混合内能
inline double p_mix(double a1, double m1, double m2, double rhoe) {
	const double a[2] = {a1, 1.0 - a1}, m[2] = {m1, m2};
	double num = rhoe, den = 0.0;
	for (int k = 0; k < 2; k++) { num -= m[k] * PH[k].q + PH[k].g * PH[k].pi * dedp_am(k, a[k], m[k]); den += dedp_am(k, a[k], m[k]); }
	return num / den;
}

// ------------------------------------------------------------------ 饱和线 Tsat(p)：g1(p,T) = g2(p,T)，取 s2 > s1 的物理分支（f = g1 - g2 对 T 递增）
struct SatTable { bool built = false; std::vector<double> lp, T; double pcrit = INF, vcrit = 0.0, pmin = 0.0; } SAT;
static double sat_f(double p, double T) { return g_pT(0, p, T) - g_pT(1, p, T); }
static double sat_newton(double p, double Tg, bool& ok) {
	double T = Tg; ok = false;
	for (int it = 0; it < 60; it++) {
		const double f = sat_f(p, T), d = s_pT(1, p, T) - s_pT(0, p, T);
		if (!(d > 0) || !std::isfinite(f)) return T;
		double Tn = T - f / d; if (Tn < 0.5 * T) Tn = 0.5 * T; if (Tn > 2.0 * T) Tn = 2.0 * T;
		if (std::fabs(Tn - T) < 1e-14 * T) { ok = true; return Tn; }
		T = Tn;
	}
	ok = std::fabs(sat_f(p, T)) < 1e-8 * std::fabs(g_pT(0, p, T)) + 1e-8; return T;
}
static void build_sat() {
	SAT = SatTable(); SAT.built = true;
	const int NP = 1600; const double l0 = std::log(1e-1), l1 = std::log(5e9);
	double Tprev = -1.0;
	for (int n = 0; n < NP; n++) {
		const double lp = l0 + (l1 - l0) * n / (NP - 1), p = std::exp(lp);
		double Troot = -1.0;
		if (Tprev < 0) {   // 第一个点：在 [20, 5000] K 上找 f 由负变正的根
			double Ta = 20.0, fa = sat_f(p, Ta);
			for (int m = 1; m <= 4000; m++) { const double Tb = 20.0 * std::pow(250.0, m / 4000.0), fb = sat_f(p, Tb);
				if (fa < 0 && fb >= 0) { double a = Ta, b = Tb; for (int it = 0; it < 200; it++) { const double c = 0.5 * (a + b); if (sat_f(p, c) < 0) a = c; else b = c; } Troot = 0.5 * (a + b); break; }
				Ta = Tb; fa = fb; }
		} else { bool ok; const double T = sat_newton(p, Tprev, ok); if (ok) Troot = T; }
		if (Troot < 0) { if (n > 0) SAT.pcrit = std::exp(SAT.lp.back()); break; }
		if (!(v_pT(1, p, Troot) > v_pT(0, p, Troot))) { SAT.pcrit = p; SAT.vcrit = v_pT(0, p, Troot); break; }   // 端点 v_l = v_v
		SAT.lp.push_back(lp); SAT.T.push_back(Troot); Tprev = Troot;
	}
	if (!SAT.lp.empty()) SAT.pmin = std::exp(SAT.lp.front());
	if (SAT.vcrit == 0.0 && !SAT.T.empty()) SAT.vcrit = v_pT(0, std::exp(SAT.lp.back()), SAT.T.back());
	// 端点处二分细化
	if (std::isfinite(SAT.pcrit) && !SAT.lp.empty()) {
		double a = std::exp(SAT.lp.back()), b = SAT.pcrit, Ta = SAT.T.back();
		for (int it = 0; it < 100; it++) { const double c = std::sqrt(a * b); bool ok; const double T = sat_newton(c, Ta, ok);
			if (ok && v_pT(1, c, T) > v_pT(0, c, T)) { a = c; Ta = T; } else b = c; }
		SAT.pcrit = a; SAT.vcrit = v_pT(0, a, Ta);
	}
	if (VERBOSE) fprintf(stderr, "Tsat table: %zu pts, p in [%.3e, %.3e], pcrit %.4e, vcrit %.4e\n", SAT.lp.size(), SAT.pmin, std::exp(SAT.lp.back()), SAT.pcrit, SAT.vcrit);
}
static double Tsat(double p) {
	if (!SAT.built) build_sat();
	if (!(p > 0) || !(p < SAT.pcrit) || SAT.lp.empty()) return NAN;
	const double lp = std::log(p); double Tg;
	if (lp <= SAT.lp.front()) Tg = SAT.T.front();
	else if (lp >= SAT.lp.back()) Tg = SAT.T.back();
	else { const double h = SAT.lp[1] - SAT.lp[0]; const int n = std::min((int)((lp - SAT.lp[0]) / h), (int)SAT.lp.size() - 2); const double w = (lp - SAT.lp[n]) / h; Tg = (1 - w) * SAT.T[n] + w * SAT.T[n + 1]; }
	bool ok; const double T = sat_newton(p, Tg, ok); return ok ? T : NAN;
}

// ==================================================================================================== 标量求根
// Brent（带区间，fa, fb 异号）
static bool brent(const std::function<double(double)>& f, double a, double b, double fa, double fb, double tol, double& root) {
	if (!(fa * fb <= 0)) return false;
	if (fa == 0) { root = a; return true; } if (fb == 0) { root = b; return true; }
	double c = a, fc = fa, d = b - a, e = d;
	for (int it = 0; it < 200; it++) {
		if (fb * fc > 0) { c = a; fc = fa; d = b - a; e = d; }
		if (std::fabs(fc) < std::fabs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
		const double tol1 = 2.0 * 1e-16 * std::fabs(b) + 0.5 * tol, xm = 0.5 * (c - b);
		if (std::fabs(xm) <= tol1 || fb == 0) { root = b; return true; }
		if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
			double s = fb / fa, p, q;
			if (a == c) { p = 2 * xm * s; q = 1 - s; }
			else { const double qq = fa / fc, r = fb / fc; p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
			if (p > 0) q = -q; p = std::fabs(p);
			if (2 * p < std::min(3 * xm * q - std::fabs(tol1 * q), std::fabs(e * q))) { e = d; d = p / q; } else { d = xm; e = d; }
		} else { d = xm; e = d; }
		a = b; fa = fb;
		b += (std::fabs(d) > tol1) ? d : (xm > 0 ? tol1 : -tol1);
		fb = f(b); if (!std::isfinite(fb)) return false;
	}
	root = b; return true;
}
// 在 (plo, inf) 上对 y = ln(p - plo) 找变号区间再 Brent；f 在 plo+ 处应为正（或 NaN 表示不可行），在 +inf 处为负
static bool root_above(const std::function<double(double)>& f, double plo, double guess, double& root) {
	auto F = [&](double y) { return f(plo + std::exp(y)); };
	double y0 = std::log(std::max(guess - plo, 1e-6 * (std::fabs(plo) + 1.0) + 1e-12));
	double f0 = F(y0);
	if (f0 == 0.0) { root = plo + std::exp(y0); return true; }   // 初猜恰为根
	double ya = y0, yb = y0, fa = f0, fb = f0, st = 0.1;
	if (!std::isfinite(f0)) {   // 往上找可行点
		for (int k = 0; k < 200 && !std::isfinite(fb); k++) { yb += 1.0; fb = F(yb); }
		if (!std::isfinite(fb)) return false; ya = yb; fa = fb;
	}
	if (fa > 0) { for (int k = 0; k < 400 && !(fb < 0); k++) { ya = yb; fa = fb; yb += st; st *= 1.6; fb = F(yb); if (!std::isfinite(fb)) return false; } }
	else { for (int k = 0; k < 400 && !(fa > 0); k++) { yb = ya; fb = fa; ya -= st; st *= 1.6; fa = F(ya); if (!std::isfinite(fa)) { ya = yb; fa = fb; for (int q = 0; q < 60; q++) { const double ym = ya - st / std::pow(2.0, q + 1); const double fm = F(ym); if (std::isfinite(fm)) { if (fm > 0) { ya = ym; fa = fm; break; } } } break; } } }
	if (fa == 0.0) { root = plo + std::exp(ya); return true; }
	if (fb == 0.0) { root = plo + std::exp(yb); return true; }
	if (!(fa > 0 && fb < 0)) return false;
	double y; if (!brent(F, ya, yb, fa, fb, 1e-15, y)) return false;
	root = plo + std::exp(y); return true;
}

// ==================================================================================================== 统一局部闭合（B5, B6）
enum ClosureType { CL_P = 0, CL_MU = 1, CL_PT = 2, CL_THETA = 3, CL_NU = 4, CL_PTG = 5 };
struct ClIn {
	double atr[2], mtr[2], etr[2];   // 输运态（未压回单位体积）
	double s;                        // 体积应变 s = 1 - sum atr
	double rhoeC;                    // 守恒的混合内能 E - 1/2 rho u^2
	double tau;                      // 级的隐式步长 a_ii dt
	double Z[2];                     // 有限 mu 时的阻抗（冻结）
	int type;
	bool nodelta;                    // MODE 0 的 SPB 原文力学松弛：不含 lambda*delta（之后单独重初始化）
	double pguess;
};
struct ClOut { double a[2], m[2], e[2], p[2], T[2], pbar; };
static const double YMIN = 1e-10;

// (B6.2)：p 闭合中 alpha_k(p) 的系数
struct PCoef { double A[2], B[2], N[2], lam[2]; };
static bool p_coef(const ClIn& in, PCoef& c) {
	double den = 0.0;
	for (int k = 0; k < 2; k++) { c.lam[k] = in.nodelta ? 0.0 : dedp_am(k, in.atr[k], in.mtr[k]); den += c.lam[k]; }
	if (!in.nodelta) { if (!(den > 0)) return false; for (int k = 0; k < 2; k++) c.lam[k] /= den; }
	const double D0 = in.nodelta ? 0.0 : in.rhoeC - in.etr[0] - in.etr[1];
	for (int k = 0; k < 2; k++) {
		c.A[k] = in.etr[k] - in.mtr[k] * PH[k].q + c.lam[k] * D0;
		c.B[k] = in.atr[k] + c.lam[k] * in.s;
		c.N[k] = c.A[k] - (c.B[k] - in.mtr[k] * PH[k].b) * PH[k].pi;
	}
	return true;
}
inline double alpha_p(const PCoef& c, const ClIn& in, int k, double p) { const Phase& P = PH[k]; return ((P.g - 1.0) * (c.A[k] + c.B[k] * p) + (p + P.g * P.pi) * in.mtr[k] * P.b) / (P.g * (p + P.pi)); }
inline double dalpha_p(const PCoef& c, int k, double p) { const Phase& P = PH[k]; return -(P.g - 1.0) * c.N[k] / (P.g * (p + P.pi) * (p + P.pi)); }

static double plo_mix() { return std::max(-PH[0].pi, -PH[1].pi); }

// B6.1  mu = inf：sum_k alpha_k(p) = 1，严格递减，安全 Newton
static bool closure_p(const ClIn& in, ClOut& o) {
	PCoef c; if (!p_coef(in, c)) return false;
	for (int k = 0; k < 2; k++) if (!(c.N[k] > 0)) return false;
	const double plo = plo_mix();
	auto Phi = [&](double p) { return alpha_p(c, in, 0, p) + alpha_p(c, in, 1, p) - 1.0; };
	// 区间: lo 处 Phi > 0（p -> plo+ 时 -> +inf），找 hi 使 Phi < 0
	double d = std::max(in.pguess - plo, 1e-3 * (std::fabs(plo) + 1.0));
	double lo = plo, hi = plo + d; double fh = Phi(hi);
	if (fh > 0) { lo = hi; for (int k = 0; k < 200 && fh > 0; k++) { d *= 4.0; hi = plo + d; fh = Phi(hi); } if (fh > 0) return false; }
	else { for (int k = 0; k < 200; k++) { const double dd = d * 0.25, pl = plo + dd, fl = Phi(pl); if (fl > 0) { lo = pl; break; } hi = pl; fh = fl; d = dd; } }
	double p = 0.5 * (lo + hi); if (in.pguess > lo && in.pguess < hi) p = in.pguess;
	for (int it = 0; it < 200; it++) {
		const double f = Phi(p);
		if (f == 0.0) break;                                   // 恰为根
		if (f > 0) lo = p; else hi = p;
		const double df = dalpha_p(c, 0, p) + dalpha_p(c, 1, p);
		double pn = p - f / df;
		if (!(pn > lo && pn < hi)) pn = (hi - plo > 1e3 * (lo - plo + 1e-300) && lo > plo) ? plo + std::sqrt((lo - plo) * (hi - plo)) : 0.5 * (lo + hi);
		if (std::fabs(pn - p) <= 2e-16 * (std::fabs(p) + std::fabs(plo)) || hi - lo <= 4e-16 * (std::fabs(hi) + std::fabs(lo))) { p = pn; break; }
		p = pn;
	}
	for (int k = 0; k < 2; k++) {
		o.a[k] = alpha_p(c, in, k, p); o.m[k] = in.mtr[k];
		if (!(o.a[k] > o.m[k] * PH[k].b)) return false;
		o.e[k] = eps_amp(k, o.a[k], o.m[k], p); o.p[k] = p; o.T[k] = T_amp(k, o.a[k], o.m[k], p);
	}
	o.pbar = p; return std::isfinite(p);
}

// 等熵版 mu = inf 闭合：(1) 各相从输运态沿等熵线 (p+pi)(alpha-mb)^gamma = C_k 到公共压力 p_s，sum alpha = 1；
//                      (2) 能量差 delta = rhoeC - sum eps(p_s) 按 lambda_k ~ d eps_k/dp|_alpha 分配（SPB 的重初始化），p = p_s + delta / sum lam
static bool closure_p_isen(const ClIn& in, ClOut& o) {
	double C[2], mb[2];
	for (int k = 0; k < 2; k++) {
		mb[k] = in.mtr[k] * PH[k].b; const double af = in.atr[k] - mb[k]; if (!(af > 0)) return false;
		const double ptr = p_ame(k, in.atr[k], in.mtr[k], in.etr[k]); if (!(ptr + PH[k].pi > 0)) return false;
		C[k] = (ptr + PH[k].pi) * std::pow(af, PH[k].g);
	}
	auto al = [&](int k, double p) { return mb[k] + std::pow(C[k] / (p + PH[k].pi), 1.0 / PH[k].g); };
	auto Phi = [&](double p) { return al(0, p) + al(1, p) - 1.0; };
	double ps; if (!root_above(Phi, plo_mix(), std::max(in.pguess, 1e-6), ps)) return false;
	double es = 0, lt = 0, a[2];
	for (int k = 0; k < 2; k++) { a[k] = al(k, ps); es += eps_amp(k, a[k], in.mtr[k], ps); lt += dedp_am(k, a[k], in.mtr[k]); }
	const double p = in.nodelta ? ps : ps + (in.rhoeC - es) / lt;
	for (int k = 0; k < 2; k++) {
		o.a[k] = a[k]; o.m[k] = in.mtr[k]; if (!(p + PH[k].pi > 0)) return false;
		o.e[k] = eps_amp(k, a[k], in.mtr[k], p); o.p[k] = p; o.T[k] = T_amp(k, a[k], in.mtr[k], p);
	}
	o.pbar = p; return std::isfinite(p);
}

// B6.2  有限 mu（theta = nu = 0）：给定 R 解 2x2 线性方程组得 (p1, p2)，再解 F(R) = R - tau mu (p1 - p2) = 0
static bool closure_mu(const ClIn& in, ClOut& o) {
	const double sa = in.atr[0] + in.atr[1];
	const double om[2] = {in.atr[0] / sa, in.atr[1] / sa};
	const double ak[2] = {in.atr[0] / sa, in.atr[1] / sa};   // a_k = atr_k + om_k s = atr_k / (1 - s)
	double lam[2], den = 0.0;
	for (int k = 0; k < 2; k++) { lam[k] = dedp_am(k, in.atr[k], in.mtr[k]); den += lam[k]; }
	for (int k = 0; k < 2; k++) lam[k] /= den;
	const double D0 = in.rhoeC - in.etr[0] - in.etr[1];
	const double w1 = in.Z[1] / (in.Z[0] + in.Z[1]), w2 = 1.0 - w1;
	auto solve = [&](double R, double& p1, double& p2) -> bool {
		const double a1 = ak[0] + R, a2 = ak[1] - R;
		const double C1 = dedp_am(0, a1, in.mtr[0]) + om[0] * in.s, C2 = dedp_am(1, a2, in.mtr[1]) + om[1] * in.s;
		const double M11 = C1 + R * w1 - lam[0] * in.s * om[0], M12 = R * w2 - lam[0] * in.s * om[1];
		const double M21 = -R * w1 - lam[1] * in.s * om[0], M22 = C2 - R * w2 - lam[1] * in.s * om[1];
		const double r1 = in.etr[0] - in.mtr[0] * PH[0].q - PH[0].g * PH[0].pi * dedp_am(0, a1, in.mtr[0]) + lam[0] * D0;
		const double r2 = in.etr[1] - in.mtr[1] * PH[1].q - PH[1].g * PH[1].pi * dedp_am(1, a2, in.mtr[1]) + lam[1] * D0;
		const double det = M11 * M22 - M12 * M21; if (!(std::fabs(det) > 0)) return false;
		p1 = (r1 * M22 - M12 * r2) / det; p2 = (M11 * r2 - M21 * r1) / det; return std::isfinite(p1) && std::isfinite(p2);
	};
	double R = 0.0, p1, p2;
	if (in.tau * mu_rate > 0) {
		const double Rlo = -(ak[0] - in.mtr[0] * PH[0].b) * (1 - 1e-13), Rhi = (ak[1] - in.mtr[1] * PH[1].b) * (1 - 1e-13);
		auto F = [&](double r) { double q1, q2; if (!solve(r, q1, q2)) return (double)NAN; return r - in.tau * mu_rate * (q1 - q2); };
		const double fl = F(Rlo), fh = F(Rhi);
		if (!(fl < 0 && fh > 0)) return false;
		if (!brent(F, Rlo, Rhi, fl, fh, 1e-16 * (Rhi - Rlo), R)) return false;
	}
	if (!solve(R, p1, p2)) return false;
	const double pp[2] = {p1, p2};
	for (int k = 0; k < 2; k++) {
		o.a[k] = ak[k] + (k == 0 ? R : -R); o.m[k] = in.mtr[k]; o.p[k] = pp[k];
		if (!(o.a[k] > o.m[k] * PH[k].b) || !(pp[k] + PH[k].pi > 0)) return false;
		o.e[k] = eps_amp(k, o.a[k], o.m[k], pp[k]); o.T[k] = T_amp(k, o.a[k], o.m[k], pp[k]);
	}
	o.pbar = o.a[0] * p1 + o.a[1] * p2; return true;
}

// B6.3  pT（mu = theta = inf, nu = 0）：T(p) 由能量显式定，sum m_k v_k(p, T) = 1
static bool closure_pT(const ClIn& in, ClOut& o) {
	const double m[2] = {in.mtr[0], in.mtr[1]};
	const double Eth = in.rhoeC - m[0] * PH[0].q - m[1] * PH[1].q;
	if (!(Eth > 0)) return false;
	auto Tof = [&](double p) { double d = 0; for (int k = 0; k < 2; k++) d += m[k] * PH[k].cv * (p + PH[k].g * PH[k].pi) / (p + PH[k].pi); return Eth / d; };
	auto Phi = [&](double p) { const double T = Tof(p); return m[0] * v_pT(0, p, T) + m[1] * v_pT(1, p, T) - 1.0; };
	double p; if (!root_above(Phi, plo_mix(), in.pguess, p)) return false;
	const double T = Tof(p);
	for (int k = 0; k < 2; k++) { o.m[k] = m[k]; o.a[k] = m[k] * v_pT(k, p, T); o.e[k] = m[k] * e_pT(k, p, T); o.p[k] = p; o.T[k] = T; }
	o.pbar = p; return T > 0;
}

// B6.4  有限 theta（mu = inf, nu = 0）：Delta Q(p) 闭式 (B6.5)
static bool closure_theta(const ClIn& in, ClOut& o) {
	PCoef c; if (!p_coef(in, c)) return false;
	for (int k = 0; k < 2; k++) if (!(c.N[k] > 0)) return false;
	const double Cp[2] = {PH[0].g * PH[0].cv * in.mtr[0], PH[1].g * PH[1].cv * in.mtr[1]};
	const double tt = in.tau * theta_rate;
	auto eval = [&](double p, double a[2], double& dQ) {
		double a0[2], T0[2];
		for (int k = 0; k < 2; k++) { a0[k] = alpha_p(c, in, k, p); T0[k] = T_amp(k, a0[k], in.mtr[k], p); }
		dQ = tt * (T0[1] - T0[0]) / (1.0 + tt * (1.0 / Cp[0] + 1.0 / Cp[1]));
		a[0] = a0[0] + (PH[0].g - 1.0) * dQ / (PH[0].g * (p + PH[0].pi));
		a[1] = a0[1] - (PH[1].g - 1.0) * dQ / (PH[1].g * (p + PH[1].pi));
	};
	auto Phi = [&](double p) { double a[2], dQ; eval(p, a, dQ); return a[0] + a[1] - 1.0; };
	double p; if (!root_above(Phi, plo_mix(), in.pguess, p)) return false;
	double a[2], dQ; eval(p, a, dQ);
	for (int k = 0; k < 2; k++) {
		o.a[k] = a[k]; o.m[k] = in.mtr[k]; if (!(o.a[k] > o.m[k] * PH[k].b)) return false;
		o.e[k] = eps_amp(k, o.a[k], o.m[k], p); o.p[k] = p; o.T[k] = T_amp(k, o.a[k], o.m[k], p);
	}
	o.pbar = p; return true;
}

// 给定 p, Y 由混合能量定 T（mu = theta = inf）
inline double T_from_energy(double p, double Y, double eC) {
	const double d = (1 - Y) * PH[0].cv * (p + PH[0].g * PH[0].pi) / (p + PH[0].pi) + Y * PH[1].cv * (p + PH[1].g * PH[1].pi) / (p + PH[1].pi);
	return (eC - (1 - Y) * PH[0].q - Y * PH[1].q) / d;
}

// B6.5  有限 nu（mu = theta = inf）：内层 R(Y) = rho Y - m2tr + tau nu (g2 - g1)(p, T(p, Y)) = 0，外层 rho sum Y_k v_k = 1
static bool closure_nu(const ClIn& in, ClOut& o) {
	const double rho = in.mtr[0] + in.mtr[1], eC = in.rhoeC / rho, tn = in.tau * nu_rate;
	auto Yof = [&](double p, double& Y, double& T) -> bool {
		auto R = [&](double y) { const double TT = T_from_energy(p, y, eC); if (!(TT > 0)) return (double)NAN; return rho * y - in.mtr[1] + tn * (g_pT(1, p, TT) - g_pT(0, p, TT)); };
		const double ya = YMIN, yb = 1.0 - YMIN; const double Ra = R(ya), Rb = R(yb);
		if (!std::isfinite(Ra) || !std::isfinite(Rb)) return false;
		if (Ra >= 0) Y = ya; else if (Rb <= 0) Y = yb; else if (!brent(R, ya, yb, Ra, Rb, 1e-16, Y)) return false;
		T = T_from_energy(p, Y, eC); return T > 0;
	};
	auto Phi = [&](double p) { double Y, T; if (!Yof(p, Y, T)) return (double)NAN; return rho * ((1 - Y) * v_pT(0, p, T) + Y * v_pT(1, p, T)) - 1.0; };
	double p; if (!root_above(Phi, plo_mix(), in.pguess, p)) return false;
	double Y, T; if (!Yof(p, Y, T)) return false;
	o.m[0] = rho * (1 - Y); o.m[1] = rho * Y;
	for (int k = 0; k < 2; k++) { o.a[k] = o.m[k] * v_pT(k, p, T); o.e[k] = o.m[k] * e_pT(k, p, T); o.p[k] = p; o.T[k] = T; }
	o.pbar = p; return true;
}

// B6.6  pTg（全部瞬时）：T = Tsat(p)，Y 由能量定，出界钳住（T 由能量定，连续）；p >= pcrit 时按 v 与 vcrit 的关系取单相
static bool pTg_state(double p, double rho, double eC, double& Y, double& T) {
	const double Ts = Tsat(p);
	if (std::isfinite(Ts)) {
		const double el = e_pT(0, p, Ts), ev = e_pT(1, p, Ts);
		Y = (eC - el) / (ev - el);
		if (Y > YMIN && Y < 1 - YMIN) { T = Ts; return true; }
		Y = std::min(std::max(Y, YMIN), 1 - YMIN);
	} else {
		if (!(p >= SAT.pcrit)) return false;
		Y = (1.0 / rho < SAT.vcrit) ? YMIN : 1 - YMIN;
	}
	T = T_from_energy(p, Y, eC); return T > 0;
}
static bool closure_pTg(const ClIn& in, ClOut& o) {
	const double rho = in.mtr[0] + in.mtr[1], eC = in.rhoeC / rho;
	auto Phi = [&](double p) { double Y, T; if (!pTg_state(p, rho, eC, Y, T)) return (double)NAN; return rho * ((1 - Y) * v_pT(0, p, T) + Y * v_pT(1, p, T)) - 1.0; };
	double p; if (!root_above(Phi, plo_mix(), in.pguess, p)) return false;
	double Y, T; if (!pTg_state(p, rho, eC, Y, T)) return false;
	o.m[0] = rho * (1 - Y); o.m[1] = rho * Y;
	for (int k = 0; k < 2; k++) { o.a[k] = o.m[k] * v_pT(k, p, T); o.e[k] = o.m[k] * e_pT(k, p, T); o.p[k] = p; o.T[k] = T; }
	o.pbar = p; return true;
}

static bool closure(const ClIn& in, ClOut& o) {
	for (int k = 0; k < 2; k++) if (!(in.atr[k] > 0) || !(in.mtr[k] > 0) || !std::isfinite(in.etr[k])) return false;
	if (!std::isfinite(in.rhoeC) || !std::isfinite(in.s)) return false;
	bool ok = false;
	switch (in.type) {
		case CL_P: ok = ISEN ? closure_p_isen(in, o) : closure_p(in, o); break;
		case CL_MU: ok = closure_mu(in, o); break;
		case CL_PT: ok = closure_pT(in, o); break;
		case CL_THETA: ok = closure_theta(in, o); break;
		case CL_NU: ok = closure_nu(in, o); break;
		case CL_PTG: ok = closure_pTg(in, o); break;
	}
	return ok && std::isfinite(o.pbar);
}

// 松弛档位 -> 每格的闭合类型（B6.7）；输入是单位体积的“内容”与比内能
//   界面判据 eps_I <= alpha1 <= 1 - eps_I；化学松弛的判据 T > Tsat(p)（P-S (10)，Pelanti 8.2.1）
static int closure_type_of(double a1, double m1, double m2, double e1, double e2, double rhoe) {
	const int mech = std::isinf(mu_rate) ? CL_P : CL_MU;
	if (RELAX_LEVEL == 0) return mech;
	if (!(a1 >= EPS_I && a1 <= 1.0 - EPS_I)) return mech;
	const int chem = std::isinf(nu_rate) ? CL_PTG : CL_NU;
	ClIn in{}; in.atr[0] = a1; in.atr[1] = 1 - a1; in.mtr[0] = m1; in.mtr[1] = m2; in.etr[0] = e1; in.etr[1] = e2; in.s = 0; in.rhoeC = rhoe; in.tau = 0; in.type = CL_PT;
	in.pguess = std::max(p_mix(a1, m1, m2, rhoe), 1.0);
	if (RELAX_LEVEL == 3) {   // p-pTg：力学平衡态下任一相 T_k > Tsat(p) 时做 pTg
		in.type = CL_P; in.nodelta = true; ClOut o; if (!closure_p(in, o)) return mech;
		const double Ts = Tsat(o.p[0]); if (std::isfinite(Ts) && (o.T[0] > Ts || o.T[1] > Ts)) return chem; return mech;
	}
	const int therm = std::isinf(theta_rate) ? CL_PT : CL_THETA;
	if (RELAX_LEVEL == 2) {
		ClOut o; if (closure_pT(in, o)) { const double Ts = Tsat(o.p[0]); if (std::isfinite(Ts) && o.T[0] > Ts) return chem; }
	}
	return therm;
}

// B6.7'：激活与投影（每步一次，P-S 的算子顺序）。瞬时松弛（速率 = inf）在激活时是状态的跳跃；若在级内激活，
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

// ==================================================================================================== 边界条件（2DEuler_TV 的 bound_cond_u）
static void bound_cond_u(Grid& U) {
	const int g = ghost_cells, nx = U.nx, ny = U.ny;
	for (int j = 0; j < ny + 2 * g; j++) for (int q = 0; q < g; q++) {
		int iL = g - 1 - q, iR = nx + g + q, sL, sR;
		if (x_bc == x_Periodic) { sL = nx + g - 1 - q; sR = g + q; }
		else if (x_bc == x_Reflective) { sL = g + q; sR = nx + g - 1 - q; }
		else { sL = g; sR = nx + g - 1; }
		for (int k = 0; k < NUM_VARS; k++) { U(iL, j, k) = U(sL, j, k); U(iR, j, k) = U(sR, j, k); }
		if (x_bc == x_Reflective) { U(iL, j, IMU) = -U(iL, j, IMU); U(iR, j, IMU) = -U(iR, j, IMU); }
	}
	if (ny > 1 || y_bc == y_Periodic) for (int i = 0; i < nx + 2 * g; i++) for (int q = 0; q < g; q++) {
		int jL = g - 1 - q, jR = ny + g + q, sL, sR;
		if (y_bc == y_Periodic) { sL = ny + g - 1 - q; sR = g + q; }
		else if (y_bc == y_Reflective) { sL = g + q; sR = ny + g - 1 - q; }
		else { sL = g; sR = ny + g - 1; }
		for (int k = 0; k < NUM_VARS; k++) { U(i, jL, k) = U(i, sL, k); U(i, jR, k) = U(i, sR, k); }
		if (y_bc == y_Reflective) { U(i, jL, IMV) = -U(i, jL, IMV); U(i, jR, IMV) = -U(i, jR, IMV); }
	}
	else for (int i = 0; i < nx + 2 * g; i++) for (int q = 0; q < g; q++) for (int k = 0; k < NUM_VARS; k++) { U(i, g - 1 - q, k) = U(i, g, k); U(i, g + 1 + q, k) = U(i, g, k); }
}
static void bound_cond_array(Array& A, int sign_x = 1, int sign_y = 1) {
	const int g = ghost_cells, nx = A.nx, ny = A.ny;
	for (int j = 0; j < ny + 2 * g; j++) for (int q = 0; q < g; q++) {
		int iL = g - 1 - q, iR = nx + g + q, sL, sR; double fs = 1.0;
		if (x_bc == x_Periodic) { sL = nx + g - 1 - q; sR = g + q; } else if (x_bc == x_Reflective) { sL = g + q; sR = nx + g - 1 - q; fs = sign_x; } else { sL = g; sR = nx + g - 1; }
		A(iL, j) = fs * A(sL, j); A(iR, j) = fs * A(sR, j);
	}
	for (int i = 0; i < nx + 2 * g; i++) for (int q = 0; q < g; q++) {
		int jL = g - 1 - q, jR = ny + g + q, sL, sR; double fs = 1.0;
		if (ny == 1 && y_bc != y_Periodic) { sL = g; sR = g; }
		else if (y_bc == y_Periodic) { sL = ny + g - 1 - q; sR = g + q; } else if (y_bc == y_Reflective) { sL = g + q; sR = ny + g - 1 - q; fs = sign_y; } else { sL = g; sR = ny + g - 1; }
		A(i, jL) = fs * A(i, sL); A(i, jR) = fs * A(i, sR);
	}
}

// ==================================================================================================== 单元工具
struct Prim { double a1, r1, r2, u, v, p1, p2, rho, cf; };
static bool cell_prim(const Grid& U, int i, int j, Prim& w) {
	const double a1 = U(i, j, IA1), a2 = 1 - a1, m1 = U(i, j, IM1), m2 = U(i, j, IM2);
	if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return false;
	w.a1 = a1; w.rho = m1 + m2; w.r1 = m1 / a1; w.r2 = m2 / a2; w.u = U(i, j, IMU) / w.rho; w.v = U(i, j, IMV) / w.rho;
	w.p1 = p_ame(0, a1, m1, U(i, j, IE1)); w.p2 = p_ame(1, a2, m2, U(i, j, IE2));
	if (!(w.p1 + PH[0].pi > 0 && w.p2 + PH[1].pi > 0)) return false;
	w.cf = std::sqrt((m1 * c2_rp(0, w.r1, w.p1) + m2 * c2_rp(1, w.r2, w.p2)) / w.rho);
	return true;
}
static void prim2cell(const Prim& w, double* Q) {
	const double a2 = 1 - w.a1, m1 = w.a1 * w.r1, m2 = a2 * w.r2, rho = m1 + m2;
	Q[IA1] = w.a1; Q[IM1] = m1; Q[IM2] = m2; Q[IMU] = rho * w.u; Q[IMV] = rho * w.v;
	Q[IE1] = eps_amp(0, w.a1, m1, w.p1); Q[IE2] = eps_amp(1, a2, m2, w.p2);
	Q[IEN] = Q[IE1] + Q[IE2] + 0.5 * rho * (w.u * w.u + w.v * w.v);
}
inline double pbar_cell(const Grid& U, int i, int j) {
	const double a1 = U(i, j, IA1);
	return a1 * p_ame(0, a1, U(i, j, IM1), U(i, j, IE1)) + (1 - a1) * p_ame(1, 1 - a1, U(i, j, IM2), U(i, j, IE2));
}

// ==================================================================================================== MODE 0：显式参考（HLLC + Godunov + 松弛）
struct HFace { double F[5], ustar, ua1, ue[2]; };   // F: m1, m2, normal mom, tangential mom, E; ustar; (u alpha1)*; (eps_k u)*
// 相 Hugoniot (II.8) 的 NASG 形式: e* - e + 1/2 (p* + p)(v* - v) = 0 对 p* 线性
static double hugoniot_p(int k, double p, double v, double vs) {
	const Phase& P = PH[k];
	const double num = (p + P.g * P.pi) * (v - P.b) - P.g * P.pi * (vs - P.b) - 0.5 * (P.g - 1.0) * p * (vs - v);
	const double den = (vs - P.b) + 0.5 * (P.g - 1.0) * (vs - v);
	return num / den;
}
// dir = 0: x 面 (法向 u)；dir = 1: y 面 (法向 v)
static HFace hllc(const Prim& L, const Prim& R, int dir) {
	auto un = [&](const Prim& w) { return dir == 0 ? w.u : w.v; };
	auto ut = [&](const Prim& w) { return dir == 0 ? w.v : w.u; };
	const double uL = un(L), uR = un(R), pL = L.a1 * L.p1 + (1 - L.a1) * L.p2, pR = R.a1 * R.p1 + (1 - R.a1) * R.p2;
	const double SL = std::min(uL - L.cf, uR - R.cf), SR = std::max(uL + L.cf, uR + R.cf);
	const double SM = (pR - pL + L.rho * uL * (SL - uL) - R.rho * uR * (SR - uR)) / (L.rho * (SL - uL) - R.rho * (SR - uR));
	auto EK = [&](const Prim& w) { const double a2 = 1 - w.a1; return eps_amp(0, w.a1, w.a1 * w.r1, w.p1) + eps_amp(1, a2, a2 * w.r2, w.p2) + 0.5 * w.rho * (w.u * w.u + w.v * w.v); };
	auto flux = [&](const Prim& w, double* F) { const double u = un(w), p = w.a1 * w.p1 + (1 - w.a1) * w.p2;
		F[0] = w.a1 * w.r1 * u; F[1] = (1 - w.a1) * w.r2 * u; F[2] = w.rho * u * u + p; F[3] = w.rho * u * ut(w); F[4] = (EK(w) + p) * u; };
	HFace f;
	const bool left = SM >= 0; const Prim& K = left ? L : R; const double uK = left ? uL : uR, SK = left ? SL : SR, pK = left ? pL : pR;
	double FK[5]; flux(K, FK);
	if ((left && SL >= 0) || (!left && SR <= 0)) {
		for (int q = 0; q < 5; q++) f.F[q] = FK[q];
		f.ustar = uK; f.ua1 = uK * K.a1;
		f.ue[0] = eps_amp(0, K.a1, K.a1 * K.r1, K.p1) * uK; f.ue[1] = eps_amp(1, 1 - K.a1, (1 - K.a1) * K.r2, K.p2) * uK; return f;
	}
	const double fac = (SK - uK) / (SK - SM), EKv = EK(K);
	const double rs = K.rho * fac;
	const double Es = fac * (EKv + (SM - uK) * (K.rho * SM + pK / (SK - uK)));
	const double UK[5] = {K.a1 * K.r1, (1 - K.a1) * K.r2, K.rho * uK, K.rho * ut(K), EKv};
	const double Us[5] = {UK[0] * fac, UK[1] * fac, rs * SM, rs * ut(K), Es};
	for (int q = 0; q < 5; q++) f.F[q] = FK[q] + SK * (Us[q] - UK[q]);
	f.ustar = SM; f.ua1 = SM * K.a1;
	const double ak[2] = {K.a1, 1 - K.a1}, rk[2] = {K.r1, K.r2}, pk[2] = {K.p1, K.p2};
	for (int k = 0; k < 2; k++) {
		const double v0 = 1.0 / rk[k], vs = 1.0 / (rk[k] * fac);
		double ps = hugoniot_p(k, pk[k], v0, vs); if (!(ps + PH[k].pi > 0) || !std::isfinite(ps)) ps = pk[k];
		f.ue[k] = eps_amp(k, ak[k], ak[k] * rk[k] * fac, ps) * SM;
	}
	return f;
}

inline double pbar_cell_q(const double* Q) {
	const double a1 = Q[IA1];
	return a1 * p_ame(0, a1, Q[IM1], Q[IE1]) + (1 - a1) * p_ame(1, 1 - a1, Q[IM2], Q[IE2]);
}

// 显式步的局部松弛：力学（SPB 原文：Gibbs 关系 + 重初始化）然后热/化学（P-S / Pelanti，用统一闭合 s = 0）
static bool relax_cell_explicit(double* Q, double dtloc) {
	const double rho = Q[IM1] + Q[IM2], rhoe = Q[IEN] - 0.5 * (Q[IMU] * Q[IMU] + Q[IMV] * Q[IMV]) / rho;
	ClIn in{}; in.atr[0] = Q[IA1]; in.atr[1] = 1 - Q[IA1]; in.mtr[0] = Q[IM1]; in.mtr[1] = Q[IM2]; in.etr[0] = Q[IE1]; in.etr[1] = Q[IE2];
	in.s = 0; in.rhoeC = rhoe; in.tau = dtloc; in.pguess = std::max(pbar_cell_q(Q), 1.0);
	ClOut o;
	const int mech = std::isinf(mu_rate) ? CL_P : CL_MU;
	in.type = mech; in.nodelta = (SPB_EXACT && mech == CL_P);
	if (!closure(in, o)) { st_cl_fail++; return false; }
	double a1 = o.a[0], m1 = o.m[0], m2 = o.m[1], e1 = o.e[0], e2 = o.e[1];
	if (in.nodelta) {   // SPB 3.4 重初始化：alpha 固定，用守恒总能由 (A2.3) 重算 p，再重置相内能
		const double p = p_mix(a1, m1, m2, rhoe);
		if (!(p + PH[0].pi > 0 && p + PH[1].pi > 0)) { st_cl_fail++; return false; }
		e1 = eps_amp(0, a1, m1, p); e2 = eps_amp(1, 1 - a1, m2, p);
	}
	const int ty = closure_type_of(a1, m1, m2, e1, e2, rhoe);
	if (ty != mech) {
		ClIn t2{}; t2.atr[0] = a1; t2.atr[1] = 1 - a1; t2.mtr[0] = m1; t2.mtr[1] = m2; t2.etr[0] = e1; t2.etr[1] = e2;
		t2.s = 0; t2.rhoeC = rhoe; t2.tau = dtloc; t2.type = ty; t2.pguess = std::max(p_mix(a1, m1, m2, rhoe), 1.0);
		ClOut o2; if (!closure(t2, o2)) { st_cl_fail++; } else { a1 = o2.a[0]; m1 = o2.m[0]; m2 = o2.m[1]; e1 = o2.e[0]; e2 = o2.e[1]; }
	}
	Q[IA1] = a1; Q[IM1] = m1; Q[IM2] = m2; Q[IE1] = e1; Q[IE2] = e2;
	return true;
}

static void reconstruct_prim_minmod(const Grid& U, std::vector<Prim>& WL, std::vector<Prim>& WR, int dir, bool& ok);

static bool explicit_step(Grid& U, double dtl) {
	const int g = ghost_cells, nx = U.nx, ny = U.ny, NXT = nx + 2 * g, NYT = ny + 2 * g;
	bound_cond_u(U);
	std::vector<Prim> W((size_t)NXT * NYT);
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) if (!cell_prim(U, i, j, W[i + (size_t)j * NXT])) return false;
	std::vector<Prim> XL, XR, YL, YR; bool ok = true;
	reconstruct_prim_minmod(U, XL, XR, 0, ok); if (ny > 1) reconstruct_prim_minmod(U, YL, YR, 1, ok);
	if (!ok) return false;
	Grid V = U; const double rx = dtl / dx, ry = dtl / dy;
	std::vector<HFace> FX((size_t)NXT * NYT), FY;
	for (int j = g; j < ny + g; j++) for (int i = g - 1; i < nx + g; i++) FX[i + (size_t)j * NXT] = hllc(XR[i + (size_t)j * NXT], XL[i + 1 + (size_t)j * NXT], 0);
	if (ny > 1) { FY.resize((size_t)NXT * NYT); for (int j = g - 1; j < ny + g; j++) for (int i = g; i < nx + g; i++) FY[i + (size_t)j * NXT] = hllc(YR[i + (size_t)j * NXT], YL[i + (size_t)(j + 1) * NXT], 1); }
	for (int j = g; j < ny + g; j++) for (int i = g; i < nx + g; i++) {
		const Prim& w = W[i + (size_t)j * NXT];
		const HFace& fp = FX[i + (size_t)j * NXT], & fm = FX[i - 1 + (size_t)j * NXT];
		double dU[NUM_VARS] = {0};
		const double ak[2] = {w.a1, 1 - w.a1}, pk[2] = {w.p1, w.p2};
		dU[IM1] = -rx * (fp.F[0] - fm.F[0]); dU[IM2] = -rx * (fp.F[1] - fm.F[1]); dU[IMU] = -rx * (fp.F[2] - fm.F[2]); dU[IMV] = -rx * (fp.F[3] - fm.F[3]); dU[IEN] = -rx * (fp.F[4] - fm.F[4]);
		dU[IA1] = -rx * ((fp.ua1 - fm.ua1) - w.a1 * (fp.ustar - fm.ustar));
		for (int k = 0; k < 2; k++) dU[IE1 + k] = -rx * ((fp.ue[k] - fm.ue[k]) + ak[k] * pk[k] * (fp.ustar - fm.ustar));
		if (ny > 1) {
			const HFace& gp = FY[i + (size_t)j * NXT], & gm = FY[i + (size_t)(j - 1) * NXT];
			dU[IM1] -= ry * (gp.F[0] - gm.F[0]); dU[IM2] -= ry * (gp.F[1] - gm.F[1]); dU[IMV] -= ry * (gp.F[2] - gm.F[2]); dU[IMU] -= ry * (gp.F[3] - gm.F[3]); dU[IEN] -= ry * (gp.F[4] - gm.F[4]);
			dU[IA1] -= ry * ((gp.ua1 - gm.ua1) - w.a1 * (gp.ustar - gm.ustar));
			for (int k = 0; k < 2; k++) dU[IE1 + k] -= ry * ((gp.ue[k] - gm.ue[k]) + ak[k] * pk[k] * (gp.ustar - gm.ustar));
		}
		for (int k = 0; k < NUM_VARS; k++) V(i, j, k) = U(i, j, k) + dU[k];
		double Q[NUM_VARS]; for (int k = 0; k < NUM_VARS; k++) Q[k] = V(i, j, k);
		if (!(Q[IA1] > 0 && Q[IA1] < 1 && Q[IM1] > 0 && Q[IM2] > 0)) return false;
		if (!relax_cell_explicit(Q, dtl)) return false;
		for (int k = 0; k < NUM_VARS; k++) V(i, j, k) = Q[k];
	}
	U = V; bound_cond_u(U); return true;
}

// 原始变量的 minmod 重构（MODE 0；SPB 附录 C 也用原始变量）
static void reconstruct_prim_minmod(const Grid& U, std::vector<Prim>& WL, std::vector<Prim>& WR, int dir, bool& ok) {
	const int g = ghost_cells, nx = U.nx, ny = U.ny, NXT = nx + 2 * g, NYT = ny + 2 * g;
	std::vector<Prim> W((size_t)NXT * NYT);
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) if (!cell_prim(U, i, j, W[i + (size_t)j * NXT])) { ok = false; return; }
	WL = W; WR = W;
	if (SpaceOrder < 2) return;
	auto mm = [](double a, double b) { return a * b > 0 ? (std::fabs(a) < std::fabs(b) ? a : b) : 0.0; };
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
		int im, ip, jm = j, jp = j, imx = i, ipx = i;
		if (dir == 0) { im = i - 1; ip = i + 1; if (im < 0 || ip >= NXT) continue; imx = im; ipx = ip; } else { jm = j - 1; jp = j + 1; if (jm < 0 || jp >= NYT) continue; }
		const Prim& c = W[i + (size_t)j * NXT], & a = W[imx + (size_t)jm * NXT], & b = W[ipx + (size_t)jp * NXT];
		const double cv[7] = {c.a1, c.r1, c.r2, c.u, c.v, c.p1, c.p2}, av[7] = {a.a1, a.r1, a.r2, a.u, a.v, a.p1, a.p2}, bv[7] = {b.a1, b.r1, b.r2, b.u, b.v, b.p1, b.p2};
		double l[7], r[7];
		for (int q = 0; q < 7; q++) { const double sl = mm(cv[q] - av[q], bv[q] - cv[q]); l[q] = cv[q] - 0.5 * sl; r[q] = cv[q] + 0.5 * sl; }
		Prim PL = c, PR = c;
		PL.a1 = l[0]; PL.r1 = l[1]; PL.r2 = l[2]; PL.u = l[3]; PL.v = l[4]; PL.p1 = l[5]; PL.p2 = l[6];
		PR.a1 = r[0]; PR.r1 = r[1]; PR.r2 = r[2]; PR.u = r[3]; PR.v = r[4]; PR.p1 = r[5]; PR.p2 = r[6];
		for (Prim* P : {&PL, &PR}) { const double m1 = P->a1 * P->r1, m2 = (1 - P->a1) * P->r2; P->rho = m1 + m2; P->cf = std::sqrt((m1 * c2_rp(0, P->r1, P->p1) + m2 * c2_rp(1, P->r2, P->p2)) / P->rho); }
		WL[i + (size_t)j * NXT] = PL; WR[i + (size_t)j * NXT] = PR;
	}
}

// ==================================================================================================== MODE 1：TV 分裂 IMEX
// ---- B3 显式算子：共用限制因子的 minmod（在守恒量上），BB 的 Rusanov（中心部分用格均值，耗散用重构值）
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
inline void face_weights(FaceVals& f, double ube) {
	if (INCDONOR) { f.wL = std::max(f.uf, 0.0) - std::max(ube, 0.0); f.wR = std::min(f.uf, 0.0) - std::min(ube, 0.0); }
	else { f.wL = f.du > 0 ? f.du : 0.0; f.wR = f.du > 0 ? 0.0 : f.du; }
}
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
	if (which == 0 && perW) PW += dn; if (which == 1 && perE) PE += dn; if (which == 2 && perS) PS += dn; if (which == 3 && perN) PN += dn;
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
	// BB 的“中心部分用格均值”形式 (CENTRAL_CELLAVG = 1) 对线性输运给出 u X_i + 1/4 u (sig_i + sig_{i+1})，含下游斜率，不是 TVD。
	auto flux_face = [&](int i0, int j0, int i1, int j1, const std::vector<double>& sig, int dir, double a, double* F, double& ubar) {
		const size_t k0 = ((size_t)i0 + (size_t)j0 * NXT) * NUM_VARS, k1 = ((size_t)i1 + (size_t)j1 * NXT) * NUM_VARS;
		double XL[NUM_VARS], XR[NUM_VARS];
		for (int k = 0; k < NUM_VARS; k++) { XL[k] = Qe(i0, j0, k) + 0.5 * sig[k0 + k]; XR[k] = Qe(i1, j1, k) - 0.5 * sig[k1 + k]; }
		const int kn = dir == 0 ? IMU : IMV;
		double uL, uR;
		if (CENTRAL_CELLAVG) { uL = Qe(i0, j0, kn) / (Qe(i0, j0, IM1) + Qe(i0, j0, IM2)); uR = Qe(i1, j1, kn) / (Qe(i1, j1, IM1) + Qe(i1, j1, IM2)); }
		else { uL = XL[kn] / (XL[IM1] + XL[IM2]); uR = XR[kn] / (XR[IM1] + XR[IM2]); }
		for (int k = 0; k < NUM_VARS; k++) {
			const double cL = CENTRAL_CELLAVG ? Qe(i0, j0, k) : XL[k], cR = CENTRAL_CELLAVG ? Qe(i1, j1, k) : XR[k];
			F[k] = 0.5 * (cL * uL + cR * uR) - 0.5 * a * (XR[k] - XL[k]);
		}
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
	// B4.1：面两侧各有质量 rho_k l_k 参与，l_k = min(c_k tau, dx/2)（声学影响域，不超过半格），阻抗 Zhat_k = rho_k l_k / tau。
	// 两侧动量平衡 Zhat_L (u_f - u*_L) = P_L - p_f，Zhat_R (u_f - u*_R) = p_f - P_R 给出
	//   u_f = (Zhat_L u*_L + Zhat_R u*_R - (P_R - P_L)) / (Zhat_L + Zhat_R)，p_f = (Zhat_R P_L + Zhat_L P_R - Zhat_L Zhat_R (u*_R - u*_L)) / (Zhat_L + Zhat_R)。
	// tau -> 0：声学 Riemann 解（Godunov）；tau 大：半格动量平衡（压力项 = tau dP/(rhobar dx)，速度跳跃项的粘性 ~ dx^2/(4 tau) = O(物质速度 dx)）。
	// c 用冻结声速（各松弛声速的上界，满足次特征条件）；X* 不可容许时取 c = inf（即 l = dx/2）。
	auto cfz = [&](int i, int j) -> double {
		const double a1 = S.Xs(i, j, IA1), a2 = S.a2s(i, j), m1 = S.Xs(i, j, IM1), m2 = S.Xs(i, j, IM2);
		if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return INF;
		const double p1 = p_ame(0, a1, m1, S.Xs(i, j, IE1)), p2 = p_ame(1, a2, m2, S.Xs(i, j, IE2));
		if (!(p1 + PH[0].pi > 0 && p2 + PH[1].pi > 0)) return INF;
		const double c21 = c2_rp(0, m1 / a1, p1), c22 = c2_rp(1, m2 / a2, p2), sa = a1 + a2;
		if (CSPEED == 1) { const double k = (a1 / sa) / (m1 / a1 * c21) + (a2 / sa) / (m2 / a2 * c22); return k > 0 ? std::sqrt(sa / ((m1 + m2) * k)) : INF; }   // Wood: 1/(rho c^2) = sum alpha_k/(rho_k c_k^2)
		const double c2 = (m1 * c21 + m2 * c22) / (m1 + m2);
		return c2 > 0 ? std::sqrt(c2) : INF; };
	std::vector<double> cz((size_t)NXT * NYT, INF);
	if (ACOUSTIC) for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) cz[S.id(i, j)] = cfz(i, j);
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
		if (SpaceOrder < 2 || ACORD < 2) return 1.0;
		const double D = xR - xL; if (D == 0.0) return 1.0;
		return 1.0 - 0.5 * (minmod(xL - xm, D) + minmod(D, xp - xR)) / D; };
	// B4.1：u_f = ubar* - tau dP/(rhobar h) - phi psi_P dP/(Z_L + Z_R)，p_f = (rho_R P_L + rho_L P_R)/(rho_L + rho_R) - phi psi_u Z_L Z_R/(Z_L + Z_R) du*，
	// phi = (1 - c tau/h)^+：后向 Euler 的声学耗散 c^2 tau/2 不足声学 Godunov 的 c h/2 时补足差额（c tau/h >= 1 时为 0，低马赫/大步长不受影响）
	auto set_face = [&](int iL, int jL, int iR, int jR, double h, bool ydir, double& ubs, double& kp, double& wL, double& pv) {
		const double rL = rhoS(iL, jL), rR = rhoS(iR, jR);
		const double uL = ydir ? S.vst(iL, jL) : S.ust(iL, jL), uR = ydir ? S.vst(iR, jR) : S.ust(iR, jR);
		ubs = 0.5 * (uL + uR); kp = tau / h / (0.5 * (rL + rR)); wL = rR / (rL + rR); pv = 0.0;
		if (!ACOUSTIC) return;
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

// ==================================================================================================== 时间步长（caldt）
static double max_speed(const Grid& U, bool acoustic, double& ux, double& uy) {
	ux = uy = 0; const int g = ghost_cells;
	for (int j = g; j < U.ny + g; j++) for (int i = g; i < U.nx + g; i++) {
		Prim w; if (!cell_prim(U, i, j, w)) continue;
		const double c = acoustic ? w.cf : 0.0;
		ux = std::max(ux, std::fabs(w.u) + c); uy = std::max(uy, std::fabs(w.v) + c);
	}
	return 0;
}
static double caldt(const Grid& U, bool acoustic) {
	double ux, uy; max_speed(U, acoustic, ux, uy);
	double rate = ux / dx + (U.ny > 1 ? uy / dy : 0.0);
	return rate > 0 ? 1.0 / rate : INF;
}

// ==================================================================================================== 算例
struct Case { std::string name; };
static double ALPHA_EPS_INIT = 1e-8;
static void set_SG(int k, double g, double pi, double q = 0, double qp = 0, double cv = 1000.0, double b = 0) { PH[k] = {g, pi, q, qp, cv, b}; }
static double rho_pT(int k, double p, double T) { return 1.0 / v_pT(k, p, T); }

static void init_case(Grid& U) {
	const int g = ghost_cells;
	auto setcell = [&](int i, int j, double a1, double r1, double r2, double u, double v, double p1, double p2) {
		Prim w; w.a1 = a1; w.r1 = r1; w.r2 = r2; w.u = u; w.v = v; w.p1 = p1; w.p2 = p2; double Q[NUM_VARS]; prim2cell(w, Q);
		for (int k = 0; k < NUM_VARS; k++) U(i, j, k) = Q[k]; };
	for (int j = g; j < U.ny + g; j++) for (int i = g; i < U.nx + g; i++) {
		const double x = x_0 + (i - g + 0.5) * dx, y = y_0 + (j - g + 0.5) * dy; (void)y;
		switch (TESTCASE) {
		case 0: { // SPB 4.5 空化管：水/空气 alpha_air = 1e-2，+-100 m/s（空气密度论文未给，取 1）
			const double u = x < 0.5 ? -100.0 : 100.0; setcell(i, j, 1 - 1e-2, 1000.0, 1.0, u, 0, 1e5, 1e5); break; }
		case 1: case 2: { // P-S 2014 水空化管 +-2 / +-500 m/s, T = 354.728 K
			const double T0 = 354.728, U0 = TESTCASE == 1 ? 2.0 : 500.0, u = x < 0.5 ? -U0 : U0;
			setcell(i, j, 1 - 1e-2, rho_pT(0, 1e5, T0), rho_pT(1, 1e5, T0), u, 0, 1e5, 1e5); break; }
		case 3: case 31: { // 十二烷液-汽激波管 (P-S 5.2 / Pelanti 8.2.2)
			const bool L = x < 0.75; setcell(i, j, L ? 1 - 1e-8 : 1e-8, 500.0, 2.0, 0, 0, L ? 1e8 : 1e5, L ? 1e8 : 1e5); break; }
		case 4: { // SPB 4.1 界面平移
			const bool L = x < 0.5; setcell(i, j, L ? 1 - 1e-8 : 1e-8, 1000.0, 10.0, 100.0, 0, 1e5, 1e5); break; }
		case 5: case 6: { // SPB 4.3.1 / 4.3.2 水-气激波管
			const bool L = x < 0.75; const double pL = TESTCASE == 5 ? 1e9 : 1e12, ra = TESTCASE == 5 ? 1.0 : 10.0;
			setcell(i, j, L ? 1 - 1e-6 : 1e-6, 1000.0, ra, 0, 0, L ? pL : 1e5, L ? pL : 1e5); break; }
		case 12: case 13: { // SPB 5.1 / 5.2 环氧-尖晶石
			const double pL = TESTCASE == 12 ? 1e10 : 2e11; const double p = x < 0.6 ? pL : 1e5;
			setcell(i, j, 0.5954, 1185.0, 3622.0, 0, 0, p, p); break; }
		case 20: case 21: { // 光滑解（时间收敛阶）
			const double gs = std::exp(-(x - 0.5) * (x - 0.5) / 0.01), u = 50.0 * gs, a1 = 0.5 + 0.2 * std::sin(2 * M_PI * x), p = 1e5 * (1 + 0.2 * gs);
			setcell(i, j, a1, TESTCASE == 20 ? 2.0 : 1000.0, 1.0, u, 0, p, p); break; }
		case 30: { // Pelanti 2022 8.2.1：NASG 水空化管 +-2 m/s，T = 353 K，p = 1e5，alpha_v = 1e-2
			const double T0 = 353.0, u = x < 0.5 ? -2.0 : 2.0; setcell(i, j, 1 - 1e-2, rho_pT(0, 1e5, T0), rho_pT(1, 1e5, T0), u, 0, 1e5, 1e5); break; }
		case 40: { // 0D 松弛：u = 0，均匀非平衡态（周期）
			const double T1 = 350.0, T2 = 380.0, p1 = 1.2e5, p2 = 0.8e5;
			setcell(i, j, 1 - 0.2, rho_pT(0, p1, T1), rho_pT(1, p2, T2), 0, 0, p1, p2); break; }
		}
	}
	bound_cond_u(U);
}
static void setup_case() {
	x_0 = 0; x_1 = 1; y_0 = 0; y_1 = 1; x_bc = x_Transmissive; y_bc = y_Transmissive;
	switch (TESTCASE) {
	case 0: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); time_stop = 1.85e-3; break;
	case 1: case 2:
		set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0);
		time_stop = TESTCASE == 1 ? 3.2e-3 : 0.58e-3; break;
	case 3: set_SG(0, 2.35, 4e8, -775.269e3, 0.0, 1077.7); set_SG(1, 1.025, 0.0, -237.547e3, -24.4e3, 1956.45); time_stop = 473e-6; if (EPS_I == 0) EPS_I = 1e-4; break;
	case 4: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); time_stop = 2.79e-3; break;
	case 5: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); time_stop = 240e-6; break;
	case 6: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); time_stop = 8.3e-6; break;
	case 12: case 13: set_SG(0, 2.43, 5.3e9); set_SG(1, 1.62, 141e9); time_stop = TESTCASE == 12 ? 80e-6 : 34e-6; break;
	case 20: set_SG(0, 1.4, 0); set_SG(1, 1.4, 0); time_stop = 5e-4; break;
	case 21: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); time_stop = 1e-4; break;
	case 30: PH[0] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4}; PH[1] = {1.467, 0.0, 2077616.0, 14317.0, 955.0, 0.0}; time_stop = 3e-3; break;
	case 31: PH[0] = {1.206, 1681e5, -996054.0, 0.0, 2532.0, 7.51e-4}; PH[1] = {1.021, 0.0, -384592.0, -4301.0, 2274.0, 0.0}; time_stop = 473e-6; if (EPS_I == 0) EPS_I = 1e-4; break;
	case 40: set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0); time_stop = 1e-4; x_bc = x_Periodic; break;
	default: fprintf(stderr, "unknown case %d\n", TESTCASE); exit(1);
	}
}

// ==================================================================================================== 输出
static void output(const Grid& U, const char* fn) {
	FILE* f = fopen(fn, "w"); if (!f) return;
	const int g = ghost_cells;
	fprintf(f, "# x y p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2\n");
	for (int j = g; j < U.ny + g; j++) for (int i = g; i < U.nx + g; i++) {
		Prim w; const bool ok = cell_prim(U, i, j, w);
		const double a1 = U(i, j, IA1), m1 = U(i, j, IM1), m2 = U(i, j, IM2), rho = m1 + m2;
		const double p = ok ? a1 * w.p1 + (1 - a1) * w.p2 : NAN;
		double cw = NAN, T1 = NAN, T2 = NAN, dg = NAN;
		if (ok) {
			cw = std::sqrt(1.0 / (rho * (a1 / (w.r1 * c2_rp(0, w.r1, w.p1)) + (1 - a1) / (w.r2 * c2_rp(1, w.r2, w.p2)))));
			T1 = T_amp(0, a1, m1, w.p1); T2 = T_amp(1, 1 - a1, m2, w.p2);
			if (T1 > 0 && T2 > 0 && PH[0].cv > 0) dg = g_pT(0, w.p1, T1) - g_pT(1, w.p2, T2);
		}
		fprintf(f, "%.10e %.10e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.8e %.8e %.10e %.10e %.12e %.12e %.8e\n",
			x_0 + (i - g + 0.5) * dx, y_0 + (j - g + 0.5) * dy, p, ok ? w.u : NAN, ok ? w.v : NAN, rho, 1 - a1, m2 / rho, ok ? w.r1 : NAN, ok ? w.r2 : NAN, ok ? w.cf : NAN, cw, T1, T2, ok ? w.p1 : NAN, ok ? w.p2 : NAN, dg);
	}
	fclose(f);
}

// ==================================================================================================== main（2DEuler_TV 的时间循环）
int main(int argc, char** argv) {
	std::string out = "out.dat"; double tend_override = -1;
	for (int a = 1; a < argc; a++) {
		char key[64]; char val[256]; if (sscanf(argv[a], "%63[^=]=%255s", key, val) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); return 1; }
		std::string k(key); const double d = atof(val);
		if (k == "case") TESTCASE = (int)d; else if (k == "N" || k == "nx") nxCells = (int)d; else if (k == "ny") nyCells = (int)d;
		else if (k == "mode") MODE = (int)d; else if (k == "cfl") CFL = d; else if (k == "cflac") CFL_ac = d; else if (k == "order") SpaceOrder = (int)d; else if (k == "rk") RKs = (int)d;
		else if (k == "tend") tend_override = d; else if (k == "out") out = val; else if (k == "ramp") RAMP = (int)d;
		else if (k == "relax") RELAX_LEVEL = (int)d; else if (k == "epsI") EPS_I = d; else if (k == "spbexact") SPB_EXACT = (int)d;
		else if (k == "mu") mu_rate = (std::string(val) == "inf") ? INF : d; else if (k == "theta") theta_rate = (std::string(val) == "inf") ? INF : d; else if (k == "nu") nu_rate = (std::string(val) == "inf") ? INF : d;
		else if (k == "tol") NEWTON_TOL = d; else if (k == "verbose") VERBOSE = (int)d; else if (k == "cavg") CENTRAL_CELLAVG = (int)d; else if (k == "isen") ISEN = (int)d; else if (k == "incdonor") INCDONOR = (int)d; else if (k == "acoustic") ACOUSTIC = (int)d; else if (k == "cspeed") CSPEED = (int)d; else if (k == "acord") ACORD = (int)d; else if (k == "stepact") STEPACT = (int)d; else if (k == "dump") DUMP = (int)d;
		else { fprintf(stderr, "unknown key %s\n", key); return 1; }
	}
	// 松弛档位的默认速率：relax>=1 时热松弛默认瞬时，relax>=2 化学默认瞬时（可由 theta=, nu= 覆盖为有限值）
	if (RELAX_LEVEL >= 1 && theta_rate == 0.0) theta_rate = INF;
	if ((RELAX_LEVEL == 2 || RELAX_LEVEL == 3) && nu_rate == 0.0) nu_rate = INF;
	setup_case(); if (tend_override > 0) time_stop = tend_override;
	if (getenv("SELFTEST")) {   // 闭合自检：p = 1e5 的状态，s = 0，应返回 1e5
		const double a1 = 1 - 1e-8, p0 = 1e5, r1 = 1000, r2 = 10;
		ClIn in{}; in.atr[0] = a1; in.atr[1] = 1 - a1; in.mtr[0] = a1 * r1; in.mtr[1] = (1 - a1) * r2;
		in.etr[0] = eps_amp(0, a1, in.mtr[0], p0); in.etr[1] = eps_amp(1, 1 - a1, in.mtr[1], p0);
		in.s = 0; in.rhoeC = in.etr[0] + in.etr[1]; in.tau = 1e-6; in.type = CL_P; in.pguess = p0;
		if (getenv("SELFTEST")[0] == '2') { in.atr[0] = 0.99999998999999995; in.atr[1] = 1.0000000050247593e-08; in.mtr[0] = 999.99998999999991; in.mtr[1] = 1.0000000049361959e-07;
			in.etr[0] = 776499992.2349999; in.etr[1] = 0.0025000000124212668; in.s = 0; in.rhoeC = 776499992.23749995; VERBOSE = 9; }
		ClOut o; const bool ok = closure(in, o);
		printf("selftest closure_p ok %d p %.12e  a %.12e %.12e\n", ok, o.pbar, o.a[0], o.a[1]);
		for (double gs : {1e4, 3e4, 9.9e4, 1e5, 1.01e5, 4e5, 1e7}) { in.pguess = gs; closure(in, o); printf("  guess %.3e -> p %.12e\n", gs, o.pbar); }
		return 0;
	}
	dx = (x_1 - x_0) / nxCells; dy = (y_1 - y_0) / nyCells;
	if (RELAX_LEVEL > 0 || TESTCASE == 40) build_sat();
	Grid U(nxCells, nyCells); init_case(U);
	double M0[3] = {0, 0, 0}; for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++) { M0[0] += U(i, j, IM1); M0[1] += U(i, j, IM2); M0[2] += U(i, j, IEN); }

	// Butcher 表（2DEuler_TV）
	const int RKmax = 3; double Ae[3][3] = {{0}}, Ai[3][3] = {{0}};
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
		if (MODE == 0) dt = CFL * caldt(U, true);
		else {
			dt = CFL * caldt(U, false);
			if (RAMP) { if (dtprev < 0) dtprev = CFL_ac * caldt(U, true) / 1.2; if (dt > 1.2 * dtprev) { dt = 1.2 * dtprev; ramping = true; } }
			if (!std::isfinite(dt)) dt = CFL_ac * caldt(U, true);
		}
		dt = std::min(dt, time_stop - t);
		bool done = false;
		for (int attempt = 0; attempt < 30 && !done; attempt++) {
			if (MODE == 0) {
				Grid V = U;
				if (SpaceOrder < 2) done = explicit_step(V, dt);
				else { Grid W = V; done = explicit_step(V, dt) && explicit_step(V, dt);
					if (done) { for (size_t q = 0; q < V.data.size(); q++) V.data[q] = 0.5 * (V.data[q] + W.data[q]);
						for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++) { double Q[NUM_VARS]; for (int k = 0; k < NUM_VARS; k++) Q[k] = V(i, j, k); if (!relax_cell_explicit(Q, dt)) { done = false; break; } for (int k = 0; k < NUM_VARS; k++) V(i, j, k) = Q[k]; }
						bound_cond_u(V); } }
				if (done) U = V;
			} else {
				const int rk = (RAMP && ramping) ? 1 : RKs; if (rk == 1) nrk1++;
				set_tab(rk);
				if (RELAX_LEVEL > 0 && STEPACT) { activate_and_project(U); bound_cond_u(U); }
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
			}
			if (!done) { st_reject++; dt *= 0.5; if (VERBOSE) fprintf(stderr, "reject: t %.4e new dt %.3e\n", t, dt); }
		}
		if (!done) { fprintf(stderr, "FAILED at t = %.6e (step %ld)\n", t, st_steps); output(U, out.c_str()); return 2; }
		t += dt; dtprev = dt; st_steps++;
		if (st_steps <= DUMP) { char fn[512]; snprintf(fn, sizeof fn, "%s.step%ld", out.c_str(), st_steps); output(U, fn); fprintf(stderr, "step %ld t %.6e dt %.4e\n", st_steps, t, dt); }
		if (VERBOSE > 1 && st_steps % 50 == 0) fprintf(stderr, "step %ld t %.4e dt %.3e newton it %ld\n", st_steps, t, dt, st_newton_it);
	}
	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double M1[3] = {0, 0, 0}; for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++) { M1[0] += U(i, j, IM1); M1[1] += U(i, j, IM2); M1[2] += U(i, j, IEN); }
	if (MODE == 1 && RELAX_LEVEL > 0 && STEPACT) { activate_and_project(U); bound_cond_u(U); }   // 输出前同样做激活投影（P-S：输出的是松弛后的态）
	output(U, out.c_str());
	printf("case %d mode %d order %d rk %d N %d CFL %g | steps %ld (rk1 %ld) rejects %ld | stages %ld newton it %ld (%.2f/stage) maxres %.2e | closure fails %ld | dM1/M %.2e dM2/M %.2e dE/E %.2e | wall %.2f s\n",
		TESTCASE, MODE, SpaceOrder, RKs, nxCells, CFL, st_steps, nrk1, st_reject, st_stage, st_newton_it, st_stage ? (double)st_newton_it / st_stage : 0.0, st_maxres, st_cl_fail,
		(M1[0] - M0[0]) / M0[0], (M1[1] - M0[1]) / M0[1], (M1[2] - M0[2]) / M0[2], wall);
	return 0;
}
