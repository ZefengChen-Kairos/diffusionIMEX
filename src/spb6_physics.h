// =====================================================================================================================
//  spb6_physics.h
//  SPB 2009 的 6+1 方程两相模型（单速度、双压力 + 冗余总能）的物理部分，显式与 IMEX 两个程序共用：
//    状态方程（NASG；b = 0 为 SG）、饱和线、标量求根、统一局部闭合（推导 B5、B6；速率 mu, theta, nu 为参数）、
//    松弛档位与激活判据、单元原始变量、算例（初值、区域、边界、终止时间）、输出用的单元诊断量。
//  两个程序都直接由 reference/code/2DEuler_TV.cpp 修改而来（网格、边界、时间步、输出、主循环在各自的文件里）：
//    SPB6_explicit.cpp  全显式参考格式（docs/derivation/A_model_EOS_explicit.md）
//    SPB6_TV_IMEX.cpp   TV 分裂 IMEX（docs/derivation/B_TV_IMEX.md）
//  所有算例用同一套格式与参数；算例只给出状态方程、初值、边界与终止时间。
//
//  守恒变量（NUM_VARS = 8）: 0 alpha1, 1 m1=(alpha rho)_1, 2 m2, 3 rho u, 4 rho v, 5 eps1=(alpha rho e)_1, 6 eps2, 7 E
//  eps_k 含生成焓 q_k；E = eps1 + eps2 + 1/2 rho |u|^2。
// =====================================================================================================================
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <limits>

#define NUM_VARS 8
enum { IA1 = 0, IM1 = 1, IM2 = 2, IMU = 3, IMV = 4, IE1 = 5, IE2 = 6, IEN = 7 };
static const double INF = std::numeric_limits<double>::infinity();

enum x_boundary_condition { x_Transmissive, x_Reflective, x_Periodic };
enum y_boundary_condition { y_Transmissive, y_Reflective, y_Periodic };

// 松弛设置（所有算例相同；命令行可改速率与档位）
double mu_rate = INF, theta_rate = -1.0, nu_rate = -1.0;   // 力学 / 热 / 化学 速率（INF = 瞬时；-1 = 未给定，按松弛档位取默认）
int RELAX_LEVEL = 0;        // 0: p;  1: p-pT;  2: p-pT-pTg;  3: p-pTg（P-S 的第四档）
double EPS_I = 1e-4;        // 热/化学松弛只在 eps_I <= alpha1 <= 1 - eps_I 的格子打开（P-S 2014、Pelanti 2022 8.2.2；所有算例相同）
int VERBOSE = 0;
long st_cl_fail = 0;

// ==================================================================================================== 状态方程（NASG；b = 0 为 SG）
// kind = 0: NASG（b = 0 为 SG），参数 g, pi, q, qp, cv, b。
// kind = 1: Cochran-Chan，kind = 2: JWL（A2.1，V12）：p = G rho (e - er(rho)) + pr(rho)，G 常数，er' = pr/rho^2。
//   CC：pr = c1 x^e1 - c2 x^e2，x = rho/rr；JWL：pr = c1 exp(-e1 rr/rho) + c2 exp(-e2 rr/rho)，er 再减 e0。cv > 0 时 T = (e - er)/cv（只作输出）
struct Phase { double g, pi, q, qp, cv, b; int kind = 0; double G = 0, c1 = 0, c2 = 0, e1 = 0, e2 = 0, rr = 0, e0 = 0; };
Phase PH[2];
inline bool is_mg(int k) { return PH[k].kind != 0; }
inline bool mg_any() { return is_mg(0) || is_mg(1); }
// MG 参考曲线：pr(rho)、pr'(rho)、er(rho)
inline void mg_ref(int k, double r, double& pr, double& dpr, double& er) {
	const Phase& P = PH[k];
	if (P.kind == 1) {
		const double x = r / P.rr, t1 = P.c1 * std::pow(x, P.e1), t2 = P.c2 * std::pow(x, P.e2);
		pr = t1 - t2; dpr = (P.e1 * t1 - P.e2 * t2) / r; er = (t1 / (P.e1 - 1.0) - t2 / (P.e2 - 1.0)) / r;
	} else {
		const double t1 = P.c1 * std::exp(-P.e1 * P.rr / r), t2 = P.c2 * std::exp(-P.e2 * P.rr / r);
		pr = t1 + t2; dpr = (P.e1 * t1 + P.e2 * t2) * P.rr / (r * r); er = t1 / (P.e1 * P.rr) + t2 / (P.e2 * P.rr) - P.e0;
	}
}

inline double eps_amp(int k, double a, double m, double p) { const Phase& P = PH[k];
	if (P.kind) { double pr, dpr, er; mg_ref(k, m / a, pr, dpr, er); return m * er + a * (p - pr) / P.G; }   // (A2.5)
	return (p + P.g * P.pi) * (a - m * P.b) / (P.g - 1.0) + m * P.q; }   // (A2.2)
inline double p_ame(int k, double a, double m, double e) { const Phase& P = PH[k];
	if (P.kind) { double pr, dpr, er; mg_ref(k, m / a, pr, dpr, er); return P.G * (e - m * er) / a + pr; }
	return (e - m * P.q) * (P.g - 1.0) / (a - m * P.b) - P.g * P.pi; }
inline double T_amp(int k, double a, double m, double p) { const Phase& P = PH[k];
	if (P.kind) { if (!(P.cv > 0)) return NAN; double pr, dpr, er; mg_ref(k, m / a, pr, dpr, er); return (a * (p - pr) / (P.G * m)) / P.cv; }
	return (a / m - P.b) * (p + P.pi) / ((P.g - 1.0) * P.cv); }
inline double v_pT(int k, double p, double T) { const Phase& P = PH[k]; return (P.g - 1.0) * P.cv * T / (p + P.pi) + P.b; }
inline double e_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.cv * T * (p + P.g * P.pi) / (p + P.pi) + P.q; }
inline double s_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.cv * (P.g * std::log(T) - (P.g - 1.0) * std::log(p + P.pi)) + P.qp; }
inline double g_pT(int k, double p, double T) { const Phase& P = PH[k]; return P.g * P.cv * T + P.b * p + P.q - T * s_pT(k, p, T); }
inline double c2_rp(int k, double r, double p) { const Phase& P = PH[k];
	if (P.kind) { double pr, dpr, er; mg_ref(k, r, pr, dpr, er); return (P.G + 1.0) * (p - pr) / r + dpr; }   // (A2.6)
	return P.g * (p + P.pi) / (r * (1.0 - r * P.b)); }
inline double dedp_am(int k, double a, double m) { const Phase& P = PH[k];
	if (P.kind) return a / P.G;
	return (a - m * P.b) / (P.g - 1.0); }    // d eps_k / d p |_{alpha,m}

// 混合压力 (A2.3)：给定 alpha_k, m_k 与混合内能
inline double p_mix(double a1, double m1, double m2, double rhoe) {
	const double a[2] = {a1, 1.0 - a1}, m[2] = {m1, m2};
	double num = rhoe, den = 0.0;
	for (int k = 0; k < 2; k++) {
		if (is_mg(k)) { double pr, dpr, er; mg_ref(k, m[k] / a[k], pr, dpr, er); num -= m[k] * er - a[k] * pr / PH[k].G; }   // (A2.7)
		else num -= m[k] * PH[k].q + PH[k].g * PH[k].pi * dedp_am(k, a[k], m[k]);
		den += dedp_am(k, a[k], m[k]);
	}
	return num / den;
}
// 相状态可容许：NASG 为 p + pi > 0；MG 为 c^2 > 0（A2.1）
inline bool p_admissible(int k, double a, double m, double p) {
	if (is_mg(k)) return std::isfinite(p) && c2_rp(k, m / a, p) > 0;
	return p + PH[k].pi > 0;
}
// 相内能 eps_k 的舍入尺度里除 |eps_k| 外的部分（B8 的 P_E）：NASG |m q|；MG |m er| + |alpha pr / G|
inline double eps_round(int k, double a, double m) {
	if (is_mg(k)) { double pr, dpr, er; mg_ref(k, m / a, pr, dpr, er); return std::fabs(m * er) + std::fabs(a * pr / PH[k].G); }
	return std::fabs(m * PH[k].q);
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
			if (p > 0) q = -q;
			p = std::fabs(p);
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
		if (!std::isfinite(fb)) return false;
		ya = yb; fa = fb;
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
	double rhoeC;                    // 守恒的混合内能 E - 1/2 rho u^2；与 sum etr 之差（能量亏损）按质量分数分给两相（B5.2）
	double tau;                      // 级的隐式步长 a_ii dt
	double Z[2];                     // 有限 mu 时的阻抗（冻结）
	int type;
	bool nodelta;                    // 不含能量亏损（只在 p-pTg 激活判据的预测里用）
	double pguess;
};
struct ClOut { double a[2], m[2], e[2], p[2], T[2], pbar; };
static const double YMIN = 1e-10;

// (B6.2)：p 闭合中 alpha_k(p) 的系数。能量亏损 delta（守恒混合内能与非守恒相内能之和的差，主要是被数值耗散的动能）
// 按质量分数 lam_k = Y_k 分给两相（B5.2）：单速度模型里第 k 相的动能是 1/2 alpha_k rho_k |u|^2（P-S 2014 §2.2 的相总能形式），
// 耗散的动能按质量分配。SPB 3.4 的重初始化（alpha 固定、整体平移 p）相当于 lam_k ∝ de_k/dp|alpha ∝ alpha_k/(gamma_k - 1)，
// 把耗散的动能主要交给可压缩性最大的相（气/汽），质量极小的轻相因此被加热到 1e5-1e6 K（devlog V10）。
struct PCoef { double A[2], B[2], N[2], lam[2]; };
static bool p_coef(const ClIn& in, PCoef& c) {
	double den = 0.0;
	for (int k = 0; k < 2; k++) { c.lam[k] = in.nodelta ? 0.0 : in.mtr[k]; den += c.lam[k]; }
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
static bool closure_p_mg(const ClIn& in, ClOut& o);
static bool closure_p(const ClIn& in, ClOut& o) {
	if (mg_any()) return closure_p_mg(in, o);
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

// B6.1 的 MG 推广（V12）：MG 相的 alpha_k(p) 是 f(alpha) = m er(m/alpha) + alpha (p - pr)/G + p alpha - A - B p 的根，
// f' = rho c^2 / G > 0；NASG 相仍用 (B6.2) 的闭式。外层 Phi(p) = sum alpha_k(p) - 1 用带区间的 Newton（区间从初猜向两侧扩张）
static bool alpha_mg(int k, double m, double A, double B, double p, double guess, double& al, double& dal) {
	const Phase& P = PH[k];
	auto F = [&](double a, double& df) { double pr, dpr, er; const double r = m / a; mg_ref(k, r, pr, dpr, er);
		df = ((P.G + 1.0) * (p - pr) + r * dpr) / P.G; return m * er + a * (p - pr) / P.G + p * a - A - B * p; };
	double lo = guess, hi = guess, d; double flo = F(lo, d), fhi = flo;
	if (!std::isfinite(flo)) return false;
	if (flo < 0) { for (int it = 0; it < 200 && fhi < 0; it++) { lo = hi; flo = fhi; hi *= 2.0; fhi = F(hi, d); if (!std::isfinite(fhi)) return false; } if (fhi < 0) return false; }
	else { for (int it = 0; it < 200 && flo > 0; it++) { hi = lo; fhi = flo; lo *= 0.5; flo = F(lo, d); if (!std::isfinite(flo) || d <= 0) return false; } if (flo > 0) return false; }
	double a = (flo == 0) ? lo : (fhi == 0 ? hi : (guess > lo && guess < hi ? guess : std::sqrt(lo * hi)));
	for (int it = 0; it < 200; it++) {
		double df; const double fa = F(a, df);
		if (fa == 0) break;
		if (fa < 0) lo = a; else hi = a;
		double an = (df > 0) ? a - fa / df : -1.0;
		if (!(an > lo && an < hi)) an = std::sqrt(lo * hi);
		if (std::fabs(an - a) <= 4e-16 * a || hi - lo <= 4e-16 * hi) { a = an; break; }
		a = an;
	}
	double df; F(a, df); if (!(df > 0)) return false;
	al = a; dal = -(a / P.G + a - B) / df; return true;
}
static bool closure_p_mg(const ClIn& in, ClOut& o) {
	PCoef c; if (!p_coef(in, c)) return false;
	double ag[2] = {in.atr[0] / (in.atr[0] + in.atr[1]), in.atr[1] / (in.atr[0] + in.atr[1])};
	auto Phi = [&](double p, double& dPhi, bool upd) -> double {   // NaN: 该 p 下某相无解（p 太低）
		double s = -1.0; dPhi = 0.0; double an[2];
		for (int k = 0; k < 2; k++) {
			double a, da;
			if (is_mg(k)) { if (!alpha_mg(k, in.mtr[k], c.A[k], c.B[k], p, ag[k], a, da)) return NAN; }   // MG 相 q = 0：A = eps^tr + lam D0
			else { if (!(p + PH[k].pi > 0) || !(c.N[k] > 0)) return NAN; a = alpha_p(c, in, k, p); da = dalpha_p(c, k, p); }
			an[k] = a; s += a; dPhi += da;
		}
		if (upd) { ag[0] = an[0]; ag[1] = an[1]; }
		return s; };
	double p = in.pguess, d;
	double f = Phi(p, d, true);
	double scale = 1e-6 * (std::fabs(p) + 1e5);
	double lo, hi, flo, fhi;
	if (!std::isfinite(f) || f > 0) {   // p 太低：向上找 Phi < 0
		lo = p; flo = f; hi = p;
		for (int it = 0; it < 300; it++) { hi = p + scale; fhi = Phi(hi, d, false); if (std::isfinite(fhi) && fhi <= 0) break; if (std::isfinite(fhi)) { lo = hi; flo = fhi; } scale *= 2.0; }
		if (!(std::isfinite(fhi) && fhi <= 0)) return false;
	} else {
		hi = p; fhi = f; lo = p;
		for (int it = 0; it < 300; it++) { lo = p - scale; flo = Phi(lo, d, false); if (std::isfinite(flo) && flo >= 0) break; if (!std::isfinite(flo)) { scale *= 0.5; if (scale < 1e-300) return false; continue; } hi = lo; fhi = flo; scale *= 2.0; }
		if (!(std::isfinite(flo) && flo >= 0)) return false;
	}
	p = (std::isfinite(f) && in.pguess >= lo && in.pguess <= hi) ? in.pguess : 0.5 * (lo + hi);
	for (int it = 0; it < 300; it++) {
		double df; const double fp = Phi(p, df, true);
		if (fp == 0.0) break;
		if (!std::isfinite(fp) || fp > 0) lo = p; else hi = p;
		double pn = (std::isfinite(fp) && df < 0) ? p - fp / df : 0.5 * (lo + hi);
		if (!(pn > lo && pn < hi)) pn = 0.5 * (lo + hi);
		if (std::fabs(pn - p) <= 2e-16 * (std::fabs(p) + std::fabs(pn)) || hi - lo <= 4e-16 * (std::fabs(hi) + std::fabs(lo))) { p = pn; break; }
		p = pn;
	}
	double df; if (!std::isfinite(Phi(p, df, true))) return false;
	for (int k = 0; k < 2; k++) {
		o.a[k] = ag[k]; o.m[k] = in.mtr[k];
		if (!p_admissible(k, o.a[k], o.m[k], p)) return false;
		o.e[k] = eps_amp(k, o.a[k], o.m[k], p); o.p[k] = p; o.T[k] = T_amp(k, o.a[k], o.m[k], p);
	}
	o.pbar = p; return std::isfinite(p);
}

// B6.2  有限 mu（theta = nu = 0）：给定 R 解 2x2 线性方程组得 (p1, p2)，再解 F(R) = R - tau mu (p1 - p2) = 0
static bool closure_mu(const ClIn& in, ClOut& o) {
	const double sa = in.atr[0] + in.atr[1];
	const double om[2] = {in.atr[0] / sa, in.atr[1] / sa};
	const double ak[2] = {in.atr[0] / sa, in.atr[1] / sa};   // a_k = atr_k + om_k s = atr_k / (1 - s)
	double lam[2], den = 0.0;
	for (int k = 0; k < 2; k++) { lam[k] = in.mtr[k]; den += lam[k]; }   // 能量亏损按质量分数分给两相（同 p_coef）
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
		// F(R) = R - tau mu (p1 - p2) 在可容许的连通区间上单调增（R 增大：相 1 膨胀降压、相 2 压缩升压）。
		// 松弛使 |p1 - p2| 减小，所以根在 0 与 R1 = tau mu (p1 - p2)(R=0) 之间。注意 (-a1, a2) 的端点不能当区间：
		// 压缩相的功按终态 p_I 计（后向 Euler），2x2 方程组在体积到 0 之前就奇异（SG 时 R = a2/gamma2），F 在端点不定号。
		// 若 R1 不可容许，奇点处 F -> 被压缩相一侧的无穷大（符号正确），向 0 二分即可得到可容许的区间。
		auto Fa = [&](double r, double& f) -> bool {
			double q1, q2; if (!solve(r, q1, q2)) return false;
			const double a1 = ak[0] + r, a2 = ak[1] - r;
			if (!(a1 > in.mtr[0] * PH[0].b && a2 > in.mtr[1] * PH[1].b && q1 + PH[0].pi > 0 && q2 + PH[1].pi > 0)) return false;
			f = r - in.tau * mu_rate * (q1 - q2); return std::isfinite(f); };
		double f0; if (!Fa(0.0, f0)) return false;
		if (f0 != 0.0) {
			double lo = 0.0, flo = f0, hi = -f0, fhi = NAN; bool ok = false;
			for (int it = 0; it < 200; it++) {
				if (Fa(hi, fhi)) { if (fhi * flo < 0) { ok = true; break; } lo = hi; flo = fhi; hi = lo + (lo - 0.0 == 0 ? -f0 : (lo)); }   // 同号：根在更远处（不应出现，防御）
				else hi = 0.5 * (lo + hi);   // 不可容许：向 lo 收缩
			}
			if (!ok) return false;
			auto F = [&](double r) { double f; return Fa(r, f) ? f : (double)NAN; };
			if (!brent(F, std::min(lo, hi), std::max(lo, hi), lo < hi ? flo : fhi, lo < hi ? fhi : flo, 1e-16 * (std::fabs(lo) + std::fabs(hi)) + 1e-300, R)) return false;
		}
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
		// 可容许区间：T(p, y) > 0 <=> eC - q1 - y (q2 - q1) > 0（分母为正）。q2 > q1 时 y 有上界、q2 < q1 时有下界。
		// T -> 0+ 时 g_k -> q_k + b_k p，R 在该端的符号由 q2 - q1 决定，与单调性（B6.5）一致，所以区间取可容许部分即可。
		double ya = YMIN, yb = 1.0 - YMIN;
		const double dq = PH[1].q - PH[0].q, yT = dq != 0.0 ? (eC - PH[0].q) / dq : INF;
		if (dq > 0) yb = std::min(yb, yT * (1.0 - 1e-12)); else if (dq < 0) ya = std::max(ya, yT * (1.0 + 1e-12));
		if (!(ya < yb)) return false;
		const double Ra = R(ya), Rb = R(yb);
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
	if (mg_any() && in.type != CL_P) return false;   // MG 相只做力学松弛（A2.1）
	switch (in.type) {
		case CL_P: ok = closure_p(in, o); break;
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

// ==================================================================================================== 单元原始变量
struct Prim { double a1, r1, r2, u, v, p1, p2, rho, cf; };
// 守恒量 Q -> 原始变量（不可容许时返回 false）
static bool prim_of(const double* Q, Prim& w) {
	const double a1 = Q[IA1], a2 = 1 - a1, m1 = Q[IM1], m2 = Q[IM2];
	if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return false;
	w.a1 = a1; w.rho = m1 + m2; w.r1 = m1 / a1; w.r2 = m2 / a2; w.u = Q[IMU] / w.rho; w.v = Q[IMV] / w.rho;
	w.p1 = p_ame(0, a1, m1, Q[IE1]); w.p2 = p_ame(1, a2, m2, Q[IE2]);
	if (!(p_admissible(0, a1, m1, w.p1) && p_admissible(1, a2, m2, w.p2))) return false;
	w.cf = std::sqrt((m1 * c2_rp(0, w.r1, w.p1) + m2 * c2_rp(1, w.r2, w.p2)) / w.rho);
	return true;
}
static void prim2cell(const Prim& w, double* Q) {
	const double a2 = 1 - w.a1, m1 = w.a1 * w.r1, m2 = a2 * w.r2, rho = m1 + m2;
	Q[IA1] = w.a1; Q[IM1] = m1; Q[IM2] = m2; Q[IMU] = rho * w.u; Q[IMV] = rho * w.v;
	Q[IE1] = eps_amp(0, w.a1, m1, w.p1); Q[IE2] = eps_amp(1, a2, m2, w.p2);
	Q[IEN] = Q[IE1] + Q[IE2] + 0.5 * rho * (w.u * w.u + w.v * w.v);
}
inline double pbar_of(const double* Q) {
	const double a1 = Q[IA1];
	return a1 * p_ame(0, a1, Q[IM1], Q[IE1]) + (1 - a1) * p_ame(1, 1 - a1, Q[IM2], Q[IE2]);
}

// ==================================================================================================== 算例（两个程序共用）
// 编号沿用测试脚本（case=N）；初值按论文，参数出处见 docs/results.md。
enum Test {
	testSPB_Cavitation = 0,       // SPB 2009 4.5   空化管，水/空气 alpha_air = 1e-2，-+100 m/s
	testPS_Cavitation2 = 1,       // P-S 2014 5.1   水空化管 -+2 m/s，T = 354.728 K
	testPS_Cavitation500 = 2,     // P-S 2014 5.1   水空化管 -+500 m/s
	testDodecaneSG = 3,           // P-S 2014 5.2   十二烷液-汽激波管（SG）
	testSPB_Interface = 4,        // SPB 2009 4.1   界面平移（WB）
	testSPB_WaterAir = 5,         // SPB 2009 4.3.1 水-气激波管 1 GPa
	testSPB_WaterAirTPa = 6,      // SPB 2009 4.3.2 水-气激波管 1 TPa
	testSPB_EpoxySpinel = 12,     // SPB 2009 5.1   环氧-尖晶石 1e10 Pa
	testSPB_EpoxySpinelExt = 13,  // SPB 2009 5.2   环氧-尖晶石 2e11 Pa
	testSmoothGas = 20,           // 光滑解（时间收敛阶），两相理想气体
	testSmoothWater = 21,         // 光滑解，水-气
	testPel_Cavitation = 30,      // Pelanti 2022 8.2.1 NASG 水空化管 -+2 m/s，T = 353 K
	testDodecaneNASG = 31,        // Pelanti 2022 8.2.2 十二烷液-汽激波管（NASG）
	testRelax0D = 40,             // 0D 松弛：均匀非平衡态（周期）
	testRelax0DChem = 41,         // 0D 化学松弛：pT 平衡、过热（周期）
	testGresho1 = 50,             // Gresho 涡（2DEuler 的 testGreshoVortex），单相极限（两相同为理想气体，alpha1 = 1/2），M = 0.1
	testGresho2 = 51,             //   M = 0.01
	testGresho3 = 52,             //   M = 0.001（低马赫 AP 检验）
	testGresho4 = 53,             //   M = 1e-4
	testGresho5 = 54,             //   M = 1e-5
	testGresho6 = 55,             //   M = 1e-6
	testSPB_Nitromethane = 60,    // SPB 2009 4.2  硝基甲烷激波管（Cochran-Chan，同一流体当两相），20 GPa / 0.2 MPa
	testPel_JWLWater = 61         // Pelanti 2022 8.1.1  爆轰气体（JWL）-水（NASG）激波管
};
struct TestSetup { double x0, x1, y0, y1, tend; x_boundary_condition xbc; y_boundary_condition ybc; };
static void set_SG(int k, double g, double pi, double q = 0, double qp = 0, double cv = 1000.0, double b = 0) { PH[k] = {g, pi, q, qp, cv, b}; }
static double rho_pT(int k, double p, double T) { return 1.0 / v_pT(k, p, T); }
// 状态方程、区域、边界、终止时间
static bool test_setup(int test, TestSetup& S) {
	S.x0 = 0; S.x1 = 1; S.y0 = 0; S.y1 = 1; S.xbc = x_Transmissive; S.ybc = y_Transmissive;
	switch (test) {
	case testSPB_Cavitation: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 1.85e-3; break;
	case testPS_Cavitation2: case testPS_Cavitation500:
		set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0);
		S.tend = test == testPS_Cavitation2 ? 3.2e-3 : 0.58e-3; break;
	case testDodecaneSG: set_SG(0, 2.35, 4e8, -775.269e3, 0.0, 1077.7); set_SG(1, 1.025, 0.0, -237.547e3, -24.4e3, 1956.45); S.tend = 473e-6; break;
	case testSPB_Interface: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 2.79e-3; break;
	case testSPB_WaterAir: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 240e-6; break;
	case testSPB_WaterAirTPa: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 8.3e-6; break;
	case testSPB_EpoxySpinel: case testSPB_EpoxySpinelExt: set_SG(0, 2.43, 5.3e9); set_SG(1, 1.62, 141e9); S.tend = test == testSPB_EpoxySpinel ? 80e-6 : 34e-6; break;
	case testSmoothGas: set_SG(0, 1.4, 0); set_SG(1, 1.4, 0); S.tend = 5e-4; break;
	case testSmoothWater: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 1e-4; break;
	case testPel_Cavitation: PH[0] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4}; PH[1] = {1.467, 0.0, 2077616.0, 14317.0, 955.0, 0.0}; S.tend = 3e-3; break;
	case testDodecaneNASG: PH[0] = {1.206, 1681e5, -996054.0, 0.0, 2532.0, 7.51e-4}; PH[1] = {1.021, 0.0, -384592.0, -4301.0, 2274.0, 0.0}; S.tend = 473e-6; break;
	case testRelax0D: case testRelax0DChem: set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0); S.tend = 1e-4; S.xbc = x_Periodic; break;
	case testGresho1: case testGresho2: case testGresho3: case testGresho4: case testGresho5: case testGresho6:
		set_SG(0, 1.4, 0); set_SG(1, 1.4, 0); S.x0 = -0.5; S.x1 = 0.5; S.y0 = -0.5; S.y1 = 0.5; S.tend = 0.4 * M_PI; S.xbc = x_Periodic; S.ybc = y_Periodic; break;
	case testSPB_Nitromethane:   // SPB 4.2：Gamma = 1.19，rho_ref = 1134，A1 = 0.819181e9，A2 = 1.50835e9，E1 = 4.52969，E2 = 1.42144；界面位置论文未写，取 0.5 m
		for (int k = 0; k < 2; k++) { PH[k] = Phase{}; PH[k].kind = 1; PH[k].G = 1.19; PH[k].rr = 1134.0; PH[k].c1 = 0.819181e9; PH[k].c2 = 1.50835e9; PH[k].e1 = 4.52969; PH[k].e2 = 1.42144; PH[k].g = 2.19; }
		S.tend = 67e-6; break;
	case testPel_JWLWater:       // Pelanti 表 2（JWL）与表 3（NASG 液态水）；[-10, 10] m，t = 0.9 ms
		PH[0] = Phase{}; PH[0].kind = 2; PH[0].G = 0.25; PH[0].rr = 1840.0; PH[0].c1 = 854.5e9; PH[0].c2 = 20.5e9; PH[0].e1 = 4.6; PH[0].e2 = 1.35; PH[0].e0 = 0.0; PH[0].g = 1.25;
		PH[1] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4};
		S.x0 = -10.0; S.x1 = 10.0; S.tend = 0.9e-3; break;
	default: return false;
	}
	return true;
}
// 点 (x, y) 处的初值（守恒量）
static void test_state(int test, double x, double y, double* Q) {
	(void)y;
	auto setcell = [&](double a1, double r1, double r2, double u, double v, double p1, double p2) {
		Prim w; w.a1 = a1; w.r1 = r1; w.r2 = r2; w.u = u; w.v = v; w.p1 = p1; w.p2 = p2; prim2cell(w, Q); };
	switch (test) {
	case testSPB_Cavitation: { const double u = x < 0.5 ? -100.0 : 100.0; setcell(1 - 1e-2, 1000.0, 1.0, u, 0, 1e5, 1e5); break; }   // 空气密度论文未给，取 1
	case testPS_Cavitation2: case testPS_Cavitation500: {
		const double T0 = 354.728, U0 = test == testPS_Cavitation2 ? 2.0 : 500.0, u = x < 0.5 ? -U0 : U0;
		setcell(1 - 1e-2, rho_pT(0, 1e5, T0), rho_pT(1, 1e5, T0), u, 0, 1e5, 1e5); break; }
	case testDodecaneSG: case testDodecaneNASG: { const bool L = x < 0.75; setcell(L ? 1 - 1e-8 : 1e-8, 500.0, 2.0, 0, 0, L ? 1e8 : 1e5, L ? 1e8 : 1e5); break; }
	case testSPB_Nitromethane: { const bool L = x < 0.5; setcell(L ? 1 - 1e-8 : 1e-8, 1134.0, 1200.0, 0, 0, L ? 20e9 : 0.2e6, L ? 20e9 : 0.2e6); break; }
	case testPel_JWLWater: { const bool L = x < 0.0; setcell(L ? 1 - 1e-8 : 1e-8, 2000.0, 1044.0, 0, 0, L ? 4.6406e10 : 1e5, L ? 4.6406e10 : 1e5); break; }
	case testSPB_Interface: { const bool L = x < 0.5; setcell(L ? 1 - 1e-8 : 1e-8, 1000.0, 10.0, 100.0, 0, 1e5, 1e5); break; }
	case testSPB_WaterAir: case testSPB_WaterAirTPa: {
		const bool L = x < 0.75; const double pL = test == testSPB_WaterAir ? 1e9 : 1e12, ra = test == testSPB_WaterAir ? 1.0 : 10.0;
		setcell(L ? 1 - 1e-6 : 1e-6, 1000.0, ra, 0, 0, L ? pL : 1e5, L ? pL : 1e5); break; }
	case testSPB_EpoxySpinel: case testSPB_EpoxySpinelExt: {
		const double pL = test == testSPB_EpoxySpinel ? 1e10 : 2e11; const double p = x < 0.6 ? pL : 1e5;
		setcell(0.5954, 1185.0, 3622.0, 0, 0, p, p); break; }
	case testSmoothGas: case testSmoothWater: {
		const double gs = std::exp(-(x - 0.5) * (x - 0.5) / 0.01), u = 50.0 * gs, a1 = 0.5 + 0.2 * std::sin(2 * M_PI * x), p = 1e5 * (1 + 0.2 * gs);
		setcell(a1, test == testSmoothGas ? 2.0 : 1000.0, 1.0, u, 0, p, p); break; }
	case testPel_Cavitation: { const double T0 = 353.0, u = x < 0.5 ? -2.0 : 2.0; setcell(1 - 1e-2, rho_pT(0, 1e5, T0), rho_pT(1, 1e5, T0), u, 0, 1e5, 1e5); break; }
	case testRelax0DChem: { const double T0 = 380.0, p0 = 0.8e5; setcell(1 - 0.2, rho_pT(0, p0, T0), rho_pT(1, p0, T0), 0, 0, p0, p0); break; }
	case testRelax0D: { const double T1 = 350.0, T2 = 380.0, p1 = 1.2e5, p2 = 0.8e5; setcell(1 - 0.2, rho_pT(0, p1, T1), rho_pT(1, p2, T2), 0, 0, p1, p2); break; }
	case testGresho1: case testGresho2: case testGresho3: case testGresho4: case testGresho5: case testGresho6: {   // 2DEuler_TV.cpp 的 testGreshoVortex：rho = 1，p0 = 1/(gamma M^2)，最大速度 1（r = 0.2）
		static const double Ms[6] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
		const double M = Ms[test - testGresho1], p0 = 1.0 / (1.4 * M * M), r = std::sqrt(x * x + y * y);
		double vel, p;
		if (r < 0.2) { vel = 5 * r; p = p0 + 25.0 / 2.0 * r * r; }
		else if (r < 0.4) { vel = 2 - 5 * r; p = p0 + 25.0 / 2.0 * r * r + 4 * (1 - 5 * r - std::log(0.2) + std::log(r)); }
		else { vel = 0; p = p0 - 2 + 4 * std::log(2.0); }
		const double ux = r < 1e-14 ? 0.0 : -vel * y / r, uy = r < 1e-14 ? 0.0 : vel * x / r;
		setcell(0.5, 1.0, 1.0, ux, uy, p, p); break; }
	}
}

// ==================================================================================================== 输出用的单元诊断量
// 列：p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2
static void cell_diag(const double* Q, double* o) {
	Prim w; const bool ok = prim_of(Q, w);
	const double a1 = Q[IA1], m1 = Q[IM1], m2 = Q[IM2], rho = m1 + m2;
	double cw = NAN, T1 = NAN, T2 = NAN, dg = NAN;
	if (ok) {
		cw = std::sqrt(1.0 / (rho * (a1 / (w.r1 * c2_rp(0, w.r1, w.p1)) + (1 - a1) / (w.r2 * c2_rp(1, w.r2, w.p2)))));
		T1 = T_amp(0, a1, m1, w.p1); T2 = T_amp(1, 1 - a1, m2, w.p2);
		if (T1 > 0 && T2 > 0 && PH[0].cv > 0 && !mg_any()) dg = g_pT(0, w.p1, T1) - g_pT(1, w.p2, T2);
	}
	o[0] = ok ? a1 * w.p1 + (1 - a1) * w.p2 : NAN; o[1] = ok ? w.u : NAN; o[2] = ok ? w.v : NAN; o[3] = rho; o[4] = 1 - a1; o[5] = m2 / rho;
	o[6] = ok ? w.r1 : NAN; o[7] = ok ? w.r2 : NAN; o[8] = ok ? w.cf : NAN; o[9] = cw; o[10] = T1; o[11] = T2; o[12] = ok ? w.p1 : NAN; o[13] = ok ? w.p2 : NAN; o[14] = dg;
}
