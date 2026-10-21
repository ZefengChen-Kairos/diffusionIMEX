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
double mu_rate = INF, theta_rate = -1.0, nu_rate = -1.0;
// 随状态变化的传质速率（Pelanti 2022 式 (84)，Downar-Zapolski / De Lorenzo 型；命令行 nulaw=1）：
//   nu = NU_CR alpha_v^0.6 ((psat(T) - p) / (NU_PC - psat(T)))^1.76（p < psat(T)，否则 0），alpha_v 为相 2（汽）的体积分数
int NU_LAW = 0; double NU_CR = 1.0, NU_PC = 22.064e6;   // 力学 / 热 / 化学 速率（INF = 瞬时；-1 = 未给定，按松弛档位取默认）
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
// 饱和压力 psat(T)：Tsat 的反函数（Tsat 在表的范围内单调增），在 ln p 上二分；T 低于表的下端返回 0，高于临界点返回 INF
static double psat(double T) {
	if (!SAT.built) build_sat();
	if (SAT.lp.empty()) return NAN;
	double a = SAT.pmin, b = SAT.pcrit * (1 - 1e-12);
	const double Ta = Tsat(a), Tb = Tsat(b);
	if (!(T > Ta)) return 0.0; if (!(T < Tb)) return INF;
	for (int it = 0; it < 80; it++) { const double c = std::sqrt(a * b); if (Tsat(c) < T) a = c; else b = c; if (b / a - 1 < 1e-14) break; }
	return std::sqrt(a * b);
}
// Pelanti 2022 式 (84)：nu = NU_CR alpha_v^0.6 ((psat(T) - p) / (NU_PC - psat(T)))^1.76，只在 p < psat(T)（过热液体，蒸发）时非零
static double nu_law(double p, double T, double av) {
	const double ps = psat(T);
	if (!std::isfinite(ps) || !(ps > p) || !(ps < NU_PC) || !(av > 0)) return 0.0;
	return NU_CR * std::pow(av, 0.6) * std::pow((ps - p) / (NU_PC - ps), 1.76);
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
	double rhoeC;                    // 守恒的混合内能 E - 1/2 rho u^2；与 sum etr 之差（能量亏损）分给两相（B5.2，defect_split）
	double tau;                      // 级的隐式步长 a_ii dt
	double Z[2];                     // 有限 mu 时的阻抗（冻结）
	int type;
	bool nodelta;                    // 不含能量亏损（只在 p-pTg 激活判据的预测里用）
	double pguess;
	double Qph[2] = {0.0, 0.0};      // 本格局部激波的物理耗散预算 Q_k（每单位体积，shock_budget）；0 = 没有物理部分（defect_split）
	double nu = -1.0;                // 本格的传质速率（>= 0 时代替 nu_rate；nulaw=1 时按步开始的状态给出）
};
struct ClOut { double a[2], m[2], e[2], p[2], T[2], pbar; };
static const double YMIN = 1e-10;

// (B5.2) 能量亏损的分配（V16）。delta = (rho e)^C - sum eps^tr + p s 是闭合里的真实亏损：激波处是物理耗散（格式的数值粘性把动能变成
// 内能），光滑区是离散误差（正负交替）。两者的分法不同：
// - 物理耗散按模型自己的激波关系分。5 方程模型的激波闭合是 SPB (II.5)（Saurel, Le Métayer, Massoni, Gavrilyuk 2007）：每相走自己的
//   Hugoniot，e_k - e_k0 = -(p0 + p1)/2 (v_k - v_k0)。剖面分析（devlog V16）：格式的极限激波满足 (II.5) 当且仅当第 k 相得到的热是
//   Q_k = m_k (v_k^pre - b_k) f_k(p_lo -> p_hi)，m_k 为本格的质量，f_k(v_pre - b) = 同一终压下相 k 每单位质量的 Hugoniot 能量减等熵能量
//   （弱激波极限 alpha_k G_k dp^3/(6 K_k^2)，G_k = (gamma_k + 1)/2，K_k = rho_k c_k^2）。激波强度 p_lo -> p_hi 取格式自身的紧致模板（±1，二维十字）上混合压力的极值，
//   波前的相比容取模板里压力最低的格。
// - 离散误差按各相可容许界之上的能量 E_k 分（V14 的负亏损规则推广到两个符号）：各相的 p_k + pi_k 同乘 (1 + delta/sum E)，两相的相对
//   误差相同，微量相既不会被加热也不会被扣成负值。质量分配（V10–V15）是格式数值粘性的分法，不是模型的：它让 SPB 5.2 的激波后密度偏 5%，
//   按 (II.5) 分后偏 0.3%。SPB 3.4 的重初始化相当于 lam_k ∝ alpha_k/(gamma_k - 1)，微量相被加热到 1e5–1e6 K（devlog V10）。
// 物理部分有预算：压缩格里 delta 中至多 Q = sum Q_k 的部分是物理耗散（一次激波通过产生的亏损之和恰为 Q），其余是误差：
//   Lambda_k(delta) = min(delta, Q) lam^(II.5)_k + (delta - min(delta, Q)) lam^E_k   (delta >= 0)，   = delta lam^E_k   (delta < 0)。
// Lambda 是 delta 的连续分段仿射函数，所以 G(P) 连续（V15 按 D0 符号分支、又乘 s 的不连续与 active-set 校正都不再需要）。
// 闭合在 delta 线性的每一段上用 A_k + B_k p = eps^tr_k - m_k q_k + p alpha^tr_k + Lambda_k(delta(p))：p 闭合在每个试探 p 处选段（closure_p），
// 其余闭合按 delta(pguess) 选段，收敛时 pguess = p。

// ------------------------------------------------------------------ 激波的物理耗散与相份额（B5.2 (b)，V22）
// 松弛极限（Kapila）的定常激波剖面：Rayleigh 线 P(v) = p0 + j^2 (v0 - v)，混合耗散 dq = (p - P) dv，按常数份额 lam_k 分给各相，
//   de_k = -p dv_k + (lam_k / Y_k) dq,   sum_k Y_k v_k = v。
// 份额 c = lam_1 是使剖面终态落在 (II.5)（每相各自的 Hugoniot）上的唯一常数；混合终态与份额无关，所以积分区间 [v0, v_end] 固定，对 c 打靶。
// 弱激波（dp < 1e-2 min_k K_k）取解析极限（定理 8）：Y_k q_k = Y_k v_k G_k dp^3 / (6 K_k^2)。
// 相 k 在 (p, v) 处的 de/dp|_v 与 de/dv|_p + p
static inline void eos_ep_ev(int k, double p, double v, double& ep, double& evp) {
	const Phase& P = PH[k];
	if (P.kind) { double pr, dpr, er; const double r = 1.0 / v; mg_ref(k, r, pr, dpr, er); ep = v / P.G; evp = (p - pr) * (1.0 + 1.0 / P.G) + r * dpr / P.G; return; }
	ep = (v - P.b) / (P.g - 1.0); evp = P.g * (p + P.pi) / (P.g - 1.0);
}
static inline bool eos_ok(int k, double p, double v) {
	if (!(v > 0) || !std::isfinite(p)) return false;
	if (is_mg(k)) return c2_rp(k, 1.0 / v, p) > 0;
	return v > PH[k].b && p + PH[k].pi > 0;
}
// 相 k 的 Hugoniot (II.5)：从 (p0, v0) 到 p1 的比容。NASG 闭式；MG 在 (0, v0] 上二分
static double hug_v(int k, double p0, double p1, double v0) {
	const Phase& P = PH[k];
	if (!P.kind) {
		const double x0 = p0 + P.pi, x1 = p1 + P.pi;
		return P.b + (v0 - P.b) * ((P.g - 1.0) * x1 + (P.g + 1.0) * x0) / ((P.g + 1.0) * x1 + (P.g - 1.0) * x0);
	}
	auto e_pv = [&](double p, double v) { double pr, dpr, er; mg_ref(k, 1.0 / v, pr, dpr, er); return er + v * (p - pr) / P.G; };
	const double e0 = e_pv(p0, v0);
	auto g = [&](double v) { return e_pv(p1, v) - e0 + 0.5 * (p0 + p1) * (v - v0); };
	double hi = v0, lo = v0, gl = g(v0); if (!(gl > 0)) return v0;
	for (int it = 0; it < 80 && gl > 0; it++) { hi = lo; lo *= 0.9; gl = g(lo); if (!std::isfinite(gl)) return NAN; }
	if (gl > 0) return NAN;
	for (int it = 0; it < 100; it++) { const double m = 0.5 * (lo + hi), gm = g(m); if (gm > 0) hi = m; else lo = m; if (hi - lo <= 1e-14 * v0) break; }
	return 0.5 * (lo + hi);
}
// 等熵体积模量 K = rho c^2 与基本导数 G = (1 + dK/dp|_s)/2（MG 沿等熵线 dv = -v dp / K 差分）
static void K_G(int k, double p, double v, double& K, double& G) {
	K = c2_rp(k, 1.0 / v, p) / v;
	if (!is_mg(k)) { G = 0.5 * (PH[k].g + 1.0); return; }
	const double dp = 1e-6 * K, K2 = c2_rp(k, 1.0 / (v * (1.0 - dp / K)), p + dp) / (v * (1.0 - dp / K));
	G = 0.5 * (1.0 + (K2 - K) / dp);
}
// 给定份额 c 积分剖面到 v_end（自适应 Cash-Karp RK45）；返回终态 v_1 与混合耗散 q = int (p - P) dv（每单位混合质量）
struct ShockProf { double Y[2], vk0[2], v0, p0, j2, vend; };
static bool prof_rhs(const ShockProf& S, double c, double v, const double* y, double* dy) {
	const double p = y[0], v1 = y[1], v2 = (v - S.Y[0] * v1) / S.Y[1], P = S.p0 + S.j2 * (S.v0 - v);
	if (!eos_ok(0, p, v1) || !eos_ok(1, p, v2)) return false;
	const double lam[2] = {c, 1.0 - c}, vk[2] = {v1, v2};
	double num = -1.0, den = 0.0, r[2], ep[2], evp[2];
	for (int k = 0; k < 2; k++) { eos_ep_ev(k, p, vk[k], ep[k], evp[k]); r[k] = lam[k] / S.Y[k] * (p - P); num += S.Y[k] * r[k] / evp[k]; den += S.Y[k] * ep[k] / evp[k]; }
	dy[0] = num / den; dy[1] = (r[0] - ep[0] * dy[0]) / evp[0]; dy[2] = p - P;
	return std::isfinite(dy[0]) && std::isfinite(dy[1]);
}
static bool prof_integrate(const ShockProf& S, double c, double& v1end, double& qend) {
	static const double A[6] = {0, 0.2, 0.3, 0.6, 1.0, 0.875};
	static const double B[6][5] = {{0}, {0.2}, {3.0 / 40, 9.0 / 40}, {0.3, -0.9, 1.2}, {-11.0 / 54, 2.5, -70.0 / 27, 35.0 / 27},
		{1631.0 / 55296, 175.0 / 512, 575.0 / 13824, 44275.0 / 110592, 253.0 / 4096}};
	static const double C5[6] = {37.0 / 378, 0, 250.0 / 621, 125.0 / 594, 0, 512.0 / 1771};
	static const double C4[6] = {2825.0 / 27648, 0, 18575.0 / 48384, 13525.0 / 55296, 277.0 / 14336, 0.25};
	double y[3] = {S.p0, S.vk0[0], 0.0}, v = S.v0;
	const double L = S.vend - S.v0, sc[3] = {std::fabs(S.p0) + S.j2 * std::fabs(L), S.vk0[0], S.j2 * L * L + 1e-300};
	double h = L / 32.0;
	for (int it = 0; it < 4000; it++) {
		const bool last = (L < 0 && v + h <= S.vend) || (L > 0 && v + h >= S.vend);
		if (last) h = S.vend - v;
		double k[6][3], yt[3];
		bool ok = true;
		for (int st = 0; st < 6 && ok; st++) {
			for (int n = 0; n < 3; n++) { yt[n] = y[n]; for (int q = 0; q < st; q++) yt[n] += h * B[st][q] * k[q][n]; }
			ok = prof_rhs(S, c, v + A[st] * h, yt, k[st]);
		}
		double err = 0.0, y5[3];
		if (ok) for (int n = 0; n < 3; n++) {
			double d5 = 0, d4 = 0; for (int st = 0; st < 6; st++) { d5 += C5[st] * k[st][n]; d4 += C4[st] * k[st][n]; }
			y5[n] = y[n] + h * d5; err = std::max(err, std::fabs(h * (d5 - d4)) / sc[n]);
		}
		if (!ok || !(err <= 1e-8)) { h *= ok ? std::max(0.1, 0.9 * std::pow(1e-8 / err, 0.25)) : 0.25; if (std::fabs(h) < 1e-14 * std::fabs(L)) return false; continue; }
		for (int n = 0; n < 3; n++) y[n] = y5[n];
		if (last) { v1end = y[1]; qend = y[2]; return true; }
		v += h;
		h *= std::min(5.0, 0.9 * std::pow(1e-8 / std::max(err, 1e-30), 0.2));
	}
	return false;
}
// 波前 (p0, v_k^pre, Y_k) 被压到 p1：相份额 lam_k 与混合耗散 q（每单位混合质量）。失败返回 false
static bool shock_share(const double Y[2], const double vpre[2], double p0, double p1, double lam[2], double& q) {
	double K[2], G[2], w[2], sw = 0.0, Kmin = INF;
	for (int k = 0; k < 2; k++) { K_G(k, p0, vpre[k], K[k], G[k]); if (!(K[k] > 0)) return false; Kmin = std::min(Kmin, K[k]); }
	const double dp = p1 - p0;
	if (dp < 1e-2 * Kmin) {    // 弱激波：定理 8
		q = 0.0;
		for (int k = 0; k < 2; k++) { w[k] = Y[k] * vpre[k] * G[k] / (K[k] * K[k]); sw += w[k]; }
		for (int k = 0; k < 2; k++) lam[k] = w[k] / sw;
		q = sw * dp * dp * dp / 6.0;
		return true;
	}
	ShockProf S; S.p0 = p0;
	double vH[2];
	for (int k = 0; k < 2; k++) { S.Y[k] = Y[k]; S.vk0[k] = vpre[k]; vH[k] = hug_v(k, p0, p1, vpre[k]); if (!std::isfinite(vH[k])) return false; }
	S.v0 = Y[0] * vpre[0] + Y[1] * vpre[1]; S.vend = Y[0] * vH[0] + Y[1] * vH[1];
	if (!(S.vend < S.v0)) return false;
	S.j2 = dp / (S.v0 - S.vend);
	// 残量 F(c) = v_1(v_end; c) - v_1H，对 c 递增；在 logit(c) 上用 Illinois 法
	auto F = [&](double t, double& qq) -> double { const double c = 1.0 / (1.0 + std::exp(-t)); double v1e;
		if (!prof_integrate(S, c, v1e, qq)) return NAN; return (v1e - vH[0]) / vpre[0]; };
	double ta = -40.0, tb = 40.0, qa, qb, fa = F(ta, qa), fb = F(tb, qb);
	if (!std::isfinite(fa) || !std::isfinite(fb)) return false;
	double t, qt;
	if (fa >= 0) { t = ta; qt = qa; }
	else if (fb <= 0) { t = tb; qt = qb; }
	else {
		int side = 0; t = ta; qt = qa;
		for (int it = 0; it < 100; it++) {
			t = (ta * fb - tb * fa) / (fb - fa);
			if (!(t > ta && t < tb)) t = 0.5 * (ta + tb);
			const double ft = F(t, qt);
			if (!std::isfinite(ft)) return false;
			if (std::fabs(ft) < 1e-10 || tb - ta < 1e-10) break;
			if (ft < 0) { ta = t; fa = ft; if (side == -1) fb *= 0.5; side = -1; }
			else { tb = t; fb = ft; if (side == 1) fa *= 0.5; side = 1; }
		}
	}
	const double c = 1.0 / (1.0 + std::exp(-t));
	lam[0] = c; lam[1] = 1.0 - c; q = qt;
	return q > 0;
}
// 守恒量 c（IA1..IEN）-> 模板状态 st = {a1, a2, m1, m2, p, u, v}；p 是混合压力 (A2.3)，不可容许时 NAN（模板里跳过）
inline void stencil_state(const double* c, double* st) {
	const double a1 = c[IA1], m1 = c[IM1], m2 = c[IM2], r = m1 + m2;
	st[0] = a1; st[1] = 1.0 - a1; st[2] = m1; st[3] = m2; st[4] = NAN;
	st[5] = r > 0 ? c[IMU] / r : 0.0; st[6] = r > 0 ? c[IMV] / r : 0.0;
	if (a1 > 0 && a1 < 1 && m1 > 0 && m2 > 0) { const double p = p_mix(a1, m1, m2, c[IE1] + c[IE2]); if (p_admissible(0, a1, m1, p) && p_admissible(1, 1.0 - a1, m2, p)) st[4] = p; }
}
// 紧致模板（±1，二维为十字；V21，半径 1–4 结果相同到 0.1%，devlog V21）上的激波预算：cell(di, dj, st) 给出模板格的状态。压缩（中心差分散度 < 0）且模板上有压力差时，本格的物质（自己的 m_k）
// 按 (II.5) 应得 Q_k = lam_k rho q（shock_share：剖面的精确常数份额与混合耗散，V22）；激波强度 p_lo -> p_hi 取模板上混合压力的极值，波前的相比容
// v_k^pre = alpha_k/m_k 取模板里压力最低的格（相不存在时取本格）。否则 Q = 0（全部按误差分）。
// 用本格自己的质量而不是波前格的 alpha_k：激波穿过接触/混合区时，被压缩的是本格的物质（均匀组分时两者的比例相同）
template <class F> static void shock_budget(F cell, bool twoD, double Q[2]) {
	Q[0] = Q[1] = 0.0;
	double c0[7], st[7], lo[7] = {0, 0, 0, 0, 0, 0, 0}, plo = 1e300, phi = -1e300, div = 0.0;
	cell(0, 0, c0);
	for (int d = 0; d < (twoD ? 2 : 1); d++) for (int q = -1; q <= 1; q++) {
		if (d == 1 && q == 0) continue;
		cell(d == 0 ? q : 0, d == 1 ? q : 0, st);
		if (q == 1) div += st[5 + d]; else if (q == -1) div -= st[5 + d];
		if (!std::isfinite(st[4])) continue;
		if (st[4] < plo) { plo = st[4]; for (int n = 0; n < 7; n++) lo[n] = st[n]; }
		if (st[4] > phi) phi = st[4];
	}
	if (!(div < 0.0) || !(phi > plo)) return;
	double Y[2], vpre[2], lam[2], q;
	const double rho = c0[2] + c0[3];
	if (!(rho > 0)) return;
	for (int k = 0; k < 2; k++) {
		Y[k] = c0[2 + k] / rho;
		vpre[k] = (lo[k] > 0 && lo[2 + k] > 0) ? lo[k] / lo[2 + k] : (c0[2 + k] > 0 ? c0[k] / c0[2 + k] : 0.0);
	}
	if (!(Y[0] > 0) || !(Y[1] > 0) || !(vpre[0] > 0) || !(vpre[1] > 0)) return;
	if (!shock_share(Y, vpre, plo, phi, lam, q)) return;
	for (int k = 0; k < 2; k++) Q[k] = lam[k] * rho * q;
}
// 可容许界之上的能量 E_k：NASG (p_k + pi_k)(alpha_k - m_k b_k)/(gamma_k - 1)；MG 冷曲线之上的热能 eps_k - m_k er
static double E_above(const ClIn& in, int k) {
	if (is_mg(k)) { if (!(in.atr[k] > 0 && in.mtr[k] > 0)) return 0.0; double pr, dpr, er; mg_ref(k, in.mtr[k] / in.atr[k], pr, dpr, er); return std::max(0.0, in.etr[k] - in.mtr[k] * er); }
	return std::max(0.0, in.etr[k] - in.mtr[k] * PH[k].q - (in.atr[k] - in.mtr[k] * PH[k].b) * PH[k].pi);
}
// 亏损分配的消融开关（只用于论文的对比图，默认 0）：0 = 按激波关系（物理份额 + 误差份额 lambda^E）；
// 1 = 没有物理份额（全部按 lambda^E）；2 = 全部按质量分数 Y_k
static int SHARE_MODE = 0;
// Lambda_k(delta) 在 delta 所在分段上的仿射系数：Lambda_k = lam_k delta + extra_k
static void defect_split(const ClIn& in, double delta, double lam[2], double extra[2]) {
	if (SHARE_MODE == 2) { const double sm = in.mtr[0] + in.mtr[1]; for (int k = 0; k < 2; k++) { lam[k] = sm > 0 ? in.mtr[k] / sm : 0.5; extra[k] = 0.0; } return; }
	double E[2], sE = 0.0; for (int k = 0; k < 2; k++) { E[k] = E_above(in, k); sE += E[k]; }
	double lE[2];
	if (sE > 0) for (int k = 0; k < 2; k++) lE[k] = E[k] / sE;
	else { const double sm = in.mtr[0] + in.mtr[1]; for (int k = 0; k < 2; k++) lE[k] = sm > 0 ? in.mtr[k] / sm : 0.5; }   // 两相都在可容许界上（退化）：按质量
	const double Q = SHARE_MODE == 1 ? 0.0 : in.Qph[0] + in.Qph[1];
	extra[0] = extra[1] = 0.0;
	if (!(delta > 0.0) || !(Q > 0.0)) { lam[0] = lE[0]; lam[1] = lE[1]; return; }
	if (delta < Q) { lam[0] = in.Qph[0] / Q; lam[1] = in.Qph[1] / Q; return; }
	for (int k = 0; k < 2; k++) { lam[k] = lE[k]; extra[k] = in.Qph[k] - Q * lE[k]; }
}
struct PCoef { double A[2], B[2], N[2], lam[2]; };
// 分段按 delta(pseg) = D0 + pseg s 取（pseg < 0：取 pguess）；在该段上 Lambda_k = lam_k delta + xi_k 对 p 仿射，A_k + B_k p 含它的 p 依赖，
// 于是对任意 p，sum_k eps_k = (rho e)^C 严格成立（能量一致性在闭合内，不依赖全局迭代收敛）
static bool p_coef(const ClIn& in, PCoef& c, double pseg = -1e300) {
	if (!in.nodelta && !(in.mtr[0] + in.mtr[1] > 0)) return false;
	const double D0 = in.nodelta ? 0.0 : in.rhoeC - in.etr[0] - in.etr[1];
	double ex[2] = {0.0, 0.0};
	if (in.nodelta) c.lam[0] = c.lam[1] = 0.0; else defect_split(in, D0 + (pseg > -1e299 ? pseg : in.pguess) * in.s, c.lam, ex);
	for (int k = 0; k < 2; k++) {
		c.A[k] = in.etr[k] - in.mtr[k] * PH[k].q + c.lam[k] * D0 + ex[k];
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
	// Lambda_k(delta(p)) 的分段在每个试探 p 处取（B5.2 (e)），Phi(p) = sum alpha_k(p) - 1 连续。单相的 alpha_k(p) 不必单调：
	// 激波处液格里的微量气体在 (II.5) 段分到 lam ~ 1，D0 ~ -p s 时其 N_k < 0（alpha_gas 随 p 增），但和 Phi 仍由液体主导而递减。
	// 所以只要求 Phi 有带区间的根、解处各相可容许，不再逐相要求 N_k > 0（那样空化 RM 在 t = 0.73 ms 无解，devlog V16 §16.8）
	PCoef c; if (!p_coef(in, c)) return false;
	const double plo = plo_mix();
	auto Phi = [&](double p) { PCoef q; p_coef(in, q, p); return alpha_p(q, in, 0, p) + alpha_p(q, in, 1, p) - 1.0; };
	auto dPhi = [&](double p) { PCoef q; p_coef(in, q, p); return dalpha_p(q, 0, p) + dalpha_p(q, 1, p); };
	// 区间：lo 处 Phi > 0，hi 处 Phi < 0（从初猜向两侧几何扩张）
	double d = std::max(in.pguess - plo, 1e-3 * (std::fabs(plo) + 1.0));
	double lo = plo, hi = plo + d; double fh = Phi(hi);
	if (fh > 0) { lo = hi; for (int k = 0; k < 200 && fh > 0; k++) { d *= 4.0; hi = plo + d; fh = Phi(hi); } if (fh > 0) return false; }
	else { bool found = false; for (int k = 0; k < 200; k++) { const double dd = d * 0.25, pl = plo + dd, fl = Phi(pl); if (fl > 0) { lo = pl; found = true; break; } hi = pl; fh = fl; d = dd; } if (!found) return false; }
	double p = 0.5 * (lo + hi); if (in.pguess > lo && in.pguess < hi) p = in.pguess;
	for (int it = 0; it < 200; it++) {
		const double f = Phi(p);
		if (f == 0.0) break;                                   // 恰为根
		if (f > 0) lo = p; else hi = p;
		const double df = dPhi(p);
		double pn = df < 0 ? p - f / df : 0.5 * (lo + hi);
		if (!(pn > lo && pn < hi)) pn = (hi - plo > 1e3 * (lo - plo + 1e-300) && lo > plo) ? plo + std::sqrt((lo - plo) * (hi - plo)) : 0.5 * (lo + hi);
		if (std::fabs(pn - p) <= 2e-16 * (std::fabs(p) + std::fabs(plo)) || hi - lo <= 4e-16 * (std::fabs(hi) + std::fabs(lo))) { p = pn; break; }
		p = pn;
	}
	p_coef(in, c, p);
	for (int k = 0; k < 2; k++) {
		o.a[k] = alpha_p(c, in, k, p); o.m[k] = in.mtr[k];
		if (!(o.a[k] > o.m[k] * PH[k].b) || !(p + PH[k].pi > 0)) return false;
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
			if (is_mg(k)) { if (!alpha_mg(k, in.mtr[k], c.A[k], c.B[k], p, ag[k], a, da)) return NAN; }   // MG 相 q = 0：A = eps^tr + lam D0 + xi
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
	const double D0 = in.rhoeC - in.etr[0] - in.etr[1];
	// 两压力时 delta = D0 + s (om1 p1 + om2 p2) 不是 D0 + s pbar，所以 Lambda_k 的 p 依赖留在 2x2 方程组里（分段按 delta(pguess) 取），保持能量严格一致
	double lam[2], ex[2]; defect_split(in, D0 + in.pguess * in.s, lam, ex);
	const double w1 = in.Z[1] / (in.Z[0] + in.Z[1]), w2 = 1.0 - w1;
	auto solve = [&](double R, double& p1, double& p2) -> bool {
		const double a1 = ak[0] + R, a2 = ak[1] - R;
		const double C1 = dedp_am(0, a1, in.mtr[0]) + om[0] * in.s, C2 = dedp_am(1, a2, in.mtr[1]) + om[1] * in.s;
		const double M11 = C1 + R * w1 - lam[0] * in.s * om[0], M12 = R * w2 - lam[0] * in.s * om[1];
		const double M21 = -R * w1 - lam[1] * in.s * om[0], M22 = C2 - R * w2 - lam[1] * in.s * om[1];
		const double r1 = in.etr[0] - in.mtr[0] * PH[0].q - PH[0].g * PH[0].pi * dedp_am(0, a1, in.mtr[0]) + lam[0] * D0 + ex[0];
		const double r2 = in.etr[1] - in.mtr[1] * PH[1].q - PH[1].g * PH[1].pi * dedp_am(1, a2, in.mtr[1]) + lam[1] * D0 + ex[1];
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
	const double rho = in.mtr[0] + in.mtr[1], eC = in.rhoeC / rho, tn = in.tau * (in.nu >= 0 ? in.nu : nu_rate);
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
// 原始变量（论文 2.3，eq:primitive）：Q = (alpha1, rho1, rho2, u, v, p1, p2, p)，p 为总能量的混合压力
// sum_k eps_k(alpha_k, m_k, p) = E - 1/2 rho |u|^2 (A2.3)；派生量：rho、冻结声速 cf
struct Prim { double a1, r1, r2, u, v, p1, p2, p, rho, cf; };
// 守恒量 U -> 原始变量（不可容许时返回 false）
static bool prim_of(const double* Q, Prim& w) {
	const double a1 = Q[IA1], a2 = 1 - a1, m1 = Q[IM1], m2 = Q[IM2];
	if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return false;
	w.a1 = a1; w.rho = m1 + m2; w.r1 = m1 / a1; w.r2 = m2 / a2; w.u = Q[IMU] / w.rho; w.v = Q[IMV] / w.rho;
	w.p = p_mix(a1, m1, m2, Q[IEN] - 0.5 * w.rho * (w.u * w.u + w.v * w.v));
	w.p1 = p_ame(0, a1, m1, Q[IE1]); w.p2 = p_ame(1, a2, m2, Q[IE2]);
	if (!(p_admissible(0, a1, m1, w.p1) && p_admissible(1, a2, m2, w.p2))) return false;
	w.cf = std::sqrt((m1 * c2_rp(0, w.r1, w.p1) + m2 * c2_rp(1, w.r2, w.p2)) / w.rho);
	return true;
}
// 原始变量 -> 守恒量 U：相能量由 (alpha_k, m_k, p_k)，E = sum_k eps_k + 1/2 rho |u|^2，于是 p = p_mix 与 E 一致
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
	testBartak = 32,              // Pelanti 2022 8.2.4 / Sirianni-Re 2026 5.2：Barták 管道快速泄压（NASG 水 350-550 K，有限速率传质 nulaw=1）
	testRelax0D = 40,             // 0D 松弛：均匀非平衡态（周期）
	testRelax0DChem = 41,         // 0D 化学松弛：pT 平衡、过热（周期）
	testVortexBB = 56,            // Battisti-Boscheri 4.1 两相等熵涡（Hu-Shu 型，周期，平移一周回到初值：精确解就是初值）；背景压力 p0 = VORTEX_P0（命令行 p0=）
	testTaylorGreen1 = 57,        // Battisti-Boscheri 4.4 低马赫 Taylor-Green 涡，单相（alpha_s = 0.5）：不可压 Euler 的定常解；p0 = VORTEX_P0，M = 1/sqrt(1.4 p0)
	testTaylorGreen2 = 58,        //   两相（alpha_s = 0.5 + 0.1 (cos 2x + cos 2y)）
	testSPB_Nitromethane = 60,    // SPB 2009 4.2  硝基甲烷激波管（Cochran-Chan，同一流体当两相），20 GPa / 0.2 MPa
	testPel_JWLWater = 61,        // Pelanti 2022 8.1.1  爆轰气体（JWL）-水（NASG）激波管
	testSPB_HeBubble = 70,        // SPB 2009 4.6  二维：空气中 Mach 1.5 激波打氦气泡（Layes & Le Métayer）
	testSPB_RM = 71,              // SPB 2009 4.7  二维：空化 Richtmyer-Meshkov（水 | 气，-200 m/s 撞左壁）
	testPS_BubbleComp = 72,       // P-S 2014 5.4  二维：活塞压缩十二烷汽泡（SG，热、化学松弛）
	testPel_UNDEX = 73,           // Pelanti 2022 8.1.2  二维：近壁水下爆炸（JWL 气泡 | NASG 水）
	// Schropff-Petitpas-Daniel 2026（JCP 547, 114545）5.1、5.1.1、5.2：[0, 1] m，透射边界。单相算例两相取同一流体、同一状态（alpha1 = 1 - 1e-8）
	testSchAirAdv = 80,           // 表 2  低马赫空气平移：rho 1000 | 0.01，u = 1，p = 1e5，界面 0.25，t = 0.5 s
	testSchWaterLM = 81,          // 表 2  低马赫水激波管：rho 1000，u 0 | 15，p 1e8 | 0.98e8，t = 1e-4 s
	testSchLax = 82,              // 表 2  Lax 激波管（空气）：(0.445, 0.698, 3.528) | (0.5, 0, 0.571)，t = 0.14
	testSchStrong = 83,           // 表 2  强激波管（空气）：(1, 0, 1e10) | (0.125, 0, 0.1)，t = 2.5e-6 s
	testSchTwoRare = 84,          // 表 2  双稀疏波（空气）：rho 1，u -1 | 1，p 0.4，t = 0.15
	testSchLiqLiq = 85,           // 5.1.1  液-液激波管：水 rho 1000，p 1e6 | 1e5，t = 0.15 ms
	testSchWaterAirAdv = 86,      // 表 3  水(1)-气(2) 平移：alpha1 1 - 1e-7 | 1e-7，rho1 1000，rho2 0.01，u = 1，p = 1e5，界面 0.25，t = 0.25 s
	testSchWaterAirST = 87,       // 表 3  水-气激波管：rho1 1000，rho2 50，p 1e9 | 1e5，界面 0.7，t = 241 us
	testSchEpoxySpinel = 88       // 表 3  环氧-尖晶石：alpha1 = 0.5954，rho 1185 / 3622，p 1e10 | 1e5，界面 0.6，t = 30 us
};
struct TestSetup { double x0, x1, y0, y1, tend; x_boundary_condition xbc; y_boundary_condition ybc;
	int bc[4] = {-1, -1, -1, -1};      // 每一边（左、右、下、上）单独给：BC_TRANS / BC_WALL；-1 表示用 xbc / ybc
	double bcu[4] = {0, 0, 0, 0};      // 固壁的法向速度（沿 +x / +y 方向；活塞）
	double bcp[4] = {0, 0, 0, 0};      // BC_PRES 的出口压力
	int sub = 1;                       // 初值的子格采样数（每方向）：圆形区域的格平均
	int nx = -1, ny = -1, relax = -1;  // 论文的网格与松弛档位（命令行 N= ny= relax= 覆盖）
	double cfl = -1;                   // IMEX 的物质 CFL（-1：程序默认 0.2）。声波本身是要解的结构时取小（B10，devlog V14）
	int rus = -1; };                   // IMEX 显式 Rusanov 的耗散速度：1 每个面局部，0 全场最大（-1：默认 1）（B3，devlog V14）
enum { BC_TRANS = 0, BC_WALL = 1, BC_PERIODIC = 2, BC_PRES = 3 };
// BC_PRES（只用于 x 两边）：给定压力出口，虚格 = 最近内部格的 alpha、m、速度，两相压力都取 BCP[side]（Pelanti 8.2.4 的大气压出口）；隐式部分的外侧压力取 BCP（不是未知量）
int BCS[4] = {0, 0, 0, 0}; double BCU[4] = {0, 0, 0, 0}; double BCP[4] = {0, 0, 0, 0}; int SUBCELL = 1;
// 输出（两个程序共用）：snap=t1,t2,... 在这些时刻写剖面（时间步截到这些时刻）；probe=x,y 每步记录最近格的压力；prec= 剖面的有效位数
double VORTEX_P0 = 1.0;   // 算例 56 的背景压力 p0（Battisti-Boscheri 表 1：1, 10, 1e2, 1e3, 1e4；命令行 p0=）
int WELLPREP = 0;         // 算例 56：1 = 离散良态初值（速度由格心流函数的中心差分给出，离散散度为零；命令行 wellprep=1，devlog V15）
int INIT_ONLY = 0;        // 1 = 只写精确格平均初值后退出（命令行 init=1；误差后处理用）
int LIMITER = 1;   // 1: minmod（所有算例）；0: 不限制的中心斜率，只用于光滑解的收敛阶试验（限制器在光滑极值处降为一阶，devlog V14）
std::vector<double> SNAP; size_t SNAPI = 0; bool HAS_PROBE = false; double PROBE_X = 0, PROBE_Y = 0; int OUTPREC = 16;
static void parse_list(const std::string& v, std::vector<double>& out) { size_t a = 0; while (a < v.size()) { size_t b = v.find(',', a); if (b == std::string::npos) b = v.size(); out.push_back(std::atof(v.substr(a, b - a).c_str())); a = b + 1; } }
inline double snap_limit(double t, double dt) { if (SNAPI < SNAP.size() && t + dt > SNAP[SNAPI]) return SNAP[SNAPI] - t; return dt; }
inline std::string snap_name(const std::string& out, size_t k) { std::string b = out; if (b.size() > 4 && b.substr(b.size() - 4) == ".dat") b = b.substr(0, b.size() - 4); return b + "_s" + std::to_string(k) + ".dat"; }
static void apply_bc_setup(const TestSetup& S) {
	const int dx_ = S.xbc == x_Periodic ? BC_PERIODIC : (S.xbc == x_Reflective ? BC_WALL : BC_TRANS);
	const int dy_ = S.ybc == y_Periodic ? BC_PERIODIC : (S.ybc == y_Reflective ? BC_WALL : BC_TRANS);
	for (int q = 0; q < 4; q++) { BCS[q] = S.bc[q] >= 0 ? S.bc[q] : (q < 2 ? dx_ : dy_); BCU[q] = S.bcu[q]; BCP[q] = S.bcp[q]; }
	SUBCELL = S.sub;
}
// 边界条件（两个程序共用）：透射复制边界格；固壁镜像，法向速度取 2 u_w - u（u_w = 0 时就是取负；总能随动能改）；周期绕回。
// side: 0 左 1 右 2 下 3 上。一维（ny = 1 且 y 非周期）时 y 向只复制
template <class GridT> static void bc_apply(GridT& u, int nx, int ny, int g) {
	auto ghost = [&](int i, int j, int si, int sj, int side) {
		for (int k = 0; k < 8; k++) u(i, j, k) = u(si, sj, k);
		if (BCS[side] == BC_PRES) {   // 两相压力取出口压力，alpha、m、速度照搬
			const double a1 = u(i, j, IA1), e1 = eps_amp(0, a1, u(i, j, IM1), BCP[side]), e2 = eps_amp(1, 1 - a1, u(i, j, IM2), BCP[side]);
			u(i, j, IEN) += (e1 + e2) - (u(i, j, IE1) + u(i, j, IE2)); u(i, j, IE1) = e1; u(i, j, IE2) = e2; return; }
		if (BCS[side] != BC_WALL) return;
		const int kn = side < 2 ? IMU : IMV; const double rho = u(i, j, IM1) + u(i, j, IM2);
		if (BCU[side] == 0.0) { u(i, j, kn) = -u(i, j, kn); return; }
		const double mn = 2.0 * rho * BCU[side] - u(si, sj, kn);
		u(i, j, IEN) += 0.5 * (mn * mn - u(si, sj, kn) * u(si, sj, kn)) / rho; u(i, j, kn) = mn; };
	for (int j = 0; j < ny + 2 * g; j++)
		for (int q = 0; q < g; q++) {
			const int iL = g - 1 - q, iR = nx + g + q;
			const int sL = BCS[0] == BC_PERIODIC ? nx + g - 1 - q : (BCS[0] == BC_WALL ? g + q : g);
			const int sR = BCS[1] == BC_PERIODIC ? g + q : (BCS[1] == BC_WALL ? nx + g - 1 - q : nx + g - 1);
			ghost(iL, j, sL, j, 0); ghost(iR, j, sR, j, 1);
		}
	if (ny > 1 || BCS[2] == BC_PERIODIC) {
		for (int i = 0; i < nx + 2 * g; i++)
			for (int q = 0; q < g; q++) {
				const int jL = g - 1 - q, jR = ny + g + q;
				const int sL = BCS[2] == BC_PERIODIC ? ny + g - 1 - q : (BCS[2] == BC_WALL ? g + q : g);
				const int sR = BCS[3] == BC_PERIODIC ? g + q : (BCS[3] == BC_WALL ? ny + g - 1 - q : ny + g - 1);
				ghost(i, jL, i, sL, 2); ghost(i, jR, i, sR, 3);
			}
	} else {
		for (int i = 0; i < nx + 2 * g; i++)
			for (int q = 0; q < g; q++)
				for (int k = 0; k < 8; k++) { u(i, g - 1 - q, k) = u(i, g, k); u(i, g + 1 + q, k) = u(i, g, k); }
	}
}
static void set_SG(int k, double g, double pi, double q = 0, double qp = 0, double cv = 1000.0, double b = 0) { PH[k] = {g, pi, q, qp, cv, b}; }
static double rho_pT(int k, double p, double T) { return 1.0 / v_pT(k, p, T); }
// 状态方程、区域、边界、终止时间
static bool test_setup(int test, TestSetup& S) {
	S.x0 = 0; S.x1 = 1; S.y0 = 0; S.y1 = 1; S.xbc = x_Transmissive; S.ybc = y_Transmissive;
	switch (test) {
	case testSPB_Cavitation: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 1.85e-3; S.nx = 1000; break;
	case testPS_Cavitation2: case testPS_Cavitation500:
		set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0);
		S.tend = test == testPS_Cavitation2 ? 3.2e-3 : 0.58e-3; S.nx = 5000; break;
	case testDodecaneSG: set_SG(0, 2.35, 4e8, -775.269e3, 0.0, 1077.7); set_SG(1, 1.025, 0.0, -237.547e3, -24.4e3, 1956.45); S.tend = 473e-6; S.nx = 2000;
		S.cfl = 0.05; break;   // 液体中 1e8 Pa 的稀疏波是要解的结构：物质 CFL 0.2 时其声学 Courant 数约 0.9，LSDIRK2 的色散使平台压力偏低（devlog V14）
	case testSPB_Interface: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 2.79e-3; S.nx = 200; break;
	case testSPB_WaterAir: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 240e-6; S.nx = 1000; break;
	case testSPB_WaterAirTPa: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 8.3e-6; S.nx = 1000; break;
	case testSPB_EpoxySpinel: case testSPB_EpoxySpinelExt: set_SG(0, 2.43, 5.3e9); set_SG(1, 1.62, 141e9); S.tend = test == testSPB_EpoxySpinel ? 80e-6 : 34e-6; S.nx = 500;
		if (test == testSPB_EpoxySpinel) S.cfl = 0.05;   // 5.1：激波与稀疏波都是声学结构（u << c），同上（devlog V14）
		else S.rus = 0;                                  // 5.2：全场最大的 Rusanov 速度更好（L1(u) 118 对 124，devlog V14）
		break;
	case testSmoothGas: set_SG(0, 1.4, 0); set_SG(1, 1.4, 0); S.tend = 5e-4; break;
	case testSmoothWater: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = 1e-4; break;
	case testPel_Cavitation: PH[0] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4}; PH[1] = {1.467, 0.0, 2077616.0, 14317.0, 955.0, 0.0}; S.tend = 3e-3; S.nx = 5000; break;
	case testDodecaneNASG: PH[0] = {1.206, 1681e5, -996054.0, 0.0, 2532.0, 7.51e-4}; PH[1] = {1.021, 0.0, -384592.0, -4301.0, 2274.0, 0.0}; S.tend = 473e-6; S.nx = 2000; break;
	case testBartak:   // Pelanti 8.2.4：表 5（NASG 水 350-550 K）；管长 1.7 m，右端封闭，左端突然打开：左边界给定大气压 1e5 Pa；1000 格，t = 17.5 ms
		PH[0] = {1.387, 8899e5, -1244191.0, 0.0, 3202.0, 4.78e-4}; PH[1] = {1.954, 0.0, 2287484.0, 6417.0, 462.0, 0.0};
		S.x0 = 0.0; S.x1 = 1.7; S.tend = 17.5e-3; S.nx = 1000; S.relax = 2; S.bc[0] = BC_PRES; S.bcp[0] = 1e5; S.bc[1] = BC_WALL; break;
	case testRelax0D: case testRelax0DChem: set_SG(0, 2.35, 1e9, -1167e3, 0.0, 1816.0); set_SG(1, 1.43, 0.0, 2030e3, -23.4e3, 1040.0); S.tend = 1e-4; S.xbc = x_Periodic; break;
	case testVortexBB:   // Battisti-Boscheri 4.1：[-5, 5]^2 周期，gamma_f = 1.35（相 1）、gamma_s = 1.4（相 2），平移 (2, 2)，t = 5 正好一周
		set_SG(0, 1.35, 0); set_SG(1, 1.4, 0); S.x0 = -5; S.x1 = 5; S.y0 = -5; S.y1 = 5; S.tend = 5.0; S.xbc = x_Periodic; S.ybc = y_Periodic;
		S.sub = 4; S.nx = 40; S.ny = 40; S.relax = 0; break;   // 初值取 4 x 4 子点的格平均（误差按格平均比较）
	case testTaylorGreen1: case testTaylorGreen2:   // BB 4.4：[0, 2 pi]^2 周期，50 x 50，两相同为理想气体 gamma = 1.4、rho = 1；t = 0.2，固定 dt = 2.5e-3（命令行 dt=）
		set_SG(0, 1.4, 0); set_SG(1, 1.4, 0); S.x0 = 0; S.x1 = 2 * M_PI; S.y0 = 0; S.y1 = 2 * M_PI; S.tend = 0.2; S.xbc = x_Periodic; S.ybc = y_Periodic;
		S.sub = 4; S.nx = 50; S.ny = 50; S.relax = 0; break;
	case testSPB_Nitromethane:   // SPB 4.2：Gamma = 1.19，rho_ref = 1134，A1 = 0.819181e9，A2 = 1.50835e9，E1 = 4.52969，E2 = 1.42144；界面位置论文未写，取 0.5 m
		for (int k = 0; k < 2; k++) { PH[k] = Phase{}; PH[k].kind = 1; PH[k].G = 1.19; PH[k].rr = 1134.0; PH[k].c1 = 0.819181e9; PH[k].c2 = 1.50835e9; PH[k].e1 = 4.52969; PH[k].e2 = 1.42144; PH[k].g = 2.19; }
		S.tend = 67e-6; S.nx = 500; break;
	case testPel_JWLWater:       // Pelanti 表 2（JWL）与表 3（NASG 液态水）；[-10, 10] m，t = 0.9 ms
		PH[0] = Phase{}; PH[0].kind = 2; PH[0].G = 0.25; PH[0].rr = 1840.0; PH[0].c1 = 854.5e9; PH[0].c2 = 20.5e9; PH[0].e1 = 4.6; PH[0].e2 = 1.35; PH[0].e0 = 0.0; PH[0].g = 1.25;
		PH[1] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4};
		S.x0 = -10.0; S.x1 = 10.0; S.tend = 0.9e-3; S.nx = 500; break;
	case testSPB_HeBubble:   // 区域 0.3 x 0.08 m，气泡直径 0.04 m（论文）；气泡中心 (0.2, 0.04)、激波初位置 x = 0.221 m（论文未给，取 t = 0 时激波刚到气泡右缘）
		set_SG(0, 1.4, 0); set_SG(1, 1.67, 0); S.x0 = 0; S.x1 = 0.3; S.y0 = 0; S.y1 = 0.08; S.tend = 600e-6;
		S.xbc = x_Transmissive; S.ybc = y_Reflective; S.sub = 4; S.nx = 600; S.ny = 160; break;
	case testSPB_RM:         // 3 m x 1 m；左、上、下为固壁，右为出流；水 SG (4.4, 6e8)，气 gamma = 1.8（SPB；P-S 用 1.4）
		set_SG(0, 4.4, 6e8); set_SG(1, 1.8, 0); S.x0 = 0; S.x1 = 3; S.y0 = 0; S.y1 = 1; S.tend = 8.6e-3;
		S.xbc = x_Transmissive; S.ybc = y_Reflective; S.bc[0] = BC_WALL; S.sub = 4; S.nx = 300; S.ny = 100;   // P-S 5.3 的网格（SPB 用 900 x 400：N=900 ny=400）
		S.rus = 0; break;   // 全场最大 Rusanov 速度（局部时 18842 次拒绝，devlog V14）
	case testPS_BubbleComp:  // 单位正方形，汽泡半径 0.2 m 居中；四周固壁，左壁为活塞 100 m/s（网格不动：活塞面固定在 x = 0，液体从该面进入）
		set_SG(0, 2.35, 4e8, -775.269e3, 0.0, 1077.7); set_SG(1, 1.025, 0.0, -237.547e3, -24.4e3, 1956.45);
		S.tend = 1.2e-3; S.xbc = x_Reflective; S.ybc = y_Reflective; S.bcu[0] = 100.0; S.sub = 4; S.nx = 200; S.ny = 200; S.relax = 2;
		S.rus = 0; break;   // 全场最大 Rusanov 速度（局部时界面旁微量汽相质量为负的舍入，devlog V14）
	case testPel_UNDEX:      // [-1.35, 1.35] x [-1.5, 0] m，y = 0 为固壁，其余透射；气泡中心 (0, -0.32) m，半径 0.05 m
		PH[0] = Phase{}; PH[0].kind = 2; PH[0].G = 0.25; PH[0].rr = 1840.0; PH[0].c1 = 854.5e9; PH[0].c2 = 20.5e9; PH[0].e1 = 4.6; PH[0].e2 = 1.35; PH[0].e0 = 0.0; PH[0].g = 1.25;
		PH[1] = {1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4};
		S.x0 = -1.35; S.x1 = 1.35; S.y0 = -1.5; S.y1 = 0.0; S.tend = 0.4e-3; S.xbc = x_Transmissive; S.ybc = y_Transmissive; S.bc[3] = BC_WALL; S.sub = 4; S.nx = 361; S.ny = 200; break;   // 论文网格 721 x 400（N=721 ny=400）IMEX 约 14 h；默认用一半分辨率（devlog V14）
	case testSchAirAdv: case testSchLax: case testSchStrong: case testSchTwoRare:
		set_SG(0, 1.4, 0); set_SG(1, 1.4, 0);
		S.tend = test == testSchAirAdv ? 0.5 : test == testSchLax ? 0.14 : test == testSchStrong ? 2.5e-6 : 0.15; S.nx = 1000; break;
	case testSchWaterLM: case testSchLiqLiq: set_SG(0, 4.4, 6e8); set_SG(1, 4.4, 6e8); S.tend = test == testSchWaterLM ? 1e-4 : 0.15e-3; S.nx = 1000; break;
	case testSchWaterAirAdv: case testSchWaterAirST: set_SG(0, 4.4, 6e8); set_SG(1, 1.4, 0); S.tend = test == testSchWaterAirAdv ? 0.25 : 241e-6; S.nx = 1000; break;
	case testSchEpoxySpinel: set_SG(0, 2.43, 5.3e9); set_SG(1, 1.62, 141e9); S.tend = 30e-6; S.nx = 1000; break;
	default: return false;
	}
	return true;
}
// 点 (x, y) 处的初值（守恒量）
// 算例 56 的流函数（扰动部分）：u = 2 + d psi/dy，v = 2 - d psi/dx，psi = eps/(2 pi) exp((1 - r^2)/2)
inline double vortex_psi(double x, double y) { const double eps = 5.0; return eps / (2 * M_PI) * std::exp(0.5 * (1 - x * x - y * y)); }
// 算例 57、58 的流函数：u = d psi/dy = sin x cos y，v = -d psi/dx = -cos x sin y
inline double tg_psi(double x, double y) { return std::sin(x) * std::sin(y); }
static void test_state(int test, double x, double y, double* Q) {
	(void)y;
	auto setcell = [&](double a1, double r1, double r2, double u, double v, double p1, double p2) {
		Prim w; w.a1 = a1; w.r1 = r1; w.r2 = r2; w.u = u; w.v = v; w.p1 = p1; w.p2 = p2; prim2cell(w, Q); };
	switch (test) {
	case testSPB_Cavitation: { const double u = x < 0.5 ? -100.0 : 100.0; setcell(1 - 1e-2, 1000.0, 1.0, u, 0, 1e5, 1e5); break; }   // 空气密度论文未给，取 1
	case testPS_Cavitation2: case testPS_Cavitation500: {
		const double T0 = 354.728, U0 = test == testPS_Cavitation2 ? 2.0 : 500.0, u = x < 0.5 ? -U0 : U0;
		setcell(1 - 1e-2, rho_pT(0, 1e5, T0), rho_pT(1, 1e5, T0), u, 0, 1e5, 1e5); break; }
	case testBartak: {   // 液态水 12.5 MPa、563.15 K，均匀的汽 alpha_v = 1e-3（两相同温）
		const double p = 12.5e6, T = 563.15; setcell(1 - 1e-3, rho_pT(0, p, T), rho_pT(1, p, T), 0, 0, p, p); break; }
	case testDodecaneSG: case testDodecaneNASG: { const bool L = x < 0.75; setcell(L ? 1 - 1e-8 : 1e-8, 500.0, 2.0, 0, 0, L ? 1e8 : 1e5, L ? 1e8 : 1e5); break; }
	case testSPB_Nitromethane: { const bool L = x < 0.5; setcell(L ? 1 - 1e-8 : 1e-8, 1134.0, 1200.0, 0, 0, L ? 20e9 : 0.2e6, L ? 20e9 : 0.2e6); break; }
	case testPel_JWLWater: { const bool L = x < 0.0; setcell(L ? 1 - 1e-8 : 1e-8, 2000.0, 1044.0, 0, 0, L ? 4.6406e10 : 1e5, L ? 4.6406e10 : 1e5); break; }
	case testSPB_HeBubble: {   // 空气 1.29 kg/m3、101325 Pa；氦 0.167 kg/m3；激波右侧为 Mach 1.5 波后状态（理想气体 Rankine-Hugoniot），向左运动
		const double g = 1.4, Ms = 1.5, p0 = 101325.0, r0 = 1.29, c0 = std::sqrt(g * p0 / r0);
		const double p1 = p0 * (1 + 2 * g / (g + 1) * (Ms * Ms - 1)), r1 = r0 * (g + 1) * Ms * Ms / ((g - 1) * Ms * Ms + 2), u1 = -Ms * c0 * (1 - r0 / r1);
		const bool in = (x - 0.2) * (x - 0.2) + (y - 0.04) * (y - 0.04) < 0.02 * 0.02;
		if (x > 0.221) setcell(1 - 1e-6, r1, 0.167, u1, 0, p1, p1);
		else setcell(in ? 1e-6 : 1 - 1e-6, r0, 0.167, 0, 0, p0, p0);
		break; }
	case testSPB_RM: {         // 水在左：x < 1.2 + sqrt(0.36 - (y - 0.5)^2)（半径 0.6 m、圆心 (1.2, 0.5) 的圆弧凸向气体）；全场 u = -200 m/s，p = 1 atm
		const double dy_ = y - 0.5, xi = 1.2 + std::sqrt(std::max(0.0, 0.36 - dy_ * dy_)); const bool W = x < xi;
		setcell(W ? 1 - 1e-6 : 1e-6, 1000.0, 100.0, -200.0, 0, 101325.0, 101325.0); break; }
	case testPS_BubbleComp: {  // p = 1e5 Pa、T = 600 K：rho_l = 458.338、rho_v = 3.408；汽泡内 alpha_l = 1e-6，液体中 alpha_v = 1e-6
		const bool in = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.04;
		setcell(in ? 1e-6 : 1 - 1e-6, rho_pT(0, 1e5, 600.0), rho_pT(1, 1e5, 600.0), 0, 0, 1e5, 1e5); break; }
	case testPel_UNDEX: {      // 与 8.1.1 同样的两个状态；气泡内 alpha_gas = 1 - 1e-6
		const bool in = x * x + (y + 0.32) * (y + 0.32) < 0.05 * 0.05;
		setcell(in ? 1 - 1e-6 : 1e-6, 2000.0, 1044.0, 0, 0, in ? 4.6406e10 : 1e5, in ? 4.6406e10 : 1e5); break; }
	case testSchAirAdv: case testSchWaterLM: case testSchLax: case testSchStrong: case testSchTwoRare: case testSchLiqLiq: {   // 单相：两相同一状态
		double xi = 0.5, rL = 0, uL = 0, pL = 0, rR = 0, uR = 0, pR = 0;
		switch (test) {
		case testSchAirAdv: xi = 0.25; rL = 1000; uL = 1; pL = 1e5; rR = 0.01; uR = 1; pR = 1e5; break;
		case testSchWaterLM: rL = 1000; uL = 0; pL = 1e8; rR = 1000; uR = 15; pR = 0.98e8; break;
		case testSchLax: rL = 0.445; uL = 0.698; pL = 3.528; rR = 0.5; uR = 0; pR = 0.571; break;
		case testSchStrong: rL = 1; uL = 0; pL = 1e10; rR = 0.125; uR = 0; pR = 0.1; break;
		case testSchTwoRare: rL = 1; uL = -1; pL = 0.4; rR = 1; uR = 1; pR = 0.4; break;
		default: rL = 1000; uL = 0; pL = 1e6; rR = 1000; uR = 0; pR = 1e5; break; }
		const bool L = x < xi; const double r = L ? rL : rR, u = L ? uL : uR, p = L ? pL : pR;
		setcell(1 - 1e-8, r, r, u, 0, p, p); break; }
	case testSchWaterAirAdv: case testSchWaterAirST: {
		const bool adv = test == testSchWaterAirAdv, L = x < (adv ? 0.25 : 0.7);
		setcell(L ? 1 - 1e-7 : 1e-7, 1000.0, adv ? 0.01 : 50.0, adv ? 1.0 : 0.0, 0, adv ? 1e5 : (L ? 1e9 : 1e5), adv ? 1e5 : (L ? 1e9 : 1e5)); break; }
	case testSchEpoxySpinel: { const double p = x < 0.6 ? 1e10 : 1e5; setcell(0.5954, 1185.0, 3622.0, 0, 0, p, p); break; }
	case testSPB_Interface: { const bool L = x < 0.5; setcell(L ? 1 - 1e-8 : 1e-8, 1000.0, 10.0, 100.0, 0, 1e5, 1e5); break; }
	case testSPB_WaterAir: case testSPB_WaterAirTPa: {
		// 4.3.1 界面 0.75 m（正文）；4.3.2 正文说“同 4.3.1”，但论文图 10 的波位置（稀疏波头 0.035、波尾 0.765、接触 0.857、激波 0.910，t = 8.3 us）
		// 反推出界面在 0.59-0.60 m，取 0.6 m（devlog V15 §15.8）
		const bool L = x < (test == testSPB_WaterAir ? 0.75 : 0.6); const double pL = test == testSPB_WaterAir ? 1e9 : 1e12, ra = test == testSPB_WaterAir ? 1.0 : 10.0;
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
	case testTaylorGreen1: case testTaylorGreen2: {   // BB (73)：rho = 1，u = sin x cos y，v = -cos x sin y，p = p0 + (cos 2x + cos 2y)/4，两相同压
		const double p = VORTEX_P0 + 0.25 * (std::cos(2 * x) + std::cos(2 * y));
		const double a2 = test == testTaylorGreen1 ? 0.5 : 0.5 + 0.1 * (std::cos(2 * x) + std::cos(2 * y));
		setcell(1 - a2, 1.0, 1.0, std::sin(x) * std::cos(y), -std::cos(x) * std::sin(y), p, p); break; }
	case testVortexBB: {   // Hu-Shu 涡（eps = 5），T = p0 + dT，rho = (T / p0)^(1/(gamma-1))，p = rho T（gamma = 1.4）：dp/dr = rho u_theta^2 / r 与 p0 无关，
		// 速度扰动与 p0 无关，密度扰动 ~ 1/p0（与论文表 1 中 L2(rho_s phi_s) 随 p0 减小、L2(动量)、L2(phi_s) 不变一致）。
		// 单速度模型要求两相同压（否则 (p1 - p2) grad alpha 破坏径向平衡）：两相取同一 rho、p，alpha2 为被动输运的光滑分布。
		// div u = 0、alpha 与熵沿流线不变，故平移解精确（devlog V15）。Dumbser 等 [67] 的 BN 原设置（两速度、两压力）在单速度模型里不可用。
		const double g = 1.4, eps = 5.0, p0 = VORTEX_P0, r2 = x * x + y * y;
		const double du = eps / (2 * M_PI) * std::exp(0.5 * (1 - r2)), dT = -(g - 1) * eps * eps / (8 * g * M_PI * M_PI) * std::exp(1 - r2);
		const double T = p0 + dT, rho = std::pow(T / p0, 1.0 / (g - 1)), p = rho * T, a2 = 0.5 + 0.25 * std::exp(0.5 * (1 - r2));
		setcell(1 - a2, rho, rho, 2.0 - du * y, 2.0 + du * x, p, p); break; }
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
