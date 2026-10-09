// =====================================================================================================================
//  Saurel-Petitpas-Berry (JCP 228 (2009) 1678) 6+1 方程松弛模型, 一维.
//    MODE 0: 论文原方法 (3.1.2 HLLC + 3.2.1 Godunov, 显式, 声学 CFL)            —— 参考解
//    MODE 1: IMEX (Toro-Vazquez / Battisti-Boscheri 型分裂):
//            显式 = 对流子系统, path-conservative Rusanov, 耗散速度 = 全局最大物质速度 max|u|
//            隐式 = 压力子系统, 混合压力的线性三对角方程 (动能 Picard 迭代)
//    两种模式共用论文的 3.3 松弛步 (III.4) 和 3.4 重初始化 (II.4)(III.5).
//  状态: alpha1, m1 = (alpha rho)_1, m2, rho u, (alpha rho e)_1, (alpha rho e)_2, rho E   (论文 (III.1) 的 7 个方程)
//  编译: g++ -O2 -std=c++17 spb_imex.cpp -o spb   运行: ./spb MODE N CFL [out.dat] [预扫遍数=20] [阶数=1|2]
// =====================================================================================================================
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

static double GAM[2] = {4.4, 1.4};        // 相 1 = 液 (stiffened gas), 相 2 = 气 (理想气体); 默认 = SPB 4.5 的水/空气
static double PIN[2] = {6.0e8, 0.0};
static double QQ[2] = {0.0, 0.0}, QP[2] = {0.0, 0.0}, CV[2] = {1000.0, 718.0};   // 生成焓 q, 熵常数 q', 定容比热 (热/化学松弛用)
static int CENTRAL = 0;   // 隐式增量通量里的面值 (相内能, 质量): 0 = 按 ubar*_f 迎风 (质量与内能同侧), 1 = 算术平均. 尖锐界面 (密度比 ~250) 上平均值不可用
static int RELAX = 0;   // 0: 只力学松弛 (SPB 原文); 1: + 热松弛 (pT); 2: + 化学松弛 (pTg, 亚稳液体蒸发)
static const int NV = 7, NG = 2;
struct Cell { double a1, m[2], mu, ei[2], E; int eq = 0; };   // eq: 上一次松弛里化学松弛 (pTg) 在本格生效  // ei[k] = (alpha rho e)_k,  E = rho E

static inline double alpha(const Cell& c, int k){ return k==0 ? c.a1 : 1.0 - c.a1; }
static inline double pk(const Cell& c, int k){ return (GAM[k]-1.0)*c.ei[k]/alpha(c,k) - GAM[k]*PIN[k]; }           // (II.3)
static inline double pmix(const Cell& c){ return alpha(c,0)*pk(c,0) + alpha(c,1)*pk(c,1); }
static inline double ck2(const Cell& c, int k){ double r = c.m[k]/alpha(c,k); return GAM[k]*std::max(pk(c,k)+PIN[k], 1e-6)/r; }
static inline double cfrozen(const Cell& c){ double r = c.m[0]+c.m[1]; return std::sqrt((c.m[0]*ck2(c,0) + c.m[1]*ck2(c,1))/r); }   // (III.2)
// 混合 EOS (II.4): rho e = p sum(alpha_k/(g_k-1)) + sum(alpha_k g_k pi_k/(g_k-1))
static inline double Acoef(double a1){ return a1/(GAM[0]-1.0) + (1.0-a1)/(GAM[1]-1.0); }
static inline double Bcoef(double a1){ return a1*GAM[0]*PIN[0]/(GAM[0]-1.0) + (1.0-a1)*GAM[1]*PIN[1]/(GAM[1]-1.0); }

static long g_relax_fail = 0, g_pfloor = 0, g_newton_steps = 0, g_newton_fail = 0;
static const double ALPHA_EPS = 1e-10;
static int SWEEPS = 20;
static int WALLL = 0; static int RANN = 1, NRAN = 0; static int UREF = 1;   /* 2026-10-07: 默认改为 1. UREF = 2 漏掉显式动量对流引起的面速度变化, 时间上只有一阶 (推导文档 12.1) */ static int QEGLOB = 1; static long g_trbad = 0; static int BNDCHK = 1;   // 二阶: 第二级的输运若给出负的相体积或质量, 本步退回一阶 (显式表的保界条件约为 CFL <= 0.19)
static int FLAT = 0, FLATA = 0; static double FK1 = 0.4, FSC = 0.0; static long g_flatcells = 0;   /* FLAT: 激波处的隐式耗散; FK1: 指示器里的 kappa_1 (Balsara 取 0.4) */
static int SLOC = 0, DONOR = 1;   // DONOR 默认开 (2026-10-07): 一维全套算例上 Newton 不收敛归零, 误差持平; SLOC 在空化中心变差, 默认关
static int AVOL = 1;    // 1: 相体积 alpha_k 作为被输运量, 与质量/相内能用同一个通量算子; 力学闭合 N_k = eps^tr_k + p alpha^tr_k (不再单列压力功项). 0: 旧形式 (path-conservative u d(alpha) + alpha* p d(u))
static int UNEW = 1;
static int PFACE = 0;   // 动量/总能通量里的面压力: 1 = 按 1/rho 加权 (加速度连续), 0 = 算术平均
static int IRF = 1;     // 面上的 1/rho: 1 = 1/(rho 的算术平均), 0 = (1/rho) 的算术平均
static int EXPLE = 1;   // 1: 相内能的对流通量 (a rho e)_k u 放在显式半步 (与质量同一个通量形式), 隐式只留增量 (a rho e)_k,f du_f 和压力功

// ---- 3.3 松弛步 (III.4) + 3.4 重初始化 (II.4)(III.5) ----
static void relax_reinit(Cell& c){
	double a[2] = {alpha(c,0), alpha(c,1)}, p0[2] = {pk(c,0), pk(c,1)}, Z[2];
	for(int k=0;k<2;k++) Z[k] = std::sqrt(c.m[k]/a[k]*GAM[k]*std::max(p0[k]+PIN[k], 1e-6)*1.0) * std::sqrt(1.0);   // Z_k = rho_k c_k
	double pI = (Z[1]*p0[0] + Z[0]*p0[1])/(Z[0]+Z[1]);                    // 界面压力 p_I, 取初值作 \hat p_I (论文 4.4: 两种估计无差别)
	// v_k(p) = v_k0 (p_k0 + g pi + (g-1) pI)/(p + g pi + (g-1) pI);  f(p) = sum alpha_k0 * 该比值 - 1, 对 p 单调减
	double num[2], off[2]; bool ok = true;
	for(int k=0;k<2;k++){ off[k] = GAM[k]*PIN[k] + (GAM[k]-1.0)*pI; num[k] = p0[k] + off[k]; if(!(num[k] > 0.0)) ok = false; }
	if(!ok){   // p_I 使某相的分子非正 (液体受拉且 p_I < 0): 改用 \hat p_I = 0 (仍与熵不等式相容的端点之外的保守取法), 再不行则跳过
		pI = 0.0; ok = true;
		for(int k=0;k<2;k++){ off[k] = GAM[k]*PIN[k]; num[k] = p0[k] + off[k]; if(!(num[k] > 0.0)) ok = false; }
	}
	if(getenv("PIHAT_FINAL")){   // \hat p_I = p (松弛后的压力): v_k/v_k0 = (p_k0 + g pi + (g-1) p)/(g (p + pi))
		auto f2 = [&](double p){ double g=-1.0; for(int k=0;k<2;k++) g += a[k]*(p0[k]+GAM[k]*PIN[k]+(GAM[k]-1.0)*p)/(GAM[k]*(p+PIN[k])); return g; };
		double lo = 1e-12, hi = std::max(std::max(p0[0],p0[1]),1.0); int n=0; while(f2(hi) > 0.0 && n++<200) hi*=2.0;
		if(f2(lo) > 0.0){ for(int it=0; it<200 && (hi-lo) > 1e-14*hi; it++){ double m=0.5*(lo+hi); if(f2(m) > 0.0) lo=m; else hi=m; }
			double p=0.5*(lo+hi), a1=a[0]*(p0[0]+GAM[0]*PIN[0]+(GAM[0]-1.0)*p)/(GAM[0]*(p+PIN[0])), a2=a[1]*(p0[1]+GAM[1]*PIN[1]+(GAM[1]-1.0)*p)/(GAM[1]*(p+PIN[1]));
			c.a1 = std::min(std::max(a1/(a1+a2), ALPHA_EPS), 1.0-ALPHA_EPS); } else g_relax_fail++;
	}else
	if(ok){
		auto f = [&](double p){ return a[0]*num[0]/(p+off[0]) + a[1]*num[1]/(p+off[1]) - 1.0; };
		double lo = std::max(-off[0], -off[1]); lo += 1e-12*(std::fabs(lo)+1.0);       // f(lo+) = +inf
		double hi = std::max(std::max(p0[0], p0[1]), lo + 1.0);
		int n = 0; while(f(hi) > 0.0 && n++ < 200) hi = lo + 2.0*(hi - lo) + 1.0;
		double p = std::min(std::max(0.5*(p0[0]+p0[1]), lo), hi);
		for(int it=0; it<200; it++){
			double fp = f(p); if(fp > 0.0) lo = p; else hi = p;
			double df = -a[0]*num[0]/((p+off[0])*(p+off[0])) - a[1]*num[1]/((p+off[1])*(p+off[1]));
			double pn = p - fp/df; if(!(pn > lo && pn < hi)) pn = 0.5*(lo+hi);
			if(std::fabs(fp) < 1e-13 || (hi-lo) <= 1e-15*(std::fabs(hi)+std::fabs(lo))){ break; } p = pn;   // 判据用残量: 根可能紧贴下界 (气相大幅膨胀), 步长判据会误停
		}
		double a1 = a[0]*num[0]/(p+off[0]), a2 = a[1]*num[1]/(p+off[1]);
		c.a1 = std::min(std::max(a1/(a1+a2), ALPHA_EPS), 1.0-ALPHA_EPS);
	}else g_relax_fail++;
	// 重初始化: 混合压力由守恒的总能给出, 再按各相 EOS 重置内能
	double rho = c.m[0]+c.m[1], rhoe = c.E - 0.5*c.mu*c.mu/rho - c.m[0]*QQ[0] - c.m[1]*QQ[1];     // 热内能 (总能里含生成焓 sum m_k q_k)
	double p = (rhoe - Bcoef(c.a1))/Acoef(c.a1);
	const double pfloor = 1.0;                                             // Pa; 仅在总能给出非正压力时触发 (计数)
	if(!(p > pfloor)){ p = pfloor; g_pfloor++; }
	for(int k=0;k<2;k++) c.ei[k] = alpha(c,k)*(p + GAM[k]*PIN[k])/(GAM[k]-1.0);
}

// ---- 热松弛 (pT) 与化学松弛 (pTg): Pelanti-Shyue 2014 的做法, 作为代数投影接在力学松弛 + 重初始化之后 ----
//  SG 带 q, q', cv:  v_k(p,T) = (g-1) cv T/(p+pi),  e_k(p,T) = cv T (p + g pi)/(p+pi) + q,  g_k(p,T) = (g cv - q') T - cv T ln(T^g/(p+pi)^(g-1)) + q
//  输入只用守恒量: 混合密度 rho, 相质量 m_k, 混合内能 rho e = E - 1/2 rho u^2 (含 q). 输出 alpha, (m_k), 再按 6+1 的规矩重置相内能.
static long g_pT_fail = 0, g_pTg_cells = 0, g_clfail = 0;
static double EPSI = 0.0; static int RAMP = 0;   // P-S (10): 热/化学松弛只在界面格 eps_I <= alpha_1 <= 1 - eps_I 打开 (0 = 处处打开)
static inline double sgv(int k,double p,double T){ return (GAM[k]-1.0)*CV[k]*T/(p+PIN[k]); }
static inline double sge(int k,double p,double T){ return CV[k]*T*(p+GAM[k]*PIN[k])/(p+PIN[k]) + QQ[k]; }
static inline double sgg(int k,double p,double T){ return (GAM[k]*CV[k]-QP[k])*T - CV[k]*T*std::log(std::pow(T,GAM[k])/std::pow(p+PIN[k],GAM[k]-1.0)) + QQ[k]; }
static inline double sgdg(int k,double p,double T){ return (GAM[k]*CV[k]-QP[k]) - CV[k]*std::log(std::pow(T,GAM[k])/std::pow(p+PIN[k],GAM[k]-1.0)) - CV[k]*GAM[k]; }   // dg/dT
static double Tsat_bis(double p){ double a=200.0, b=1500.0; auto f=[&](double T){ return sgg(0,p,T)-sgg(1,p,T); }; double fa=f(a), fb=f(b); if(!(fa*fb < 0.0)) return NAN;
	for(int it=0; it<80; it++){ double m=0.5*(a+b), fm=f(m); if(fa*fm <= 0.0){ b=m; } else { a=m; fa=fm; } } return 0.5*(a+b); }
static double PCRIT = -1.0, VCRIT = 0.0, VLC = 0.0, VVC = 0.0;   // 饱和线端点 (推导文档 11.3), 与二维 phys.h 相同
static double Tsat_raw(double p){ static double Tg = 370.0; double T = Tg;   // Newton (热启动), 不收敛退回二分
	for(int it=0; it<30; it++){ double f = sgg(0,p,T)-sgg(1,p,T), d = sgdg(0,p,T)-sgdg(1,p,T), Tn = T - f/d;
		if(!(Tn > 200.0 && Tn < 1500.0) || !std::isfinite(Tn)) break;
		if(std::fabs(Tn-T) < 1e-13*Tn){ Tg = Tn; return Tn; } T = Tn; }
	double Tb = Tsat_bis(p); if(std::isfinite(Tb)) Tg = Tb; return Tb; }
static void init_crit(){ auto f = [](double p){ const double T = Tsat_raw(p); return std::isfinite(T) ? sgv(1,p,T)-sgv(0,p,T) : (double)NAN; };
	PCRIT = 1e300; double a = 1e3, fa = f(a); for(double b = a*1.05; b < 1e12; b *= 1.05){ const double fb = f(b); if(!std::isfinite(fb)){ if(std::isfinite(fa)){ PCRIT = a; const double Tc = Tsat_raw(a); VLC = sgv(0,a,Tc); VVC = sgv(1,a,Tc); VCRIT = VLC; } break; }
		if(std::isfinite(fa) && fa > 0.0 && fb <= 0.0){ for(int it=0; it<100; it++){ const double m = std::sqrt(a*b); if(f(m) > 0.0) a = m; else b = m; } PCRIT = a; VCRIT = VLC = VVC = sgv(0,a,Tsat_raw(a)); break; }
		a = b; fa = fb; } }
static double Tsat(double p){ if(PCRIT < 0.0) init_crit(); return p < PCRIT ? Tsat_raw(p) : (double)NAN; }
template<class Fn> static bool solve_logp(Fn f, double pg, double& root){   // 在 ln p 上从猜值向外找变号再二分
	double a=std::log(std::max(pg,1e-3)), b=a, fa=f(std::exp(a)), fb=fa; if(!std::isfinite(fa)) return false; double st=0.05; bool found=false;
	for(int k=0;k<120 && !found;k++){ double an=a-st, bn=b+st, fan=f(std::exp(an)), fbn=f(std::exp(bn));
		if(std::isfinite(fan) && fan*fa <= 0.0){ b=a; fb=fa; a=an; fa=fan; found=true; break; }
		if(std::isfinite(fbn) && fbn*fb <= 0.0){ a=b; fa=fb; b=bn; fb=fbn; found=true; break; }
		if(std::isfinite(fan)){ a=an; fa=fan; } if(std::isfinite(fbn)){ b=bn; fb=fbn; } st*=1.5; }
	if(!found) return false;
	for(int it=0; it<80 && (b-a) > 1e-14; it++){ double m=0.5*(a+b), fm=f(std::exp(m)); if(!std::isfinite(fm)) return false; if(fa*fm <= 0.0){ b=m; } else { a=m; fa=fm; } }
	root=std::exp(0.5*(a+b)); return true; }
static void thermo_relax(Cell& c){
	if(RELAX < 1) return;
	if(!(c.a1 >= EPSI && c.a1 <= 1.0-EPSI)){ c.eq = 0; return; }
	const double rho=c.m[0]+c.m[1], v=1.0/rho, e=(c.E - 0.5*c.mu*c.mu/rho)/rho, p0=pmix(c);
	double Y[2]={c.m[0]/rho, c.m[1]/rho};
	// pT: T(p) 由混合能量给出, 对 p 解 sum Y_k v_k(p,T) = v
	auto Tof=[&](double p){ return (e - Y[0]*QQ[0] - Y[1]*QQ[1])/(Y[0]*CV[0]*(p+GAM[0]*PIN[0])/(p+PIN[0]) + Y[1]*CV[1]*(p+GAM[1]*PIN[1])/(p+PIN[1])); };
	double p; if(!solve_logp([&](double q){ double T=Tof(q); return Y[0]*sgv(0,q,T)+Y[1]*sgv(1,q,T)-v; }, p0, p)){ g_pT_fail++; return; }
	double T=Tof(p); if(!(T > 0.0)){ g_pT_fail++; return; }
	if(RELAX == 3){   // p-pTg (P-S 的第四档): 跳过 pT; 判据 = 力学平衡态下任一相温度 T_k > T_sat(p_eq). 不满足则保持力学松弛态不动
		const double Tk[2] = { (p0+PIN[0])*(c.a1/c.m[0])/((GAM[0]-1.0)*CV[0]), (p0+PIN[1])*((1.0-c.a1)/c.m[1])/((GAM[1]-1.0)*CV[1]) };
		p = p0; T = std::max(Tk[0], Tk[1]); }
	if(RELAX >= 2){   // pTg: 只在液体亚稳 (T > T_sat(p)) 时做 (P-S 的判据); T = T_sat(p), Y_v 由体积定, 对 p 解混合能量
		double Ts=Tsat(p);
		const bool was = c.eq; c.eq = 0;
		if(std::isfinite(Ts) && (T > Ts || (was && T > Ts*(1.0-1e-9)))){
			auto Yv=[&](double q,double Tq){ double vl=sgv(0,q,Tq), vv=sgv(1,q,Tq); return (v-vl)/(vv-vl); };
			double pe; bool ok = solve_logp([&](double q){ double Tq=Tsat(q); if(!std::isfinite(Tq)) return (double)NAN; double y=Yv(q,Tq); return y*sge(1,q,Tq)+(1.0-y)*sge(0,q,Tq)-e; }, p, pe);
			if(ok){ double Te=Tsat(pe), y=Yv(pe,Te); if(y > 0.0 && y < 1.0){ p=pe; T=Te; Y[1]=y; Y[0]=1.0-y; c.m[0]=Y[0]*rho; c.m[1]=Y[1]*rho; g_pTg_cells++; c.eq = 1; } }
		}
	}
	if(RELAX == 3 && !c.eq) return;
	double a1 = c.m[0]*sgv(0,p,T); c.a1 = std::min(std::max(a1, ALPHA_EPS), 1.0-ALPHA_EPS);
	for(int k=0;k<2;k++) c.ei[k] = alpha(c,k)*(p + GAM[k]*PIN[k])/(GAM[k]-1.0);
}

// 饱和线上的平衡混合内能 (含 q): rho e_eq(rho, p),  T = T_sat(p),  Y_v = (v - v_l)/(v_v - v_l). 无两相根 -> false
static bool Eeq(double rho, double p, double& E, double* Yv = nullptr, double* Tout = nullptr){
	//  两相根不存在时 (y 落在 (0,1) 之外) 不换闭合, 而是沿同一个平衡 EOS 连续延拓: Y_v 钳在界上 (全凝结 / 全蒸发), 温度由 v = sum Y_k v_k(p,T) 定 (pT 平衡).
	//  在 y = 界处 T = T_sat, 两支的 (v, e) 相同, 所以 rho e_eq(rho, p) 对 (rho, p) 连续.
	if(!(p > 0.0)) return false; if(PCRIT < 0.0) init_crit(); const double YMIN = 1e-9; double T, y;
	if(p < PCRIT){ T = Tsat(p); if(!std::isfinite(T)) return false; const double vl = sgv(0,p,T), vv = sgv(1,p,T); y = (1.0/rho - vl)/(vv - vl); if(!std::isfinite(y)) return false; }
	else y = (VVC > VLC*(1.0+1e-12)) ? (1.0/rho - VLC)/(VVC - VLC) : ((1.0/rho < VCRIT) ? -1.0 : 2.0);
	if(p >= PCRIT || !(y > YMIN && y < 1.0-YMIN)){ y = std::min(std::max(y, YMIN), 1.0-YMIN);
		T = (1.0/rho)/((1.0-y)*(GAM[0]-1.0)*CV[0]/(p+PIN[0]) + y*(GAM[1]-1.0)*CV[1]/(p+PIN[1])); if(!(T > 0.0)) return false; }
	E = rho*(y*sge(1,p,T) + (1.0-y)*sge(0,p,T)); if(Yv) *Yv = y; if(Tout) *Tout = T; return std::isfinite(E); }

// pT 平衡混合内能 (含 q), 质量分数 y 给定: T 由 v = sum Y_k v_k(p,T) 定
static bool EpT(double rho, double p, double y, double& E, double* Tout = nullptr){
	if(!(p > 0.0) || !(y > 0.0 && y < 1.0)) return false;
	const double T = (1.0/rho)/((1.0-y)*(GAM[0]-1.0)*CV[0]/(p+PIN[0]) + y*(GAM[1]-1.0)*CV[1]/(p+PIN[1])); if(!(T > 0.0)) return false;
	E = rho*(y*sge(1,p,T) + (1.0-y)*sge(0,p,T)); if(Tout) *Tout = T; return std::isfinite(E); }

// ---- 二阶: 分片线性重构 (minmod), 变量 = 状态向量各分量; 面态不可容许时退回格均值 ----
static int ORDER = 1; static long g_rk1 = 0, g_qelim = 0;
static void reconstruct(const std::vector<Cell>& U, int N, std::vector<Cell>& UL, std::vector<Cell>& UR){
	// 原始变量 W = (alpha1, rho1, rho2, u, p1, p2) 上做 minmod 重构 (论文附录 C 也用原始变量), 再由面上的原始量组装状态.
	// 不能对守恒量逐分量重构: 空穴里 p = (rho e - B(alpha))/A(alpha) 是两个大数之差, alpha 与能量的斜率各自限幅后不再相容, 压力会掉到零以下.
	UL = U; UR = U; if(ORDER < 2) return;
	const int M = N+2*NG; std::vector<double> W((size_t)M*6);
	for(int i=0;i<M;i++){ const Cell& c = U[i]; double* w = &W[(size_t)i*6]; double r = c.m[0]+c.m[1];
		w[0]=c.a1; w[1]=c.m[0]/c.a1; w[2]=c.m[1]/(1.0-c.a1); w[3]=c.mu/r; w[4]=pk(c,0); w[5]=pk(c,1); }
	auto build = [](const double* w)->Cell{ Cell c; c.a1=w[0]; c.m[0]=w[0]*w[1]; c.m[1]=(1.0-w[0])*w[2]; double r=c.m[0]+c.m[1]; c.mu=r*w[3];
		c.ei[0]=w[0]*(w[4]+GAM[0]*PIN[0])/(GAM[0]-1.0); c.ei[1]=(1.0-w[0])*(w[5]+GAM[1]*PIN[1])/(GAM[1]-1.0); c.E=c.ei[0]+c.ei[1]+0.5*r*w[3]*w[3]+c.m[0]*QQ[0]+c.m[1]*QQ[1]; return c; };
	for(int i=1;i<M-1;i++){ double l[6], r[6];
		for(int q=0;q<6;q++){ double c0=W[(size_t)i*6+q], a=c0-W[(size_t)(i-1)*6+q], b=W[(size_t)(i+1)*6+q]-c0;
			double sl = (a*b > 0.0) ? (std::fabs(a) < std::fabs(b) ? a : b) : 0.0; l[q]=c0-0.5*sl; r[q]=c0+0.5*sl; }
		UL[i] = build(l); UR[i] = build(r); }
}

// ---- MODE 0: 论文的 HLLC + Godunov ----
struct Face { double F[4], ua, u, aei_u[2]; };   // 守恒通量 (m1,m2,rho u,rho E), (u alpha1)*, u*, (alpha rho e u)*_k
static Face hllc(const Cell& L, const Cell& R){
	auto prim = [&](const Cell& c, double& rho, double& u, double& p, double& cs){ rho=c.m[0]+c.m[1]; u=c.mu/rho; p=pmix(c); cs=cfrozen(c); };
	double rL,uL,pL,cL,rR,uR,pR,cR; prim(L,rL,uL,pL,cL); prim(R,rR,uR,pR,cR);
	double SL = std::min(uL-cL, uR-cR), SR = std::max(uL+cL, uR+cR);
	double SM = ((rR*uR*uR+pR) - (rL*uL*uL+pL) - SR*rR*uR + SL*rL*uL)/(rR*uR - rL*uL - SR*rR + SL*rL);   // 论文式 (符号整理为标准 HLL 接触速度)
	auto flux = [&](const Cell& c, double u, double p, double* F){ F[0]=c.m[0]*u; F[1]=c.m[1]*u; F[2]=c.mu*u+p; F[3]=(c.E+p)*u; };
	Face f; const Cell& K = (SM >= 0.0) ? L : R; double rK = (SM>=0)?rL:rR, uK = (SM>=0)?uL:uR, pKm = (SM>=0)?pL:pR, SK = (SM>=0)?SL:SR;
	double FK[4]; flux(K, uK, pKm, FK);
	if((SM >= 0.0 && SL >= 0.0) || (SM < 0.0 && SR <= 0.0)){   // 超声速: 迎风态
		for(int q=0;q<4;q++) f.F[q]=FK[q]; f.u=uK; f.ua=uK*K.a1; for(int k=0;k<2;k++) f.aei_u[k]=K.ei[k]*uK; return f; }
	double fac = (SK-uK)/(SK-SM), rs = rK*fac;
	double ps = pKm + rK*(uK-SK)*(uK-SM);
	double Es = (K.E*(uK-SK) + pKm*uK - ps*SM)/(rs*(SM-SK))*rs;           // (rho E)*
	double Us[4] = {K.m[0]*fac, K.m[1]*fac, rs*SM, Es}, UK[4] = {K.m[0], K.m[1], K.mu, K.E};
	for(int q=0;q<4;q++) f.F[q] = FK[q] + SK*(Us[q]-UK[q]);
	f.u = SM; f.ua = SM*K.a1;
	for(int k=0;k<2;k++){   // 相内能沿各相 Hugoniot (II.8): p_k*(rho_k*) 再由 EOS 得 e_k*
		double a = alpha(K,k), r0 = K.m[k]/a, r1 = r0*fac, p0 = pk(K,k), g = GAM[k];
		double den = (g-1.0)*r1 - (g+1.0)*r0, p1 = (std::fabs(den) > 1e-300) ? (p0+PIN[k])*((g-1.0)*r0 - (g+1.0)*r1)/den - PIN[k] : p0;
		if(!(p1 + PIN[k] > 0.0) || !std::isfinite(p1)) p1 = p0;
		f.aei_u[k] = a*(p1 + g*PIN[k])/(g-1.0)*SM;
	}
	return f;
}
static void step_explicit(std::vector<Cell>& U, int N, double dt, double dx){
	std::vector<Face> F(N+2*NG);
	std::vector<Cell> UL, UR; reconstruct(U, N, UL, UR);
	for(int i=NG-1;i<N+NG;i++) F[i] = hllc(UR[i], UL[i+1]);               // 面 i = (i, i+1)
	std::vector<Cell> V = U;
	for(int i=NG;i<N+NG;i++){
		const Face &fp = F[i], &fm = F[i-1]; const Cell& c = U[i]; double r = dt/dx;
		V[i].m[0] = c.m[0] - r*(fp.F[0]-fm.F[0]); V[i].m[1] = c.m[1] - r*(fp.F[1]-fm.F[1]);
		V[i].mu = c.mu - r*(fp.F[2]-fm.F[2]); V[i].E = c.E - r*(fp.F[3]-fm.F[3]);
		V[i].a1 = c.a1 - r*((fp.ua-fm.ua) - c.a1*(fp.u-fm.u));
		for(int k=0;k<2;k++) V[i].ei[k] = c.ei[k] - r*((fp.aei_u[k]-fm.aei_u[k]) + alpha(c,k)*pk(c,k)*(fp.u-fm.u));
		V[i].a1 = std::min(std::max(V[i].a1, ALPHA_EPS), 1.0-ALPHA_EPS);
	}
	U = V;
}

// ---- MODE 1: IMEX ----
//  分裂 (Toro-Vazquez; 内能方程按同一原则分):
//    对流 (显式):  alpha_t + u alpha_x = 0;  (m_k)_t + (m_k u)_x = 0;  (rho u)_t + (rho u^2)_x = 0;  (rho E)_t + (1/2 rho u^2 u)_x = 0;  相内能无对流通量
//    压力 (隐式):  (rho u)_t + p_x = 0;  (alpha rho e)_k,t + ((alpha rho e)_k u)_x + (alpha p)_k u_x = 0;  (rho E)_t + ((rho e + p) u)_x = 0
//  显式: path-conservative Rusanov, D^{-/+} = 1/2 (F_R - F_L + ubar (alpha_R - alpha_L) e_alpha) -/+ 1/2 s (Q_R - Q_L), s = max|u| (全域)
//  隐式: 由混合 EOS (II.4) rho E = A(alpha) p + B(alpha) + 1/2 rho u^2 与 u^{n+1} = u* - dt p_x / rho 得
//        A_i p_i - dt^2/dx^2 d_i[ H_f (1/rho)_f d_f p ] = (rho E)*_i - B_i - K_i - dt/dx d_i[ H_f ubar*_f ],   H = rho e + p (迎风),  K 用 Picard 迭代
//  级算子 (半隐式 IMEX-RK 的一级): 显式通量在 Qe 上算, 更新的基底是 Qs, 步长 dt = a_ii * Dt. 结果写回 U (已松弛).
static void step_imex(std::vector<Cell>& U, const std::vector<Cell>& Qs, int N, double dt, double dx, int picard){   // 进入时 U = Qe
	const int M = N+2*NG; double sglob = 0.0;
	for(int i=NG;i<N+NG;i++) sglob = std::max(sglob, std::fabs(U[i].mu/(U[i].m[0]+U[i].m[1])));
	// ---------- 显式对流 ----------
	std::vector<Cell> V = Qs; std::vector<double> axr[2], af[2]; for(int k=0;k<2;k++){ axr[k].assign(M,0.0); af[k].assign(M,0.0); }
	auto toQ = [&](const Cell& c, double* Q){ Q[0]=c.a1; Q[1]=c.m[0]; Q[2]=c.m[1]; Q[3]=c.mu; Q[4]=c.ei[0]; Q[5]=c.ei[1]; Q[6]=c.E; };
	auto conv = [&](const Cell& c, double* F){ double r=c.m[0]+c.m[1], u=c.mu/r; F[0]=AVOL*c.a1*u; F[1]=c.m[0]*u; F[2]=c.m[1]*u; F[3]=c.mu*u; F[4]=EXPLE*c.ei[0]*u; F[5]=EXPLE*c.ei[1]*u; F[6]=0.5*c.mu*u*u + (c.m[0]*QQ[0]+c.m[1]*QQ[1])*u + EXPLE*(c.ei[0]+c.ei[1])*u; };   // 生成焓随显式质量通量走
	std::vector<Cell> UL_, UR_; reconstruct(U, N, UL_, UR_);
	std::vector<double> ue(M, 0.0);
	std::vector<double> Dm((size_t)M*NV, 0.0), Dp((size_t)M*NV, 0.0);       // 面 i 的 D^- (给格 i) 与 D^+ (给格 i+1)
	for(int i=NG-1;i<N+NG;i++){
		const Cell &L = UR_[i], &R = UL_[i+1]; double QL[NV], QR[NV], FL[NV], FR[NV]; toQ(L,QL); toQ(R,QR); conv(L,FL); conv(R,FR);
		double ub = 0.5*(L.mu/(L.m[0]+L.m[1]) + R.mu/(R.m[0]+R.m[1])); ue[i] = ub;
		const double s = SLOC ? std::max(std::max(std::fabs(L.mu/(L.m[0]+L.m[1])), std::fabs(R.mu/(R.m[0]+R.m[1]))), std::max(std::fabs(U[i].mu/(U[i].m[0]+U[i].m[1])), std::fabs(U[i+1].mu/(U[i+1].m[0]+U[i+1].m[1])))) : sglob;   // 耗散速度: 面局部 (SLOC=1) 或全场最大
		for(int q=0;q<NV;q++){ double jump = FR[q]-FL[q] + ((q==0 && !AVOL) ? ub*(QR[0]-QL[0]) : 0.0);
			Dm[(size_t)i*NV+q] = 0.5*jump - 0.5*s*(QR[q]-QL[q]); Dp[(size_t)i*NV+q] = 0.5*jump + 0.5*s*(QR[q]-QL[q]); }
	}
	for(int i=NG;i<N+NG;i++){ double Q[NV]; toQ(Qs[i],Q);
		for(int q=0;q<NV;q++) Q[q] -= dt/dx*(Dm[(size_t)i*NV+q] + Dp[(size_t)(i-1)*NV+q]);
		{ // 格内光滑部分: 守恒通量差 F(Q_iR) - F(Q_iL), 加上 alpha 的非守恒项 u_i (alpha_iR - alpha_iL)
		  double FLc[NV], FRc[NV]; conv(UL_[i],FLc); conv(UR_[i],FRc); double ui = U[i].mu/(U[i].m[0]+U[i].m[1]);
		  for(int q=0;q<NV;q++) Q[q] -= dt/dx*(FRc[q]-FLc[q]); if(!AVOL) Q[0] -= dt/dx*ui*(UR_[i].a1-UL_[i].a1); }
		if(AVOL){   // 相体积按与质量、内能同一个守恒通量输运 (不钳位): a1x + a2x = 1 - r d_i[ubar^e_f], 差额就是本格的体积应变
			auto uo = [&](const Cell& z){ return z.mu/(z.m[0]+z.m[1]); };
			const double divu = 0.5*(uo(UR_[i])+uo(UL_[i+1])) - 0.5*(uo(UR_[i-1])+uo(UL_[i]));
			axr[0][i] = Q[0]; axr[1][i] = 1.0 - dt/dx*divu - Q[0];
			if(!(axr[0][i] > 0.0 && axr[1][i] > 0.0 && Q[1] > 0.0 && Q[2] > 0.0)) g_trbad++;   /* 输运出了负的相体积/质量: 本级不保界 */
			 Q[0] = axr[0][i]/(axr[0][i]+axr[1][i]); }
		Cell& c = V[i]; c.a1=std::min(std::max(Q[0],ALPHA_EPS),1.0-ALPHA_EPS); c.m[0]=Q[1]; c.m[1]=Q[2]; c.mu=Q[3]; c.ei[0]=Q[4]; c.ei[1]=Q[5]; c.E=Q[6]; }
	for(int g=0;g<NG;g++){ V[g] = V[NG]; V[N+NG+g] = V[N+NG-1]; }
	// ---------- 隐式压力 + 松弛 (耦合) ----------
	//  物质 CFL 下一步之内气相体积分数可以变 O(1) (空化), 所以松弛不能放在一个 "alpha 冻结" 的声学隐式步之后:
	//  冻结 alpha 的混合 EOS 给的是液体的刚度, 而真实响应是 Wood (平衡) 刚度. 因此把论文 3.3 的松弛关系放进隐式方程:
	//    面速度      u_f = ubar*_f - dt (1/rho)_f (p_{i+1} - p_i)/dx
	//    相内能      (a rho e)^0_k = (a rho e)^x_k - dt/dx d_i[ (a rho e)_k,up u_f ] - dt/dx alpha*_k p_i d_i[u_f]        (论文 3.2.1, 压力功取 p^{n+1})
	//    松弛 (3.3)  e_k - e_k^0 + \hat p_I (v_k - v_k^0) = 0, \hat p_I = p (论文允许的两种估计之一)
	//                => alpha_k [ (p + g_k pi_k)/(g_k - 1) + p ] = (a rho e)^0_k + p alpha*_k  ==> alpha_k(p)
	//    约束        G_i(p) = alpha_1 + alpha_2 - 1 = 0        (每格一个未知量 p_i; Jacobian 三对角; Newton)
	std::vector<double> rho(M), us(M), p(M), irf(M), ubf(M), uf(M), ef[2];
	ef[0].assign(M,0.0); ef[1].assign(M,0.0);
	for(int i=0;i<M;i++){ const Cell& c = V[i]; rho[i]=c.m[0]+c.m[1]; us[i]=c.mu/rho[i]; const Cell& b = EXPLE ? Qs[i] : c; p[i] = (b.ei[0]+b.ei[1] - Bcoef(b.a1))/Acoef(b.a1); }   // 初猜: 显式半步输运了内能但没有压力功, V 的压力没有意义, 取级基底的压力
	std::vector<Cell> VLr, VRr; reconstruct(V, N, VLr, VRr);
	for(int i=NG-1;i<N+NG;i++){ ubf[i] = 0.5*(us[i]+us[i+1]); irf[i] = IRF ? 2.0/(rho[i]+rho[i+1]) : 0.5*(1.0/rho[i]+1.0/rho[i+1]);
		// 被输运的相内能取面两侧 (重构值) 的算术平均: 显式半步已经对它加了 Rusanov 耗散 -1/2 s d(a rho e), 两者合起来才是一个完整的 Rusanov 通量,
		// 且与显式质量通量的中心部分 1/2 (m_L u_L + m_R u_R) 同一种面平均. 这里若取迎风, 质量按平均走而内能按迎风走, 比内能 (温度) 在 alpha 梯度大的格里系统性偏高.
		for(int k=0;k<2;k++) ef[k][i] = CENTRAL ? 0.5*(VRr[i].ei[k]+VLr[i+1].ei[k]) : (ubf[i] > 0.0 ? VRr[i].ei[k] : (ubf[i] < 0.0 ? VLr[i+1].ei[k] : 0.5*(VRr[i].ei[k]+VLr[i+1].ei[k]))); }
	// 质量跟随隐式速度增量 du_f = u_f^{n+1} - ubar*_f (与相内能同一迎风侧): 否则热内能按 u^{n+1} 走而质量按 u^n 走,
	//  比内能 (温度) 每步错 O(du dt/dx); 力学松弛下 SG 的压力与密度无关所以看不出来, 热/化学松弛下温度直接定饱和压力.
	std::vector<double> mf[2], uup(M,0.0); mf[0].assign(M,0.0); mf[1].assign(M,0.0);
	for(int i=NG-1;i<N+NG;i++){ for(int k=0;k<2;k++) mf[k][i] = CENTRAL ? 0.5*(VRr[i].m[k]+VLr[i+1].m[k]) : (ubf[i] > 0.0 ? VRr[i].m[k] : (ubf[i] < 0.0 ? VLr[i+1].m[k] : 0.5*(VRr[i].m[k]+VLr[i+1].m[k])));
		for(int k=0;k<2;k++){ const double aL = alpha(VRr[i],k), aR = alpha(VLr[i+1],k); af[k][i] = CENTRAL ? 0.5*(aL+aR) : (ubf[i] > 0.0 ? aL : (ubf[i] < 0.0 ? aR : 0.5*(aL+aR))); }
		uup[i] = CENTRAL ? ubf[i] : (ubf[i] > 0.0 ? us[i] : (ubf[i] < 0.0 ? us[i+1] : 0.5*(us[i]+us[i+1]))); }
	// 增量通量的参考速度: 必须就是显式中心通量里用的那个面速度 ubar^e_f = 1/2 (u_L + u_R)|_{Qe}, 这样总体积通量恰为 u_f^{n+1}:
	//   sum_k alpha^tr_k = 1 - r d_i[u_f^{n+1}]  (散度完全隐式). 若用 ubar*_f (来自 Qs), 二阶的第二级里 ubar^e - ubar* 含有上一级的压力梯度加速度,
	//   等于把刚性的声学压缩项显式外推, 近纯液体里声学 CFL 大时不稳定.
	//   UREF = 2 (默认): 只把那一块刚性的差去掉, uR = ubar* - [ubar(Qs) - ubar(Qe)]; 显式动量对流引起的 ubar* 变化 (非刚性, 激波处可达 O(u)) 仍留在 ubar* 里,
	//   不进增量通量 (否则增量通量可以搬空一格, 破坏正性). 一阶 (Qs = Qe) 时 uR = ubar*.   UREF = 1: uR = ubar^e (散度严格隐式, 但激波处 Newton 不收敛).
	auto ucell = [&](const Cell& z){ return z.mu/(z.m[0]+z.m[1]); };
	std::vector<double> uR(M, 0.0); for(int i=NG-1;i<N+NG;i++) uR[i] = !AVOL || UREF == 0 ? ubf[i] : (UREF == 1 ? ue[i] : ubf[i] - 0.5*(ucell(Qs[i])+ucell(Qs[i+1])) + 0.5*(ucell(U[i])+ucell(U[i+1])));
	// DONOR = 1: 增量通量 Xhat_f du_f 的面值按增量自己的方向取供体格 (du_f > 0 取左, 否则取右), 而不是按 ubar*_f 的符号
	auto setface = [&](int f, double d){ if(!DONOR || CENTRAL) return; const Cell& z = d > 0.0 ? VRr[f] : VLr[f+1];
		for(int k=0;k<2;k++){ ef[k][f] = z.ei[k]; mf[k][f] = z.m[k]; af[k][f] = alpha(z,k); } uup[f] = d > 0.0 ? us[f] : us[f+1]; };
	const double r = dt/dx;
	// ---- 激波处的隐式耗散 (FLAT=1): 指示器 = Balsara 的 Euler flattener (声速), 耗散 = 压力方程里的隐式 Laplace 项 (eps-MHD 文稿 (43)-(45)) ----
	//   chi~_i = min[1, max(0, -(dx div u + K1 c_min)/(K1 c_min))],  dx div u = u_{i+1} - u_{i-1},  c_min = 三格模板内的最小 (冻结) 声速;  chi_i = 模板内 chi~ 的最大值
	//   各相热内能的耗散通量 F^d_k,f = -kappa_f V'_k,f (p_R - p_L)/dx,  kappa_f = chi_f |u|_f dx,  V'_k = alpha_k/(gamma_k - 1)  (sum_k V'_k = d(rho e)/dp|_alpha)
	//   各相热内能的耗散通量 F^d_k,f = -kappa_f V'_k,f (p_R - p_L)/dx, V'_k,f = alpha_k,f/(gamma_k-1), 各相各自守恒; 总能的通量是两相之和. 质量、动量不动; p 均匀时为零.
	//   kapR[k][i] = 格 i 右面上相 k 的系数 (dt/dx^2) kappa V'_k, kapL[k][i] = 左面上的.
	std::vector<double> KA(M,0.0), kapR[2], kapL[2]; for(int k=0;k<2;k++){ kapR[k].assign(M,0.0); kapL[k].assign(M,0.0); }   // KA[f] = (dt/dx^2) kappa_f V'_f;  kapR/kapL[k][i] = 本格相 k 的份额 * KA[i] / KA[i-1]
	if(FLAT){ std::vector<double> ch(M,0.0), chi(M,0.0), cs(M), uq(M);
		for(int i=0;i<M;i++){ cs[i] = cfrozen(U[i]); uq[i] = U[i].mu/(U[i].m[0]+U[i].m[1]); }
		for(int i=1;i<M-1;i++){ const double cm = std::min(cs[i], std::min(cs[i-1], cs[i+1])); ch[i] = std::min(1.0, std::max(0.0, -((uq[i+1]-uq[i-1]) + FK1*cm)/(FK1*cm))); }
		for(int i=1;i<M-1;i++) chi[i] = std::max(ch[i], std::max(ch[i-1], ch[i+1]));
		for(int f=NG;f<N+NG-1;f++){ const double kf = r*std::max(chi[f],chi[f+1])*(std::max(std::fabs(uq[f]),std::fabs(uq[f+1])) + FSC*std::max(cs[f],cs[f+1]));   // 耗散速度: |u| (FSC=0, 对流量级) 或 |u| + c (FSC=1, 声学量级)
			double ka = 0.0; for(int k=0;k<2;k++){ const double aL = alpha(V[f],k), aR = alpha(V[f+1],k), af_ = FLATA == 0 ? std::min(aL,aR) : (FLATA == 1 ? 0.5*(aL+aR) : 2.0*aL*aR/(aL+aR));
				kapR[k][f] = kf*af_/(GAM[k]-1.0); kapL[k][f+1] = kapR[k][f]; ka += kapR[k][f]; } KA[f] = ka; }   // 各相各自守恒: 相 k 的内能只在相 k 内扩散, 面系数取两侧 alpha_k 的较小者 (一侧没有这一相就不通)
		for(int i=NG;i<N+NG;i++) if(chi[i] > 0.0) g_flatcells++; }
	std::vector<double> G(M), lo(M), di(M), up(M), dp(M), al[2]; al[0].assign(M,0.0); al[1].assign(M,0.0); std::vector<int> useq(M,0);
	// 热力学闭合的单格残量 (只依赖本格与左右邻格的压力):  G = [ rho e^0 - rho e_cl(rho, p, Y) ] / S
	//   rho e^0 = sum_k (a rho e)^0_k + sum_k q_k m_k,  m_k = m*_k - dt/dx d_i[ m_k,f du_f ]
	//   cl = 1: pT 平衡 (Y 冻结 = m_2/rho);  cl = 2: pTg 平衡 (饱和线, Y_v 由体积定, 出界则钳住退化为 pT)
	static const bool noeq = getenv("NOEQ") != nullptr, nopt = getenv("NOPT") != nullptr;
	auto clos = [&](int i)->int{ if(noeq) return 0; if(!(V[i].a1 >= EPSI && V[i].a1 <= 1.0-EPSI)) return 0; if(RELAX >= 2 && V[i].eq) return 2; if((RELAX == 1 || RELAX == 2) && !nopt) return 1; return 0; };
	auto Gth = [&](int i, int cl, double pl, double pc, double pr, double& g)->bool{
		const double ufr = (i < N+NG-1) ? ubf[i] - r*irf[i]*(pr - pc) : ubf[i], ufl = (i > NG) ? ubf[i-1] - r*irf[i-1]*(pc - pl) : ubf[i-1];
		const double du = ufr-ufl, duR = ufr-uR[i], duL = ufl-uR[i-1]; setface(i, duR); setface(i-1, duL);
		const double m0 = V[i].m[0] - r*(mf[0][i]*duR - mf[0][i-1]*duL), m1 = V[i].m[1] - r*(mf[1][i]*duR - mf[1][i-1]*duL), rn = m0+m1;
		if(!(m0 > 0.0 && m1 > 0.0)) return false;
		const double E0 = V[i].ei[0]+V[i].ei[1] - r*((ef[0][i]+ef[1][i])*(ufr-EXPLE*uR[i]) - (ef[0][i-1]+ef[1][i-1])*(ufl-EXPLE*uR[i-1])) + (AVOL ? pc*(axr[0][i]+axr[1][i]-1.0 - r*((af[0][i]+af[1][i])*duR - (af[0][i-1]+af[1][i-1])*duL)) : -r*pc*du) + m0*QQ[0] + m1*QQ[1]
		                  + KA[i]*(pr-pc) - KA[i-1]*(pc-pl);
		double Ee; if(!(cl == 2 ? Eeq(rn, pc, Ee) : EpT(rn, pc, m1/rn, Ee))) return false;
		g = (E0 - Ee)/(std::fabs(V[i].ei[0]+V[i].ei[1]) + 1.0); return std::isfinite(g); };
	auto eval = [&](const std::vector<double>& q, bool jac)->bool{   // 残量 (和 Jacobian); 返回可容许性
		for(int i=NG-1;i<N+NG;i++){ int a = std::max(i,NG), b = std::min(i+1,N+NG-1); uf[i] = ubf[i] - r*irf[i]*(q[b]-q[a]); setface(i, uf[i]-uR[i]); }
		for(int i=NG;i<N+NG;i++){
			double g = -1.0, dc = 0.0, dl = 0.0, dr = 0.0;
			useq[i] = 0;
			if(const int cl = clos(i)){   // 热力学闭合; Jacobian 用差分 (残量只依赖三个压力)
				const double pl = q[std::max(i-1,NG)], pc = q[i], pr = q[std::min(i+1,N+NG-1)]; double g0, g1;
				if(Gth(i, cl, pl, pc, pr, g0)){ G[i] = g0; useq[i] = cl; bool okj = true;
					if(jac){ const double h = 1e-6*pc + 1e-9;
						okj = Gth(i,cl,pl,pc+h,pr,g1); if(okj) di[i] = (g1-g0)/h;
						if(okj && i < N+NG-1){ const double hr = 1e-6*std::fabs(pr)+1e-9; okj = Gth(i,cl,pl,pc,pr+hr,g1); if(okj) up[i] = (g1-g0)/hr; } else up[i] = 0.0;
						if(okj && i > NG){ const double hl = 1e-6*std::fabs(pl)+1e-9; okj = Gth(i,cl,pl+hl,pc,pr,g1); if(okj) lo[i] = (g1-g0)/hl; } else if(i == NG) lo[i] = 0.0; }
					setface(i, uf[i]-uR[i]); setface(i-1, uf[i-1]-uR[i-1]); if(okj) continue; useq[i] = 0; }
				setface(i, uf[i]-uR[i]); setface(i-1, uf[i-1]-uR[i-1]); g_clfail++;
			}
			for(int k=0;k<2;k++){
				const double as = alpha(V[i],k), du = uf[i]-uf[i-1];
				const double atr = AVOL ? axr[k][i] - r*(af[k][i]*(uf[i]-uR[i]) - af[k][i-1]*(uf[i-1]-uR[i-1])) : 0.0;   // 输运后的相体积 (未压回 1)
				const double Nk = kapR[k][i]*(q[std::min(i+1,N+NG-1)]-q[i]) - kapL[k][i]*(q[i]-q[std::max(i-1,NG)]) + (AVOL ? V[i].ei[k] - r*(ef[k][i]*(uf[i]-uR[i]) - ef[k][i-1]*(uf[i-1]-uR[i-1])) + q[i]*atr
				                       : V[i].ei[k] - r*(ef[k][i]*(uf[i]-EXPLE*uR[i]) - ef[k][i-1]*(uf[i-1]-EXPLE*uR[i-1])) - r*as*q[i]*du + q[i]*as);
				if(!(Nk > 0.0) || !(q[i] + PIN[k] > 0.0)) return false;
				const double cf = (GAM[k]-1.0)/(GAM[k]*(q[i]+PIN[k])), ak = cf*Nk; al[k][i] = ak; g += ak;
				if(jac){
					const double wR = (i < N+NG-1) ? r*irf[i] : 0.0, wL = (i > NG) ? r*irf[i-1] : 0.0;      // 透射边界: 边界面压差恒为 0
					dc -= cf*(kapR[k][i]+kapL[k][i]); dr += cf*kapR[k][i]; dl += cf*kapL[k][i];
					if(AVOL){ dc += cf*(-r*((ef[k][i]+q[i]*af[k][i])*wR + (ef[k][i-1]+q[i]*af[k][i-1])*wL) + atr) - ak/(q[i]+PIN[k]);
						dr += cf*r*wR*(ef[k][i] + q[i]*af[k][i]); dl += cf*r*wL*(ef[k][i-1] + q[i]*af[k][i-1]); }
					else{
					const double dNc = -r*(ef[k][i]*wR + ef[k][i-1]*wL) - r*as*du - r*as*q[i]*(wR+wL) + as;
					dc += cf*dNc - ak/(q[i]+PIN[k]); dr += cf*r*wR*(ef[k][i] + as*q[i]); dl += cf*r*wL*(ef[k][i-1] + as*q[i]); }
				}
			}
			G[i] = g; if(jac){ di[i] = dc; up[i] = dr; lo[i] = dl; }
		}
		return true; };
	bool cvg = false; std::vector<double> q(M);
	if(getenv("JAC") && g_newton_steps==0){ eval(p,true); std::vector<double> G0=G, d0=di, u0=up, l0=lo; int ic=NG+N/2;
		for(int j=ic-1;j<=ic+1;j++){ std::vector<double> pp=p; double h_=1e-3*pp[j]; pp[j]+=h_; eval(pp,false);
			fprintf(stderr,"dG[%d]/dp[%d]: fd %.6e  an %.6e\n",ic,j,(G[ic]-G0[ic])/h_, j==ic?d0[ic]:(j>ic?u0[ic]:l0[ic])); } eval(p,true); }
	// 预扫 (非线性 Gauss-Seidel): 邻格压力取当前值, 每格解自己的 G_i(q) = 0. G_i 对 q 单调减 (q -> 0+: 气相体积 -> inf; q 大: 体积 -> 0),
	//  带括号的割线/二分一定有解. 作用: 空化格里压力要一步降两个量级, 1/p 型的 G 让 Newton 从高压侧起步时发散; 预扫把每格先放到正确量级.
	auto Gloc = [&](int i, double val)->double{
		const double pl = (i > NG) ? p[i-1] : val, pr = (i < N+NG-1) ? p[i+1] : val;
		const double ufr = ubf[i] - r*irf[i]*(pr - val), ufl = ubf[i-1] - r*irf[i-1]*(val - pl), du = ufr - ufl; double g = -1.0;
		if(const int cl = clos(i)){ double gg; if(Gth(i, cl, pl, val, pr, gg)) return gg; }
		setface(i, ufr-uR[i]); setface(i-1, ufl-uR[i-1]);
		for(int k=0;k<2;k++){ const double as = alpha(V[i],k);
			const double atr = AVOL ? axr[k][i] - r*(af[k][i]*(ufr-uR[i]) - af[k][i-1]*(ufl-uR[i-1])) : 0.0;
			const double Nk = kapR[k][i]*(pr-val) - kapL[k][i]*(val-pl) + (AVOL ? V[i].ei[k] - r*(ef[k][i]*(ufr-uR[i]) - ef[k][i-1]*(ufl-uR[i-1])) + val*atr
			                       : V[i].ei[k] - r*(ef[k][i]*(ufr-EXPLE*uR[i]) - ef[k][i-1]*(ufl-EXPLE*uR[i-1])) - r*as*val*du + val*as);
			if(!(Nk > 0.0)) return -1.0;                                    // 能量被搬空: 体积为 0, 记为 G < 0 (压力取高了)
			g += (GAM[k]-1.0)/(GAM[k]*(val+PIN[k]))*Nk; }
		return g; };
	for(int sw=0; sw<SWEEPS; sw++){ double chg = 0.0;
		for(int i=NG;i<N+NG;i++){
			double a = p[i] > 0.0 ? p[i] : 1.0, ga = Gloc(i,a), b = a, gb = ga; int n = 0;
			if(ga > 0.0){ while(gb > 0.0 && n++ < 200 && b < 1e13){ a = b; ga = gb; b *= 2.0; gb = Gloc(i,b); } }
			else        { while(ga < 0.0 && n++ < 200 && a > 1e-9){ b = a; gb = ga; a *= 0.5; ga = Gloc(i,a); } }
			if(!(ga >= 0.0 && gb <= 0.0)) continue;
			for(int it=0; it<200 && (b-a) > 1e-12*b; it++){ double m = 0.5*(a+b), gm = Gloc(i,m); if(gm > 0.0){ a = m; } else { b = m; } }
			double pn = 0.5*(a+b); chg = std::max(chg, std::fabs(pn-p[i])/(pn+PIN[0])); p[i] = pn;
		}
		if(chg < 1e-10) break; }
	if(!eval(p, true)){ g_relax_fail++; if(getenv("DBG") && g_newton_steps<2) fprintf(stderr,"initial point inadmissible\n"); }
	else for(int it=0, itmax = getenv("NITER") ? atoi(getenv("NITER")) : 60; it<=itmax; it++){ if(it == itmax){ if(getenv("NITER")) cvg = true; break; }
		double gmax = 0.0; for(int i=NG;i<N+NG;i++) gmax = std::max(gmax, std::fabs(G[i]));
		if(gmax < 1e-12){ cvg = true; break; }
		if(getenv("DBGN") && g_newton_steps>=atoi(getenv("DBGN"))){ int iw=NG; for(int i=NG;i<N+NG;i++) if(std::fabs(G[i])>=gmax) iw=i; fprintf(stderr,"  call %ld dt %.3e it %d gmax %.3e @%d p %.4e eq %d useq %d a1 %.4f\n",g_newton_steps,dt,it,gmax,iw-NG,p[iw],V[iw].eq,useq[iw],V[iw].a1); }
		for(int i=NG;i<N+NG;i++) dp[i] = -G[i];
		for(int i=NG+1;i<N+NG;i++){ double w = lo[i]/di[i-1]; di[i] -= w*up[i-1]; dp[i] -= w*dp[i-1]; }     // Thomas
		dp[N+NG-1] /= di[N+NG-1]; for(int i=N+NG-2;i>=NG;i--) dp[i] = (dp[i] - up[i]*dp[i+1])/di[i];
		double lam = 1.0; bool ok = false;
		for(int bt=0; bt<40; bt++){   // 可容许域内回缩 (p > 0, N_k > 0)
			for(int i=NG;i<N+NG;i++) q[i] = p[i] + lam*dp[i];
			bool pos = true; for(int i=NG;i<N+NG;i++) if(!(q[i] > 0.0)){ pos = false; break; }
			if(pos && eval(q, true)){ ok = true; break; } lam *= 0.5; }
		if(getenv("DBGN") && g_newton_steps>=atoi(getenv("DBGN"))) fprintf(stderr,"     lam %.3e ok %d\n",lam,(int)ok);
		if(!ok){ eval(p, false); break; }
		p = q;
	}
	g_newton_steps++; if(!cvg) g_newton_fail++;
	for(int g=0;g<NG;g++){ p[g] = p[NG]; p[N+NG+g] = p[N+NG-1]; }
	for(int i=NG-1;i<N+NG;i++){ int a = std::max(i,NG), b = std::min(i+1,N+NG-1); uf[i] = ubf[i] - r*irf[i]*(p[b]-p[a]); setface(i, uf[i]-uR[i]); }
	std::vector<Cell> W = V;
	for(int i=NG;i<N+NG;i++){
		Cell& c = W[i]; const Cell& v = V[i];
		// 面压力: 由面两侧加速度连续 (p_f - p_i)/rho_i = (p_{i+1} - p_f)/rho_{i+1} 得 p_f = (rho_{i+1} p_i + rho_i p_{i+1})/(rho_i + rho_{i+1}); 密度均匀时即算术平均.
		//  算术平均的面压力在密度梯度大的弥散界面里把重侧的压降也摊给轻格, 轻格的加速度被放大 ~rho_重/rho_轻 (水-气激波管里气侧速度过冲).
		auto pfc = [&](int f)->double{ return PFACE ? (rho[f+1]*p[f] + rho[f]*p[f+1])/(rho[f]+rho[f+1]) : 0.5*(p[f]+p[f+1]); };
		c.mu = v.mu - r*(pfc(i)-pfc(i-1));
		auto Hup = [&](int f)->double{ return (ef[0][f] + ef[1][f])*(uf[f]-EXPLE*uR[f]) + pfc(f)*uf[f]; };                       // (rho e)_up + p_f
		c.E  = v.E - r*(Hup(i) - Hup(i-1)) + KA[i]*(p[i+1]-p[i]) - KA[i-1]*(p[i]-p[i-1]);                    // 守恒的总能更新 (B-B (32b) 型)
		{ const double duR = uf[i]-uR[i], duL = uf[i-1]-uR[i-1];             // 质量增量通量 m_k,up du_f, 连同它带走的动量、动能、生成焓 (面通量, 守恒)
		  for(int k=0;k<2;k++){ const double FR = mf[k][i]*duR, FL = mf[k][i-1]*duL;
			const double tR = UNEW ? uf[i] : uup[i], tL = UNEW ? uf[i-1] : uup[i-1];   // 增量质量通量携带的速度: 1 = 新时间层的面速度 u_f (均匀加速度场下 u^{n+1} 保持均匀的充要条件), 0 = 旧速度
			c.m[k] = v.m[k] - r*(FR - FL); c.mu -= r*(tR*FR - tL*FL);
			c.E -= r*((QQ[k] + 0.5*tR*tR)*FR - (QQ[k] + 0.5*tL*tL)*FL); } }
		if(cvg && useq[i]){ double Ee, y = c.m[1]/(c.m[0]+c.m[1]), T; const double rn = c.m[0]+c.m[1]; if(useq[i] == 2 ? Eeq(rn, p[i], Ee, &y, &T) : EpT(rn, p[i], y, Ee, &T)){ c.m[1] = y*rn; c.m[0] = (1.0-y)*rn;
			c.a1 = std::min(std::max(c.m[0]*sgv(0,p[i],T), ALPHA_EPS), 1.0-ALPHA_EPS); for(int k=0;k<2;k++) c.ei[k] = alpha(c,k)*(p[i] + GAM[k]*PIN[k])/(GAM[k]-1.0); } }
		else if(cvg){ double a1 = al[0][i]/(al[0][i]+al[1][i]); c.a1 = std::min(std::max(a1, ALPHA_EPS), 1.0-ALPHA_EPS);
			for(int k=0;k<2;k++) c.ei[k] = alpha(c,k)*(p[i] + GAM[k]*PIN[k])/(GAM[k]-1.0); }
	}
	U = W;
	for(int i=NG;i<N+NG;i++){ relax_reinit(U[i]); thermo_relax(U[i]); }
}

int main(int argc, char** argv){
	const int mode = argc > 1 ? atoi(argv[1]) : 1; const int N = argc > 2 ? atoi(argv[2]) : 1000;
	const double CFL = argc > 3 ? atof(argv[3]) : (mode ? 0.4 : 0.5); const char* out = argc > 4 ? argv[4] : "prof.dat";
	const int picard = 0; if(argc > 5) SWEEPS = atoi(argv[5]); if(argc > 6) ORDER = atoi(argv[6]); if(getenv("CENTRAL")) CENTRAL = atoi(getenv("CENTRAL")); if(getenv("EXPLE")) EXPLE = atoi(getenv("EXPLE")); if(getenv("IRF")) IRF = atoi(getenv("IRF")); if(getenv("PFACE")) PFACE = atoi(getenv("PFACE")); if(getenv("UNEW")) UNEW = atoi(getenv("UNEW")); if(getenv("AVOL")) AVOL = atoi(getenv("AVOL")); if(getenv("BNDCHK")) BNDCHK = atoi(getenv("BNDCHK")); if(getenv("QEGLOB")) QEGLOB = atoi(getenv("QEGLOB")); if(getenv("UREF")) UREF = atoi(getenv("UREF")); if(getenv("RANN")) RANN = atoi(getenv("RANN")); if(getenv("FLAT")) FLAT = atoi(getenv("FLAT")); if(getenv("FLATA")) FLATA = atoi(getenv("FLATA")); if(getenv("FSC")) FSC = atof(getenv("FSC")); if(getenv("K1")) FK1 = atof(getenv("K1")); if(getenv("SLOC")) SLOC = atoi(getenv("SLOC")); if(getenv("DONOR")) DONOR = atoi(getenv("DONOR")); if(getenv("NRAN")) NRAN = atoi(getenv("NRAN"));
	// 论文 4.5 空化算例: 1 m 管, 水 rho = 1000, p = 1e5 Pa, 处处 alpha_air = 1e-2, x < 0.5: u = -100, x > 0.5: u = +100, t = 1.85 ms
	double Lx = 1.0, tend = 1.85e-3, p0 = 1.0e5, rw = 1000.0, ra = 1.0, aair = 1.0e-2, U0 = 100.0;     // 空气密度论文 4.5 未给出, 取 1 kg/m^3
	const int icase = getenv("CASE") ? atoi(getenv("CASE")) : 0;
	if(icase == 1 || icase == 2){   // Pelanti-Shyue 2014 水空化管 (只做力学松弛这一档): 液态水/水蒸气 SG, alpha_v = 1e-2, p = 1e5, rho_l = 1150, rho_v = 0.63
		GAM[0] = 2.35; PIN[0] = 1.0e9; QQ[0] = -1167e3; CV[0] = 1816.0; QP[0] = 0.0;          // 液态水 (P-S 2014 Table 1)
		GAM[1] = 1.43; PIN[1] = 0.0;   QQ[1] = 2030e3;  CV[1] = 1040.0; QP[1] = -23.4e3;      // 水蒸气
		{ const double T0 = 354.728; rw = (p0+PIN[0])/((GAM[0]-1.0)*CV[0]*T0); ra = (p0+PIN[1])/((GAM[1]-1.0)*CV[1]*T0); }
		if(getenv("RELAX")) RELAX = atoi(getenv("RELAX"));
		if(icase == 1){ U0 = 2.0; tend = 3.2e-3; } else { U0 = 500.0; tend = 0.58e-3; } }
	if(icase == 3){   // P-S 2014 5.2 十二烷液-汽激波管: x0 = 0.75, 左: 液体 p = 1e8 (alpha_v = 1e-8), 右: 蒸汽 p = 1e5 (alpha_l = 1e-8), rho_v = 2, rho_l = 500, t = 473 mus
		GAM[0] = 2.35;  PIN[0] = 4.0e8; QQ[0] = -775.269e3; CV[0] = 1077.7;  QP[0] = 0.0;
		GAM[1] = 1.025; PIN[1] = 0.0;   QQ[1] = -237.547e3; CV[1] = 1956.45; QP[1] = -24.4e3;
		if(getenv("RELAX")) RELAX = atoi(getenv("RELAX")); tend = 473e-6; EPSI = 1e-4; }
	if(getenv("EPSI")) EPSI = atof(getenv("EPSI")); RAMP = (icase == 3 || icase == 5 || icase == 6); if(icase == 7 || icase == 8){ GAM[0] = GAM[1] = 1.4; PIN[0] = PIN[1] = 0.0; tend = icase == 7 ? 0.2 : 5e-4; RAMP = 1; } if(icase == 9){ Lx = 3.0; tend = 2e-3; WALLL = 1; } if(icase == 12 || icase == 13){ GAM[0] = 2.43; PIN[0] = 5.3e9; GAM[1] = 1.62; PIN[1] = 141e9; tend = icase == 12 ? 80e-6 : 34e-6; RAMP = 1; }   /* SPB 5.1 / 5.2 环氧-尖晶石混合物激波管 */ if(icase == 20){ GAM[0]=GAM[1]=1.4; PIN[0]=PIN[1]=0.0; tend=5e-4; RAMP=0; } if(icase == 21){ GAM[0]=4.4; PIN[0]=6e8; GAM[1]=1.4; PIN[1]=0.0; tend=1e-4; RAMP=0; }   /* 20, 21: 光滑算例 (时间阶检验) */ if(icase == 4) tend = 2.79e-3; if(icase == 5) tend = 240e-6; if(icase == 6) tend = 8.3e-6; if(getenv("TEND")) tend = atof(getenv("TEND")); if(getenv("RAMP")) RAMP = atoi(getenv("RAMP"));
	if(getenv("TEND")) tend = atof(getenv("TEND"));
	const double dx = Lx/N; std::vector<Cell> U(N+2*NG);
	for(int i=0;i<N+2*NG;i++){ double x = (i-NG+0.5)*dx; Cell& c = U[i]; double u = x < 0.5 ? -U0 : U0;
		double pc0 = p0;
		c.a1 = 1.0-aair; c.m[0] = c.a1*rw; c.m[1] = aair*ra;
		if(icase == 7 || icase == 8){ u = 0.0; const bool L = x < 0.5; pc0 = L ? (icase == 7 ? 1.0 : 1.0e6) : (icase == 7 ? 0.1 : 1.0e5); const double rr = L ? (icase == 7 ? 1.0 : 10.0) : (icase == 7 ? 0.125 : 1.0);   // 7: Sod; 8: 压比 10 的强一点的激波管 (单一理想气体, 两"相"同 EOS, alpha = 0.5)
			c.a1 = 0.5; c.m[0] = 0.5*rr; c.m[1] = 0.5*rr; }
		if(icase == 9){ u = -200.0; const bool L = x < 0.9; c.a1 = L ? 1.0-1e-6 : 1e-6; c.m[0] = c.a1*1000.0; c.m[1] = (1.0-c.a1)*100.0; }   // 空化 RM 的平面类比: 左壁, 水 (x<0.9) 和气都以 -200 m/s 撞壁
		if(icase == 4){ u = 100.0; const bool L = x < 0.5; c.a1 = L ? 1.0-1e-8 : (getenv("C4R") ? atof(getenv("C4R")) : 1e-8); c.m[0] = c.a1*1000.0; c.m[1] = (1.0-c.a1)*10.0; }                 // SPB 4.1 界面平移
		if(icase == 5 || icase == 6){ u = 0.0; const bool L = x < 0.75; pc0 = L ? (icase == 5 ? 1.0e9 : 1.0e12) : 1.0e5;                                   // SPB 4.3.1 / 4.3.2 水-气激波管
			c.a1 = L ? 1.0-1e-6 : 1e-6; c.m[0] = c.a1*1000.0; c.m[1] = (1.0-c.a1)*(icase == 5 ? 1.0 : 10.0); }
		if(icase == 12 || icase == 13){ u = 0.0; pc0 = x < 0.6 ? (icase == 12 ? 1.0e10 : 2.0e11) : 1.0e5; c.a1 = 0.5954; c.m[0] = c.a1*1185.0; c.m[1] = (1.0-c.a1)*3622.0; }
		if(icase == 3){ u = 0.0; const bool L = x < 0.75; pc0 = L ? 1.0e8 : 1.0e5; c.a1 = L ? 1.0-1e-8 : 1e-8; c.m[0] = c.a1*500.0; c.m[1] = (1.0-c.a1)*2.0; }
		if(icase == 20 || icase == 21){ const double gs=std::exp(-(x-0.5)*(x-0.5)/0.01); u = 50.0*gs; c.a1 = 0.5+0.2*std::sin(6.283185307179586*x); c.m[0]=c.a1*(icase==20?2.0:1000.0); c.m[1]=(1.0-c.a1)*1.0; pc0 = 1e5*(1.0+0.2*gs); } double rho = c.m[0]+c.m[1]; c.mu = rho*u;
		for(int k=0;k<2;k++) c.ei[k] = alpha(c,k)*(pc0 + GAM[k]*PIN[k])/(GAM[k]-1.0);
		c.E = c.ei[0]+c.ei[1] + 0.5*rho*u*u + c.m[0]*QQ[0]+c.m[1]*QQ[1]; }
	auto totals = [&](double& M1, double& M2, double& P, double& E){ M1=M2=P=E=0; for(int i=NG;i<N+NG;i++){ M1+=U[i].m[0]*dx; M2+=U[i].m[1]*dx; P+=U[i].mu*dx; E+=U[i].E*dx; } };
	double M10,M20,P0,E0; totals(M10,M20,P0,E0);
	double t = 0.0; long nstep = 0; double bM1=0,bM2=0,bE=0;                  // 边界通量累计 (透射边界有流出, 守恒检查要扣掉)
	while(t < tend*(1.0-1e-12)){
		double smax = 0.0, amax = 0.0;
		for(int i=NG;i<N+NG;i++){ double u = std::fabs(U[i].mu/(U[i].m[0]+U[i].m[1])); smax = std::max(smax,u); amax = std::max(amax, u + cfrozen(U[i])); }
		double dt = CFL*dx/(mode ? std::max(smax,1e-12) : amax);
		bool ramping = false;
		if(mode && RAMP){ static double dtprev = -1.0; if(dtprev < 0.0) dtprev = CFL*dx/amax/1.2; if(1.2*dtprev < dt){ dt = 1.2*dtprev; ramping = true; } dtprev = dt; }   // 从静止起算: 首步取声学步长, 之后每步至多放大 1.2 倍直到物质 CFL
		dt = std::min(dt, tend - t);
		auto bc = [&](std::vector<Cell>& X){ for(int g=0;g<NG;g++){ X[NG-1-g] = WALLL ? X[NG+g] : X[NG]; if(WALLL) X[NG-1-g].mu = -X[NG-1-g].mu; X[N+NG+g] = X[N+NG-1]; } };   // WALLL: 左端固壁 (镜像), 只给显式参考用
		bc(U);
		auto admissible = [&](const std::vector<Cell>& X)->bool{ for(int i=NG;i<N+NG;i++){ const Cell& c = X[i];
			if(!(c.a1 >= ALPHA_EPS && c.a1 <= 1.0-ALPHA_EPS && c.m[0] > 0.0 && c.m[1] > 0.0)) return false;   /* 与显式半步里对 alpha 的钳位同一个界: 否则钳位只改 alpha 不改相内能, 破坏均匀 p */
			for(int k=0;k<2;k++) if(!(pk(c,k) + PIN[k] > 0.0)) return false; } return true; };
		bc(U);
		if(!mode){   // 论文显式方法: 一阶 = 单步; 二阶 = MUSCL + Heun
			auto one = [&](){ bc(U); step_explicit(U, N, dt, dx); for(int i=NG;i<N+NG;i++){ relax_reinit(U[i]); thermo_relax(U[i]); } };
			if(ORDER < 2) one();
			else{ std::vector<Cell> U0 = U; one(); one();
				for(int i=NG;i<N+NG;i++){ Cell& c = U[i]; const Cell& o = U0[i]; c.a1=0.5*(c.a1+o.a1); c.mu=0.5*(c.mu+o.mu); c.E=0.5*(c.E+o.E);
					for(int k=0;k<2;k++){ c.m[k]=0.5*(c.m[k]+o.m[k]); c.ei[k]=0.5*(c.ei[k]+o.ei[k]); } relax_reinit(c); thermo_relax(c); } }
		}else if(ORDER < 2){ std::vector<Cell> Qs = U; step_imex(U, Qs, N, dt, dx, picard); }
		else{
			// 半隐式 IMEX-RK, 刚性精确的 (2,2,2) 表 (B-B 用的那张): gamma = 1 - 1/sqrt(2)
			//   隐式 a = [[g,0],[1-g,g]], 显式 at = [[0,0],[1/(2g),0]], b = bt = (1-g, g)  => 末级即 U^{n+1}
			//   级 i:  Qe = U^n + Dt sum_{j<i} at_ij K_j,  Qs = U^n + Dt sum_{j<i} a_ij K_j,  U_i = 级算子(Qe, Qs, a_ii Dt),  K_i = (U_i - Qs)/(a_ii Dt)
			const double g = 1.0 - 1.0/std::sqrt(2.0), at21 = 1.0/(2.0*g), a21 = 1.0 - g;
			const std::vector<Cell> Un = U;
			std::vector<Cell> U1 = Un; step_imex(U1, Un, N, g*dt, dx, picard);                       // 级 1: Qe = Qs = U^n
			std::vector<Cell> Qe = Un, Qs = Un;
			auto axpy = [&](std::vector<Cell>& X, double w){   // X = U^n + w * (U1 - U^n)   (K_1 Dt = (U1 - U^n)/g)
				for(int i=NG;i<N+NG;i++){ Cell& c = X[i]; const Cell& o = Un[i]; const Cell& n1 = U1[i];
					c.a1=o.a1+w*(n1.a1-o.a1); c.mu=o.mu+w*(n1.mu-o.mu); c.E=o.E+w*(n1.E-o.E);
					for(int k=0;k<2;k++){ c.m[k]=o.m[k]+w*(n1.m[k]-o.m[k]); c.ei[k]=o.ei[k]+w*(n1.ei[k]-o.ei[k]); } } };
			axpy(Qe, at21/g); axpy(Qs, a21/g); bc(Qe); bc(Qs);
			// Qe 只用来算显式通量: 逐格把外推量限到可容许 (theta in [0,1]); 不影响守恒. Qs 是更新基底, 不能逐格限 -> 不可容许则本步退回一阶.
			double thmin = 1.0;
			for(int i=NG;i<N+NG;i++){ auto okc=[&](const Cell& c)->bool{ if(!(c.a1 >= ALPHA_EPS && c.a1 <= 1.0-ALPHA_EPS && c.m[0] > 0.0 && c.m[1] > 0.0)) return false;   /* 与显式半步里对 alpha 的钳位同一个界: 否则钳位只改 alpha 不改相内能, 破坏均匀 p */ for(int k=0;k<2;k++) if(!(pk(c,k) + PIN[k] > 0.0)) return false; return true; };
				double th = 1.0; int n = 0; while(!okc(Qe[i]) && n++ < 30){ th *= 0.5; if(n==30) th = 0.0; const double w = th*at21/g; Cell& c = Qe[i]; const Cell& o = Un[i]; const Cell& n1 = U1[i];
					c.a1=o.a1+w*(n1.a1-o.a1); c.mu=o.mu+w*(n1.mu-o.mu); c.E=o.E+w*(n1.E-o.E); for(int k=0;k<2;k++){ c.m[k]=o.m[k]+w*(n1.m[k]-o.m[k]); c.ei[k]=o.ei[k]+w*(n1.ei[k]-o.ei[k]); } }
				if(th < 1.0) g_qelim++; thmin = std::min(thmin, th); }
			if(QEGLOB && thmin < 1.0){ axpy(Qe, thmin*at21/g); }   // 全场同一个 theta: 等价于本步把显式表的 c2 缩小, 通量仍在一个自洽的中间时刻状态上算
			bc(Qe);
			bool ok2 = !getenv("FORCE_RK1") && admissible(Qs);
			if(getenv("FO_UNTIL") && t < atof(getenv("FO_UNTIL"))) ok2 = false;
			if(RANN && (ramping || nstep < NRAN)) ok2 = false;   // Rannacher 式起步: 间断初值的起始阶段用后向 Euler (L 稳定, 单调), 之后才切到二阶
			if(ok2){ U = Qe; g_trbad = 0; step_imex(U, Qs, N, g*dt, dx, picard); if(BNDCHK && g_trbad) ok2 = false; }
			if(!ok2){ g_rk1++; U = Un; const int o = ORDER; if(!getenv("FB_MUSCL")) ORDER = 1; step_imex(U, Un, N, dt, dx, picard); ORDER = o; }   // 退回步在空间上也取一阶 (MUSCL + 后向 Euler 在大 CFL 下不稳)     // 外推级态不可容许 (at21/g = 5.8 倍外推): 本步退回一阶
		}
		t += dt; nstep++;
		if(getenv("PDEV")){ int iw=NG; double dv=0; for(int i=NG;i<N+NG;i++){ double e=std::fabs(pmix(U[i])-1e5); if(e>dv){dv=e;iw=i;} }
			if(dv > atof(getenv("PDEV"))){ fprintf(stderr,"step %ld worst cell %d dp %.3e | ",nstep,iw-NG,dv); for(int j=iw-3;j<=iw+3;j++) fprintf(stderr,"[%d a1 %.3e p1 %.6e p2 %.6e u %.6f] ",j-NG,U[j].a1,pk(U[j],0),pk(U[j],1),U[j].mu/(U[j].m[0]+U[j].m[1])); fprintf(stderr,"rk1 %ld qelim %ld\n",g_rk1,g_qelim); } }
		if(getenv("TRACE")){ double um=0,pmn=1e99,pmx=0; int im=0; for(int i=NG;i<N+NG;i++){ double u=std::fabs(U[i].mu/(U[i].m[0]+U[i].m[1])); if(u>um){um=u;im=i-NG;} double pp=pmix(U[i]); pmn=std::min(pmn,pp); pmx=std::max(pmx,pp);} 
			fprintf(stderr,"step %ld t %.4e dt %.3e umax %.4f @%d pmin %.4e pmax %.4e rk1 %ld qelim %ld nwfail %ld relaxfail %ld\n",nstep,t,dt,um,im,pmn,pmx,g_rk1,g_qelim,g_newton_fail,g_relax_fail); }
		if(getenv("PROG") && nstep%500==0){ fprintf(stderr,"step %ld t %.4e dt %.3e pTg cells %ld\n",nstep,t,dt,g_pTg_cells); }
	}
	(void)bM1; (void)bM2; (void)bE;
	double M1,M2,P,E; totals(M1,M2,P,E);
	if(RELAX) printf("closure fallbacks %ld | relax level %d: pT failures %ld, pTg applied in %ld cell-updates, Tsat(1e5) = %.2f K\n", g_clfail, RELAX, g_pT_fail, g_pTg_cells, Tsat(1e5));
	FILE* f = fopen(out,"w");
	for(int i=NG;i<N+NG;i++){ const Cell& c = U[i]; double rho = c.m[0]+c.m[1];
		const double p = pmix(c), cw = std::sqrt(1.0/(rho*(c.a1/(GAM[0]*(p+PIN[0])) + (1.0-c.a1)/(GAM[1]*(p+PIN[1])))));   // 列 7-10: rho_1, rho_2, 冻结声速, Wood 声速
		fprintf(f, "%.8e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", (i-NG+0.5)*dx, p, c.mu/rho, rho, 1.0-c.a1, c.m[1]/rho, c.m[0]/c.a1, c.m[1]/(1.0-c.a1), cfrozen(c), cw); }
	fclose(f);
	printf("mode %d order %d N %d CFL %g steps %ld (RK1 fallback %ld, Qe limited cells %ld) | Newton fail %ld/%ld  relax fail %ld  p-floor hits %ld | d(mass_w)/M %.2e d(mass_a)/M %.2e dE/E %.2e (含透射边界流出)\n",
	       mode, ORDER, N, CFL, nstep, g_rk1, g_qelim, g_newton_fail, g_newton_steps, g_relax_fail, g_pfloor, (M1-M10)/M10, (M2-M20)/M20, (E-E0)/E0);
	return 0;
}
