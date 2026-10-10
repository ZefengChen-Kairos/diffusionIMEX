// =====================================================================================================================
//  SPB6_explicit.cpp
//  全显式参考格式（docs/derivation/A_model_EOS_explicit.md）：SPB 2009 的 HLLC + Godunov 更新（SPB 3.1-3.2），
//  力学松弛 + 用总能重初始化（SPB 3.3-3.4），P-S 2014 / Pelanti 2022 的热、化学松弛（瞬时或有限速率）。
//  二阶：原始变量 MUSCL（minmod）+ Heun；瞬时松弛（投影）每级之后都做，有限速率松弛每步只在平均之后做一次（A3）。
//
//  编译: g++ -O2 -std=c++17 SPB6_explicit.cpp -o spb6_ex
//  运行: ./spb6_ex key=value ...   键见 spb6_common.h 的 parse_common（显式默认声学 CFL 0.4）
// =====================================================================================================================
#include "spb6_common.h"

// ==================================================================================================== 显式参考（HLLC + Godunov + 松弛）
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

// 显式步的局部松弛（A3）：力学（SPB 原文：Gibbs 关系 + 重初始化）然后热/化学（P-S / Pelanti，用统一闭合 s = 0）。
// 瞬时松弛是投影（幂等），Heun 的每一级之后都做；有限速率（mu, theta, nu 有限）是对时间的积分，每个时间步只做一次（finite = true，
// tau = dt）。若每级都做，Heun 的一步对 y' = -lambda y 给出 R(1/2 (y + R(R(y)))) = (1 - 2 lambda dt) y + O(dt^2)：速率放大约 2 倍。
static bool relax_cell_explicit(double* Q, double dtloc, bool finite) {
	const double rho = Q[IM1] + Q[IM2], rhoe = Q[IEN] - 0.5 * (Q[IMU] * Q[IMU] + Q[IMV] * Q[IMV]) / rho;
	double a1 = Q[IA1], m1 = Q[IM1], m2 = Q[IM2], e1 = Q[IE1], e2 = Q[IE2];
	const int mech = std::isinf(mu_rate) ? CL_P : CL_MU;
	if (mech == CL_P || finite) {
		ClIn in{}; in.atr[0] = a1; in.atr[1] = 1 - a1; in.mtr[0] = m1; in.mtr[1] = m2; in.etr[0] = e1; in.etr[1] = e2;
		in.s = 0; in.rhoeC = rhoe; in.tau = dtloc; in.pguess = std::max(pbar_cell_q(Q), 1.0);
		in.type = mech; in.nodelta = (mech == CL_P);
		ClOut o; if (!closure(in, o)) { st_cl_fail++; return false; }
		a1 = o.a[0]; m1 = o.m[0]; m2 = o.m[1]; e1 = o.e[0]; e2 = o.e[1];
		if (in.nodelta) {   // SPB 3.4 重初始化：alpha 固定，用守恒总能由 (A2.3) 重算 p，再重置相内能
			const double p = p_mix(a1, m1, m2, rhoe);
			if (!(p + PH[0].pi > 0 && p + PH[1].pi > 0)) { st_cl_fail++; return false; }
			e1 = eps_amp(0, a1, m1, p); e2 = eps_amp(1, 1 - a1, m2, p);
		}
	}
	int ty = closure_type_of(a1, m1, m2, e1, e2, rhoe);
	if (!finite) { if (ty == CL_NU) ty = CL_PT; else if (ty == CL_THETA) ty = mech; }   // 级内只做瞬时部分（nu 有限时 theta = inf）
	if (ty != mech) {
		ClIn t2{}; t2.atr[0] = a1; t2.atr[1] = 1 - a1; t2.mtr[0] = m1; t2.mtr[1] = m2; t2.etr[0] = e1; t2.etr[1] = e2;
		t2.s = 0; t2.rhoeC = rhoe; t2.tau = dtloc; t2.type = ty; t2.pguess = std::max(p_mix(a1, m1, m2, rhoe), 1.0);
		ClOut o2; if (!closure(t2, o2)) { st_cl_fail++; } else { a1 = o2.a[0]; m1 = o2.m[0]; m2 = o2.m[1]; e1 = o2.e[0]; e2 = o2.e[1]; }
	}
	Q[IA1] = a1; Q[IM1] = m1; Q[IM2] = m2; Q[IE1] = e1; Q[IE2] = e2;
	return true;
}

static void reconstruct_prim_minmod(const Grid& U, std::vector<Prim>& WL, std::vector<Prim>& WR, int dir, bool& ok);

static bool explicit_step(Grid& U, double dtl, bool finite) {
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
		if (!relax_cell_explicit(Q, dtl, finite)) return false;
		for (int k = 0; k < NUM_VARS; k++) V(i, j, k) = Q[k];
	}
	U = V; bound_cond_u(U); return true;
}

// 原始变量的 minmod 重构（SPB 附录 C 也用原始变量）
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

// ==================================================================================================== main
int main(int argc, char** argv) {
	CFL = 0.4;   // 声学 CFL
	CommonArgs A;
	for (int a = 1; a < argc; a++) {
		char key[64]; char val[256]; if (sscanf(argv[a], "%63[^=]=%255s", key, val) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); return 1; }
		if (!parse_common(key, val, atof(val), A)) { fprintf(stderr, "unknown key %s\n", key); return 1; }
	}
	setup_common(A);
	Grid U(nxCells, nyCells); init_case(U);
	double M0[3]; conserved_totals(U, M0);
	auto t0 = std::chrono::steady_clock::now();
	t = 0;
	while (t < time_stop * (1 - 1e-12)) {
		dt = std::min(CFL * caldt(U, true), time_stop - t);
		if (!(dt > 1e-12 * time_stop)) { fprintf(stderr, "FAILED: time step underflow (dt = %.3e) at t = %.6e (step %ld)\n", dt, t, st_steps); output(U, A.out.c_str()); return 2; }
		bool done = false;
		for (int attempt = 0; attempt < 30 && !done; attempt++) {
			Grid V = U;
			if (SpaceOrder < 2) done = explicit_step(V, dt, true);
			else { Grid W = V; done = explicit_step(V, dt, false) && explicit_step(V, dt, false);   // Heun；有限速率松弛只在平均之后做一次（A3）
				if (done) { for (size_t q = 0; q < V.data.size(); q++) V.data[q] = 0.5 * (V.data[q] + W.data[q]);
					for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++) { double Q[NUM_VARS]; for (int k = 0; k < NUM_VARS; k++) Q[k] = V(i, j, k); if (!relax_cell_explicit(Q, dt, true)) { done = false; break; } for (int k = 0; k < NUM_VARS; k++) V(i, j, k) = Q[k]; }
					bound_cond_u(V); } }
			if (done) U = V;
			if (!done) { st_reject++; dt *= 0.5; if (VERBOSE) fprintf(stderr, "reject: t %.4e new dt %.3e\n", t, dt); }
		}
		if (!done) { fprintf(stderr, "FAILED at t = %.6e (step %ld)\n", t, st_steps); output(U, A.out.c_str()); return 2; }
		t += dt; st_steps++;
		if (VERBOSE > 1 && st_steps % 50 == 0) fprintf(stderr, "step %ld t %.4e dt %.3e\n", st_steps, t, dt);
	}
	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double M1[3]; conserved_totals(U, M1);
	output(U, A.out.c_str());
	printf("case %d mode 0 order %d N %d CFL %g | steps %ld rejects %ld | closure fails %ld | dM1/M %.2e dM2/M %.2e dE/E %.2e | wall %.2f s\n",
		TESTCASE, SpaceOrder, nxCells, CFL, st_steps, st_reject, st_cl_fail, (M1[0] - M0[0]) / M0[0], (M1[1] - M0[1]) / M0[1], (M1[2] - M0[2]) / M0[2], wall);
	return 0;
}
