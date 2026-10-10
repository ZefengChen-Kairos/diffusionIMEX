// =====================================================================================================================
//  SPB6_explicit.cpp
//  全显式参考格式（docs/derivation/A_model_EOS_explicit.md）：SPB 2009 的 HLLC + Godunov 更新（SPB 3.1-3.2），
//  力学松弛（SPB 3.3，p_I 取松弛后的 p），守恒总能与相内能之和的差（能量亏损）在松弛闭合里按质量分数分给两相
//  （代替 SPB 3.4 的重初始化；P-S 2014 §2.2 的相总能形式、§3.1 的一致性 (ii)），
//  P-S 2014 / Pelanti 2022 的热、化学松弛（瞬时或有限速率）。
//
//  与 SPB6_TV_IMEX.cpp 一样由 reference/code/2DEuler_TV.cpp 修改而来（全局设置、Grid/Array、TVDReconstruction、Update、
//  u_initialization、bound_cond_u、caldt、bound_cond_array、VTK 输出、main 的时间循环），状态方程、闭合、算例在 spb6_physics.h（两个程序共用）。
//  与 IMEX 的不同：RusanovC_x/y 换成 HLLC_x/y（SPB 3.1.2），重构在原始变量上（SPB 附录 C），PressureSubsystem 换成逐格的
//  RelaxationSubsystem；时间为二阶 Heun：瞬时松弛（投影）每级之后都做，有限速率松弛每步只在平均之后做一次（A3）。
//
//  编译: g++ -O2 -std=c++17 -I/usr/include/eigen3 SPB6_explicit.cpp -o spb6_ex
//  运行: ./spb6_ex key=value ...   键：case N ny cfl order tend out relax epsI mu theta nu verbose vtk
// =====================================================================================================================
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <functional>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <sstream>
#include "spb6_physics.h"


enum Reconstruction{NoRecon, Minmod};

Test test = testSPB_WaterAir;
x_boundary_condition x_bc;
y_boundary_condition y_bc;
int nxCells = 1000;
int nyCells = 1;
const int ghost_cells = 2;

double CFL = 0.4;                     // 声学 CFL
int SpaceOrder = 2;
Reconstruction recon = Minmod;

double t, time_stop, dt;
int timestep = 0;

double x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1;
double dx, dy;

std::string out_file = "out.dat", vtk_dir = "";
long st_steps = 0, st_reject = 0;


struct Grid
{
	std::vector<double> data;
	int nxCells;
	int nyCells;

	Grid(int xcell, int ycell) : data((size_t)(xcell + 2 * ghost_cells) * (ycell + 2 * ghost_cells) * NUM_VARS, 0.0){
		nxCells = xcell;
		nyCells = ycell;
	}

	double operator()(int i, int j, int k)const
	{
		return data[NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int i, int j, int k){
		return data[NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)) + k];
	}
	const double* cell(int i, int j) const { return &data[NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells))]; }
	double* cell(int i, int j) { return &data[NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells))]; }
};

struct Array
{
	std::vector<double> data;
	int nxCells;
	int nyCells;

	Array(int xcell, int ycell, double v = 0.0) : data((size_t)(xcell + 2 * ghost_cells) * (ycell + 2 * ghost_cells), v){
		nxCells = xcell;
		nyCells = ycell;
	}

	double operator()(int i, int j)const
	{
		return data[(size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)];
	}
	double& operator()(int i, int j){
		return data[(size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)];
	}
};

// 面上的原始变量（重构结果；HLLC 需要原始变量与声速）
struct PrimGrid
{
	std::vector<Prim> data;
	int nxCells;
	int nyCells;
	PrimGrid(int xcell, int ycell) : data((size_t)(xcell + 2 * ghost_cells) * (ycell + 2 * ghost_cells)), nxCells(xcell), nyCells(ycell) {}
	const Prim& operator()(int i, int j) const { return data[(size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)]; }
	Prim& operator()(int i, int j) { return data[(size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)]; }
};


using namespace std;
void copy(Grid &state_1, Grid &state_0);
void ucon2upri(Grid &u_pri,Grid &u_con);
void u_initialization(Grid &u);
void bound_cond_u(Grid &u);
double caldt(Grid &u_con, bool acoustic);
bool TVDReconstruction(Grid &u_con, PrimGrid &w, PrimGrid &w_xL, PrimGrid &w_xR, PrimGrid &w_yL, PrimGrid &w_yR);
bool Update(Grid &u_new, Grid &u_con, PrimGrid &w, PrimGrid &w_xL, PrimGrid &w_xR, PrimGrid &w_yL, PrimGrid &w_yR, double h);
bool RelaxationSubsystem(Grid &u_con, double h, bool finite);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri, const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime);


void copy(Grid &state_1, Grid &state_0){
	state_1.data = state_0.data;
}

// 输出用的“原始变量”：0 rho, 1 u, 2 v, 3 p（与 2DEuler 相同），4 alpha2, 5 Y2, 6 T1, 7 T2
void ucon2upri(Grid &u_pri,Grid &u_con){
	for(int i = 0; i < u_con.nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < u_con.nyCells + 2 * ghost_cells; j++){
			double o[15]; cell_diag(u_con.cell(i, j), o);
			u_pri(i,j,0) = o[3]; u_pri(i,j,1) = o[1]; u_pri(i,j,2) = o[2]; u_pri(i,j,3) = o[0];
			u_pri(i,j,4) = o[4]; u_pri(i,j,5) = o[5]; u_pri(i,j,6) = o[10]; u_pri(i,j,7) = o[11];
		}
	}
}


// ==================================================================================================== HLLC（SPB 3.1.2；代替 2DEuler 的 RusanovC_x/y）
struct HFace { double F[5], ustar, ua1, ue[2]; };   // F: m1, m2, normal mom, tangential mom, E; ustar; (u alpha1)*; (eps_k u)*
// 相 Hugoniot (II.8) 的 NASG 形式: e* - e + 1/2 (p* + p)(v* - v) = 0 对 p* 线性
static double hugoniot_p(int k, double p, double v, double vs) {
	const Phase& P = PH[k];
	const double num = (p + P.g * P.pi) * (v - P.b) - P.g * P.pi * (vs - P.b) - 0.5 * (P.g - 1.0) * p * (vs - v);
	const double den = (vs - P.b) + 0.5 * (P.g - 1.0) * (vs - v);
	return num / den;
}
// dir = 0: x 面 (法向 u)；dir = 1: y 面 (法向 v)
static HFace HLLC(const Prim& L, const Prim& R, int dir) {
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


// Godunov 更新（SPB 3.2）：u_new = u_con - h div F；alpha_1 与相内能方程的非守恒项用格内的 alpha_1、alpha_k p_k 与面上的 u*
bool Update(Grid &u_new, Grid &u_con, PrimGrid &w, PrimGrid &w_xL, PrimGrid &w_xR, PrimGrid &w_yL, PrimGrid &w_yR, double h){
	const int g = ghost_cells, nx = nxCells, ny = nyCells, NXT = nx + 2 * g, NYT = ny + 2 * g;
	const double rx = h / dx, ry = h / dy;
	std::vector<HFace> FX((size_t)NXT * NYT), FY;
	for (int j = g; j < ny + g; j++) for (int i = g - 1; i < nx + g; i++) FX[i + (size_t)j * NXT] = HLLC(w_xR(i, j), w_xL(i + 1, j), 0);
	if (ny > 1) { FY.resize((size_t)NXT * NYT); for (int j = g - 1; j < ny + g; j++) for (int i = g; i < nx + g; i++) FY[i + (size_t)j * NXT] = HLLC(w_yR(i, j), w_yL(i, j + 1), 1); }
	copy(u_new, u_con);
	for (int j = g; j < ny + g; j++) for (int i = g; i < nx + g; i++) {
		const Prim& wc = w(i, j);
		const HFace& fp = FX[i + (size_t)j * NXT], & fm = FX[i - 1 + (size_t)j * NXT];
		double dU[NUM_VARS] = {0};
		const double ak[2] = {wc.a1, 1 - wc.a1}, pk[2] = {wc.p1, wc.p2};
		dU[IM1] = -rx * (fp.F[0] - fm.F[0]); dU[IM2] = -rx * (fp.F[1] - fm.F[1]); dU[IMU] = -rx * (fp.F[2] - fm.F[2]); dU[IMV] = -rx * (fp.F[3] - fm.F[3]); dU[IEN] = -rx * (fp.F[4] - fm.F[4]);
		dU[IA1] = -rx * ((fp.ua1 - fm.ua1) - wc.a1 * (fp.ustar - fm.ustar));
		for (int k = 0; k < 2; k++) dU[IE1 + k] = -rx * ((fp.ue[k] - fm.ue[k]) + ak[k] * pk[k] * (fp.ustar - fm.ustar));
		if (ny > 1) {
			const HFace& gp = FY[i + (size_t)j * NXT], & gm = FY[i + (size_t)(j - 1) * NXT];
			dU[IM1] -= ry * (gp.F[0] - gm.F[0]); dU[IM2] -= ry * (gp.F[1] - gm.F[1]); dU[IMV] -= ry * (gp.F[2] - gm.F[2]); dU[IMU] -= ry * (gp.F[3] - gm.F[3]); dU[IEN] -= ry * (gp.F[4] - gm.F[4]);
			dU[IA1] -= ry * ((gp.ua1 - gm.ua1) - wc.a1 * (gp.ustar - gm.ustar));
			for (int k = 0; k < 2; k++) dU[IE1 + k] -= ry * ((gp.ue[k] - gm.ue[k]) + ak[k] * pk[k] * (gp.ustar - gm.ustar));
		}
		for (int k = 0; k < NUM_VARS; k++) u_new(i, j, k) = u_con(i, j, k) + dU[k];
		if (!(u_new(i, j, IA1) > 0 && u_new(i, j, IA1) < 1 && u_new(i, j, IM1) > 0 && u_new(i, j, IM2) > 0)) return false;
	}
	return true;
}


// 原始变量 (alpha1, rho1, rho2, u, v, p1, p2) 的 minmod 重构（SPB 附录 C）；w 为格平均的原始变量
bool TVDReconstruction(Grid &u_con, PrimGrid &w, PrimGrid &w_xL, PrimGrid &w_xR, PrimGrid &w_yL, PrimGrid &w_yR){
	const int g = ghost_cells, NXT = nxCells + 2 * g, NYT = nyCells + 2 * g;
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) if (!prim_of(u_con.cell(i, j), w(i, j))) return false;
	w_xL.data = w.data; w_xR.data = w.data; w_yL.data = w.data; w_yR.data = w.data;
	SpaceOrder = recon == NoRecon ? 1 : 2;
	if (recon == NoRecon) return true;
	auto mm = [](double a, double b) { return a * b > 0 ? (std::fabs(a) < std::fabs(b) ? a : b) : 0.0; };
	for (int dir = 0; dir < (nyCells > 1 ? 2 : 1); dir++) {
		PrimGrid &WL = dir == 0 ? w_xL : w_yL, &WR = dir == 0 ? w_xR : w_yR;
		for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) {
			int im = i, ip = i, jm = j, jp = j;
			if (dir == 0) { im = i - 1; ip = i + 1; if (im < 0 || ip >= NXT) continue; } else { jm = j - 1; jp = j + 1; if (jm < 0 || jp >= NYT) continue; }
			const Prim& c = w(i, j), & a = w(im, jm), & b = w(ip, jp);
			const double cv[7] = {c.a1, c.r1, c.r2, c.u, c.v, c.p1, c.p2}, av[7] = {a.a1, a.r1, a.r2, a.u, a.v, a.p1, a.p2}, bv[7] = {b.a1, b.r1, b.r2, b.u, b.v, b.p1, b.p2};
			double l[7], r[7];
			for (int q = 0; q < 7; q++) { const double sl = mm(cv[q] - av[q], bv[q] - cv[q]); l[q] = cv[q] - 0.5 * sl; r[q] = cv[q] + 0.5 * sl; }
			Prim PL = c, PR = c;
			PL.a1 = l[0]; PL.r1 = l[1]; PL.r2 = l[2]; PL.u = l[3]; PL.v = l[4]; PL.p1 = l[5]; PL.p2 = l[6];
			PR.a1 = r[0]; PR.r1 = r[1]; PR.r2 = r[2]; PR.u = r[3]; PR.v = r[4]; PR.p1 = r[5]; PR.p2 = r[6];
			for (Prim* P : {&PL, &PR}) { const double m1 = P->a1 * P->r1, m2 = (1 - P->a1) * P->r2; P->rho = m1 + m2; P->cf = std::sqrt((m1 * c2_rp(0, P->r1, P->p1) + m2 * c2_rp(1, P->r2, P->p2)) / P->rho); }
			WL(i, j) = PL; WR(i, j) = PR;
		}
	}
	return true;
}


// 初值：算例的定义（区域、状态方程、初值）在 spb6_physics.h（显式与 IMEX 共用）
void u_initialization(Grid &u){
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			double Q[NUM_VARS];
			test_state(test, x, y, Q);
			for(int k = 0; k < NUM_VARS; k++) u(i,j,k) = Q[k];
		}
	}
}


// 透射：复制边界格；固壁：镜像，只把法向动量取负；周期：绕回。一维时 y 向只复制。
void bound_cond_u(Grid &u){
	const int g = ghost_cells, nx = u.nxCells, ny = u.nyCells;
	for(int j = 0; j < ny + 2 * g; j++){
		for(int q = 0; q < g; q++){
			int iL = g - 1 - q, iR = nx + g + q, sL, sR;
			switch(x_bc){
				case x_Periodic: sL = nx + g - 1 - q; sR = g + q; break;
				case x_Reflective: sL = g + q; sR = nx + g - 1 - q; break;
				default: sL = g; sR = nx + g - 1; break;
			}
			for(int k = 0; k < NUM_VARS; k++){ u(iL, j, k) = u(sL, j, k); u(iR, j, k) = u(sR, j, k); }
			if(x_bc == x_Reflective){ u(iL, j, IMU) = -u(iL, j, IMU); u(iR, j, IMU) = -u(iR, j, IMU); }
		}
	}
	if(ny > 1 || y_bc == y_Periodic){
		for(int i = 0; i < nx + 2 * g; i++){
			for(int q = 0; q < g; q++){
				int jL = g - 1 - q, jR = ny + g + q, sL, sR;
				switch(y_bc){
					case y_Periodic: sL = ny + g - 1 - q; sR = g + q; break;
					case y_Reflective: sL = g + q; sR = ny + g - 1 - q; break;
					default: sL = g; sR = ny + g - 1; break;
				}
				for(int k = 0; k < NUM_VARS; k++){ u(i, jL, k) = u(i, sL, k); u(i, jR, k) = u(i, sR, k); }
				if(y_bc == y_Reflective){ u(i, jL, IMV) = -u(i, jL, IMV); u(i, jR, IMV) = -u(i, jR, IMV); }
			}
		}
	}else{
		for(int i = 0; i < nx + 2 * g; i++){
			for(int q = 0; q < g; q++){
				for(int k = 0; k < NUM_VARS; k++){ u(i, g - 1 - q, k) = u(i, g, k); u(i, g + 1 + q, k) = u(i, g, k); }
			}
		}
	}
}


// 时间步：声学速度 |u| + c_f（acoustic = true）。返回 1 / (a_x/dx + a_y/dy)；dt = CFL * 返回值
double caldt(Grid &u_con, bool acoustic){
	double ux = 0, uy = 0; const int g = ghost_cells;
	for(int j = g; j < u_con.nyCells + g; j++){
		for(int i = g; i < u_con.nxCells + g; i++){
			Prim w; if(!prim_of(u_con.cell(i, j), w)) continue;
			const double c = acoustic ? w.cf : 0.0;
			ux = std::max(ux, std::fabs(w.u) + c); uy = std::max(uy, std::fabs(w.v) + c);
		}
	}
	double rate = ux / dx + (u_con.nyCells > 1 ? uy / dy : 0.0);
	return rate > 0 ? 1.0 / rate : INF;
}


// ==================================================================================================== 松弛（代替 2DEuler 的 PressureSubsystem）
// 局部松弛（A3）：力学（SPB 3.3 的 Gibbs 关系，能量亏损按质量分数分给两相）然后热/化学（P-S / Pelanti，用统一闭合 s = 0）。
// 瞬时松弛是投影（幂等），Heun 的每一级之后都做；有限速率（mu, theta, nu 有限）是对时间的积分，每个时间步只做一次（finite = true，
// tau = dt）。若每级都做，Heun 的一步对 y' = -lambda y 给出 R(1/2 (y + R(R(y)))) = (1 - 2 lambda dt) y + O(dt^2)：速率放大约 2 倍。
static bool relax_cell_explicit(double* Q, double dtloc, bool finite) {
	const double rho = Q[IM1] + Q[IM2], rhoe = Q[IEN] - 0.5 * (Q[IMU] * Q[IMU] + Q[IMV] * Q[IMV]) / rho;
	double a1 = Q[IA1], m1 = Q[IM1], m2 = Q[IM2], e1 = Q[IE1], e2 = Q[IE2];
	const int mech = std::isinf(mu_rate) ? CL_P : CL_MU;
	if (mech == CL_P || finite) {
		ClIn in{}; in.atr[0] = a1; in.atr[1] = 1 - a1; in.mtr[0] = m1; in.mtr[1] = m2; in.etr[0] = e1; in.etr[1] = e2;
		in.s = 0; in.rhoeC = rhoe; in.tau = dtloc; in.pguess = std::max(pbar_of(Q), 1.0);
		in.type = mech; in.nodelta = false;   // 守恒总能（"+1"）在闭合里：能量亏损按质量分数分给两相，代替 SPB 3.4 的重初始化（A3）
		ClOut o; if (!closure(in, o)) { st_cl_fail++; return false; }
		a1 = o.a[0]; m1 = o.m[0]; m2 = o.m[1]; e1 = o.e[0]; e2 = o.e[1];
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
bool RelaxationSubsystem(Grid &u_con, double h, bool finite){
	for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++)
		if (!relax_cell_explicit(u_con.cell(i, j), h, finite)) return false;
	bound_cond_u(u_con);
	return true;
}

// 一个显式 Euler 级：重构、HLLC、Godunov 更新、松弛（瞬时；finite = true 时含有限速率）
static bool explicit_stage(Grid &u_con, double h, bool finite){
	bound_cond_u(u_con);
	PrimGrid w(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells), w_yL(nxCells, nyCells), w_yR(nxCells, nyCells);
	if (!TVDReconstruction(u_con, w, w_xL, w_xR, w_yL, w_yR)) return false;
	Grid u_new(nxCells, nyCells);
	if (!Update(u_new, u_con, w, w_xL, w_xR, w_yL, w_yR, h)) return false;
	if (!RelaxationSubsystem(u_new, h, finite)) return false;
	copy(u_con, u_new);
	return true;
}

double getsoundspeed(double *con){
	Prim w; return prim_of(con, w) ? w.cf : NAN;
}


// ==================================================================================================== 输出
// 一维剖面（测试脚本读的格式）：x y p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2
void write_profile(Grid &u_con, const std::string &fn){
	FILE* f = fopen(fn.c_str(), "w"); if (!f) return;
	const int g = ghost_cells;
	fprintf(f, "# x y p u v rho alpha2 Y2 rho1 rho2 c_frozen c_wood T1 T2 p1 p2 g1-g2\n");
	for (int j = g; j < nyCells + g; j++) for (int i = g; i < nxCells + g; i++) {
		double o[15]; cell_diag(u_con.cell(i, j), o);
		fprintf(f, "%.10e %.10e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.8e %.8e %.10e %.10e %.16e %.16e %.8e\n",
			x_0 + (i - g + 0.5) * dx, y_0 + (j - g + 0.5) * dy, o[0], o[1], o[2], o[3], o[4], o[5], o[6], o[7], o[8], o[9], o[10], o[11], o[12], o[13], o[14]);
	}
	fclose(f);
}

void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime){
	if(vtk_dir.empty()) return;
	bool isPlotTime = false;
	std::string testname = "case" + std::to_string((int)test);
	for(int k=0; k<(int)IntermediatePlotTimes.size(); k++){
		double Tdiff = std::abs(IntermediatePlotTimes[k]-t);
		if(Tdiff <= dt/2.0){isPlotTime = true; break;}
	}
	if(iteration == 0) isPlotTime = true;
	if(plotINTERVAL > 0 && iteration % plotINTERVAL == 0) isPlotTime = true;
	if(isPlotTime || isFinalTime){
		ucon2upri(u_pri, u_con);
		printAllVariablesToVTK(u_con,u_pri, testname,iteration,isFinalTime);
	}
}

/* VISIT OUTPUT FILE CREATION */

void printVTKHeader(std::ofstream& outFile, int& nstep){
	outFile << "# vtk DataFile Version 3.0" << std::endl;
	outFile << "vtk output" << std::endl;
	outFile << "ASCII" << std::endl;
	outFile << "DATASET RECTILINEAR_GRID" << std::endl;
	outFile << "FIELD FieldData 2" << std::endl;
	outFile << "CYCLE 1 1 int" << std::endl;
	outFile << nstep << std::endl;
	outFile << "TIME 1 1 double" << std::endl;
	outFile << t << std::endl;

	outFile << "DIMENSIONS "<< nxCells + 1 << " " << nyCells + 1 << " 1" << std::endl;

	outFile << "X_COORDINATES " << nxCells + 1 <<" double" << std::endl;
	for(int i=0; i<nxCells + 1 ; i++){
		outFile << x_0 + ((double)i)*dx << " ";
	}
	outFile <<  std::endl;

	outFile << "Y_COORDINATES " << nyCells + 1 <<" double" << std::endl;
	for(int j=0; j<nyCells+1 ; j++){
		outFile << y_0 + ((double)j)*dy << " ";
	}
	outFile << std::endl;

	outFile << "Z_COORDINATES 1 double" << std::endl;
	outFile << "0" << std::endl;
	outFile << "CELL_DATA " << nxCells*nyCells << std::endl;
	return;
}

/** Prints a scalar variable into a VTK file. var < NUM_VARS: 守恒量；var >= NUM_VARS: u_pri 的第 var - NUM_VARS 个 */
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName){
	outFile << variableName << " 1 " << nxCells * nyCells << " double" << std::endl;
	int printWidth = 24;
	for(int j = ghost_cells; j < ghost_cells+ nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			double varPrint = var < NUM_VARS ? u_con(i,j,var) : u_pri(i,j,var - NUM_VARS);
			outFile << std::setw(printWidth) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << varPrint << " ";
		}
		outFile << std::endl;
	}
	outFile << std::endl;
}

/**Prints as many variables as desired to a single VTK file. */
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime){
	(void)isFinalTime;
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(0) << 100*CFL;
	std::string Path = vtk_dir + "/" + testname + "/explicit_" + std::to_string(nxCells) + "x"+ std::to_string(nyCells) + "_Order" + std::to_string(SpaceOrder) + "_CFL0" + ss.str();
	filesystem::create_directories(Path);
	std::string fileString = Path + "/plt_" + testname + "_" + to_string(nstep) + ".vtk";
	std::ofstream outFile(fileString.c_str());
	printVTKHeader(outFile,nstep);
	std::vector< std::tuple<int,std::string> > variablesToPrint;
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 0, std::string("rho")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 1, std::string("velX")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 2, std::string("velY")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 3, std::string("p")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 4, std::string("alpha2")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 5, std::string("Y2")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 6, std::string("T1")));
	variablesToPrint.push_back(std::make_tuple(NUM_VARS + 7, std::string("T2")));
	outFile << "FIELD FieldData " << variablesToPrint.size() << std::endl;
	for(unsigned int i = 0;i <variablesToPrint.size(); i++){
		printVTKVariable(outFile,u_con,u_pri,std::get<0>(variablesToPrint[i]),std::get<1>(variablesToPrint[i]));
	}
	outFile << "POINT_DATA " << (nxCells+1)*(nyCells+1) << std::endl;
}

static void conserved_totals(Grid &U, double M[3]) {
	M[0] = M[1] = M[2] = 0;
	for (int j = ghost_cells; j < nyCells + ghost_cells; j++) for (int i = ghost_cells; i < nxCells + ghost_cells; i++) { M[0] += U(i, j, IM1); M[1] += U(i, j, IM2); M[2] += U(i, j, IEN); }
}


int main(int argc, char *argv[]){
	double tend_override = -1;
	for (int a = 1; a < argc; a++) {
		char key[64]; char val[256]; if (sscanf(argv[a], "%63[^=]=%255s", key, val) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); return 1; }
		std::string k(key); const double d = atof(val); const bool isinf_val = std::string(val) == "inf";
		if (k == "case") test = (Test)(int)d; else if (k == "N" || k == "nx") nxCells = (int)d; else if (k == "ny") nyCells = (int)d;
		else if (k == "cfl") CFL = d; else if (k == "order") recon = (int)d < 2 ? NoRecon : Minmod; else if (k == "tend") tend_override = d; else if (k == "out") out_file = val;
		else if (k == "relax") RELAX_LEVEL = (int)d; else if (k == "epsI") EPS_I = d;
		else if (k == "mu") mu_rate = isinf_val ? INF : d; else if (k == "theta") theta_rate = isinf_val ? INF : d; else if (k == "nu") nu_rate = isinf_val ? INF : d;
		else if (k == "verbose") VERBOSE = (int)d;
		else if (k == "vtk") vtk_dir = val;
		else { fprintf(stderr, "unknown key %s\n", key); return 1; }
	}
	SpaceOrder = recon == NoRecon ? 1 : 2;
	if (theta_rate < 0) theta_rate = RELAX_LEVEL >= 1 ? INF : 0.0;
	if (nu_rate < 0) nu_rate = (RELAX_LEVEL == 2 || RELAX_LEVEL == 3) ? INF : 0.0;

	TestSetup TS; if (!test_setup(test, TS)) { fprintf(stderr, "unknown case %d\n", (int)test); return 1; }
	x_0 = TS.x0; x_1 = TS.x1; y_0 = TS.y0; y_1 = TS.y1; time_stop = TS.tend; x_bc = TS.xbc; y_bc = TS.ybc;
	if (tend_override > 0) time_stop = tend_override;
	dx = (x_1 - x_0) / nxCells;
	dy = (y_1 - y_0) / nyCells;
	if (RELAX_LEVEL > 0) build_sat();

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells);
	u_initialization(u_con);
	bound_cond_u(u_con);
	double M0[3]; conserved_totals(u_con, M0);

	const int plotINTERVAL = 0;
	std::vector<double> IntermediatePlotTimes = {0.25 * time_stop, 0.5 * time_stop, 0.75 * time_stop, time_stop};
	IntermediatePlotIntermediatePlot(u_con, u_pri, timestep, plotINTERVAL, IntermediatePlotTimes, false);

	auto t0 = std::chrono::steady_clock::now();
	do{
		dt = std::min(CFL * caldt(u_con, true), time_stop - t);
		if (!(dt > 1e-12 * time_stop)) { fprintf(stderr, "FAILED: time step underflow (dt = %.3e) at t = %.6e (step %ld)\n", dt, t, st_steps); write_profile(u_con, out_file); return 2; }
		bool done = false;
		for (int attempt = 0; attempt < 30 && !done; attempt++) {
			Grid V(nxCells, nyCells); copy(V, u_con);
			if (SpaceOrder < 2) done = explicit_stage(V, dt, true);
			else {   // Heun：U1 = R(U + dt L(U))，U2 = R(U1 + dt L(U1))，U^{n+1} = R_f((U + U2)/2)；有限速率松弛只在平均之后做一次（A3）
				Grid W(nxCells, nyCells); copy(W, V);
				done = explicit_stage(V, dt, false) && explicit_stage(V, dt, false);
				if (done) {
					for (size_t q = 0; q < V.data.size(); q++) V.data[q] = 0.5 * (V.data[q] + W.data[q]);
					done = RelaxationSubsystem(V, dt, true);
				}
			}
			if (done) copy(u_con, V);
			if (!done) { st_reject++; dt *= 0.5; if (VERBOSE) fprintf(stderr, "reject: t %.4e new dt %.3e\n", t, dt); }
		}
		if (!done) { fprintf(stderr, "FAILED at t = %.6e (step %ld)\n", t, st_steps); write_profile(u_con, out_file); return 2; }
		t = t + dt; st_steps++; timestep++;
		if (VERBOSE > 1 && st_steps % 50 == 0) fprintf(stderr, "step %ld t %.4e dt %.3e\n", st_steps, t, dt);
		IntermediatePlotIntermediatePlot(u_con, u_pri, timestep, plotINTERVAL, IntermediatePlotTimes, t >= time_stop * (1 - 1e-12));
	}while(t < time_stop * (1 - 1e-12));

	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double M1[3]; conserved_totals(u_con, M1);
	write_profile(u_con, out_file);
	printf("case %d mode 0 order %d N %d CFL %g | steps %ld rejects %ld | closure fails %ld | dM1/M %.2e dM2/M %.2e dE/E %.2e | wall %.2f s\n",
		(int)test, SpaceOrder, nxCells, CFL, st_steps, st_reject, st_cl_fail, (M1[0] - M0[0]) / M0[0], (M1[1] - M0[1]) / M0[1], (M1[2] - M0[2]) / M0[2], wall);
	return 0;
}
