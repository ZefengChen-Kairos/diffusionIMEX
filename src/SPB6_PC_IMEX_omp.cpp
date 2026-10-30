// =====================================================================================================================
//  SPB6_PC_IMEX_omp.cpp（与 spb6_physics_omp.h 一起）：SPB6_PC_IMEX.cpp 的计算效率优化版，格式、闭合与判据不变（docs/code/openmp_optimisation.md）。
//    - 求根：p 闭合与 MG 闭合的 Newton 从初猜起步并先判 Newton 步、与 p 无关的量每次闭合只算一次；激波份额的 Illinois 从弱激波份额起步
//    - 内存：级数据、Newton/BiCGStab 工作向量、RK 级张量在各级之间复用；mallopt 使大块不走 mmap
//    - OpenMP：显式重构与更新、caldt、级数据准备、Newton 的标量循环、RK 级组合；激波预算动态分配。结果与线程数无关、逐位可重复
//    - 可选 jreuse=1：同一级内复用差分 Jacobian（修正 Newton）
//  编译: g++ -O2 -std=c++17 -fopenmp -I/usr/include/eigen3 SPB6_PC_IMEX_omp.cpp -o spb6_imex_omp
// ---------------------------------------------------------------------------------------------------------------------
//  原文件头：
//  SPB6_TV_IMEX.cpp
//  SPB 2009 的 6+1 方程两相模型（单速度、双压力 + 冗余总能）的 TV 分裂 IMEX，推导见 docs/derivation/B_TV_IMEX.md。
//
//  本文件由 reference/code/2DEuler_TV.cpp 直接修改而来：全局设置、Tensor/Grid/Array、f_convective/g_convective、
//  RusanovC_x/y、Update、TVDReconstruction、u_initialization、bound_cond_u、caldt、bound_cond_array、M_P、
//  Pressure_Subsystem_AX、GMRES_Solver、PressureSubsystem、VTK 输出和 main 的 Qe/Qs/Hs 时间循环，名字、位置与作用都保留。
//  改动（理由见推导各节）：
//    - 守恒量从 (rho, rho u, rho v, E) 换成 8 个（spb6_physics.h）；状态方程、松弛闭合、算例在 spb6_physics.h，与显式程序共用。
//    - 显式（对流）部分：每个被输运量 X 的通量都是 X u（体积输运形式，B1、B2）；重构在守恒量上做，各量共用一个限制因子（B3，WB 所需）；
//      Rusanov 的耗散速度取物质速度，且不小于面态速度（B3 的保正条件）。
//    - 隐式（压力）部分：未知量仍是每格一个压力 P，但方程 G(P) = pbar_cl(P) - P = 0 是非线性的（闭合含状态方程与三种松弛，B5、B6），
//      所以用 Newton（差分 Jacobian）代替 Picard；面速度只用面两侧的压力（紧致面梯度，B4.1），每行 5 个元素；
//      线性方程组用 GMRES_Solver（与 2DEuler 一样实际是 BiCGStab，这里加了 Jacobi 预条件与相对容差）或一维 Thomas。
//      2DEuler 的 ShockDetector 与压力 Laplacian 去掉，作用由 B4.1 的声学迎风项承担。
//    - 时间：2DEuler 的半隐式 IMEX-RK（Ae/Ai/Hs）不变；起步改为声学步长后每步至多放大 1.2 倍（代替 full_step 个声学步，B10）；
//      不可容许时步长减半重算（B8、B10）；瞬时热/化学松弛每步开始时激活并投影到平衡流形（B6.7）。
//    - 能量亏损 delta 的分配（B5.2，V16）：物理耗散按 SPB (II.5) 的相 Hugoniot 能隙份额（每级先算紧致模板（±1，V21）上的激波预算），离散误差按
//      可容许能量份额；对 delta 连续，Newton 无分支校正。逐步说明见 docs/derivation/E_code_path.md。
//
//  编译: g++ -O2 -std=c++17 -fopenmp -I/usr/include/eigen3 SPB6_TV_IMEX.cpp -o spb6_imex
//  运行: ./spb6_imex key=value ...   键：case N ny cfl rus order limiter tend out relax epsI mu theta nu rk cflac ramp tol dt init nulaw nucr
//        snap probe prec p0 wellprep vtk verbose share acdt relaxform unify unifyproj（含义见 E 页 E0）
// =====================================================================================================================
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>
#include <functional>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <memory>
#include <sstream>
#include "spb6_physics_omp.h"
#ifdef __GLIBC__
#include <malloc.h>
#endif


enum Reconstruction{NoRecon, Minmod};
enum dtCalculation{Fullydt, IMEXdt};

dtCalculation caltimestep = IMEXdt;   // IMEXdt: 物质 CFL（起步除外）；Fullydt: 一直用声学 CFL
Test test = testSPB_WaterAir;         // 切换算例：改这里（或命令行 case=N）。每个算例的区域、状态方程、各边边界条件、终止时间、论文网格、松弛档位
                                      // 都在 spb6_physics.h 的 test_setup / test_state 里，物质 CFL 与 Rusanov 耗散速度（局部/全场）也按算例给（TestSetup::cfl、rus）；命令行 N= ny= relax= tend= cfl= rus= 可覆盖。
                                      // 一维：0 SPB 4.5 空化管 | 1, 2 P-S ±2 / ±500 m/s 水 | 3 十二烷 SG | 4 SPB 4.1 界面平移 | 5, 6 SPB 4.3.1 / 4.3.2 水-气
                                      //       12, 13 SPB 5.1 / 5.2 环氧-尖晶石 | 20, 21 光滑解 | 30 Pelanti 8.2.1 | 31 十二烷 NASG | 40, 41 0D 松弛
                                      //       60 SPB 4.2 硝基甲烷（Cochran-Chan）| 61 Pelanti 8.1.1 JWL-水
                                      // 二维：56 两相等熵涡（Battisti-Boscheri 4.1，p0=）| 70 SPB 4.6 激波-氦泡 | 71 SPB 4.7 空化 RM
                                      //       72 P-S 5.4 十二烷汽泡压缩（relax 2）| 73 Pelanti 8.1.2 近壁水下爆炸
x_boundary_condition x_bc;
y_boundary_condition y_bc;
int nxCells = 1000;
int nyCells = 1;
const int ghost_cells = 2;            // 2DEuler 为 5（WENO5）；这里只用二阶 minmod

int RKs = 2;                          // 1: 后向 Euler；2: LSDIRK2(2,2,2)；3: SSP2(3,3,2)

double CFL = 0.2;                     // 物质 CFL
double CFL_ac = 0.4;                  // 起步的声学 CFL
int SpaceOrder = 2;
Reconstruction recon = Minmod;
int RAMP = 1;                         // 起步：首步声学步长，之后每步至多放大 ramp_factor 倍；被限制的步用一阶（2DEuler 的 full_step）
const double ramp_factor = 1.2;
double DTFIX = 0.0;                   // > 0：固定时间步（时间收敛阶试验用）

double NewtonTol = 1e-13;             // 相对 |P| + rho c_W^2 + P_E（B8）
double P_REF = 0.0;                   // 压力未知量的参考值（B8）：Newton 的未知量是 P - P_REF。面压力差、面压力、动量通量都只含扰动量，
                                      // 低马赫时 p0 ~ 1/M^2 的舍入不进入它们（否则经隐式算子放大为约 nu^2 eps p0 的压力噪声，devlog V11）
int NewtonIter = 60;
int JREUSE = 0;                       // 1：同一级内 Newton 迭代复用上一次的差分 Jacobian（收缩率 <= 0.25 且上一步未阻尼时；见 Newton_Solver）
double GMREStol = 1e-13;              // 相对 |b|
int GMRESiter = 4000;

double t, time_stop, dt, adt;
int timestep = 0;

double x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1;
double dx, dy;
double vx_max, vy_max;                // Rusanov 耗散速度 a_Omega（x、y 方向，全场最大）
int RUS_LOCAL = 1;                    // 1：每个面取 max(|u^L|, |u^R|)；0：全场最大 vx_max（按算例，TestSetup::rus；B3）

std::string out_file = "out.dat", vtk_dir = "";
double st_Tmax[2] = {0.0, 0.0};
std::vector<double> DIVF;   // 最后一级的面散度 delta . u_f = s / tau（每格，诊断：算例 57、58 的低马赫一致性）   // 各相温度的最大值（诊断，日志行；alpha_k > 1e-9 的格）
long st_steps = 0, st_reject = 0, st_newton_it = 0, st_stage = 0, st_rk1 = 0, st_floor = 0;   // st_floor：Newton 停在闭合舍入底而被接受的级数（devlog 16.9）
double st_maxres = 0;


struct Tensor
{
	std::vector<double> data;
	int nxCells;
	int nyCells;
	int nRKS;

	Tensor(int RKS,int xcell, int ycell) : data((size_t)RKS * (xcell + 2 * ghost_cells) * (ycell + 2 * ghost_cells) * NUM_VARS, 0.0){
		nxCells = xcell;
		nyCells = ycell;
		nRKS = RKS;
	}

	double operator()(int r,int i, int j, int k)const
	{
		return data[(size_t)r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int r,int i, int j, int k){
		return data[(size_t)r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * ((size_t)i + (size_t)j * (nxCells + 2 * ghost_cells)) + k];
	}
};

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

// 显式部分与隐式部分之间传递的量（main 里按网格分配）
Array *phi_x, *phi_y;     // B3 的共用限制因子（B4.1 的二阶 psi 也用它）
Array *a2_star;           // alpha2*（显式更新后，不等于 1 - alpha1*：两者之和 = 1 - r delta[ubar^e]）
Array *ube_x, *ube_y;     // 显式中心通量用的面速度 ubar^e（x 面 (i,j) 表示 i+1/2）
Grid *wxL_g, *wxR_g, *wyL_g, *wyR_g;   // 显式级态 Q^e 的重构面值（格 i 的左/右面）：隐式增量的供体（B4.2，V21）

// 压力未知量（内部格），也是下一级 Newton 的初猜
struct PField { int nx, ny; std::vector<double> P; inline double& at(int i, int j) { return P[(size_t)i + (size_t)j * nx]; } inline double at(int i, int j) const { return P[(size_t)i + (size_t)j * nx]; } };
// 每步开始时确定的各格闭合类型（B6.7）；空 = 各级自行判断
std::vector<int> STEP_TYPE;
std::vector<double> STEP_NU;   // nulaw=1：每格的传质速率，按步开始的状态（pT 平衡态的 p、T 与 alpha_v）由 nu_law 给出，本步各级冻结


using namespace std;
void copystate(double *state_1, double *state_0);
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void g_convective(double *g, double *con);
void ucon2upri(Grid &u_pri,Grid &u_con);
void u_initialization(Grid &u);
void bound_cond_u(Grid &u);
double caldt(Grid &u_con, bool acoustic);
void bound_cond_array(Array &u, int sign_x, int sign_y);
bool Update(Grid &u_con, Grid &Qe, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR);
void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
bool PressureSubsystem(Grid &u_con, Grid &Qs, Grid &Qe, PField &P, double h, const int rks);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri, const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime);

inline int IG(int i) { return i + ghost_cells; }   // 内部格 i (0-based) -> 带虚格下标


void copystate(double *state_1, double *state_0){
	for(int i = 0; i < NUM_VARS; i++){
		state_1[i] = state_0[i];
	}
}

void copy(Grid &state_1, Grid &state_0){
	state_1.data = state_0.data;
}

// TV 分裂的对流（显式）部分：每个被输运量 X（含 alpha_1）的通量是 X u（B1）。压力部分全在 PressureSubsystem 里，所以没有 f_pressure。
void f_convective(double *f, double *con){
	double vx = con[IMU] / (con[IM1] + con[IM2]);
	for(int k = 0; k < NUM_VARS; k++){
		f[k] = con[k] * vx;
	}
}

void g_convective(double *g, double *con){
	double vy = con[IMV] / (con[IM1] + con[IM2]);
	for(int k = 0; k < NUM_VARS; k++){
		g[k] = con[k] * vy;
	}
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

// B3.1：中心部分与耗散都用重构面值（MUSCL-Rusanov）。耗散速度是物质速度：RUS_LOCAL = 1 时取本面两侧面态速度的最大值
// （B3 的保正条件是逐面的 a_f >= |u^L_f|, |u^R_f|），= 0 时取全场最大 vx_max（devlog V14）
void RusanovC_x(double *flux, double *QL, double *QR){
	double fL_c[NUM_VARS], fR_c[NUM_VARS];

	f_convective(fL_c, QL);
	f_convective(fR_c, QR);
	const double a = RUS_LOCAL ? std::max(std::fabs(QL[IMU] / (QL[IM1] + QL[IM2])), std::fabs(QR[IMU] / (QR[IM1] + QR[IM2]))) : vx_max;

	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL_c[k] + fR_c[k]) - 0.5 * a * (QR[k] - QL[k]);
	}
}

void RusanovC_y(double *flux, double *QL, double *QR){
	double gL_c[NUM_VARS], gR_c[NUM_VARS];

	g_convective(gL_c, QL);
	g_convective(gR_c, QR);
	const double a = RUS_LOCAL ? std::max(std::fabs(QL[IMV] / (QL[IM1] + QL[IM2])), std::fabs(QR[IMV] / (QR[IM1] + QR[IM2]))) : vy_max;

	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (gL_c[k] + gR_c[k]) - 0.5 * a * (QR[k] - QL[k]);
	}
}


// 显式更新 X* = Qs - adt div F^e(Qe)（B3）。另给出 alpha2*（用 alpha2 的通量 ubar^e - F[alpha1]）与面速度 ubar^e。
// 耗散速度 a_Omega 先对格速度、再对所有面态速度取最大：MUSCL 的面速度不是格速度的凸组合，a_Omega 必须不小于它（B3 的保正条件）。
bool Update(Grid &u_con, Grid &Qe, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR){
	const int g = ghost_cells, nx = nxCells, ny = nyCells;
	// 各格、各面独立：OpenMP 并行；max 归约与次序无关，结果逐位不变
	double vxm = 0, vym = 0; bool rok = true;
	#pragma omp parallel for collapse(2) schedule(static) reduction(max:vxm,vym) reduction(&&:rok)
	for(int j = 0; j < ny + 2 * g; j++){
		for(int i = 0; i < nx + 2 * g; i++){
			const double rho = Qe(i,j,IM1) + Qe(i,j,IM2);
			if(!(rho > 0)) { rok = false; continue; }
			if(i >= g && i < nx + g && j >= g && j < ny + g){
				vxm = std::max(vxm, std::fabs(Qe(i,j,IMU) / rho));
				vym = std::max(vym, std::fabs(Qe(i,j,IMV) / rho));
			}
		}
	}
	if(!rok) return false;
	#pragma omp parallel for collapse(2) schedule(static) reduction(max:vxm)
	for(int j = g; j < ny + g; j++){
		for(int i = g - 1; i < nx + g; i++){
			const double uL = w_xR(i,j,IMU) / (w_xR(i,j,IM1) + w_xR(i,j,IM2)), uR = w_xL(i + 1,j,IMU) / (w_xL(i + 1,j,IM1) + w_xL(i + 1,j,IM2));
			vxm = std::max(vxm, std::max(std::fabs(uL), std::fabs(uR)));
			(*ube_x)(i,j) = 0.5 * (uL + uR);   // = sum_k F[alpha_k]（alpha2 的面值 = 1 - alpha1 的面值）
		}
	}
	if(ny > 1){
		#pragma omp parallel for collapse(2) schedule(static) reduction(max:vym)
		for(int j = g - 1; j < ny + g; j++){
			for(int i = g; i < nx + g; i++){
				const double vL = w_yR(i,j,IMV) / (w_yR(i,j,IM1) + w_yR(i,j,IM2)), vR = w_yL(i,j + 1,IMV) / (w_yL(i,j + 1,IM1) + w_yL(i,j + 1,IM2));
				vym = std::max(vym, std::max(std::fabs(vL), std::fabs(vR)));
				(*ube_y)(i,j) = 0.5 * (vL + vR);
			}
		}
	}
	vx_max = vxm; vy_max = vym;

	#pragma omp parallel for collapse(2) schedule(static)
	for(int j = g; j < ny + g; j++){
		for(int i = g; i < nx + g; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
			}
			RusanovC_x(fp,Q_xL, Q_xR);
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			RusanovC_x(fm,Q_xL, Q_xR);

			for(int k = 0; k < NUM_VARS; k++){
				u_con(i,j,k) = Qs(i,j,k) - adt / dx * (fp[k] - fm[k]);
			}
			// alpha1：u d(alpha1)/dx 的 path-conservative 离散。涨落 D^-_f = F_f - alpha_i ubar_f、D^+_f = alpha_{i+1} ubar_f - F_f，
			// 两者之和 = ubar_f (alpha_{i+1} - alpha_i)（线性路径的积分）；即体积通量的 Rusanov 差减去本格 alpha_i 乘面速度的散度
			// （Saurel-Abgrall 1999 的形式）。体积应变用本格的 alpha，与 m、eps 的通量相容（界面处速度跳跃时微量相不失真）
			u_con(i,j,IA1) += adt / dx * Qe(i,j,IA1) * ((*ube_x)(i,j) - (*ube_x)(i - 1,j));

			if(ny > 1){
				double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					Q_yL[k] = w_yR(i,j,k);
					Q_yR[k] = w_yL(i,j + 1,k);
				}
				RusanovC_y(gp,Q_yL, Q_yR);
				for(int k = 0; k < NUM_VARS; k++){
					Q_yL[k] = w_yR(i,j - 1,k);
					Q_yR[k] = w_yL(i,j,k);
				}
				RusanovC_y(gm,Q_yL, Q_yR);
				for(int k = 0; k < NUM_VARS; k++){
					u_con(i,j,k) -= adt / dy * (gp[k] - gm[k]);
				}
				u_con(i,j,IA1) += adt / dy * Qe(i,j,IA1) * ((*ube_y)(i,j) - (*ube_y)(i,j - 1));
			}
			(*a2_star)(i,j) = 1.0 - u_con(i,j,IA1);
		}
	}
	bound_cond_u(u_con);
	bound_cond_array(*a2_star, 1, 1);
	return true;
}


// B3：体积输运量 (alpha_1, m_1, m_2, eps_1, eps_2) 共用一个限制因子 phi_i = min_X minmod(dX_L, dX_R)/(中心差分)，斜率 = phi_i * 中心差分。
// 共用一个系数，量之间的线性关系（均匀 p 时 eps_k = c_k(p) alpha_k + q_k m_k）在面值上原样成立，这是 WB 所需的。
// 速度 (u, v) 是单速度模型里独立的原始变量，各分量单独 minmod（2DEuler 的 TVDReconstruction 对原始变量逐个限制）；
// 面上的动量与总能由面值组装：(rho u)^f = rho^f u^f，E^f = sum eps_k^f + 1/2 rho^f |u^f|^2。均匀 u 时速度斜率为 0，面速度就是 u（WB 不变）。
// 若速度也共用 phi，任一量的极值会使所有量降为一阶：二维涡里 rho u 沿 x = 0 取极值，v 在其梯度最大处被降为一阶（devlog V11）。
// “模板内为常数（到舍入）”的判据用该量的舍入尺度 S_X（B3 的表；devlog V8、V11）：
//   alpha_1：1；m_k：自身大小；eps_k：sum |eps_k| + 1/2 rho |u|^2；u, v：|u| + c_f。
// phi_x/phi_y 存 min(体积输运量的 phi, 速度的光滑比)：B4.1 的二阶 psi 用它判断“所有量都光滑”。
// 诊断（facediag=1，不改变计算）：重构面态是否违反 NASG/SG 的可容许约束 alpha>0, m>0, alpha-mb>0, eps-mq-(alpha-mb)pi>0（论文引理 Face states）。
// FD_CUR[0] 面值总数，[1] 违反数，[2] 模板三格都可容许时的违反数；级成功后累加到 FD_TOT[级是否 Qe=Qs][...]
int FACEDIAG = 0; long FD_CUR[3] = {0, 0, 0}; long FD_TOT[2][3] = {{0, 0, 0}, {0, 0, 0}};
static bool fd_adm(double a1, double m1, double m2, double e1, double e2) {
	const double a[2] = {a1, 1.0 - a1}, m[2] = {m1, m2}, e[2] = {e1, e2};
	for (int k = 0; k < 2; k++) {
		if (is_mg(k)) continue;
		const double L1 = a[k] - m[k] * PH[k].b, L2 = e[k] - m[k] * PH[k].q - L1 * PH[k].pi;
		if (!(a[k] > 0 && m[k] > 0 && L1 > 0 && L2 > 0)) return false;
	}
	return true;
}
void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR){
	if (FACEDIAG) { FD_CUR[0] = FD_CUR[1] = FD_CUR[2] = 0; }
	const int g = ghost_cells, NXT = nxCells + 2 * g, NYT = nyCells + 2 * g;
	w_xL.data = u_con.data; w_xR.data = u_con.data; w_yL.data = u_con.data; w_yR.data = u_con.data;
	std::fill(phi_x->data.begin(), phi_x->data.end(), 0.0); std::fill(phi_y->data.begin(), phi_y->data.end(), 0.0);
	SpaceOrder = recon == NoRecon ? 1 : 2;
	if(recon == NoRecon) return;
	static std::vector<double> vsc, esc; vsc.assign((size_t)NXT * NYT, 0.0); esc.assign((size_t)NXT * NYT, 0.0);
	// 一维（ny = 1）时只有内部那一行的面值被用到（虚行是它的副本，结果不被读取）：只算这一行，工作量为原来的 1/5
	const int j0 = nyCells == 1 ? g : 0, j1 = nyCells == 1 ? g + 1 : NYT;
	#pragma omp parallel for collapse(2) schedule(static)
	for(int j = j0; j < j1; j++){
		for(int i = 0; i < NXT; i++){
			Prim w; const double rho = u_con(i,j,IM1) + u_con(i,j,IM2), mm = std::hypot(u_con(i,j,IMU), u_con(i,j,IMV)); const size_t c = (size_t)i + (size_t)j * NXT;
			vsc[c] = (rho > 0 ? mm / rho : 0.0) + (prim_of(u_con.cell(i, j), w) ? w.cf : 0.0);
			esc[c] = std::fabs(u_con(i,j,IE1)) + std::fabs(u_con(i,j,IE2)) + (rho > 0 ? 0.5 * mm * mm / rho : 0.0);
		}
	}
	const int VOL[5] = {IA1, IM1, IM2, IE1, IE2};
	for(int dir = 0; dir < (nyCells > 1 ? 2 : 1); dir++){
		Grid &wL = dir == 0 ? w_xL : w_yL, &wR = dir == 0 ? w_xR : w_yR;
		Array &phis = dir == 0 ? *phi_x : *phi_y;
		// 每格只写自己的面值与 phi：并行（facediag=1 的诊断计数是共享的，那时串行）
		#pragma omp parallel for collapse(2) schedule(static) if(!FACEDIAG)
		for(int j = j0; j < j1; j++){
			for(int i = 0; i < NXT; i++){
				int ia = i, ja = j, ib = i, jb = j;
				if(dir == 0){ ia = i - 1; ib = i + 1; if(ia < 0 || ib >= NXT) continue; }
				else{ ja = j - 1; jb = j + 1; if(ja < 0 || jb >= NYT) continue; }
				const size_t ca = (size_t)ia + (size_t)ja * NXT, cc = (size_t)i + (size_t)j * NXT, cb = (size_t)ib + (size_t)jb * NXT;
				double phi = 1.0;
				for(int q = 0; q < 5; q++){
					const int k = VOL[q];
					const double xa = u_con(ia,ja,k), xc = u_con(i,j,k), xb = u_con(ib,jb,k);
					double sc;
					if(k == IA1) sc = 3.0;
					else if(k == IM1 || k == IM2) sc = std::fabs(xa) + std::fabs(xc) + std::fabs(xb);
					else sc = esc[ca] + esc[cc] + esc[cb];
					const double a = xc - xa, b = xb - xc, tiny = 1e-12 * sc + 1e-300;
					if(!LIMITER) break;                                                   // 不限制（收敛阶试验）：phi = 1，斜率为中心差分
					if(std::fabs(a) <= tiny && std::fabs(b) <= tiny) continue;          // 该量在模板内（到舍入）为常数：不约束
					if(a * b <= 0){ phi = 0.0; break; }                                   // 极值点：minmod 斜率为 0
					const double c = 0.5 * (a + b), m = std::fabs(a) < std::fabs(b) ? a : b;
					phi = std::min(phi, m / c);
				}
				// 面态可容许（引理 lem:faces，V22）：NASG/SG 相的线性泛函 L1 = alpha - m b、L2 = eps - m q - L1 pi 的面值为
				// L_i -/+ (phi/4)(L_{i+1} - L_{i-1})；phi <= 2 L_i / |L_{i+1} - L_{i-1}| 使面值不低于 L_i / 2（与 minmod 对单个量的界相同）
				if(LIMITER && phi > 0.0){
					const int cs[3][2] = {{ia, ja}, {i, j}, {ib, jb}};
					for(int ph = 0; ph < 2; ph++){
						if(is_mg(ph)) continue;
						double L1[3], L2[3];
						for(int n = 0; n < 3; n++){
							const int ci = cs[n][0], cj = cs[n][1];
							const double a = ph == 0 ? u_con(ci,cj,IA1) : 1.0 - u_con(ci,cj,IA1), m = u_con(ci,cj,ph == 0 ? IM1 : IM2), e = u_con(ci,cj,ph == 0 ? IE1 : IE2);
							L1[n] = a - m * PH[ph].b; L2[n] = e - m * PH[ph].q - L1[n] * PH[ph].pi;
						}
						const double* Ls[2] = {L1, L2};
						for(int l = 0; l < 2; l++){
							const double dL = std::fabs(Ls[l][2] - Ls[l][0]);
							if(Ls[l][1] > 0.0 && dL > 0.0) phi = std::min(phi, 2.0 * Ls[l][1] / dL);
						}
					}
				}
				for(int q = 0; q < 5; q++){
					const int k = VOL[q];
					const double sig = phi * 0.5 * (u_con(ib,jb,k) - u_con(ia,ja,k));
					wL(i,j,k) = u_con(i,j,k) - 0.5 * sig;
					wR(i,j,k) = u_con(i,j,k) + 0.5 * sig;
				}
				if (FACEDIAG && i >= g && i < nxCells + g && j >= g && j < nyCells + g) {
					const bool stencil_ok = fd_adm(u_con(ia,ja,IA1), u_con(ia,ja,IM1), u_con(ia,ja,IM2), u_con(ia,ja,IE1), u_con(ia,ja,IE2))
						&& fd_adm(u_con(i,j,IA1), u_con(i,j,IM1), u_con(i,j,IM2), u_con(i,j,IE1), u_con(i,j,IE2))
						&& fd_adm(u_con(ib,jb,IA1), u_con(ib,jb,IM1), u_con(ib,jb,IM2), u_con(ib,jb,IE1), u_con(ib,jb,IE2));
					for (int side = 0; side < 2; side++) {
						Grid &w = side == 0 ? wL : wR;
						FD_CUR[0]++;
						if (!fd_adm(w(i,j,IA1), w(i,j,IM1), w(i,j,IM2), w(i,j,IE1), w(i,j,IE2))) { FD_CUR[1]++; if (stencil_ok) FD_CUR[2]++; }
					}
				}
				// 速度：各分量单独 minmod；sm 为速度的光滑比（极值处 0，模板内为常数时不约束）
				const double ra = u_con(ia,ja,IM1) + u_con(ia,ja,IM2), rc = u_con(i,j,IM1) + u_con(i,j,IM2), rb = u_con(ib,jb,IM1) + u_con(ib,jb,IM2);
				const double tiny = 1e-12 * (vsc[ca] + vsc[cc] + vsc[cb]) + 1e-300;
				double vL[2], vR[2], sm = 1.0;
				for(int d = 0; d < 2; d++){
					const int k = d == 0 ? IMU : IMV;
					const double va = u_con(ia,ja,k) / ra, vc = u_con(i,j,k) / rc, vb = u_con(ib,jb,k) / rb, a = vc - va, b = vb - vc;
					double sl = 0.0;
					if(!LIMITER) sl = 0.5 * (a + b);
					else if(!(std::fabs(a) <= tiny && std::fabs(b) <= tiny)){
						if(a * b > 0){ sl = std::fabs(a) < std::fabs(b) ? a : b; sm = std::min(sm, sl / (0.5 * (a + b))); }
						else sm = 0.0;
					}
					vL[d] = vc - 0.5 * sl; vR[d] = vc + 0.5 * sl;
				}
				const double rL = wL(i,j,IM1) + wL(i,j,IM2), rR = wR(i,j,IM1) + wR(i,j,IM2);
				wL(i,j,IMU) = rL * vL[0]; wL(i,j,IMV) = rL * vL[1]; wL(i,j,IEN) = wL(i,j,IE1) + wL(i,j,IE2) + 0.5 * rL * (vL[0] * vL[0] + vL[1] * vL[1]);
				wR(i,j,IMU) = rR * vR[0]; wR(i,j,IMV) = rR * vR[1]; wR(i,j,IEN) = wR(i,j,IE1) + wR(i,j,IE2) + 0.5 * rR * (vR[0] * vR[0] + vR[1] * vR[1]);
				phis(i,j) = std::min(phi, sm);
			}
		}
	}
}


// 初值：算例的定义（区域、状态方程、初值）在 spb6_physics.h（显式与 IMEX 共用）
void u_initialization(Grid &u){
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			double Q[NUM_VARS];
			if (SUBCELL <= 1) test_state(test, x, y, Q);
			else {   // 格平均：SUBCELL x SUBCELL 个子点上守恒量的平均（圆形区域，二维算例）
				for (int k = 0; k < NUM_VARS; k++) Q[k] = 0.0;
				for (int a = 0; a < SUBCELL; a++) for (int b = 0; b < SUBCELL; b++) {
					double Qs[NUM_VARS]; test_state(test, x + ((a + 0.5) / SUBCELL - 0.5) * dx, y + ((b + 0.5) / SUBCELL - 0.5) * dy, Qs);
					for (int k = 0; k < NUM_VARS; k++) Q[k] += Qs[k] / (SUBCELL * SUBCELL);
				}
			}
			for(int k = 0; k < NUM_VARS; k++) u(i,j,k) = Q[k];
		}
	}
}


// 透射：复制边界格；固壁：镜像，只把法向动量取负（2DEuler 原来两个动量分量都取负）；周期：绕回。一维时 y 向只复制。
void bound_cond_u(Grid &u){
	bc_apply(u, u.nxCells, u.nyCells, ghost_cells);   // spb6_physics.h：每一边单独的透射 / 固壁（含活塞速度）/ 周期
}


// 时间步：acoustic = true 时用 |u| + c_f（起步、Fullydt），否则用物质速度 |u|（IMEXdt）。返回 1 / (a_x/dx + a_y/dy)；dt = CFL * 返回值
double caldt(Grid &u_con, bool acoustic){
	double ux = 0, uy = 0; const int g = ghost_cells;
	#pragma omp parallel for collapse(2) schedule(static) reduction(max:ux,uy)
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

void bound_cond_array(Array &u, int sign_x = 1, int sign_y = 1){
	const int g = ghost_cells, nx = u.nxCells, ny = u.nyCells;
	auto src = [&](int side, int q, int n, double sg, int& sL, double& fs) {   // 每一边单独（BCS）
		fs = 1.0;
		if (BCS[side] == BC_PERIODIC) sL = (side % 2 == 0) ? n + g - 1 - q : g + q;
		else if (BCS[side] == BC_WALL) { sL = (side % 2 == 0) ? g + q : n + g - 1 - q; fs = sg; }
		else sL = (side % 2 == 0) ? g : n + g - 1; };
	for(int j = 0; j < ny + 2 * g; j++){
		for(int q = 0; q < g; q++){
			int sL, sR; double fL, fR; src(0, q, nx, sign_x, sL, fL); src(1, q, nx, sign_x, sR, fR);
			u(g - 1 - q, j) = fL * u(sL, j); u(nx + g + q, j) = fR * u(sR, j);
		}
	}
	for(int i = 0; i < nx + 2 * g; i++){
		for(int q = 0; q < g; q++){
			int sL, sR; double fL, fR;
			if(ny == 1 && BCS[2] != BC_PERIODIC){ sL = g; sR = g; fL = fR = 1.0; }
			else { src(2, q, ny, sign_y, sL, fL); src(3, q, ny, sign_y, sR, fR); }
			u(i, g - 1 - q) = fL * u(i, sL); u(i, ny + g + q) = fR * u(i, sR);
		}
	}
}


// ==================================================================================================== B6.7：激活与投影
// 瞬时松弛（速率 = inf）在激活时是状态的跳跃；若在级内激活，IMEX-RK 的级组合（LSDIRK2 的 Q^s_2 = U^n + ((1-g)/g)(U_1 - U^n)，系数 2.41）
// 会外推这个跳跃，微量相的含量变负。刚性精确 IMEX-RK 的 AP 性质要求初值在流形上。所以每步开始时按 U^n 判断激活（P-S 在双曲步之后判断），
// 把 U^n 投影到该闭合的流形上（s = 0 的同一个闭合，守恒 m 或 rho、rho u、E），本步各级冻结这些类型。有限速率的类型不投影。
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
	ClOut o; if (!closure(in, o)) return false;   // 也是可行性判据：无解的格退回力学类型
	if (!UNIFY_PROJ) return true;                 // unifyproj=0：类型照常判断与冻结，U^n 不投影
	U(i, j, IA1) = o.a[0]; U(i, j, IM1) = o.m[0]; U(i, j, IM2) = o.m[1]; U(i, j, IE1) = o.e[0]; U(i, j, IE2) = o.e[1];
	return true;
}
// 对整个场：判断类型、投影；投影失败的格退回力学类型（不投影）
void RelaxationActivation(Grid& U) {
	const int nx = U.nxCells, ny = U.nyCells;
	STEP_TYPE.assign((size_t)nx * ny, CL_P); if (NU_LAW) STEP_NU.assign((size_t)nx * ny, 0.0);
	#pragma omp parallel for schedule(static)
	for (int q = 0; q < nx * ny; q++) {
		const int i = q % nx, j = q / nx;
		int ty = cell_type_U(U, IG(i), IG(j));
		if (!project_cell(U, IG(i), IG(j), ty)) ty = std::isinf(mu_rate) ? CL_P : CL_MU;
		STEP_TYPE[(size_t)i + (size_t)j * nx] = ty;
		if (NU_LAW) {   // 速率取本格 pT 平衡态（mu = theta = inf，不传质）的 p、T 与汽相体积分数
			double nu = 0.0;
			if (ty == CL_NU) {
				const double rho = U(IG(i), IG(j), IM1) + U(IG(i), IG(j), IM2);
				ClIn in{}; in.atr[0] = U(IG(i), IG(j), IA1); in.atr[1] = 1.0 - in.atr[0]; in.mtr[0] = U(IG(i), IG(j), IM1); in.mtr[1] = U(IG(i), IG(j), IM2);
				in.etr[0] = U(IG(i), IG(j), IE1); in.etr[1] = U(IG(i), IG(j), IE2); in.s = 0.0; in.tau = 0.0; in.type = CL_PT; in.Z[0] = in.Z[1] = 1.0;
				in.rhoeC = U(IG(i), IG(j), IEN) - 0.5 * (U(IG(i), IG(j), IMU) * U(IG(i), IG(j), IMU) + U(IG(i), IG(j), IMV) * U(IG(i), IG(j), IMV)) / rho;
				in.pguess = std::max(p_mix(in.atr[0], in.mtr[0], in.mtr[1], in.rhoeC), 1.0);
				ClOut o; if (closure(in, o)) nu = nu_law(o.p[0], o.T[0], o.a[1]);
			}
			STEP_NU[(size_t)i + (size_t)j * nx] = nu;
		}
	}
	bound_cond_u(U);
}


// ==================================================================================================== 隐式（压力）部分
// 一级的面数据（B4）：显式预测态 X*、单位体积内容、供体速度、面速度/面压力的系数、闭合类型、残量尺度
struct StageData {
	int nx, ny, NXT, NYT; double rx, ry, tau;
	Grid Xs;                       // 显式更新后的 X*
	Array a2s;                     // alpha2*
	std::vector<double> cont;      // 每格单位体积内容 (alpha1, alpha2, m1, m2, e1, e2) / sum alpha*（闭合类型判断用）
	std::vector<double> dxL, dxR, dyL, dyR;   // 增量供体：面 (i,j) 两侧的显式重构面值 (alpha1, 1 - alpha1, m1, m2, e1, e2)（B4.2，V21）
	Array ust, vst;                // u*, v*（供体速度）
	Array ubsx, kpx, ubex, ubsy, kpy, ubey;   // 面量（x 面 (i,j) = i+1/2）：u_f = ubs - kp (P_R - P_L)，ube = 显式中心通量的面速度
	Array wLx, wLy, pvx, pvy;      // 面压力：p_f = wL P_L + (1 - wL) P_R - pv（B4.1 的半格动量平衡 + 声学 Riemann 解）
	std::vector<int> type; std::vector<double> Z0, Z1;
	std::vector<double> Qph;       // 局部激波的物理耗散预算 Q_k（每格两相，B5.2 的 shock_budget，按级冻结）
	std::vector<double> Kw;        // 残量尺度 K + P_E（X* 处，B8）
	StageData(int x, int y) : nx(x), ny(y), NXT(x + 2 * ghost_cells), NYT(y + 2 * ghost_cells), Xs(x, y), a2s(x, y), ust(x, y), vst(x, y),
		ubsx(x, y), kpx(x, y), ubex(x, y), ubsy(x, y), kpy(x, y), ubey(x, y), wLx(x, y, 0.5), wLy(x, y, 0.5), pvx(x, y), pvy(x, y) {
		cont.assign((size_t)NXT * NYT * 6, 0.0); dxL.assign((size_t)NXT * NYT * 6, 0.0); dxR = dxL; dyL = dxL; dyR = dxL; Qph.assign((size_t)NXT * NYT * 2, 0.0); Kw.assign((size_t)NXT * NYT, 0.0); type.assign((size_t)NXT * NYT, CL_P); Z0.assign((size_t)NXT * NYT, 1.0); Z1.assign((size_t)NXT * NYT, 1.0);
	}
	inline size_t id(int i, int j) const { return (size_t)i + (size_t)j * NXT; }
};

// 面上的隐式增量（B4.2'）：总体积通量对 u_f 迎风。显式部分已按 ubar^e 从其上游格扫过 ubar^e tau，增量补足到 u_f tau：
//   I_f[X] = X^L_f wL + X^R_f wR,  wL = (u_f)^+ - (ubar^e)^+,  wR = (u_f)^- - (ubar^e)^-
// 供体 X^{L,R}_f 是显式 Rusanov 用的同一组重构面值（V21）：显式通量 + 增量 = X^L [(u^L + a)/2 + wL] + X^R [(u^R - a)/2 + wR]，
// a >= |u^{L,R}| 时邻格系数自动非负，退回的量不超过显式送来的量；本格系数非负只需 r (a + |u_f|) <= 1/2。于是 Q^e = Q^s 的级里
// 每一相 (alpha_k, m_k, eps_k) 是可容许面态的非负组合，与含量比无关（微量相保正）。按格内含量退回时余量与 Xhat_nb / Xhat_i 成比例（devlog 20.5）。
// 同号时就是 Xhat_up du（du < 0 退回的是显式扫过来的同一种物质）；反向时退回显式扫过的全部、再从另一侧扫 u_f。
// wL + wR = du，故 sum_k alphahat 的增量 = du（B4.4 不变）；I_f 对 P 连续（分段线性）。
struct FaceVals { double uf, du, pf, wL, wR, utL, utR; const double* cL; const double* cR; };
inline void face_weights(FaceVals& f, double ube) { f.wL = std::max(f.uf, 0.0) - std::max(ube, 0.0); f.wR = std::min(f.uf, 0.0) - std::min(ube, 0.0); }
// 固壁面：返回边号（0 左 1 右 2 下 3 上），否则 -1
inline int wall_x(const StageData& S, int i) { if (i == ghost_cells - 1 && BCS[0] == BC_WALL) return 0; if (i == S.nx + ghost_cells - 1 && BCS[1] == BC_WALL) return 1; return -1; }
inline int wall_y(const StageData& S, int j) { if (j == ghost_cells - 1 && BCS[2] == BC_WALL) return 2; if (j == S.ny + ghost_cells - 1 && BCS[3] == BC_WALL) return 3; return -1; }
// x 面 (i+1/2, j)，左格压力 PL，右格压力 PR（B4.1, B4.2）
inline void face_x(const StageData& S, int i, int j, double PL, double PR, FaceVals& f) {
	f.pf = S.wLx(i, j) * PL + (1.0 - S.wLx(i, j)) * PR - S.pvx(i, j);
	f.cL = &S.dxL[S.id(i, j) * 6]; f.cR = &S.dxR[S.id(i, j) * 6]; f.utL = S.vst(i, j); f.utR = S.vst(i + 1, j);
	if (const int w = wall_x(S, i); w >= 0) {   // 固壁：面速度 = 壁速度（静止壁 0；活塞时体积通量从虚格进入）
		if (BCU[w] == 0.0) { f.uf = 0.0; f.du = 0.0; f.wL = f.wR = 0.0; return; }
		f.uf = BCU[w]; f.du = f.uf - S.ubex(i, j); face_weights(f, S.ubex(i, j)); return; }
	f.uf = S.ubsx(i, j) - S.kpx(i, j) * (PR - PL);
	f.du = f.uf - S.ubex(i, j);
	face_weights(f, S.ubex(i, j));
}
inline void face_y(const StageData& S, int i, int j, double PL, double PR, FaceVals& f) {
	f.pf = S.wLy(i, j) * PL + (1.0 - S.wLy(i, j)) * PR - S.pvy(i, j);
	f.cL = &S.dyL[S.id(i, j) * 6]; f.cR = &S.dyR[S.id(i, j) * 6]; f.utL = S.ust(i, j); f.utR = S.ust(i, j + 1);
	if (const int w = wall_y(S, j); w >= 0) {
		if (BCU[w] == 0.0) { f.uf = 0.0; f.du = 0.0; f.wL = f.wR = 0.0; return; }
		f.uf = BCU[w]; f.du = f.uf - S.ubey(i, j); face_weights(f, S.ubey(i, j)); return; }
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
	F.E = inc(f, 4) + inc(f, 5) + 0.5 * f.uf * f.uf * M + 0.5 * (f.utL * f.utL * ML + f.utR * f.utR * MR) + (P_REF + f.pf) * f.uf;
	return F;
}

struct CellOut { ClOut cl; double mx, my, E, s, sig; };
// 闭合压力的舍入底（devlog 16.9）：能量亏损 D0 = (rho e)^C - sum eps^tr 是大数之差，带 eps_m (|(rho e)^C| + sum |eps^tr| + sum |m q|) 的舍入；
// 它按 lambda_k 进入相 k 的 A_k，经 d alpha_k / d A_k = (gamma_k - 1)/(gamma_k (p + pi_k)) 与 Phi'(p) = sum d alpha_k/dp 变成压力误差。
// 亏损分给软相（低压下的微量气体）时这个底远高于混合物刚度给出的底。只对 NASG/SG 的 p 闭合计算，其余返回 0。
static double closure_roundoff(const ClIn& in, double p) {
	if (in.type != CL_P || in.nodelta || is_mg(0) || is_mg(1)) return 0.0;
	PCoef c; if (!p_coef(in, c, p)) return 0.0;
	double dphi = 0.0, sens = 0.0;
	for (int k = 0; k < 2; k++) { dphi += dalpha_p(c, k, p); sens += std::fabs(c.lam[k]) * (PH[k].g - 1.0) / (PH[k].g * (p + PH[k].pi)); }
	const double en = 2.220446049250313e-16 * (std::fabs(in.rhoeC) + std::fabs(in.etr[0]) + std::fabs(in.etr[1]) + std::fabs(in.mtr[0] * PH[0].q) + std::fabs(in.mtr[1] * PH[1].q));
	return dphi < 0 ? en * sens / -dphi : 0.0;
}
// 单格残量 G = pbar_cl - Pc（B4.3-B4.6, B5.4）。邻格压力由调用者给（边界虚格按边界条件取值）
static bool cell_residual(const StageData& S, int i, int j, double Pc, double PW, double PE, double PS, double PN, double& G, CellOut* out) {
	FaceVals fe, fw, fn, fs;
	face_x(S, i, j, Pc, PE, fe); face_x(S, i - 1, j, PW, Pc, fw);
	const bool twoD = S.ny > 1;
	if (twoD) { face_y(S, i, j, Pc, PN, fn); face_y(S, i, j - 1, PS, Pc, fs); }
	double X[6] = {S.Xs(i, j, IA1), S.a2s(i, j), S.Xs(i, j, IM1), S.Xs(i, j, IM2), S.Xs(i, j, IE1), S.Xs(i, j, IE2)};
	for (int k = 2; k < 6; k++) {
		X[k] -= S.rx * (inc(fe, k) - inc(fw, k));
		if (twoD) X[k] -= S.ry * (inc(fn, k) - inc(fs, k));
	}
	// alpha 只在显式步由 path-conservative Rusanov 更新（u^e），隐式步不再输运 alpha；alpha2 = 1 - alpha1
	X[1] = 1.0 - X[0];
	// 体积应变 s = tau div(u_f)（面速度的离散散度）；压力功 alpha_k^tr p_k s。闭合的输入写成 alpha_k^tr (1 - s)，
	// 于是 omega_k = alpha_k^tr、alpha_k = alpha_k^tr + chi_k R，闭合的公式不变
	const double s = S.rx * (fe.uf - fw.uf) + (twoD ? S.ry * (fn.uf - fs.uf) : 0.0);
	X[0] *= 1.0 - s; X[1] *= 1.0 - s;
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
	in.Z[0] = S.Z0[S.id(i, j)]; in.Z[1] = S.Z1[S.id(i, j)]; in.pguess = P_REF + Pc; in.nodelta = false; in.Qph[0] = S.Qph[S.id(i, j) * 2]; in.Qph[1] = S.Qph[S.id(i, j) * 2 + 1];
	if (NU_LAW && !STEP_NU.empty()) in.nu = STEP_NU[(size_t)(i - ghost_cells) + (size_t)(j - ghost_cells) * S.nx];
	ClOut o; if (!closure(in, o)) {
		if (VERBOSE >= 3) fprintf(stderr, "    closure fail cell (%d,%d) type %d atr %.6e %.6e mtr %.6e %.6e etr %.6e %.6e s %.3e rhoeC %.6e Pc %.6e\n", i - ghost_cells, j - ghost_cells, in.type, in.atr[0], in.atr[1], in.mtr[0], in.mtr[1], in.etr[0], in.etr[1], s, in.rhoeC, Pc);
		return false; }
	G = (o.pbar - P_REF) - Pc;
	if (out) { out->cl = o; out->mx = mx; out->my = my; out->E = E; out->s = s; out->sig = closure_roundoff(in, o.pbar); }
	return std::isfinite(G);
}

// 2DEuler 的 M_P：带虚格的下标 (i, j) -> 压力未知量的编号。透射/固壁的虚格对应边界格本身，周期绕回。
// 边界虚格的压力是本格压力的一个副本，在 Jacobian 里是一个单独的元素（列号为本格），不并进对角元：
// 这样边界行与内部行的算术逐位相同，y 向均匀的问题各行逐位相同（devlog V8）。
int M_P(int i, int j){
	int I, J;
	if(i < ghost_cells) I = (x_bc == x_Periodic) ? nxCells + i : ghost_cells;
	else if(i > nxCells + ghost_cells - 1) I = (x_bc == x_Periodic) ? i - nxCells : nxCells + ghost_cells - 1;
	else I = i;
	if(j < ghost_cells) J = (y_bc == y_Periodic && nyCells > 1) ? nyCells + j : ghost_cells;
	else if(j > nyCells + ghost_cells - 1) J = (y_bc == y_Periodic && nyCells > 1) ? j - nyCells : nyCells + ghost_cells - 1;
	else J = j;
	return (I - ghost_cells) + nxCells * (J - ghost_cells);
}
inline double Pnb(const PField& P, int i, int j, int di, int dj) { return P.P[M_P(IG(i) + di, IG(j) + dj)]; }
static bool residual_at(const StageData& S, const PField& P, int i, int j, double dc, int which, double dn, double& G, CellOut* out = nullptr) {
	// dc 扰动本格；which = 0 W, 1 E, 2 S, 3 N 只扰动该邻格的压力 dn（-1 不扰动）
	const double Pc = P.at(i, j) + dc;
	double PW = Pnb(P, i, j, -1, 0), PE = Pnb(P, i, j, 1, 0), PS = Pnb(P, i, j, 0, -1), PN = Pnb(P, i, j, 0, 1);
	if (which == 0) PW += dn;
	if (which == 1) PE += dn;
	if (which == 2) PS += dn;
	if (which == 3) PN += dn;
	if (i == 0 && BCS[0] == BC_PRES) PW = BCP[0] - P_REF;              // 给定压力出口：外侧压力不是未知量（Jacobian 的该列为 0）
	if (i == S.nx - 1 && BCS[1] == BC_PRES) PE = BCP[1] - P_REF;
	return cell_residual(S, IG(i), IG(j), Pc, PW, PE, PS, PN, G, out);
}

// Jacobian 的一行：本格 + W, E, S, N（列号可重复，不合并）
struct SpRow { int n = 0; int col[5]; double val[5]; };

// 2DEuler 的 Pressure_Subsystem_AX：算子乘向量。这里的算子是 G(P) 的差分 Jacobian（逐行固定次序求和）
void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, const Eigen::VectorXd& x, const std::vector<SpRow>& A){
	const int n = (int)A.size();
	#pragma omp parallel for schedule(static)
	for(int r = 0; r < n; r++){ double s = 0; for(int q = 0; q < A[r].n; q++) s += A[r].val[q] * x[A[r].col[q]]; Ax[r] = s; }
}
// 求和按固定大小的块（块内顺序、块间顺序都固定）：结果与线程数无关、逐位可重复
// 块大小固定（与线程数无关），1024 使 n = 4000–100000 的二维问题也能分到全部线程（4096 时 n < 8192 只有 1–2 块，内积实际是串行的）
template <class F> static double block_sum(int n, F term) {
	const int B = 1024, nb = (n + B - 1) / B; std::vector<double> part(nb, 0.0);
	#pragma omp parallel for schedule(static)
	for (int c = 0; c < nb; c++) { double s = 0; const int e = std::min(n, (c + 1) * B); for (int r = c * B; r < e; r++) s += term(r); part[c] = s; }
	double s = 0; for (int c = 0; c < nb; c++) s += part[c]; return s;
}

// 2DEuler 的 GMRES_Solver（原函数实际是 BiCGStab）：加两级预条件（Jacobi 用本格扰动的那个元素，粗空间为常数向量）与相对容差
bool GMRES_Solver(Eigen::VectorXd& x, const Eigen::VectorXd& b, const std::vector<SpRow>& A){
	const int n = (int)b.size();
	// 工作向量在各次调用之间复用（每次 Newton 迭代都调用一次）
	static Eigen::VectorXd dinv, r, r0, p, v, s, t, ph, sh, Ax;
	for (Eigen::VectorXd* w : {&dinv, &r, &r0, &p, &v, &s, &t, &ph, &sh, &Ax}) w->resize(n);
	#pragma omp parallel for schedule(static)
	for (int q = 0; q < n; q++) dinv[q] = (A[q].n > 0 && A[q].val[0] != 0) ? 1.0 / A[q].val[0] : 1.0;
	auto dot = [&](const Eigen::VectorXd& a, const Eigen::VectorXd& c) { return block_sum(n, [&](int r) { return a[r] * c[r]; }); };
	// 两级预条件（B8）：Jacobi + 常数向量的粗空间校正 M^{-1} r = D^{-1} r + 1 (1^T r)/(1^T A 1)（Nicolaides）。
	// 面速度与动量只依赖压力差，P 的常数平移只使 G 变化 O(1)，而对角元约 1 + 4 nu^2：Jacobi 预条件后常数模态的特征值约 1/(4 nu^2)，
	// M = 1e-6 时约 1e-10，BiCGStab 发散（devlog V11）。粗空间校正把它移到 1 附近，其余模态不变
	const double E = block_sum(n, [&](int r) { double e = 0; for (int q = 0; q < A[r].n; q++) e += A[r].val[q]; return e; });
	const double Einv = (E != 0.0 && std::isfinite(E)) ? 1.0 / E : 0.0;
	auto prec = [&](const Eigen::VectorXd& in, Eigen::VectorXd& out) { const double cc = block_sum(n, [&](int r) { return in[r]; }) * Einv;
		#pragma omp parallel for schedule(static)
		for (int r = 0; r < n; r++) out[r] = dinv[r] * in[r] + cc; };
	x.resize(n); prec(b, x);
	p.setZero(); v.setZero();
	Pressure_Subsystem_AX(Ax, x, A);
	#pragma omp parallel for schedule(static)
	for (int q = 0; q < n; q++) { r[q] = b[q] - Ax[q]; r0[q] = r[q]; }
	const double bn = std::sqrt(dot(b, b)) + 1e-300; double rho = 1, alpha = 1, omega = 1;
	for (int it = 0; it < GMRESiter; it++) {
		if (std::sqrt(dot(r, r)) <= GMREStol * bn) return true;
		const double rho1 = dot(r0, r); if (rho1 == 0) return false;
		const double beta = (rho1 / rho) * (alpha / omega); rho = rho1;
		#pragma omp parallel for schedule(static)
		for (int q = 0; q < n; q++) p[q] = r[q] + beta * (p[q] - omega * v[q]);
		prec(p, ph);
		Pressure_Subsystem_AX(v, ph, A); alpha = rho / dot(r0, v);
		#pragma omp parallel for schedule(static)
		for (int q = 0; q < n; q++) s[q] = r[q] - alpha * v[q];
		if (std::sqrt(dot(s, s)) <= GMREStol * bn) {
			#pragma omp parallel for schedule(static)
			for (int q = 0; q < n; q++) x[q] += alpha * ph[q];
			return true; }
		prec(s, sh);
		Pressure_Subsystem_AX(t, sh, A); omega = dot(t, s) / dot(t, t);
		#pragma omp parallel for schedule(static)
		for (int q = 0; q < n; q++) { x[q] += alpha * ph[q] + omega * sh[q]; r[q] = s[q] - omega * t[q]; }
	}
	return std::sqrt(dot(r, r)) <= 1e3 * GMREStol * bn;
}

// 修正 Newton 的重建判据：上一步完整（未阻尼）且残量至少缩小到 1/4
static inline bool fresh_ok(double ratio, double lamprev) { return lamprev == 1.0 && ratio <= 0.25; }
// B8：全局 Newton，G_i(P) = pbar_cl,i - P_i。Jacobian 用逐格差分；一维非周期用 Thomas，其余用 GMRES_Solver；线搜索保证可容许且残量下降
static bool Newton_Solver(const StageData& S, PField& P, int& iters, double& resfinal) {
	const int nx = P.nx, ny = P.ny, n = nx * ny;
	// 工作数组在各级之间复用（只改尺寸，不重新分配）
	static std::vector<double> G, Gn, lo, di, up, sc, cp, dp; static Eigen::VectorXd dP, rhs; static PField Pn; static std::vector<SpRow> A;
	for (std::vector<double>* w : {&G, &Gn, &lo, &di, &up, &sc, &cp, &dp}) w->resize(n);
	dP.resize(n); rhs.resize(n); A.resize(n); Pn = P;
	auto evalAll = [&](const PField& Q, std::vector<double>& R) -> bool {
		bool ok = true;   // 各格独立（OpenMP；每格的算术与串行相同，结果逐位不变）
		#pragma omp parallel for schedule(static) reduction(&&:ok)
		for (int q = 0; q < n; q++) { double g; if (!residual_at(S, Q, q % nx, q / nx, 0, -1, 0, g)) ok = false; else R[q] = g; }
		return ok; };
	// 残量范数：|G| / (|P| + K + P_E)；按 |P| 归一在 p << pi 时达不到（B8）
	auto norm = [&](const std::vector<double>& R, const PField& Q) { double m = 0;
		#pragma omp parallel for schedule(static) reduction(max:m)
		for (int q = 0; q < n; q++) { const int i = q % nx, j = q / nx; m = std::max(m, std::fabs(R[q]) / (std::fabs(P_REF + Q.P[q]) + S.Kw[S.id(IG(i), IG(j))] + 1.0)); }
		return m; };
	if (!evalAll(P, G)) { iters = 0; resfinal = INF; return false; }
	double res = norm(G, P);
	const bool tri = (ny == 1 && x_bc != x_Periodic);
	// P 的自然尺度（B8）：s_q = N_q / D_q。N = |P| + K + P_E + 1 是残量尺度；D = 1 + (K + P_E) sum_f (tau/h_f) kp_f 是 dG/dP 对角元的尺度（约 1 + nu^2），
	// s 就是残量所能分辨的压力变化。差分步长取 1e-7 s：步长取 1e-7 |P| 时扰动在本格引起的体积变化约 1e-7 nu^2/gamma，
	// 低马赫时 nu ~ 1/M，M = 1e-5 时已是 O(1)，差分 Jacobian 失效（devlog V11）
	// jreuse=1（修正 Newton）：Jacobian 只在级的第一次迭代、上一步收缩率 > 0.25、上一步用了阻尼、或复用的 Jacobian 让线搜索失败时重建。
	// 收敛判据不变（残量 < tol，或步长 < tol：用旧 Jacobian 时收敛是线性的，误差约为收缩率乘最后一步，仍在容差内）
	bool fresh = false; double ratio = 1.0, lamprev = 1.0;
	for (int it = 0; it < NewtonIter; it++) {
		iters = it; resfinal = res;
		if (VERBOSE >= 3) { int qm = 0; for (int q = 0; q < n; q++) if (std::fabs(G[q]) > std::fabs(G[qm])) qm = q; fprintf(stderr, "    newton it %d res %.3e  worst cell %d G %.4e P %.8e\n", it, res, qm, G[qm], P_REF + P.P[qm]); }
		if (res < NewtonTol) return true;
		const bool build = !JREUSE || it == 0 || !fresh_ok(ratio, lamprev);
		if (build) {
		#pragma omp parallel for schedule(static)
		for (int q = 0; q < n; q++) {
			const int I = IG(q % nx), J = IG(q / nx); const double K = S.Kw[S.id(I, J)];
			double D = S.rx * (S.kpx(I, J) + S.kpx(I - 1, J));
			if (ny > 1) D += S.ry * (S.kpy(I, J) + S.kpy(I, J - 1));
			sc[q] = (std::fabs(P_REF + P.P[q]) + K + 1.0) / (1.0 + K * D);
		}
		bool jok = true;
		#pragma omp parallel for schedule(static) reduction(&&:jok)
		for (int r = 0; r < n; r++) {
			const int i = r % nx, j = r / nx; const double h = 1e-7 * sc[r];
			double g1; if (!residual_at(S, P, i, j, h, -1, 0, g1)) { jok = false; continue; }
			SpRow& row = A[r]; row.n = 0; row.col[row.n] = r; row.val[row.n++] = (g1 - G[r]) / h;
			// 邻格元素按 W, E, S, N 的固定次序存放；同列（边界虚格 = 本格，或周期小网格上重合的邻格）也分开存，不合并
			const int dd[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
			for (int w = 0; w < 4; w++) {
				if (ny == 1 && w >= 2) break;
				const int c = M_P(IG(i) + dd[w][0], IG(j) + dd[w][1]);
				const double hn = 1e-7 * sc[c];
				double g2; if (!residual_at(S, P, i, j, 0, w, hn, g2)) { jok = false; break; }
				row.col[row.n] = c; row.val[row.n++] = (g2 - G[r]) / hn;
			}
		}
		if (!jok) return false;
		}
		fresh = build;
		#pragma omp parallel for schedule(static)
		for (int q = 0; q < n; q++) rhs[q] = -G[q];
		if (tri) {
			#pragma omp parallel for schedule(static)
			for (int r = 0; r < n; r++) { lo[r] = up[r] = di[r] = 0; for (int q = 0; q < A[r].n; q++) { const int c = A[r].col[q]; if (c == r) di[r] += A[r].val[q]; else if (c == r - 1) lo[r] = A[r].val[q]; else if (c == r + 1) up[r] = A[r].val[q]; } }
			cp[0] = up[0] / di[0]; dp[0] = rhs[0] / di[0];
			for (int r = 1; r < n; r++) { const double m = di[r] - lo[r] * cp[r - 1]; cp[r] = up[r] / m; dp[r] = (rhs[r] - lo[r] * dp[r - 1]) / m; }
			dP[n - 1] = dp[n - 1]; for (int r = n - 2; r >= 0; r--) dP[r] = dp[r] - cp[r] * dP[r + 1];
		} else {
			if (!GMRES_Solver(dP, rhs, A)) return false;
		}
		// 步长判据（B8，devlog V11）：Newton 步相对未知量的大小 |P - P_REF| 与残量的分辨率 s 已在容差内时，P + dP 就是解。
		// 残量可能有高于容差的舍入底（P 的舍入经面速度、体积通量与闭合放大约 nu^2 倍，nu = c tau/h），这时残量判据达不到，而这个判据可以。
		// 尺度不用 |P| + K + P_E：低马赫时它比 s 大约 nu^2 倍，M = 1e-6 时第一步就被接受，残量停在 6e-12（devlog V11）
		double sn = 0.0;
		#pragma omp parallel for schedule(static) reduction(max:sn)
		for (int q = 0; q < n; q++) sn = std::max(sn, std::fabs(dP[q]) / (std::fabs(P.P[q]) + sc[q]));
		if (sn < NewtonTol) {
			#pragma omp parallel for schedule(static)
			for (int q = 0; q < n; q++) P.P[q] += dP[q];
			st_newton_it++; iters = it + 1;
			if (!evalAll(P, G)) return false;
			resfinal = norm(G, P); return true;
		}
		double lam = 1.0; bool acc = false; const double res0 = res;
		for (int ls = 0; ls < 40; ls++) {
			#pragma omp parallel for schedule(static)
			for (int q = 0; q < n; q++) Pn.P[q] = P.P[q] + lam * dP[q];
			if (evalAll(Pn, Gn)) { const double rn = norm(Gn, Pn); if (rn < res * (1.0 - 1e-4 * lam) || rn < NewtonTol) { acc = true; res = rn; break; } }
			lam *= 0.5;
		}
		if (!acc && !fresh) { ratio = 1.0; lamprev = 0.0; continue; }   // 复用的 Jacobian 不下降：重建后重算这一步（P 不变）
		ratio = res / res0; lamprev = lam;
		// 线搜索无法再降低残量：到了舍入底。残量已在 10 倍容差之内（相对 1e-12，远小于离散误差）就接受当前解，否则失败（devlog 16.9）
		if (!acc) {   // 线搜索无法再降低残量：逐格检查残量是否已在 tol N_i + sigma_i 之内（sigma_i 为闭合压力的舍入底，devlog 16.9）；是则收敛，否则失败
			resfinal = res; bool floor_ok = true;
			#pragma omp parallel for schedule(static) reduction(&&:floor_ok)
			for (int q = 0; q < n; q++) {
				const int i = q % nx, j = q / nx; double g; CellOut o;
				if (!residual_at(S, P, i, j, 0, -1, 0, g, &o)) { floor_ok = false; continue; }
				if (std::fabs(g) > NewtonTol * (std::fabs(P_REF + P.P[q]) + S.Kw[S.id(IG(i), IG(j))] + 1.0) + o.sig) floor_ok = false;
			}
			if (floor_ok) st_floor++;
			return floor_ok; }
		std::swap(P.P, Pn.P); std::swap(G, Gn);
		st_newton_it++;
	}
	iters = NewtonIter; resfinal = res; return res < NewtonTol;
}


// 一级的隐式部分（B4-B8）：输入 u_con = X*（Update 的结果）、Qs、Qe；P 进入时为 Newton 初猜，返回时为解；输出 u_con = 本级的终态
bool PressureSubsystem(Grid &u_con, Grid &Qs, Grid &Qe, PField &P, double h, const int rks){
	(void)Qs; (void)Qe; (void)rks;
	const int g = ghost_cells, nx = nxCells, ny = nyCells;
	const double tau = h;
	// 级数据在各级之间复用（只在网格尺寸变化时重新分配）：每级重新分配约 60 个 double/格并清零，大块内存走 mmap，缺页的代价与整个闭合相当。
	// 每级都会重写的量不必清零；只在 CL_MU 格写的阻抗 Z0、Z1 每级重置为 1（与新分配时相同）
	static std::unique_ptr<StageData> SBUF;
	if (!SBUF || SBUF->nx != nx || SBUF->ny != ny) SBUF.reset(new StageData(nx, ny));
	StageData& S = *SBUF;
	std::fill(S.Z0.begin(), S.Z0.end(), 1.0); std::fill(S.Z1.begin(), S.Z1.end(), 1.0);
	S.tau = tau; S.rx = tau / dx; S.ry = tau / dy;
	const int NXT = S.NXT, NYT = S.NYT;
	S.Xs = u_con; S.a2s = *a2_star; S.ubex = *ube_x; S.ubey = *ube_y;
	const Array &phx = *phi_x, &phy = *phi_y;
	// 单位体积内容、供体速度（以下各循环都是逐格或逐面独立的，OpenMP 并行，结果逐位不变）
	bool sok = true;
	// 一维时 y 向虚行的格量与 x 面量都不被读取（面、残量、尺度只在 j = g 行上取）：只算这一行
	const int jr0 = ny == 1 ? g : 0, jr1 = ny == 1 ? g + 1 : NYT;
	#pragma omp parallel for collapse(2) schedule(static) reduction(&&:sok)
	for (int j = jr0; j < jr1; j++) for (int i = 0; i < NXT; i++) {
		const double sa = S.Xs(i, j, IA1) + S.a2s(i, j), rho = S.Xs(i, j, IM1) + S.Xs(i, j, IM2);
		if (!(sa > 0 && rho > 0)) { sok = false; continue; }
		double* c = &S.cont[S.id(i, j) * 6];
		c[0] = S.Xs(i, j, IA1) / sa; c[1] = S.a2s(i, j) / sa; c[2] = S.Xs(i, j, IM1) / sa; c[3] = S.Xs(i, j, IM2) / sa; c[4] = S.Xs(i, j, IE1) / sa; c[5] = S.Xs(i, j, IE2) / sa;
		S.ust(i, j) = S.Xs(i, j, IMU) / rho; S.vst(i, j) = S.Xs(i, j, IMV) / rho;
	}
	if (!sok) return false;
	// 增量供体（V21）：x 面 (i,j) = i+1/2 的左态 = 格 i 的右面值，右态 = 格 i+1 的左面值；alpha2 的面值 = 1 - alpha1 的面值（与显式通量一致）
	{
		const int VK[5] = {IA1, IM1, IM2, IE1, IE2};
		auto put = [&](double* d, const Grid& w, int i, int j) { d[0] = w(i, j, VK[0]); d[1] = 1.0 - w(i, j, VK[0]); for (int q = 1; q < 5; q++) d[q + 1] = w(i, j, VK[q]); };
		#pragma omp parallel for collapse(2) schedule(static)
		for (int j = jr0; j < jr1; j++) for (int i = 0; i < NXT - 1; i++) { put(&S.dxL[S.id(i, j) * 6], *wxR_g, i, j); put(&S.dxR[S.id(i, j) * 6], *wxL_g, i + 1, j); }
		if (ny > 1) {
		#pragma omp parallel for collapse(2) schedule(static)
		for (int j = 0; j < NYT - 1; j++) for (int i = 0; i < NXT; i++) { put(&S.dyL[S.id(i, j) * 6], *wyR_g, i, j); put(&S.dyR[S.id(i, j) * 6], *wyL_g, i, j + 1); }
		}
	}
	auto rhoS = [&](int i, int j) { return S.Xs(i, j, IM1) + S.Xs(i, j, IM2); };
	// B4.1 声学项用冻结声速 c_f（各松弛闭合声速的上界，次特征条件）。在 X* 上计算；X* 不可容许时取 c = inf（该面不加声学项）。
	auto cfz = [&](int i, int j) -> double {
		const double a1 = S.Xs(i, j, IA1), a2 = S.a2s(i, j), m1 = S.Xs(i, j, IM1), m2 = S.Xs(i, j, IM2);
		if (!(a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0)) return INF;
		const double p1 = p_ame(0, a1, m1, S.Xs(i, j, IE1)), p2 = p_ame(1, a2, m2, S.Xs(i, j, IE2));
		if (!(p_admissible(0, a1, m1, p1) && p_admissible(1, a2, m2, p2))) return INF;
		const double c21 = c2_rp(0, m1 / a1, p1), c22 = c2_rp(1, m2 / a2, p2);
		const double c2 = (m1 * c21 + m2 * c22) / (m1 + m2);
		return c2 > 0 ? std::sqrt(c2) : INF; };
	static std::vector<double> cz; cz.resize((size_t)NXT * NYT);
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = jr0; j < jr1; j++) for (int i = 0; i < NXT; i++) cz[S.id(i, j)] = cfz(i, j);
	// 残量尺度（B8）：Wood 体积模量 K = rho c_W^2 + 能量的压力尺度 P_E = sum(|eps_k| + |m_k q_k|) / sum de_k/dp
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = jr0; j < jr1; j++) for (int i = 0; i < NXT; i++) {
		const double a1 = S.Xs(i, j, IA1), a2 = S.a2s(i, j), m1 = S.Xs(i, j, IM1), m2 = S.Xs(i, j, IM2), sa = a1 + a2;
		double K = std::max(PH[0].pi, PH[1].pi), PE = 0.0;
		if (a1 > 0 && a2 > 0 && m1 > 0 && m2 > 0) {
			const double p1 = p_ame(0, a1, m1, S.Xs(i, j, IE1)), p2 = p_ame(1, a2, m2, S.Xs(i, j, IE2));
			if (p_admissible(0, a1, m1, p1) && p_admissible(1, a2, m2, p2)) K = 1.0 / ((a1 / sa) / (m1 / a1 * c2_rp(0, m1 / a1, p1)) + (a2 / sa) / (m2 / a2 * c2_rp(1, m2 / a2, p2)));
			const double dd = dedp_am(0, a1, m1) + dedp_am(1, a2, m2);
			if (dd > 0) PE = (std::fabs(S.Xs(i, j, IE1)) + eps_round(0, a1, m1) + std::fabs(S.Xs(i, j, IE2)) + eps_round(1, a2, m2)) / dd;
		}
		S.Kw[S.id(i, j)] = K + PE;
	}
	// 上一次的压力场（Newton 初猜）与 u*：只用来算二阶限制因子 psi（显式量，残量对 P 仍是线性的）
	static Array Pold(0, 0); if (Pold.nxCells != nx || Pold.nyCells != ny) Pold = Array(nx, ny);
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = 0; j < NYT; j++) for (int i = 0; i < NXT; i++) Pold(i, j) = P.P[M_P(i, j)];
	auto minmod = [](double a, double b) { return a * b <= 0 ? 0.0 : (std::fabs(a) < std::fabs(b) ? a : b); };
	// psi = 1 - (sig_L + sig_R)/(2 D)：与 D 同号且 |sig| <= |D|，故 psi in [0,1]；光滑区 psi = O(dx)，间断处 psi -> 1
	auto psi_face = [&](double xm, double xL, double xR, double xp) {
		if (SpaceOrder < 2) return 1.0;
		const double D = xR - xL; if (D == 0.0) return 1.0;
		if (!LIMITER) return std::min(1.0, std::max(0.0, 1.0 - 0.25 * ((xR - xm) + (xp - xL)) / D));   // 中心斜率：psi D = X^R - X^L
		return 1.0 - 0.5 * (minmod(xL - xm, D) + minmod(D, xp - xR)) / D; };
	// B4.1：u_f = ubar* - tau dP/(rhobar h) - phi psi_P dP/(Z_L + Z_R)，p_f = (rho_R P_L + rho_L P_R)/(rho_L + rho_R) - phi psi_u Z_L Z_R/(Z_L + Z_R) D*，
	// phi = (1 - c tau/h)^+：时间步分辨声波时把声学耗散补足到 Godunov（c tau/h >= 1 时为 0，低马赫/大步长不受影响）
	// B4.1(c)：速度跳跃 D* 取面两端两个顶点散度的平均乘 h（体粘性形式）。x 面 (i+1/2, j)：
	//   D = [2 (u_{i+1,j} - u_{i,j}) + (u_{i+1,j+1} - u_{i,j+1}) + (u_{i+1,j-1} - u_{i,j-1})]/4 + (dx/dy) [(v_{i,j+1} - v_{i,j-1}) + (v_{i+1,j+1} - v_{i+1,j-1})]/4
	// 声学子系统只通过散度看见速度：无散场在这一项的核里，低马赫时它不作用在涡上（逐维的法向跳跃 u_{i+1} - u_i 会，devlog V11）。
	// 一维与 y 向均匀时 D = u_{i+1} - u_i（逐位：0.25 (2d + (d + d)) = d）。
	auto dil_jump = [&](int i, int j, bool ydir) -> double {
		if (!ydir) {
			const double d = S.ust(i + 1, j) - S.ust(i, j);
			if (ny == 1 || j < 1 || j + 1 >= NYT) return d;
			return 0.25 * (2.0 * d + ((S.ust(i + 1, j + 1) - S.ust(i, j + 1)) + (S.ust(i + 1, j - 1) - S.ust(i, j - 1))))
				+ 0.25 * (dx / dy) * ((S.vst(i, j + 1) - S.vst(i, j - 1)) + (S.vst(i + 1, j + 1) - S.vst(i + 1, j - 1)));
		}
		const double d = S.vst(i, j + 1) - S.vst(i, j);
		if (i < 1 || i + 1 >= NXT) return d;
		return 0.25 * (2.0 * d + ((S.vst(i + 1, j + 1) - S.vst(i + 1, j)) + (S.vst(i - 1, j + 1) - S.vst(i - 1, j))))
			+ 0.25 * (dy / dx) * ((S.ust(i + 1, j) - S.ust(i - 1, j)) + (S.ust(i + 1, j + 1) - S.ust(i - 1, j + 1)));
	};
	auto set_face = [&](int iL, int jL, int iR, int jR, double hh, bool ydir, double& ubs, double& kp, double& wL, double& pv) {
		const double rL = rhoS(iL, jL), rR = rhoS(iR, jR);
		const double uL = ydir ? S.vst(iL, jL) : S.ust(iL, jL), uR = ydir ? S.vst(iR, jR) : S.ust(iR, jR);
		ubs = 0.5 * (uL + uR); kp = tau / hh / (0.5 * (rL + rR)); wL = rR / (rL + rR); pv = 0.0;
		const double cL = cz[S.id(iL, jL)], cR = cz[S.id(iR, jR)];
		if (!(std::isfinite(cL) && std::isfinite(cR))) return;
		const double phi = std::max(0.0, 1.0 - tau * std::min(cL, cR) / hh);
		if (phi == 0.0) return;
		const int di = iR - iL, dj = jR - jL, im = iL - di, jm = jL - dj, ip = iR + di, jp = jR + dj;
		double psiP = 1.0, psiU = 1.0;
		if (im >= 0 && jm >= 0 && ip < NXT && jp < NYT) {
			// 压力项：只在所有被输运量都光滑处（共用限制因子 phi > 0，B3）才降到二阶：接触间断处 p、u 连续但阻抗跳跃，phi = 0 保持一阶耗散。
			// 速度项（体粘性，作用在散度上）：只看 u 自身的光滑性。乘共用 phi 时，空化扇区里 alpha、m、eps 的极值使它变成一阶体粘性，
			// 把膨胀的动能变成内能（P-S ±500 m/s 的起步加热、SPB 4.5 的扇区误差，devlog V14）
			const Array& ph = ydir ? phy : phx;
			const double phis = std::min(ph(iL, jL), ph(iR, jR));
			psiP = 1.0 - (1.0 - psi_face(Pold(im, jm), Pold(iL, jL), Pold(iR, jR), Pold(ip, jp))) * phis;
			const Array& V = ydir ? S.vst : S.ust;
			psiU = psi_face(V(im, jm), uL, uR, V(ip, jp));
		}
		const double ZL = rL * cL, ZR = rR * cR;
		kp += phi * psiP / (ZL + ZR); pv = phi * psiU * ZL * ZR / (ZL + ZR) * dil_jump(iL, jL, ydir);
	};
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = jr0; j < jr1; j++) for (int i = 0; i < NXT - 1; i++) set_face(i, j, i + 1, j, dx, false, S.ubsx(i, j), S.kpx(i, j), S.wLx(i, j), S.pvx(i, j));
	if (ny > 1) {
		#pragma omp parallel for collapse(2) schedule(static)
		for (int j = 0; j < NYT - 1; j++) for (int i = 0; i < NXT; i++) set_face(i, j, i, j + 1, dy, true, S.ubsy(i, j), S.kpy(i, j), S.wLy(i, j), S.pvy(i, j));
	}
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = g; j < ny + g; j++) for (int i = g; i < nx + g; i++) {
		const double* c = &S.cont[S.id(i, j) * 6];
		const double rho = S.Xs(i, j, IM1) + S.Xs(i, j, IM2), sa = S.Xs(i, j, IA1) + S.a2s(i, j);
		const double rhoe = (S.Xs(i, j, IEN) - 0.5 * (S.Xs(i, j, IMU) * S.Xs(i, j, IMU) + S.Xs(i, j, IMV) * S.Xs(i, j, IMV)) / rho) / sa;
		S.type[S.id(i, j)] = STEP_TYPE.empty() ? closure_type_of(c[0], c[2], c[3], c[4], c[5], rhoe) : STEP_TYPE[(size_t)(i - g) + (size_t)(j - g) * nx];
		if (S.type[S.id(i, j)] == CL_MU) {
			for (int k = 0; k < 2; k++) { const double a = c[k], m = c[2 + k], pk = p_ame(k, a, m, c[4 + k]), r = m / a; const double Z = r * std::sqrt(std::max(c2_rp(k, r, pk), 1e-30)); (k == 0 ? S.Z0 : S.Z1)[S.id(i, j)] = Z; }
		}
	}
	// 激波预算（B5.2，V16；V21 起为 ±1 的紧致模板）：X* 混合压力在模板上的极值给出局部激波强度，压缩格的物理耗散 Q_k 由 (II.5) 的相 Hugoniot 定；按级冻结
	{
		const bool twoD = ny > 1;
		// 只有压缩且有压力差的格要积分激波剖面，代价集中在少数格：动态分配（每格结果与线程无关）
		#pragma omp parallel for schedule(dynamic, 8)
		for (int q = 0; q < nx * ny; q++) {
			const int i = IG(q % nx), j = IG(q / nx);
			auto cell = [&](int di, int dj, double* st) { stencil_state(S.Xs.cell(i + di, j + dj), st); };
			shock_budget(cell, twoD, &S.Qph[S.id(i, j) * 2]);
		}
	}
	// 隐式：全局 Newton。B5.2 的 Lambda_k(delta) 是 delta 的连续分段仿射函数（分段按 pguess 处的 delta 取），G(P) 连续、分段光滑
	int iters; double res; const bool ok = Newton_Solver(S, P, iters, res);
	if (ok) st_maxres = std::max(st_maxres, res);
	if (!ok) { if (VERBOSE) fprintf(stderr, "  Newton failed: it %d res %.3e tau %.3e\n", iters, res, tau); return false; }
	// 终态
	DIVF.assign((size_t)nx * ny, 0.0);
	bool fok = true; double Tm0 = 0.0, Tm1 = 0.0;
	#pragma omp parallel for schedule(static) reduction(&&:fok) reduction(max:Tm0,Tm1)
	for (int q = 0; q < nx * ny; q++) {
		const int i = q % nx, j = q / nx;
		double Gv; CellOut o; if (!residual_at(S, P, i, j, 0, -1, 0, Gv, &o)) { fok = false; continue; }
		const int I = IG(i), J = IG(j);
		u_con(I, J, IA1) = o.cl.a[0]; u_con(I, J, IM1) = o.cl.m[0]; u_con(I, J, IM2) = o.cl.m[1]; u_con(I, J, IMU) = o.mx; u_con(I, J, IMV) = o.my;
		u_con(I, J, IE1) = o.cl.e[0]; u_con(I, J, IE2) = o.cl.e[1]; u_con(I, J, IEN) = o.E;
		if ((size_t)q < DIVF.size()) DIVF[q] = o.s / tau;
		if (o.cl.a[0] > 1e-9 && std::isfinite(o.cl.T[0])) Tm0 = std::max(Tm0, o.cl.T[0]);
		if (o.cl.a[1] > 1e-9 && std::isfinite(o.cl.T[1])) Tm1 = std::max(Tm1, o.cl.T[1]);
	}
	if (!fok) return false;
	st_Tmax[0] = std::max(st_Tmax[0], Tm0); st_Tmax[1] = std::max(st_Tmax[1], Tm1);
	bound_cond_u(u_con);
	st_stage++;
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
		const int P = OUTPREC, Q = std::min(OUTPREC, 8), R = std::min(OUTPREC, 10);
		fprintf(f, "%.10e %.10e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e\n",
			x_0 + (i - g + 0.5) * dx, y_0 + (j - g + 0.5) * dy, P, o[0], P, o[1], P, o[2], P, o[3], P, o[4], P, o[5], P, o[6], P, o[7], Q, o[8], Q, o[9], R, o[10], R, o[11], P, o[12], P, o[13], Q, o[14]);
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
	std::string Path = vtk_dir + "/" + testname + "/IMEX_" + std::to_string(nxCells) + "x"+ std::to_string(nyCells) + "_RK" + std::to_string(RKs) + "_Order" + std::to_string(SpaceOrder) + "_CFL0" + ss.str();
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


// 压力探针：最近格的混合压力随时间（Pelanti 8.1.2 图 3）
void probe_write(Grid& u_con, bool first){
	const int i = std::min(nxCells - 1, std::max(0, (int)std::floor((PROBE_X - x_0) / dx))), j = std::min(nyCells - 1, std::max(0, (int)std::floor((PROBE_Y - y_0) / dy)));
	double o[15]; cell_diag(u_con.cell(i + ghost_cells, j + ghost_cells), o);
	FILE* f = fopen((out_file + ".probe").c_str(), first ? "w" : "a"); if (!f) return;
	fprintf(f, "%.10e %.10e %.10e %.10e %.10e\n", t, o[0], o[5], o[4], o[10]); fclose(f);   // t p Y2 alpha2 T1
}

int main(int argc, char *argv[]){
#ifdef __GLIBC__
	// 大块内存留在堆里（不走 mmap/munmap，也不归还系统）：每级的临时数组不再每次缺页
	mallopt(M_MMAP_THRESHOLD, 1 << 30); mallopt(M_TRIM_THRESHOLD, 1 << 30); mallopt(M_TOP_PAD, 64 << 20);
#endif
	// ---- 命令行（2DEuler 用全局变量；这里允许 key=value 覆盖）
	double tend_override = -1; bool gotN = false, gotNy = false, gotRelax = false, gotCFL = false, gotRus = false;
	for (int a = 1; a < argc; a++) {
		char key[64]; char val[256]; if (sscanf(argv[a], "%63[^=]=%255s", key, val) != 2) { fprintf(stderr, "bad arg %s\n", argv[a]); return 1; }
		std::string k(key); const double d = atof(val); const bool isinf_val = std::string(val) == "inf";
		if (k == "case") test = (Test)(int)d; else if (k == "N" || k == "nx") { nxCells = (int)d; gotN = true; } else if (k == "ny") { nyCells = (int)d; gotNy = true; }
		else if (k == "cfl") { CFL = d; gotCFL = true; } else if (k == "order") recon = (int)d < 2 ? NoRecon : Minmod; else if (k == "tend") tend_override = d; else if (k == "out") out_file = val;
		else if (k == "relax") { RELAX_LEVEL = (int)d; gotRelax = true; } else if (k == "epsI") EPS_I = d;
		else if (k == "mu") mu_rate = isinf_val ? INF : d; else if (k == "theta") theta_rate = isinf_val ? INF : d; else if (k == "nu") nu_rate = isinf_val ? INF : d;
		else if (k == "verbose") VERBOSE = (int)d;
		else if (k == "rk") RKs = (int)d; else if (k == "cflac") CFL_ac = d; else if (k == "ramp") RAMP = (int)d; else if (k == "tol") NewtonTol = d; else if (k == "dt") DTFIX = d;
		else if (k == "vtk") vtk_dir = val; else if (k == "limiter") LIMITER = (int)d; else if (k == "rus") { RUS_LOCAL = (int)d; gotRus = true; }
		else if (k == "p0") VORTEX_P0 = d; else if (k == "wellprep") WELLPREP = (int)d; else if (k == "init") INIT_ONLY = (int)d; else if (k == "nulaw") NU_LAW = (int)d; else if (k == "nucr") NU_CR = d; else if (k == "facediag") FACEDIAG = (int)d; else if (k == "share") SHARE_MODE = (int)d; else if (k == "acdt") caltimestep = (int)d ? Fullydt : IMEXdt;
		else if (k == "jreuse") JREUSE = (int)d; else if (k == "relaxform") RELAX_FORM = (int)d; else if (k == "unify") UNIFY = (int)d; else if (k == "unifyproj") UNIFY_PROJ = (int)d;
		else if (k == "snap") parse_list(val, SNAP); else if (k == "prec") OUTPREC = (int)d;
		else if (k == "probe") { std::vector<double> v; parse_list(val, v); if (v.size() == 2) { HAS_PROBE = true; PROBE_X = v[0]; PROBE_Y = v[1]; } }
		else { fprintf(stderr, "unknown key %s\n", key); return 1; }
	}
	SpaceOrder = recon == NoRecon ? 1 : 2;

	// ---- 算例：区域、边界、终止时间、状态方程（2DEuler main 里的 switch(test)，定义在 spb6_physics.h）
	TestSetup TS; if (!test_setup(test, TS)) { fprintf(stderr, "unknown case %d\n", (int)test); return 1; }
	// 算例自带的网格与松弛档位（论文的设置）；命令行 N= ny= relax= 覆盖
	if (!gotN && TS.nx > 0) nxCells = TS.nx;
	if (!gotNy && TS.ny > 0) nyCells = TS.ny;
	if (!gotRelax && TS.relax >= 0) RELAX_LEVEL = TS.relax;
	if (!gotCFL && TS.cfl > 0) CFL = TS.cfl;
	if (!gotRus && TS.rus >= 0) RUS_LOCAL = TS.rus;
	// 松弛档位的默认速率（relax >= 1 时热松弛默认瞬时，relax = 2, 3 化学默认瞬时；可由 theta=, nu= 覆盖为有限值）
	if (theta_rate < 0) theta_rate = RELAX_LEVEL >= 1 ? INF : 0.0;
	if (nu_rate < 0) nu_rate = (RELAX_LEVEL == 2 || RELAX_LEVEL == 3) ? INF : 0.0;
	if (NU_LAW) nu_rate = 0.0;   // 有限速率的化学类型；每格的速率由 STEP_NU 给（nu_rate 只是不用的占位）
	if (RELAX_LEVEL > 0 && !mg_any()) build_sat();   // 饱和线表在并行循环之前建好（Tsat 惰性构造不是线程安全的）
	x_0 = TS.x0; x_1 = TS.x1; y_0 = TS.y0; y_1 = TS.y1; time_stop = TS.tend; x_bc = TS.xbc; y_bc = TS.ybc;
	apply_bc_setup(TS);
	if (tend_override > 0) time_stop = tend_override;
	dx = (x_1 - x_0) / nxCells;
	dy = (y_1 - y_0) / nyCells;
	if (RELAX_LEVEL > 0) build_sat();

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells), u_con_n(nxCells, nyCells), Qe(nxCells, nyCells), Qs(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells),w_yL(nxCells, nyCells), w_yR(nxCells, nyCells);
	Array phx(nxCells, nyCells), phy(nxCells, nyCells), a2s(nxCells, nyCells), ubx(nxCells, nyCells), uby(nxCells, nyCells);
	phi_x = &phx; phi_y = &phy; a2_star = &a2s; ube_x = &ubx; ube_y = &uby;
	wxL_g = &w_xL; wxR_g = &w_xR; wyL_g = &w_yL; wyR_g = &w_yR;
	u_initialization(u_con);
	bound_cond_u(u_con);
	double M0[3]; conserved_totals(u_con, M0);
	Grid u_init(nxCells, nyCells); copy(u_init, u_con);   // 算例 56：平移一周后的精确格平均就是初值
	if (INIT_ONLY) { write_profile(u_init, out_file); return 0; }   // init=1：只写精确格平均初值（误差后处理 tests/allmach_table.py 的参考解）
	// 算例 56 的离散良态初值（wellprep=1，devlog V15）：格平均速度的离散散度是 O(h^2)，低马赫时乘 rho c^2 产生长波声波，
	// 粗网格上压力的表观阶降到 1.2-1.3。改用格心流函数 psi_ij 的中心差分 u = 2 + D_y psi，v = 2 - D_x psi：
	// 面速度取两侧格速度的平均（B4.1），格的散度是宽中心差分 D_x u + D_y v = D_x D_y psi - D_y D_x psi = 0（差分算子可交换），逐位为零到舍入。
	// 密度、alpha、相内能仍取格平均；误差仍与精确格平均（u_init）比
	// 算例 57、58（Taylor-Green）同理，流函数 sin x sin y，无平移
	const bool isTG = test == testTaylorGreen1 || test == testTaylorGreen2;
	if ((test == testVortexBB || isTG) && WELLPREP) {
		auto xc = [&](int i) { i = (i + nxCells) % nxCells; return x_0 + (i + 0.5) * dx; };
		auto yc = [&](int j) { j = (j + nyCells) % nyCells; return y_0 + (j + 0.5) * dy; };
		auto psi = [&](double x, double y) { return isTG ? tg_psi(x, y) : vortex_psi(x, y); };
		const double U0 = isTG ? 0.0 : 2.0;
		for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) {
			const double u = U0 + (psi(xc(i), yc(j + 1)) - psi(xc(i), yc(j - 1))) / (2 * dy);
			const double v = U0 - (psi(xc(i + 1), yc(j)) - psi(xc(i - 1), yc(j))) / (2 * dx);
			double* Q = &u_con(IG(i), IG(j), 0); const double rho = Q[IM1] + Q[IM2];
			Q[IMU] = rho * u; Q[IMV] = rho * v; Q[IEN] = Q[IE1] + Q[IE2] + 0.5 * rho * (u * u + v * v);
		}
		bound_cond_u(u_con); conserved_totals(u_con, M0);
	}

	const int plotINTERVAL = 0;
	std::vector<double> IntermediatePlotTimes = {0.25 * time_stop, 0.5 * time_stop, 0.75 * time_stop, time_stop};
	bool isFinalTime = false;
	IntermediatePlotIntermediatePlot(u_con, u_pri, timestep, plotINTERVAL, IntermediatePlotTimes, isFinalTime);

	// Butcher 表（2DEuler）：1 = 后向 Euler，2 = LSDIRK2(2,2,2)，3 = SSP2(3,3,2)；被起步斜坡限制的步用 1
	double Ae[3][3], Ai[3][3];
	auto set_tab = [&](int rk) {
		for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) Ae[a][b] = Ai[a][b] = 0;
		if(rk == 1){
			Ai[0][0] = 1.0;
		}else if(rk == 2){
			double bb = 1.0 - 1.0 / sqrt(2.0);
			double cc = 1.0 / (2.0 * bb);
			Ae[1][0] = cc;
			Ai[0][0] = bb;
			Ai[1][0] = 1 - bb;
			Ai[1][1] = bb;
		}else if(rk == 3){
			double r12 = 0.5, r13=1.0/3.0, r14=1.0/4.0;
			Ae[1][0] = r12;
			Ae[2][0] = r12;
			Ae[2][1] = r12;
			Ai[0][0] = r14;
			Ai[1][1] = r14;
			Ai[2][0] = r13;
			Ai[2][1] = r13;
			Ai[2][2] = r13;
		}
	};

	auto t0 = std::chrono::steady_clock::now();
	if (HAS_PROBE) probe_write(u_con, true);
	PField P{nxCells, nyCells, std::vector<double>((size_t)nxCells * nyCells, 0.0)};
	// P_REF：初始混合压力的最小值（B8）。取最小值而不是平均值：P - P_REF 的量级不超过当地压力，残量的舍入底不因参考值变大
	P_REF = INF;
	for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) P_REF = std::min(P_REF, pbar_of(u_con.cell(IG(i), IG(j))));
	if (!std::isfinite(P_REF)) P_REF = 0.0;
	for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) P.at(i, j) = pbar_of(u_con.cell(IG(i), IG(j))) - P_REF;
	double dtprev = -1;
	bool startup = true;
	Tensor Hs(3, nxCells, nyCells); Grid Ust(nxCells, nyCells); PField Pst = P;   // 级的工作量：每次尝试都整体重写，在循环外分配一次   // 起步斜坡：首次达到物质 CFL 步长之前的步用一阶；之后被 1.2 倍限制的步（拒绝减半之后）仍用二阶
	do{
		// ---- 时间步（caldt）：物质 CFL；起步时首步为声学步长，之后每步至多放大 ramp_factor 倍（2DEuler 的 full_step），被限制的步用一阶
		bool ramping = false;
		if (DTFIX > 0) dt = DTFIX;
		else {
			dt = CFL * caldt(u_con, caltimestep == Fullydt);
			if (RAMP) { if (dtprev < 0) dtprev = CFL_ac * caldt(u_con, true) / ramp_factor; if (dt > ramp_factor * dtprev) { dt = ramp_factor * dtprev; ramping = true; } }
			if (!std::isfinite(dt)) dt = CFL_ac * caldt(u_con, true);
		}
		dt = std::min(dt, time_stop - t); dt = snap_limit(t, dt);
		if (!(dt > 1e-12 * time_stop)) { fprintf(stderr, "FAILED: time step underflow (dt = %.3e) at t = %.6e (step %ld)\n", dt, t, st_steps); write_profile(u_con, out_file); return 2; }

		bool done = false;
		for (int attempt = 0; attempt < 30 && !done; attempt++) {
			const bool first_order = RAMP && ramping && startup;
			const int rk = first_order ? 1 : RKs; if (rk == 1) st_rk1++;
			set_tab(rk);
			if (RELAX_LEVEL > 0) RelaxationActivation(u_con);   // B6.7
			copy(u_con_n,u_con);
			Pst = P;
			bool ok = true;
			for(int rks = 0; rks < rk && ok; rks++){
				copy(Qe,u_con_n);
				copy(Qs,u_con_n);
				// 只算内部格：虚格紧接着由 bound_cond_u 整体重写（Hs 的虚格也只在这里被读）
				for(int r = 0; r < rks; r++){
					#pragma omp parallel for collapse(2) schedule(static)
					for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
						for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
							for(int k = 0; k < NUM_VARS; k++){
								Qe(i,j,k) = Qe(i,j,k) + Ae[rks][r] * dt * Hs(r,i,j,k);
								Qs(i,j,k) = Qs(i,j,k) + Ai[rks][r] * dt * Hs(r,i,j,k);
							}
						}
					}
				}
				bound_cond_u(Qe);
				bound_cond_u(Qs);
				adt =  Ai[rks][rks] * dt;

				TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR);

				ok = Update(Ust, Qe, Qs, w_xL,w_xR,w_yL,w_yR);                      // 显式（对流）部分：X*

				ok = ok && PressureSubsystem(Ust, Qs, Qe, Pst, adt, rks);           // 隐式（压力 + 闭合）部分
				if (ok && FACEDIAG) for (int q = 0; q < 3; q++) FD_TOT[rks == 0 ? 0 : 1][q] += FD_CUR[q];

				if (ok) {
					#pragma omp parallel for collapse(2) schedule(static)
					for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
						for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
							for(int k = 0; k < NUM_VARS; k++){
								Hs(rks,i,j,k) = (Ust(i,j,k) - Qs(i,j,k)) / adt;
							}
						}
					}
				}
			}
			if (ok) { copy(u_con, Ust); P = Pst; done = true; }
			// 不可容许的处理（B10）：步长减半重算，仍用 LSDIRK2（V23）
			if (!done) { st_reject++; dt *= 0.5; if (VERBOSE) fprintf(stderr, "reject: t %.4e new dt %.3e\n", t, dt); }
		}
		if (!done) { fprintf(stderr, "FAILED at t = %.6e (step %ld)\n", t, st_steps); write_profile(u_con, out_file); return 2; }
		t = t + dt; dtprev = dt; st_steps++; timestep++; if (!ramping) startup = false;
		if (SNAPI < SNAP.size() && t >= SNAP[SNAPI] * (1 - 1e-12)) { write_profile(u_con, snap_name(out_file, SNAPI)); SNAPI++; }
		if (HAS_PROBE) probe_write(u_con, false);
		if (VERBOSE > 1 && st_steps % 50 == 0) fprintf(stderr, "step %ld t %.4e dt %.3e newton it %ld\n", st_steps, t, dt, st_newton_it);
		IntermediatePlotIntermediatePlot(u_con, u_pri, timestep, plotINTERVAL, IntermediatePlotTimes, t >= time_stop * (1 - 1e-12));
	}while(t < time_stop * (1 - 1e-12));

	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double M1[3]; conserved_totals(u_con, M1);
	if (RELAX_LEVEL > 0) RelaxationActivation(u_con);   // 输出前同样做激活投影（P-S：输出的是松弛后的态）
	write_profile(u_con, out_file);
	if (test == testVortexBB) {   // L2 = sqrt(sum e^2 dx dy)（论文表 1 的量：alpha2 rho2、alpha2 rho2 u、alpha2；另加混合密度 rho、速度矢量 |u - u0|、混合压力 p）
		double e[6] = {0, 0, 0, 0, 0, 0};
		for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) {
			const double* Q = u_con.cell(IG(i), IG(j)); const double* Q0 = u_init.cell(IG(i), IG(j));
			const double u1 = Q[IMU] / (Q[IM1] + Q[IM2]), u0 = Q0[IMU] / (Q0[IM1] + Q0[IM2]);
			const double r = Q[IM1] + Q[IM2], r0 = Q0[IM1] + Q0[IM2];
			const double vel = std::hypot(u1 - u0, Q[IMV] / r - Q0[IMV] / r0);
			const double d[6] = {Q[IM2] - Q0[IM2], Q[IM2] * u1 - Q0[IM2] * u0, Q0[IA1] - Q[IA1], pbar_of(Q) - pbar_of(Q0), r - r0, vel};
			for (int k = 0; k < 6; k++) e[k] += d[k] * d[k] * dx * dy;
		}
		printf("L2 p0 %g N %d | a2r2 %.4e a2r2u %.4e a2 %.4e | rho %.4e vel %.4e p %.4e\n", VORTEX_P0, nxCells, std::sqrt(e[0]), std::sqrt(e[1]), std::sqrt(e[2]), std::sqrt(e[4]), std::sqrt(e[5]), std::sqrt(e[3]));
	}
	if (isTG) {   // BB 4.4：混合速度对不可压定常解（格平均）的误差 L1 = sum |e| dx dy，L2 = sqrt(sum |e|^2 dx dy)，Linf；e = (u, v) 之差的模。
		// 另给压力的二阶脉动 p - p0 对 (cos 2x + cos 2y)/4 的 L2、密度偏差 max |rho - 1|（低马赫一致性：O(M^2)）
		double e1 = 0, e2 = 0, einf = 0, ep = 0, drho = 0, er = 0;
		for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) {
			const double* Q = u_con.cell(IG(i), IG(j)); const double* Q0 = u_init.cell(IG(i), IG(j));
			const double r = Q[IM1] + Q[IM2], r0 = Q0[IM1] + Q0[IM2];
			const double du = Q[IMU] / r - Q0[IMU] / r0, dv = Q[IMV] / r - Q0[IMV] / r0, e = std::hypot(du, dv);
			e1 += e * dx * dy; e2 += e * e * dx * dy; einf = std::max(einf, e);
			const double d = (pbar_of(Q) - VORTEX_P0) - (pbar_of(Q0) - VORTEX_P0); ep += d * d * dx * dy; er += (r - r0) * (r - r0) * dx * dy; drho = std::max(drho, std::fabs(r - 1.0));
		}
		// BB 表 4 的量：混合速度的宽中心差分散度 div_h u = (u_{i+1} - u_{i-1})/(2 dx) + (v_{j+1} - v_{j-1})/(2 dy)（周期），应为 O(M^2)
		double d1 = 0, d2 = 0, dinf = 0;
		auto vel = [&](int i, int j, int k) { i = (i + nxCells) % nxCells; j = (j + nyCells) % nyCells; const double* Q = u_con.cell(IG(i), IG(j)); return Q[k] / (Q[IM1] + Q[IM2]); };
		for (int j = 0; j < nyCells; j++) for (int i = 0; i < nxCells; i++) {
			const double dv = (vel(i + 1, j, IMU) - vel(i - 1, j, IMU)) / (2 * dx) + (vel(i, j + 1, IMV) - vel(i, j - 1, IMV)) / (2 * dy);
			d1 += std::fabs(dv) * dx * dy; d2 += dv * dv * dx * dy; dinf = std::max(dinf, std::fabs(dv));
		}
		double f1 = 0, f2 = 0, finf = 0;   // 格式自己的面散度（紧致，投影所约束的量）
		for (double dv : DIVF) { f1 += std::fabs(dv) * dx * dy; f2 += dv * dv * dx * dy; finf = std::max(finf, std::fabs(dv)); }
		printf("TG face-div L1 %.4e L2 %.4e Linf %.4e\n", f1, std::sqrt(f2), finf);
		printf("TG p0 %g M %.3e N %d | div L1 %.4e L2 %.4e Linf %.4e | u-err L1 %.4e L2 %.4e Linf %.4e | p2 L2 %.4e | rho L2 %.4e max|rho-1| %.3e\n", VORTEX_P0, 1.0 / std::sqrt(1.4 * VORTEX_P0), nxCells,
			d1, std::sqrt(d2), dinf, e1, std::sqrt(e2), einf, std::sqrt(ep), std::sqrt(er), drho);
	}
	printf("case %d mode 1 order %d rk %d N %d CFL %g | steps %ld (rk1 %ld) rejects %ld floor %ld | stages %ld newton it %ld (%.2f/stage) maxres %.2e Tmax %.4g %.4g | closure fails %ld | dM1/M %.2e dM2/M %.2e dE/E %.2e | wall %.2f s\n",
		(int)test, SpaceOrder, RKs, nxCells, CFL, st_steps, st_rk1, st_reject, st_floor, st_stage, st_newton_it, st_stage ? (double)st_newton_it / st_stage : 0.0, st_maxres, st_Tmax[0], st_Tmax[1], st_cl_fail,
		(M1[0] - M0[0]) / M0[0], (M1[1] - M0[1]) / M0[1], (M1[2] - M0[2]) / M0[2], wall);
	if (UNIFY == 2) printf("unify2 | calls / rate-inf closure failed / max difference: p %ld %ld %.2e | pT %ld %ld %.2e | pTg %ld %ld %.2e\n", st_un_n[0], st_un_fail[0], st_un_dmax[0], st_un_n[1], st_un_fail[1], st_un_dmax[1], st_un_n[2], st_un_fail[2], st_un_dmax[2]);
	if (FACEDIAG) printf("facediag | stage Qe=Qs: faces %ld violations %ld (admissible stencil %ld) | later stages: faces %ld violations %ld (admissible stencil %ld)\n", FD_TOT[0][0], FD_TOT[0][1], FD_TOT[0][2], FD_TOT[1][0], FD_TOT[1][1], FD_TOT[1][2]);
	return 0;
}
