#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <functional>
#include <iomanip>
#include <random>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Core>
#include <filesystem>
#include <sstream>

#define NUM_VARS 8
enum Formula{EnergyFormula, PressureFormula};
enum MagField{StaggerCT, UnStaggerCT, PotentialA, NoDivFree};
enum Test{test0, test1, test2, test3, test4, test5, test6, test7,test8,test9,test10, shocktest, convergencetest,test2d0, test1y, testOrszag_Tang, testblastwave, testblastwave1D,testRotor, testFieldLoop, testKHInstability, testdoubleshear, testRK_Riemann, testPlasmaEquilibrium, testAdvectedPlasmaEquilibrium, testMHD_vortex, testHallReconnection, testHarrisHall};
enum Reconstruction{NoRecon, Minmod, VanLeer, ENO3, WENO5};
enum shock_threshold_Left{one, val_left};
enum shock_threshold_Right{zero, val_right};
enum shock_detector{Yes_add, No_First, No_Second};
enum SolvePressure{no_correction, once_correction, all_correction};
enum dtCalculation{Fullydt, Alfvendt, IMEXdt};
enum x_boundary_condition{x_Transmissive, x_Reflective, x_Periodic, x_Dirichlet};
enum y_boundary_condition{y_Transmissive, y_Reflective, y_Periodic, y_Dirichlet};
enum Add_Pvisc{MaterialSpeedp, AlfvenSpeedp, Noaddvisp};
enum Energy_Update{ConservationofU, DefinitionofU};
enum Flux_Diss{Local_Diss, Global_Diss};

dtCalculation caltimestep = IMEXdt;
Formula formula = PressureFormula;
MagField magfield = PotentialA;
Add_Pvisc add_pvisc = MaterialSpeedp;
Flux_Diss flux_diss = Global_Diss;
bool GLF = true;
bool add_bvisc = true;    // 磁 A 子系统数值粘性:论文一直需要;系数=面平均快磁声速(对齐 CASTOR add_bvisc=.TRUE.)
Energy_Update U_update = DefinitionofU;
Test test = testHarrisHall;
const int nxCells = 128;
const int nyCells = 64;
const int nzCells = 1;
const int ghost_cells = 3;
const int RKs = 2;
double CFL = 0.1;
const int full_step = 10;     // 对齐 CASTOR:特征值 magnetosonic 只在 step 1,之后材料速度(硬编码 timestep==1)
double A0 = 1e-3;
Reconstruction recon = Minmod;

int SpaceOrder;
x_boundary_condition x_bc;
y_boundary_condition y_bc;
SolvePressure PressureSolver = all_correction;

bool add_RHSpvisc = false;
bool add_RHSBvisc = false;
bool Energy_Modify = false;
bool   add_resistive     = true;   // 物理电阻 η 开关(区别于 GLF/nu_x,nu_y 的数值粘性)
double eta_resistive      = 0;    // 电阻系数 η
bool   add_hyperresistive = true;  // 超电阻 η_H 开关(可选,见下面的边界警告)
double eta_hyper          = 6.25 * 1e-5;    // 超电阻系数 η_H
double Hall_di = 1.0;

// ===================== Hall MHD (IMEX) =====================================================
// Non-ideal induction terms follow the Dematte et al. VRMHD pattern (Eqs. 59/61): they enter
// BOTH the implicit magnetic subsystem (so the B feeding the momentum/pressure solve is
// non-ideal-consistent) AND the constrained-transport electric field (which produces the
// final, divergence-free B).  This is not double counting: the CT update restarts from the
// stage accumulator Bx_p0_Qs / By_0p_Qs, not from the cell-centred solve.
//
// Induction equation actually advanced (mu0 = 1, J = curl B), = report Eq. (7.3):
//   dB/dt = curl( u x B - (d_i/rho) J x B ) + eta lap B - eta_H lap^2 B
// i.e. the out-of-plane electric field is
//   E_z = (v_y B_x - v_x B_y) + eta J_z + (d_i/rho)(J x B)_z - eta_H lap J_z .
// NOTE (paper erratum): report Eq. (3.31) prints E_H = -(J x B)/(ne) while Eq. (3.32) uses the
// same "B* - dt curl E" form as the resistive Eq. (3.29).  Those two are inconsistent; the
// correct Hall electric field is E_H = +(J x B)/(ne).  The sign used here is the correct one
// (and agrees with the sign already used in this code's explicit A** Hall term).
bool   add_hall       = true;    // master switch for the Hall term
bool   hall_implicit  = true;    // true: implicit in the magnetic subsystem (recommended)
                                 // false: legacy explicit source operator, comparison only
bool   hall_in_emf    = true;    // put E_H (+ eta J + hyper) into the CT corner EMF.
                                 // Without this the Hall field NEVER reaches the in-plane
                                 // staggered components and no Hall reconnection is possible.
double hall_omega_cap = 8.0;     // ACCURACY cap on h*omega_whistler(grid); 0 disables.
                                 // Hall is implicit so this is not a stability limit; the
                                 // explicit RK4 bound for reference is 2.83 (report Sec 7.1.2).
double hall_kick_amp  = 1.125;   // delta A_z amplitude psi_0.  The survey (report Sec 7.2.4,
                                 // reference code case 27) uses the SINGLE-KICK first-peak
                                 // protocol psi_0 = 0.75 * B_eff, B_eff = 2 B1 B2/(B1+B2):
                                 // a deliberately large kick that drives reconnection at once
                                 // so the first rate maximum can be sampled, NOT the 0.01*B0
                                 // convention of Sec 7.3, which belongs to the production
                                 // Double-Harris runs.  For b = 3: B_eff = 1.5 -> psi_0 = 1.125,
                                 // giving |delta B|max = psi_0 * k_x = 0.5625.
bool   hall_diag_on   = true;    // reconnection diagnostics of report Sec 7.4
bool   hall_snap_on   = true;    // append 2D fields (Az, Bz, Jz, vx, rho) at each diagnostic
                                 // output to hall_fields.bin, for plotting
double hall_diag_dt   = 0.5;     // analysis-output cadence Delta t_out (Table 7.2)
double hall_xjump_cells = 5.0;   // mask a rate sample if the X cell moved more than this many
                                 // cells between neighbouring outputs (report Sec 8.4.2)

// Upstream asymptotic states of the testHallReconnection initial condition, used by the
// Cassak-Shay normalisation of report Eqs. (7.34)-(7.37).  Side 1 = weak field side.
//   Bx_bg(y) = -1 + 2 tanh(y/lambda):  y -> +inf gives |B1| = 1,  y -> -inf gives |B2| = 3
//   rho(y)   = 1.5 + 0.5 tanh(y/lambda): y -> +inf gives rho1 = 2, y -> -inf gives rho2 = 1
// i.e. B1/B2 = 1/3, rho1/rho2 = 2/1, b = 3  ==  case R3 of report Table 7.1.
double hall_B1 = 1.0, hall_B2 = 3.0, hall_rho1 = 2.0, hall_rho2 = 1.0;
// ---- Test 12 of the report (Table 4.6): symmetric Harris sheet, Hall MHD, GEM setup ------
// Exact Harris kinetic equilibrium, uniform temperature T = p_inf/rho_inf = 0.5:
//   Bx(y)  = B0 tanh(y/lam)
//   rho(y) = sech^2(y/lam) + rho_inf          (rho_inf = 0.2)
//   p(y)   = p_inf + (B0^2/2) sech^2(y/lam)   (p_inf  = 0.1, so p + Bx^2/2 = 0.6 = const)
// GEM single-mode seed (Birn et al. 2001):  dAz = psi0 cos(kx x) cos(ky y),
//   kx = 2 pi / Lx = 0.5, ky = pi / Ly = 0.5, psi0 = 0.1 -> max|dBy| = psi0 kx = 0.05.
// eta = 5e-3 (S ~ 2500), d_i = 1, NO hyper-resistivity, rho_floor = 0.02.
double h12_B0 = 1.0, h12_lam = 0.5, h12_rho_inf = 0.2, h12_p_inf = 0.1, h12_psi0 = 0.1;
double h12_kx = 0.5, h12_ky = 0.5;
inline double h12_Az_bg(double y){ return h12_B0*h12_lam*log(cosh(y/h12_lam)); }
inline double h12_Az(double x, double y){
	return h12_Az_bg(y) + h12_psi0*cos(h12_kx*x)*cos(h12_ky*y);
}

double hall_rho_floor = 0.02;    // density floor used when forming corner 1/rho (Table 7.2)
double hall_ptot   = 5.5;        // total (thermal + magnetic) pressure of the equilibrium
double hall_lambda = 0.5;        // sheet half-thickness (= d_i/2 for the survey config)
double hall_yNL    = 0.2747;     // lambda*atanh((B2-B1)/(B2+B1))

long   g_gmres_matvec = 0;       // matrix-vector products in the last GMRES solve
double g_dt_hall = 0.0, g_dt_ideal = 0.0;   // diagnostics from caldt
// ============================================================================================

double GMREStol = 1e-14;
double GMRESiter = 1e32;
int GMRESrest = 30;
bool add_restart = false;


shock_detector add_shock_detector = No_Second;
double shock_kappa1 = 0.001;
double shock_threshold = 0.9;
shock_threshold_Left threshold_Left = val_left;
shock_threshold_Right threshold_Right = zero;

int nxPoints = nxCells + 1, nyPoints = nyCells + 1, nzPoints = nzCells + 1;
double t, time_stop, dt, adt, adt2, fully_dt, imex_dt;
int timestep = 0;
double x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1, B_max_x = 0, B_max_y = 0, U_max_x = 0,U_max_y = 0,a_max_x = 0,a_max_y = 0, a_max = 0;
std::string Path;

double dx = (x_1 - x_0)/(nxPoints-1);
double dy = (y_1 - y_0)/(nyPoints-1);
double dz = (z_1 - z_0)/(nzPoints-1);
double ch_x,ch_y, cp_x, cp_y, cA_x, cA_y, cA, vx_max, vy_max, ch, v_max;


// ===================== 压力保正:Casulli 嵌套 Newton(piecewise-linear V) =====================
// 思路:把 pressure subsystem 的对角项 V(p)=p/(gamma-1) 在 p<p_floor 处改成大斜率 Mbig 的直线,
//       于是 V 全局非线性 -> 用 Newton 外层 + GMRES 内层求解;大斜率把 p 弹回 p_floor 之上。
bool   floor_on   = true;          // 总开关:false 则完全退回你原来的线性单次 GMRES
double p_floor    = 1e-7;          // 压力地板
double Mbig_fac   = 1e6;           // 大斜率倍率(相对 1/(gamma-1));GMRES 病态就调小
int    Pnewton_max = 20;           // Newton 最大迭代
double Pnewton_tol = 1e-10;        // Newton 残差容差
Eigen::VectorXd p_freeze;          // 当前 Newton 迭代值,供 Jacobian(=AX)读取斜率

// ============================================================================================
// =====================  EOS 抽象层(Ideal Gas / Redlich-Kwong)  =============================
// 论文核心:EOS 的全部非线性只通过对角项 V(p)=rho*e(rho,p) 及其斜率 dV/dp 进入隐式系统,
// 耦合矩阵 T(焓 Laplacian)与 EOS 无关。因此移植只需把 V(p)、声速、con<->pri 改成 EOS-aware。
// 理想气体路径保持逐字节不变(Vbase_slope==1/(gamma-1), 声速==sqrt(gamma p/rho))。
enum EOSType { IDEAL_GAS, REDLICH_KWONG };
EOSType eos_type = IDEAL_GAS;                 // 顶层开关;RK 测试在 main 里自动置 REDLICH_KWONG

// ---- Ideal Gas 参数----
double gamma_plasma = 5.0/3.0;
// ---- Redlich-Kwong 参数----
double m_alfa = 0.5, m_b = 0.5, m_R = 0.4, m_cv = 1.0;

// ---- RK 热/卡 EOS 的 Newton 残差(由 1D 代码移植)----
//   热 EOS:  p = m_R T/(v-b) - alfa/sqrt(T)/v/(v+b),   v = 1/rho
//   卡 EOS:  e = m_cv T + (3/2) alfa/b * ln(v/(v+b)) / sqrt(T)
inline double RK_fNR_T (double T,double v,double p){ return m_R*T/(v-m_b) - m_alfa/sqrt(T)/v/(v+m_b) - p; }
inline double RK_dfNR_T(double T,double v,double /*p*/){ return m_R/(v-m_b) + 0.5*m_alfa/T/sqrt(T)/v/(v+m_b); }
inline double RK_fNR_p (double T,double v,double e){ double U=log(v/(v+m_b)); return m_cv*T + 1.5*m_alfa/m_b*U/sqrt(T) - e; }
inline double RK_dfNR_p(double T,double v,double /*e*/){ double U=log(v/(v+m_b)); return m_cv - 0.75*m_alfa/m_b*U/T/sqrt(T); }


// 由 (rho,p) 解温度 T
inline double REDKWONG_computeTfromRhoAndP(double rho, double p){
    double v = 1.0/rho;
    double T_old, T_new = (1.0/rho - m_b)*(p + m_alfa*rho*rho)/m_R;
    const double eps = 1e-14;
    do{ T_old = T_new; T_new = T_old - RK_fNR_T(T_old,v,p)/RK_dfNR_T(T_old,v,p);
    }while(fabs(T_new - T_old)/fabs(T_old) > eps);
    return T_new;
}
// e(rho,p)
inline double REDKWONG_e_from_rho_p(double rho, double p){
    double v = 1.0/rho;
    double T = REDKWONG_computeTfromRhoAndP(rho,p);
    double U = log(v/(v+m_b));
    return m_cv*T + 1.5*m_alfa/m_b*U/sqrt(T);
}
// p(rho,e):先由卡 EOS 解 T,再代入热 EOS
inline double REDKWONG_p_from_rho_e(double rho, double e){
    double v = 1.0/rho;
    double T_old, T_new = (e + rho*m_alfa)/m_cv;
    const double eps = 1e-14;
    do{ T_old = T_new; T_new = T_old - RK_fNR_p(T_old,v,e)/RK_dfNR_p(T_old,v,e);
    }while(fabs(T_new - T_old)/fabs(T_old) > eps);
    double a = m_alfa/sqrt(T_new);
    return m_R*T_new/(v-m_b) - a/v/(v+m_b);
}
// RK 声速 a = sqrt(dp/drho|_s)(由 1D REDKWONG_computeSoundSpeedFromEoS 移植)
inline double REDKWONG_soundspeed(double rho, double p){
    double v = 1.0/rho;
    double T = REDKWONG_computeTfromRhoAndP(rho,p);
    double U = log(v/(v+m_b));
    double p_T = m_R/(v-m_b) + m_alfa/(2.0*T*sqrt(T)*v*(v+m_b));
    double e_T = m_cv - 0.75*m_alfa/m_b*U/(T*sqrt(T));
    double p_v = -m_R*T/pow((v-m_b),2) + m_alfa*(2.0*v+m_b)/(sqrt(T)*v*v*pow((v+m_b),2));
    double e_v = 1.5*m_alfa/m_b*(1.0/v - 1.0/(v+m_b))/sqrt(T);
    double a2 = p/(rho*rho)*(p_T/e_T) + v*v*(e_v - e_T*(p_v/p_T));
    return sqrt(a2);
}

// ---- 统一 EOS 接口:其余代码全部经由这些函数,EOS 切换只改 eos_type ----
inline double eos_p_from_rho_e(double rho, double e){
    return (eos_type==IDEAL_GAS) ? (gamma_plasma - 1.0)*rho*e
                                 : REDKWONG_p_from_rho_e(rho, e);
}
inline double eos_e_from_rho_p(double rho, double p){
    return (eos_type==IDEAL_GAS) ? p/((gamma_plasma - 1.0)*rho)
                                 : REDKWONG_e_from_rho_p(rho, p);
}
inline double eos_sound_a(double rho, double p){
    return (eos_type==IDEAL_GAS) ? sqrt(gamma_plasma*p/rho)
                                 : REDKWONG_soundspeed(rho, p);
}
inline double eos_sound_a2(double rho, double p){
    double a = eos_sound_a(rho, p); return a*a;
}
// 快磁声速 cf:用 a^2 通用表达,Bn 为法向磁场分量
inline double eos_fast_speed(double rho, double p, double Bn, double B2){
    double a2 = eos_sound_a2(rho, p);
    double t  = a2 + B2/rho;
    double disc = t*t - 4.0*a2*Bn*Bn/rho;     // 4 gamma p Bn^2/rho^2 = 4 a^2 Bn^2/rho
    if(disc < 0) disc = 0;
    return sqrt(0.5*(t + sqrt(disc)));
}

// ===================== 压力保正:Casulli 嵌套 Newton(piecewise V) =====================
// V(p) = rho*e(rho,p) 为隐式压力系统的对角项。理想气体 V(p)=p/(gamma-1) 不依赖 rho;
// RK 的 V(p)=rho*e(rho,p) 依赖 rho,故 V_mod/V_slope 改为 rho-aware。
// 保正:在物理斜率之下设一个大斜率 M(p<p_floor 切换),把 p 弹回 p_floor 之上(论文 (58)-(59))。

// 物理 V(p)=rho*e(rho,p)
inline double Vbase(double p, double rho){
    return (eos_type==IDEAL_GAS) ? p/(gamma_plasma - 1.0)
                                 : rho*REDKWONG_e_from_rho_p(rho, p);
}
// 物理 dV/dp = rho*de/dp:理想气体解析,RK 用中心差分(对齐 1D, small=1e-8)
inline double Vbase_slope(double p, double rho){
    if(eos_type==IDEAL_GAS) return 1.0/(gamma_plasma - 1.0);
    const double d = 1e-8;
    double ep = REDKWONG_e_from_rho_p(rho, p + d);
    double em = REDKWONG_e_from_rho_p(rho, p - d);
    return rho*(ep - em)/(2.0*d);
}
inline double V_mod(double p, double rho){                  // 带保正地板的 V(p)
    if(!floor_on) return Vbase(p, rho);
    double M = Mbig_fac/(gamma_plasma - 1.0);               // 大斜率(与理想气体一致的量纲)
    return (p >= p_floor) ? Vbase(p, rho)
                          : Vbase(p_floor, rho) + M*(p - p_floor);
}
inline double V_slope(double p, double rho){                // dV/dp(Newton Jacobian 对角)
    if(!floor_on) return Vbase_slope(p, rho);
    double M = Mbig_fac/(gamma_plasma - 1.0);
    return (p >= p_floor) ? Vbase_slope(p, rho) : M;
}
// AX 算子模式开关:false=线性载体算子 L(对角=1/(gamma-1));true=Newton Jacobian(对角=V_slope)
bool g_AX_jacobian = false;

// --- Dirichlet ambient state (blast wave); E/A computed with gamma_plasma at runtime ---
double rho_amb = 1.0, p_amb = 0.1, Bx_amb = 70.0/sqrt(4.0*M_PI);
bool g_bc_prim = false;   // true while bound_cond_u is applied to a PRIMITIVE grid (slot 4 = p, not E)
double g_fixed_dt = 0.0;   // 0 = 用真实 IMEX dt(生产);>0 仅为对 CASTOR 强制固定 dt
inline double E_amb_val(){ return p_amb/(gamma_plasma - 1.0) + 0.5*Bx_amb*Bx_amb; } // ambient total energy
// ============================================================================================


struct Tensor
{
	double* data;
	int nxCells;
	int nyCells;
	int nRKS;
	
	// Access to data as l−value
	
	Tensor(int RKS,int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		nRKS = RKS;
		data = new double[nRKS * (nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells) * NUM_VARS];
	}
	
	double operator()(int r,int i, int j, int k)const
	{
		return data[r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int r,int i, int j, int k){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	~Tensor() {
		delete[] data;
	}
};

struct Grid
{
	double* data;
	int nxCells;
	int nyCells;
	
	// Access to data as l−value
	
	Grid(int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		data = new double[(nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells) * NUM_VARS];
	}
	
	double operator()(int i, int j, int k)const
	{
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int i, int j, int k){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	~Grid() {
		delete[] data;
	}
};

struct Array
{
	double* data;
	int nxCells;
	int nyCells;
	
	// Access to data as l−value
	
	Array(int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		data = new double[(nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells)];
	}
	
	double operator()(int i, int j)const
	{
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[i + j * (nxCells + 2 * ghost_cells)];
	}
	double& operator()(int i, int j){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[i + j * (nxCells + 2 * ghost_cells)];
	}
	~Array() {
		delete[] data;
	}
};

class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherit its traits:
		template <>
		struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> > {};
	}  // namespace internal
}  // namespace Eigen

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int StorageIndex;
	enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };
	
	Index rows() const { return mp_mat->rows(); }
	Index cols() const { return mp_mat->cols(); }
	
	template <typename Rhs>
	Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}
	
	// Custom API:
	MatrixReplacement() : mp_mat(0) {}
	
	void attachMyMatrix(const SparseMatrix<double>& mat) { mp_mat = &mat; }
	const SparseMatrix<double> my_matrix() const { return *mp_mat; }
	
private:
	const SparseMatrix<double>* mp_mat;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
	namespace internal {
		
		template <typename Rhs>
		struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape,
		GemvProduct>  // GEMV stands for matrix-vector
		: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> > {
			typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;
			
			template <typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha) {
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				eigen_assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);
				
				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):
				for (Index i = 0; i < lhs.cols(); ++i) dst += rhs(i) * lhs.my_matrix().col(i);
			}
		};
		
	}  // namespace internal
}  // namespace Eigen

using namespace std;
double getAlfvenSpeed(double *pri);
void copystate(double *state_1, double *state_0);
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void f_pressure(double *f, double *con);
void f_magnetic(double *f, double *con);
void g_convective(double *g, double *con);
void g_pressure(double *g, double *con);
void g_magnetic(double *g, double *con);
void f_full_con(double *f, double *con);
void g_full_con(double *g, double *con);
void ini_con2pri(double *pri, double *con);
void ini_pri2con(double *con, double *pri);
void upri2ucon(Grid &u_con, Grid &u_pri);
void ucon2upri(Grid &u_pri,Grid &u_con);
void half_time(double *halfL, double *halfR, double *L, double *R);
double shock_filter(double val);
void MagneticFlux_x(double *Q, double *flux);
void MagneticFlux_y(double *Q, double *flux);
void f_KineticEnergyFlux(double *Q, double *flux);
void g_KineticEnergyFlux(double *Q, double *flux);
void MagneticFlux_Update_Energy(Grid &u_con, Grid &u_test);
void MagneticFlux_Pressure_x(double *Q, double *flux);
void MagneticFlux_Pressure_y(double *Q, double *flux);
void MagneticFlux_Update_Pressure(Grid &u_con, Grid &u_test, Grid &Qe_pri);
void RusanovK_x(double *QL, double *QR, double *flux);
void RusanovK_y(double *QL, double *QR, double *flux);
void RusanovK_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
void RusanovK_Pressure_x(double *QL, double *QR, double *flux);
void RusanovK_Pressure_y(double *QL, double *QR, double *flux);
void RusanovK_Update_Pressure(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
double ini_calch_x(double *pri);
double ini_calch_y(double *pri);
void RusanovC_x(double *QL, double *QR, double *flux);
void RusanovC_y(double *QL, double *QR, double *flux);
void RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h);
void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
double calecon(double *con);
double calepri(double *pri);
double calcf(double* pri);
vector<double> ConvergenceTestFunction(double x);
vector<double> ConvergenceTestExact(double x);
vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<vector<double>(double)> func, int nGaussPoints);
void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p);
void bound_cond_u(Grid &u);
void bound_cond_A(Grid &A);
double caldt(Grid &u_pri);
void update(Grid &u_con,Grid &flux, char dir);
void calch(Grid &u_pri);
void calca(Grid &u_pri);
double ini_calch(double *u_pri);
void calcp(Grid &u_pri);
void store(Grid &u_con, Array& vxn,Array& vyn, Array& Bxn, Array& Byn);
void bound_cond_array(Array &u);
double calenthalpy(double *con);
double calkinetic(double *con);
double calmagenergy(double *con);
void MagneticEigenvalues(Grid &u_pri);
void EnergyEigenvalues(Grid &u_pri);
void Staggered_CT(Grid& u_con, Array &Bx_p0, Array &By_0p, double h);
void A_initial(Grid &A);
inline double minmod(double a, double b);
void Unstaggered_CT(Grid& u_con, Grid& A,Grid& A_s, double h);
int M_P(int i, int j);
int M_B(int i, int j, int k);
void Magnetic_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void Magnetic_Subsystem_potential_AX(Eigen::VectorXd& Ax_out, Eigen::VectorXd& x,Grid& u_con, Grid& Qe, Grid& Qs, double h);
void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void Pressure_NestedNewton(Eigen::VectorXd& p, Eigen::VectorXd& b, Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim);
void Energy_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void GMRES_Solver2(Eigen::VectorXd&x,Eigen::VectorXd& b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&,double)> matrix_vector_product, Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim);
void GMRES_Solver(Eigen::VectorXd &x, Eigen::VectorXd &b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, double)> Ax_func,Grid &u_con, Grid &Qe, Grid &Qs, double h,const int dim, std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond = {});
double htilde(double hL, double hR, double mL, double mR);
void PreMagSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, double h, const int rks);
void MagSubsystem_UpdateA(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,Grid &A_s,double h, const int rks);
void PreSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks);
void PreMagSubsystem_Energy_potential(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Grid &A,double h, const int rks);
void PreSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks);
void PreMagSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Array &A,double h, const int rks);
void PreMagSubsystem_Pressure_potential(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Grid &A,double h, const int rks);
void MagSubsystem_UpdateB(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,double h, const int rks);
double cal_divergence_error(Grid &u_con, Array &Bx_0p, Array &By_p0);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime);
void bound_cond_staggered(Array &u, int axis);
double getAlfvenSpeed(double *pri){
	double rho = pri[0];
	double cA = pri[5]  / sqrt(rho);
	return cA;
}

void copystate(double *state_1, double *state_0){
	for(int i = 0; i < NUM_VARS; i++){
		state_1[i] = state_0[i];
	}
}

void copy(Grid &state_1, Grid &state_0){
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				state_1(i,j,k) = state_0(i,j,k);
			}
		}
	}
}

void f_convective(double *f, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[8];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	double k = calkinetic(con);
	f[0] = rho * vx;
	f[1] = rho * vx * vx;
	f[2] = rho * vy * vx;
	f[3] = rho * vz * vx;
	f[4] = k * vx;
	f[5] = 0;
	f[6] = 0;
	f[7] = 0;
}

void f_pressure(double *f, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[8];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	f[0] = 0;
	f[1] = p;
	f[2] = 0;
	f[3] = 0;
	f[4] = (U + p) * vx;
	f[5] = 0;
	f[6] = 0;
	f[7] = 0;
}

void f_magnetic(double *f, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[8];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	f[0] = 0;
	f[1] = 0.5 * B_square - Bx * Bx;
	f[2] = -By * Bx;
	f[3] = -Bz * Bx;
	f[4] = B_square * vx - vB * Bx;
	f[5] = 0;
	f[6] = By * vx - vy * Bx;
	f[7] = Bz * vx - vz * Bx;
}


void g_convective(double *g, double *con){
	// calculate f
	//	// calculate g
	//	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	double k = calkinetic(con);
	g[0] = rho * vy;
	g[1] = rho * vx * vy;
	g[2] = rho * vy * vy;
	g[3] = rho * vz * vy;
	g[4] = k * vy;
	g[5] = 0;
	g[6] = 0;
	g[7] = 0;
}

void g_pressure(double *g, double *con){
	// calculate f
	//	// calculate g
	//	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	g[0] = 0;
	g[1] = 0;
	g[2] = p;
	g[3] = 0;
	g[4] = (U + p) * vy;
	g[5] = 0;
	g[6] = 0;
	g[7] = 0;
}

void g_magnetic(double *g, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	
	g[0] = 0;
	g[1] = -Bx * By;
	g[2] = 0.5 * B_square - By * By;
	g[3] = -Bz * By;
	g[4] = B_square * vy - vB * By;
	g[5] = Bx * vy - vx * By;
	g[6] = 0;
	g[7] = Bz * vy - vz * By;
}

void f_full_con(double *f, double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p + 0.5 * B_square - Bx * Bx;
	f[2] = rho * vy * vx - By * Bx;
	f[3] = rho * vz * vx - Bz * Bx;
	f[4] = (kinetic_energy + B_square) * vx - vB * Bx + enthalpy * rho * vx;
	f[5] = 0;
	f[6] = By * vx - vy * Bx;
	f[7] = Bz * vx - vz * Bx;
}

void g_full_con(double *g, double *con){
	// calculate g
	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	g[0] = rho * vy;
	g[1] = rho * vx * vy - Bx * By;
	g[2] = rho * vy * vy + p + 0.5 * B_square - By * By;
	g[3] = rho * vz * vy - Bz * By;
	g[4] = (U + p + 0.5 * B_square) * vy - vB * By;
	g[5] = Bx * vy - vx * By;
	g[6] = 0;
	g[7] = Bz * vy - vz * By;
}


void ini_con2pri(double *pri, double *con){
	// convert conservative (rho rhov E) to primitive (rho v p)
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = con[3] / con[0];
	pri[4] = eos_p_from_rho_e(con[0], e);
	pri[5] = con[5];
	pri[6] = con[6];
	pri[7] = con[7];
}

void ini_pri2con(double *con, double *pri){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	double e = calepri(pri);
	// rho rhov E
	double vx = pri[1];
	double vy = pri[2];
	double vz = pri[3];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * pri[2];
	con[3] = pri[0] * pri[3];
	con[4] = pri[0] * e + 0.5 * pri[0] * v_square + 0.5 * B_square;
	con[5] = pri[5];
	con[6] = pri[6];
	con[7] = pri[7];
}

void upri2ucon(Grid &u_con, Grid &u_pri){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_pri[k] = u_pri(i,j,k);
			}
			ini_pri2con(ini_u_con, ini_u_pri);
			for(int k = 0; k < NUM_VARS; k++){
				u_con(i, j, k) = ini_u_con[k];
			}
		}
	}
}

void ucon2upri(Grid &u_pri,Grid &u_con){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_con[k] = u_con(i, j, k);
			}
			ini_con2pri(ini_u_pri, ini_u_con);
			for(int k = 0; k < NUM_VARS; k++){
				u_pri(i, j, k) = ini_u_pri[k];
			}
		}
	}
}

void MagneticFlux_x(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = (2 - gamma_plasma) * B2 * 0.5 - Q[5] * Q[5];
	flux[2] = -Q[5] * Q[6];
	flux[3] = -Q[5] * Q[7];
	flux[4] = Q[1] / Q[0] * B2 * 0.5 - Q[5] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_y(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = -Q[5] * Q[6];
	flux[2] = (2 - gamma_plasma) * B2 * 0.5 - Q[6] * Q[6];
	flux[3] = -Q[6] * Q[7];
	flux[4] = Q[2] / Q[0] * B2 * 0.5 - Q[6] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_Pressure_x(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = B2 * 0.5 - Q[5] * Q[5];
	flux[2] = -Q[5] * Q[6];
	flux[3] = -Q[5] * Q[7];
	flux[4] = Q[1] / Q[0] * B2 - Q[5] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_Pressure_y(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = -Q[5] * Q[6];
	flux[2] = B2 * 0.5 - Q[6] * Q[6];
	flux[3] = -Q[6] * Q[7];
	flux[4] = Q[2] / Q[0] * B2 - Q[6] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}


void f_KineticEnergyFlux(double *Q, double *flux){
	flux[0] = 0;
	flux[1] = -(gamma_plasma - 1) * 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) / Q[0];
	flux[2] = 0;
	flux[3] = 0;
	flux[4] = 0;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void g_KineticEnergyFlux(double *Q, double *flux){
	flux[0] = 0;
	flux[1] = 0;
	flux[2] = -(gamma_plasma - 1) * 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) / Q[0];
	flux[3] = 0;
	flux[4] = 0;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

double ini_calch_x(double *pri){
	double rho = pri[0];
	double p = pri[4];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B2 = Bx * Bx + By * By + Bz * Bz;
	
	double cf_x = eos_fast_speed(rho, p, Bx, B2);
	
	return fabs(pri[1]) + cf_x;
}

double ini_calch_y(double *pri){
	double rho = pri[0];
	double p = pri[4];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B2 = Bx * Bx + By * By + Bz * Bz;
	
	
	double cf_y = eos_fast_speed(rho, p, By, B2);
	
	return fabs(pri[2]) + cf_y;
}

void ShockDetector(Array &shock_eta, Grid &u_con) {
	// 转换守恒变量到原始变量 (rho, u, v, w, p, Bx, By, Bz)
	Grid u_pri(nxCells, nyCells);
	ucon2upri(u_pri, u_con);
	
	// 中间数组:先存 chi_tilde (公式 4.40),再做邻居扩展 (公式 4.41)
	Array shock_eta_tilde(nxCells, nyCells);
	
	// ---- 第一步:计算 chi_tilde (公式 4.40) ----
	for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			
			// 3x3 邻域内的快磁声速 c_f = sqrt(a^2 + c_A^2),取最小值
			//   a   = sqrt(gamma p / rho)        (声速)
			//   c_A = |B| / sqrt(rho)            (Alfven 速)
			//   u_pri: 0=rho, 4=p, 5=Bx, 6=By, 7=Bz
			double cf_min_nbr = 1.0e300;
			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					double rho_nbr = u_pri(i + di, j + dj, 0);
					double p_nbr   = u_pri(i + di, j + dj, 4);
					double Bx_nbr  = u_pri(i + di, j + dj, 5);
					double By_nbr  = u_pri(i + di, j + dj, 6);
					double Bz_nbr  = u_pri(i + di, j + dj, 7);
					
					double a2  = eos_sound_a2(rho_nbr, p_nbr);
					double cA2 = (Bx_nbr * Bx_nbr + By_nbr * By_nbr + Bz_nbr * Bz_nbr)
					/ rho_nbr;
					double cf  = std::sqrt(a2 + cA2);
					
					if (cf < cf_min_nbr) cf_min_nbr = cf;
				}
			}
			
			// 速度散度 div.u = du/dx + dv/dy (中心差分; w/vz 不进 2D 散度)
			//   u_pri(.,.,1)=u, u_pri(.,.,2)=v
			double div_u = (u_pri(i + 1, j, 1) - u_pri(i - 1, j, 1)) / (2.0 * dx)
			+ (u_pri(i, j + 1, 2) - u_pri(i, j - 1, 2)) / (2.0 * dy);
			
			// 公式 (4.40):
			//   chi_tilde = min[1, max(0, -(h*div.u + k1*cf_min) / (k1*cf_min))]
			double h           = std::min(dx, dy);
			double numerator   = h * div_u + shock_kappa1 * cf_min_nbr;
			double denominator = shock_kappa1 * cf_min_nbr;
			
			shock_eta_tilde(i, j) = std::min(1.0,
				std::max(0.0, -numerator / denominator));
		}
	}
	
	// 边界条件 (让邻域扩展可以访问 ghost cell)
	bound_cond_array(shock_eta_tilde);
	
	// ---- 第二步:邻居扩展 (公式 4.41) ----
	//   chi_{i,j} = max over 3x3 of chi_tilde
	for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			double chi_max = 0.0;
			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					if (shock_eta_tilde(i + di, j + dj) > chi_max) {
						chi_max = shock_eta_tilde(i + di, j + dj);
					}
				}
			}
			shock_eta(i, j) = chi_max;
//			shock_eta(i, j) = 1.0;
		}
	}
	bound_cond_array(shock_eta);
}

void Energy_flux_x(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	f_convective(fCL, QL);
	f_convective(fCR, QR);
	fCL[4] = 0;
	fCR[4] = 0;
	double fKL[NUM_VARS], fKR[NUM_VARS];
	
	f_KineticEnergyFlux(QL, fKL);
	f_KineticEnergyFlux(QR, fKR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_x(QL, fBL);
	MagneticFlux_x(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k] + fKL[k];
		fR[k] = fCR[k] + fBR[k] + fKR[k];
	}
	
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Energy_flux_y(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	g_convective(fCL, QL);
	g_convective(fCR, QR);
	fCL[4] = 0;
	fCR[4] = 0;
	double fKL[NUM_VARS], fKR[NUM_VARS];
	
	g_KineticEnergyFlux(QL, fKL);
	g_KineticEnergyFlux(QR, fKR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_y(QL, fBL);
	MagneticFlux_y(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k] + fKL[k];
		fR[k] = fCR[k] + fBR[k] + fKR[k];
	}
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Energy_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
			}
			Energy_flux_x(Q_xL, Q_xR, fp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			Energy_flux_y(Q_yL, Q_yR, gp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			Energy_flux_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			Energy_flux_y(Q_yL, Q_yR, gm);
			
			for(int k = 0; k < 5; k++){
				u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
			}
		}
	}
}

void Pressure_flux_x(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	f_convective(fCL, QL);
	f_convective(fCR, QR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_Pressure_x(QL, fBL);
	MagneticFlux_Pressure_x(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k];
		fR[k] = fCR[k] + fBR[k];
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
		
//		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Pressure_flux_y(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	g_convective(fCL, QL);
	g_convective(fCR, QR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_Pressure_y(QL, fBL);
	MagneticFlux_Pressure_y(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k];
		fR[k] = fCR[k] + fBR[k];
	}
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
		
//		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}


void Pressure_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
			}
			Pressure_flux_x(Q_xL, Q_xR, fp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			Pressure_flux_y(Q_yL, Q_yR, gp);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			Pressure_flux_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			Pressure_flux_y(Q_yL, Q_yR, gm);
			
			for(int k = 0; k < 5; k++){
				u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
			}
		}
	}
}


void RusanovC_x(double *QL, double *QR, double *flux){
	double fL[NUM_VARS], fR[NUM_VARS];
	
	f_convective(fL, QL);
	f_convective(fR, QR);
	
	switch(formula){
		case EnergyFormula:
			fL[4] = 0;
			fR[4] = 0;
			break;
		case PressureFormula:
			break;
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
	}
}

void RusanovC_y(double *QL, double *QR, double *flux){
	double gL[NUM_VARS], gR[NUM_VARS];
	g_convective(gL, QL);
	g_convective(gR, QR);
	
	switch(formula){
		case EnergyFormula:
			gL[4] = 0;
			gR[4] = 0;
			break;
		case PressureFormula:
			break;
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
	}
}


void RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
				
			}
			RusanovC_x(Q_xL, Q_xR, fp);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			RusanovC_y(Q_yL, Q_yR, gp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			RusanovC_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			RusanovC_y(Q_yL, Q_yR, gm);
			
			switch(formula){
				case EnergyFormula:
					for(int k = 0; k < 5; k++){
						u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
					}
					break;
				case PressureFormula:
					for(int k = 0; k < 5; k++){
						u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
					}
					break;
			}
		}
	}
}

void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR){
	double c;
	double q0, q1, q2,q0L, q1L, q2L,q0R, q1R, q2R, wr0, wr1, wr2, wr0_bar, wr1_bar, wr2_bar, beta0, beta1, beta2, d0, d1, d2;
	double diff1, diff2, diff3, diff4;
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
			for(int k = 0; k < NUM_VARS; k++){
				double slope_xR = u_con(i + 1,j,k) - u_con(i,j,k);
				double slope_xL = u_con(i,j,k) - u_con(i - 1,j,k);
				
				double slope_yR = u_con(i,j + 1,k) - u_con(i,j,k);
				double slope_yL = u_con(i,j,k) - u_con(i,j - 1,k);
				switch (recon) {
					case NoRecon:
						SpaceOrder = 1;
						w_xL(i,j,k) = u_con(i,j,k);
						w_xR(i,j,k) = u_con(i,j,k);
						w_yL(i,j,k) = u_con(i,j,k);
						w_yR(i,j,k) = u_con(i,j,k);
						break;
					case Minmod:
						SpaceOrder = 2;
						if(slope_xL * slope_xR <= 0){
							c = 0;
						}else if(fabs(slope_xR) <= fabs(slope_xL)){
							c = slope_xR;
						}else if(fabs(slope_xL) < fabs(slope_xR)){
							c = slope_xL;
						}
						w_xL(i,j,k) = u_con(i,j,k) - 0.5 * c;
						w_xR(i,j,k) = u_con(i,j,k) + 0.5 * c;
						
						if(slope_yL * slope_yR <= 0){
							c = 0;
						}else if(fabs(slope_yR) <= fabs(slope_yL)){
							c = slope_yR;
						}else if(fabs(slope_yL) < fabs(slope_yR)){
							c = slope_yL;
						}
						w_yL(i,j,k) = u_con(i,j,k) - 0.5 * c;
						w_yR(i,j,k) = u_con(i,j,k) + 0.5 * c;
						break;
					case VanLeer:
						SpaceOrder = 2;
						if(fabs(slope_xR) < 1e-16){
							c = 0;
						}else{
							double r = slope_xL/slope_xR;
							if(r < 0){
								c = 0;
							}else{
								c = fmin(2/(1 + r), (2 * r)/(1 + r));
							}
						}
						w_xL(i,j,k) = u_con(i,j,k) - 0.5 * c * slope_xL;
						w_xR(i,j,k) = u_con(i,j,k) + 0.5 * c * slope_xR;
						
						if(fabs(slope_yR) < 1e-16){
							c = 0;
						}else{
							double r = slope_yL/slope_yR;
							if(r < 0){
								c = 0;
							}else{
								c = fmin(2/(1 + r), (2 * r)/(1 + r));
							}
						}
						w_yL(i,j,k) = u_con(i,j,k) - 0.5 * c * slope_yL;
						w_yR(i,j,k) = u_con(i,j,k) + 0.5 * c * slope_yR;
						break;
					case ENO3:
						SpaceOrder = 3;
						q0L = -1.0/6.0 * u_con(i - 2, j, k) + 5.0/6.0 * u_con(i - 1, j, k) + 1.0/3.0 * u_con(i, j, k);
						q1L = 1.0/3.0 * u_con(i - 1, j, k) + 5.0/6.0 * u_con(i, j, k) - 1.0/6.0 * u_con(i + 1, j, k);
						q2L = 11.0/6.0 * u_con(i, j, k) - 7.0/6.0 * u_con(i + 1, j, k) + 1.0/3.0 * u_con(i + 2, j, k);
						
						
						q0R = 1.0/3.0 * u_con(i - 2, j, k) - 7.0/6.0 * u_con(i - 1, j, k) + 11.0/6.0 * u_con(i, j, k);
						q1R = -1.0/6.0 * u_con(i - 1, j, k) + 5.0/6.0 * u_con(i, j, k) + 1.0/3.0 * u_con(i + 1, j, k);
						q2R = 1.0/3.0 * u_con(i, j, k) + 5.0/6.0 * u_con(i + 1, j, k) - 1.0/6.0 * u_con(i + 2, j, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u_con(i, j, k) - u_con(i - 1, j, k));
						diff2 = fabs(u_con(i + 1, j, k) - u_con(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u_con(i + 2, j, k) - u_con(i + 1, j, k)) - fabs(u_con(i + 1, j, k) - u_con(i, j, k)));
							diff4 = fabs(abs(u_con(i + 1, j, k) - u_con(i, j, k)) - fabs(u_con(i, j, k) - u_con(i - 1, j, k)));
							if (diff3 >= diff4){
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}else{
								w_xL(i, j, k) = q2L;
								w_xR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u_con(i, j, k) - u_con(i - 1, j, k)) - fabs(u_con(i - 1, j, k) - u_con(i - 2, j, k)));
							diff4 = fabs(abs(u_con(i + 1, j, k) - u_con(i, j, k)) - fabs(u_con(i, j, k) - u_con(i - 1, j, k)));
							if (diff4 >= diff3){
								w_xL(i, j, k) = q0L;
								w_xR(i, j, k) = q0R;
							}else{
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}
						}
						
						q0L = -1.0/6.0 * u_con(i, j - 2, k) + 5.0/6.0 * u_con(i, j - 1, k) + 1.0/3.0 * u_con(i, j, k);
						q1L = 1.0/3.0 * u_con(i, j - 1, k) + 5.0/6.0 * u_con(i, j, k) - 1.0/6.0 * u_con(i, j + 1, k);
						q2L = 11.0/6.0 * u_con(i, j, k) - 7.0/6.0 * u_con(i, j + 1, k) + 1.0/3.0 * u_con(i, j + 2, k);
						
						
						q0R = 1.0/3.0 * u_con(i, j - 2, k) - 7.0/6.0 * u_con(i, j - 1, k) + 11.0/6.0 * u_con(i, j, k);
						q1R = -1.0/6.0 * u_con(i, j - 1, k) + 5.0/6.0 * u_con(i, j, k) + 1.0/3.0 * u_con(i, j + 1, k);
						q2R = 1.0/3.0 * u_con(i, j, k) + 5.0/6.0 * u_con(i, j + 1, k) - 1.0/6.0 * u_con(i, j + 2, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u_con(i, j, k) - u_con(i, j - 1, k));
						diff2 = fabs(u_con(i, j + 1, k) - u_con(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u_con(i, j + 2, k) - u_con(i, j + 1, k)) - fabs(u_con(i, j + 1, k) - u_con(i, j, k)));
							diff4 = fabs(abs(u_con(i, j + 1, k) - u_con(i, j, k)) - fabs(u_con(i, j, k) - u_con(i, j - 1, k)));
							if (diff3 >= diff4){
								w_yL(i, j, k) = q1L;
								w_yR(i, j, k) = q1R;
							}else{
								w_yL(i, j, k) = q2L;
								w_yR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u_con(i, j, k) - u_con(i, j - 1, k)) - fabs(u_con(i, j - 1, k) - u_con(i, j - 2, k)));
							diff4 = fabs(abs(u_con(i, j + 1, k) - u_con(i, j, k)) - fabs(u_con(i, j, k) - u_con(i, j - 1, k)));
							if (diff4 >= diff3){
								w_yL(i, j, k) = q0L;
								w_yR(i, j, k) = q0R;
							}else{
								w_yL(i, j, k) = q1L;
								w_yR(i, j, k) = q1R;
							}
						}
						break;
					case WENO5:
						SpaceOrder = 4;
						q0 = -1.0/6.0 * u_con(i - 2, j, k) + 5.0/6.0 * u_con(i - 1, j, k) + 1.0/3.0 * u_con(i, j, k);
						q1 = 1.0/3.0 * u_con(i - 1, j, k) + 5.0/6.0 * u_con(i, j, k) - 1.0/6.0 * u_con(i + 1, j, k);
						q2 = 11.0/6.0 * u_con(i, j, k) - 7.0/6.0 * u_con(i + 1, j, k) + 1.0/3.0 * u_con(i + 2, j, k);
						
						beta0 = 13.0/12.0 * pow((u_con(i - 2, j, k) - 2 * u_con(i - 1, j, k) + u_con(i, j, k)),2.0) + 1.0/4.0 * pow((u_con(i - 2, j, k) - 4 * u_con(i - 1, j, k) + 3 * u_con(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u_con(i - 1, j, k) - 2 * u_con(i, j, k) + u_con(i + 1, j, k)),2.0) + 1.0/4.0 * pow((u_con(i - 1, j, k) - u_con(i + 1, j, k)),2.0);
						beta2 = 13.0/12.0 * pow((u_con(i, j, k) - 2 * u_con(i + 1, j, k) + u_con(i + 2, j, k)),2.0) + 1.0/4.0 * pow((3 * u_con(i, j, k) - 4 * u_con(i + 1, j, k) + u_con(i + 2, j, k)),2.0);
						
						d0 = 3.0/10.0;
						d1 = 6.0/10.0;
						d2 = 1.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_xL(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						q0 = 1.0/3.0 * u_con(i - 2, j, k) - 7.0/6.0 * u_con(i - 1, j, k) + 11.0/6.0 * u_con(i, j, k);
						q1 = -1.0/6.0 * u_con(i - 1, j, k) + 5.0/6.0 * u_con(i, j, k) + 1.0/3.0 * u_con(i + 1, j, k);
						q2 = 1.0/3.0 * u_con(i, j, k) + 5.0/6.0 * u_con(i + 1, j, k) - 1.0/6.0 * u_con(i + 2, j, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_xR(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						
						q0 = -1.0/6.0 * u_con(i, j - 2, k) + 5.0/6.0 * u_con(i, j - 1, k) + 1.0/3.0 * u_con(i, j, k);
						q1 = 1.0/3.0 * u_con(i, j - 1, k) + 5.0/6.0 * u_con(i, j, k) - 1.0/6.0 * u_con(i, j + 1, k);
						q2 = 11.0/6.0 * u_con(i, j, k) - 7.0/6.0 * u_con(i, j + 1, k) + 1.0/3.0 * u_con(i, j + 2, k);
						
						beta0 = 13.0/12.0 * pow((u_con(i, j - 2, k) - 2 * u_con(i, j - 1, k) + u_con(i, j, k)),2.0) + 1.0/4.0 * pow((u_con(i, j - 2, k) - 4 * u_con(i, j - 1, k) + 3 * u_con(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u_con(i, j - 1, k) - 2 * u_con(i, j, k) + u_con(i, j + 1, k)),2.0) + 1.0/4.0 * pow((u_con(i, j - 1, k) - u_con(i, j + 1, k)),2.0);
						beta2 = 13.0/12.0 * pow((u_con(i, j, k) - 2 * u_con(i, j + 1, k) + u_con(i, j + 2, k)),2.0) + 1.0/4.0 * pow((3 * u_con(i, j, k) - 4 * u_con(i, j + 1, k) + u_con(i, j + 2, k)),2.0);
						
						d0 = 3.0/10.0;
						d1 = 6.0/10.0;
						d2 = 1.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_yL(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						q0 = 1.0/3.0 * u_con(i, j - 2, k) - 7.0/6.0 * u_con(i, j - 1, k) + 11.0/6.0 * u_con(i, j, k);
						q1 = -1.0/6.0 * u_con(i, j - 1, k) + 5.0/6.0 * u_con(i, j, k) + 1.0/3.0 * u_con(i, j + 1, k);
						q2 = 1.0/3.0 * u_con(i, j, k) + 5.0/6.0 * u_con(i, j + 1, k) - 1.0/6.0 * u_con(i, j + 2, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_yR(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						break;
				}
			}
		}
	}
	bound_cond_u(w_xL);
	bound_cond_u(w_xR);
	bound_cond_u(w_yL);
	bound_cond_u(w_yR);
}



double calecon(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	double e = (U - 0.5 * rho * v_square - 0.5 * B_square)/rho;
	return e;
}

double calepri(double *pri){
	double e = eos_e_from_rho_p(pri[0], pri[4]);
	return e;
}

double calcf(double* pri){
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double cs = eos_sound_a(pri[0], pri[4]);
	double ca_x = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double ca_y = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double ca_z = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double cf_x = sqrt(0.5 * ((cs * cs + ca_x * ca_x) + sqrt((cs * cs + ca_x * ca_x) * (cs * cs + ca_x * ca_x) - ((4 * (cs * cs * Bx * Bx))/pri[0]))));
	double cf_y = sqrt(0.5 * ((cs * cs + ca_y * ca_y) + sqrt((cs * cs + ca_y * ca_y) * (cs * cs + ca_y * ca_y) - ((4 * (cs * cs * By * By))/pri[0]))));
	double cf_z = sqrt(0.5 * ((cs * cs + ca_z * ca_z) + sqrt((cs * cs + ca_z * ca_z) * (cs * cs + ca_z * ca_z) - ((4 * (cs * cs * Bz * Bz))/pri[0]))));
	double cf = max({cf_x, cf_y, cf_z});
	return cf;
}

vector<double> ConvergenceTestFunction(double x){
	vector<double> f(8);
	f[0] = 1;
	f[1] = 0;
	f[2] = -0.2 * cos(2 * M_PI * x);
	f[3] = 0.2 * sin(2 * M_PI * x);
	f[4] = 1;
	f[5] = 1;
	f[6] = 0.2 * cos(2 * M_PI * x);
	f[7] = -0.2 * sin(2 * M_PI * x);
	return f;
}

vector<double> ConvergenceTestExact(double x){
	vector<double> f(NUM_VARS);
	f[0] = 1;
	f[1] = 0;
	f[2] = -0.2 * cos(2 * M_PI * x);
	f[3] = 0.2 * sin(2 * M_PI * x);
	f[4] = 1;
	f[5] = 1;
	f[6] = 0.2 * cos(2 * M_PI * x);
	f[7] = -0.2 * sin(2 * M_PI * x);
	return f;
}

vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<vector<double>(double)> func, int nGaussPoints){
	vector<double> Weights(nGaussPoints);
	vector<double> Location(nGaussPoints);
	
	// Select Gauss-Legendre nodes & weights
	switch(nGaussPoints){
		case 1:
			Location = {0.0};
			Weights  = {1.0};
			break;
		case 2:
			Location = {-dx/2.0/sqrt(3.0), dx/2.0/sqrt(3.0)};
			Weights  = {0.5, 0.5};
			break;
		case 3:
			Location = {-dx/2.0*sqrt(3.0/5.0), 0.0, dx/2.0*sqrt(3.0/5.0)};
			Weights  = {5.0/18.0, 8.0/18.0, 5.0/18.0};
			break;
		case 4:
			Location = {-dx/2.0*0.8611363115940526,
				-dx/2.0*0.3399810435848563,
				+dx/2.0*0.3399810435848563,
				+dx/2.0*0.8611363115940526};
			
			Weights  = {0.5*0.3478548451374538,
				0.5*0.6521451548625461,
				0.5*0.6521451548625461,
				0.5*0.3478548451374538};
			break;
		case 5:
			Location = {-dx/2.0*0.9061798459386640,
				-dx/2.0*0.5384693101056831,
				0.0,
				+dx/2.0*0.5384693101056831,
				+dx/2.0*0.9061798459386640};
			
			Weights  = {0.5*0.2369268850561891,
				0.5*0.4786286704993665,
				0.5*0.5688888888888889,
				0.5*0.4786286704993665,
				0.5*0.2369268850561891};
			break;
		case 6:
			Location = {-dx/2.0*0.9324695142031521,
				-dx/2.0*0.6612093864662645,
				-dx/2.0*0.2386191860831969,
				+dx/2.0*0.2386191860831969,
				+dx/2.0*0.6612093864662645,
				+dx/2.0*0.9324695142031521};
			
			Weights  = {0.5*0.1713244923791704,
				0.5*0.3607615730481386,
				0.5*0.4679139345726910,
				0.5*0.4679139345726910,
				0.5*0.3607615730481386,
				0.5*0.1713244923791704};
			break;
		default:
			std::cerr << "Number of Gauss-Legendre nodes not implemented\n";
			std::exit(EXIT_FAILURE);
	}
	vector<double> integral;
	integral.resize(func(xc).size(), 0.0);
	
	for(int i = 0; i < nGaussPoints; i++){
		double xGauss = xc + Location[i];
		double w      = Weights[i];
		vector<double> f = func(xGauss);
		for(size_t k = 0; k < f.size(); k++){
			integral[k] += w * f[k];
		}
	}
	
	return integral;
}

void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p){
	vector<double> result;
	double uL[NUM_VARS], uR[NUM_VARS], uIN[NUM_VARS], uOUT[NUM_VARS], Uavg[NUM_VARS];;
	double r;
	double rd;                 // 分界半径：你需要指定
	double delta;     // 平滑厚度：你需要指定
	double s_r;
	double beta_s;
	double R;
	double eta;
	vector<double> Weights(2);
	vector<double> Location(2);
	double xGauss, w;
	const double a = 1.0 / std::sqrt(3.0);
	double xG[2];double yG[2];double xg;double yg;
	double Mx;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			switch (test) {
				case test0:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case test1:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0.75;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0.75;
						u(i, j, 6) = -1;
						u(i, j, 7) = 0;
					}
					break;
				case test1y:
					if(y < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 1;
						u(i, j, 6) = 0.75;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = -1;
						u(i, j, 6) = 0.75;
						u(i, j, 7) = 0;
					}
					break;
				case test2:
					if(x < -0.1){
						u(i,j,0) = 1.08;
						u(i,j,1) = 1.2;
						u(i,j,2) = 0.01;
						u(i,j,3) = 0.5;
						u(i,j,4) = 0.95;
						u(i, j, 5) = 0.564190;
						u(i, j, 6) = 1.015541;
						u(i, j, 7) = 0.564190;
					}else{
						u(i,j,0) = 0.9891;
						u(i,j,1) = -0.0131;
						u(i,j,2) = 0.0269;
						u(i,j,3) = 0.010037;
						u(i,j,4) =  0.97159;
						u(i, j, 5) = 0.564190;
						u(i, j, 6) = 1.135262;
						u(i, j, 7) = 0.564923;
					}
					break;
				case test3:
					if(x < -0.1){
						u(i,j,0) = 1.7;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.7;
						u(i, j, 5) = 1.1;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.2;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = -1.496891;
						u(i,j,4) =  0.2;
						u(i, j, 5) = 1.1;
						u(i, j, 6) = 0.785887;
						u(i, j, 7) = 0.618370;
					}
					break;
				case test4:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 1.3;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.4;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) =  0.4;
						u(i, j, 5) = 1.3;
						u(i, j, 6) = -1;
						u(i, j, 7) = 0;
					}
					break;
				case test5:
					if(x < 0){
						u(i,j,0) = 0.15;
						u(i,j,1) = 21.55;
						u(i,j,2) = 1.0;
						u(i,j,3) = 1.0;
						u(i,j,4) = 0.28;
						u(i, j, 5) = 0.05 / sqrt(4 * M_PI);
						u(i, j, 6) = -2.0 / sqrt(4 * M_PI);
						u(i, j, 7) = -1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 0.1;
						u(i,j,1) = -26.45;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) =  0.1;
						u(i, j, 5) = 0.05 / sqrt(4 * M_PI);
						u(i, j, 6) = 2.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}
					break;
				case test6:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 36.87;
						u(i,j,2) = -0.115;
						u(i,j,3) = -0.0386;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 4.0 / sqrt(4 * M_PI);
						u(i, j, 6) = 4.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = -36.87;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 4 / sqrt(4 * M_PI);
						u(i, j, 6) = 4 / sqrt(4 * M_PI);
						u(i, j, 7) = 1 / sqrt(4 * M_PI);
					}
					break;
				case test7:
					if(x < 0){
						u(i,j,0) = 1.0 / (4 * M_PI);
						u(i,j,1) = -1;
						u(i,j,2) = 1;
						u(i,j,3) = -1;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 1.0 / sqrt(4 * M_PI);
						u(i, j, 6) = -1.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 1 / (4 * M_PI);
						u(i,j,1) = -1;
						u(i,j,2) = -1;
						u(i,j,3) = -1;
						u(i,j,4) = 1;
						u(i, j, 5) = 1 / sqrt(4 * M_PI);
						u(i, j, 6) = 1 / sqrt(4 * M_PI);
						u(i, j, 7) = 1 / sqrt(4 * M_PI);
					}
					break;
				case test8:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = -2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0.5;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = 2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0.5;
						u(i, j, 7) = 0;
					}
					{
						double dlt = 4.0*dx;              // 抹宽约 4 个 cell(δ 越大越平滑)
						double vx  = 1 * sqrt(dx*dx);     // -2 -> +2 平滑过渡,去掉 x=0 速度间断
						u(i,j,0) = 1.0;
						u(i,j,1) = vx;
						u(i,j,2) = 0; u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i,j,5) = 0; u(i,j,6) = 0.5; u(i,j,7) = 0;
						break;
					}
					break;
				case test9:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 10;
						u(i,j,2) = 1;
						u(i,j,3) = -1;
						u(i,j,4) = 20.0;
						u(i, j, 5) = 5.0 / sqrt(4 * M_PI);
						u(i, j, 6) = 5.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -10;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 5 / sqrt(4 * M_PI);
						u(i, j, 6) = 5 / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}
					break;
				case test10:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 0.7;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.3;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 1.0;
						u(i,j,4) = 0.2;
						u(i, j, 5) = 0.7;
						u(i, j, 6) = 1.0;
						u(i, j, 7) = 0;
					}
					break;
				case shocktest:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = -2;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case convergencetest:
					//					u(i,j,0) = 1;
					//					u(i,j,1) = 0;
					//					u(i,j,2) = -0.2 * cos(2 * M_PI * x);
					//					u(i,j,3) = 0.2 * sin(2 * M_PI * x);
					//					u(i,j,4) = 1;
					//					u(i, j, 5) = 1;
					//					u(i, j, 6) = 0.2 * cos(2 * M_PI * x);
					//					u(i, j, 7) = -0.2 * sin(2 * M_PI * x);
					
					result = GaussLegendreNumericalIntegration(x, dx, ConvergenceTestFunction, 2 );
					u(i,j,0) = result[0];
					u(i,j,1) = result[1];
					u(i,j,2) = result[2];
					u(i,j,3) = result[3];
					u(i,j,4) = result[4];
					u(i, j, 5) = result[5];
					u(i, j, 6) = result[6];
					u(i, j, 7) = result[7];
					break;
				case test2d0:
					if(sqrt(x * x + y * y) <= 0.4){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case testOrszag_Tang:
					u(i, j, 0) = gamma_plasma * gamma_plasma;
					u(i, j, 1) = -sin(2 * M_PI * y);
					u(i, j, 2) = sin(2 * M_PI * x);
					u(i, j, 3) = 0;
					u(i, j, 4) = gamma_plasma;
					u(i, j, 5) = -sin(2 * M_PI * y);
					u(i, j, 6) = sin(4 * M_PI * x);
					u(i, j, 7) = 0;
					break;
				case testblastwave: {
					double rad = sqrt(x*x + y*y);
					double rD = 0.1;
					double prsL = 1000.0, prsR = 0.1;
					double eps = 1 * sqrt(dx*dx + dy*dy);   // smoothing width
					
					u(i,j,0) = 1.0;
					u(i,j,1) = 0; u(i,j,2) = 0; u(i,j,3) = 0;
					u(i,j,4) = 0.5*(prsR+prsL) + 0.5*(prsR-prsL) * erf((rad - rD)/eps);
					u(i,j,5) = 70.0 / sqrt(4.0 * M_PI);
					u(i,j,6) = 0.0;
					u(i,j,7) = 0.0;
					break;
				}
				case testblastwave1D:					
					if(fabs(x) <= 0.1){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1000;
						u(i, j, 5) = 1 / sqrt(4 * M_PI) * 70;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 1 / sqrt(4 * M_PI) * 70;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
//					
					break;
				case testRK_Riemann:
					// Redlich-Kwong EOS MHD shock tube(对应 1D 代码 test8 / 论文 4.3.1),x 方向 Riemann
					if(x <= 0){
						u(i,j,0) = 1;
						u(i,j,1) = 1;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 2;
						u(i,j,5) = 0;
						u(i,j,6) = 0;
						u(i,j,7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -1.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 1.0;
						u(i,j,5) = 0;
						u(i,j,6) = 0;
						u(i,j,7) = 0;
					}
					break;
				case testRotor:
					if(sqrt((x) * (x) + (y) * (y)) <= 0.1){
						u(i,j,0) = 10;
						u(i,j,1) = -10 * y;
						u(i,j,2) = 10 * x;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 2.5 / sqrt(4 * M_PI);
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 2.5 / sqrt(4 * M_PI);
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					
					break;
				case testFieldLoop:
					if(sqrt(x * x + y * y) < 0.3){
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 1;
						u(i,j,3) = 0;
						u(i,j,4) = 1e5;
						u(i, j, 5) = (-A0 * y / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
						u(i, j, 6) = (A0 * x / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 1;
						u(i,j,3) = 0;
						u(i,j,4) = 1e5;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case testKHInstability:
					Mx = A0;
					if (y >= -9.0/32.0 && y < -7.0/32.0) {
						eta = 0.5 * (1.0 + sin(16.0 * M_PI * (y + 0.25)));
					}
					else if (y >= -7.0/32.0 && y < 7.0/32.0) {
						eta = 1.0;
					}
					else if (y >= 7.0/32.0 && y < 9.0/32.0) {
						eta = 0.5 * (1.0 - sin(16.0 * M_PI * (y - 0.25)));
					}
					else {
						eta = 0.0;
					}
					u(i,j,0) = 1.4;                         // rho = gamma = 1.4
					u(i,j,1) = Mx * (1.0 - 2.0 * eta);        // u = Mx(1 - 2 eta(y))
					u(i,j,2) = 0.1 * Mx * sin(2.0 * M_PI * x);// v = 0.1 Mx sin(2 pi x)
					u(i,j,3) = 0.0;                           // w = 0
					u(i,j,4) = 1.0;                           // p = 1
					u(i,j,5) = 0.1 * Mx;                      // Bx = 0.1 Mx (uniform)
					u(i,j,6) = 0.0;                           // By = 0
					u(i,j,7) = 0.0;   
					break;
				case testdoubleshear:
					u(i,j,0) = 1;
					if(y < 0.5){
						u(i,j,1) = tanh(30 * (y - 0.25));
					}else{
						u(i,j,1) = tanh(30 * (0.75 - y));
					}
					u(i,j,2) = 0.05 * sin(2 * M_PI * x);
					u(i,j,4) = 1e5 / gamma_plasma;
					break;
				case testPlasmaEquilibrium:
					{
						double B0 = 8.0, R = 0.9;
						double r  = sqrt(x*x + y*y);
						double Bz = (r <= R) ? B0 * (r*r)/(R*R) : B0;     // 纵向场 Bz(r)
						double p  = (r <= R) ? 0.5*(B0*B0 - Bz*Bz) + 0.1  // 力平衡压力
						: 0.1;
						
						u(i,j,0) = 1.0;          // rho
						u(i,j,1) = 0;         // u  (vx)
						u(i,j,2) = 0;         // v  (vy)
						u(i,j,3) = 0.0;          // w
						u(i,j,4) = p;            // p
						u(i,j,5) = 0.0;          // Bx
						u(i,j,6) = 0.0;          // By
						u(i,j,7) = Bz;           // Bz  ← 场全在 z 方向
					}
					break;
				case testAdvectedPlasmaEquilibrium:
					{
						double B0 = 8.0, R = 0.9;
						double r  = sqrt(x*x + y*y);
						double Bz = (r <= R) ? B0 * (r*r)/(R*R) : B0;     // 纵向场 Bz(r)
						double p  = (r <= R) ? 0.5*(B0*B0 - Bz*Bz) + 0.1  // 力平衡压力
						: 0.1;
						
						u(i,j,0) = 1.0;          // rho
						u(i,j,1) = 0.04;         // u  (vx)
						u(i,j,2) = 0.04;         // v  (vy)
						u(i,j,3) = 0.0;          // w
						u(i,j,4) = p;            // p
						u(i,j,5) = 0.0;          // Bx
						u(i,j,6) = 0.0;          // By
						u(i,j,7) = Bz;           // Bz  ← 场全在 z 方向
					}
					break;
				case testMHD_vortex:{
					const double eps      = 1.0;              // 论文里定义的 ε
					const double mu_tilde = sqrt(4.0 * M_PI);  // 论文里定义的 μ̃
					
					double r2 = x*x + y*y;
					
					double expHalf = exp(0.5 * (1.0 - r2));  // e^{0.5(1-r^2)}，用于 u, v, Bx, By
					double expFull = exp(1.0 - r2);          // e^{(1-r^2)}，用于压强磁场修正项
					
					double rho0 = 1.0;
					double u0   = -eps / (2.0*M_PI) * expHalf * y;
					double v0   =  eps / (2.0*M_PI) * expHalf * x;
					double Bx0  = -1.0 / (2.0*M_PI) * expHalf * y;
					double By0  =  1.0 / (2.0*M_PI) * expHalf * x;
					
					double p0 = 1.0
					+ (1.0/(8.0*M_PI)) * (mu_tilde/(2.0*M_PI)) * (mu_tilde/(2.0*M_PI)) * (1.0 - r2) * expFull
					- 0.5 * (eps/(2.0*M_PI)) * (eps/(2.0*M_PI)) * expHalf;
					
					u(i,j,0) = rho0;
					u(i,j,1) = u0;
					u(i,j,2) = v0;
					u(i,j,3) = 0.0;
					u(i,j,4) = p0;
					u(i,j,5) = Bx0;
					u(i,j,6) = By0;
					u(i,j,7) = 0.0;
					break;
				}
				case testHarrisHall: {
					// Table 4.6 Test 12 — symmetric Harris sheet, GEM perturbation
					double th = tanh(y/h12_lam);
					double sech2 = 1.0 - th*th;
					double Bx_bg = h12_B0*th;
					// dAz = psi0 cos(kx x) cos(ky y)
					double dAz_dy = -h12_psi0*h12_ky * cos(h12_kx*x) * sin(h12_ky*y);
					double dAz_dx = -h12_psi0*h12_kx * sin(h12_kx*x) * cos(h12_ky*y);
					u(i,j,0) = sech2 + h12_rho_inf;
					u(i,j,1) = 0.0;
					u(i,j,2) = 0.0;
					u(i,j,3) = 0.0;
					u(i,j,4) = h12_p_inf + 0.5*h12_B0*h12_B0*sech2;
					u(i,j,5) = Bx_bg + dAz_dy;    // Bx =  dAz/dy
					u(i,j,6) = -dAz_dx;           // By = -dAz/dx
					u(i,j,7) = 0.0;
					break;
				}
				case testHallReconnection: {
					// Harris-sheet 背景 Bx(y),By=Bz=0
					double Bx_bg = -1.0 + 2.0*tanh(y/hall_lambda);
					// 扰动 δAz 的解析偏导数, 必须跟下面 A_initial 里的 Az_pert 严格对应
					// δAz = amp*cos(0.5x)*cos(0.5(y-y_NL)),  |δB|max = 0.5*amp
					double kk = 0.5;
					double dAz_dy_pert = -hall_kick_amp*kk * cos(kk*x) * sin(kk*(y - hall_yNL));
					double dAz_dx_pert = -hall_kick_amp*kk * sin(kk*x) * cos(kk*(y - hall_yNL));
					double Bx_tot = Bx_bg + dAz_dy_pert;      // Bx = dAz/dy
					double By_tot = -dAz_dx_pert;             // By = -dAz/dx
					double rho = 1.5 + 0.5*tanh(y/hall_lambda);
					double p   = hall_ptot - 0.5*(Bx_tot*Bx_tot + By_tot*By_tot);  // 总压处处成立(含扰动)
					u(i,j,0) = rho;
					u(i,j,1) = 0.0;
					u(i,j,2) = 0.0;
					u(i,j,3) = 0.0;
					u(i,j,4) = p;
					u(i,j,5) = Bx_tot;
					u(i,j,6) = By_tot;
					u(i,j,7) = 0.0;
					break;
				}
			}
		}
	}
	bound_cond_u(u);
	
	if(test == testHarrisHall){
		// Same Toth (2000) Sec 4.1 line-integral initialisation as the asymmetric case below:
		// faces are exact differences of the corner flux function, so div_h B = 0 at round-off.
		auto Azc = [&](int ii, int jj)->double{
			return h12_Az(x_0 + (ii - ghost_cells + 1.0)*dx,
			              y_0 + (jj - ghost_cells + 1.0)*dy);
		};
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
				Bx_p0(i,j) =  (Azc(i,j) - Azc(i,j-1)) / dy;
				By_0p(i,j) = -(Azc(i,j) - Azc(i-1,j)) / dx;
			}
		}
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				u(i,j,5) = 0.5*(Bx_p0(i-1,j) + Bx_p0(i,j));
				u(i,j,6) = 0.5*(By_0p(i,j-1) + By_0p(i,j));
				// keep total pressure exactly constant for the field actually carried
				u(i,j,4) = h12_p_inf + 0.5*h12_B0*h12_B0
				         - 0.5*(u(i,j,5)*u(i,j,5) + u(i,j,6)*u(i,j,6) + u(i,j,7)*u(i,j,7));
			}
		}
		bound_cond_u(u);
	}else if(test == testHallReconnection){
		// The staggered faces must be a DISCRETE CURL of the corner flux function, not a cell
		// average of the analytic field.  Averaging leaves a truncation-level divergence that
		// constrained transport then conserves forever: it measured 3.9e-3 at 64x32, whereas
		// the reference calculation reports ~1e-11 (Table 7.2).  Taking exact differences of
		// A_z at the corners makes div_h B vanish to round-off at t = 0.
		auto Az_corner = [&](int ii, int jj)->double{
			double xc = x_0 + (ii - ghost_cells + 1.0)*dx;   // x_{i+1/2}
			double yc = y_0 + (jj - ghost_cells + 1.0)*dy;   // y_{j+1/2}
			return -yc + 2.0*hall_lambda*log(cosh(yc/hall_lambda))
			     + hall_kick_amp * cos(0.5*xc) * cos(0.5*(yc - hall_yNL));
		};
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
				Bx_p0(i,j) =  (Az_corner(i,j) - Az_corner(i,j-1)) / dy;
				By_0p(i,j) = -(Az_corner(i,j) - Az_corner(i-1,j)) / dx;
			}
		}
		// Re-derive the collocated field from those faces, and re-impose total-pressure
		// balance on it, so the equilibrium is exact for the field the scheme actually carries.
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				u(i,j,5) = 0.5*(Bx_p0(i-1,j) + Bx_p0(i,j));
				u(i,j,6) = 0.5*(By_0p(i,j-1) + By_0p(i,j));
				u(i,j,4) = hall_ptot - 0.5*(u(i,j,5)*u(i,j,5) + u(i,j,6)*u(i,j,6)
				                          + u(i,j,7)*u(i,j,7));
			}
		}
		bound_cond_u(u);
	}else{
		for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++){
			for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
				Bx_p0(i,j) = 0.5 * (u(i,j,5) + u(i + 1,j,5));
				By_0p(i,j) = 0.5 * (u(i,j,6) + u(i,j + 1,6));
			}
		}
	}
	
	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1); 
}


void bound_cond_u(Grid &u){
	switch (x_bc) {
		case x_Transmissive:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int i = 0; i < ghost_cells; i++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i, j, k) = u(ghost_cells, j, k);
						u(nxCells + ghost_cells + i, j, k) = u(nxCells + ghost_cells - 1, j, k);
					}
				}
			}
			break;
		case x_Reflective:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j,0) = u(2 * ghost_cells - 1 - ghost,j,0);
					u(nxCells + ghost_cells + ghost,j,0) = u(nxCells + ghost_cells - 1 - ghost,j,0);
					
					u(ghost,j,1) = -u(2 * ghost_cells - 1 - ghost,j,1);
					u(nxCells + ghost_cells + ghost,j,1) = -u(nxCells + ghost_cells - 1 - ghost,j,1);
					
					u(ghost,j,2) = -u(2 * ghost_cells - 1 - ghost,j,2);
					u(nxCells + ghost_cells + ghost,j,2) = -u(nxCells + ghost_cells - 1 - ghost,j,2);
					
					u(ghost,j,3) = -u(2 * ghost_cells - 1 - ghost,j,3);
					u(nxCells + ghost_cells + ghost,j,3) = -u(nxCells + ghost_cells - 1 - ghost,j,3);
					
					u(ghost,j,4) = u(2 * ghost_cells - 1 - ghost,j,4);
					u(nxCells + ghost_cells + ghost,j,4) = u(nxCells + ghost_cells - 1 - ghost,j,4);
					
					u(ghost,j,5) = -u(2 * ghost_cells - 1 - ghost,j,5);
					u(nxCells + ghost_cells + ghost,j,5) = -u(nxCells + ghost_cells - 1 - ghost,j,5);
					
					u(ghost,j,6) = -u(2 * ghost_cells - 1 - ghost,j,6);
					u(nxCells + ghost_cells + ghost,j,6) = -u(nxCells + ghost_cells - 1 - ghost,j,6);
					
					u(ghost,j,7) = -u(2 * ghost_cells - 1 - ghost,j,7);
					u(nxCells + ghost_cells + ghost,j,7) = -u(nxCells + ghost_cells - 1 - ghost,j,7);
				}
			}
			break;
		case x_Periodic:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < NUM_VARS; k++){
						u(ghost,j,k) = u(nxCells + ghost,j,k);
						u(nxCells + ghost_cells + ghost,j,k) = u(ghost_cells + ghost,j,k);
					}
				}
			}
			break;
		case x_Dirichlet:{
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int g = 0; g < ghost_cells; g++){
					int L = g, R = nxCells + ghost_cells + g;
					u(L,j,0)=rho_amb; u(L,j,1)=0; u(L,j,2)=0; u(L,j,3)=0;   // rho,v: Dirichlet ambient (v=0)
					u(R,j,0)=rho_amb; u(R,j,1)=0; u(R,j,2)=0; u(R,j,3)=0;
					for(int k = 5; k < NUM_VARS; k++){      // B: transmissive (copy interior) -> match CASTOR
						u(L,j,k) = u(ghost_cells, j, k);
						u(R,j,k) = u(nxCells + ghost_cells - 1, j, k);
					}
					// slot4: 钉压强 p_amb;守恒模式 E = p_amb/(g-1)+1/2 B_ghost^2(E 跟透射 B 浮动,内能恒=ambient)
					if(g_bc_prim){ u(L,j,4)=p_amb; u(R,j,4)=p_amb; }
					else{
						u(L,j,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(L,j,5)*u(L,j,5)+u(L,j,6)*u(L,j,6)+u(L,j,7)*u(L,j,7));
						u(R,j,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(R,j,5)*u(R,j,5)+u(R,j,6)*u(R,j,6)+u(R,j,7)*u(R,j,7));
					}
				}
			}
			break;
		}
	}

	switch(y_bc){
		case y_Transmissive:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int j = 0; j < ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i, j, k) = u(i, ghost_cells, k);
						u(i,nyCells + ghost_cells + j,k) = u(i, nyCells + ghost_cells - 1, k);
					}
				}
			}
			break;
		case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost,0) = u(i,2 * ghost_cells - 1 - ghost,0);
					u(i,nyCells + ghost_cells + ghost,0) = u(i,nyCells + ghost_cells - 1 - ghost,0);
					
					u(i,ghost,1) = -u(i,2 * ghost_cells - 1 - ghost,1);
					u(i,nyCells + ghost_cells + ghost,1) = -u(i,nyCells + ghost_cells - 1 - ghost,1);
					
					u(i,ghost,2) = -u(i,2 * ghost_cells - 1 - ghost,2);
					u(i,nyCells + ghost_cells + ghost,2) = -u(i,nyCells + ghost_cells - 1 - ghost,2);
					
					u(i,ghost,3) = -u(i,2 * ghost_cells - 1 - ghost,3);
					u(i,nyCells + ghost_cells + ghost,3) = -u(i,nyCells + ghost_cells - 1 - ghost,3);
					
					u(i,ghost,4) = u(i,2 * ghost_cells - 1 - ghost,4);
					u(i,nyCells + ghost_cells + ghost,4) = u(i,nyCells + ghost_cells - 1 - ghost,4);
					
					u(i,ghost,5) = -u(i,2 * ghost_cells - 1 - ghost,5);
					u(i,nyCells + ghost_cells + ghost,5) = -u(i,nyCells + ghost_cells - 1 - ghost,5);
					
					u(i,ghost,6) = -u(i,2 * ghost_cells - 1 - ghost,6);
					u(i,nyCells + ghost_cells + ghost,6) = -u(i,nyCells + ghost_cells - 1 - ghost,6);
					
					u(i,ghost,7) = -u(i,2 * ghost_cells - 1 - ghost,7);
					u(i,nyCells + ghost_cells + ghost,7) = -u(i,nyCells + ghost_cells - 1 - ghost,7);
				}
			}
			break;
		case y_Periodic:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i,ghost,k) = u(i,nyCells + ghost,k);
						u(i,nyCells + ghost_cells + ghost,k) = u(i,ghost_cells + ghost,k);
					}
				}
			}
			break;
		case y_Dirichlet:{
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int g = 0; g < ghost_cells; g++){
					int B = g, T = nyCells + ghost_cells + g;
					u(i,B,0)=rho_amb; u(i,B,1)=0; u(i,B,2)=0; u(i,B,3)=0;
					u(i,T,0)=rho_amb; u(i,T,1)=0; u(i,T,2)=0; u(i,T,3)=0;
					for(int k = 5; k < NUM_VARS; k++){      // B: transmissive (copy interior) -> match CASTOR
						u(i,B,k) = u(i, ghost_cells, k);
						u(i,T,k) = u(i, nyCells + ghost_cells - 1, k);
					}
					// slot4: 钉压强 p_amb;守恒模式 E = p_amb/(g-1)+1/2 B_ghost^2
					if(g_bc_prim){ u(i,B,4)=p_amb; u(i,T,4)=p_amb; }
					else{
						u(i,B,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(i,B,5)*u(i,B,5)+u(i,B,6)*u(i,B,6)+u(i,B,7)*u(i,B,7));
						u(i,T,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(i,T,5)*u(i,T,5)+u(i,T,6)*u(i,T,6)+u(i,T,7)*u(i,T,7));
					}
				}
			}
			break;
		}
	}
}

void bound_cond_A(Grid &A) {
	// ===== x方向 =====
	switch (x_bc) {
		case x_Transmissive: case x_Dirichlet: {
			// 线性外推 (论文 Eq.63); 对线性环境 A_z=Bx*y 精确，等价于 Dirichlet
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < 3; k++){
					// 左边界: 斜率从最内两个cell取
					double slope_L = A(ghost_cells + 1, j, k) - A(ghost_cells, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(ghost_cells - 1 - g, j, k) = A(ghost_cells, j, k) - (g + 1) * slope_L;
					}
					// 右边界
					int iR = nxCells + ghost_cells - 1;
					double slope_R = A(iR, j, k) - A(iR - 1, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(iR + 1 + g, j, k) = A(iR, j, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case x_Reflective: {
			// A 在反射边界的行为需要特别处理
			// 对于 Dirichlet 类型的边界, 也用线性外推
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < 3; k++){
					double slope_L = A(ghost_cells + 1, j, k) - A(ghost_cells, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(ghost_cells - 1 - g, j, k) = A(ghost_cells, j, k) - (g + 1) * slope_L;
					}
					int iR = nxCells + ghost_cells - 1;
					double slope_R = A(iR, j, k) - A(iR - 1, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(iR + 1 + g, j, k) = A(iR, j, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case x_Periodic: {
			// 周期BC: 正常复制
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						A(ghost, j, k) = A(nxCells + ghost, j, k);
						A(nxCells + ghost_cells + ghost, j, k) = A(ghost_cells + ghost, j, k);
					}
				}
			}
			break;
		}
	}
	// ===== y方向 =====
	switch (y_bc) {
		case y_Transmissive: case y_Dirichlet: {
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int k = 0; k < 3; k++){
					// 下边界
					double slope_L = A(i, ghost_cells + 1, k) - A(i, ghost_cells, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, ghost_cells - 1 - g, k) = A(i, ghost_cells, k) - (g + 1) * slope_L;
					}
					// 上边界
					int jR = nyCells + ghost_cells - 1;
					double slope_R = A(i, jR, k) - A(i, jR - 1, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, jR + 1 + g, k) = A(i, jR, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case y_Reflective: {
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int k = 0; k < 3; k++){
					double slope_L = A(i, ghost_cells + 1, k) - A(i, ghost_cells, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, ghost_cells - 1 - g, k) = A(i, ghost_cells, k) - (g + 1) * slope_L;
					}
					int jR = nyCells + ghost_cells - 1;
					double slope_R = A(i, jR, k) - A(i, jR - 1, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, jR + 1 + g, k) = A(i, jR, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case y_Periodic: {
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						A(i, ghost, k) = A(i, nyCells + ghost, k);
						A(i, nyCells + ghost_cells + ghost, k) = A(i, ghost_cells + ghost, k);
					}
				}
			}
			break;
		}
	}
}

double caldt(Grid &u_pri){
	vx_max = 0, vy_max = 0, v_max = 0;
	double dt;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			if(fabs(u_pri(i,j,1)) > fabs(vx_max)){
				vx_max = fabs(u_pri(i,j,1));
			}
			if(fabs(u_pri(i,j,2)) > fabs(vy_max)){
				vy_max = fabs(u_pri(i,j,2));
			}
			if(sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3)) > fabs(v_max)){
				v_max = sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			}
		}
	}
	
	//timestep == 1
	switch(caltimestep){
		case Fullydt:
			a_max_x = ch_x;
			a_max_y = ch_y;
			a_max = ch;
			break;
		case Alfvendt:
			if(timestep < full_step + 1){
				a_max_x = ch_x;
				a_max_y = ch_y;
				a_max = ch;
				//std::cout << "ch: "<< ch << endl;
			}else{
				a_max_x = cA_x;
				a_max_y = cA_y;
				a_max = cA;
			}
			break;
		case IMEXdt:
			if(timestep < full_step + 1){
				a_max_x = ch_x;
				a_max_y = ch_y;
				a_max = ch;
				//std::cout << "ch: "<< ch << endl;
			}else{
				a_max_x = vx_max;
				a_max_y = vy_max;
				a_max = v_max;
				//dt = CFL * 1 / (vx_max / dx);
			}
			break;
	}
	dt = CFL * 1 / (2 * a_max / fmin(dx,dy));
	fully_dt = CFL * 1 / (2 * ch / fmin(dx,dy));
	imex_dt = CFL * 1 / (2 * a_max / fmin(dx,dy));
	g_dt_ideal = dt;
	g_dt_hall  = 0.0;
	// ---- Hall accuracy cap -----------------------------------------------------------------
	// With the Hall term implicit this is NOT a stability constraint; it bounds the phase error
	// of the whistler branch, whose grid-scale frequency is omega_w = d_i |B| k^2 / rho.
	// This is the THEORETICAL rho^{-1} form of report Eq. (7.4), not the implemented
	// d_i sqrt(max(B^2/rho)) criterion of Table 7.2, which scales as rho^{-1/2} and therefore
	// under-estimates the frequency in density cavities (report Sec 7.1.2 documents this gap).
	// For reference, the explicit RK4 bound would be |omega dt| <= 2.83.
	if(add_hall && hall_omega_cap > 0.0){
		double k2max = M_PI*M_PI*(1.0/(dx*dx) + 1.0/(dy*dy));
		double om_max = 0.0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double B2 = u_pri(i,j,5)*u_pri(i,j,5) + u_pri(i,j,6)*u_pri(i,j,6)
				          + u_pri(i,j,7)*u_pri(i,j,7);
				double r  = fmax(u_pri(i,j,0), 1e-12);
				double om = Hall_di * sqrt(B2) / r * k2max;
				if(om > om_max) om_max = om;
			}
		}
		if(om_max > 0.0){
			g_dt_hall = hall_omega_cap / om_max;
			dt = fmin(dt, g_dt_hall);
		}
	}
	// ----------------------------------------------------------------------------------------
	static double prev_dt = 0.0;
	if(timestep > 1 && prev_dt > 0.0) dt = fmin(dt, 2.0 * prev_dt);  // CASTOR enlarge_dt = 2.0
	prev_dt = dt;
	if(timestep > full_step){
		a_max_x = vx_max;
		a_max_y = vy_max;
		a_max = v_max;
	}
	return dt;
}

void update(Grid &u_con,Grid &flux, char dir){
	if(dir == 'x'){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				for(int k = 0; k < NUM_VARS; k++){
					u_con(i, j, k) = u_con(i, j, k) - (adt/dx) * (flux(i, j, k) - flux(i - 1, j, k));
				}
			}
		}
	}else if(dir == 'y'){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				for(int k = 0; k < NUM_VARS; k++){
					u_con(i, j, k) = u_con(i, j, k) - (adt/dy) * (flux(i, j, k) - flux(i, j - 1, k));
				}
			}
		}
	}
}

void calch(Grid &u_pri){
	double ch_grid_x = 1e-8, ch_grid_y = 1e-8, ch_grid = 1e-8;
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			double rho = u_pri(i,j,0);
			double p = u_pri(i,j,4);
			double Bx = u_pri(i,j,5);
			double By = u_pri(i,j,6);
			double Bz = u_pri(i,j,7);
			double B2 = Bx * Bx + By * By + Bz * Bz;
			double B = sqrt(B2);
			double v = sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			double cf_x = eos_fast_speed(rho, p, Bx, B2);
			double cf_y = eos_fast_speed(rho, p, By, B2);
			
			double cf = eos_fast_speed(rho, p, B, B2);
			
			if(ch_grid_x < u_pri(i,j,1) + cf_x){
				ch_grid_x = u_pri(i,j,1) + cf_x;
			}
			if(ch_grid_y < u_pri(i,j,2) + cf_y){
				ch_grid_y = u_pri(i,j,2) + cf_y;
			}
			if(ch_grid < v + cf){
				ch_grid = v + cf;
			}
		}
	}
	ch_x = ch_grid_x;
	ch_y = ch_grid_y;
	ch = ch_grid;
}

void calca(Grid &u_pri){
	double ca_grid_x = 1e-8;
	double ca_grid_y = 1e-8;
	double ca_grid = 1e-8;
	double pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			
			for(int k = 0; k < NUM_VARS; k++){
				pri[k] = u_pri(i,j,k);
			}
			double B = sqrt(pri[5] * pri[5] + pri[6] * pri[6] + pri[7] * pri[7]);
			if(ca_grid_x < pri[5]  / sqrt(pri[0])){
				ca_grid_x = pri[5]  / sqrt(pri[0]);
			}
			if(ca_grid_y < pri[6]  / sqrt(pri[0])){
				ca_grid_y = pri[6]  / sqrt(pri[0]);
			}
			if(ca_grid < B / sqrt(pri[0])){
				ca_grid = B  / sqrt(pri[0]);
			}
		}
	}
	cA_x = ca_grid_x;
	cA_y = ca_grid_y;
	cA = ca_grid;
}

void calcp(Grid &u_pri){
	cp_x = sqrt(0.18 * ch_x);
	cp_y = sqrt(0.18 * ch_y);
}

void store(Grid &u_con, Array& vxn,Array& vyn, Array& Bxn, Array& Byn){
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			vxn(i,j) = u_con(i,j,1) / u_con(i,j,0);
			vyn(i,j) = u_con(i,j,2) / u_con(i,j,0);
			Bxn(i,j) = u_con(i,j,5);
			Byn(i,j) = u_con(i,j,6);
		}
	}
}

void bound_cond_array(Array &u){
	switch (x_bc) {
		case x_Transmissive:case x_Dirichlet:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < ghost_cells; k++){
					u(k,j) = u(ghost_cells,j);
					u(nxCells + ghost_cells + k,j) = u(nxCells + ghost_cells - 1,j);
				}
			}
			break;
		case x_Reflective:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j) = -u(2 * ghost_cells - 1 - ghost,j);
					u(nxCells + ghost_cells + ghost,j) = -u(nxCells + ghost_cells - 1 - ghost,j);
				}
			}
			break;
		case x_Periodic:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j) = u(nxCells + ghost,j);
					u(nxCells + ghost_cells + ghost,j) = u(ghost_cells + ghost,j);
				}
			}
			break;
	}
	
	switch(y_bc){
		case y_Transmissive: case y_Dirichlet:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = u(i,ghost_cells);
					u(i,nyCells + ghost_cells + ghost) = u(i,nyCells + ghost_cells - 1);
				}
			}
			break;
		case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = -u(i,2 * ghost_cells - 1 - ghost);
					u(i,nyCells + ghost_cells + ghost) = -u(i,nyCells + ghost_cells - 1 - ghost);
				}
			}
			break;
		case y_Periodic:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = u(i,nyCells + ghost);
					u(i,nyCells + ghost_cells + ghost) = u(i,ghost_cells + ghost);
				}
			}
			break;
	}
}

double calenthalpy(double *con){
	double e = calecon(con);
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double h = e + pri[4] / pri[0];
	return h;
}

double caltotalenthalpy(double *con){
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double h = (con[4] + pri[4]) / pri[0];
	return h;
}

double calkinetic(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double kinetic_energy = 0.5 * rho*(vx * vx + vy * vy + vz * vz);
	//double kinetic_energy = 0.5 * (con[1] * con[1] + con[2] * con[2] + con[3] * con[3]) / (con[0]);
	return kinetic_energy;
}

double calmagenergy(double *con){
	double mag_energy = 0.5 * (con[5] * con[5] + con[6] * con[6] + con[7] * con[7]);
	return mag_energy;
}

void MagneticEigenvalues(Grid &u_pri){
	double bmax_x = 0;
	double bmax_y = 0;
	double b_x,b_y;
	for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			b_x = 0.5 * (fabs(u_pri(i,j,1)) + sqrt(u_pri(i,j,1) * u_pri(i,j,1) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			
			b_y = 0.5 * (fabs(u_pri(i,j,2)) + sqrt(u_pri(i,j,2) * u_pri(i,j,2) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			if(fabs(bmax_x) < fabs(b_x)){
				bmax_x = b_x;
			}
			if(fabs(bmax_y) < fabs(b_y)){
				bmax_y = b_y;
			}
		}
	}
	B_max_x = bmax_x;
	B_max_y = bmax_y;
}

void EnergyEigenvalues(Grid &u_pri){
	double umax_x = 0, umax_y = 0;
	double u_x, u_y;
	for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			double c = eos_sound_a(u_pri(i,j,0), u_pri(i,j,4));
			double m = 0.5 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7));
			double k = 0.5 * (u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			u_x = 0.5 * (u_pri(i,j,1) + sqrt(u_pri(i,j,1) * u_pri(i,j,1) + 4 * (c * c - (gamma_plasma - 1) * (m + k + u_pri(i,j,1) * u_pri(i,j,1)) )));
			u_y = 0.5 * (u_pri(i,j,2) + sqrt(u_pri(i,j,2) * u_pri(i,j,2) + 4 * (c * c - (gamma_plasma - 1) * (m + k + u_pri(i,j,2) * u_pri(i,j,2)) )));
			if(umax_x < u_x){
				umax_x = u_x;
			}
			if(umax_y < u_y){
				umax_y = u_y;
			}
		}
	}
	U_max_x = umax_x;
	U_max_y = umax_y;
}

// CT-consistent BC for a staggered NORMAL-face field.
//   axis=0 -> Bx_p0 (face at i+1/2, normal = x);  axis=1 -> By_0p (face at j+1/2, normal = y).
// 法向(normal):边界法向面 index ghost_cells-1 / N+ghost-1 是 CT 的真实 DOF,必须保留;
//               只把外侧 ghost 面由边界面零阶外推(否则破坏边界 cell 的 div-free)。
// 切向(tangential):普通 transmissive 复制。Periodic 方向:周期复制。
// 通用 bound_cond_array(cell-centered transmissive)会把边界法向面也冲掉 -> 边界 divB 指数增长。
void bound_cond_staggered(Array &u, int axis){
	// ---------- x-direction ----------
	if(x_bc == x_Periodic){
		for(int j = 0; j < nyCells + ghost_cells*2; j++)
			for(int g = 0; g < ghost_cells; g++){
				u(g,j) = u(nxCells+g,j);
				u(nxCells+ghost_cells+g,j) = u(ghost_cells+g,j);
			}
	}else if(axis==0){ // x is NORMAL (Bx_p0): preserve boundary face ghost_cells-1
		for(int j = 0; j < nyCells + ghost_cells*2; j++){
			for(int k = 0; k < ghost_cells-1; k++) u(k,j) = u(ghost_cells-1,j);
			for(int k = 0; k < ghost_cells;   k++) u(nxCells+ghost_cells+k,j) = u(nxCells+ghost_cells-1,j);
		}
	}else{ // x is TANGENTIAL (By_0p): plain transmissive
		for(int j = 0; j < nyCells + ghost_cells*2; j++)
			for(int k = 0; k < ghost_cells; k++){
				u(k,j) = u(ghost_cells,j);
				u(nxCells+ghost_cells+k,j) = u(nxCells+ghost_cells-1,j);
			}
	}
	// ---------- y-direction ----------
	if(y_bc == y_Periodic){
		for(int i = 0; i < nxCells + ghost_cells*2; i++)
			for(int g = 0; g < ghost_cells; g++){
				u(i,g) = u(i,nyCells+g);
				u(i,nyCells+ghost_cells+g) = u(i,ghost_cells+g);
			}
	}else if(axis==1){ // y is NORMAL (By_0p): preserve boundary face ghost_cells-1
		for(int i = 0; i < nxCells + ghost_cells*2; i++){
			for(int k = 0; k < ghost_cells-1; k++) u(i,k) = u(i,ghost_cells-1);
			for(int k = 0; k < ghost_cells;   k++) u(i,nyCells+ghost_cells+k) = u(i,nyCells+ghost_cells-1);
		}
	}else{ // y is TANGENTIAL (Bx_p0): plain transmissive
		for(int i = 0; i < nxCells + ghost_cells*2; i++)
			for(int k = 0; k < ghost_cells; k++){
				u(i,k) = u(i,ghost_cells);
				u(i,nyCells+ghost_cells+k) = u(i,nyCells+ghost_cells-1);
			}
	}
}

// Cell-centred current density from a collocated magnetic field, second-order centred
// differences with mu0 = 1.  This is exactly report Eqs. (3.25)-(3.27).
inline void CellCurrent(Grid& B, int i, int j, double& Jx, double& Jy, double& Jz){
	Jx =  (B(i,j+1,7) - B(i,j-1,7)) / (2.0*dy);
	Jy = -(B(i+1,j,7) - B(i-1,j,7)) / (2.0*dx);
	Jz =  (B(i+1,j,6) - B(i-1,j,6)) / (2.0*dx)
	    - (B(i,j+1,5) - B(i,j-1,5)) / (2.0*dy);
}

void Staggered_CT(Grid& u_con, Array &Bx_p0, Array &By_0p, Array &Bx_p0_Qs, Array &By_0p_Qs,double h){
	Array Ez_0p(nxCells,nyCells), Ez_p0(nxCells,nyCells), Ez_bar_pp(nxCells,nyCells),Ez_pp(nxCells,nyCells),Ez_tilde(nxCells,nyCells);
	double Q_L[NUM_VARS], Q_R[NUM_VARS], fL[NUM_VARS], fR[NUM_VARS], gL[NUM_VARS], gR[NUM_VARS];
	Grid w_xL(nxCells,nyCells), w_xR(nxCells,nyCells),w_yL(nxCells,nyCells), w_yR(nxCells,nyCells);
	
	TVDReconstruction(u_con,w_xL,w_xR,w_yL,w_yR);
	
	for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_xR(i,j,k);
				Q_R[k] = w_xL(i + 1,j,k);
			}
			
			f_full_con(fL,Q_L);
			f_full_con(fR,Q_R);
			
			Ez_p0(i,j) = -(0.5 * (fL[6] + fR[6]) - 0.5 * fmax(fabs(Q_L[1]/Q_L[0]),fabs(Q_R[1]/Q_R[0])) * (Q_R[6] - Q_L[6]));
			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_yR(i,j,k);
				Q_R[k] = w_yL(i,j + 1,k);
			}
			
			g_full_con(gL, Q_L);
			g_full_con(gR, Q_R);
		
			Ez_0p(i,j) = (0.5 * (gL[5] + gR[5]) - 0.5 * fmax(fabs(Q_L[2]/Q_L[0]),fabs(Q_R[2]/Q_R[0])) * (Q_R[5] - Q_L[5]));
		}
	}
	
	bound_cond_array(Ez_p0);
	bound_cond_array(Ez_0p);
	
	for(int j = 0; j < nyCells + 2 * ghost_cells - 2; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells - 2; i++){
			Ez_bar_pp(i,j) = 0.25 * (Ez_p0(i,j) + Ez_p0(i,j + 1) + Ez_0p(i,j) + Ez_0p(i + 1,j));
		}
	}
	
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			Ez_tilde(i,j) = (u_con(i,j,2) / u_con(i,j,0) * u_con(i,j,5) - u_con(i,j,1) / u_con(i,j,0) * u_con(i,j,6));
		}
	}
	bound_cond_array(Ez_bar_pp);
	
	
	for(int j = 1; j < nyCells + 2 * ghost_cells - 2; j++){
		for(int i = 1; i < nxCells + 2 * ghost_cells - 2; i++){
			double term_1 = 0.5 * ((Ez_p0(i,j) - Ez_tilde(i,j)) / (dx / 2.0) + (Ez_p0(i,j + 1) - Ez_tilde(i,j + 1)) / (dx / 2.0));
			double term_2 = 0.5 * ((Ez_tilde(i + 1,j) - Ez_p0(i,j)) / (dx / 2.0) + (Ez_tilde(i + 1,j + 1) - Ez_p0(i,j + 1)) / (dx / 2.0));
			double term_3 = 0.5 * ((Ez_0p(i,j) - Ez_tilde(i,j)) / (dy / 2.0) + (Ez_0p(i + 1,j) - Ez_tilde(i + 1,j)) / (dy / 2.0));
			double term_4 = 0.5 * ((Ez_tilde(i,j + 1) - Ez_0p(i,j)) / (dy / 2.0) + (Ez_tilde(i + 1,j + 1) - Ez_0p(i + 1,j)) / (dy / 2.0));
			
			Ez_pp(i,j) = Ez_bar_pp(i,j) + dx / 8.0 * (term_1 - term_2) + dy / 8.0 * (term_3 - term_4);
		}
	}
	
	// ---- ideal Faraday step --------------------------------------------------------------
	for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
			Bx_p0(i,j) = Bx_p0_Qs(i,j) - h / dy * (Ez_pp(i,j) - Ez_pp(i,j - 1));
			By_0p(i,j) = By_0p_Qs(i,j) + h / dx * (Ez_pp(i,j) - Ez_pp(i - 1,j));
		}
	}
	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1);   // y-normal

	// ===================== non-ideal corner electric field ==================================
	// Second Faraday step for the resistive / Hall / hyper-resistive electric field, applied
	// to the faces the ideal step has just produced.  This mirrors the reference explicit
	// HLLD/CT code, where the non-ideal block runs in post_step on the already-updated face
	// field, and it fixes two things at once relative to the earlier single-pass version:
	//
	//   * the corner current is now built DIRECTLY from the staggered faces,
	//         Jz = (By(i+1,j) - By(i,j))/dx - (Bx(i,j+1) - Bx(i,j))/dy ,
	//     the minimal one-cell stencil natural to constrained transport.  Building it from the
	//     collocated field instead (as before) silently widened the Hall stencil, because each
	//     collocated component is itself a two-face average.
	//   * the faces entering it are post-ideal, i.e. the same time level as u_con, so the two
	//     halves of the electric field are no longer sampled at different levels.
	//
	// Sign convention: E_z^Hall = +(d_i/rho)(J x B)_z, identical to the reference implementation
	// (report Eq. (3.31) prints a minus sign, which is inconsistent with its own Eq. (3.32);
	// the code, here and there, uses the correct sign).
	if((add_hall && hall_in_emf) || add_resistive || add_hyperresistive){
		Array Jz_c(nxCells,nyCells), Ez_non(nxCells,nyCells);
		const int NI = nxCells + 2 * ghost_cells, NJ = nyCells + 2 * ghost_cells;

		for(int i = 0; i < NI - 1; i++){
			for(int j = 0; j < NJ - 1; j++){
				Jz_c(i,j) = (By_0p(i+1,j) - By_0p(i,j)) / dx
				          - (Bx_p0(i,j+1) - Bx_p0(i,j)) / dy;
			}
		}

		for(int i = 1; i < NI - 2; i++){
			for(int j = 1; j < NJ - 2; j++){
				double Enon = 0.0;
				if(add_resistive) Enon += eta_resistive * Jz_c(i,j);

				if(add_hall && hall_in_emf){
					// corner in-plane field: average the two faces that meet at the corner
					double Bx_c = 0.5*(Bx_p0(i,j) + Bx_p0(i,j+1));
					double By_c = 0.5*(By_0p(i,j) + By_0p(i+1,j));
					// Jx = dBz/dy, Jy = -dBz/dx from the four surrounding cell-centred Bz
					double Bz_up    = 0.5*(u_con(i,j+1,7) + u_con(i+1,j+1,7));
					double Bz_down  = 0.5*(u_con(i,j,  7) + u_con(i+1,j,  7));
					double Bz_right = 0.5*(u_con(i+1,j,7) + u_con(i+1,j+1,7));
					double Bz_left  = 0.5*(u_con(i,j,  7) + u_con(i,  j+1,7));
					double Jx_c =  (Bz_up    - Bz_down) / dy;
					double Jy_c = -(Bz_right - Bz_left) / dx;
					double rho_c = 0.25*(u_con(i,j,0) + u_con(i+1,j,0)
					                   + u_con(i,j+1,0) + u_con(i+1,j+1,0));
					rho_c = fmax(rho_c, hall_rho_floor);
					Enon += (Hall_di / rho_c) * (Jx_c*By_c - Jy_c*Bx_c);
				}

				if(add_hyperresistive){
					// E_hyp = -eta_H lap J  reproduces  dB/dt = -eta_H lap^2 B exactly
					double lapJz = (Jz_c(i+1,j) - 2.0*Jz_c(i,j) + Jz_c(i-1,j))/(dx*dx)
					             + (Jz_c(i,j+1) - 2.0*Jz_c(i,j) + Jz_c(i,j-1))/(dy*dy);
					Enon += -eta_hyper * lapJz;
				}
				Ez_non(i,j) = Enon;
			}
		}

		for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
				Bx_p0(i,j) -= h / dy * (Ez_non(i,j) - Ez_non(i,j - 1));
				By_0p(i,j) += h / dx * (Ez_non(i,j) - Ez_non(i - 1,j));
			}
		}
		bound_cond_staggered(Bx_p0, 0);
		bound_cond_staggered(By_0p, 1);
	}
	// ========================================================================================

	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			u_con(i,j,5) = 0.5 * (Bx_p0(i - 1,j) + Bx_p0(i,j));
			u_con(i,j,6) = 0.5 * (By_0p(i,j - 1) + By_0p(i,j));
		}
	}
	bound_cond_u(u_con);
}

void A_initial(Grid &A){
	double C;
	double R;
	double r;
	// 注意: i 和 j 都要覆盖全部网格(包括ghost cells)
	for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
		for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
			A(i,j,0) = 0;
			A(i,j,1) = 0;
			double Ax_val = 0.0, Ay_val = 0.0, Az_val = 0.0;
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			switch(test) {
				case test0: 
				case shocktest:
				case convergencetest: 
				case test2d0:
				case testdoubleshear:
					break;
				// ============================================================
				// 1D x方向 Riemann Problems
				// Az = Bx*y - By*(x - x_d),  Ay = Bz*(x - x_d),  Ax = 0
				// Bx 在两侧相同 (div B = 0 的1D约束)
				// ============================================================
				
				case test1: {  // x_d = 0, Bx=0.75
					double Bx0 = 0.75;
					double By_here = (x < 0) ? 1.0 : -1.0;
					Az_val = Bx0 * y - By_here * (x - 0.0);
					Ay_val = 0.0;  // Bz = 0
					break;
				}
					
				case test2: {  // x_d = -0.1, Bx = 2/sqrt(4pi) ≈ 0.564190
					double xd = -0.1;
					double Bx0 = 0.564190;
					if (x < xd) {
						Az_val = Bx0 * y - 1.015541 * (x - xd);
						Ay_val = 0.564190 * (x - xd);
					} else {
						Az_val = Bx0 * y - 1.135262 * (x - xd);
						Ay_val = 0.564923 * (x - xd);
					}
					break;
				}
					
				case test3: {  // x_d = -0.1, Bx = 1.1
					double xd = -0.1;
					double Bx0 = 1.1;
					if (x < xd) {
						Az_val = Bx0 * y - 1.0 * (x - xd);
						Ay_val = 0.0 * (x - xd);
					} else {
						Az_val = Bx0 * y - 0.785887 * (x - xd);
						Ay_val = 0.618370 * (x - xd);
					}
					break;
				}
					
				case test4: {  // x_d = 0, Bx = 1.3
					double Bx0 = 1.3;
					double By_here = (x < 0) ? 1.0 : -1.0;
					Az_val = Bx0 * y - By_here * (x - 0.0);
					Ay_val = 0.0;  // Bz = 0
					break;
				}
					
				case test5: {  // x_d = 0
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 0.05 / s4p;
					if (x < 0) {
						Az_val = Bx0 * y - (-2.0 / s4p) * x;
						Ay_val = (-1.0 / s4p) * x;
					} else {
						Az_val = Bx0 * y - (2.0 / s4p) * x;
						Ay_val = (1.0 / s4p) * x;
					}
					break;
				}
					
				case test6: {  // x_d = 0, B相同两侧 → A连续
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 4.0 / s4p;
					double By0 = 4.0 / s4p;
					double Bz0 = 1.0 / s4p;
					Az_val = Bx0 * y - By0 * x;
					Ay_val = Bz0 * x;
					break;
				}
					
				case test7: {  // x_d = 0
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 1.0 / s4p;
					double Bz0 = 1.0 / s4p;  // 两侧相同
					if (x < 0) {
						Az_val = Bx0 * y - (-1.0 / s4p) * x;  // By_L = -1/s4p
					} else {
						Az_val = Bx0 * y - (1.0 / s4p) * x;   // By_R = +1/s4p
					}
					Ay_val = Bz0 * x;
					break;
				}
					
					// ============================================================
					// 1D y方向 Riemann Problem
					// Az = -By*x + Bx*(y - y_d),  Ax = -Bz*(y - y_d),  Ay = 0
					// By 在两侧相同 (div B = 0 的1D y约束)
					// ============================================================
					
				case test1y: {  // y_d = 0, By = 0.75
					double By0 = 0.75;
					double Bx_here = (y < 0) ? 1.0 : -1.0;
					Az_val = -By0 * x + Bx_here * (y - 0.0);
					Ay_val = 0.0;
					break;
				}
				case testOrszag_Tang:
					Az_val = cos(2.0*M_PI*y)/(2.0*M_PI) + cos(4.0*M_PI*x)/(4.0*M_PI);
					break;
				
				case testblastwave:
					Az_val = (1.0/sqrt(4.0*M_PI) * 70.0) * y;
					break;
				
				case testRotor:
					Az_val = (2.5/sqrt(4.0*M_PI)) * y;
					break;
				
				case testFieldLoop: {
					double C = A0 / sqrt(4.0*M_PI);
					double R = 0.3;
					double r = sqrt(x*x + y*y);
					Az_val = (r < R) ? C*(R-r) : 0.0;
					break;
				}
					
				case testKHInstability:
					Az_val = 0.1 * A0 * y;
					break;
				
				case testRK_Riemann: {  // x-Riemann,By 跨 x=0 跳变(3.6->4)/s4p;Bx,Bz 两侧相同
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 2.0 / s4p;
					double Bz0 = 2.0 / s4p;
					if (x <= 0) Az_val = Bx0 * y - (3.6 / s4p) * x;
					else        Az_val = Bx0 * y - (4.0  / s4p) * x;
					Ay_val = Bz0 * x;
					break;
				}
				case testPlasmaEquilibrium:
				case testAdvectedPlasmaEquilibrium:
					Az_val = 0.0;     // 面内 Bx=By=0 -> Az=0;Bz 作为独立分量演化,不经 Az
					break;
				case testMHD_vortex: {  // 换成你实际的 enum 名字
					double r2 = x*x + y*y;
					Az_val = (1.0 / (2.0*M_PI)) * exp(0.5 * (1.0 - r2));
					// Bz = 0，所以 Ax = Ay = 0，不用改
					break;
				}
				case testHarrisHall: {
					Az_val = h12_Az(x, y);
					Ay_val = 0.0;
					break;
				}
				case testHallReconnection: {
					// Az_bg(y) 满足 dAz_bg/dy = Bx_bg(y) = -1+2tanh(2y):
					//   ∫(-1+2tanh(2y))dy = -y + ln(cosh(2y))  (常数项不影响 B, 取0)
					double Az_bg   = -y + 2.0*hall_lambda*log(cosh(y/hall_lambda));
					double Az_pert = hall_kick_amp * cos(0.5*x) * cos(0.5*(y - hall_yNL));
					Az_val = Az_bg + Az_pert;
					Ay_val = 0.0;   // Bz≡0, 不需要 Ay
					break;
				}
				default:
					break;
				
			}
			
			A(i, j, 0) = Ax_val;
			A(i, j, 1) = Ay_val;
			A(i, j, 2) = Az_val;
		}
	}
	bound_cond_A(A);
}

inline double minmod(double a, double b) {
	if (a * b <= 0.0) return 0.0;
	return (fabs(a) < fabs(b)) ? a : b;
}

void Unstaggered_CT(Grid& u_con, Grid& A,Grid& A_s, double h) {
	// ===== Step 1: 计算 Godunov 角点电场 (复用 Staggered CT 的逻辑) =====
	Array Ez_0p(nxCells,nyCells), Ez_p0(nxCells,nyCells);
	Array Ez_bar_pp(nxCells,nyCells), Ez_pp(nxCells,nyCells), Ez_tilde(nxCells,nyCells);
	double Q_L[NUM_VARS], Q_R[NUM_VARS], fL[NUM_VARS], fR[NUM_VARS], gL[NUM_VARS], gR[NUM_VARS];
	Grid w_xL(nxCells,nyCells), w_xR(nxCells,nyCells), w_yL(nxCells,nyCells), w_yR(nxCells,nyCells);
	
	TVDReconstruction(u_con, w_xL, w_xR, w_yL, w_yR);
	
	for(int i = 0; i < nxCells + 2*ghost_cells - 1; i++){
		for(int j = 0; j < nyCells + 2*ghost_cells - 1; j++){
			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_xR(i,j,k);
				Q_R[k] = w_xL(i+1,j,k);
			}
			f_full_con(fL, Q_L);
			f_full_con(fR, Q_R);
			
			Ez_p0(i,j) = -(0.5 * (fL[6] + fR[6]) - 0.5 * fmax(fabs(Q_L[1]/Q_L[0]),fabs(Q_R[1]/Q_R[0])) * (Q_R[6] - Q_L[6]));
						
			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_yR(i,j,k);
				Q_R[k] = w_yL(i,j+1,k);
			}
			g_full_con(gL, Q_L);
			g_full_con(gR, Q_R);
			
			Ez_0p(i,j) = (0.5 * (gL[5] + gR[5]) - 0.5 * fmax(fabs(Q_L[2]/Q_L[0]),fabs(Q_R[2]/Q_R[0])) * (Q_R[5] - Q_L[5]));
		}
	}
	bound_cond_array(Ez_p0);
	bound_cond_array(Ez_0p);
	
	// 角点 Ez (Gardiner-Stone 重构)
	for(int j = 0; j < nyCells + 2*ghost_cells - 2; j++){
		for(int i = 0; i < nxCells + 2*ghost_cells - 2; i++){
			Ez_bar_pp(i,j) = 0.25*(Ez_p0(i,j)+Ez_p0(i,j+1)+Ez_0p(i,j)+Ez_0p(i+1,j));
		}
	}
	for(int j = 0; j < nyCells + 2*ghost_cells; j++){
		for(int i = 0; i < nxCells + 2*ghost_cells; i++){
			Ez_tilde(i,j) = u_con(i,j,2)/u_con(i,j,0)*u_con(i,j,5) - u_con(i,j,1)/u_con(i,j,0)*u_con(i,j,6);
		}
	}
	bound_cond_array(Ez_bar_pp);
	
	for(int j = 1; j < nyCells + 2*ghost_cells - 2; j++){
		for(int i = 1; i < nxCells + 2*ghost_cells - 2; i++){
			double t1 = 0.5*((Ez_p0(i,j)-Ez_tilde(i,j))/(dx/2.0) + (Ez_p0(i,j+1)-Ez_tilde(i,j+1))/(dx/2.0));
			double t2 = 0.5*((Ez_tilde(i+1,j)-Ez_p0(i,j))/(dx/2.0) + (Ez_tilde(i+1,j+1)-Ez_p0(i,j+1))/(dx/2.0));
			double t3 = 0.5*((Ez_0p(i,j)-Ez_tilde(i,j))/(dy/2.0) + (Ez_0p(i+1,j)-Ez_tilde(i+1,j))/(dy/2.0));
			double t4 = 0.5*((Ez_tilde(i,j+1)-Ez_0p(i,j))/(dy/2.0) + (Ez_tilde(i+1,j+1)-Ez_0p(i+1,j))/(dy/2.0));
			Ez_pp(i,j) = Ez_bar_pp(i,j) + dx/8.0*(t1-t2) + dy/8.0*(t3-t4);
		}
	}
	
	// ===== Step 2: 用角点 Ez 平均到 cell center 更新 Az =====
	// Ez_pp(i,j) 定义在角点 (i+1/2, j+1/2)
	// cell (i,j) 的四个角: (i-1,j-1), (i,j-1), (i-1,j), (i,j)
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Ez_cell = 0.25*(Ez_pp(i,j) + Ez_pp(i-1,j) + Ez_pp(i,j-1) + Ez_pp(i-1,j-1));
			A(i, j, 2) = A_s(i, j, 2) - h * Ez_cell;
		}
	}
	bound_cond_A(A);
	
	// ===== Step 3: B = curl A =====
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i, j, 5) =  (A(i,j+1,2) - A(i,j-1,2)) / (2.0*dy);
			u_con(i, j, 6) = -(A(i+1,j,2) - A(i-1,j,2)) / (2.0*dx);
		}
	}
	
	bound_cond_u(u_con);
}

int M_P(int i, int j){
	int coord = 0;
	int xL, xR, yL, yR;
	int I, J;
	
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = ghost_cells;
				break;
			case x_Reflective:
				I = 2 * ghost_cells - 1 - i;
				break;
			case x_Periodic:
				I = nxCells + i;
				break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = nxCells + ghost_cells - 1;
				break;
			case x_Reflective:
				I = 2 * (nxCells + ghost_cells) - i - 1;
				break;
			case x_Periodic:
				I = i - nxCells;
				break;
		}
	}else{
		I = i;
	}
	
	if(j < ghost_cells){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet:
				J = ghost_cells;
				break;
			case y_Reflective:
				J = 2 * ghost_cells - 1 - j;
				break;
			case y_Periodic:
				J = nyCells + j;
				break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet:
				J = nyCells + ghost_cells - 1;
				break;
			case y_Reflective:
				J = 2 * (nyCells + ghost_cells) - j - 1;
				break;
			case y_Periodic:
				J = j - nyCells;
				break;
		}
	}else{
		J = j;
	}
	coord = (I - ghost_cells) + nxCells * (J - ghost_cells);
	return coord;
}

int M_B(int i, int j, int k){
	int coord = 0;
	int xL, xR, yL, yR;
	int I, J;
	
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = ghost_cells;
				break;
			case x_Reflective:
				I = 2 * ghost_cells - 1 - i;
				break;
			case x_Periodic:
				I = nxCells + i;
				break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = nxCells + ghost_cells - 1;
				break;
			case x_Reflective:
				I = 2 * (nxCells + ghost_cells) - i - 1;
				break;
			case x_Periodic:
				I = i - nxCells;
				break;
		}
	}else{
		I = i;
	}
	
	if(j < ghost_cells){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet:
				J = ghost_cells;
				break;
			case y_Reflective:
				J = 2 * ghost_cells - 1 - j;
				break;
			case y_Periodic:
				J = nyCells + j;
				break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet:
				J = nyCells + ghost_cells - 1;
				break;
			case y_Reflective:
				J = 2 * (nyCells + ghost_cells) - j - 1;
				break;
			case y_Periodic:
				J = j - nyCells;
				break;
		}
	}else{
		J = j;
	}
	
	if(k == 5){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells);
	}else if(k == 6){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells) + nxCells * nyCells;
	}else if(k == 7){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells) + 2 * nxCells * nyCells;
	}else{
		std::cout << "Error" << endl;
	}
	
	return coord;
}

void Magnetic_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h){
	int reflct_xL, reflct_yL, reflct_xR, reflct_yR;
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			reflct_xL = 1, reflct_yL = 1,reflct_xR = 1,reflct_yR = 1;
			int xL, xR, yL, yR;
			switch(x_bc){
				case x_Transmissive: case x_Dirichlet:
					xL = max(ghost_cells,i-1);
					xR = min(nxCells + ghost_cells - 1,i+1);
					break;
				case x_Reflective:
					xL = max(ghost_cells,i-1);
					xR = min(nxCells + ghost_cells - 1,i+1);
					
					if(xL == ghost_cells){
						reflct_xL = -1;
					}
					
					if(xR == nxCells + ghost_cells - 1){
						reflct_xR = -1;
					}
					break;
				case x_Periodic:
					if(i == ghost_cells){
						xL = nxCells + ghost_cells - 1;
					}else{
						xL = i - 1;
					}
					
					if(i == nxCells + ghost_cells - 1){
						xR = ghost_cells;
					}else{
						xR = i + 1;
					}
					break;
			}
			
			switch(y_bc){
				case y_Transmissive: case y_Dirichlet:
					yL = max(ghost_cells,j-1);
					yR = min(nyCells + ghost_cells - 1,j+1);
					break;
				case y_Reflective:
					yL = max(ghost_cells,j-1);
					yR = min(nyCells + ghost_cells - 1,j+1);
					
					if(yL == ghost_cells){
						reflct_yL = -1;
					}
					
					if(yR == nyCells + ghost_cells - 1){
						reflct_yR = -1;
					}
					break;
				case y_Periodic:
					if(j == ghost_cells){
						yL = nyCells + ghost_cells - 1;
					}else{
						yL = j - 1;
					}
					
					if(j == nyCells + ghost_cells - 1){
						yR = ghost_cells;
					}else{
						yR = j + 1;
					}
					break;
			}
			
			double rho_p0 = 0.5 * (u_con(i,j,0) + u_con(i + 1,j,0));
			double rho_m0 = 0.5 * (u_con(i,j,0) + u_con(i - 1,j,0));
			double rho_0p = 0.5 * (u_con(i,j,0) + u_con(i,j + 1,0));
			double rho_0m = 0.5 * (u_con(i,j,0) + u_con(i,j - 1,0));
			
			double Bxn_p0 = 0.5 * (Qe(i,j,5) + Qe(i + 1,j,5));
			double Bxn_m0 = 0.5 * (Qe(i,j,5) + Qe(i - 1,j,5));
			double Bxn_0p = 0.5 * (Qe(i,j,5) + Qe(i,j + 1,5));
			double Bxn_0m = 0.5 * (Qe(i,j,5) + Qe(i,j - 1,5));
			double Bxn_pp = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j + 1,5) + Qe(i + 1,j + 1,5));
			double Bxn_pm = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j - 1,5) + Qe(i + 1,j - 1,5));
			double Bxn_mp = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j + 1,5) + Qe(i - 1,j + 1,5));
			double Bxn_mm = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j - 1,5) + Qe(i - 1,j - 1,5));
			
			double Byn_p0 = 0.5 * (Qe(i,j,6) + Qe(i + 1,j,6));
			double Byn_m0 = 0.5 * (Qe(i,j,6) + Qe(i - 1,j,6));
			double Byn_0p = 0.5 * (Qe(i,j,6) + Qe(i,j + 1,6));
			double Byn_0m = 0.5 * (Qe(i,j,6) + Qe(i,j - 1,6));
			double Byn_pp = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j + 1,6) + Qe(i + 1,j + 1,6));
			double Byn_pm = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j - 1,6) + Qe(i + 1,j - 1,6));
			double Byn_mp = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j + 1,6) + Qe(i - 1,j + 1,6));
			double Byn_mm = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j - 1,6) + Qe(i - 1,j - 1,6));
			
			double Bzn_p0 = 0.5 * (Qe(i,j,7) + Qe(i + 1,j,7));
			double Bzn_m0 = 0.5 * (Qe(i,j,7) + Qe(i - 1,j,7));
			double Bzn_0p = 0.5 * (Qe(i,j,7) + Qe(i,j + 1,7));
			double Bzn_0m = 0.5 * (Qe(i,j,7) + Qe(i,j - 1,7));
			double Bzn_pp = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j + 1,7) + Qe(i + 1,j + 1,7));
			double Bzn_pm = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j - 1,7) + Qe(i + 1,j - 1,7));
			double Bzn_mp = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j + 1,7) + Qe(i - 1,j + 1,7));
			double Bzn_mm = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j - 1,7) + Qe(i - 1,j - 1,7));
			
			double Bx_p0 = 0.5 * (x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR);
			double Bx_m0 = 0.5 * (x(M_B(i,j,5)) + x(M_B(i-1,j,5)) * reflct_xL);
			double Bx_0p = 0.5 * (x(M_B(i,j,5)) + x(M_B(i,j+1,5)) * reflct_yR);
			double Bx_0m = 0.5 * (x(M_B(i,j,5)) + x(M_B(i,j-1,5)) * reflct_yL);
			double Bx_pp = 0.25 * (x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR + x(M_B(i,j+1,5)) * reflct_yR + x(M_B(i+1,j+1,5)) * reflct_xR * reflct_yR);
			double Bx_pm = 0.25 * (x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR + x(M_B(i,j-1,5)) * reflct_yL + x(M_B(i+1,j-1,5)) * reflct_xR * reflct_yL);
			double Bx_mp = 0.25 * (x(M_B(i,j,5)) + x(M_B(i-1,j,5)) * reflct_xL + x(M_B(i,j+1,5)) * reflct_yR + x(M_B(i-1,j+1,5)) * reflct_xL * reflct_yR);
			double Bx_mm = 0.25 * (x(M_B(i,j,5)) + x(M_B(i-1,j,5)) * reflct_xL + x(M_B(i,j-1,5)) * reflct_yL + x(M_B(i-1,j-1,5)) * reflct_xL * reflct_yL);
			
			double By_p0 = 0.5 * (x(M_B(i,j,6)) + x(M_B(i+1,j,6)) * reflct_xR);
			double By_m0 = 0.5 * (x(M_B(i,j,6)) + x(M_B(i-1,j,6)) * reflct_xL);
			double By_0p = 0.5 * (x(M_B(i,j,6)) + x(M_B(i,j+1,6)) * reflct_yR);
			double By_0m = 0.5 * (x(M_B(i,j,6)) + x(M_B(i,j-1,6)) * reflct_yL);
			double By_pp = 0.25 * (x(M_B(i,j,6)) + x(M_B(i+1,j,6)) * reflct_xR + x(M_B(i,j+1,6)) * reflct_yR + x(M_B(i+1,j+1,6)) * reflct_xR * reflct_yR);
			double By_pm = 0.25 * (x(M_B(i,j,6)) + x(M_B(i+1,j,6)) * reflct_xR + x(M_B(i,j-1,6)) * reflct_yL + x(M_B(i+1,j-1,6)) * reflct_xR * reflct_yL);
			double By_mp = 0.25 * (x(M_B(i,j,6)) + x(M_B(i-1,j,6)) * reflct_xL + x(M_B(i,j+1,6)) * reflct_yR + x(M_B(i-1,j+1,6)) * reflct_xL * reflct_yR);
			double By_mm = 0.25 * (x(M_B(i,j,6)) + x(M_B(i-1,j,6)) * reflct_xL + x(M_B(i,j-1,6)) * reflct_yL + x(M_B(i-1,j-1,6)) * reflct_xL * reflct_yL);
			
			double Bz_p0 = 0.5 * (x(M_B(i,j,7)) + x(M_B(i+1,j,7)) * reflct_xR);
			double Bz_m0 = 0.5 * (x(M_B(i,j,7)) + x(M_B(i-1,j,7)) * reflct_xL);
			double Bz_0p = 0.5 * (x(M_B(i,j,7)) + x(M_B(i,j+1,7)) * reflct_yR);
			double Bz_0m = 0.5 * (x(M_B(i,j,7)) + x(M_B(i,j-1,7)) * reflct_yL);
			double Bz_pp = 0.25 * (x(M_B(i,j,7)) + x(M_B(i+1,j,7)) * reflct_xR + x(M_B(i,j+1,7)) * reflct_yR + x(M_B(i+1,j+1,7)) * reflct_xR * reflct_yR);
			double Bz_pm = 0.25 * (x(M_B(i,j,7)) + x(M_B(i+1,j,7)) * reflct_xR + x(M_B(i,j-1,7)) * reflct_yL + x(M_B(i+1,j-1,7)) * reflct_xR * reflct_yL);
			double Bz_mp = 0.25 * (x(M_B(i,j,7)) + x(M_B(i-1,j,7)) * reflct_xL + x(M_B(i,j+1,7)) * reflct_yR + x(M_B(i-1,j+1,7)) * reflct_xL * reflct_yR);
			double Bz_mm = 0.25 * (x(M_B(i,j,7)) + x(M_B(i-1,j,7)) * reflct_xL + x(M_B(i,j-1,7)) * reflct_yL + x(M_B(i-1,j-1,7)) * reflct_xL * reflct_yL);
			
			//Laplacian Operator
//			double nu_x = B_max_x * dx;
//			double nu_y = B_max_y * dy;
			
			
			double Bx_1 = Bxn_0p / rho_0p * (1 / dx * (Bxn_pp * By_pp - Bxn_mp * By_mp) - 1 / dy * ((0.5 * (Qe(i, j + 1, 5) * x(M_B(i,j + 1,5)) * reflct_yR + Qe(i, j + 1, 6) * x(M_B(i,j + 1,6)) * reflct_yR + Qe(i, j + 1, 7) * x(M_B(i,j + 1,7)) * reflct_yR ) - Qe(i, j + 1, 6) * x(M_B(i,j + 1,6)) * reflct_yR) - (0.5 * (Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i,j,7))) - Qe(i, j, 6) * x(M_B(i,j,6))) ) );
			
			double Bx_2 = Byn_0p / rho_0p * (1 / dx * ((0.5 * (Bxn_pp * Bx_pp + Byn_pp * By_pp + Bzn_pp * Bz_pp) - Bxn_pp * Bx_pp) - (0.5 * (Bxn_mp * Bx_mp + Byn_mp * By_mp + Bzn_mp * Bz_mp) - Bxn_mp * Bx_mp )) - 1 / dy * (Qe(i,j + 1,6) * x(M_B(i,j + 1,5)) * reflct_yR - Qe(i,j,6) * x(M_B(i,j,5)))  );
			
			double Bx_3 = Bxn_0m / rho_0m * (1 / dx * (Bxn_pm * By_pm - Bxn_mm * By_mm) - 1 / dy * ((0.5 *(Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i,j,7))) - Qe(i, j, 6) * x(M_B(i, j, 6))) - (0.5 * (Qe(i, j - 1, 5) * x(M_B(i,j - 1,5)) * reflct_yL + Qe(i, j - 1, 6) * x(M_B(i,j - 1,6))  * reflct_yL + Qe(i, j - 1, 7) * x(M_B(i,j - 1,7))  * reflct_yL ) - Qe(i, j - 1, 6) * x(M_B(i,j - 1,6))  * reflct_yL) ) );
			
			double Bx_4 =  Byn_0m / rho_0m * (1 / dx * ((0.5 * (Bxn_pm * Bx_pm + Byn_pm * By_pm + Bzn_pm * Bz_pm) - Bxn_pm * Bx_pm) - (0.5 * (Bxn_mm * Bx_mm + Byn_mm * By_mm + Bzn_mm * Bz_mm) - Bxn_mm * Bx_mm)) - 1 / dy * (Qe(i,j,6) * x(M_B(i,j,5)) - Qe(i,j - 1,6) * x(M_B(i,j - 1,5)) * reflct_yL)  );
			
			
			Ax(M_B(i,j,5)) = x(M_B(i,j,5)) + (h * h/ dy) * ((Bx_1 + Bx_2) - (Bx_3 + Bx_4));

			double By_1 = Byn_p0 / rho_p0 * (-1 / dx * ((0.5 * ((Qe(i + 1, j, 5) * x(M_B(i + 1,j,5)) * reflct_xR + Qe(i + 1, j, 6) * x(M_B(i+1,j,6)) * reflct_xR + Qe(i + 1, j, 7) * x(M_B(i + 1, j, 7)) * reflct_xR )) - Qe(i + 1, j, 5) * x(M_B(i + 1, j, 5)) * reflct_xR ) - (0.5 * ((Qe(i, j, 5) * x(M_B(i, j, 5)) + Qe(i, j, 6) * x(M_B(i, j, 6)) + Qe(i, j, 7) * x(M_B(i, j, 7)))) - Qe(i, j, 5) * x(M_B(i, j, 5)) )) + 1 / dy * (Bxn_pp * By_pp - Bxn_pm * By_pm));
			
			double By_2 = Bxn_p0 / rho_p0 * (-1 / dx * (Qe(i + 1, j, 5) * x(M_B(i + 1, j, 6)) * reflct_xR - Qe(i, j, 5) * x(M_B(i, j, 6))) + 1 / dy * ((0.5 * (Bxn_pp * Bx_pp + Byn_pp * By_pp + Bzn_pp * Bz_pp) - Byn_pp * By_pp) - (0.5 * (Bxn_pm * Bx_pm + Byn_pm * By_pm + Bzn_pm * Bz_pm ) - Byn_pm * By_pm )));
			
			double By_3 = Byn_m0 / rho_m0 * (-1 / dx * ((0.5 * (Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i, j, 7)) ) - Qe(i, j, 5) * x(M_B(i, j, 5))) - (0.5 * ((Qe(i - 1, j, 5) * x(M_B(i - 1, j, 5)) * reflct_xL + Qe(i - 1, j, 6) * x(M_B(i - 1, j, 6)) * reflct_xL + Qe(i - 1, j, 7) * x(M_B(i - 1, j, 7)) * reflct_xL)) - Qe(i - 1, j, 5) * x(M_B(i - 1, j, 5)) * reflct_xL )) + 1 / dy * (Bxn_mp * By_mp - Bxn_mm * By_mm));
			
			double By_4 = Bxn_m0 / rho_m0 * (-1 / dx * (Qe(i, j, 5) * x(M_B(i, j, 6)) - Qe(i - 1, j, 5) * x(M_B(i - 1, j, 6)) * reflct_xL) + 1 / dy * ((0.5 * (Bxn_mp * Bx_mp + Byn_mp * By_mp + Bzn_mp * Bz_mp) - Byn_mp * By_mp)  - (0.5 * (Bxn_mm * Bx_mm + Byn_mm * By_mm + Bzn_mm * Bz_mm ) - Byn_mm * By_mm )));
			
			Ax(M_B(i,j,6)) = x(M_B(i,j,6)) + (h * h / dx) * ((By_1 + By_2) - (By_3 + By_4));
			

			double Bz_1 = Bzn_p0 / rho_p0 * (-1/dx * ((0.5 * (Qe(i + 1,j,5) * x(M_B(i + 1,j,5)) * reflct_xR + Qe(i + 1,j,6) * x(M_B(i + 1,j,6))* reflct_xR + Qe(i + 1,j,7) * x(M_B(i + 1,j,7))* reflct_xR) - Qe(i + 1,j,5) * x(M_B(i + 1,j,5))* reflct_xR ) - (0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,5) * x(M_B(i,j,5)))) + 1 / dy * (Byn_pp * Bx_pp - Byn_pm * Bx_pm ) ); 
			
			double Bz_2 =  Bxn_p0 / rho_p0 * (-1 / dx * (Qe(i + 1,j,5) * x(M_B(i + 1, j, 7)) * reflct_xR - Qe(i,j,5) * x(M_B(i, j, 7))) - 1 / dy * (Byn_pp * Bz_pp - Byn_pm * Bz_pm));
			
			double Bz_3 = Bzn_m0 / rho_m0 * (-1/dx * ((0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,5) * x(M_B(i,j,5))) - (0.5 * (Qe(i - 1,j,5) * x(M_B(i - 1,j,5)) * reflct_xL + Qe(i - 1,j,6) * x(M_B(i - 1,j,6)) * reflct_xL + Qe(i - 1,j,7) * x(M_B(i - 1,j,7))* reflct_xL) - Qe(i-1,j,5) * x(M_B(i-1,j,5)) * reflct_xL )) + 1 / dy * (Byn_mp * Bx_mp - Byn_mm * Bx_mm ) ); 
			
			double Bz_4 =  Bxn_m0 / rho_m0 * (-1 / dx * (Qe(i,j,5) * x(M_B(i, j, 7)) - Qe(i - 1,j,5) * x(M_B(i - 1, j, 7)) * reflct_xL ) - 1 / dy * (Byn_mp * Bz_mp - Byn_mm * Bz_mm));
			
			double Bz_5 = Bzn_0p / rho_0p * (1 / dx * (Bxn_pp * By_pp - Bxn_mp * By_mp) - 1 / dy * ((0.5 * (Qe(i,j+1,5) * x(M_B(i,j+1,5)) * reflct_yR + Qe(i,j+1,6) * x(M_B(i,j+1,6))* reflct_yR + Qe(i,j+1,7) * x(M_B(i,j+1,7)) * reflct_yR ) - Qe(i,j+1,6) * x(M_B(i,j+1,6))* reflct_yR)  - (0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7)) ) - Qe(i,j,6) * x(M_B(i,j,6)))) );
			
			double Bz_6 = Byn_0p / rho_0p * (-1 / dx * (Bxn_pp * Bz_pp - Bxn_mp * Bz_mp) - 1 / dy * (Qe(i,j + 1,6) * x(M_B(i, j + 1, 7)) * reflct_yR - Qe(i,j,6) * x(M_B(i, j, 7))) );
			
			double Bz_7 = Bzn_0m / rho_0m * (1 / dx * (Bxn_pm * By_pm - Bxn_mm * By_mm) - 1 / dy * ((0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,6) * x(M_B(i,j,6))) - (0.5 * (Qe(i,j - 1,5) * x(M_B(i,j - 1,5)) * reflct_yL + Qe(i,j - 1,6) * x(M_B(i,j - 1,6)) * reflct_yL + Qe(i,j - 1,7) * x(M_B(i,j - 1,7)) * reflct_yL ) - Qe(i,j - 1,6) * x(M_B(i,j - 1,6))* reflct_yL)) );
			
			double Bz_8 = Byn_0m / rho_0m * (-1 / dx * (Bxn_pm * Bz_pm - Bxn_mm * Bz_mm) - 1 / dy * (Qe(i,j,6) * x(M_B(i, j, 7)) - Qe(i,j - 1,6) * x(M_B(i, j - 1, 7))* reflct_yL) );
			

			
			Ax(M_B(i,j,7)) = x(M_B(i,j,7)) + h *h  / dx * ((Bz_1 + Bz_2) - (Bz_3 + Bz_4)) + h * h/ dy * ((Bz_5 + Bz_6) - (Bz_7 + Bz_8));
			
			// ---- implicit Hall term (Dematte Eq. 59 pattern, extended to the Hall EMF) ----
			// The Hall term is quadratic in B, so it is linearised the same way the ideal
			// tension already is in this operator: the current is taken from the implicit
			// unknown x, the multiplying field is frozen at the explicit stage state Qe.
			//   B^{n+1} + h * curl[ (d_i/rho) (curl x) x B^n ]  =  B*
			if(add_hall && hall_implicit){
				auto xB = [&](int ii,int jj,int c)->double{ return x(M_B(ii,jj,c)); };
				auto Jx_x = [&](int ii,int jj)->double{
					return  (xB(ii,jj+1,7) - xB(ii,jj-1,7)) / (2.0*dy); };
				auto Jy_x = [&](int ii,int jj)->double{
					return -(xB(ii+1,jj,7) - xB(ii-1,jj,7)) / (2.0*dx); };
				auto Jz_x = [&](int ii,int jj)->double{
					return  (xB(ii+1,jj,6) - xB(ii-1,jj,6)) / (2.0*dx)
					      - (xB(ii,jj+1,5) - xB(ii,jj-1,5)) / (2.0*dy); };
				auto EHx = [&](int ii,int jj)->double{
					double r = fmax(u_con(ii,jj,0), 1e-12);
					return (Hall_di/r)*(Jy_x(ii,jj)*Qe(ii,jj,7) - Jz_x(ii,jj)*Qe(ii,jj,6)); };
				auto EHy = [&](int ii,int jj)->double{
					double r = fmax(u_con(ii,jj,0), 1e-12);
					return (Hall_di/r)*(Jz_x(ii,jj)*Qe(ii,jj,5) - Jx_x(ii,jj)*Qe(ii,jj,7)); };
				auto EHz = [&](int ii,int jj)->double{
					double r = fmax(u_con(ii,jj,0), 1e-12);
					return (Hall_di/r)*(Jx_x(ii,jj)*Qe(ii,jj,6) - Jy_x(ii,jj)*Qe(ii,jj,5)); };
				// dB/dt = -curl E  =>  LHS gains +h*curl E_H
				Ax(M_B(i,j,5)) += h * ( (EHz(i,j+1) - EHz(i,j-1)) / (2.0*dy) );
				Ax(M_B(i,j,6)) += h * (-(EHz(i+1,j) - EHz(i-1,j)) / (2.0*dx) );
				Ax(M_B(i,j,7)) += h * ( (EHy(i+1,j) - EHy(i-1,j)) / (2.0*dx)
				                      - (EHx(i,j+1) - EHx(i,j-1)) / (2.0*dy) );
			}
			
			double lap_Bx, lap_By, lap_Bz;
			double nu_x, nu_y;
			
			if(add_resistive){
				double lapBx_res = (x(M_B(i+1,j,5))*reflct_xR - 2.0*x(M_B(i,j,5)) + x(M_B(i-1,j,5))*reflct_xL)/(dx*dx)
				+ (x(M_B(i,j+1,5))*reflct_yR - 2.0*x(M_B(i,j,5)) + x(M_B(i,j-1,5))*reflct_yL)/(dy*dy);
				double lapBy_res = (x(M_B(i+1,j,6))*reflct_xR - 2.0*x(M_B(i,j,6)) + x(M_B(i-1,j,6))*reflct_xL)/(dx*dx)
				+ (x(M_B(i,j+1,6))*reflct_yR - 2.0*x(M_B(i,j,6)) + x(M_B(i,j-1,6))*reflct_yL)/(dy*dy);
				double lapBz_res = (x(M_B(i+1,j,7))*reflct_xR - 2.0*x(M_B(i,j,7)) + x(M_B(i-1,j,7))*reflct_xL)/(dx*dx)
				+ (x(M_B(i,j+1,7))*reflct_yR - 2.0*x(M_B(i,j,7)) + x(M_B(i,j-1,7))*reflct_yL)/(dy*dy);
				Ax(M_B(i,j,5)) -= h * eta_resistive * lapBx_res;
				Ax(M_B(i,j,6)) -= h * eta_resistive * lapBy_res;
				Ax(M_B(i,j,7)) -= h * eta_resistive * lapBz_res;
			}
			
			// ---- 超电阻项(可选): dB/dt ⊃ -η_H∇⁴B, 用同一个各向同性 Laplacian 作用两次。
			// 警告: 内层 Laplacian 在 ±2 格处没有重新套用 reflct_xL/xR/yL/yR 镜像修正,
			// 在 Periodic 边界下没问题;如果是 Reflective/Dirichlet 强场问题,先关掉
			// add_hyperresistive 或者照 reflct_* 那套逻辑自己把镜像补全再开。----
			if(add_hyperresistive){
				auto lapB_at = [&](int ii, int jj, int comp) -> double {
					return (x(M_B(ii+1,jj,comp)) - 2.0*x(M_B(ii,jj,comp)) + x(M_B(ii-1,jj,comp)))/(dx*dx)
					+ (x(M_B(ii,jj+1,comp)) - 2.0*x(M_B(ii,jj,comp)) + x(M_B(ii,jj-1,comp)))/(dy*dy);
				};
				for(int comp = 5; comp <= 7; comp++){
					double biharm = (lapB_at(i+1,j,comp) - 2.0*lapB_at(i,j,comp) + lapB_at(i-1,j,comp))/(dx*dx)
					+ (lapB_at(i,j+1,comp) - 2.0*lapB_at(i,j,comp) + lapB_at(i,j-1,comp))/(dy*dy);
					Ax(M_B(i,j,comp)) += h * eta_hyper * biharm;
				}
			}
			
			if(GLF){
				nu_x = B_max_x * dx;
				nu_y = B_max_y * dy;
				lap_Bx = nu_x/(dx * dx) * (x(M_B(i-1,j,5)) * reflct_xL - 2.0*x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,5)) * reflct_yL - 2.0*x(M_B(i,j,5)) + x(M_B(i,j + 1,5)) * reflct_yR);
				Ax(M_B(i,j,5)) = Ax(M_B(i,j,5)) - h * lap_Bx;
				
				lap_By = nu_x/(dx * dx) * (x(M_B(i - 1,j,6)) * reflct_xL - 2.0*x(M_B(i,j,6)) + x(M_B(i + 1,j,6)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,6)) * reflct_yL - 2.0*x(M_B(i,j,6)) + x(M_B(i,j + 1,6)) * reflct_yR);
				
				Ax(M_B(i,j,6)) = Ax(M_B(i,j,6)) - h * lap_By;
				
				lap_Bz = nu_x/(dx * dx) * (x(M_B(i - 1,j,7)) * reflct_xL - 2.0*x(M_B(i,j,7)) + x(M_B(i + 1,j,7)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,7)) * reflct_yL - 2.0*x(M_B(i,j,7)) + x(M_B(i,j + 1,7)) * reflct_yR);
				Ax(M_B(i,j,7)) = Ax(M_B(i,j,7)) - h * lap_Bz;
			}else{
				double b_x = 0.5 * (fabs(Qe_pri(i,j,1)) + sqrt(Qe_pri(i,j,1) * Qe_pri(i,j,1) + 4 * (Qe_pri(i,j,5) * Qe_pri(i,j,5) + Qe_pri(i,j,6) * Qe_pri(i,j,6) + Qe_pri(i,j,7) * Qe_pri(i,j,7)) / Qe_pri(i,j,0)  ));
				
				double b_y = 0.5 * (fabs(Qe_pri(i,j,2)) + sqrt(Qe_pri(i,j,2) * Qe_pri(i,j,2) + 4 * (Qe_pri(i,j,5) * Qe_pri(i,j,5) + Qe_pri(i,j,6) * Qe_pri(i,j,6) + Qe_pri(i,j,7) * Qe_pri(i,j,7)) / Qe_pri(i,j,0)  ));
				
				double b_x_m = 0.5 * (fabs(Qe_pri(i - 1,j,1)) + sqrt(Qe_pri(i - 1,j,1) * Qe_pri(i - 1,j,1) + 4 * (Qe_pri(i - 1,j,5) * Qe_pri(i - 1,j,5) + Qe_pri(i - 1,j,6) * Qe_pri(i - 1,j,6) + Qe_pri(i - 1,j,7) * Qe_pri(i - 1,j,7)) / Qe_pri(i - 1,j,0)  ));
				
				double b_y_m = 0.5 * (fabs(Qe_pri(i,j - 1,2)) + sqrt(Qe_pri(i,j - 1,2) * Qe_pri(i,j - 1,2) + 4 * (Qe_pri(i,j - 1,5) * Qe_pri(i,j - 1,5) + Qe_pri(i,j - 1,6) * Qe_pri(i,j - 1,6) + Qe_pri(i,j - 1,7) * Qe_pri(i,j - 1,7)) / Qe_pri(i,j - 1,0)  ));
				
				double b_x_p = 0.5 * (fabs(Qe_pri(i + 1,j,1)) + sqrt(Qe_pri(i + 1,j,1) * Qe_pri(i + 1,j,1) + 4 * (Qe_pri(i + 1,j,5) * Qe_pri(i + 1,j,5) + Qe_pri(i + 1,j,6) * Qe_pri(i + 1,j,6) + Qe_pri(i + 1,j,7) * Qe_pri(i + 1,j,7)) / Qe_pri(i + 1,j,0)  ));
				
				double b_y_p = 0.5 * (fabs(Qe_pri(i,j + 1,2)) + sqrt(Qe_pri(i,j + 1,2) * Qe_pri(i,j + 1,2) + 4 * (Qe_pri(i,j + 1,5) * Qe_pri(i,j + 1,5) + Qe_pri(i,j + 1,6) * Qe_pri(i,j + 1,6) + Qe_pri(i,j + 1,7) * Qe_pri(i,j + 1,7)) / Qe_pri(i,j + 1,0)  ));
				
				double nu_x_p = fmax(fabs(b_x), fabs(b_x_p)) * dx;
				double nu_x_m = fmax(fabs(b_x), fabs(b_x_m)) * dx;
				double nu_y_p = fmax(fabs(b_y), fabs(b_y_p)) * dy;
				double nu_y_m = fmax(fabs(b_y), fabs(b_y_m)) * dy;
				
				lap_Bx = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,5)) * reflct_xR - x(M_B(i,j,5))) - nu_x_m * (x(M_B(i,j,5)) - x(M_B(i-1,j,5))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,5)) * reflct_yR -  x(M_B(i,j,5)))  - nu_y_m * (x(M_B(i,j,5)) - x(M_B(i,j - 1,5)) * reflct_yL));
				Ax(M_B(i,j,5)) = Ax(M_B(i,j,5)) - h * lap_Bx;
				
				lap_By = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,6)) * reflct_xR - x(M_B(i,j,6))) - nu_x_m * (x(M_B(i,j,6)) - x(M_B(i-1,j,6))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,6)) * reflct_yR -  x(M_B(i,j,6)))  - nu_y_m * (x(M_B(i,j,6)) - x(M_B(i,j - 1,6)) * reflct_yL));
				
				Ax(M_B(i,j,6)) = Ax(M_B(i,j,6)) - h * lap_By;
				
				lap_Bz = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,7)) * reflct_xR - x(M_B(i,j,7))) - nu_x_m * (x(M_B(i,j,7)) - x(M_B(i-1,j,7))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,7)) * reflct_yR -  x(M_B(i,j,7)))  - nu_y_m * (x(M_B(i,j,7)) - x(M_B(i,j - 1,7)) * reflct_yL));
				Ax(M_B(i,j,7)) = Ax(M_B(i,j,7)) - h * lap_Bz;
			}
		}
	}
}

void Magnetic_Subsystem_potential_AX(
	Eigen::VectorXd& Ax_out, Eigen::VectorXd& x,
	Grid& u_con, Grid& Qe, Grid& Qs, double h
) {
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qs);   // A-viscosity 特征值用 Qs(隐式累加器),对齐 CASTOR avisc=MagneticEigenvalues(Qs)

	// GLF: A-viscosity 系数全局化 —— 对齐 CASTOR (SI_MHD.f90 avisc(iDim,:)=MAXVAL(avisc(iDim,1:nElem)))
	// 每级从 Qs(=Qe_pri)在 real cell 上取方向全局 max,而非主循环里步首 u_pri 的 B_max_x/y。
	double gmax_bx = 0.0, gmax_by = 0.0;
	if(GLF){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Bsq = Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7);
			double bxv = 0.5*(fabs(Qe_pri(i,j,1))+sqrt(Qe_pri(i,j,1)*Qe_pri(i,j,1)+4.0*Bsq/Qe_pri(i,j,0)));
			double byv = 0.5*(fabs(Qe_pri(i,j,2))+sqrt(Qe_pri(i,j,2)*Qe_pri(i,j,2)+4.0*Bsq/Qe_pri(i,j,0)));
			if(bxv > gmax_bx) gmax_bx = bxv;
			if(byv > gmax_by) gmax_by = byv;
		}
	}

	// ★ 把GMRES向量复制到Grid, 做线性外推BC
	Grid AG(nxCells, nyCells);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			AG(i,j,0) = x((i-ghost_cells) + nxCells*(j-ghost_cells));
			AG(i,j,1) = x((i-ghost_cells) + nxCells*(j-ghost_cells) + nxCells*nyCells);
			AG(i,j,2) = x((i-ghost_cells) + nxCells*(j-ghost_cells) + 2*nxCells*nyCells);
		}
	bound_cond_A(AG);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Bxn_00 = Qe(i,j,5),   Byn_00 = Qe(i,j,6),   Bzn_00 = Qe(i,j,7);
			double Bxn_P0 = Qe(i+1,j,5), Byn_P0 = Qe(i+1,j,6), Bzn_P0 = Qe(i+1,j,7);
			double Bxn_M0 = Qe(i-1,j,5), Byn_M0 = Qe(i-1,j,6), Bzn_M0 = Qe(i-1,j,7);
			double Bxn_0P = Qe(i,j+1,5), Byn_0P = Qe(i,j+1,6), Bzn_0P = Qe(i,j+1,7);
			double Bxn_0M = Qe(i,j-1,5), Byn_0M = Qe(i,j-1,6), Bzn_0M = Qe(i,j-1,7);
			double rho_00 = u_con(i,j,0);
			
			// B^{n+1} = curl(A^{n+1}) — 全部从AG读
			double bx_00 = (AG(i,j+1,2) - AG(i,j-1,2)) / (2.0*dy);
			double by_00 = -(AG(i+1,j,2) - AG(i-1,j,2)) / (2.0*dx);
			double bz_00 = (AG(i+1,j,1) - AG(i-1,j,1)) / (2.0*dx)
			- (AG(i,j+1,0) - AG(i,j-1,0)) / (2.0*dy);
			
			double bx_P0 = (AG(i+1,j+1,2) - AG(i+1,j-1,2)) / (2.0*dy);
			double by_P0 = -(AG(i+2,j,2) - AG(i,j,2)) / (2.0*dx);
			double bz_P0 = (AG(i+2,j,1) - AG(i,j,1)) / (2.0*dx)
			- (AG(i+1,j+1,0) - AG(i+1,j-1,0)) / (2.0*dy);
			
			double bx_M0 = (AG(i-1,j+1,2) - AG(i-1,j-1,2)) / (2.0*dy);
			double by_M0 = -(AG(i,j,2) - AG(i-2,j,2)) / (2.0*dx);
			double bz_M0 = (AG(i,j,1) - AG(i-2,j,1)) / (2.0*dx)
			- (AG(i-1,j+1,0) - AG(i-1,j-1,0)) / (2.0*dy);
			
			double bx_0P = (AG(i,j+2,2) - AG(i,j,2)) / (2.0*dy);
			double by_0P = -(AG(i+1,j+1,2) - AG(i-1,j+1,2)) / (2.0*dx);
			double bz_0P = (AG(i+1,j+1,1) - AG(i-1,j+1,1)) / (2.0*dx)
			- (AG(i,j+2,0) - AG(i,j,0)) / (2.0*dy);
			
			double bx_0M = (AG(i,j,2) - AG(i,j-2,2)) / (2.0*dy);
			double by_0M = -(AG(i+1,j-1,2) - AG(i-1,j-1,2)) / (2.0*dx);
			double bz_0M = (AG(i+1,j-1,1) - AG(i-1,j-1,1)) / (2.0*dx)
			- (AG(i,j,0) - AG(i,j-2,0)) / (2.0*dy);
			
			// ==============================================================
			// 张量 T^B(α,β):  α 是动量分量(用 B_new=b), β 是通量方向(用 B_n)
			// 对应 Fortran MagneticFluxSI_momentum:
			//   FB(α,β) = B28Pe·δ_{αβ} - B_n^β · B_new^α     (μ0=1)
			//   其中 B28Pe = 0.5 · (B_n · B_new)
			// ==============================================================
			
			// —— x 方向通量 (β=x)，用在 P0, M0 位置 ——
			// 需要：T(2,1)=Txx, T(3,1)=Tyx, T(4,1)=Tzx
			double bdot_P0 = Bxn_P0*bx_P0 + Byn_P0*by_P0 + Bzn_P0*bz_P0;
			double Txx_P0 = bdot_P0/2.0 - Bxn_P0*bx_P0;   // 对角：α=β=x, 对称不变
			double Tyx_P0 = -Bxn_P0*by_P0;                 // α=y, β=x:  -B_n^x · B_new^y
			double Tzx_P0 = -Bxn_P0*bz_P0;                 // α=z, β=x:  -B_n^x · B_new^z
			
			double bdot_M0 = Bxn_M0*bx_M0 + Byn_M0*by_M0 + Bzn_M0*bz_M0;
			double Txx_M0 = bdot_M0/2.0 - Bxn_M0*bx_M0;
			double Tyx_M0 = -Bxn_M0*by_M0;
			double Tzx_M0 = -Bxn_M0*bz_M0;
			
			// —— y 方向通量 (β=y)，用在 0P, 0M 位置 ——
			// 需要：T(2,2)=Txy, T(3,2)=Tyy, T(4,2)=Tzy
			double bdot_0P = Bxn_0P*bx_0P + Byn_0P*by_0P + Bzn_0P*bz_0P;
			double Txy_0P = -Byn_0P*bx_0P;                 // α=x, β=y:  -B_n^y · B_new^x
			double Tyy_0P = bdot_0P/2.0 - Byn_0P*by_0P;   // 对角：α=β=y, 对称不变
			double Tzy_0P = -Byn_0P*bz_0P;                 // α=z, β=y:  -B_n^y · B_new^z
			
			double bdot_0M = Bxn_0M*bx_0M + Byn_0M*by_0M + Bzn_0M*bz_0M;
			double Txy_0M = -Byn_0M*bx_0M;
			double Tyy_0M = bdot_0M/2.0 - Byn_0M*by_0M;
			double Tzy_0M = -Byn_0M*bz_0M;
			
			// 散度: divT_α = ∂_x T(α,x) + ∂_y T(α,y)
			double divT_x = (Txx_P0-Txx_M0)/(2.0*dx) + (Txy_0P-Txy_0M)/(2.0*dy);
			double divT_y = (Tyx_P0-Tyx_M0)/(2.0*dx) + (Tyy_0P-Tyy_0M)/(2.0*dy);
			double divT_z = (Tzx_P0-Tzx_M0)/(2.0*dx) + (Tzy_0P-Tzy_0M)/(2.0*dy);
			
			double Fx = h/rho_00 * divT_x;
			double Fy = h/rho_00 * divT_y;
			double Fz = h/rho_00 * divT_z;
			
			double cross_x = Byn_00*Fz - Bzn_00*Fy;
			double cross_y = Bzn_00*Fx - Bxn_00*Fz;
			double cross_z = Bxn_00*Fy - Byn_00*Fx;
			
			// ---- implicit non-ideal terms on the vector potential ----
			// dA/dt = u x B - eta J - (d_i/rho) J x B  (Weyl gauge),  plus -eta_H lap^2 A,
			// which reproduces -eta_H lap^2 B exactly because curl commutes with lap^2.
			// The resistive part is written as +eta lap A, which likewise gives +eta lap B.
			// Here J^{n+1} = curl_h b with b = curl_h A^{n+1} already assembled above; note
			// this composition is a 2*Delta operator applied twice, so the A path carries a
			// wider Hall stencil (and a softer modified wavenumber) than the CT path.
			double hall_x = 0.0, hall_y = 0.0, hall_z = 0.0;
			if(add_hall && hall_implicit){
				double Jx_n = (bz_0P - bz_0M) / (2.0*dy);
				double Jy_n = -(bz_P0 - bz_M0) / (2.0*dx);
				double Jz_n = (by_P0 - by_M0) / (2.0*dx) - (bx_0P - bx_0M) / (2.0*dy);
				double rr = fmax(rho_00, 1e-12);
				hall_x = (Hall_di/rr) * (Jy_n*Bzn_00 - Jz_n*Byn_00);
				hall_y = (Hall_di/rr) * (Jz_n*Bxn_00 - Jx_n*Bzn_00);
				hall_z = (Hall_di/rr) * (Jx_n*Byn_00 - Jy_n*Bxn_00);
			}
			double res_x = 0.0, res_y = 0.0, res_z = 0.0;
			if(add_resistive && eta_resistive != 0.0){
				for(int k = 0; k < 3; k++){
					double lk = (AG(i+1,j,k) - 2.0*AG(i,j,k) + AG(i-1,j,k))/(dx*dx)
					          + (AG(i,j+1,k) - 2.0*AG(i,j,k) + AG(i,j-1,k))/(dy*dy);
					if(k==0) res_x = eta_resistive*lk;
					if(k==1) res_y = eta_resistive*lk;
					if(k==2) res_z = eta_resistive*lk;
				}
			}
			double hyp_x = 0.0, hyp_y = 0.0, hyp_z = 0.0;
			if(add_hyperresistive){
				auto lapA = [&](int ii,int jj,int k)->double{
					return (AG(ii+1,jj,k) - 2.0*AG(ii,jj,k) + AG(ii-1,jj,k))/(dx*dx)
					     + (AG(ii,jj+1,k) - 2.0*AG(ii,jj,k) + AG(ii,jj-1,k))/(dy*dy); };
				for(int k = 0; k < 3; k++){
					double bh = (lapA(i+1,j,k) - 2.0*lapA(i,j,k) + lapA(i-1,j,k))/(dx*dx)
					          + (lapA(i,j+1,k) - 2.0*lapA(i,j,k) + lapA(i,j-1,k))/(dy*dy);
					if(k==0) hyp_x = eta_hyper*bh;
					if(k==1) hyp_y = eta_hyper*bh;
					if(k==2) hyp_z = eta_hyper*bh;
				}
			}
			
			// Laplacian — 也从AG读
			double lap_Ax = 0.0, lap_Ay = 0.0, lap_Az = 0.0;
			if(GLF){
				double nu_x = gmax_bx * dx;   // 每级 Qs real-cell 全局 max(对齐 CASTOR GLF avisc)
				double nu_y = gmax_by * dy;
				lap_Ax = nu_x/(dx*dx)*(AG(i-1,j,0)-2.0*AG(i,j,0)+AG(i+1,j,0))
				+ nu_y/(dy*dy)*(AG(i,j-1,0)-2.0*AG(i,j,0)+AG(i,j+1,0));
				lap_Ay = nu_x/(dx*dx)*(AG(i-1,j,1)-2.0*AG(i,j,1)+AG(i+1,j,1))
				+ nu_y/(dy*dy)*(AG(i,j-1,1)-2.0*AG(i,j,1)+AG(i,j+1,1));
				lap_Az = nu_x/(dx*dx)*(AG(i-1,j,2)-2.0*AG(i,j,2)+AG(i+1,j,2))
				+ nu_y/(dy*dy)*(AG(i,j-1,2)-2.0*AG(i,j,2)+AG(i,j+1,2));
			}else{
				double Bsq = Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7);
				double b_x   = 0.5*(fabs(Qe_pri(i,j,1))+sqrt(Qe_pri(i,j,1)*Qe_pri(i,j,1)+4.0*Bsq/Qe_pri(i,j,0)));
				double b_x_p = 0.5*(fabs(Qe_pri(i+1,j,1))+sqrt(Qe_pri(i+1,j,1)*Qe_pri(i+1,j,1)+4.0*(Qe_pri(i+1,j,5)*Qe_pri(i+1,j,5)+Qe_pri(i+1,j,6)*Qe_pri(i+1,j,6)+Qe_pri(i+1,j,7)*Qe_pri(i+1,j,7))/Qe_pri(i+1,j,0)));
				double b_x_m = 0.5*(fabs(Qe_pri(i-1,j,1))+sqrt(Qe_pri(i-1,j,1)*Qe_pri(i-1,j,1)+4.0*(Qe_pri(i-1,j,5)*Qe_pri(i-1,j,5)+Qe_pri(i-1,j,6)*Qe_pri(i-1,j,6)+Qe_pri(i-1,j,7)*Qe_pri(i-1,j,7))/Qe_pri(i-1,j,0)));
				double b_y   = 0.5*(fabs(Qe_pri(i,j,2))+sqrt(Qe_pri(i,j,2)*Qe_pri(i,j,2)+4.0*Bsq/Qe_pri(i,j,0)));
				double b_y_p = 0.5*(fabs(Qe_pri(i,j+1,2))+sqrt(Qe_pri(i,j+1,2)*Qe_pri(i,j+1,2)+4.0*(Qe_pri(i,j+1,5)*Qe_pri(i,j+1,5)+Qe_pri(i,j+1,6)*Qe_pri(i,j+1,6)+Qe_pri(i,j+1,7)*Qe_pri(i,j+1,7))/Qe_pri(i,j+1,0)));
				double b_y_m = 0.5*(fabs(Qe_pri(i,j-1,2))+sqrt(Qe_pri(i,j-1,2)*Qe_pri(i,j-1,2)+4.0*(Qe_pri(i,j-1,5)*Qe_pri(i,j-1,5)+Qe_pri(i,j-1,6)*Qe_pri(i,j-1,6)+Qe_pri(i,j-1,7)*Qe_pri(i,j-1,7))/Qe_pri(i,j-1,0)));
				// 面平均(对齐 CASTOR ComputeCellLaplace: nup=0.5(avisc_i+avisc_{i+1}));原来用 max 不对
				double nu_xp = 0.5*(fabs(b_x)+fabs(b_x_p))*dx;
				double nu_xm = 0.5*(fabs(b_x)+fabs(b_x_m))*dx;
				double nu_yp = 0.5*(fabs(b_y)+fabs(b_y_p))*dy;
				double nu_ym = 0.5*(fabs(b_y)+fabs(b_y_m))*dy;
				for(int k = 0; k < 3; k++){
					double lk = 1.0/(dx*dx)*(nu_xp*(AG(i+1,j,k)-AG(i,j,k))-nu_xm*(AG(i,j,k)-AG(i-1,j,k)))
					+ 1.0/(dy*dy)*(nu_yp*(AG(i,j+1,k)-AG(i,j,k))-nu_ym*(AG(i,j,k)-AG(i,j-1,k)));
					if(k==0) lap_Ax=lk; if(k==1) lap_Ay=lk; if(k==2) lap_Az=lk;
				}
			}
			// 输出 — 直接用平坦索引, 不用M_B
			// A-solve numerical viscosity: OFF by default to match CASTOR (add_bvisc=.FALSE.).
			// 在低 β 强爆炸前沿,这个粘性会把 A 解搅坏 -> B 偏高 -> e=E-B^2/2 负压。
			if(!add_bvisc){ lap_Ax = 0.0; lap_Ay = 0.0; lap_Az = 0.0; }
			int idx = (i-ghost_cells) + nxCells*(j-ghost_cells);
			Ax_out(idx)                         = AG(i,j,0) - h*cross_x - h*lap_Ax
			                                      + h*hall_x - h*res_x + h*hyp_x;
			Ax_out(idx + nxCells*nyCells)        = AG(i,j,1) - h*cross_y - h*lap_Ay
			                                      + h*hall_y - h*res_y + h*hyp_y;
			Ax_out(idx + 2*nxCells*nyCells)      = AG(i,j,2) - h*cross_z - h*lap_Az
			                                      + h*hall_z - h*res_z + h*hyp_z;
		}
	}
}

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	Array shock_eta(nxCells,nyCells);
	ShockDetector(shock_eta, Qe);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);

			// Dirichlet: homogeneous part (ghost = -mirror) — 照抄 epsilon_BC 的 get_x_unknown
			// Dirichlet: -镜像(对齐参照 2DEuler_epsilon_BC get_x_unknown);lifting 用 1×(get_x_const)
			auto getp = [&](int ii, int jj) -> double {
				int II = ii, JJ = jj; double s = 1.0;
				if(x_bc == x_Dirichlet){
					if(ii < ghost_cells){ II = 2*ghost_cells - 1 - ii; s = -s; }
					else if(ii >= nxCells + ghost_cells){ II = 2*(nxCells+ghost_cells) - 1 - ii; s = -s; }
				}
				if(y_bc == y_Dirichlet){
					if(jj < ghost_cells){ JJ = 2*ghost_cells - 1 - jj; s = -s; }
					else if(jj >= nyCells + ghost_cells){ JJ = 2*(nyCells+ghost_cells) - 1 - jj; s = -s; }
				}
				return s * x(M_P(II, JJ));
			};
			double pii = x(M_P(i,j));

			// 强 Dirichlet:边界 cell 在 LHS 钉死 p=p_amb(恒等行),配 RHS b_p_vect=p_amb
			if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ){ Ax(M_P(i,j)) = pii; continue; }
			// eq (52): Lagrange degree-2 three-point stencil for d/dx(h dp/dx)
			// V(p) 对角项:floor_on 时用冻结斜率 V_slope(p_freeze)(此函数即作 Newton Jacobian);
			// 否则退回原线性斜率 1/(gamma-1)。耦合项(焓 Laplacian + 黏性)完全不变。
			double Vdiag = g_AX_jacobian ? V_slope( p_freeze(M_P(i,j)), u_con(i,j,0) ) : (1.0/(gamma_plasma - 1.0));
//			Ax(M_P(i,j)) = Vdiag * pii - h * h * (
//				1 / (dx * dx) * ( (0.75*hm0 + 0.25*hp0)*getp(i-1,j) - (hm0 + hp0)*pii + (0.25*hm0 + 0.75*hp0)*getp(i+1,j) )
//			  + 1 / (dy * dy) * ( (0.75*h0m + 0.25*h0p)*getp(i,j-1) - (h0m + h0p)*pii + (0.25*h0m + 0.75*h0p)*getp(i,j+1) ) );
			
			Ax(M_P(i,j)) = Vdiag * pii - h * h * (1 / (dx * dx) * (0.5 * (hp0 + h00) * (x(M_P(i + 1, j)) - x(M_P(i,j))) - 0.5 * (h00 + hm0) * (x(M_P(i,j)) - x(M_P(i-1,j)))) + 1 / (dy * dy) * (0.5 * (h0p + h00) * (x(M_P(i,j+1)) - x(M_P(i,j))) - 0.5 * (h0m + h00) * (x(M_P(i,j)) - x(M_P(i,j-1)))));
			
			//Laplacian Operator
			double nu_x, nu_y;
			double nu_x_p, nu_y_p;
			double nu_x_m, nu_y_m;
			
			double lap_p;
			if(GLF){
				switch(add_pvisc){
					case MaterialSpeedp:
						nu_x = vx_max * dx;
						nu_y = vy_max * dy;
						break;
					case AlfvenSpeedp:
						nu_x = cA_x * dx;
						nu_y = cA_y * dy;
						break;
					case Noaddvisp:
						nu_x = 0;
						nu_y = 0;
						break;
				}
				lap_p = nu_x/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * (getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * (pii - getp(i - 1,j))) + nu_y/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * (pii - getp(i,j-1)));
			}else{
				switch(add_pvisc){
					case MaterialSpeedp:
						nu_x_p = fmax(fabs(Qe_pri(i,j,1)), fabs(Qe_pri(i + 1,j,1))) * dx;
						nu_x_m = fmax(fabs(Qe_pri(i,j,1)), fabs(Qe_pri(i - 1,j,1))) * dx;
						nu_y_p = fmax(fabs(Qe_pri(i,j,2)), fabs(Qe_pri(i,j + 1,2))) * dy;
						nu_y_m = fmax(fabs(Qe_pri(i,j,2)), fabs(Qe_pri(i,j - 1,2))) * dy;
						break;
					case AlfvenSpeedp:
						nu_x_p = fmax(fabs(Qe_pri(i,j,5)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i + 1,j,5)  / sqrt(Qe_pri(i + 1,j,0)))) * dx;
						nu_x_m = fmax(fabs(Qe_pri(i,j,5)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i - 1,j,5)  / sqrt(Qe_pri(i - 1,j,0)))) * dx;
						nu_y_p = fmax(fabs(Qe_pri(i,j,6)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i,j + 1,6)  / sqrt(Qe_pri(i,j + 1,0)))) * dy;
						nu_y_m = fmax(fabs(Qe_pri(i,j,6)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i,j - 1,6)  / sqrt(Qe_pri(i,j - 1,0)))) * dy;
						break;
					case Noaddvisp:
						nu_x_p = 0;
						nu_x_m = 0;
						nu_y_p = 0;
						nu_y_m = 0;
						break;
				}
				lap_p =  1/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * nu_x_p *(getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * nu_x_m *(pii - getp(i - 1,j))) + 1/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * nu_y_p * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * nu_y_m * (pii - getp(i,j-1)));
			}
			// 黏性项(仅激波处激活的耗散)统一按局部 dV/dp = d(rho*e)/dp 缩放:
			// 理想气体 Vbase_slope == 1/(gamma-1)(逐字节不变),RK == rho*de/dp(冻结于 p_freeze,保持算子线性)。
			double visc_scale = Vbase_slope(Qe_pri(i,j,4), u_con(i,j,0));
			Ax(M_P(i,j)) = Ax(M_P(i,j)) - visc_scale * h * lap_p;
//			Ax(M_P(i,j)) = Ax(M_P(i,j)) - h * lap_p;
		}
	}
}

// ===================== 嵌套 Newton:保正压力求解 =====================
// 解非线性系统  g(p) = V_mod(p) + T p - b = 0,  其中 T = 焓 Laplacian + 黏性(线性,Qe 冻结)。
// 复用 Pressure_Subsystem_AX,由 g_AX_jacobian 选择算子形态:
//   g_AX_jacobian=false 时返回线性载体算子 L p = (1/(gamma-1)) p + T p;
//   g_AX_jacobian=true  时返回 Newton Jacobian  J x = V_slope(p_freeze,rho) x + T x。
// 残差 g = L(p) - (1/(gamma-1)) p + V_mod(p,rho) - b。载体常数 1/(gamma-1) 精确抵消,
// 提取出 EOS 无关的耦合 T;故对理想气体与 RK 同一套机制(论文 Eq.(52)-(53))。
void Pressure_NestedNewton(Eigen::VectorXd& p, Eigen::VectorXd& b,
                           Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim){
	// 仅在 "理想气体 + 关保正" 时退回单次线性 GMRES;RK 是非线性 EOS,必须走 Newton。
	if(!floor_on && eos_type==IDEAL_GAS){
		g_AX_jacobian = false;           // 线性算子 L(对角=1/(gamma-1))
		GMRES_Solver(p, b, Pressure_Subsystem_AX, u_con, Qe, Qs, h, dim);
		return;
	}
	const double a = 1.0/(gamma_plasma - 1.0);          // 载体常数(精确抵消)
	Eigen::VectorXd g(dim), dp(dim), neg_g(dim), Lp(dim);
	p_freeze = p;                        // Newton 初值 = 当前压力
	int it = 0;
	for(it = 0; it < Pnewton_max; it++){
		// ---- 残差 g(p_freeze):用线性载体算子 L ----
		g_AX_jacobian = false;
		Pressure_Subsystem_AX(Lp, p_freeze, u_con, Qe, Qs, h);
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				int idx = M_P(i,j);
				bool bnd = ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1))
				         || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1));
				if(bnd) g(idx) = Lp(idx) - b(idx);                          // 恒等行: p - p_amb
				else    g(idx) = Lp(idx) - a*p_freeze(idx) + V_mod(p_freeze(idx), u_con(i,j,0)) - b(idx);
			}
		}
		double res = g.norm();
		if(res < Pnewton_tol){
//			std::cout << "  [P-Newton] converged it=" << it << " res=" << res << std::endl;
			break;
		}
		// ---- Jacobian 解:J dp = -g(g_AX_jacobian=true,AX 即 Jacobian) ----
		neg_g = -g;
		dp.setZero();
		g_AX_jacobian = true;
		GMRES_Solver(dp, neg_g, Pressure_Subsystem_AX, u_con, Qe, Qs, h, dim);
		p_freeze += dp;
//		std::cout << "  [P-Newton] it=" << it << " res=" << res << " |dp|=" << dp.norm() << std::endl;
	}
	g_AX_jacobian = false;               // 复位,避免影响后续(理想气体)线性求解
	if(it == Pnewton_max)
		std::cout << "  [P-Newton] WARNING: hit max iters (" << Pnewton_max << ")" << std::endl;
	p = p_freeze;
}
// ====================================================================

void Energy_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			double h_00 = caltotalenthalpy(con00_n) *  con00_n[0]/ u_con(i,j,0);
			double h_p0 = caltotalenthalpy(conp0_n) *  conp0_n[0]/ u_con(i+1,j,0);
			double h_m0 = caltotalenthalpy(conm0_n) *  conm0_n[0]/ u_con(i-1,j,0);
			double h_0p = caltotalenthalpy(con0p_n) *  con0p_n[0]/ u_con(i,j+1,0);
			double h_0m = caltotalenthalpy(con0m_n) *  con0m_n[0]/ u_con(i,j-1,0);
			
			double hhp_x = 0.5 * (h_00 + h_p0);
			double hhm_x = 0.5 * (h_00 + h_m0);
			double hhp_y = 0.5 * (h_00 + h_0p);
			double hhm_y = 0.5 * (h_00 + h_0m);

			// Dirichlet: homogeneous part (ghost = -mirror) — 照抄 epsilon_BC 的 get_x_unknown
			auto getE = [&](int ii, int jj) -> double {   // CASTOR ghost-cell Dirichlet: 齐次 ghost=0(同 getp,避免 -镜像×2 角点漏)
				bool xg = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
				bool yg = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
				if(x_bc == x_Dirichlet && xg) return 0.0;
				if(y_bc == y_Dirichlet && yg) return 0.0;
				return x(M_P(ii, jj));
			};
			double Ei = x(M_P(i,j));
			Ax(M_P(i,j)) = Ei
			- (gamma_plasma - 1) * h * h * (
				1.0/(dx*dx) * ( hhp_x * (getE(i+1,j) - Ei) - hhm_x * (Ei - getE(i-1,j)) )
				+ 1.0/(dy*dy) * ( hhp_y * (getE(i,j+1) - Ei) - hhm_y * (Ei - getE(i,j-1)) )
			);

			Ax(M_P(i,j)) = Ax(M_P(i,j)) * dx * dy;
		}
	}
}


void GMRES_Solver(
	Eigen::VectorXd &x,
	Eigen::VectorXd &b,
	std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, double)> Ax_func,
	Grid &u_con, Grid &Qe, Grid &Qs, double h,
	const int dim,
	std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond)
{
	const int JMAX = 25;
	const int outmax = dim;
	// left preconditioner: applyM(out,in) = M^{-1} in  (identity if no precond)
	auto applyM = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		if(precond) precond(out, in); else out = in;
	};

	std::vector<Eigen::VectorXd> v(JMAX + 1, Eigen::VectorXd(dim));
	std::vector<Eigen::VectorXd> Av(JMAX + 1, Eigen::VectorXd(dim));
	Eigen::VectorXd w(dim), Ax(dim), r0(dim), tmp(dim);
	
	std::vector<std::vector<double>> H(JMAX + 2, std::vector<double>(JMAX + 2, 0.0));
	std::vector<double> alpha(JMAX + 2, 0.0);
	std::vector<double> gam(JMAX + 2, 0.0);
	std::vector<double> s(JMAX + 2, 0.0);
	std::vector<double> c(JMAX + 2, 0.0);
	
	double beta, pv;
	int j = 0;
	
	for (int outer = 0; outer < outmax; outer++) {
		
		Ax.setZero();
		Ax_func(Ax, x, u_con, Qe, Qs, h);  g_gmres_matvec++;
		tmp = b - Ax;
		applyM(r0, tmp);            // r0 = M^{-1}(b - A x)

		gam[0] = r0.norm();

		if (gam[0] <= GMREStol) {
//			std::cout << "GMRES converged: outer=" << outer
//			<< " residual=" << gam[0] << std::endl;
			return;
		}

		v[0] = r0 / gam[0];
		Ax.setZero();
		Ax_func(Ax, v[0], u_con, Qe, Qs, h);  g_gmres_matvec++;
		applyM(Av[0], Ax);          // Av[0] = M^{-1} A v[0]
		
		// Reset H
		for (int ii = 0; ii < JMAX + 2; ii++)
			std::fill(H[ii].begin(), H[ii].end(), 0.0);
		
		for (j = 0; j < JMAX; j++) {
			
			// Arnoldi: modified Gram-Schmidt
			for (int i = 0; i <= j; i++) {
				H[i][j] = v[i].dot(Av[j]);
			}
			w = Av[j];
			for (int i = 0; i <= j; i++) {
				w -= H[i][j] * v[i];
			}
			H[j + 1][j] = w.norm();
			
			// Apply previous Givens rotations
			for (int i = 0; i < j; i++) {
				double tmp1 = c[i + 1] * H[i][j] + s[i + 1] * H[i + 1][j];
				double tmp2 = s[i + 1] * H[i][j] - c[i + 1] * H[i + 1][j];
				H[i][j]     = tmp1;
				H[i + 1][j] = tmp2;
			}
			
			// New Givens rotation
			beta = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
			s[j + 1] = H[j + 1][j] / beta;
			c[j + 1] = H[j][j] / beta;
			H[j][j]  = beta;
			gam[j + 1] = s[j + 1] * gam[j];
			gam[j]     = c[j + 1] * gam[j];
			
			if (gam[j + 1] > GMREStol) {
				v[j + 1] = w / H[j + 1][j];
				Ax.setZero();
				Ax_func(Ax, v[j + 1], u_con, Qe, Qs, h);  g_gmres_matvec++;
				applyM(Av[j + 1], Ax);   // Av[j+1] = M^{-1} A v[j+1]
			} else {
				// Converged: back-substitute and update x
				for (int i = j; i >= 0; i--) {
					alpha[i] = gam[i];
					for (int k = i + 1; k <= j; k++) {
						alpha[i] -= H[i][k] * alpha[k];
					}
					alpha[i] /= H[i][i];
				}
				for (int i = 0; i <= j; i++) {
					x += alpha[i] * v[i];
				}
//				std::cout << "GMRES converged: outer=" << outer
//				<< " inner=" << j + 1
//				<< " residual=" << gam[j + 1] << std::endl;
				return;
			}
		}
		
		// Not converged: back-substitute and restart
		std::fill(alpha.begin(), alpha.end(), 0.0);
		for (int i = JMAX - 1; i >= 0; i--) {
			alpha[i] = gam[i];
			for (int k = i + 1; k < JMAX; k++) {
				alpha[i] -= H[i][k] * alpha[k];
			}
			alpha[i] /= H[i][i];
		}
		for (int i = 0; i < JMAX; i++) {
			x += alpha[i] * v[i];
		}
		
		if (outer % 10 == 0 || outer < 3) {
//			std::cout << "GMRES restart " << outer
//			<< " residual=" << gam[JMAX] << std::endl;
		}
	}
	
	std::cout << "WARNING: GMRES reached max outer iterations, residual=" << gam[0] << std::endl;
}

void MagSubsystem_UpdateA(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,Grid &A_s,double h, const int rks){
	Array shock_eta(nxCells,nyCells);
	double Ax_star, Ay_star, Az_star,U_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_U_vect(nxCells*nyCells), b_vect(3 * nxCells * nyCells), A_vect(3 * nxCells * nyCells), U_vect(nxCells * nyCells);
	
	
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			
			// B^n at cell (i,j)
			double Bxn = Qe(i,j,5);
			double Byn = Qe(i,j,6);
			double Bzn = Qe(i,j,7);
			
			// ρ^{n+1}
			double rho = u_con(i,j,0);
			
			// S = v^{**} - h/ρ * ∇p^n  (逐点, 用中心差分梯度G)
			double Sx = u_con(i,j,1) / rho
			- h / rho * (Qe_pri(i+1,j,4) - Qe_pri(i-1,j,4)) / (2.0 * dx);
			double Sy = u_con(i,j,2) / rho
			- h / rho * (Qe_pri(i,j+1,4) - Qe_pri(i,j-1,4)) / (2.0 * dy);
			double Sz = u_con(i,j,3) / rho;
			
			// A^{**} = A^n - h * (B^n × S)
			double Ax_star = A_s(i,j,0) - h * (Byn * Sz - Bzn * Sy);
			double Ay_star = A_s(i,j,1) - h * (Bzn * Sx - Bxn * Sz);
			double Az_star = A_s(i,j,2) - h * (Bxn * Sy - Byn * Sx);
			
			// 在 Ax_star/Ay_star/Az_star 那一段,紧跟着现有的 Sx,Sy,Sz 之后加:
			double Jx_h =  (Qe(i,j+1,7) - Qe(i,j-1,7)) / (2.0*dy);
			double Jy_h = -(Qe(i+1,j,7) - Qe(i-1,j,7)) / (2.0*dx);
			double Jz_h =  (Qe(i+1,j,6) - Qe(i-1,j,6)) / (2.0*dx) - (Qe(i,j+1,5) - Qe(i,j-1,5)) / (2.0*dy);
			
			double JxB_x = Jy_h*Bzn - Jz_h*Byn;
			double JxB_y = Jz_h*Bxn - Jx_h*Bzn;
			double JxB_z = Jx_h*Byn - Jy_h*Bxn;
			
			// When hall_implicit the term lives on the LHS of Magnetic_Subsystem_potential_AX
			// instead, so it must NOT also appear in the right-hand side here.
			if(add_hall && !hall_implicit){
				Ax_star -= h * (Hall_di/rho) * JxB_x;
				Ay_star -= h * (Hall_di/rho) * JxB_y;
				Az_star -= h * (Hall_di/rho) * JxB_z;
			}
			
			// 存入RHS向量
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells))                         = Ax_star;
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells)     = Ay_star;
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = Az_star;
			
			// 初始猜测 = A^n
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells))                         = A(i,j,0);
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells)     = A(i,j,1);
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = A(i,j,2);
		}
	}

	GMRES_Solver(A_vect, b_vect, Magnetic_Subsystem_potential_AX, u_con, Qe, Qs, h, 3*nxCells*nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			A(i,j,0) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells));
			A(i,j,1) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells);
			A(i,j,2) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells);
		}
	}
	
	//	bound_cond_u(A, 3);
	bound_cond_A(A);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) = 1 / (2.0 * dy) * (A(i,j + 1,2) - A(i,j - 1,2));
			u_con(i,j,6) = -1 / (2.0 * dx) * (A(i + 1,j,2) - A(i - 1,j,2));
			u_con(i,j,7) = 1 / (2.0 * dx) * (A(i + 1,j,1) - A(i - 1,j,1)) - 1 / (2.0 * dy) * (A(i,j + 1,0) - A(i,j - 1,0));
			
		}
	}
	bound_cond_u(u_con);
}

void PreSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks){
	Array shock_eta(nxCells,nyCells);
	double Bx_star, By_star, Bz_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_p_vect(nxCells*nyCells), b_Bt_vect(3 * nxCells * nyCells), B_vect(3 * nxCells * nyCells), p_vect(nxCells * nyCells);
	
	
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	
	
	
	bound_cond_u(u_con);
	
	Array rhovx_star(nxCells, nyCells);
	Array rhovy_star(nxCells, nyCells);
	Array rhovz_star(nxCells, nyCells);
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			rhovx_star(i,j) = u_con(i,j,1);
			rhovy_star(i,j) = u_con(i,j,2);
			rhovz_star(i,j) = u_con(i,j,3);
			
		}
	}
	bound_cond_array(rhovx_star);
	bound_cond_array(rhovy_star);
	bound_cond_array(rhovz_star);
	
	
//	std::cout << "rk = " <<rks + 1 <<", First Internal: "<<endl;
	vector<T> triplets_p;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double B2_00 = pow(u_con(i,j,5), 2.0) + pow(u_con(i,j,6), 2.0) + pow(u_con(i,j,7), 2.0);
			//			double B2_00 = pow(Qe(i,j,5), 2.0) + pow(Qe(i,j,6), 2.0) + pow(Qe(i,j,7), 2.0);
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
			
			double p_star = u_con(i,j,4) - 0.5 * B2_00;
			// ↓↓↓ 在这里插入 ↓↓↓
			double Jx_00 =  (u_con(i,j+1,7) - u_con(i,j-1,7)) / (2.0*dy);
			double Jy_00 = -(u_con(i+1,j,7) - u_con(i-1,j,7)) / (2.0*dx);
			double Jz_00 =  (u_con(i+1,j,6) - u_con(i-1,j,6)) / (2.0*dx) - (u_con(i,j+1,5) - u_con(i,j-1,5)) / (2.0*dy);
			double Joule_heating = add_resistive ? eta_resistive * (Jx_00*Jx_00 + Jy_00*Jy_00 + Jz_00*Jz_00) : 0.0;
			double internal = p_star - 0.5*(rhovx_star(i,j)*rhovx_star(i,j) + rhovy_star(i,j)*rhovy_star(i,j) + rhovz_star(i,j)*rhovz_star(i,j))/u_con(i,j,0);
			double RHS;

			RHS = p_star - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j) + rhovz_star(i,j) * rhovz_star(i,j)) / u_con(i,j,0) - h * (1 / dx * (0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) - 0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * (0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));
			RHS += h * Joule_heating;
//			if(internal < 1e-6){
//				std::cout << internal << " ";
//			}

			b_p_vect(M_P(i, j)) = RHS;

			// 强 Dirichlet:边界 cell 的 RHS = p_amb(配算子恒等行 Ax=pii -> 解出 p=p_amb)
			if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ) b_p_vect(M_P(i, j)) = p_amb;

			p_vect(M_P(i, j)) = Qe_pri(i,j,4);
		}
	}

//	std::cout << endl;


	Pressure_NestedNewton(p_vect, b_p_vect, u_con, Qe, Qs, h, nxCells*nyCells);
	
//	std::cout << "rk = " <<rks + 1<< ", First Solution: "<<std::endl;
//	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
//		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
//			if(p_vect(M_P(i,j))  < 1e-6){
//				std::cout << p_vect(M_P(i,j))  << " ";
//			}
//		}
//	}
//	std::cout << endl;
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			double fL, fR,	flux_xm, flux_xp;
			fL = p_vect(M_P(i - 1,j));
			fR = p_vect(M_P(i,j));

			flux_xm = 0.5 * (fL + fR);

			fL = p_vect(M_P(i,j));
			fR = p_vect(M_P(i + 1,j));

			flux_xp = 0.5 * (fL + fR);

			u_con(i,j,1) = rhovx_star(i,j) - h / dx * (flux_xp - flux_xm);

			double gL, gR, flux_ym, flux_yp;
			gL = p_vect(M_P(i,j - 1));
			gR = p_vect(M_P(i,j));

			flux_ym = 0.5 * (gL + gR);

			gL = p_vect(M_P(i,j));
			gR = p_vect(M_P(i,j + 1));

			flux_yp = 0.5 * (gL + gR);

			u_con(i,j,2) = rhovy_star(i,j) - h / dy * (flux_yp - flux_ym);
		}
	}
	
	
	bound_cond_u(u_con);
	
	
	
//	std::cout << "rk = " << rks + 1<< ", Second Internal: "<<endl;
	if(PressureSolver == all_correction || (PressureSolver == once_correction && rks == RKs - 1)){
		//std::cout <<"Second Pressure Solved "<< "rk: " << rks + 1 << endl;
		vector<T> triplets_p_second;
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
				double B2_00 = pow(u_con(i,j,5), 2.0) + pow(u_con(i,j,6), 2.0) + pow(u_con(i,j,7), 2.0);
				
				double con00_n[NUM_VARS];
				double conm0_n[NUM_VARS];
				double conp0_n[NUM_VARS];
				double con0m_n[NUM_VARS];
				double con0p_n[NUM_VARS];
				
				
				for(int k = 0; k < NUM_VARS; k++){
					con00_n[k] = Qe(i,j,k);
					conm0_n[k] = Qe(i - 1,j,k);
					conp0_n[k] = Qe(i + 1,j,k);
					con0m_n[k] = Qe(i,j - 1,k);
					con0p_n[k] = Qe(i,j + 1,k);
				}
				
				
				double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
				double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
				double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
				double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
				double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
				
				double p_star = u_con(i,j,4) - 0.5 * B2_00;
				// ↓↓↓ 在这里插入 ↓↓↓
				double Jx_00 =  (u_con(i,j+1,7) - u_con(i,j-1,7)) / (2.0*dy);
				double Jy_00 = -(u_con(i+1,j,7) - u_con(i-1,j,7)) / (2.0*dx);
				double Jz_00 =  (u_con(i+1,j,6) - u_con(i-1,j,6)) / (2.0*dx) - (u_con(i,j+1,5) - u_con(i,j-1,5)) / (2.0*dy);
				double Joule_heating = add_resistive ? eta_resistive * (Jx_00*Jx_00 + Jy_00*Jy_00 + Jz_00*Jz_00) : 0.0;
				
				double internal = p_star - 0.5*(rhovx_star(i,j)*rhovx_star(i,j) + rhovy_star(i,j)*rhovy_star(i,j) + rhovz_star(i,j)*rhovz_star(i,j))/u_con(i,j,0);
				
				
				double RHS;
				
				RHS = p_star - 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) - h * (1 / dx * (0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) - 0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * (0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));
				RHS += h * Joule_heating;
//				if(internal < 1e-8){
//					std::cout << internal << " ";
//				}
				b_p_vect(M_P(i, j)) = RHS;
				// 强 Dirichlet:边界 cell RHS = p_amb(配算子恒等行)
				if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ) b_p_vect(M_P(i, j)) = p_amb;

			}
		}
//		std::cout<<endl;

		Pressure_NestedNewton(p_vect, b_p_vect, u_con, Qe, Qs, h, nxCells*nyCells);
		
//		std::cout << "p Solution: " << std::endl;
//		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
//			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
//				if(p_vect(M_P(i,j))  < 0){
//					std::cout << p_vect(M_P(i,j))  << " ";
//				}
//			}
//		}
//		std::cout << endl;
	}
//	std::cout << "rk = " <<rks + 1<< ", Internal Energy: "<< endl;
	switch(U_update){
		case DefinitionofU:
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					u_con(i,j,4) = Vbase(p_vect(M_P(i,j)), u_con(i,j,0)) + 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) + 0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7));
					
					double internal_energy = u_con(i,j,4) - (0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7))) - (0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) );
					
//					if(internal_energy < 0){
//						std::cout << "E_total=" << u_con(i,j,4) << " ";
//						std::cout << "KE=" << (0.5 * (u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0)) << " ";
//						std::cout << "ME=" << (0.5 * (u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))) << " ";
//						std::cout << "IE=" << internal_energy;
//					}
				}
			}
//			std::cout<<endl;
			break;
		case ConservationofU:
			Array final_pressure(nxCells, nyCells);
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					final_pressure(i,j) = p_vect(M_P(i,j));
				}
			}
			bound_cond_array(final_pressure);
			double s_max;
			double priL[NUM_VARS], priR[NUM_VARS];
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					double fL, fR,	flux_xm, flux_xp;
					double gL, gR, flux_ym, flux_yp;
					
					double e_00 = eos_e_from_rho_p(u_con(i,j,0),     final_pressure(i,j));
					double e_p0 = eos_e_from_rho_p(u_con(i + 1,j,0), final_pressure(i + 1,j));
					double e_m0 = eos_e_from_rho_p(u_con(i - 1,j,0), final_pressure(i - 1,j));
					double e_0p = eos_e_from_rho_p(u_con(i,j + 1,0), final_pressure(i,j + 1));
					double e_0m = eos_e_from_rho_p(u_con(i,j - 1,0), final_pressure(i,j - 1));
					
					double h_00 = e_00 + final_pressure(i,j) / u_con(i,j,0);
					double h_p0 = e_p0 + final_pressure(i + 1,j) / u_con(i + 1,j,0);
					double h_m0 = e_m0 + final_pressure(i - 1,j) / u_con(i - 1,j,0);
					double h_0p = e_0p + final_pressure(i,j + 1) / u_con(i,j + 1,0);
					double h_0m = e_0m + final_pressure(i,j - 1) / u_con(i,j - 1,0);
					
					
					priL[0] = u_con(i - 1,j,0); priL[4] = final_pressure(i - 1,j);
					priL[1] = u_con(i - 1,j,1) / u_con(i - 1,j,0); priL[2] = u_con(i - 1,j,2) / u_con(i - 1,j,0); priL[3] = u_con(i - 1,j,3) / u_con(i - 1,j,0);
					priL[5] = u_con(i - 1,j,5); priL[6] = u_con(i - 1,j,6); priL[7] = u_con(i - 1,j,7);
					
					priR[0] = u_con(i,j,0); priR[4] = final_pressure(i,j);
					priR[1] = u_con(i,j,1) / u_con(i,j,0); priR[2] = u_con(i,j,2) / u_con(i,j,0); priR[3] = u_con(i,j,3) / u_con(i,j,0);
					priR[5] = u_con(i,j,5); priR[6] = u_con(i,j,6); priR[7] = u_con(i,j,7);
					fL = h_m0 * u_con(i-1,j,1);
					fR = h_00 * u_con(i,j,1);
					
					flux_xm = 0.5 * (fL + fR);
					
					priL[0] = u_con(i,j,0); priL[4] = final_pressure(i,j);
					priL[1] = u_con(i,j,1) / u_con(i,j,0); priL[2] = u_con(i,j,2) / u_con(i,j,0); priL[3] = u_con(i,j,3) / u_con(i,j,0);
					priL[5] = u_con(i,j,5); priL[6] = u_con(i,j,6); priL[7] = u_con(i,j,7);
					
					priR[0] = u_con(i + 1,j,0); priR[4] = final_pressure(i + 1,j);
					priR[1] = u_con(i + 1,j,1) / u_con(i + 1,j,0); priR[2] = u_con(i + 1,j,2) / u_con(i + 1,j,0); priR[3] = u_con(i + 1,j,3) / u_con(i + 1,j,0);
					priR[5] = u_con(i + 1,j,5); priR[6] = u_con(i + 1,j,6); priR[7] = u_con(i + 1,j,7);
					fL = h_00 * u_con(i,j,1);
					fR = h_p0 * u_con(i+1,j,1);
					
					flux_xp = 0.5 * (fL + fR);
					
					priL[0] = u_con(i,j - 1,0); priL[4] = final_pressure(i,j - 1);
					priL[1] = u_con(i,j - 1,1) / u_con(i,j - 1,0); priL[2] = u_con(i,j - 1,2) / u_con(i,j - 1,0); priL[3] = u_con(i,j - 1,3) / u_con(i,j - 1,0);
					priL[5] = u_con(i,j - 1,5); priL[6] = u_con(i,j - 1,6); priL[7] = u_con(i,j - 1,7);
					
					priR[0] = u_con(i,j,0); priR[4] = final_pressure(i,j);
					priR[1] = u_con(i,j,1) / u_con(i,j,0); priR[2] = u_con(i,j,2) / u_con(i,j,0); priR[3] = u_con(i,j,3) / u_con(i,j,0);
					priR[5] = u_con(i,j,5); priR[6] = u_con(i,j,6); priR[7] = u_con(i,j,7);
					
					gL = h_0m * u_con(i,j - 1,2);
					gR = h_00 * u_con(i,j,2);
					
					flux_ym = 0.5 * (gL + gR);
					
					priL[0] = u_con(i,j,0); priL[4] = final_pressure(i,j);
					priL[1] = u_con(i,j,1) / u_con(i,j,0); priL[2] = u_con(i,j,2) / u_con(i,j,0); priL[3] = u_con(i,j,3) / u_con(i,j,0);
					priL[5] = u_con(i,j,5); priL[6] = u_con(i,j,6); priL[7] = u_con(i,j,7);
					
					priR[0] = u_con(i,j + 1,0); priR[4] = final_pressure(i,j + 1);
					priR[1] = u_con(i,j + 1,1) / u_con(i,j + 1,0); priR[2] = u_con(i,j + 1,2) / u_con(i,j + 1,0); priR[3] = u_con(i,j + 1,3) / u_con(i,j + 1,0);
					priR[5] = u_con(i,j + 1,5); priR[6] = u_con(i,j + 1,6); priR[7] = u_con(i,j + 1,7);
					
					gL = h_00 * u_con(i,j,2);
					gR = h_0p * u_con(i,j + 1,2);
					
					flux_yp = 0.5 * (gL + gR);
					
					u_con(i,j,4) = u_con(i,j,4) - h / (dx) * (flux_xp - flux_xm) - h / (dy) * (flux_yp - flux_ym);
					
					
					double internal_energy = u_con(i,j,4) - (0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7))) - (0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) );
					
					if(internal_energy < 0){
						std::cout << "E_total=" << u_con(i,j,4) << " ";
						std::cout << "KE=" << (0.5 * (u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0)) << " ";
						std::cout << "ME=" << (0.5 * (u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))) << " ";
						std::cout << "IE=" << internal_energy;
					}
				}
			}
			std::cout<<endl;
			break;
	}
	
	bound_cond_u(u_con);
	
	
	b_p_vect.resize(0);
	b_Bt_vect.resize(0);
	B_vect.resize(0);
	p_vect.resize(0);
	bound_cond_u(u_con);
}

void PreSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	double U_star;
	Eigen::VectorXd b_U_vect(nxCells*nyCells), U_vect(nxCells * nyCells);
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	std::cout << "rk = " <<rks + 1 <<", RHS: "<<endl;
	vector<T> triplets_U;
	SpMat M_U(nxCells * nyCells, nxCells * nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			double h_00 = caltotalenthalpy(con00_n) *  con00_n[0]/ u_con(i,j,0);
			double h_p0 = caltotalenthalpy(conp0_n) *  conp0_n[0]/ u_con(i+1,j,0);
			double h_m0 = caltotalenthalpy(conm0_n) *  conm0_n[0]/ u_con(i-1,j,0);
			double h_0p = caltotalenthalpy(con0p_n) *  con0p_n[0]/ u_con(i,j+1,0);
			double h_0m = caltotalenthalpy(con0m_n) *  con0m_n[0]/ u_con(i,j-1,0);
			
			U_star = u_con(i,j,4) - 0.5 * h / dx * (h_p0 * u_con(i + 1,j,1) - h_m0 * u_con(i - 1,j,1)) - 0.5 * h / dy * (h_0p * u_con(i,j + 1,2) - h_0m * u_con(i,j - 1,2));


			b_U_vect(M_P(i, j)) = U_star * dx * dy;

			// Dirichlet lifting: subtract operator applied to (interior=0, ghost=E_bc)
			if(x_bc == x_Dirichlet || y_bc == y_Dirichlet){
				double hhp_x = 0.5*(h_00+h_p0), hhm_x = 0.5*(h_00+h_m0);
				double hhp_y = 0.5*(h_00+h_0p), hhm_y = 0.5*(h_00+h_0m);
				// ghost-cell Dirichlet 1×(配 getE ghost=0);钉内能=p_amb -> E = p_amb/(g-1) + 1/2 B_ghost^2
				// (用透射 B,使内能恒=ambient,B 偏离 ambient 时也一致)
				auto getEc = [&](int ii, int jj) -> double {
					bool xg = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
					bool yg = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
					bool g = (x_bc == x_Dirichlet && xg) || (y_bc == y_Dirichlet && yg);
					if(!g) return 0.0;
					return p_amb/(gamma_plasma-1.0)
					     + 0.5*(u_con(ii,jj,5)*u_con(ii,jj,5)+u_con(ii,jj,6)*u_con(ii,jj,6)+u_con(ii,jj,7)*u_con(ii,jj,7));
				};
				double L_const = - (gamma_plasma - 1) * h * h * (
					1.0/(dx*dx) * ( hhp_x * getEc(i+1,j) + hhm_x * getEc(i-1,j) )
					+ 1.0/(dy*dy) * ( hhp_y * getEc(i,j+1) + hhm_y * getEc(i,j-1) )
				);
				L_const *= dx * dy;
				b_U_vect(M_P(i, j)) -= L_const;
			}
			U_vect(M_P(i,j)) = u_con(i,j,4); // 初始猜测
			
			if(U_star < 1e-6){
				std::cout << U_star << " ";
			}
		}
		
	}
	std::cout << endl;
	// Diagonal (Jacobi) preconditioner for the energy solve (matches CASTOR precopDiag_SIMHD_energy)
	auto precE = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		out.resize(in.size());
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double c00[NUM_VARS], cp0[NUM_VARS], cm0[NUM_VARS], c0p[NUM_VARS], c0m[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					c00[k]=Qe(i,j,k); cp0[k]=Qe(i+1,j,k); cm0[k]=Qe(i-1,j,k); c0p[k]=Qe(i,j+1,k); c0m[k]=Qe(i,j-1,k);
				}
				double h00 = caltotalenthalpy(c00)*c00[0]/u_con(i,j,0);
				double hp0 = caltotalenthalpy(cp0)*cp0[0]/u_con(i+1,j,0);
				double hm0 = caltotalenthalpy(cm0)*cm0[0]/u_con(i-1,j,0);
				double h0p = caltotalenthalpy(c0p)*c0p[0]/u_con(i,j+1,0);
				double h0m = caltotalenthalpy(c0m)*c0m[0]/u_con(i,j-1,0);
				double hhp_x=0.5*(h00+hp0), hhm_x=0.5*(h00+hm0);
				double hhp_y=0.5*(h00+h0p), hhm_y=0.5*(h00+h0m);
				double diag = 1.0 + (gamma_plasma-1.0)*h*h*( (hhp_x+hhm_x)/(dx*dx) + (hhp_y+hhm_y)/(dy*dy) );
				out(M_P(i,j)) = in(M_P(i,j)) / (diag * dx * dy);
			}
		}
	};
	GMRES_Solver(U_vect, b_U_vect, Energy_Subsystem_AX, u_con, Qe, Qs, h, nxCells*nyCells, precE);

	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			u_con(i,j,4) = U_vect(M_P(i, j));
		}
	}
	bound_cond_u(u_con);
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			u_con(i,j,1) -= h * (1/(2*dx)) * (gamma_plasma-1) * (u_con(i+1,j,4) - u_con(i-1,j,4));
			u_con(i,j,2) -= h * (1/(2*dy)) * (gamma_plasma-1) * (u_con(i,j+1,4) - u_con(i,j-1,4));
		}
	}
	
	bound_cond_u(u_con);
	
	std::cout << "rk = " <<rks + 1 <<", Internal Energy: "<<endl;
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double internal = u_con(i,j,4) - 0.5*(u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0) -  0.5*(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7));
			if(internal < 1e-6){
				std::cout << internal << " ";
			}
		}
	}
	std::cout << endl;
}


void MagSubsystem_UpdateB(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,double h, const int rks){
	double Bx_star, By_star, Bz_star, U_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_U_vect(nxCells*nyCells), b_Bt_vect(3 * nxCells * nyCells), B_vect(3 * nxCells * nyCells), U_vect(nxCells * nyCells);
	
	
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double rho_p0 = 0.5 * (u_con(i,j,0) + u_con(i + 1,j,0));
			double rho_m0 = 0.5 * (u_con(i,j,0) + u_con(i - 1,j,0));
			double rho_0p = 0.5 * (u_con(i,j,0) + u_con(i,j + 1,0));
			double rho_0m = 0.5 * (u_con(i,j,0) + u_con(i,j - 1,0));
			double rho_pp = 0.25 * (u_con(i,j,0) + u_con(i + 1,j,0) + u_con(i,j + 1,0) + u_con(i + 1,j + 1,0));
			double rho_pm = 0.25 * (u_con(i,j,0) + u_con(i + 1,j,0) + u_con(i,j - 1,0) + u_con(i + 1,j - 1,0));
			double rho_mp = 0.25 * (u_con(i,j,0) + u_con(i - 1,j,0) + u_con(i,j + 1,0) + u_con(i - 1,j + 1,0));
			double rho_mm = 0.25 * (u_con(i,j,0) + u_con(i - 1,j,0) + u_con(i,j - 1,0) + u_con(i - 1,j - 1,0));
			
			double rhovxstar_p0 = 0.5 * (u_con(i,j,1) + u_con(i + 1,j,1));
			double rhovxstar_m0 = 0.5 * (u_con(i,j,1) + u_con(i - 1,j,1));
			double rhovxstar_0p = 0.5 * (u_con(i,j,1) + u_con(i,j + 1,1));
			double rhovxstar_0m = 0.5 * (u_con(i,j,1) + u_con(i,j - 1,1));
			double rhovxstar_pp = 0.25 * (u_con(i,j,1) + u_con(i + 1,j,1) + u_con(i,j + 1,1) + u_con(i + 1,j + 1,1));
			double rhovxstar_pm = 0.25 * (u_con(i,j,1) + u_con(i + 1,j,1) + u_con(i,j - 1,1) + u_con(i + 1,j - 1,1));
			double rhovxstar_mp = 0.25 * (u_con(i,j,1) + u_con(i - 1,j,1) + u_con(i,j + 1,1) + u_con(i - 1,j + 1,1));
			double rhovxstar_mm = 0.25 * (u_con(i,j,1) + u_con(i - 1,j,1) + u_con(i,j - 1,1) + u_con(i - 1,j - 1,1));
			
			double rhovystar_p0 = 0.5 * (u_con(i,j,2) + u_con(i + 1,j,2));
			double rhovystar_m0 = 0.5 * (u_con(i,j,2) + u_con(i - 1,j,2));
			double rhovystar_0p = 0.5 * (u_con(i,j,2) + u_con(i,j + 1,2));
			double rhovystar_0m = 0.5 * (u_con(i,j,2) + u_con(i,j - 1,2));
			double rhovystar_pp = 0.25 * (u_con(i,j,2) + u_con(i + 1,j,2) + u_con(i,j + 1,2) + u_con(i + 1,j + 1,2));
			double rhovystar_pm = 0.25 * (u_con(i,j,2) + u_con(i + 1,j,2) + u_con(i,j - 1,2) + u_con(i + 1,j - 1,2));
			double rhovystar_mp = 0.25 * (u_con(i,j,2) + u_con(i - 1,j,2) + u_con(i,j + 1,2) + u_con(i - 1,j + 1,2));
			double rhovystar_mm = 0.25 * (u_con(i,j,2) + u_con(i - 1,j,2) + u_con(i,j - 1,2) + u_con(i - 1,j - 1,2));
			
			double rhovzstar_p0 = 0.5 * (u_con(i,j,3) + u_con(i + 1,j,3));
			double rhovzstar_m0 = 0.5 * (u_con(i,j,3) + u_con(i - 1,j,3));
			double rhovzstar_0p = 0.5 * (u_con(i,j,3) + u_con(i,j + 1,3));
			double rhovzstar_0m = 0.5 * (u_con(i,j,3) + u_con(i,j - 1,3));
			double rhovzstar_pp = 0.25 * (u_con(i,j,3) + u_con(i + 1,j,3) + u_con(i,j + 1,3) + u_con(i + 1,j + 1,3));
			double rhovzstar_pm = 0.25 * (u_con(i,j,3) + u_con(i + 1,j,3) + u_con(i,j - 1,3) + u_con(i + 1,j - 1,3));
			double rhovzstar_mp = 0.25 * (u_con(i,j,3) + u_con(i - 1,j,3) + u_con(i,j + 1,3) + u_con(i - 1,j + 1,3));
			double rhovzstar_mm = 0.25 * (u_con(i,j,3) + u_con(i - 1,j,3) + u_con(i,j - 1,3) + u_con(i - 1,j - 1,3));
			
			
			double p_p0 = 0.5 * (Qe_pri(i,j,4) + Qe_pri(i + 1,j,4));
			double p_m0 = 0.5 * (Qe_pri(i,j,4) + Qe_pri(i - 1,j,4));
			double p_0p = 0.5 * (Qe_pri(i,j,4) + Qe_pri(i,j + 1,4));
			double p_0m = 0.5 * (Qe_pri(i,j,4) + Qe_pri(i,j - 1,4));
			double p_pp = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i + 1,j,4) + Qe_pri(i,j + 1,4) + Qe_pri(i + 1,j + 1,4));
			double p_pm = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i + 1,j,4) + Qe_pri(i,j - 1,4) + Qe_pri(i + 1,j - 1,4));
			double p_mp = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i - 1,j,4) + Qe_pri(i,j + 1,4) + Qe_pri(i - 1,j + 1,4));
			double p_mm = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i - 1,j,4) + Qe_pri(i,j - 1,4) + Qe_pri(i - 1,j - 1,4));
			
			double Bxn_p0 = 0.5 * (Qe(i,j,5) + Qe(i + 1,j,5));
			double Bxn_m0 = 0.5 * (Qe(i,j,5) + Qe(i - 1,j,5));
			double Bxn_0p = 0.5 * (Qe(i,j,5) + Qe(i,j + 1,5));
			double Bxn_0m = 0.5 * (Qe(i,j,5) + Qe(i,j - 1,5));
			double Bxn_pp = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j + 1,5) + Qe(i + 1,j + 1,5));
			double Bxn_pm = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j - 1,5) + Qe(i + 1,j - 1,5));
			double Bxn_mp = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j + 1,5) + Qe(i - 1,j + 1,5));
			double Bxn_mm = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j - 1,5) + Qe(i - 1,j - 1,5));
			
			double Byn_p0 = 0.5 * (Qe(i,j,6) + Qe(i + 1,j,6));
			double Byn_m0 = 0.5 * (Qe(i,j,6) + Qe(i - 1,j,6));
			double Byn_0p = 0.5 * (Qe(i,j,6) + Qe(i,j + 1,6));
			double Byn_0m = 0.5 * (Qe(i,j,6) + Qe(i,j - 1,6));
			double Byn_pp = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j + 1,6) + Qe(i + 1,j + 1,6));
			double Byn_pm = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j - 1,6) + Qe(i + 1,j - 1,6));
			double Byn_mp = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j + 1,6) + Qe(i - 1,j + 1,6));
			double Byn_mm = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j - 1,6) + Qe(i - 1,j - 1,6));
			
			double Bzn_p0 = 0.5 * (Qe(i,j,7) + Qe(i + 1,j,7));
			double Bzn_m0 = 0.5 * (Qe(i,j,7) + Qe(i - 1,j,7));
			double Bzn_0p = 0.5 * (Qe(i,j,7) + Qe(i,j + 1,7));
			double Bzn_0m = 0.5 * (Qe(i,j,7) + Qe(i,j - 1,7));
			double Bzn_pp = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j + 1,7) + Qe(i + 1,j + 1,7));
			double Bzn_pm = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j - 1,7) + Qe(i + 1,j - 1,7));
			double Bzn_mp = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j + 1,7) + Qe(i - 1,j + 1,7));
			double Bzn_mm = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j - 1,7) + Qe(i - 1,j - 1,7));
			
			
			
			Bx_star = Qs(i,j,5) - (h / dy) * ((Bxn_0p / rho_0p * (rhovystar_0p - (h / dy) * (Qe_pri(i,j + 1,4) - Qe_pri(i,j,4))) - Byn_0p / rho_0p * (rhovxstar_0p - h / dx * (p_pp - p_mp)  )) - (Bxn_0m / rho_0m * (rhovystar_0m - (h / dy) * (Qe_pri(i,j,4) - Qe_pri(i,j - 1,4))) - Byn_0m / rho_0m * (rhovxstar_0m - h / dx * (p_pm - p_mm)  )) ); 
			
			By_star = Qs(i,j,6) - (h / dx) * ((Byn_p0 / rho_p0 * (rhovxstar_p0 - (h / dx) * (Qe_pri(i + 1,j,4) - Qe_pri(i,j,4))) - Bxn_p0 / rho_p0 * (rhovystar_p0 - h / dy * (p_pp - p_pm) ) ) - (Byn_m0/rho_m0 * (rhovxstar_m0 - (h / dx) * (Qe_pri(i,j,4) - Qe_pri(i - 1,j,4))) - Bxn_m0 / rho_m0 * (rhovystar_m0 - h / dy * (p_mp - p_mm)  ) ));
			
			double Bz_star_1 = ((Bzn_p0 / rho_p0 * (rhovxstar_p0 - h / dx * (Qe_pri(i + 1,j,4) - Qe_pri(i,j,4))) - Bxn_p0 / rho_p0 * (rhovzstar_p0) ) - (Bzn_m0 / rho_m0 * (rhovxstar_m0 - h / dx * (Qe_pri(i,j,4) - Qe_pri(i - 1,j,4))) - Bxn_m0 / rho_m0 * (rhovzstar_m0)) );
			
			double Bz_star_2 = ((Bzn_0p / rho_0p * (rhovystar_0p - h / dy * (Qe_pri(i,j + 1,4) - Qe_pri(i,j,4))) - Byn_0p / rho_0p * (rhovzstar_0p) ) - (Bzn_0m / rho_0m * (rhovystar_0m - h / dy * (Qe_pri(i,j,4) - Qe_pri(i,j - 1,4))) - Byn_0m / rho_0m * (rhovzstar_0m)) );
			
			Bz_star = Qs(i,j,7) - (h / dx) * Bz_star_1 - (h / dy) * Bz_star_2;
			
			
			//Laplacian Operator
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) =  (Bx_star);
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells) = (By_star);
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = (Bz_star);
			
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = Qe(i,j,5);
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells) = Qe(i,j,6);
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = Qe(i,j,7);
			
		}
	}
	GMRES_Solver(B_vect, b_Bt_vect, Magnetic_Subsystem_AX, u_con, Qe, Qs, h, 3*nxCells*nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells));
			u_con(i,j,6) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells);
			u_con(i,j,7) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells);
		}
	}
	bound_cond_u(u_con);
}

// Legacy EXPLICIT Hall operator, used only when hall_implicit == false, as a reference for
// the implicit path.  The previous version of this routine was wrong in two ways: it rotated
// the in-plane field with a pointwise Cayley transform of dB/dt = -(d_i/rho) J x B, which is
// the vector-POTENTIAL equation (E_H itself) rather than the induction equation (-curl E_H),
// and whatever it did write into Bx,By was overwritten immediately afterwards by the CT face
// average.  It now applies -h * curl_h E_H consistently to all three components.
//
// Note this remains forward Euler on a purely dispersive term, so it has no stability region
// on the imaginary axis for the two-stage explicit tableau used here; it will need
// hall_omega_cap well below the RK4 value 2.83 to survive, and is not the production path.
void HallMHD_SourceRelax(Grid& u_con, double h){
	Grid Bnew(nxCells, nyCells);

	auto EH = [&](int ii, int jj, int comp)->double{
		double Jx, Jy, Jz;
		CellCurrent(u_con, ii, jj, Jx, Jy, Jz);
		double r  = fmax(u_con(ii,jj,0), 1e-12);
		double Bx = u_con(ii,jj,5), By = u_con(ii,jj,6), Bz = u_con(ii,jj,7);
		if(comp == 0) return (Hall_di/r) * (Jy*Bz - Jz*By);
		if(comp == 1) return (Hall_di/r) * (Jz*Bx - Jx*Bz);
		return                (Hall_di/r) * (Jx*By - Jy*Bx);
	};

	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double curl_x =  (EH(i,j+1,2) - EH(i,j-1,2)) / (2.0*dy);
			double curl_y = -(EH(i+1,j,2) - EH(i-1,j,2)) / (2.0*dx);
			double curl_z =  (EH(i+1,j,1) - EH(i-1,j,1)) / (2.0*dx)
			              -  (EH(i,j+1,0) - EH(i,j-1,0)) / (2.0*dy);
			Bnew(i,j,5) = u_con(i,j,5) - h * curl_x;
			Bnew(i,j,6) = u_con(i,j,6) - h * curl_y;
			Bnew(i,j,7) = u_con(i,j,7) - h * curl_z;
		}
	}
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) = Bnew(i,j,5);
			u_con(i,j,6) = Bnew(i,j,6);
			u_con(i,j,7) = Bnew(i,j,7);
		}
	}
	bound_cond_u(u_con);
}

double cal_divergence_error(Grid &u_con, Array &Bx_0p, Array &By_p0){
	// 每种 div-free 方案测它自己保证为零的散度定义,且测当前真实 B(不是旧的交错场):
	//   StaggerCT                        -> 交错面场 (Bx_p0,By_0p) 的有限体积散度
	//   PotentialA / UnStaggerCT / NoDivFree -> cell-centered B = curl A 的散度,用 CASTOR
	//        的 ComputeCellGradientA 同款梯度算子:内部中心差分,域边界线性外推=单侧差分。
	// 注意:边界**不能**用中心差分去碰 bound_cond_u 的透射 ghost B(≠curl A),否则
	//       边界 div(curl)≠0 且随时间增长(透射 ghost 与 curl A 不一致)。CASTOR 用同一个
	//       外推梯度算子构造 B 和测 divB,故边界也保持一致。
	// (旧版固定用 Bx_p0/By_0p,对后两种方案是没更新的初始值 -> 假的 0)
	double error = 0;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double d;
			if(magfield == StaggerCT){
				d = (Bx_0p(i,j) - Bx_0p(i-1,j)) / dx + (By_p0(i,j) - By_p0(i,j-1)) / dy;
			}else{
				// CASTOR ComputeCellGradientA: 内部中心差分;域边界线性外推 -> 单侧差分
				double dBxdx, dBydy;
				if(i == ghost_cells)                    dBxdx = (u_con(i+1,j,5) - u_con(i,j,5)) / dx;
				else if(i == nxCells + ghost_cells - 1) dBxdx = (u_con(i,j,5) - u_con(i-1,j,5)) / dx;
				else                                    dBxdx = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx);
				if(j == ghost_cells)                    dBydy = (u_con(i,j+1,6) - u_con(i,j,6)) / dy;
				else if(j == nyCells + ghost_cells - 1) dBydy = (u_con(i,j,6) - u_con(i,j-1,6)) / dy;
				else                                    dBydy = (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
				d = dBxdx + dBydy;
			}
			error += dx * dy * d * d;
		}
	}
	return sqrt(error);
}

void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime){
	bool isPlotTime = false;
	std::string testname;
	switch(test){
		case test0: 
			testname = "test0";
			break;
		case test1: 
			testname = "test1";
			break;
		case test2: 
			testname = "test2";
			break;
		case test3: 
			testname = "test3";
			break;
		case test4: 
			testname = "test4";
			break;
		case test5:
			testname = "test5";
			break;
		case test6:
			testname = "test6";
			break;
		case test7:
			testname = "test7";
			break;
		case test8:
			testname = "test8";
			break;
		case test9:
			testname = "test9";
			break;
		case test10:
			testname = "test10";
			break;
		
		case shocktest:
			testname = "shocktest";
			break;
		case convergencetest:
			testname = "convergencetest";
			break;
		case test2d0:
			testname = "test2d0";
			break;
		case test1y:
			testname = "test1y";
			break;
		case testOrszag_Tang:
			testname = "Orszag_Tang";
			break;
		case testblastwave:
			testname = "blast_wave";
			break;
		case testblastwave1D:
			testname = "blast_wave1D";
			break;
		case testRotor:
			testname = "Rotor";
			break;
		case testFieldLoop:
			testname = "FieldLoop" + std::to_string(A0) ;
			break;
		case testKHInstability:
			testname = "KHInstability" + std::to_string(A0);
			break;
		case testdoubleshear:
			testname = "doubleshear";
			break;
		case testRK_Riemann:
			testname = "RK_Riemann";
			break;
		case testPlasmaEquilibrium:
			testname = "PlasmaEquilibrium";
			break;
		case testAdvectedPlasmaEquilibrium:
			testname = "AdvectedPlasmaEquilibrium";
			break;
		case testMHD_vortex:
			testname = "MHD_vortex";
			break;
		case testHallReconnection:
			testname = "testHallReconnection";
			break;
		case testHarrisHall:
			testname = "testHarrisHall";
			break;
	}
	
	for(int k=0; k<(int)IntermediatePlotTimes.size(); k++){
		double Tdiff = std::abs(IntermediatePlotTimes[k]-t);
		if(Tdiff <= dt/2.0){isPlotTime = true; break;}
	}
	if(iteration == 0) isPlotTime = true;
	if(iteration % plotINTERVAL == 0) isPlotTime = true;
	if(isPlotTime || isFinalTime){
		std::cout<<"Intermediate plot"<< std::endl;
		printAllVariablesToVTK(u_con,u_pri, testname,iteration,isFinalTime);
	}
}

/* VISIT OUTPUT FILE CREATION */

/*Prints the header section for a VTK file. */


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
	
	for(int i=1; i<nxCells + 2 ; i++){
		outFile << (((double)i)-0.5)*dx << " ";
	}
	
	outFile <<  std::endl;
	
	outFile << "Y_COORDINATES " << nyCells + 1 <<" double" << std::endl;
	
	for(int j=1; j<nyCells+2 ; j++){
		outFile << (((double)j)-0.5)*dy << " ";
	}
	
	outFile << std::endl;
	
	outFile << "Z_COORDINATES 1 double" << std::endl;
	
	outFile << "0" << std::endl;
	
	outFile << "CELL_DATA " << nxCells*nyCells << std::endl;
	
	return;
}

/** Prints a scalar variable into a VTK file.
*/
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName){
	outFile << variableName << " 1 " << nxCells * nyCells << " double" << std::endl;
	
	int printWidth = 24;
	
	for(int j = 0; j < nyCells; j++){
		for(int i = 0; i < nxCells; i++){
			double varPrint;
			if(var < NUM_VARS){
				varPrint = u_con(i,j,var);
			}else if(var < NUM_VARS + 4){
				varPrint = u_pri(i,j,var - NUM_VARS + 1);
			}else{
				varPrint = sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6)  + u_con(i,j,7) * u_con(i,j,7) );
			}
			outFile << std::setw(printWidth) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << varPrint << " ";
		}
		outFile << std::endl;
	}
	
	outFile << std::endl;
}

/**Prints as many variables as desired to a single VTK file.
*/
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime){
	string formula_type;
	string div_type;
	string update_type;
	string dt_type;
	string eos;
	string positive_preserving;
	switch (recon) {
		case NoRecon:
			SpaceOrder = 1;
			break;
		case Minmod: case VanLeer:
			SpaceOrder = 2;
			break;
		case ENO3:
			SpaceOrder = 3;
			break;
		case WENO5:
			SpaceOrder = 4;
			break;
	}
	
	switch (formula) {
		case EnergyFormula:
			formula_type = "EnergyFormula";
			update_type = "";
			break;
		case PressureFormula:
			formula_type = "PressureFormula";
			switch(U_update){
				case ConservationofU:
					update_type = "conservation";
					break;
				case DefinitionofU:
					update_type = "definition";
					break;
			}
			break;
	}
	
	switch(magfield){
		case StaggerCT:
			div_type = "StaggerCT";
			break;
		case UnStaggerCT:
			div_type = "UnStaggerCT";
			break;
		case PotentialA:
			div_type = "Potential";
			break;
		case NoDivFree:
			div_type = "NoDivFree";
			break;
	}
	
	switch(caltimestep){
		case Fullydt:
			dt_type = "Fullydt";
			break;
		case Alfvendt:
			dt_type = "Alfvendt";
			break;
		case IMEXdt:
			dt_type = "IMEXdt" + std::to_string(full_step);
			break;
	}
	
	switch(eos_type){
		case IDEAL_GAS:
			eos = "Ideal Gas";
			break;
		case REDLICH_KWONG:
			eos = "Redlich Kwong";
			break;
	}
	
	if(floor_on){
		positive_preserving = "PreservingON";
	}else{
		positive_preserving = "PreservingOFF";
	}
	if(nxCells == 1 || nyCells == 1){
		div_type = "";
	}
	
	
	
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(0) << 100*CFL;
	std::string s = ss.str();  
	
	Path = "/Users/chenzefeng/Desktop/Study/CamSC/phd/research/HallMHD/Plot/" + eos  +"/" +positive_preserving +"/"+ testname + "/" + formula_type + "/"+ div_type + "/"+ update_type + "/" + std::to_string(nxCells) + "x"+ std::to_string(nyCells) + "_"+ dt_type+"_RK" + std::to_string(RKs) + "_Order" + std::to_string(SpaceOrder) +"_CFL0" + s;
	
	filesystem::create_directories(Path);
	
	std::string fileString = Path + "/" + "plt_" + testname; // "./Plot/" + name;
	
	fileString.append("_");
	fileString.append(to_string(nstep).substr(0,6));
	//if(isFinalTime) fileString.append("_FINAL");
	fileString.append(".vtk");
	
	std::ofstream outFile(fileString.c_str());
		
	
	printVTKHeader(outFile,nstep);
	
	std::vector< std::tuple<int,std::string,int,int,int> > variablesToPrint;
	
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(0,std::string("rho"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(1,std::string("mom_X"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(2,std::string("mom_Y"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(3,std::string("mom_Z"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(4,"U",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(5,std::string("B_X"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(6,std::string("B_Y"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(7,std::string("B_Z"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(8,"velX",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(9,"velY",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(10,"velZ",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(11,"p",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(12,"B",0,0,0));
	
	
	outFile << "FIELD FieldData " << variablesToPrint.size() << std::endl;
	
	for(unsigned int i = 0;i <variablesToPrint.size(); i++){
		printVTKVariable(outFile,u_con,u_pri,std::get<0>(variablesToPrint[i]),std::get<1>(variablesToPrint[i]));
	}
	
	
	outFile << "POINT_DATA " << (nxCells+1)*(nyCells+1) << std::endl;
	
	return;
}


// ============================================================================================
// Reconnection diagnostics, report Sec 7.4.  Reconstructs the flux function, tracks the X/O
// pair, and reports the reconnection rate against the Cassak-Shay scale.  Method A (d psi/dt)
// is the reported rate; Method B (|E_z| sampled in the X cell) is the consistency cross-check.
// ============================================================================================
struct HallDiagRec {
	double t, psi, AzX, xX, yX, Ez_X, vout, Bz_absmax, q1, q2, q3, q4, rho_min, p_min;
	double vy_max, Jz_absmax, By_sheet;
	double Jz_signed_max;   // max(J_z) WITH sign.  The reference plots this, not max|J_z|:
	                        // a Harris sheet has J_z < 0 throughout, so max(J_z) starts near
	                        // zero while max|J_z| starts at B0/lambda.  They are different
	                        // diagnostics and must not be compared with each other.
	int    xo_ok;           // 1 = X/O identification passed its consistency checks
};
static std::vector<HallDiagRec> g_hall_diag;

void HallReconnectionDiagnostics(Grid& u_con, double tnow){
	const int i0 = ghost_cells, i1 = nxCells + ghost_cells - 1;
	const int j0 = ghost_cells, j1 = nyCells + ghost_cells - 1;

	// ---- flux function A_z along the two-segment path of Sec 7.4.1 ----
	Array Az(nxCells, nyCells);
	Az(i0, j0) = 0.0;
	for(int j = j0 + 1; j <= j1; j++)
		Az(i0,j) = Az(i0,j-1) + 0.5*(u_con(i0,j-1,5) + u_con(i0,j,5)) * dy;
	for(int j = j0; j <= j1; j++)
		for(int i = i0 + 1; i <= i1; i++)
			Az(i,j) = Az(i-1,j) - 0.5*(u_con(i-1,j,6) + u_con(i,j,6)) * dx;
	// The reconstruction fixes an arbitrary additive constant.  Anchor it on the WEAK-field
	// upstream corner rather than the strong-field one: the strong-field boundary is where the
	// transmissive condition misbehaves first, and anchoring there injects that drift straight
	// into the X-point rate below.
	{
		const double anchor = Az(i0, j1);
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++) Az(i,j) -= anchor;
	}

	// ---- sheet centre: row of largest x-averaged |J_z| within 5 d_i of the neutral line ----
	int jX = j0; double best = -1.0;
	for(int j = j0; j <= j1; j++){
		double y = y_0 + (j - ghost_cells + 0.5)*dy;
		double yNL = (test == testHarrisHall) ? 0.0 : hall_yNL;
		if(fabs(y - yNL) > 5.0*Hall_di) continue;
		double acc = 0.0;
		for(int i = i0; i <= i1; i++){
			double Jx, Jy, Jz; CellCurrent(u_con, i, j, Jx, Jy, Jz);
			acc += fabs(Jz);
		}
		if(acc > best){ best = acc; jX = j; }
	}

	// ---- O/X candidates on that row: global max and min of A_z ----
	int iMax = i0, iMin = i0;
	for(int i = i0; i <= i1; i++){
		if(Az(i,jX) > Az(iMax,jX)) iMax = i;
		if(Az(i,jX) < Az(iMin,jX)) iMin = i;
	}
	double psi = Az(iMax,jX) - Az(iMin,jX);          // Eq. (7.21), label independent

	// ---- X / O identification ------------------------------------------------------------
	// The X point is the SADDLE of the flux function and the O point an extremum, so the
	// signs of the two second derivatives separate them unambiguously:
	//     A_xx * A_yy < 0  -> saddle  -> X
	//     A_xx * A_yy > 0  -> extremum -> O
	// The earlier rule "the candidate with the larger |J_z| is X" only works when current has
	// already concentrated at the reconnection site.  In a SYMMETRIC Harris sheet at early
	// times the current is essentially uniform along the layer, the two candidates are
	// indistinguishable by |J_z|, and the rule silently selected the O point at the domain
	// edge -- which then corrupted every quantity sampled "at the X point".
	auto curv = [&](int i, int j, double& axx, double& ayy){
		int im = (i > i0) ? i-1 : i0, ip = (i < i1) ? i+1 : i1;
		int jm = (j > j0) ? j-1 : j0, jp = (j < j1) ? j+1 : j1;
		axx = (Az(ip,j) - 2.0*Az(i,j) + Az(im,j)) / (dx*dx);
		ayy = (Az(i,jp) - 2.0*Az(i,j) + Az(i,jm)) / (dy*dy);
	};
	double axxA, ayyA, axxB, ayyB;
	curv(iMax, jX, axxA, ayyA);
	curv(iMin, jX, axxB, ayyB);
	bool maxIsSaddle = (axxA*ayyA < 0.0), minIsSaddle = (axxB*ayyB < 0.0);
	int iX; int xo_ok = 1;
	if(maxIsSaddle && !minIsSaddle)      iX = iMax;
	else if(minIsSaddle && !maxIsSaddle) iX = iMin;
	else{
		// degenerate: both or neither look like a saddle.  Fall back to the current-based
		// rule and flag the sample so downstream statistics can drop it.
		double JxA, JyA, JzA, JxB2, JyB2, JzB2;
		CellCurrent(u_con, iMax, jX, JxA, JyA, JzA);
		CellCurrent(u_con, iMin, jX, JxB2, JyB2, JzB2);
		iX = (fabs(JzA) > fabs(JzB2)) ? iMax : iMin;
		xo_ok = 0;
	}
	double xX = x_0 + (iX - ghost_cells + 0.5)*dx;
	double yX = y_0 + (jX - ghost_cells + 0.5)*dy;

	// ---- Method B: total out-of-plane electric field sampled in the X cell, Eq. (7.28) ----
	double Jx, Jy, Jz; CellCurrent(u_con, iX, jX, Jx, Jy, Jz);
	double rX  = fmax(u_con(iX,jX,0), 1e-12);
	double vxX = u_con(iX,jX,1)/rX, vyX = u_con(iX,jX,2)/rX;
	double Ez  = vyX*u_con(iX,jX,5) - vxX*u_con(iX,jX,6);
	if(add_resistive) Ez += eta_resistive * Jz;
	if(add_hall)      Ez += (Hall_di/rX) * (Jx*u_con(iX,jX,6) - Jy*u_con(iX,jX,5));
	if(add_hyperresistive){
		double jzp, jzm, jzu, jzd, a, b, c;
		CellCurrent(u_con, iX+1, jX, a,b,jzp); CellCurrent(u_con, iX-1, jX, a,b,jzm);
		CellCurrent(u_con, iX, jX+1, a,b,jzu); CellCurrent(u_con, iX, jX-1, a,b,jzd);
		Ez += -eta_hyper * ((jzp - 2.0*Jz + jzm)/(dx*dx) + (jzu - 2.0*Jz + jzd)/(dy*dy));
	}

	// ---- outflow: max |v_x| within |y - y_X| <= 2.0 = 4 lambda ----
	double vout = 0.0, rho_min = 1e30, p_min = 1e30, Bz_absmax = 0.0;
	double vy_max = 0.0, Jz_absmax = 0.0, By_sheet = 0.0;
	double Jz_signed_max = -1e30;
	double q[4] = {0,0,0,0}; int qn[4] = {0,0,0,0};
	for(int i = i0; i <= i1; i++){
		for(int j = j0; j <= j1; j++){
			double xx = x_0 + (i - ghost_cells + 0.5)*dx;
			double yy = y_0 + (j - ghost_cells + 0.5)*dy;
			double r  = u_con(i,j,0);
			if(r < rho_min) rho_min = r;
			double ke = 0.5*(u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2)
			               + u_con(i,j,3)*u_con(i,j,3))/fmax(r,1e-12);
			double me = 0.5*(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6)
			               + u_con(i,j,7)*u_con(i,j,7));
			double pp = (gamma_plasma - 1.0)*(u_con(i,j,4) - ke - me);
			if(pp < p_min) p_min = pp;
			if(fabs(u_con(i,j,7)) > Bz_absmax) Bz_absmax = fabs(u_con(i,j,7));
			if(fabs(yy - yX) <= 2.0){
				vout   = fmax(vout,   fabs(u_con(i,j,1)/fmax(r,1e-12)));
				vy_max = fmax(vy_max, fabs(u_con(i,j,2)/fmax(r,1e-12)));
			}
			{	double ja, jb, jc; CellCurrent(u_con, i, j, ja, jb, jc);
				Jz_absmax     = fmax(Jz_absmax, fabs(jc));
				Jz_signed_max = fmax(Jz_signed_max, jc); }
			if(j == jX) By_sheet = fmax(By_sheet, fabs(u_con(i,j,6)));
			// quadrupolar signature: mean B_z per quadrant in a 2 d_i box around the X point
			if(fabs(xx - xX) <= 2.0*Hall_di && fabs(yy - yX) <= 2.0*Hall_di){
				int k = (xx > xX ? 0 : 1) + (yy > yX ? 0 : 2);
				q[k] += u_con(i,j,7); qn[k]++;
			}
		}
	}
	// A quadrant box that runs off the domain returns a partial (or empty) average, which
	// would silently read as a small or exactly-zero lobe.  Require all four boxes to be
	// populated and comparable in size before trusting the quadrupole signature.
	int qn_min = qn[0], qn_max = qn[0];
	for(int k = 1; k < 4; k++){ qn_min = std::min(qn_min, qn[k]); qn_max = std::max(qn_max, qn[k]); }
	if(qn_min == 0 || qn_max > 2*std::max(qn_min,1)) xo_ok = 0;
	for(int k = 0; k < 4; k++) if(qn[k]) q[k] /= qn[k];

	HallDiagRec rec = {tnow, psi, Az(iX,jX), xX, yX, fabs(Ez), vout, Bz_absmax,
	                   q[0], q[1], q[2], q[3], rho_min, p_min,
	                   vy_max, Jz_absmax, By_sheet, Jz_signed_max, xo_ok};
	g_hall_diag.push_back(rec);

	// ---- Cassak-Shay normalisation, Eqs. (7.34)-(7.37) ----
	// Upstream states: symmetric for Test 12, asymmetric for the survey case.  The
	// Cassak-Shay expressions reduce to B_eff = B0 and v_out,mix = vA_inf when B1 = B2.
	double B1u = hall_B1, B2u = hall_B2, r1u = hall_rho1, r2u = hall_rho2;
	if(test == testHarrisHall){
		B1u = B2u = h12_B0; r1u = r2u = h12_rho_inf;
	}
	double Beff = 2.0*B1u*B2u/(B1u + B2u);
	double vmix = sqrt(B1u*B2u*(B1u + B2u) / (r1u*B2u + r2u*B1u));
	double Ecs = Beff*vmix, Epred = 0.1*Ecs;

	// ---- rate estimators ------------------------------------------------------------------
	// The two protocols in the report use DIFFERENT estimators and must not be mixed:
	//
	//   survey (Sec 7.2.4, reference code case 27) : E = |-dA_z/dt| sampled AT the X point,
	//        normalised by E_CS = B_eff v_out,mix, with the first-peak selection rule below.
	//   Double-Harris sequence (Eq. 7.21/7.26)     : E = |d psi/dt| with psi = max_x A_z - min_x A_z.
	//
	// Both are written out.  Note that in the Weyl gauge dA_z/dt = -E_z exactly, so the survey
	// estimator and the Method-B sample |E_z(x_X)| should agree; their difference measures the
	// gauge drift of the reconstruction plus the motion of the X point across cells, and is a
	// useful consistency check rather than an independent physical quantity.
	// A centred difference of A_z(X) is only a rate if the X CELL did not move meanwhile: the
	// identification can swap X and O, or hop to a competing island, and then the difference
	// records the displacement of the sample point rather than any reconnected flux.  Samples
	// whose X cell moved by more than hall_xjump_cells between the two neighbours are masked,
	// which is the in-code version of the offline masking described in report Sec 8.4.2.
	const double t_transit = 0.5*(x_1 - x_0) / vmix;
	std::vector<double> tt, EsurvV, EpsiV;
	std::vector<int> okV;
	for(size_t n = 1; n + 1 < g_hall_diag.size(); n++){
		const HallDiagRec& a2 = g_hall_diag[n-1];
		const HallDiagRec& b2 = g_hall_diag[n+1];
		tt.push_back(g_hall_diag[n].t);
		EsurvV.push_back(fabs((b2.AzX - a2.AzX)/(b2.t - a2.t)));
		EpsiV .push_back(fabs((b2.psi - a2.psi)/(b2.t - a2.t)));
		bool ok = (fabs(b2.xX - a2.xX) <= hall_xjump_cells*dx)
		       && (fabs(b2.yX - a2.yX) <= hall_xjump_cells*dy)
		       && a2.xo_ok && b2.xo_ok && g_hall_diag[n].xo_ok;
		okV.push_back(ok ? 1 : 0);
	}

	// first-peak rule: the first interior local maximum before t_transit that reaches at least
	// half the in-window maximum; otherwise the window maximum is used.
	double Ewin = 0.0; int nmask = 0;
	for(size_t n = 0; n < tt.size(); n++){
		if(!okV[n]){ nmask++; continue; }
		if(tt[n] <= t_transit) Ewin = fmax(Ewin, EsurvV[n]);
	}
	double Efirst = Ewin; double t_first = -1.0; bool got = false;
	for(size_t n = 1; n + 1 < tt.size() && !got; n++){
		if(tt[n] > t_transit) break;
		if(!okV[n] || !okV[n-1] || !okV[n+1]) continue;
		if(EsurvV[n] >= EsurvV[n-1] && EsurvV[n] >= EsurvV[n+1] && EsurvV[n] >= 0.5*Ewin){
			Efirst = EsurvV[n]; t_first = tt[n]; got = true;
		}
	}

	std::ofstream f("hall_reconnection.dat");
	f << "# report Sec 7.4 diagnostics for testHallReconnection\n";
	f << "# B1=" << hall_B1 << " B2=" << hall_B2 << " rho1=" << hall_rho1
	  << " rho2=" << hall_rho2 << "  b=" << hall_B2/hall_B1 << "\n";
	f << "# B_eff=" << Beff << " v_out,mix=" << vmix
	  << " E_CS=" << Ecs << " E_pred=0.1*E_CS=" << Epred
	  << " t_transit=" << t_transit << "\n";
	f << "# first peak (survey rule): E=" << Efirst << " E/E_CS=" << Efirst/Ecs
	  << " at t=" << t_first << (got ? " (interior local max)" : " (window max fallback)")
	  << "   masked samples (X-cell jump): " << nmask << "/" << tt.size() << "\n";
	// Machine-readable column names.  Plotting scripts MUST parse this line rather than
	// hard-coding indices; inserting a column has silently broken downstream plots before.
	f << "#cols t E_surv E_surv_over_ECS E_X E_psi E_psi_over_Epred Phi AzX valid"
	     " xX yX v_out vy_max Jz_signed_max Jz_absmax By_sheet Bz_absmax"
	     " q_pp q_mp q_pm q_mm rho_min p_min\n";
	f << "# Phi = max_x Az - min_x Az on the sheet row (Eq. 7.21).\n";
	f << "# Jz_signed_max = max(Jz) WITH sign  <-- this is the quantity the reference plots.\n";
	f << "# Jz_absmax     = max|Jz|            <-- NOT the same curve; do not cross-compare.\n";
	f << "# valid = 1 when the X/O identification passed (A_z saddle test + unclipped"
	     " quadrant boxes) and the X cell did not jump between neighbouring outputs.\n";
	f << "# q_pp/q_mp/q_pm/q_mm = mean Bz in the (+x,+y)/(-x,+y)/(+x,-y)/(-x,-y) quadrants"
	     " of a 2 d_i box centred on the X point.\n";
	f << std::scientific << std::setprecision(8);
	for(size_t n = 0; n < tt.size(); n++){
		const HallDiagRec& m = g_hall_diag[n+1];
		f << m.t << " " << EsurvV[n] << " " << EsurvV[n]/Ecs << " " << m.Ez_X << " "
		  << EpsiV[n] << " " << EpsiV[n]/Epred << " "
		  << m.psi << " " << m.AzX << " "
		  << okV[n] << " "
		  << m.xX << " " << m.yX << " " << m.vout << " "
		  << m.vy_max << " " << m.Jz_signed_max << " " << m.Jz_absmax << " "
		  << m.By_sheet << " " << m.Bz_absmax << " "
		  << m.q1 << " " << m.q2 << " " << m.q3 << " " << m.q4 << " "
		  << m.rho_min << " " << m.p_min << "\n";
	}
	f.close();

	// ---- optional binary field snapshot: header (t, nx, ny, x0, y0, dx, dy) then arrays ----
	if(hall_snap_on){
		static bool snap_init = false;
		std::ofstream fb("hall_fields.bin",
		    std::ios::binary | (snap_init ? std::ios::app : std::ios::trunc));
		snap_init = true;
		double hdr[7] = {tnow, (double)nxCells, (double)nyCells, x_0, y_0, dx, dy};
		fb.write((char*)hdr, sizeof(hdr));
		std::vector<double> buf(nxCells*nyCells);
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++)
			buf[(i-i0) + nxCells*(j-j0)] = Az(i,j);
		fb.write((char*)buf.data(), buf.size()*sizeof(double));
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++)
			buf[(i-i0) + nxCells*(j-j0)] = u_con(i,j,7);              // Bz
		fb.write((char*)buf.data(), buf.size()*sizeof(double));
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++){
			double a,b,c; CellCurrent(u_con, i, j, a, b, c);
			buf[(i-i0) + nxCells*(j-j0)] = c;                          // Jz
		}
		fb.write((char*)buf.data(), buf.size()*sizeof(double));
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++)
			buf[(i-i0) + nxCells*(j-j0)] = u_con(i,j,1)/fmax(u_con(i,j,0),1e-12);  // vx
		fb.write((char*)buf.data(), buf.size()*sizeof(double));
		for(int j = j0; j <= j1; j++) for(int i = i0; i <= i1; i++)
			buf[(i-i0) + nxCells*(j-j0)] = u_con(i,j,0);              // rho
		fb.write((char*)buf.data(), buf.size()*sizeof(double));
		fb.close();
	}

	double Enow = 0.0;
	size_t N = g_hall_diag.size();
	if(N >= 3) Enow = fabs((g_hall_diag[N-1].AzX - g_hall_diag[N-3].AzX)
	                     / (g_hall_diag[N-1].t   - g_hall_diag[N-3].t));
	std::cout << "HALL t=" << tnow << " E_surv/E_CS=" << (Ecs > 0 ? Enow/Ecs : 0.0)
	          << " E_X/E_CS=" << (Ecs > 0 ? fabs(Ez)/Ecs : 0.0)
	          << " psi=" << psi
	          << " X=(" << xX << "," << yX << ")"
	          << " vout=" << vout << " vy=" << vy_max
	          << " maxJz=" << Jz_signed_max << " |Jz|max=" << Jz_absmax
	          << " By=" << By_sheet << " xoOK=" << xo_ok
	          << " Phi=" << psi << " |Bz|max=" << Bz_absmax
	          << " quad=[" << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << "]"
	          << "  [firstpeak E/E_CS=" << Efirst/Ecs << " @t=" << t_first << "]"
	          << std::endl;
}

int main(int argc, char *argv[]){
	switch (test) {
		case test0:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test1: case test1y:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test2:
			x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test3:
			//x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test4:
			//x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test5:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.04;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test6:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.03;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test7:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test8:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test9:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.08;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test10:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case shocktest:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.8;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case convergencetest:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 1.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case test2d0:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testOrszag_Tang:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 1.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			gamma_plasma = 1.4;
			break;
		case testblastwave:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.01;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;   // match CASTOR: rho/v/p Dirichlet(ambient) + B transmissive
			y_bc = y_Transmissive;
			break;
		case testblastwave1D:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.01;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testRotor:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testFieldLoop:
			x_0 = -1, x_1 = 1, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 1;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testKHInstability:
			x_0 = 0, x_1 = 2, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = (4.8 / A0 ) / 6.0;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testdoubleshear:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 1.8;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testRK_Riemann:
			x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1, z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testPlasmaEquilibrium:
			x_0 = -3.0, x_1 = 3.0, y_0 = -3.0, y_1 = 3.0, z_0 = -1, z_1 = 1;
			time_stop = 2.0;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Transmissive;     // 静态平衡,零流;透射即可
			y_bc = y_Transmissive;
			break;
		case testAdvectedPlasmaEquilibrium:
			x_0 = -3.0, x_1 = 3.0, y_0 = -3.0, y_1 = 3.0, z_0 = -1, z_1 = 1;
			time_stop = 400.0;         // 平流绕域一圈
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;         // 平流版用周期(让环绕回来)
			y_bc = y_Periodic;
			break;
		case testMHD_vortex:
			x_0 = -10.0, x_1 = 10.0, y_0 = -10.0, y_1 = 10.0, z_0 = -1, z_1 = 1;
			time_stop = 1000.0;         // 平流绕域一圈
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;         // 平流版用周期(让环绕回来)
			y_bc = y_Periodic;
			break;
		case testHarrisHall:
			// Report Table 4.6, Test 12: Hall MHD Harris sheet (GEM challenge geometry).
			// Lx x Ly = 4pi x 2pi on 128 x 64, periodic in x, transmissive in y,
			// eta = 5e-3, d_i = 1, psi0 = 0.1, rho_floor = 0.02, t_end = 12.
			// The reference EXPLICIT calculation terminates near t = 12: X-line density
			// cavitation drives |JxB|/rho up while its Hall CFL uses a density floor of 0.1
			// against an evolved floor of 0.02, so the whistler speed is underestimated by
			// sqrt(0.1/0.02) ~ 2.2 and grid-scale modes outrun the explicit stability limit.
			// Here the Hall term is implicit and the accuracy cap in caldt uses the ACTUAL
			// density, with no floor mismatch, so surviving past t ~ 9-12 is the point of
			// running this case.
			x_0 = -2.0*M_PI, x_1 = 2.0*M_PI, y_0 = -M_PI, y_1 = M_PI, z_0 = -1, z_1 = 1;
			time_stop = 12.0;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;
			y_bc = y_Transmissive;
			// non-ideal switches for Table 4.6 Test 12
			add_resistive     = true;  eta_resistive = 5.0e-3;   // S ~ 2500
			add_hyperresistive = false; eta_hyper    = 0.0;      // not used in Test 12
			Hall_di = 1.0;
			break;
		case testHallReconnection:
			// Report Sec 7.2.4 first-peak survey configuration: Lx x Ly = 4pi x 2pi,
			// 256 x 128, dx = dy = 0.049087 = d_i/20, lambda = 0.5 = d_i/2, eta = 0,
			// eta_H = 6.25e-5, transmissive boundaries, duration t = 20.
			x_0 = -2.0*M_PI, x_1 = 2.0*M_PI, y_0 = -M_PI, y_1 = M_PI, z_0 = -1, z_1 = 1;
			time_stop = 20.0;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
	}
	
	dx = (x_1 - x_0)/(nxPoints-1);
	dy = (y_1 - y_0)/(nyPoints-1);
	dz = (z_1 - z_0)/(nzPoints-1);

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells), u_con_n(nxCells, nyCells), Qe(nxCells, nyCells),Qe_pri(nxCells, nyCells), Qs(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells),w_yL(nxCells, nyCells), w_yR(nxCells, nyCells),A(nxCells,nyCells);
	
	// 新增A的IMEX相关变量
	Grid A_n(nxCells, nyCells);           // 时间步开头的A
	Tensor HsA(RKs, nxCells, nyCells), HsB(RKs, nxCells, nyCells);    // A的时间导数历史，需要3分量
	Array Bx_p0(nxCells,nyCells), By_0p(nxCells,nyCells),Bx_p0_n(nxCells,nyCells), By_0p_n(nxCells,nyCells), Az(nxCells, nyCells), Bx(nxCells,nyCells),By(nxCells,nyCells), Bz(nxCells,nyCells);
	Array Bx_p0_Qe(nxCells,nyCells), By_0p_Qe(nxCells,nyCells),Bx_p0_Qs(nxCells,nyCells), By_0p_Qs(nxCells,nyCells);
	u_initialization(u_pri, Bx_p0, By_0p);
	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1); 
	upri2ucon(u_con, u_pri);
	A_initial(A);
	
	
	std::vector<double> bRK(RKs);
	Array Ae(RKs, RKs), Ai(RKs, RKs);
	Tensor Hs(RKs, nxCells, nyCells);
	Array shock_eta(nxCells, nyCells);
	std::vector<std::pair<int, double>> div_error;
	std::vector<std::pair<int, double>> B_inf;
	std::vector<std::pair<int, double>> dt_comp;
	const int plotINTERVAL = 100;
	std::vector<double> IntermediatePlotTimes = {
		0.25 * time_stop,
		0.5 * time_stop,
		0.75 * time_stop
	};
	bool isFinalTime = false;
	
	IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);

	if(RKs == 1){
		Ae(0,0) = 0.0;
		Ai(0,0) = 1.0;
		bRK[0]  = 1.0;
	}else if(RKs == 2){
		double bb = 1.0 - 1.0 / sqrt(2.0);
		double cc = 1.0 / (2.0 * bb);
		Ae(0,0) = 0.0;
		Ae(0,1) = 0.0;
		
		Ae(1,0) = cc;
		Ae(1,1) = 0.0;
		
		Ai(0,0) = bb;
		Ai(0,1) = 0.0;
		
		Ai(1,0) = 1 - bb;
		Ai(1,1) = bb;
		
		bRK[0] = 1.0 - bb;
		bRK[1] = bb;
	}else if(RKs == 3){
		double r12 = 0.5, r13=1.0/3.0, r14=1.0/4.0;
		Ae(0,0) = 0.0;
		Ae(0,1) = 0.0;
		Ae(0,2) = 0.0;
		
		Ae(1,0) = r12;
		Ae(1,1) = 0.0;
		Ae(1,2) = 0.0;
		
		Ae(2,0) = r12;
		Ae(2,1) = r12;
		Ae(2,2) = 0.0;
		
		Ai(0,0) = r14;
		Ai(0,1) = 0.0;
		Ai(0,2) = 0.0;
		
		Ai(1,0) = 0.0;
		Ai(1,1) = r14;
		Ai(1,2) = 0.0;
		
		Ai(2,0) = r13;
		Ai(2,1) = r13;
		Ai(2,2) = r13;
		
		bRK[0] = r13;
		bRK[1] = r13;
		bRK[2] = r13;
	}
	div_error.emplace_back(timestep, cal_divergence_error(u_con, Bx_p0, By_0p));
	double B_inf_m = 0;
	double divB = 0;
	for(int i = 0; i < nxCells + 2*ghost_cells; i++){
		for(int j = 0; j < nyCells + 2*ghost_cells; j++){
			divB = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx) + (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
			if(fabs(divB) > B_inf_m){
				B_inf_m = fabs(divB);
			}
		}
	}
	B_inf.emplace_back(timestep, B_inf_m);
	
	// Report Sec 7.2.4 duration; the diagnostics are emitted on the Delta t_out grid.
	if((test == testHallReconnection || test == testHarrisHall) && hall_diag_on) HallReconnectionDiagnostics(u_con, t);
	double next_hall_diag = hall_diag_dt;

	do{
		timestep++;
		g_gmres_matvec = 0;
		g_bc_prim = true;  bound_cond_u(u_pri);  g_bc_prim = false;  // primitive ghosts (slot 4 = p)
		calch(u_pri);
		calcp(u_pri);
		calca(u_pri);
		
		dt = fmin(caldt(u_pri), time_stop - t);
		if(g_fixed_dt > 0.0) dt = fmin(g_fixed_dt, time_stop - t);   // 对比用:强制固定 dt(=CASTOR 的 FIXDT)
		if(fabs(dt - (time_stop - t)) < 1e-8){
			isFinalTime = true;
		}
		
		Grid A_n_snapshot(nxCells, nyCells);
		copy(A_n_snapshot, A);
		
		
		MagneticEigenvalues(u_pri);
		EnergyEigenvalues(u_pri);
		t = t + dt;
		//std::cout << B_max << " " << U_max << endl;
		
		upri2ucon(u_con, u_pri);
		
		copy(u_con_n,u_con);
		
		for(int i = 0; i < nxCells + 2*ghost_cells; i++){
			for(int j = 0; j < nyCells + 2*ghost_cells; j++){
				Bx_p0_n(i,j) = Bx_p0(i,j); 
				By_0p_n(i,j) = By_0p(i,j);
			}
		}
		
		for(int rks = 0; rks < RKs; rks++){
			copy(Qe,u_con_n);
			copy(Qs,u_con_n);
			
			Grid A_s(nxCells, nyCells), A_e(nxCells, nyCells);
			copy(A_s, A_n_snapshot);
			copy(A_e, A_n_snapshot);
			
			for(int i = 0; i < nxCells + 2*ghost_cells; i++){
				for(int j = 0; j < nyCells + 2*ghost_cells; j++){
					Bx_p0_Qe(i,j) = Bx_p0_n(i,j); 
					By_0p_Qe(i,j) = By_0p_n(i,j);
					Bx_p0_Qs(i,j) = Bx_p0_n(i,j);
					By_0p_Qs(i,j) = By_0p_n(i,j);
				}
			}
						
			ucon2upri(Qe_pri, Qe);
			
			for(int r = 0; r < rks; r++){
				for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
					for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
						for(int k = 0; k < NUM_VARS; k++){
							Qe(i,j,k) = Qe(i,j,k) + Ae(rks,r) * dt * Hs(r,i,j,k);
							Qs(i,j,k) = Qs(i,j,k) + Ai(rks,r) * dt * Hs(r,i,j,k);
						}
					}
				}
			}
			for(int r = 0; r < rks; r++){
				for(int i = 0; i < nxCells + 2*ghost_cells; i++){
					for(int j = 0; j < nyCells + 2*ghost_cells; j++){
						for(int k = 0; k < 3; k++){   // A 只有 3 分量
							A_e(i,j,k) += Ae(rks,r) * dt * HsA(r,i,j,k);
							A_s(i,j,k) += Ai(rks,r) * dt * HsA(r,i,j,k);
						}
						Bx_p0_Qe(i,j) += Ae(rks,r) * dt * HsB(r,i,j,0); 
						By_0p_Qe(i,j) += Ae(rks,r) * dt * HsB(r,i,j,1);
						Bx_p0_Qs(i,j) += Ai(rks,r) * dt * HsB(r,i,j,0);
						By_0p_Qs(i,j) += Ai(rks,r) * dt * HsB(r,i,j,1);
					}
				}
			}
			
						
			copy(A, A_s);
			
			adt =  Ai(rks,rks) * dt;
			
			
			TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR);
			
			// Convective Subsystem F^c
			RusanovC_Update(u_con, Qs, w_xL,w_xR,w_yL,w_yR, adt);
			bound_cond_u(u_con);

			// Magnetic Subsystem F^B
			if(nxCells == 1 || nyCells == 1){
				MagSubsystem_UpdateB(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A, adt,rks);
			}else{
				switch(magfield){
					case StaggerCT: case UnStaggerCT: case NoDivFree:
						MagSubsystem_UpdateB(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A, adt, rks);
						break;
					case PotentialA:
						MagSubsystem_UpdateA(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A,A_s, adt, rks);
						break;
				}
			}
			bound_cond_u(u_con);
			
			// Store Magnetic Field
			for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
					Bx(i,j) = u_con(i,j,5);
					By(i,j) = u_con(i,j,6);
					Bz(i,j) = u_con(i,j,7);
				}
			}
			
			Grid u_test(nxCells, nyCells);
			for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
					w_xL(i,j,5) = u_con(i,j,5);
					w_xL(i,j,6) = u_con(i,j,6);
					w_xL(i,j,7) = u_con(i,j,7);
					
					w_xR(i,j,5) = u_con(i,j,5);
					w_xR(i,j,6) = u_con(i,j,6);
					w_xR(i,j,7) = u_con(i,j,7);
					
					w_yL(i,j,5) = u_con(i,j,5);
					w_yL(i,j,6) = u_con(i,j,6);
					w_yL(i,j,7) = u_con(i,j,7);
					
					w_yR(i,j,5) = u_con(i,j,5);
					w_yR(i,j,6) = u_con(i,j,6);
					w_yR(i,j,7) = u_con(i,j,7);
					
				}
			}
			
			switch(formula){
				case EnergyFormula:
					Energy_flux_Update(u_con, Qs, w_xL, w_xR, w_yL, w_yR, adt);
					break;
				case PressureFormula:
					Pressure_flux_Update(u_con, Qs, w_xL, w_xR, w_yL, w_yR, adt);
					break;
			}
			bound_cond_u(u_con); 
			
			// Pressure Subsystem F^p + Momentum Update
			switch(formula){
				case EnergyFormula:
					PreSubsystem_Energy(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, adt,rks);
					break;
				case PressureFormula:
					PreSubsystem_Pressure(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, adt,rks);
					break;
			}
			bound_cond_u(u_con); 
			
			if(add_hall && !hall_implicit) HallMHD_SourceRelax(u_con, adt);
			
			// divergence-free
			if(nxCells != 1 && nyCells != 1){
				switch(magfield){
					case StaggerCT:
						Staggered_CT(u_con, Bx_p0, By_0p, Bx_p0_Qs,By_0p_Qs,adt);
						break;
					case UnStaggerCT:
						Unstaggered_CT(u_con, A, A_s,adt);
						break;
					case PotentialA: case NoDivFree:
						break;
				}
			}
			
			
			bound_cond_u(u_con);
			
			// Modify Energy
			if(Energy_Modify){
				for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
					for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
						u_con(i,j,4) = u_con(i,j,4) + ((u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7)) - (Bx(i,j) * Bx(i,j) + By(i,j) * By(i,j) + Bz(i,j) * Bz(i,j))) / 2.0;
					}
				}
			}
			
			bound_cond_u(u_con);
			ucon2upri(u_pri, u_con);
			
			for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
				for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						Hs(rks,i,j,k) = (u_con(i,j,k) - Qs(i,j,k)) / adt;
					}
				}
			}
			
			for(int i = 0; i < nxCells + 2*ghost_cells; i++){
				for(int j = 0; j < nyCells + 2*ghost_cells; j++){
					for(int k = 0; k < 3; k++){
						HsA(rks,i,j,k) = (A(i,j,k) - A_s(i,j,k)) / adt;
					}
					HsB(rks,i,j,0) = (Bx_p0(i,j) - Bx_p0_Qs(i,j)) / adt;
					HsB(rks,i,j,1) = (By_0p(i,j) - By_0p_Qs(i,j)) / adt;
				}
			}
			
				
			bound_cond_u(u_con);
		}
				
		ucon2upri(u_pri, u_con);
		
		std::cout<<"[timestep, dt, time] = "  << timestep << " "<<dt << " " << t << ", Divergence Error: "<< cal_divergence_error(u_con, Bx_p0, By_0p)<< endl;
		{   // ===== 跨代码对比诊断(归一化无关量) =====
			double emin=1e30, mass=0, intE=0, mmax=0; int im=-1, jm=-1;
			for(int i=ghost_cells;i<nxCells+ghost_cells;i++) for(int j=ghost_cells;j<nyCells+ghost_cells;j++){
				double rho=u_con(i,j,0);
				double k=0.5*(u_con(i,j,1)*u_con(i,j,1)+u_con(i,j,2)*u_con(i,j,2)+u_con(i,j,3)*u_con(i,j,3))/rho;
				double m=0.5*(u_con(i,j,5)*u_con(i,j,5)+u_con(i,j,6)*u_con(i,j,6)+u_con(i,j,7)*u_con(i,j,7));
				double e=u_con(i,j,4)-k-m;
				mass+=rho; intE+=e; if(m>mmax)mmax=m; if(e<emin){emin=e;im=i;jm=j;}
			}
			double xe=x_0+(im-ghost_cells+0.5)*dx, ye=y_0+(jm-ghost_cells+0.5)*dy;
			std::cout<<"DIAG ts="<<timestep<<" t="<<t<<" emin="<<emin<<" at("<<xe<<","<<ye<<")"
			<<" mass="<<mass*dx*dy<<" intE="<<intE*dx*dy<<" magE_max="<<mmax
			<<" gmresMV="<<g_gmres_matvec
			<<" dt_ideal="<<g_dt_ideal<<" dt_hall="<<g_dt_hall<<std::endl;
		}
		if((test == testHallReconnection || test == testHarrisHall) && hall_diag_on
		   && t >= next_hall_diag - 1e-12){
			HallReconnectionDiagnostics(u_con, t);
			while(next_hall_diag <= t + 1e-12) next_hall_diag += hall_diag_dt;
		}
		IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);
		div_error.emplace_back(timestep, cal_divergence_error(u_con, Bx_p0, By_0p));
		
		
		double B_inf_m = 0;
		double divB = 0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				divB = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx) + (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
				if(fabs(divB) > B_inf_m){
					B_inf_m = fabs(divB);
				}
			}
		}
		B_inf.emplace_back(timestep, B_inf_m);
		dt_comp.emplace_back(timestep, imex_dt/fully_dt);
		
	}while(t < time_stop);
	
	ofstream outFilerho(Path + "/RK_rho.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerho << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,0) << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx(Path + "/RK_vx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,1) << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy(Path + "/RK_vy.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,2) << endl;
		}
	}
	outFilevy.close();
	
	ofstream outFileE(Path + "/RK_E.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileE << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,4) << endl;
		}
	}
	outFileE.close();
	
	
	ofstream outFilep(Path + "/RK_p.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilep << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,4)) << endl;
		}
	}
	outFilep.close();
	
	ofstream outFileBx(Path + "/RK_Bx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,5)) << endl;
		}
	}
	outFileBx.close();

	ofstream outFileBy(Path + "/RK_By.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,6)) << endl;
		}
	}
	outFileBy.close();
	
	ofstream outFileBz(Path + "/RK_Bz.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBz << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,7)) << endl;
		}
	}
	outFileBz.close();
	
	
	ofstream outFilev(Path + "/RK_v_mag.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilev << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2)) << endl;
		}
	}
	outFilev.close();
	
	
	ofstream outFiledt_comp(Path + "/RK_dt_comp.dat");
	for(const auto& p : dt_comp)
		outFiledt_comp << p.first << " " << p.second << "\n";
	outFiledt_comp.close();
	
	ofstream outFilediv(Path + "/RK_div_error.dat");
	for(const auto& p : div_error)
		outFilediv << p.first << " " << p.second << "\n";
	outFilediv.close();
	
	ofstream outFileB_inf(Path + "/RK_B_inf.dat");
	for(const auto& p : B_inf)
		outFileB_inf << p.first << " " << p.second << "\n";
	outFileB_inf.close();
}