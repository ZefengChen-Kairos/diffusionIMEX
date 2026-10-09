#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
// boost removed
#include <functional>
#include <iomanip>
#include <random>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Core>



#define NUM_VARS 7
enum Formula{EnergyFormula, PressureFormula};
enum Test{testshockbubble,test2dahshockbubble,testDodecaneShockTube};
enum Reconstruction{First_Order,Second_Order};
enum shock_threshold_Left{one, val_left};
enum shock_threshold_Right{zero, val_right};
enum shock_detector{Yes_add, No_First, No_Second};
enum SolvePressure{no_correction, once_correction, all_correction};
enum dtCalculation{Fullydt, IMEXdt};
enum x_boundary_condition{x_Transmissive, x_Reflective, x_Periodic};
enum y_boundary_condition{y_Transmissive, y_Reflective, y_Periodic};

dtCalculation caltimestep = IMEXdt;
Formula formula = PressureFormula;
x_boundary_condition x_bc = x_Transmissive;
y_boundary_condition y_bc = y_Transmissive;
SolvePressure PressureSolver = all_correction;
Test test = testDodecaneShockTube;

const int nxCells = 1000;
const int nyCells = 1;
const int nzCells = 1;
const int ghost_cells = 3;

const int RKs = 2;
double CFL = 0.1;
Reconstruction recon = Second_Order;

bool add_Evisc = true;    // shock-activated chi-Laplacian on the mixture pressure (eq 53-54)

const double ALPHA_MIN = 1e-8;
const double mu = 1e60;
double GMREStol = 1e-8;
double GMRESiter = 2000;
int GMRESrest = 400;
bool add_restart = false;

//pressure formula
shock_detector add_shock_detector = No_Second;
double shock_kappa1 = 0.001;
double shock_threshold = 0.9;
shock_threshold_Left threshold_Left = val_left;
shock_threshold_Right threshold_Right = zero;


int nxPoints = nxCells + 1, nyPoints = nyCells + 1, nzPoints = nzCells + 1;
double t, time_stop, dt, adt, adt2;
int timestep = 0;
double x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1, B_max_x = 0, B_max_y = 0, U_max_x = 0,U_max_y = 0,a_max_x = 0,a_max_y = 0, vx_max, vy_max;
//double t, dt, time_stop = 0.2,x_0 = -1.0,y_0 = -1.0,x_1 = 1.0, y_1 = 1.0, CFL= 0.8;

double dx = (x_1 - x_0)/(nxPoints-1);
double dy = (y_1 - y_0)/(nyPoints-1);
double dz = (z_1 - z_0)/(nzPoints-1);
double ch_x,ch_y, cp_x, cp_y, ca_x, ca_y;
Eigen::VectorXd g_chi;   // Balsara flattener shock indicator (eq 52), per interior cell



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

class EquationOfState {
public:
	virtual ~EquationOfState() = default;
	virtual double getGamma() const = 0;
	virtual double getPInf() const = 0;
	virtual double getP_fromE(double E, double rho, double u, double v) const = 0;
	virtual double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,
		double u, double v, double alpha) const = 0;
	virtual double getAlphaE_fromAlphaP(double alphaP, double alphaRho,
		double u, double v, double alpha) const = 0;
	// 仅保留带 rho 的版本
	virtual double getAlphaP(double alphaRhoe, double rho, double alpha) const = 0;
	virtual double internalEnergy(double rho, double p) const = 0;
	virtual double pressure(double rho, double e) const = 0;
	virtual double calcs(double rho, double p) const = 0;
	virtual void e_linear_coeff(double rho, double& c, double& d) const = 0;
	virtual EquationOfState* clone() const = 0;
};

class MieGruneisenEOS : public EquationOfState {
public:
	// 通用压力公式
	double pressure(double rho, double e) const override {
		return p_ref(rho) + Gamma(rho) * rho * (e - e_ref(rho));
	}
	// 通用内能公式
	double internalEnergy(double rho, double p) const override {
		return e_ref(rho) + (p - p_ref(rho)) / (Gamma(rho) * rho);
	}
	// 线性化系数
	void e_linear_coeff(double rho, double& c, double& d) const override {
		double Gamma_val = Gamma(rho);
		c = 1.0 / (Gamma_val * rho);
		d = e_ref(rho) - p_ref(rho) / (Gamma_val * rho);
	}
	// 通用 getAlphaP 实现（利用 pressure 函数）
	double getAlphaP(double alphaRhoe, double rho, double alpha) const override {
		double e = alphaRhoe / (alpha * rho);
		double p = pressure(rho, e);
		return alpha * p;
	}
	
	// 以下纯虚函数由派生类实现
	virtual double p_ref(double rho) const = 0;
	virtual double e_ref(double rho) const = 0;
	virtual double Gamma(double rho) const = 0;
	virtual double calcs(double rho, double p) const override = 0;
};

class StiffenedGasEOS : public MieGruneisenEOS {
	double gamma, p_inf;
public:
	StiffenedGasEOS(double gamma_, double p_inf_) : gamma(gamma_), p_inf(p_inf_) {}
	
	double getGamma() const override { return gamma; }
	double getPInf() const override { return p_inf; }
	double getP_fromE(double E, double rho, double u, double v) const override {
		double ke = 0.5*rho*(u*u+v*v);
		return (E-ke)*(gamma-1) - gamma*p_inf;
	}
	double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,
		double u, double v, double alpha) const override {
			double ke = 0.5*alphaRho*(u*u+v*v);
			double alphaRhoe = alphaRhoE - ke;
			return alphaRhoe*(gamma-1) - alpha*gamma*p_inf;
		}
	double getAlphaE_fromAlphaP(double alphaP, double alphaRho,
		double u, double v, double alpha) const override {
			double alphaRhoe = (alphaP + alpha*gamma*p_inf)/(gamma-1);
			double ke = 0.5*alphaRho*(u*u+v*v);
			return alphaRhoe + ke;
		}
	// 实现 getAlphaP 带 rho 参数（虽然 rho 未使用）
	double getAlphaP(double alphaRhoe, double rho, double alpha) const override {
		return alphaRhoe*(gamma-1) - alpha*gamma*p_inf;
	}
	// Mie-Grüneisen 接口
	double p_ref(double /*rho*/) const override { return -gamma * p_inf; }
	double e_ref(double /*rho*/) const override { return 0.0; }
	double Gamma(double /*rho*/) const override { return gamma - 1.0; }
	double calcs(double rho, double p) const override {
		// 刚性 EOS 重构过冲可能给出 p+pinf<0 或 rho<=0 -> sqrt(负)=NaN.
		// floor 到小正值, 保证声速有限非负 (正性保护).
		double num = fmax(gamma * (p + p_inf), 1e-12);
		double den = fmax(rho, 1e-12);
		return std::sqrt(num / den);
	}
	EquationOfState* clone() const override { return new StiffenedGasEOS(*this); }
};

class IdealGasEOS : public StiffenedGasEOS {
public:
	IdealGasEOS(double gamma_) : StiffenedGasEOS(gamma_, 0.0) {}
	// 所有函数自动正确，除了可能需要重写 getPInf 返回 0（基类已处理）
};

class JWLEOS : public MieGruneisenEOS {
	double A, B, R1, R2, omega, rho0;
public:
	JWLEOS(double A_, double B_, double R1_, double R2_, double omega_, double rho0_)
	: A(A_), B(B_), R1(R1_), R2(R2_), omega(omega_), rho0(rho0_) {}
	
	double getGamma() const override { return 1.0 + omega; }
	double getPInf() const override { return 0.0; }
	double getP_fromE(double E, double rho, double u, double v) const override {
		double ke = 0.5*rho*(u*u+v*v);
		return pressure(rho, (E-ke)/rho);
	}
	double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,
		double u, double v, double alpha) const override {
			double ke = 0.5*alphaRho*(u*u+v*v);
			double alphaRhoe = alphaRhoE - ke;
			double rho = alphaRho / alpha;
			double e = alphaRhoe / alphaRho;
			double p = pressure(rho, e);
			return alpha * p;
		}
	double getAlphaE_fromAlphaP(double alphaP, double alphaRho,
		double u, double v, double alpha) const override {
			double p = alphaP / alpha;
			double rho = alphaRho / alpha;
			double e = internalEnergy(rho, p);
			double alphaRhoe = alphaRho * e;
			double ke = 0.5*alphaRho*(u*u+v*v);
			return alphaRhoe + ke;
		}
	// 实现 getAlphaP 需要用到 rho
	double getAlphaP(double alphaRhoe, double rho, double alpha) const override {
		double e = alphaRhoe / (alpha * rho);
		double p = pressure(rho, e);
		return alpha * p;
	}
	// Mie-Grüneisen 函数
	double p_ref(double rho) const override {
		double v = rho0 / rho;
		return A * std::exp(-R1*v) + B * std::exp(-R2*v);
	}
	double e_ref(double rho) const override {
		double v = rho0 / rho;
		return A/(R1*rho0)*std::exp(-R1*v) + B/(R2*rho0)*std::exp(-R2*v);
	}
	double Gamma(double) const override { return omega; }
	double calcs(double rho, double p) const override {
		double v = rho0 / rho;
		double pref = p_ref(rho);
		double dp_ref_drho = (rho0/(rho*rho)) * (A*R1*std::exp(-R1*v) + B*R2*std::exp(-R2*v));
		double c2 = dp_ref_drho + (1.0 + omega)/rho * (p - pref);   // MG: c² = p_ref' + (1+Γ)(p-p_ref)/ρ
		return std::sqrt(std::max(c2, 0.0));
	}
	EquationOfState* clone() const override { return new JWLEOS(*this); }
};

class CochranChanEOS : public EquationOfState {
private:
	const double A;      // Pa
	const double B;      // Pa
	const double E1;     // 无量纲
	const double E2;     // 无量纲
	const double Gamma0; // 无量纲
	const double rho0;   // 参考密度 (kg/m³)
	
public:
	CochranChanEOS(double A_, double B_, double E1_, double E2_, double Gamma0_, double rho0_)
	: A(A_), B(B_), E1(E1_), E2(E2_), Gamma0(Gamma0_), rho0(rho0_) {}
	
	// 以下实现所有纯虚函数
	
	double getGamma() const override {
		// 返回有效比热比（近似值，用于需要 γ 的场合）
		return 1.0 + Gamma0;
	}
	
	double getPInf() const override {
		return 0.0;
	}
	
	double getP_fromE(double E, double rho, double u, double v) const override {
		double ke = 0.5 * rho * (u*u + v*v);
		double e = (E - ke) / rho;
		return pressure(rho, e);
	}
	
	double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,
		double u, double v, double alpha) const override {
			double ke = 0.5 * alphaRho * (u*u + v*v);
			double alphaRhoe = alphaRhoE - ke;
			double rho = alphaRho / alpha;          // 相密度
			double e = alphaRhoe / alphaRho;        // 比内能
			double p = pressure(rho, e);
			return alpha * p;
		}
	
	double getAlphaE_fromAlphaP(double alphaP, double alphaRho,
		double u, double v, double alpha) const override {
			double p = alphaP / alpha;
			double rho = alphaRho / alpha;
			double e = internalEnergy(rho, p);
			double alphaRhoe = alphaRho * e;
			double ke = 0.5 * alphaRho * (u*u + v*v);
			return alphaRhoe + ke;
		}
	
	// 核心函数：由相内能密度、相密度、体积分数直接求相压力密度
	double getAlphaP(double alphaRhoe, double rho, double alpha) const override {
		double e = alphaRhoe / (alpha * rho);
		double p = pressure(rho, e);
		return alpha * p;
	}
	
	// 比内能 e(ρ, p)
	double internalEnergy(double rho, double p) const override {
		double eta = rho / rho0;
		// 计算参考压力
		double pref = A * std::pow(eta, E1) - B * std::pow(eta, E2);
		// 计算参考内能（通过积分 p_ref/ρ^2 dρ 得到）
		auto term = [&](double coeff, double expo) -> double {
			if (std::abs(expo - 1.0) < 1e-12) {
				return -coeff / rho0 * std::log(eta);
			} else {
				return -coeff / rho0 * (std::pow(eta, expo - 1.0) - 1.0) / (1.0 - expo);
			}
		};
		double eref = term(A, E1) + term(B, E2);
		return eref + (p - pref) / (Gamma0 * rho);
	}
	
	// 压力 p(ρ, e)
	double pressure(double rho, double e) const override {
		double eta = rho / rho0;
		double pref = A * std::pow(eta, E1) - B * std::pow(eta, E2);
		auto term = [&](double coeff, double expo) -> double {
			if (std::abs(expo - 1.0) < 1e-12) {
				return -coeff / rho0 * std::log(eta);
			} else {
				return -coeff / rho0 * (std::pow(eta, expo - 1.0) - 1.0) / (1.0 - expo);
			}
		};
		double eref = term(A, E1) + term(B, E2);
		return pref + Gamma0 * rho * (e - eref);
	}
	
	// 声速 c = sqrt( (∂p/∂ρ)_s )
	double calcs(double rho, double p) const override {
		double eta = rho / rho0;
		double pref = A * std::pow(eta, E1) - B * std::pow(eta, E2);
		double dp_ref_drho = (A * E1 / rho0) * std::pow(eta, E1 - 1.0)
		- (B * E2 / rho0) * std::pow(eta, E2 - 1.0);
		double c2 = dp_ref_drho + (1.0 + Gamma0)/rho * (p - pref);   // MG: c² = p_ref' + (1+Γ)(p-p_ref)/ρ
		return std::sqrt(std::max(c2, 0.0));
	}
	
	// 线性化系数 e = c * p + d （固定 ρ）
	void e_linear_coeff(double rho, double& c, double& d) const override {
		double eta = rho / rho0;
		double pref = A * std::pow(eta, E1) - B * std::pow(eta, E2);
		auto term = [&](double coeff, double expo) -> double {
			if (std::abs(expo - 1.0) < 1e-12) {
				return -coeff / rho0 * std::log(eta);
			} else {
				return -coeff / rho0 * (std::pow(eta, expo - 1.0) - 1.0) / (1.0 - expo);
			}
		};
		double eref = term(A, E1) + term(B, E2);
		c = 1.0 / (Gamma0 * rho);
		d = eref - pref / (Gamma0 * rho);
	}
	
	EquationOfState* clone() const override {
		return new CochranChanEOS(*this);
	}
};

using namespace std;
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void g_convective(double *g, double *con);
void f_full_con(double *f, double *con,const EquationOfState& eos1,const EquationOfState& eos2);
void g_full_con(double *g, double *con,const EquationOfState& eos1,const EquationOfState& eos2);
void ini_con2pri(double *pri, double *con, const EquationOfState& eos1,const EquationOfState& eos2);
void ini_pri2con(double *con, double *pri,const EquationOfState& eos1, const EquationOfState& eos2);
void upri2ucon(Grid &u_con, Grid &u_pri, const EquationOfState& eos1,const EquationOfState& eos2);
void ucon2upri(Grid &u_pri,Grid &u_con, const EquationOfState& eos1,const EquationOfState& eos2);
void RusanovC_x(double *QL, double *QR, double *flux, const EquationOfState& eos1,const EquationOfState& eos2);
void RusanovC_y(double *QL, double *QR, double *flux, const EquationOfState& eos1,const EquationOfState& eos2);
void RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR, const EquationOfState& eos1,const EquationOfState& eos2);
void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p);
void bound_cond_u(Grid &u);
double caldt(Grid &u_pri);
void bound_cond_array(Array &u);
double calenthalpy(double *con, const EquationOfState& eos1,const EquationOfState& eos2);
void ShockDetector(Array &shock_eta, Grid &u_con, const EquationOfState& eos1,const EquationOfState& eos2);
int M_P(int i, int j);
void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,Array& chi,double h, const EquationOfState& eos1, const EquationOfState& eos2);
void GMRES_Solver(Eigen::VectorXd& Ax, Eigen::VectorXd&x,Eigen::VectorXd& b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&,double, const EquationOfState&, const EquationOfState&)> matrix_vector_product, Grid& u_con, Grid& Qe, Grid& Qs, Array& chi,double h, const EquationOfState& eos1, const EquationOfState& eos2, const int dim);
void PressureSubsystem(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, double h, const int rks, const EquationOfState& eos1, const EquationOfState& eos2);
double calkinetic(double *con);



// ===== Battisti-Boscheri (12d)/(23): 隐式中心 fluctuation 形式解 α^{n+1} =====
// α_i + Δt/(2Δx)[u_{i+1/2}(α_{i+1}-α_i) - u_{i-1/2}(α_{i-1}-α_i)] + (y) = α^n + Δt·D[α^n]
// D = 显式 Rusanov 型耗散 (36)-(37), ν± = ½(σ_i+σ_{i±1}), σ = ½Δx·aΩ
void copy(Grid &state_1, Grid &state_0){
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				state_1(i,j,k) = state_0(i,j,k);
			}
		}
	}
}

void ini_con2pri(double *pri, double *con, const EquationOfState& eos1,const EquationOfState& eos2){
	// convert conservative (rho rhov E) to primitive (rho v p)
	// rho v p
	double alpha1 = con[0];
	double alpha2 = 1 - alpha1;
	double rho1 = con[1] / alpha1;
	double rho2 = con[2] / alpha2;
	double rho = con[1] + con[2];
	double vx = con[3] / rho;
	double vy = con[4] / rho;
	
	double alpha1p1 = eos1.getAlphaP_fromAlphaE(con[5], con[1], vx, vy, alpha1);
	double alpha2p2 = eos2.getAlphaP_fromAlphaE(con[6], con[2], vx, vy, alpha2);

	
	double c1, c2, d1, d2;
	
	pri[0] = con[0];
	pri[1] = con[1] / alpha1;
	pri[2] = con[2] / alpha2;
	pri[3] = con[3] / rho;
	pri[4] = con[4] / rho;
	pri[5] = alpha1p1;
	pri[6] = alpha2p2;
}

void ini_pri2con(double *con, double *pri,const EquationOfState& eos1, const EquationOfState& eos2){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	// rho rhov E
	double alpha1 = pri[0];
	double alpha2 = 1 - alpha1;
	double rho1 = pri[1];
	double rho2 = pri[2];
	double rho = alpha1 * rho1 + alpha2 * rho2;
	double vx = pri[3];
	double vy = pri[4];

	con[0] = pri[0];
	con[1] = alpha1 * rho1;
	con[2] = alpha2 * rho2;
	con[3] = rho * pri[3];
	con[4] = rho * pri[4];
	con[5] = eos1.getAlphaE_fromAlphaP(pri[5], con[1], vx, vy, alpha1);
	con[6] = eos2.getAlphaE_fromAlphaP(pri[6], con[2], vx, vy, alpha2);
}

void upri2ucon(Grid &u_con, Grid &u_pri,const EquationOfState& eos1,const EquationOfState& eos2){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_pri[k] = u_pri(i,j,k);
			}
			ini_pri2con(ini_u_con, ini_u_pri,eos1, eos2);
			for(int k = 0; k < NUM_VARS; k++){
				u_con(i, j, k) = ini_u_con[k];
			}
		}
	}
}

void ucon2upri(Grid &u_pri,Grid &u_con,const EquationOfState& eos1,const EquationOfState& eos2){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_con[k] = u_con(i, j, k);
			}
			ini_con2pri(ini_u_pri, ini_u_con,eos1,eos2);
			for(int k = 0; k < NUM_VARS; k++){
				u_pri(i, j, k) = ini_u_pri[k];
			}
		}
	}
}

static void buildA_advX(const double* Q, double A[NUM_VARS][NUM_VARS]){
	double rho = Q[1] + Q[2];
	double u   = Q[3]/rho;
	double v   = Q[4]/rho;
	double Y1  = Q[1]/rho, Y2 = Q[2]/rho;
	double K   = 0.5*(u*u + v*v);                      // 比动能
	for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) A[r][c]=0.0;
	A[0][0]=u;                                                            //
	A[1][1]=u*(1.0-Y1); A[1][2]=-u*Y1;      A[1][3]=Y1;                   //
	A[2][1]=-u*Y2;      A[2][2]=u*(1.0-Y2); A[2][3]=Y2;                   //
	A[3][1]=-u*u;       A[3][2]=-u*u;       A[3][3]=2.0*u;                //
	A[4][1]=-u*v;       A[4][2]=-u*v;       A[4][3]=v;  A[4][4]=u;        //
	A[5][1]=K*u*(1.0-3.0*Y1); A[5][2]=-3.0*Y1*K*u; A[5][3]=Y1*(u*u+K);
	A[5][4]=Y1*u*v; // α1E1
	A[6][1]=-3.0*Y2*K*u;      A[6][2]=K*u*(1.0-3.0*Y2); A[6][3]=Y2*(u*u+K);
	A[6][4]=Y2*u*v; // α2E2
	// A[5][5]=A[6][6]=0 : 动能通量不依赖能量变量本身
}

// ---- y 方向拟线性矩阵 A_y（u↔v 对称）--------------------------------------
static void buildA_advY(const double* Q, double A[NUM_VARS][NUM_VARS]){
	double rho = Q[1] + Q[2];
	double u   = Q[3]/rho;
	double v   = Q[4]/rho;
	double Y1  = Q[1]/rho, Y2 = Q[2]/rho;
	double K   = 0.5*(u*u + v*v);
	for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) A[r][c]=0.0;
	A[0][0]=v;                                                            //
	A[1][1]=v*(1.0-Y1); A[1][2]=-v*Y1;      A[1][4]=Y1;                   //
	A[2][1]=-v*Y2;      A[2][2]=v*(1.0-Y2); A[2][4]=Y2;                   //
	A[3][1]=-u*v;       A[3][2]=-u*v;       A[3][3]=v;  A[3][4]=u;        //
	A[4][1]=-v*v;       A[4][2]=-v*v;       A[4][4]=2.0*v;                //
	A[5][1]=K*v*(1.0-3.0*Y1); A[5][2]=-3.0*Y1*K*v; A[5][3]=Y1*u*v;
	A[5][4]=Y1*(v*v+K); // α1E1
	A[6][1]=-3.0*Y2*K*v;      A[6][2]=K*v*(1.0-3.0*Y2); A[6][3]=Y2*u*v;
	A[6][4]=Y2*(v*v+K); // α2E2
}

// ---- path 平均（3 点 Gauss–Legendre）-------------------------------------
static void pathAvgAx(const double* QL, const double* QR, double
	Ahat[NUM_VARS][NUM_VARS]){
		const double sq = sqrt(15.0)/10.0;
		const double s[3] = { 0.5 - sq, 0.5, 0.5 + sq };
		const double w[3] = { 5.0/18.0, 8.0/18.0, 5.0/18.0 };
		for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) Ahat[r][c]=0.0;
		for(int g=0; g<3; ++g){
			double Qs[NUM_VARS];
			for(int m=0;m<NUM_VARS;++m) Qs[m] = QL[m] + s[g]*(QR[m]-QL[m]);
			double Ag[NUM_VARS][NUM_VARS]; buildA_advX(Qs, Ag);
			for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) Ahat[r][c] +=
			w[g]*Ag[r][c];
		}
	}

static void pathAvgAy(const double* QL, const double* QR, double
	Ahat[NUM_VARS][NUM_VARS]){
		const double sq = sqrt(15.0)/10.0;
		const double s[3] = { 0.5 - sq, 0.5, 0.5 + sq };
		const double w[3] = { 5.0/18.0, 8.0/18.0, 5.0/18.0 };
		for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) Ahat[r][c]=0.0;
		for(int g=0; g<3; ++g){
			double Qs[NUM_VARS];
			for(int m=0;m<NUM_VARS;++m) Qs[m] = QL[m] + s[g]*(QR[m]-QL[m]);
			double Ag[NUM_VARS][NUM_VARS]; buildA_advY(Qs, Ag);
			for(int r=0;r<NUM_VARS;++r) for(int c=0;c<NUM_VARS;++c) Ahat[r][c] +=
			w[g]*Ag[r][c];
		}
	}

// ---- 路径平均速度 = Ahat[0][0] (α 行唯一的非零元) --------------------
static double pathAvgVel(const double* QL, const double* QR, int comp){
	const double sq = sqrt(15.0)/10.0;
	const double s[3] = { 0.5 - sq, 0.5, 0.5 + sq };
	const double w[3] = { 5.0/18.0, 8.0/18.0, 5.0/18.0 };
	double a = 0.0;
	for(int g=0; g<3; ++g){
		double m1 = QL[1] + s[g]*(QR[1]-QL[1]);
		double m2 = QL[2] + s[g]*(QR[2]-QL[2]);
		double mv = QL[comp] + s[g]*(QR[comp]-QL[comp]);
		a += w[g]*mv/(m1+m2);
	}
	return a;
}

// ---- 界面 fluctuation D^- / D^+（Rusanov: S = 物质速度）-------------------
void Path_RusanovC_x(double *QL, double *QR, double *Dm, double *Dp,
	const EquationOfState& eos1, const EquationOfState& eos2){
		double Ah[NUM_VARS][NUM_VARS];
		pathAvgAx(QL, QR, Ah);
		double dQ[NUM_VARS];
		for(int k = 0; k < NUM_VARS; k++) dQ[k] = QR[k] - QL[k];
		double rhoL = QL[1] + QL[2], rhoR = QR[1] + QR[2];
		double vxL = QL[3]/rhoL,     vxR = QR[3]/rhoR;
		double S = fmax(fabs(vxL), fabs(vxR));
		S = vx_max;
		// >>> 能量行耗散只取动能跳 (Battisti-Boscheri §3.3), 不含内能
		double kjmp[NUM_VARS];
		for(int k = 0; k < NUM_VARS; k++) kjmp[k] = dQ[k];
		{
			double keL = 0.5*(QL[3]*QL[3]+QL[4]*QL[4])/(rhoL*rhoL);
			double keR = 0.5*(QR[3]*QR[3]+QR[4]*QR[4])/(rhoR*rhoR);
			kjmp[5] = QR[1]*keR - QL[1]*keL;
			kjmp[6] = QR[2]*keR - QL[2]*keL;
		}
		for(int r = 0; r < NUM_VARS; r++){
			double Av = 0.0;
			for(int c = 0; c < NUM_VARS; c++) Av += Ah[r][c]*dQ[c];
			Dm[r] = 0.5*Av - 0.5*S*kjmp[r];
			Dp[r] = 0.5*Av + 0.5*S*kjmp[r];
		}
	}

void Path_RusanovC_y(double *QL, double *QR, double *Dm, double *Dp,
	const EquationOfState& eos1, const EquationOfState& eos2){
		double Ah[NUM_VARS][NUM_VARS];
		pathAvgAy(QL, QR, Ah);
		double dQ[NUM_VARS];
		for(int k = 0; k < NUM_VARS; k++) dQ[k] = QR[k] - QL[k];
		double rhoL = QL[1] + QL[2], rhoR = QR[1] + QR[2];
		double vyL = QL[4]/rhoL,     vyR = QR[4]/rhoR;
		double S = fmax(fabs(vyL), fabs(vyR));
		S = vy_max;
		// >>> 能量行耗散只取动能跳 (Battisti-Boscheri §3.3), 不含内能
		double kjmp[NUM_VARS];
		for(int k = 0; k < NUM_VARS; k++) kjmp[k] = dQ[k];
		{
			double keL = 0.5*(QL[3]*QL[3]+QL[4]*QL[4])/(rhoL*rhoL);
			double keR = 0.5*(QR[3]*QR[3]+QR[4]*QR[4])/(rhoR*rhoR);
			kjmp[5] = QR[1]*keR - QL[1]*keL;
			kjmp[6] = QR[2]*keR - QL[2]*keL;
		}
		for(int r = 0; r < NUM_VARS; r++){
			double Av = 0.0;
			for(int c = 0; c < NUM_VARS; c++) Av += Ah[r][c]*dQ[c];
			Dm[r] = 0.5*Av - 0.5*S*kjmp[r];
			Dp[r] = 0.5*Av + 0.5*S*kjmp[r];
		}
	}

// ---- 显式更新：u_con = Qs + adt * L --------------------------------------
//  w_xL/w_xR/w_yL/w_yR 由外部重构好后传入（cell 的左右/上下边界外推值）
void Path_RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,
	Grid &w_yL, Grid &w_yR,
	const EquationOfState& eos1, const EquationOfState& eos2){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				
				double QL[NUM_VARS], QR[NUM_VARS];
				double Dm[NUM_VARS], Dp[NUM_VARS];
				double L[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++) L[k] = 0.0;
				
				// ---- x 界面----------------------------------------------------
				// 右界面 i+1/2 : cell i 收 D^-
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_xR(i,j,k);
					QR[k]=w_xL(i+1,j,k); }
				Path_RusanovC_x(QL, QR, Dm, Dp, eos1, eos2);
				for(int k=0;k<NUM_VARS;k++) L[k] -= Dm[k]/dx;
				// 左界面 i-1/2 : cell i 收 D^+
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_xR(i-1,j,k);
					QR[k]=w_xL(i,j,k); }
				Path_RusanovC_x(QL, QR, Dm, Dp, eos1, eos2);
				for(int k=0;k<NUM_VARS;k++) L[k] -= Dp[k]/dx;
				
				// ---- y 界面----------------------------------------------------
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_yR(i,j,k);
					QR[k]=w_yL(i,j+1,k); }
				Path_RusanovC_y(QL, QR, Dm, Dp, eos1, eos2);
				for(int k=0;k<NUM_VARS;k++) L[k] -= Dm[k]/dy;
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_yR(i,j-1,k);
					QR[k]=w_yL(i,j,k); }
				Path_RusanovC_y(QL, QR, Dm, Dp, eos1, eos2);
				for(int k=0;k<NUM_VARS;k++) L[k] -= Dp[k]/dy;
				
				// ---- in-cell B 项 (x): path w_xL(i) -> w_xR(i)-----------------
				double Ah[NUM_VARS][NUM_VARS], dQc[NUM_VARS];
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_xL(i,j,k);
					QR[k]=w_xR(i,j,k); }
				pathAvgAx(QL, QR, Ah);
				for(int k=0;k<NUM_VARS;k++) dQc[k] = QR[k] - QL[k];
				for(int r=0;r<NUM_VARS;r++){
					double Bv = 0.0;
					for(int c=0;c<NUM_VARS;c++) Bv += Ah[r][c]*dQc[c];
					L[r] -= Bv/dx;
				}
				// ---- in-cell B 项 (y): path w_yL(i) -> w_yR(i)-----------------
				for(int k=0;k<NUM_VARS;k++){ QL[k]=w_yL(i,j,k);
					QR[k]=w_yR(i,j,k); }
				pathAvgAy(QL, QR, Ah);
				for(int k=0;k<NUM_VARS;k++) dQc[k] = QR[k] - QL[k];
				for(int r=0;r<NUM_VARS;r++){
					double Bv = 0.0;
					for(int c=0;c<NUM_VARS;c++) Bv += Ah[r][c]*dQc[c];
					L[r] -= Bv/dy;
				}
				
				// ---- 更新全部 7 个变量----------------------------------------
				// 能量行 (5,6) 已是动能通量 →自动「只更新动能部分」，enthalpy 留给 implicit 步
				for(int k=0;k<NUM_VARS;k++)
					u_con(i,j,k) = Qs(i,j,k) + adt * L[k];
				
				// α 安全 clamp
				u_con(i,j,0) = fmax(ALPHA_MIN, fmin(u_con(i,j,0), 1.0 -ALPHA_MIN));
			}
		}
	}


void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,
	const EquationOfState& eos1, const EquationOfState& eos2){
		double c;
		// reconstruct in PRIMITIVE variables so p,u stay continuous across
		// material interfaces (no spurious pressure overshoot)
		Grid u_pri(nxCells, nyCells);
		ucon2upri(u_pri, u_con, eos1, eos2);
		
		Grid wpri_xL(nxCells, nyCells), wpri_xR(nxCells, nyCells);
		Grid wpri_yL(nxCells, nyCells), wpri_yR(nxCells, nyCells);

		const double pinf1 = eos1.getPInf(), pinf2 = eos2.getPInf();
		const double RHO_FLOOR = 1e-10, P_EPS = 1e-6;   // 正性下限
		// 正性保护 (Zhang-Shu): 给定 cell 均值 ubar[] 与一维(已 minmod)斜率 c[],
		// 求最大缩放 θ∈[0,1] 使两侧面值 ubar±0.5θc 满足 α∈[αmin,1-αmin]、ρk>0、
		// pk+p∞,k>0. 所有约束对 θ 线性 -> 解析取最小. 接触/激波处自动退化到一阶,
		// 光滑区 θ=1 (二阶不变).
		auto posTheta = [&](const double* ubar, const double* c)->double{
			double th = 1.0;
			auto lim=[&](double base,double half,double floor){   // base ± half*θ ≥ floor (双侧)
				double a=fabs(half); if(a<=1e-300) return;
				double t=(base-floor)/a;                // 来自 base - a*θ ≥ floor
				if(t<th) th=(t<0.0?0.0:t);
			};
			lim(ubar[0], 0.5*c[0], ALPHA_MIN);             // α ≥ αmin
			lim(1.0-ubar[0], -0.5*c[0], ALPHA_MIN);        // α ≤ 1-αmin
			lim(ubar[1], 0.5*c[1], RHO_FLOOR);             // ρ1 > 0
			lim(ubar[2], 0.5*c[2], RHO_FLOOR);             // ρ2 > 0
			// p1+pinf1>0 <=> (α1p1)+α1*pinf1 > 0 (α1>0); 线性组合 5 与 0
			lim(ubar[5]+ubar[0]*pinf1, 0.5*(c[5]+c[0]*pinf1), P_EPS);
			// p2+pinf2>0 <=> (α2p2)+α2*pinf2 > 0, α2=1-α1
			lim(ubar[6]+(1.0-ubar[0])*pinf2, 0.5*(c[6]-c[0]*pinf2), P_EPS);
			return th;
		};

		for(int j = ghost_cells - 1; j < nyCells + ghost_cells + 1; j++){   // >>> FIX: 覆盖 ghost 行
			for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
				if(recon == First_Order){
					for(int k = 0; k < NUM_VARS; k++){
						wpri_xL(i,j,k)=wpri_xR(i,j,k)=wpri_yL(i,j,k)=wpri_yR(i,j,k)=u_pri(i,j,k);
					}
					continue;
				}
				// Second_Order: 先算各分量 minmod 斜率 (x,y 分开)
				double ubar[NUM_VARS], cx[NUM_VARS], cy[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					ubar[k]=u_pri(i,j,k);
					double sR = u_pri(i+1,j,k)-u_pri(i,j,k), sL = u_pri(i,j,k)-u_pri(i-1,j,k);
					if(sL*sR<=0) cx[k]=0; else cx[k]=(fabs(sR)<=fabs(sL))?sR:sL;
					double tR = u_pri(i,j+1,k)-u_pri(i,j,k), tL = u_pri(i,j,k)-u_pri(i,j-1,k);
					if(tL*tR<=0) cy[k]=0; else cy[k]=(fabs(tR)<=fabs(tL))?tR:tL;
				}
				double thx = posTheta(ubar, cx), thy = posTheta(ubar, cy);
				for(int k = 0; k < NUM_VARS; k++){
					wpri_xL(i,j,k)=ubar[k]-0.5*thx*cx[k];
					wpri_xR(i,j,k)=ubar[k]+0.5*thx*cx[k];
					wpri_yL(i,j,k)=ubar[k]-0.5*thy*cy[k];
					wpri_yR(i,j,k)=ubar[k]+0.5*thy*cy[k];
				}
			}
		}
		
		// map primitive edge states back to conservative edge states
		double pe[NUM_VARS], ce[NUM_VARS];
		for(int j = ghost_cells - 1; j < nyCells + ghost_cells + 1; j++){
			for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
				for(int k = 0; k < NUM_VARS; k++) pe[k] = wpri_xL(i,j,k);
				ini_pri2con(ce, pe, eos1, eos2); for(int k=0;k<NUM_VARS;k++) w_xL(i,j,k)=ce[k];
				
				for(int k = 0; k < NUM_VARS; k++) pe[k] = wpri_xR(i,j,k);
				ini_pri2con(ce, pe, eos1, eos2); for(int k=0;k<NUM_VARS;k++) w_xR(i,j,k)=ce[k];
				
				for(int k = 0; k < NUM_VARS; k++) pe[k] = wpri_yL(i,j,k);
				ini_pri2con(ce, pe, eos1, eos2); for(int k=0;k<NUM_VARS;k++) w_yL(i,j,k)=ce[k];
				
				for(int k = 0; k < NUM_VARS; k++) pe[k] = wpri_yR(i,j,k);
				ini_pri2con(ce, pe, eos1, eos2); for(int k=0;k<NUM_VARS;k++) w_yR(i,j,k)=ce[k];
			}
		}
		bound_cond_u(w_xL); bound_cond_u(w_xR);
		bound_cond_u(w_yL); bound_cond_u(w_yR);
	}

void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p){
	vector<double> result;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			double a1, rho1, rho2, vx, vy, p;
			switch (test) {
				case testshockbubble:
					if(x <= 0.05){            // 激波后 phase-1 气体
						a1 = 1.0 - 1e-6; rho1 = 1.3333; rho2 = 1e3;
						vx = 0.3535*sqrt(1e5); vy = 0; p = 1.5e5;
					}else if(x <= 0.4){       // 激波前 phase-1 气体
						a1 = 1.0 - 1e-6; rho1 = 1.0;    rho2 = 1e3;
						vx = 0; vy = 0; p = 1e5;
					}else if(x <= 0.6){       // phase-2 气泡(轻)
						a1 = 1e-6;       rho1 = 1e3;    rho2 = 0.1379;
						vx = 0; vy = 0; p = 1e5;
					}else{                    // phase-1 气体
						a1 = 1.0 - 1e-6; rho1 = 1.0;    rho2 = 1e3;
						vx = 0; vy = 0; p = 1e5;
					}
					u(i,j,0) = a1;              // α1
					u(i,j,1) = rho1;           // ρ1   (相密度,不是分质量)
					u(i,j,2) = rho2;           // ρ2
					u(i,j,3) = vx;             // vx   (速度,不是动量)
					u(i,j,4) = vy;             // vy
					u(i,j,5) = a1       * p;   // α1 p1   ← α 加权!
					u(i,j,6) = (1.0-a1) * p;
					break;
				case test2dahshockbubble:{
					const double SHOCK_POS = 50.0;       // 激波位置 (mm)
					const double BUBBLE_CX = 175.0;      // 气泡圆心 x (mm)
					const double BUBBLE_CY = 0.0;        // 气泡圆心 y (mm)
					const double BUBBLE_R  = 25.0;       // 气泡半径 (mm)
					const double ALPHA_MIN = 1e-6;       // 极小体积分数
					
					double x = x_0 + (i - ghost_cells + 0.5) * dx;
					double y = y_0 + (j - ghost_cells + 0.5) * dy;
					
					// 判断是否在气泡内
					double dxb = x - BUBBLE_CX;
					double dyb = y - BUBBLE_CY;
					bool in_bubble = (dxb*dxb + dyb*dyb) < (BUBBLE_R * BUBBLE_R);
					
					double a1, rho1, rho2, vx, vy, p;
					
					if (in_bubble) {
						// 氦气泡内部
						a1   = ALPHA_MIN;                 // 空气体积分数极小
						rho1 = 1.3764;                      // 微量空气密度（任意，不会影响）
						rho2 = 0.138;                    // 氦气密度
						vx   = 0.0;
						vy   = 0.0;
						p    = 1.0;                      // 气泡内压力与周围激波前空气一致
					} else if (x < SHOCK_POS) {
						// 激波后空气区域
						a1   = 1.0 - ALPHA_MIN;
						rho1 = 1.3764;
						rho2 = 0.138;                    // 微量氦气（体积分数极小）
						vx   = 0.394;
						vy   = 0.0;
						p    = 1.5698;
					} else {
						// 激波前空气区域
						a1   = 1.0 - ALPHA_MIN;
						rho1 = 1.0;
						rho2 = 0.138;                    // 微量氦气
						vx   = 0.0;
						vy   = 0.0;
						p    = 1.0;
					}
					
					// 六方程存储：α₁, ρ₁, ρ₂, u, v, α₁p₁, α₂p₂
					u(i,j,0) = a1;
					u(i,j,1) = rho1;
					u(i,j,2) = rho2;
					u(i,j,3) = vx;
					u(i,j,4) = vy;
					u(i,j,5) = a1 * p;          // α₁p₁
					u(i,j,6) = (1.0 - a1) * p;  // α₂p₂
					break;
				}
				case testDodecaneShockTube:{
					const double x_disc = 0.75;
					const double ALPHA_MIN = 1e-8;
					
					// 左侧状态（液体为主）
					const double alpha1_L = ALPHA_MIN;          // 蒸气体积分数
					const double rho1_L   = 2.0;               // 蒸气密度
					const double rho2_L   = 500.0;             // 液体密度
					const double u_L      = 0.0;
					const double p_L      = 1.0e8;
					
					// 右侧状态（蒸气为主）
					const double alpha1_R = 1.0 - ALPHA_MIN;   // 蒸气体积分数
					const double rho1_R   = 2.0;
					const double rho2_R   = 500.0;
					const double u_R      = 0.0;
					const double p_R      = 1.0e5;
					
					double x = x_0 + (i - ghost_cells + 0.5) * dx;
					bool left = (x < x_disc);
					
					double alpha1 = left ? alpha1_L : alpha1_R;
					double rho1   = left ? rho1_L   : rho1_R;
					double rho2   = left ? rho2_L   : rho2_R;
					double vx      = left ? u_L      : u_R;
					double p      = left ? p_L      : p_R;
					double vy      = 0.0;  // 1D
					
					// 存储原始变量（假设数组顺序: α, ρ₁, ρ₂, u, v, p）
					u(i,j,0) = alpha1;
					u(i,j,1) = rho1;
					u(i,j,2) = rho2;
					u(i,j,3) = vx;
					u(i,j,4) = vy;
					u(i,j,5) = alpha1 * p;      // α₁p₁
					u(i,j,6) = (1.0 - alpha1) * p;  // α₂p₂
					break;
				}
					
			}
		}
	}
	bound_cond_u(u);
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
					
					u(ghost,j,1) = u(2 * ghost_cells - 1 - ghost,j,1);
					u(nxCells + ghost_cells + ghost,j,1) = u(nxCells + ghost_cells - 1 - ghost,j,1);
					
					u(ghost,j,2) = u(2 * ghost_cells - 1 - ghost,j,2);
					u(nxCells + ghost_cells + ghost,j,2) = u(nxCells + ghost_cells - 1 - ghost,j,2);
					
					u(ghost,j,3) = -u(2 * ghost_cells - 1 - ghost,j,3);
					u(nxCells + ghost_cells + ghost,j,3) = -u(nxCells + ghost_cells - 1 - ghost,j,3);
					
					u(ghost,j,4) = -u(2 * ghost_cells - 1 - ghost,j,4);
					u(nxCells + ghost_cells + ghost,j,4) = -u(nxCells + ghost_cells - 1 - ghost,j,4);
					
					u(ghost,j,5) = u(2 * ghost_cells - 1 - ghost,j,5);
					u(nxCells + ghost_cells + ghost,j,5) = u(nxCells + ghost_cells - 1 - ghost,j,5);
					
					u(ghost,j,6) = u(2 * ghost_cells - 1 - ghost,j,6);
					u(nxCells + ghost_cells + ghost,j,6) = u(nxCells + ghost_cells - 1 - ghost,j,6);
					
					u(ghost,j,7) = u(2 * ghost_cells - 1 - ghost,j,7);
					u(nxCells + ghost_cells + ghost,j,7) = u(nxCells + ghost_cells - 1 - ghost,j,7);
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
					
					u(i,ghost,1) = u(i,2 * ghost_cells - 1 - ghost,1);
					u(i,nyCells + ghost_cells + ghost,1) = u(i,nyCells + ghost_cells - 1 - ghost,1);
					
					u(i,ghost,2) = u(i,2 * ghost_cells - 1 - ghost,2);
					u(i,nyCells + ghost_cells + ghost,2) = u(i,nyCells + ghost_cells - 1 - ghost,2);
					
					u(i,ghost,3) = -u(i,2 * ghost_cells - 1 - ghost,3);
					u(i,nyCells + ghost_cells + ghost,3) = -u(i,nyCells + ghost_cells - 1 - ghost,3);
					
					u(i,ghost,4) = -u(i,2 * ghost_cells - 1 - ghost,4);
					u(i,nyCells + ghost_cells + ghost,4) = -u(i,nyCells + ghost_cells - 1 - ghost,4);
					
					u(i,ghost,5) = u(i,2 * ghost_cells - 1 - ghost,5);
					u(i,nyCells + ghost_cells + ghost,5) = u(i,nyCells + ghost_cells - 1 - ghost,5);
					
					u(i,ghost,6) = u(i,2 * ghost_cells - 1 - ghost,6);
					u(i,nyCells + ghost_cells + ghost,6) = u(i,nyCells + ghost_cells - 1 - ghost,6);
					
					u(i,ghost,7) = u(i,2 * ghost_cells - 1 - ghost,7);
					u(i,nyCells + ghost_cells + ghost,7) = u(i,nyCells + ghost_cells - 1 - ghost,7);
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
	}
}

// frozen mixture sound speed at one cell:  cf^2 = sum_k Y_k c_k^2
double calMixSoundSpeed(double *con, const EquationOfState& eos1,
	const EquationOfState& eos2) {
		double alpha1 = con[0];
		double alpha2 = 1.0 - alpha1;
		double rho1 = con[1] / alpha1;
		double rho2 = con[2] / alpha2;
		double rho_mix = con[1] + con[2];
		
		double vx = con[3] / rho_mix;
		double vy = con[4] / rho_mix;
		
		// 从相总能计算相压力
		double alpha1p1 = eos1.getAlphaP_fromAlphaE(con[5], con[1], vx, vy, alpha1);
		double alpha2p2 = eos2.getAlphaP_fromAlphaE(con[6], con[2], vx, vy, alpha2);
		double p1 = alpha1p1 / alpha1;
		double p2 = alpha2p2 / alpha2;
		
		// 使用 EOS 的 calcs 方法获得相声速
		double c1 = eos1.calcs(rho1, p1);
		double c2 = eos2.calcs(rho2, p2);
		
		// Wood 公式
		double inv_rho_c2 = alpha1 / (rho1 * c1 * c1) + alpha2 / (rho2 * c2 * c2);
		double c_mix = 1.0 / sqrt(rho_mix * inv_rho_c2);
		
		return c_mix;
	}

double caldt(Grid &u_pri,const EquationOfState& eos1,const EquationOfState& eos2){
	vx_max = 0, vy_max = 0;
	double dt;
	double cf_max = 0.0;
	double con[NUM_VARS],pri[NUM_VARS];
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			if(fabs(u_pri(i,j,3)) > fabs(vx_max)) vx_max = fabs(u_pri(i,j,3));
			if(fabs(u_pri(i,j,4)) > fabs(vy_max)) vy_max = fabs(u_pri(i,j,4));
			
			for(int k = 0; k < NUM_VARS; k++) pri[k] = u_pri(i,j,k);
			ini_pri2con(con, pri, eos1, eos2);
			double cf = calMixSoundSpeed(con, eos1, eos2);
			if(cf > cf_max) cf_max = cf;
		}
	}
	
	//timestep == 1
	switch(caltimestep){
		case Fullydt:
			a_max_x = vx_max + cf_max;
			a_max_y = vy_max + cf_max;
			break;
		case IMEXdt:
			if(timestep < 5){
				a_max_x = vx_max + cf_max;
				a_max_y = vy_max + cf_max;
				//std::cout << "ch: "<< ch << endl;
			}else{
				a_max_x = vx_max;
				a_max_y = vy_max;
				//dt = CFL * 1 / (vx_max / dx);
			}
			break;
	}
	dt = CFL * fmin(dx,dy) / fmax(a_max_x,a_max_y);
	
	static double dt_prev = 0.0;
	if(timestep > 1 && dt_prev > 0) dt = fmin(dt, 1.2 * dt_prev);
	dt_prev = dt;
	return dt;
}

void bound_cond_array(Array &u){
	switch (x_bc) {
		case x_Transmissive:
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
		case y_Transmissive:
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

// ---- cell 混合 frozen 声速（α 解析约掉，只除以总密度 ρ）-------------------
//  c² = [ γ1(α1p1+α1π1) + γ2(α2p2+α2π2) ] / ρ
static double cellSound(const double* Q,
	const EquationOfState& eos1, const EquationOfState& eos2) {
		double alpha1 = Q[0];
		double alpha2 = 1.0 - alpha1;
		double rho1 = Q[1] / alpha1;
		double rho2 = Q[2] / alpha2;
		double rho_mix = Q[1] + Q[2];
		double u = Q[3] / rho_mix;
		double v = Q[4] / rho_mix;
		
		// 从相总能计算相压力
		double alpha1p1 = eos1.getAlphaP_fromAlphaE(Q[5], Q[1], u, v, alpha1);
		double alpha2p2 = eos2.getAlphaP_fromAlphaE(Q[6], Q[2], u, v, alpha2);
		double p1 = alpha1p1 / alpha1;
		double p2 = alpha2p2 / alpha2;
		
		// 使用各相 EOS 的 calcs 方法获得声速
		double c1 = eos1.calcs(rho1, p1);
		double c2 = eos2.calcs(rho2, p2);
		
		// Wood 公式计算混合声速
		double inv_rho_c2 = alpha1 / (rho1 * c1 * c1) + alpha2 / (rho2 * c2 * c2);
		double c_mix = 1.0 / sqrt(rho_mix * inv_rho_c2);
		
		return c_mix;
	}

// ---- Balsara flattener shock detector χ ∈ [0,1]，存进 Array chi -----------
void ShockDetector(Array &chi,Grid &u_con, 
	const EquationOfState& eos1, const EquationOfState& eos2){
		const int NTx = nxCells + 2*ghost_cells;
		const int NTy = nyCells + 2*ghost_cells;
		const double k1    = 1e-3;
		const double delta = fmin(dx, dy);
		
		Array U(nxCells,nyCells), V(nxCells,nyCells), cs(nxCells,nyCells),
		chit(nxCells,nyCells);
		
		// 1) 速度场 + 声速（整个含 ghost 网格）
		for(int j=0;j<NTy;j++) for(int i=0;i<NTx;i++){
			double Q[NUM_VARS]; for(int k=0;k<NUM_VARS;k++) Q[k]=u_con(i,j,k);
			double rho = Q[1] + Q[2];
			U(i,j)  = Q[3]/rho;
			V(i,j)  = Q[4]/rho;
			cs(i,j) = cellSound(Q, eos1, eos2);
		}
		
		// 2) χ̃ : 中心差分散度 + von Neumann 邻域最小声速（内部 +1 层）
		for(int j=1;j<NTy-1;j++) for(int i=1;i<NTx-1;i++){
			double divu = (U(i+1,j)-U(i-1,j))/(2.0*dx)
			+ (V(i,j+1)-V(i,j-1))/(2.0*dy);
			double amin = cs(i,j);
			amin = fmin(amin, cs(i-1,j)); amin = fmin(amin, cs(i+1,j));
			amin = fmin(amin, cs(i,j-1)); amin = fmin(amin, cs(i,j+1));
			double arg = -(delta*divu + k1*amin)/(k1*amin);
			chit(i,j)  = fmin(1.0, fmax(0.0, arg));
		}
		
		// 3) χ = 3×3 stencil 取 max（内部格）
		for(int j=ghost_cells;j<nyCells+ghost_cells;j++)
			for(int i=ghost_cells;i<nxCells+ghost_cells;i++){
				double m = 0.0;
				for(int dj=-1;dj<=1;dj++) for(int di=-1;di<=1;di++)
					m = fmax(m, chit(i+di,j+dj));
				chi(i,j) = m;
			}
		Array P(nxCells, nyCells);
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
				double con[NUM_VARS], pri[NUM_VARS];
				for (int k = 0; k < NUM_VARS; k++) con[k] = u_con(i, j, k);
				ini_con2pri(pri, con, eos1, eos2);
				P(i, j) = pri[5] + pri[6];   // 混合压力
			}
		}
		
		// 2. 多次扫描（2~3次）以确保传播链式进行
		const int num_iter = 2;   // 通常2次足够，若激波极强可增至3
		for (int iter = 0; iter < num_iter; iter++) {
			
			// x 方向传播
			for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
				for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
					if (chi(i, j) > 0.0) {   // 当前单元是激波区
						// 向右传播（i+1）
						if (i + 1 < nxCells + ghost_cells && chi(i+1, j) == 0.0 && P(i, j) > P(i+1, j))
							chi(i+1, j) = chi(i, j);
						// 向左传播（i-1）
						if (i - 1 >= ghost_cells && chi(i-1, j) == 0.0 && P(i, j) > P(i-1, j))
							chi(i-1, j) = chi(i, j);
					}
				}
			}
			
			// y 方向传播
			for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
				for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
					if (chi(i, j) > 0.0) {
						// 向上传播（j+1）
						if (j + 1 < nyCells + ghost_cells && chi(i, j+1) == 0.0 && P(i, j) > P(i, j+1))
							chi(i, j+1) = chi(i, j);
						// 向下传播（j-1）
						if (j - 1 >= ghost_cells && chi(i, j-1) == 0.0 && P(i, j) > P(i, j-1))
							chi(i, j-1) = chi(i, j);
					}
				}
			}
		}
		bound_cond_array(chi);
	}

void calEnthalpy_phasic(double *con, double &h1, double &h2,const EquationOfState& eos1, const EquationOfState& eos2){
		double alpha1 = con[0];
		double alpha2 = 1.0 - alpha1;
		double rho    = con[1] + con[2];        // 混合密度
		double u      = con[3]/rho;             // 速度(动量/密度)
		double v      = con[4]/rho;
		double ke     = 0.5*(u*u + v*v);        // 比动能
		
		double m1 = con[1], m2 = con[2];        // 相质量
		
		// 相内能密度 = 相总能 - 相动能(动能用相质量)
		double arhoe1 = con[5] - m1*ke;         // α₁ρ₁e₁
		double arhoe2 = con[6] - m2*ke;         // α₂ρ₂e₂
		
		double ap1 = eos1.getAlphaP(arhoe1, con[1] / alpha1,alpha1);   // α₁p₁
		double ap2 = eos2.getAlphaP(arhoe2, con[2] / alpha2,alpha2);   // α₂p₂
		
		// 焓 h_k = (α_k ρ_k e_k + α_k p_k) / 相质量
		//  = (arhoe + ap) / m_k,只除相质量不除 α
		h1 = (m1 > 1e-30) ? (arhoe1 + ap1)/m1 : 0.0;
		h2 = (m2 > 1e-30) ? (arhoe2 + ap2)/m2 : 0.0;
	}

double calkinetic(double *con){
	double rho = con[1] + con[2];
	double vx = con[3] / rho;
	double vy = con[4] / rho;
	double vz = con[5] / rho;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy);
	//double kinetic_energy = 0.5 * (con[1] * con[1] + con[2] * con[2] + con[3] * con[3]) / (con[0]);
	return kinetic_energy;
}


int M_P(int i, int j, int k){
	int coord = 0;
	int I, J;
	
	// ---- x 方向边界处理(算出有效索引 I)----
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive: I = ghost_cells;                       break;
			case x_Reflective:   I = 2 * ghost_cells - 1 - i;           break;
			case x_Periodic:     I = nxCells + i;                       break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: I = nxCells + ghost_cells - 1;         break;
			case x_Reflective:   I = 2 * (nxCells + ghost_cells) - i - 1; break;
			case x_Periodic:     I = i - nxCells;                       break;
		}
	}else{
		I = i;
	}
	
	// ---- y 方向边界处理(算出有效索引 J)----
	if(j < ghost_cells){
		switch(y_bc){
			case y_Transmissive: J = ghost_cells;                       break;
			case y_Reflective:   J = 2 * ghost_cells - 1 - j;           break;
			case y_Periodic:     J = nyCells + j;                       break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive: J = nyCells + ghost_cells - 1;         break;
			case y_Reflective:   J = 2 * (nyCells + ghost_cells) - j - 1; break;
			case y_Periodic:     J = j - nyCells;                       break;
		}
	}else{
		J = j;
	}
	
	// ---- interleaved 布局:同一 cell 的两相压力相邻 ----
	// cell_index = (I-gc) + nx*(J-gc),范围 [0, N)
	// 每个 cell 占 2 个未知量:相 k=0 在 2*cell,相 k=1 在 2*cell+1
	int cell_index = (I - ghost_cells) + nxCells * (J - ghost_cells);
	coord = 2 * cell_index + k;     // k = 0 -> p1, k = 1 -> p2
	
	return coord;
}

void get_phase_linear_coeffs(double rho, const EquationOfState& eos,double& A, double& B) {
		double c, d;
		eos.e_linear_coeff(rho, c, d);
		A = rho * c;
		B = rho * d;
	}

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x,Grid& u_con, Grid& Qe, Grid& Qs, Array &chi,double h,const EquationOfState& eos1, const EquationOfState& eos2){
		const double gamma1 = eos1.getGamma();
		const double gamma2 = eos2.getGamma();
		const double pinf_1 = eos1.getPInf();
		const double pinf_2 = eos2.getPInf();
		
		// A_k = 1/Γ_k(ρ_k). 需要逐格取值: 下面把 A_k 放进散度里面.
		auto Acoef = [&](int ii,int jj,int phase)->double{
			double an = phase==0 ? u_con(ii,jj,0) : 1.0-u_con(ii,jj,0);
			double Aa,Bb;
			get_phase_linear_coeffs(u_con(ii,jj,1+phase)/fmax(an,ALPHA_MIN),
				phase==0 ? eos1 : eos2, Aa, Bb);
			return Aa;
		};

		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				
				// ===== 0. 取邻居 con,转 pri =====
				double con00_n[NUM_VARS], pri00_n[NUM_VARS];
				double conm0_n[NUM_VARS], prim0_n[NUM_VARS];
				double conp0_n[NUM_VARS], prip0_n[NUM_VARS];
				double con0m_n[NUM_VARS], pri0m_n[NUM_VARS];
				double con0p_n[NUM_VARS], pri0p_n[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					con00_n[k]=Qe(i,j,k);   conm0_n[k]=Qe(i-1,j,k); conp0_n[k]=Qe(i+1,j,k);
					con0m_n[k]=Qe(i,j-1,k); con0p_n[k]=Qe(i,j+1,k);
				}
				ini_con2pri(pri00_n, con00_n, eos1, eos2);
				ini_con2pri(prim0_n, conm0_n, eos1, eos2);
				ini_con2pri(prip0_n, conp0_n, eos1, eos2);
				ini_con2pri(pri0m_n, con0m_n, eos1, eos2);
				ini_con2pri(pri0p_n, con0p_n, eos1, eos2);
				
				// ===================== Σ : path-consistent 离散 =====================
				//  用 (6g) 的展开形式, 而不是 -u·(Y2∇(α1p1) - Y1∇(α2p2)):
				//      Σ = -u·[ (Y2 p1 + Y1 p2) ∇α1 + α1 Y2 ∇p1 - α2 Y1 ∇p2 ]
				//  这样 α 梯度和 p 梯度被显式分开:
				//   * ∇α1 是显式 path-conservative 子系统刚更新过的量, 必须用同一个
				//     PC 算子. α 只被显式步改变, 所以该算子就是增量本身:
				//         (u·∇α1)^PC ≡ ( α1^{Qs} - α1^{n+1} ) / h
				//     逐位继承重构 / 路径积分 / Rusanov 耗散 / clamp, 无需重新差分,
				//     也不需要任何失配修正项 (原 Rcorr 因此被删除).
				//   * ∇p_k 是压力梯度, 在接触间断上恒为 0, 保持中心差分 ->
				//     Rusanov 耗散不会进入椭圆算子 (避免把隐式系统变成抛物型).
				double rho = u_con(i,j,1)+u_con(i,j,2);
				double Y1 = u_con(i,j,1)/rho;
				double Y2 = u_con(i,j,2)/rho;
				double alpha1 = u_con(i,j,0);
				double alpha2 = 1 - alpha1;
				double vxs = u_con(i,j,3)/rho;      // u*  (显式步后的速度)
				double vys = u_con(i,j,4)/rho;      // v*

				const double Aadv = (Qs(i,j,0) - u_con(i,j,0)) / h;   // = (u·∇α1)^PC

				double a1c = fmax(u_con(i,j,0), ALPHA_MIN), a2c = fmax(1.0-u_con(i,j,0), ALPHA_MIN);

				double p_mix = x(M_P(i,j,0)) + x(M_P(i,j,1));
				double p_mix_p0 = x(M_P(i+1,j,0)) + x(M_P(i+1,j,1));
				double p_mix_m0 = x(M_P(i-1,j,0)) + x(M_P(i-1,j,1));
				double p_mix_0p = x(M_P(i,j+1,0)) + x(M_P(i,j+1,1));
				double p_mix_0m = x(M_P(i,j-1,0)) + x(M_P(i,j-1,1));

				// 用 α_k∇p_k = ∇q_k - p_k∇α_k 把 1/α 从"邻格"里消掉, 只留对角上的 p_k:
				//   Σ = -[ u*·(Y2 ∇q1 - Y1 ∇q2)|_中心差分
				//        + (Y2 p1 + Y1 p2)·( (u·∇α1)^PC - u*·∇α1|_中心差分 ) ]
				// 第一行与原来的 (6g) 左式离散逐字相同 -> 矩阵稀疏性/条件数不变;
				// 第二行是把 α 梯度从中心换成 PC 的那一项, 系数是 (6g) 自己给的
				// (Y2p1+Y1p2), 而不是外挂的修正.
				double Gq1x = ( x(M_P(i+1,j,0)) - x(M_P(i-1,j,0)) )/(2*dx);
				double Gq2x = ( x(M_P(i+1,j,1)) - x(M_P(i-1,j,1)) )/(2*dx);
				double Gq1y = ( x(M_P(i,j+1,0)) - x(M_P(i,j-1,0)) )/(2*dy);
				double Gq2y = ( x(M_P(i,j+1,1)) - x(M_P(i,j-1,1)) )/(2*dy);
				double Ga1x = ( u_con(i+1,j,0) - u_con(i-1,j,0) )/(2*dx);
				double Ga1y = ( u_con(i,j+1,0) - u_con(i,j-1,0) )/(2*dy);

				// (Y2 p1 + Y1 p2): 只用本格未知量 -> 纯对角项 (无邻格 1/α 放大)
				double Pw = Y2 * x(M_P(i,j,0))/a1c + Y1 * x(M_P(i,j,1))/a2c;

				double S_impl = vxs * ( Y2*Gq1x - Y1*Gq2x ) + vys * ( Y2*Gq1y - Y1*Gq2y )
				              + Pw * ( Aadv - vxs*Ga1x - vys*Ga1y );

				// (b) Picard 速度修正 u^{n+1} = u* - (h/ρ)∇p_mix^{n+1}, 方括号冻结在 n.
				//     接触间断上 ∇p_mix^{n+1} = 0, 该项对 WB 无贡献, 用中心差分即可.
				//     冻结括号保持原来的 (6g) 左式写法, 与原代码逐字相同.
				double Bx_n = Y2*(prip0_n[5]-prim0_n[5])/(2*dx) - Y1*(prip0_n[6]-prim0_n[6])/(2*dx);
				double By_n = Y2*(pri0p_n[5]-pri0m_n[5])/(2*dy) - Y1*(pri0p_n[6]-pri0m_n[6])/(2*dy);

				double sigma = -S_impl
				             + (h/rho) * ( (p_mix_p0-p_mix_m0)/(2*dx)*Bx_n + (p_mix_0p-p_mix_0m)/(2*dy)*By_n );
				
				
				// 算phase1 LHS
				double h1_00, h2_00, h1_p0, h2_p0, h1_m0, h2_m0, h1_0p, h2_0p,h1_0m, h2_0m;
				calEnthalpy_phasic(con00_n, h1_00,h2_00,eos1, eos2);
				calEnthalpy_phasic(conp0_n, h1_p0,h2_p0,eos1, eos2);
				calEnthalpy_phasic(conm0_n, h1_m0,h2_m0,eos1, eos2);
				calEnthalpy_phasic(con0p_n, h1_0p,h2_0p,eos1, eos2);
				calEnthalpy_phasic(con0m_n, h1_0m,h2_0m,eos1, eos2);
				h1_p0 = 0.5 * (h1_00 + h1_p0);
				h1_m0 = 0.5 * (h1_00 + h1_m0);
				h1_0p = 0.5 * (h1_00 + h1_0p);
				h1_0m = 0.5 * (h1_00 + h1_0m);
				h2_p0 = 0.5 * (h2_00 + h2_p0);
				h2_m0 = 0.5 * (h2_00 + h2_m0);
				h2_0p = 0.5 * (h2_00 + h2_0p);
				h2_0m = 0.5 * (h2_00 + h2_0m);
				auto Ck = [&](int ii,int jj,int phase)->double{
					// Boscheri (18): C_k = α_k^{n+1} · (ρ_k h_k)^n / ρ^{n+1},  (ρ_k h_k)^n = ((αρe)^n+(αp)^n)/α^n
					double cc[NUM_VARS], pp[NUM_VARS];
					for(int k=0;k<NUM_VARS;k++) cc[k]=Qe(ii,jj,k);
					ini_con2pri(pp,cc,eos1,eos2);
					double rr  = cc[1]+cc[2];
					double uu = cc[3]/rr, vv = cc[4]/rr;
					double ke = 0.5*(uu*uu+vv*vv);
					double arhoe = cc[5+phase] - cc[1+phase]*ke;        // (α_kρ_ke_k)^n
					double an    = phase==0 ? cc[0] : 1.0-cc[0];
					double rhk   = (arhoe + pp[5+phase]) / fmax(an,ALPHA_MIN);   // (ρ_k h_k)^n
					double an1   = phase==0 ? u_con(ii,jj,0) : 1.0-u_con(ii,jj,0);
					double rho_n1= u_con(ii,jj,1)+u_con(ii,jj,2);
					return an1 * rhk / rho_n1;
				};
				double C1_00=Ck(i,j,0), C2_00=Ck(i,j,1);
				double W1_p0=0.5*(C1_00+Ck(i+1,j,0)), W1_m0=0.5*(C1_00+Ck(i-1,j,0));
				double W1_0p=0.5*(C1_00+Ck(i,j+1,0)), W1_0m=0.5*(C1_00+Ck(i,j-1,0));
				double W2_p0=0.5*(C2_00+Ck(i+1,j,1)), W2_m0=0.5*(C2_00+Ck(i-1,j,1));
				double W2_0p=0.5*(C2_00+Ck(i,j+1,1)), W2_0m=0.5*(C2_00+Ck(i,j-1,1));
				
				double A1,B1,A2,B2;
				get_phase_linear_coeffs(u_con(i,j,1) / alpha1, eos1,A1, B1);
				get_phase_linear_coeffs(u_con(i,j,2) / alpha2, eos2,A2, B2);
				
				Ax(M_P(i,j,0)) = x(M_P(i,j,0)) * A1 - h * h * ( 1 / (dx * dx) * (W1_p0 * (p_mix_p0 - p_mix) - W1_m0 * (p_mix - p_mix_m0)) + 1 / (dy * dy) * (W1_0p * (p_mix_0p - p_mix) - W1_0m * (p_mix - p_mix_0m)) ) + h * sigma ;
				Ax(M_P(i,j,1)) = x(M_P(i,j,1)) * A2 - h * h * ( 1 / (dx * dx) * (W2_p0 * (p_mix_p0 - p_mix) - W2_m0 * (p_mix - p_mix_m0)) + 1 / (dy * dy) * (W2_0p * (p_mix_0p - p_mix) - W2_0m * (p_mix - p_mix_0m)) ) - h * sigma ;
				
				// ===== shock-detector Laplacian: -Δt/(γ_k-1)·∇·(κ ∇(α_k p_k)) =====
				// χ 已由 ShockDetector 预算 (Array chi, 含一层 ghost); κ 用对流速度
				
				// 各格 |速度|（分母用各格混合密度）
				
				
				// κ at cells (x 用 |u|·dx, y 用 |v|·dy)
				
				double kx_00 = chi(i,  j  ) * vx_max * dx;
				double kx_p0 = chi(i+1,j  ) * vx_max * dx;
				double kx_m0 = chi(i-1,j  ) * vx_max * dx;
				double ky_00 = chi(i,  j  ) * vy_max * dy;
				double ky_0p = chi(i,  j+1) * vy_max * dy;
				double ky_0m = chi(i,  j-1) * vy_max * dy;
				
//				double kx_00 = 0;
//				double kx_p0 = 0;
//				double kx_m0 = 0;
//				double ky_00 = 0;
//				double ky_0p = 0;
//				double ky_0m = 0;
				
//				double kx_00 = vx_max * dx;
//				double kx_p0 = vx_max * dx;
//				double kx_m0 = vx_max * dx;
//				double ky_00 = vy_max * dy;
//				double ky_0p = vy_max * dy;
//				double ky_0m = vy_max * dy;
				
//				// 各格混合密度
//				double r00 = u_con(i,  j,  1) + u_con(i,  j,  2);
//				double rp0 = u_con(i+1,j,  1) + u_con(i+1,j,  2);
//				double rm0 = u_con(i-1,j,  1) + u_con(i-1,j,  2);
//				double r0p = u_con(i,  j+1,1) + u_con(i,  j+1,2);
//				double r0m = u_con(i,  j-1,1) + u_con(i,  j-1,2);
//				
//				// 局部速度版 κ = |u_local|·dx (x),  |v_local|·dy (y)
//				double kx_00 = fabs(u_con(i,  j,  3)/r00) * dx;
//				double kx_p0 = fabs(u_con(i+1,j,  3)/rp0) * dx;
//				double kx_m0 = fabs(u_con(i-1,j,  3)/rm0) * dx;
//				double ky_00 = fabs(u_con(i,  j,  4)/r00) * dy;
//				double ky_0p = fabs(u_con(i,  j+1,4)/r0p) * dy;
//				double ky_0m = fabs(u_con(i,  j-1,4)/r0m) * dy;
				
				// face-average κ
				double kxR = 0.5*(kx_00 + kx_p0), kxL = 0.5*(kx_m0 + kx_00);
				double kyT = 0.5*(ky_00 + ky_0p), kyB = 0.5*(ky_0m + ky_00);
				
				// q_k = α_k p_k 邻居值（未知量）
				double q1_00 = x(M_P(i,  j,  0)), q1_p0 = x(M_P(i+1,j,  0)), q1_m0 = x(M_P(i-1,j,  0));
				double q1_0p = x(M_P(i,  j+1,0)), q1_0m = x(M_P(i,  j-1,0));
				double q2_00 = x(M_P(i,  j,  1)), q2_p0 = x(M_P(i+1,j,  1)), q2_m0 = x(M_P(i-1,j,  1));
				double q2_0p = x(M_P(i,  j+1,1)), q2_0m = x(M_P(i,  j-1,1));
				
				// ∇·(κ∇(A_k q_k)) —— A_k 写在散度**里面**, 不是乘在外面.
				// A_k q_k = α_kρ_k e_k - α_k B_k, 所以两相相加后恒为
				//   ∇·(κ∇(A1q1 + A2q2))
				// 的面差 -> 逐面对消, 对任意 Γ(ρ) 都严格守恒.
				// (写在外面时 A1≠A2 且逐格变化 -> 混合方程上不再是通量差;
				//  当前两个 EOS 类 Γ 都是常数, 两种写法数值上等价.)
				double A1_00=Acoef(i,j,0),   A1_p0=Acoef(i+1,j,0), A1_m0=Acoef(i-1,j,0);
				double A1_0p=Acoef(i,j+1,0), A1_0m=Acoef(i,j-1,0);
				double A2_00=Acoef(i,j,1),   A2_p0=Acoef(i+1,j,1), A2_m0=Acoef(i-1,j,1);
				double A2_0p=Acoef(i,j+1,1), A2_0m=Acoef(i,j-1,1);

				double lap_shock1 = 1.0/(dx*dx)*( kxR*(A1_p0*q1_p0-A1_00*q1_00) - kxL*(A1_00*q1_00-A1_m0*q1_m0) )
				                  + 1.0/(dy*dy)*( kyT*(A1_0p*q1_0p-A1_00*q1_00) - kyB*(A1_00*q1_00-A1_0m*q1_0m) );
				double lap_shock2 = 1.0/(dx*dx)*( kxR*(A2_p0*q2_p0-A2_00*q2_00) - kxL*(A2_00*q2_00-A2_m0*q2_m0) )
				                  + 1.0/(dy*dy)*( kyT*(A2_0p*q2_0p-A2_00*q2_00) - kyB*(A2_00*q2_00-A2_0m*q2_0m) );

				Ax(M_P(i,j,0)) += -h * lap_shock1;   // A_k 已在散度内
				Ax(M_P(i,j,1)) += -h * lap_shock2;
			}
		}
	}


//============================================================================
//  Boscheri 有限速率压力弛豫(SG,半隐线性化,一个除法,cell-local) —— 终版
//  输入:u_con = 压力系统解完、回代后的守恒量(相总能布局)
//  输出:更新 α(con[0])、相总能(con[5],con[6])重分配
//  守恒不变:相质量 con[1,2]、动量 con[3,4]、混合总能 con[5]+con[6]
//  斜率 K_k = dp_k/dα_k 含界面压力做功(de=-p_I dv);da 限幅防过冲/病态
//
// ============================================================================
void PressureRelaxation_Boscheri(Grid& u_con, double mu,const EquationOfState& eos1,const EquationOfState& eos2){
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			
			// ---- 弛豫前状态(上标 H = 压力系统解完的值)----
			double a1H = u_con(i,j,0);
			double a2H = 1.0 - a1H;
			
			// 近纯相:该相基本不存在,弛豫无意义且 p_k 噪声大 -> 跳过
			if(a1H < 10.0*ALPHA_MIN || a2H < 10.0*ALPHA_MIN) continue;
			
			double m1  = u_con(i,j,1);          // α1ρ1(守恒,弛豫不变)
			double m2  = u_con(i,j,2);          // α2ρ2
			double rho = m1 + m2;
			double u   = u_con(i,j,3)/rho;
			double v   = u_con(i,j,4)/rho;
			double ke  = 0.5*(u*u + v*v);       // 比动能(冻结)
			
			// 相密度
			double rho1H = m1/a1H;
			double rho2H = m2/a2H;
			
			// 相内能密度(相总能 - 相动能)
			double arhoe1 = u_con(i,j,5) - m1*ke;   // α1ρ1 e1
			double arhoe2 = u_con(i,j,6) - m2*ke;   // α2ρ2 e2
			
			// ---- 通用 EOS: 相压力 / 声速 / Grüneisen ----
			double e1H = arhoe1/m1, e2H = arhoe2/m2;        // 比内能
			double p1H = eos1.pressure(rho1H, e1H);
			double p2H = eos2.pressure(rho2H, e2H);
			double c1  = eos1.calcs(rho1H, p1H);
			double c2  = eos2.calcs(rho2H, p2H);
			double cl1,dl1,cl2,dl2;
			eos1.e_linear_coeff(rho1H, cl1, dl1);           // e = c_lin·p + d
			eos2.e_linear_coeff(rho2H, cl2, dl2);
			double Gam1 = 1.0/(cl1*rho1H), Gam2 = 1.0/(cl2*rho2H);   // Γ = (∂p/∂e)_ρ/ρ = 1/(c_lin·ρ)
			double Z1 = rho1H*c1, Z2 = rho2H*c2;
			double pI = (Z2*p1H + Z1*p2H)/(Z1 + Z2);

			// ---- 通用弛豫斜率 K_k = -(1/α)[ρc² + Γ(pI - p)] ----
			double K1 = -(1.0/a1H)*( rho1H*c1*c1 + Gam1*(pI - p1H) );
			double K2 = -(1.0/a2H)*( rho2H*c2*c2 + Gam2*(pI - p2H) );
			
			// ---- 一个除法解 Δα1,带 da_eq 限幅 ----
			//  健康格: denom=1-dtμ(K1+K2)=1+dtμ·S >1 (S=-(K1+K2)>0),半隐da<da_eq 自然成立
			//  病态格(denom≤0.5,即 K1+K2>0): 不放大,退化显式;再用 da_eq夹住不过冲到对面
			double denom = 1.0 - dt*mu*(K1 + K2);
			double da    = dt*mu*(p1H - p2H);          // 显式基值(有界)
			if(denom > 0.5) da /= denom;               // 健康:半隐(放大 ≤2×)
			
			// 线性化平衡位移(朝 p1=p2);S>0 时有意义,S≤0 时da_eq→大,限幅不绑定
			double S     = -(K1 + K2);
			double da_eq = (p1H - p2H)/fmax(S, 1e-12);
			if(fabs(da) > fabs(da_eq)) da = da_eq;      //绝不过冲平衡(也治大μ过冲)
			
			double a1_new = a1H + da;
			
			// clamp
			a1_new = std::max(ALPHA_MIN, std::min(a1_new, 1.0 - ALPHA_MIN));
			double a2_new = 1.0 - a1_new;
			
			// ---- 更新相密度(相质量守恒)----
			double rho1_new = m1/a1_new;
			double rho2_new = m2/a2_new;
			
			// ---- 更新相内能(弛豫路径 de=-p_I dv,两相同 p_I 保守恒)----
			double v1H = 1.0/rho1H,       v2H = 1.0/rho2H;
			double v1_new = 1.0/rho1_new, v2_new = 1.0/rho2_new;
			double e1_new = arhoe1/m1 - pI*(v1_new - v1H);   // arhoe/m =比内能
			double e2_new = arhoe2/m2 - pI*(v2_new - v2H);
			
			// ---- 写回守恒量 ----
			u_con(i,j,0) = a1_new;                  // α1
			u_con(i,j,5) = m1*e1_new + m1*ke;       // α1ρ1E1 = α1ρ1 e1 + 动能
			u_con(i,j,6) = m2*e2_new + m2*ke;       // α2ρ2E2
			// 相质量 con[1,2]、动量 con[3,4] 不变
		}
	}
}

static inline double safe_div(double num, double den, double eps = 1e-14) {
	if (std::abs(den) < eps) den = (den >= 0.0 ? eps : -eps);
	return num / den;
}


void GMRES_Solver(Eigen::VectorXd& x, Eigen::VectorXd& b,
	std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, Array&, double,
		const EquationOfState&, const EquationOfState&)> matrix_vector_product,
	Grid& u_con, Grid& Qe, Grid& Qs, Array& chi, double h,
	const EquationOfState& eos1, const EquationOfState& eos2, const int dim){
		
		// ---- restarted GMRES(m), matrix-free ----
		const int m = std::min(GMRESrest > 0 ? GMRESrest : 30, dim);   // Krylov 维数
		double norm_b = b.norm();
		double tol = GMREStol * std::max(norm_b, 1e-30);
		
		Eigen::VectorXd Ax(dim), r(dim), w(dim);
		Eigen::MatrixXd V(dim, m + 1);
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m + 1, m);
		Eigen::VectorXd cs(m), sn(m), g(m + 1), y(m);
		
		for (int outer = 0; outer < GMRESiter; ++outer){
			matrix_vector_product(Ax, x, u_con, Qe, Qs, chi, h, eos1, eos2);
			r = b - Ax;
			double beta = r.norm();
			if (beta < tol) return;
			
			V.col(0) = r / beta;
			g.setZero();
			g(0) = beta;
			
			int k = 0;
			for (; k < m; ++k){
				Eigen::VectorXd vk = V.col(k);
				matrix_vector_product(w, vk, u_con, Qe, Qs, chi, h, eos1, eos2);
				// modified Gram-Schmidt
				for (int i = 0; i <= k; ++i){
					H(i, k) = w.dot(V.col(i));
					w.noalias() -= H(i, k) * V.col(i);
				}
				H(k + 1, k) = w.norm();
				if (H(k + 1, k) > 1e-14) V.col(k + 1) = w / H(k + 1, k);   // ★ near-breakdown 保护
				// apply previous Givens rotations to column k
				for (int i = 0; i < k; ++i){
					double t0 =  cs(i) * H(i, k) + sn(i) * H(i + 1, k);
					H(i + 1, k) = -sn(i) * H(i, k) + cs(i) * H(i + 1, k);
					H(i, k) = t0;
				}
				// new Givens rotation
				double dnm = std::hypot(H(k, k), H(k + 1, k));
				cs(k) = H(k, k) / dnm;
				sn(k) = H(k + 1, k) / dnm;
				H(k, k) = cs(k) * H(k, k) + sn(k) * H(k + 1, k);
				H(k + 1, k) = 0.0;
				g(k + 1) = -sn(k) * g(k);
				g(k)     =  cs(k) * g(k);
				if (std::abs(g(k + 1)) < tol){ ++k; break; }
			}
			// back-substitution: H(0:k,0:k) y = g(0:k)
			for (int i = k - 1; i >= 0; --i){
				double s0 = g(i);
				for (int j = i + 1; j < k; ++j) s0 -= H(i, j) * y(j);
				y(i) = s0 / H(i, i);
			}
			for (int i = 0; i < k; ++i) x.noalias() += y(i) * V.col(i);
			
			if (std::abs(g(k)) < tol) return;   // converged
		}
		
		matrix_vector_product(Ax, x, u_con, Qe, Qs, chi, h, eos1, eos2);
		std::cerr << "Warning: GMRES not converged, |r|/|b| = "
		<< ((b - Ax).norm() / std::max(norm_b, 1e-30)) << std::endl;
	}

// ============================================================================
//  直接解压力系统: 用"着色探测"组装稀疏矩阵 (调用现有 matvec 几次, 冻结 EOS 系数
//  只算一遍), 再 SparseLU 直接解. 替代刚性 EOS 下又慢又难收敛的 matrix-free GMRES.
//  算子是线性的 (A(0)=0), 所以探测即得矩阵列. 带组装正确性自检.
// ============================================================================
void PressureSolve_Direct(Eigen::VectorXd& x, Eigen::VectorXd& b,
    Grid& u_con, Grid& Qe, Grid& Qs, Array& chi, double h,
    const EquationOfState& eos1, const EquationOfState& eos2, const int dim){
	const int gc = ghost_cells;
	const int S = 7;   // 着色步长 (>= 2*D+1)
	const int D = 2;   // 提取半径 (格)
	auto idx2cell = [&](int q, int& ci, int& cj, int& p){
		p = q % 2; int cell = q / 2;
		ci = cell % nxCells + gc;
		cj = cell / nxCells + gc;
	};
	std::vector<Eigen::Triplet<double>> trip; trip.reserve((size_t)dim*12);
	Eigen::VectorXd v(dim), w(dim);
	int ncol = 2*S*S;
	for(int color = 0; color < ncol; ++color){
		v.setZero(); int cnt=0;
		for(int q=0; q<dim; ++q){ int ci,cj,p; idx2cell(q,ci,cj,p);
			int cc = p + 2*((ci-gc)%S) + 2*S*((cj-gc)%S);
			if(cc==color){ v[q]=1.0; ++cnt; } }
		if(cnt==0) continue;                  // 空探测跳过 (1D 时大多数空)
		Pressure_Subsystem_AX(w, v, u_con, Qe, Qs, chi, h, eos1, eos2);
		for(int q=0; q<dim; ++q){ if(v[q]==0.0) continue;
			int ci,cj,p; idx2cell(q,ci,cj,p);
			for(int dj=-D; dj<=D; ++dj) for(int di=-D; di<=D; ++di){
				int ri=ci+di, rj=cj+dj;
				if(ri<gc||ri>=nxCells+gc||rj<gc||rj>=nyCells+gc) continue;
				for(int rp=0; rp<2; ++rp){ int r=M_P(ri,rj,rp);
					double val=w[r]; if(val!=0.0) trip.emplace_back(r,q,val); }
			}
		}
	}
	Eigen::SparseMatrix<double> A(dim,dim);
	A.setFromTriplets(trip.begin(), trip.end());
	// 自检: A·x_rand 应 == matvec(x_rand)
	static bool checked=false;
	if(!checked){
		Eigen::VectorXd xr=Eigen::VectorXd::Random(dim), w1(dim);
		Pressure_Subsystem_AX(w1, xr, u_con, Qe, Qs, chi, h, eos1, eos2);
		double err=(A*xr - w1).norm()/std::max(w1.norm(),1e-30);
		std::cerr<<"[direct] assembly rel-err="<<err<<" nnz="<<A.nonZeros()<<std::endl;
		checked=true;
	}
	// RHS 健康检查: vanishing-phase 除法可能产生 NaN/Inf
	static int rhswarn=0;
	if((!b.allFinite()) && rhswarn<3){
		std::cerr<<"[direct] RHS b has NaN/Inf!"<<std::endl; ++rhswarn;
	}
	// 对称 Jacobi 均衡化: As = D^{-1/2} A D^{-1/2}, 把刚性液相/气相数量级差
	// 巨大的行/列拉回 O(1), 使缩放后系统良态 -> SparseLU 稳定可解.
	Eigen::VectorXd d(dim);
	for(int k=0;k<dim;++k){ double akk=std::abs(A.coeff(k,k)); d[k]=(akk>1e-300)?1.0/std::sqrt(akk):1.0; }
	Eigen::SparseMatrix<double> As = A;
	for(int k=0;k<As.outerSize();++k)
		for(Eigen::SparseMatrix<double>::InnerIterator it(As,k); it; ++it)
			it.valueRef() *= d[it.row()]*d[it.col()];
	Eigen::VectorXd bs(dim); for(int k=0;k<dim;++k) bs[k]=b[k]*d[k];
	As.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
	lu.compute(As);
	if(lu.info()==Eigen::Success){
		Eigen::VectorXd ys = lu.solve(bs);
		if(lu.info()==Eigen::Success && ys.allFinite()){
			for(int k=0;k<dim;++k) x[k]=ys[k]*d[k];     // 反缩放
			return;
		}
	}
	// 退路: 缩放后系统上 GMRES + 对角预条件
	Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
	solver.set_restart(200); solver.setMaxIterations(5000); solver.setTolerance(1e-9);
	solver.compute(As);
	Eigen::VectorXd xg(dim); for(int k=0;k<dim;++k) xg[k]=x[k]/d[k];
	Eigen::VectorXd ys = solver.solveWithGuess(bs, xg);
	for(int k=0;k<dim;++k) x[k]=ys[k]*d[k];
	static int warn=0;
	if(solver.info()!=Eigen::Success && warn<3){
		std::cerr<<"[direct] SparseLU+GMRES fallback iters="<<solver.iterations()<<" err="<<solver.error()<<std::endl; ++warn;
	}
}

void PressureSubsystem(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, double h, const int rks, const EquationOfState& eos1, const EquationOfState& eos2){
	Array shock_eta(nxCells,nyCells);
	double Bx_star, By_star, Bz_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe,eos1, eos2);
	ucon2upri(Qs_pri, Qs,eos1, eos2);
	Eigen::VectorXd b_p_vect(2*nxCells*nyCells), p_vect(2 * nxCells * nyCells);
	
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	double p_inf_1 = eos1.getPInf();
	double p_inf_2 = eos2.getPInf();
	double gamma_1 = eos1.getGamma();
	double gamma_2 = eos2.getGamma();
	
	ShockDetector(shock_eta, Qe,eos1, eos2);

	bound_cond_u(u_con);

	Array rhovx_star(nxCells, nyCells);
	Array rhovy_star(nxCells, nyCells);
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			rhovx_star(i,j) = u_con(i,j,3);
			rhovy_star(i,j) = u_con(i,j,4);
		}
	}
	bound_cond_array(rhovx_star);
	bound_cond_array(rhovy_star);
	
	// (ρ_k h_k)^n —— 相焓, 不含 α^{n+1}/ρ^{n+1} 因子
	auto Hk = [&](int ii,int jj,int phase)->double{
		double cc[NUM_VARS], pp[NUM_VARS];
		for(int k=0;k<NUM_VARS;k++) cc[k]=Qe(ii,jj,k);
		ini_con2pri(pp,cc,eos1,eos2);
		double rr  = cc[1]+cc[2];
		double uu = cc[3]/rr, vv = cc[4]/rr;
		double ke = 0.5*(uu*uu+vv*vv);
		double arhoe = cc[5+phase] - cc[1+phase]*ke;        // (α_kρ_ke_k)^n
		double an    = phase==0 ? cc[0] : 1.0-cc[0];
		return (arhoe + pp[5+phase]) / fmax(an,ALPHA_MIN);   // (ρ_kh_k)^n
	};

	// ---------------- 焓通量散度的 path-consistent 离散 -----------------------
	//   ∇·(α_k^{n+1} H_k u*) = H_k (u*·∇α_k) + α_k ∇·(H_k u*)         (乘积法则)
	// 第一项含 α 梯度 -> 必须和显式 path-conservative 步用同一算子, 即增量
	//   (u·∇α1)^PC = (α1^{Qs} - α1^{n+1})/h,  且 (u·∇α2)^PC = -(u·∇α1)^PC
	// 第二项不含 α 梯度 -> 中心差分 (接触上 H 均匀且 ∇·u=0, 恒为 0).
	// 与 (6g) 展开式的 Σ 配对后, 接触间断上 p^{n+1} = p^n 精确成立 —— 这就是
	// 原来靠 Rcorr 硬凑出来的 WB, 现在由离散本身给出.
	// ------------- 焞通量散度: PC 一致的**守恒面通量** -----------------
	//   F^E_{k,i+1/2} = Hhat_{k,i+1/2} * Phi_{k,i+1/2}
	// Phi_k = 显式 Path_RusanovC 用的同一个 α_k Rusanov 面通量: 同重构态
	// (w_xR(i), w_xL(i+1))、同路径平均速度 uhat、同 S = vx_max.
	//   Phi_1 = 0.5*uhat*(a1L+a1R) - 0.5*S*(a1R-a1L),   Phi_2 = uhat - Phi_1
	// 两相相加后仍是 (Hhat_1 Phi_1 + Hhat_2 Phi_2) 的面差 -> 逐面对消, 混合总能
	// 严格守恒. 而当 uhat 均匀(= 接触间断)时恒等式
	//   (Phi_p - Phi_m)/dx == (D^- + D^+ + B)_0/dx == (α1^{Qs} - α1^{n+1})/h
	// 成立 (差为 (uhat_i-uhat_p)a_R - (uhat_i-uhat_m)a_L), 所以它与 Σ 里那个
	// 非守恒乘积算子在接触上逐位重合 -> WB 精确, 同时不丢守恒性.
	Array PhiX(nxCells,nyCells), UhX(nxCells,nyCells);
	Array PhiY(nxCells,nyCells), UhY(nxCells,nyCells);
	for(int i = ghost_cells-1; i <= nxCells + ghost_cells - 1; i++){
		for(int j = ghost_cells-1; j <= nyCells + ghost_cells - 1; j++){
			double QL[NUM_VARS], QR[NUM_VARS];
			for(int k=0;k<NUM_VARS;k++){ QL[k]=w_xR(i,j,k); QR[k]=w_xL(i+1,j,k); }
			double uh = pathAvgVel(QL,QR,3);
			UhX(i,j)  = uh;
			PhiX(i,j) = 0.5*uh*(QL[0]+QR[0]) - 0.5*vx_max*(QR[0]-QL[0]);
			for(int k=0;k<NUM_VARS;k++){ QL[k]=w_yR(i,j,k); QR[k]=w_yL(i,j+1,k); }
			double vh = pathAvgVel(QL,QR,4);
			UhY(i,j)  = vh;
			PhiY(i,j) = 0.5*vh*(QL[0]+QR[0]) - 0.5*vy_max*(QR[0]-QL[0]);
		}
	}
	// 速度修正: 焞通量里该用的是 u* (显式步后), 不是路径平均的 uhat^n.
	// 差额以**面量**形式补上: Psi_k = Phi_k + ahat_k*(u*_face - uhat_face),
	// 仍然是面通量 -> 仍然逐面对消; 而接触间断上 u* == uhat 精确成立
	// (row 1+2 给 u*dRho, row 3 给 u*u*dRho, 所以 (rho u)^{n+1} = u*rho^{n+1}),
	// 所以修正项在那里恰好为零, WB 不受影响.
	Array AhX(nxCells,nyCells), AhY(nxCells,nyCells);
	Array dUX(nxCells,nyCells), dUY(nxCells,nyCells);
	for(int i = ghost_cells-1; i <= nxCells + ghost_cells - 1; i++){
		for(int j = ghost_cells-1; j <= nyCells + ghost_cells - 1; j++){
			AhX(i,j) = 0.5*( w_xR(i,j,0) + w_xL(i+1,j,0) );
			AhY(i,j) = 0.5*( w_yR(i,j,0) + w_yL(i,j+1,0) );
			double us_i  = rhovx_star(i,j)  /(u_con(i,j,1)  +u_con(i,j,2));
			double us_ip = rhovx_star(i+1,j)/(u_con(i+1,j,1)+u_con(i+1,j,2));
			double vs_j  = rhovy_star(i,j)  /(u_con(i,j,1)  +u_con(i,j,2));
			double vs_jp = rhovy_star(i,j+1)/(u_con(i,j+1,1)+u_con(i,j+1,2));
			dUX(i,j) = 0.5*(us_i+us_ip) - UhX(i,j);
			dUY(i,j) = 0.5*(vs_j+vs_jp) - UhY(i,j);
		}
	}
	auto PhiXk = [&](int ii,int jj,int phase)->double{
		double a = phase==0 ? AhX(ii,jj) : 1.0-AhX(ii,jj);
		double f = phase==0 ? PhiX(ii,jj) : UhX(ii,jj)-PhiX(ii,jj);
		return f + a*dUX(ii,jj);
	};
	auto PhiYk = [&](int ii,int jj,int phase)->double{
		double a = phase==0 ? AhY(ii,jj) : 1.0-AhY(ii,jj);
		double f = phase==0 ? PhiY(ii,jj) : UhY(ii,jj)-PhiY(ii,jj);
		return f + a*dUY(ii,jj);
	};

	auto EnthDiv = [&](int ii,int jj,int phase)->double{
		double Hc = Hk(ii,jj,phase);
		double HxP = 0.5*(Hc + Hk(ii+1,jj,phase)), HxM = 0.5*(Hc + Hk(ii-1,jj,phase));
		double HyP = 0.5*(Hc + Hk(ii,jj+1,phase)), HyM = 0.5*(Hc + Hk(ii,jj-1,phase));
		return ( HxP*PhiXk(ii,jj,phase) - HxM*PhiXk(ii-1,jj,phase) )/dx
		     + ( HyP*PhiYk(ii,jj,phase) - HyM*PhiYk(ii,jj-1,phase) )/dy;
	};
	vector<T> triplets_p;
	SpMat M_p(nxCells * nyCells, nxCells * nyCells);
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
			
			double alpha1 = u_con(i,j,0);
			double alpha2 = 1 - alpha1;
						
			double p_star_1 = u_con(i,j,5);   // = (α1E1)*; 常数由 -α1·B1 在 RHS 处理(通用)
			double p_star_2 = u_con(i,j,6);   // = (α2E2)*
			
			double rho = u_con(i,j,1) + u_con(i,j,2);
			double rho_p0 = u_con(i+1,j,1) + u_con(i+1,j,2);
			double rho_m0 = u_con(i-1,j,1) + u_con(i-1,j,2);
			double rho_0p = u_con(i,j+1,1) + u_con(i,j+1,2);
			double rho_0m = u_con(i,j-1,1) + u_con(i,j-1,2);
			
			double h1_00, h2_00, h1_p0, h2_p0, h1_m0, h2_m0, h1_0p, h2_0p,h1_0m, h2_0m;
			calEnthalpy_phasic(con00_n, h1_00,h2_00,eos1, eos2);
			calEnthalpy_phasic(conp0_n, h1_p0,h2_p0,eos1, eos2);
			calEnthalpy_phasic(conm0_n, h1_m0,h2_m0,eos1, eos2);
			calEnthalpy_phasic(con0p_n, h1_0p,h2_0p,eos1, eos2);
			calEnthalpy_phasic(con0m_n, h1_0m,h2_0m,eos1, eos2);
			
			double A1,B1,A2,B2;
			get_phase_linear_coeffs(u_con(i,j,1) / alpha1, eos1,A1, B1);
			get_phase_linear_coeffs(u_con(i,j,2) / alpha2, eos2,A2, B2);
			
			double RHS_1 = p_star_1 - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j)) * u_con(i,j,1) / (rho * rho) - h * EnthDiv(i,j,0) - alpha1 * B1;   // const = α1·B1 (通用)
			
			double RHS_2 = p_star_2 - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j)) * u_con(i,j,2) / (rho * rho) - h * EnthDiv(i,j,1) - alpha2 * B2;   // const = α2·B2 (通用)
			
			
			b_p_vect(M_P(i,j,0)) = RHS_1;
			b_p_vect(M_P(i,j,1)) = RHS_2;

			p_vect(M_P(i,j,0)) = Qe_pri(i,j,5);
			p_vect(M_P(i,j,1)) = Qe_pri(i,j,6);
		}
	}

	PressureSolve_Direct(p_vect, b_p_vect, u_con, Qe, Qs, shock_eta, h, eos1, eos2, 2 * nxCells*nyCells);
	
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			double pmix_00 = p_vect(M_P(i,  j,0)) + p_vect(M_P(i,  j,1));
			double pmix_m0 = p_vect(M_P(i-1,j,0)) + p_vect(M_P(i-1,j,1));
			double pmix_p0 = p_vect(M_P(i+1,j,0)) + p_vect(M_P(i+1,j,1));
			double pmix_0m = p_vect(M_P(i,j - 1,0)) + p_vect(M_P(i,j-1,1));
			double pmix_0p = p_vect(M_P(i,j+1,0)) + p_vect(M_P(i,j+1,1));
			
			double fL, fR,	flux_xm, flux_xp;
			fL = pmix_m0;
			fR = pmix_00;
			flux_xm = 0.5 * (fL + fR);
			
			fL = pmix_00;
			fR = pmix_p0;
			flux_xp = 0.5 * (fL + fR);
			
			u_con(i,j,3) = u_con(i,j,3) - h / dx * (flux_xp - flux_xm);
			
			double gL, gR, flux_ym, flux_yp;
			gL = pmix_0m;
			gR = pmix_00;
			flux_ym = 0.5 * (gL + gR);
			
			gL = pmix_00;
			gR = pmix_0p;
			flux_yp = 0.5 * (gL + gR);
			
			u_con(i,j,4) = u_con(i,j,4) - h / dy * (flux_yp - flux_ym);

		}
	}

	bound_cond_u(u_con);
	
	//std::cout <<"Second Pressure Solved "<< "rk: " << rks + 1 << endl;
	vector<T> triplets_p_second;
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
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
			
			double alpha1 = u_con(i,j,0);
			double alpha2 = 1 - alpha1;
			
			double p_star_1 = u_con(i,j,5);   // = (α1E1)*; 常数由 -α1·B1 在 RHS 处理(通用)
			double p_star_2 = u_con(i,j,6);   // = (α2E2)*
			
			double rho = u_con(i,j,1) + u_con(i,j,2);
			double rho_p0 = u_con(i+1,j,1) + u_con(i+1,j,2);
			double rho_m0 = u_con(i-1,j,1) + u_con(i-1,j,2);
			double rho_0p = u_con(i,j+1,1) + u_con(i,j+1,2);
			double rho_0m = u_con(i,j-1,1) + u_con(i,j-1,2);
			
			double h1_00, h2_00, h1_p0, h2_p0, h1_m0, h2_m0, h1_0p, h2_0p,h1_0m, h2_0m;
			calEnthalpy_phasic(con00_n, h1_00,h2_00,eos1, eos2);
			calEnthalpy_phasic(conp0_n, h1_p0,h2_p0,eos1, eos2);
			calEnthalpy_phasic(conm0_n, h1_m0,h2_m0,eos1, eos2);
			calEnthalpy_phasic(con0p_n, h1_0p,h2_0p,eos1, eos2);
			calEnthalpy_phasic(con0m_n, h1_0m,h2_0m,eos1, eos2);
			
			double A1,B1,A2,B2;
			get_phase_linear_coeffs(u_con(i,j,1) / alpha1, eos1,A1, B1);
			get_phase_linear_coeffs(u_con(i,j,2) / alpha2, eos2,A2, B2);
			
			double RHS_1 = p_star_1 - 0.5 * (u_con(i,j,3) * u_con(i,j,3) + u_con(i,j,4) * u_con(i,j,4)) * u_con(i,j,1) / (rho * rho) - h * EnthDiv(i,j,0) - alpha1 * B1;   // const = α1·B1 (通用)
			
			double RHS_2 = p_star_2 - 0.5 * (u_con(i,j,3) * u_con(i,j,3) + u_con(i,j,4) * u_con(i,j,4)) * u_con(i,j,2) / (rho * rho) - h * EnthDiv(i,j,1) - alpha2 * B2;   // const = α2·B2 (通用)
			
			
			b_p_vect(M_P(i,j,0)) = RHS_1;
			b_p_vect(M_P(i,j,1)) = RHS_2;
			
			p_vect(M_P(i,j,0)) = Qe_pri(i,j,5);
			p_vect(M_P(i,j,1)) = Qe_pri(i,j,6);
			
		}
	}
	PressureSolve_Direct(p_vect, b_p_vect, u_con, Qe, Qs, shock_eta, h, eos1, eos2, 2 * nxCells*nyCells);
	
	for(int i=ghost_cells;i<ghost_cells+nxCells;i++){
		for(int j=ghost_cells;j<ghost_cells+nyCells;j++){
			double alpha1 = u_con(i,j,0);
			double alpha2 = 1 - alpha1;
			double rho = u_con(i,j,1)+u_con(i,j,2);
			double uu = u_con(i,j,3)/rho, vv = u_con(i,j,4)/rho;
			double ke = 0.5*(uu*uu+vv*vv);
			double a1 = u_con(i,j,0), a2 = 1-a1;
			double a1p1 = p_vect(M_P(i,j,0)), a2p2 = p_vect(M_P(i,j,1));
			// α_kρ_ke_k = (α_kp_k + α_kγπ)/(γ-1);  α_kE_k = 内能 + 动能
			double A1,B1,A2,B2;
			get_phase_linear_coeffs(u_con(i,j,1) / alpha1, eos1,A1, B1);
			get_phase_linear_coeffs(u_con(i,j,2) / alpha2, eos2,A2, B2);
			u_con(i,j,5) = eos1.getAlphaE_fromAlphaP(a1p1, u_con(i,j,1), uu, vv, alpha1);
			u_con(i,j,6) = eos2.getAlphaE_fromAlphaP(a2p2, u_con(i,j,2), uu, vv, alpha2);
		}
	}
}

// 幻影相重置: alpha_k < eps_v 的格子,弱相压力强制等于主相压力
// 守恒: 相质量、动量、混合内能(总能)不变
void reset_vanishing_phase(Grid& u_con, const EquationOfState& eos1, const EquationOfState& eos2){
	const double eps_v = 1e-4;
	double g1 = eos1.getGamma(), g2 = eos2.getGamma();
	double pi1 = eos1.getPInf(), pi2 = eos2.getPInf();
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double a1 = u_con(i,j,0), a2 = 1.0 - a1;
			if(a1 >= eps_v && a2 >= eps_v) continue;
			double m1 = u_con(i,j,1), m2 = u_con(i,j,2);
			double rho = m1 + m2;
			double uu = u_con(i,j,3)/rho, vv = u_con(i,j,4)/rho;
			double ke = 0.5*(uu*uu + vv*vv);
			double E1 = u_con(i,j,5) - m1*ke;   // alpha1 rho1 e1
			double E2 = u_con(i,j,6) - m2*ke;
			double Eint = E1 + E2;
			if(a1 < eps_v){
				// 主相 = 2: 先用主相内能定 p2, 再令 p1 = p2
				double p2 = ((g2-1.0)*E2 - a2*g2*pi2)/a2;
				double p1 = p2;
				double E1n = (a1*p1 + a1*g1*pi1)/(g1-1.0);
				double E2n = Eint - E1n;          // 混合内能守恒
				u_con(i,j,5) = E1n + m1*ke;
				u_con(i,j,6) = E2n + m2*ke;
			}else{
				double p1 = ((g1-1.0)*E1 - a1*g1*pi1)/a1;
				double p2 = p1;
				double E2n = (a2*p2 + a2*g2*pi2)/(g2-1.0);
				double E1n = Eint - E2n;
				u_con(i,j,5) = E1n + m1*ke;
				u_con(i,j,6) = E2n + m2*ke;
			}
		}
	}
}

// 压力正性 floor: 刚性 EOS 的稀疏波可把相压力拉穿内聚极限 p_k < -p∞,k
// (p+p∞<0 -> 声速虚 -> 失稳). 无相变的力学模型必须把每相压力 clamp 到
// p_k+p∞,k ≥ eps. 通用: 用 getAlphaP_fromAlphaE / getAlphaE_fromAlphaP 往返.
void clamp_pressure_floor(Grid& u_con, const EquationOfState& eos1, const EquationOfState& eos2){
	double pi1=eos1.getPInf(), pi2=eos2.getPInf();
	double eps1=fmax(1e2, 1e-6*pi1), eps2=fmax(1e2, 1e-6*pi2);
	for(int i=ghost_cells;i<nxCells+ghost_cells;i++) for(int j=ghost_cells;j<nyCells+ghost_cells;j++){
		double a1=u_con(i,j,0), a2=1.0-a1;
		double m1=u_con(i,j,1), m2=u_con(i,j,2), rho=m1+m2;
		double u=u_con(i,j,3)/rho, v=u_con(i,j,4)/rho;
		double ap1=eos1.getAlphaP_fromAlphaE(u_con(i,j,5),m1,u,v,a1);
		double ap2=eos2.getAlphaP_fromAlphaE(u_con(i,j,6),m2,u,v,a2);
		double p1=ap1/a1, p2=ap2/a2;
		if(!(p1+pi1>=eps1)){ p1=-pi1+eps1; u_con(i,j,5)=eos1.getAlphaE_fromAlphaP(a1*p1,m1,u,v,a1); }
		if(!(p2+pi2>=eps2)){ p2=-pi2+eps2; u_con(i,j,6)=eos2.getAlphaE_fromAlphaP(a2*p2,m2,u,v,a2); }
	}
}

// ============================================================================
//  精确压力松弛 (Pelanti-Shyue 6.pdf Appendix C.1)
//  瞬时力学平衡 p1=p2=p*，mixture-energy-consistent
//  守恒不变: 相质量 con[1,2]、动量 con[3,4]、混合总能 con[5]+con[6]
//  布局: con[5]=α1E1, con[6]=α2E2 (相总能)
// ============================================================================
void PressureRelaxation_Exact(Grid& u_con,
	const EquationOfState& eos1, const EquationOfState& eos2){
		double g1 = eos1.getGamma(),  g2 = eos2.getGamma();
		double pi1 = eos1.getPInf(),  pi2 = eos2.getPInf();
		
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				
				double a1 = u_con(i,j,0);
				double a2 = 1.0 - a1;
				if(a1 < 1e-4 || a2 < 1e-4) continue;  // 近纯相由 reset_vanishing_phase 处理
				
				double m1 = u_con(i,j,1), m2 = u_con(i,j,2);   // 相质量(冻结)
				double rho = m1 + m2;
				double u = u_con(i,j,3)/rho, v = u_con(i,j,4)/rho;
				double ke = 0.5*(u*u + v*v);                   // 比动能(冻结)
				
				// 相内能密度 (相总能 - 相动能)
				double E1 = u_con(i,j,5) - m1*ke;   // α1ρ1e1
				double E2 = u_con(i,j,6) - m2*ke;   // α2ρ2e2
				double Eint = E1 + E2;              // 混合内能(松弛中守恒)
				
				// 旧相压力 (SG: αp = (γ-1)αρe - αγπ)
				double p1 = ((g1-1.0)*E1 - a1*g1*pi1)/a1;
				double p2 = ((g2-1.0)*E2 - a2*g2*pi2)/a2;
				
				// ---- C.1 二次方程 a p*² + b p* + d = 0 ----
				double C1 = 2.0*g1*pi1 + (g1-1.0)*p1;
				double C2 = 2.0*g2*pi2 + (g2-1.0)*p2;
				double A = 1.0 + g2*a1 + g1*a2;
				double B = C1*a2 + C2*a1 - (1.0+g2)*a1*p1 - (1.0+g1)*a2*p2;
				double D = -(C2*a1*p1 + C1*a2*p2);
				double disc = B*B - 4.0*A*D;
				double pstar0 = (-B + sqrt(fmax(disc, 0.0))) / (2.0*A);   // 物理正根
				
				// ---- C.2 平衡体积分数 α1* ----
				double a1s = a1 * ((g1-1.0)*pstar0 + 2.0*p1 + C1)
				/ ((g1+1.0)*pstar0 + C1);
				a1s = std::max(ALPHA_MIN, std::min(a1s, 1.0 - ALPHA_MIN));
				double a2s = 1.0 - a1s;
				
				// ---- 用混合 EOS 重算 p* (保证内能精确守恒) ----
				//  Eint = p* κ + Π,  κ=Σα_k/(γ_k-1), Π=Σα_k γ_k π_k/(γ_k-1)
				double kap = a1s/(g1-1.0) + a2s/(g2-1.0);
				double Pi  = a1s*g1*pi1/(g1-1.0) + a2s*g2*pi2/(g2-1.0);
				double pstar = (Eint - Pi)/kap;
				
				// ---- 重分配相内能 (两相同 p*，和=Eint 精确) ----
				double E1n = a1s*(pstar + g1*pi1)/(g1-1.0);   // α1ρ1e1*
				double E2n = a2s*(pstar + g2*pi2)/(g2-1.0);   // α2ρ2e2*
				
				// ---- 写回 ----
				u_con(i,j,0) = a1s;
				u_con(i,j,5) = E1n + m1*ke;   // α1E1
				u_con(i,j,6) = E2n + m2*ke;   // α2E2
				// 相质量 con[1,2]、动量 con[3,4] 不变
			}
		}
	}

void dump_profile(const char* fn, Grid& u_con, const EquationOfState& eos1, const EquationOfState& eos2){
	FILE* f=fopen(fn,"w");
	double pri[NUM_VARS],con[NUM_VARS]; int jj=ghost_cells;
	for(int i=ghost_cells;i<nxCells+ghost_cells;i++){
		for(int k=0;k<NUM_VARS;k++)con[k]=u_con(i,jj,k);
		ini_con2pri(pri,con,eos1,eos2);
		fprintf(f,"%g %g %g %g %g %g %g %g\n", x_0+dx*(i-ghost_cells+0.5),
			pri[0],pri[1],pri[2],pri[3],pri[5],pri[6],pri[5]+pri[6]);
	}
	fclose(f);
}
void diag_print(const char* tag, Grid& u_con, const EquationOfState& eos1, const EquationOfState& eos2){
	double pmin=1e300,pmax=-1e300,vmin=1e300,vmax=-1e300,p1min=1e300,p1max=-1e300,p2min=1e300,p2max=-1e300;
	double pri[NUM_VARS],con[NUM_VARS];
	int jj=ghost_cells;
	for(int i=ghost_cells;i<nxCells+ghost_cells;i++){
		for(int k=0;k<NUM_VARS;k++)con[k]=u_con(i,jj,k);
		ini_con2pri(pri,con,eos1,eos2);
		double pm=pri[5]+pri[6];
		double p1=pri[5]/pri[0], p2=pri[6]/(1-pri[0]);
		double v=pri[3];
		pmin=fmin(pmin,pm);pmax=fmax(pmax,pm);
		vmin=fmin(vmin,v);vmax=fmax(vmax,v);
		p1min=fmin(p1min,p1);p1max=fmax(p1max,p1);
		p2min=fmin(p2min,p2);p2max=fmax(p2max,p2);
	}
	printf("%-12s pmix:[%.6e,%.6e] dv:%.3e  p1:[%.4e,%.4e] p2:[%.4e,%.4e]\n",tag,pmin,pmax,vmax-vmin,p1min,p1max,p2min,p2max);
}

int main(int argc, char *argv[]){
	double gamma1,gamma2,pinf1,pinf2,a1,b1;
	EquationOfState* eos1_ptr = nullptr;
	EquationOfState* eos2_ptr = nullptr;
	
	switch (test) {
		case testshockbubble:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.0014;
			gamma1   = 1.4;
			gamma2   = 1.67;
			pinf1    = 0.0;
			pinf2    = 0.0;
			eos1_ptr = new StiffenedGasEOS(gamma1, pinf1);
			eos2_ptr = new StiffenedGasEOS(gamma2, pinf2);
			break;
		case test2dahshockbubble:
			// 物理域范围（mm）
			x_0 = 0.0;   x_1 = 325.0;
			y_0 = 0.0;   y_1 = 44.5;
			z_0 = -0.5;  z_1 = 0.5;   // 2D 问题，z 方向薄层
			time_stop = 280;        // 根据物理时间调整，可能需要更长时间（例如 0.002？）
			gamma1 = 1.4;
			gamma2 = 1.67;
			pinf1 = 0.0;
			pinf2 = 0.0;
			eos1_ptr = new StiffenedGasEOS(gamma1, pinf1);
			eos2_ptr = new StiffenedGasEOS(gamma2, pinf2);
			break;
		case testDodecaneShockTube:
			// 设置域范围
			x_0 = 0.0; x_1 = 1.0;
			y_0 = 0.0; y_1 = 1.0;  // 1D
			z_0 = 0.0; z_1 = 1.0;
			time_stop = 473e-6;     // 473 微秒
			// 创建EOS
			eos1_ptr = new StiffenedGasEOS(1.025, 0.0);       // 蒸气
			eos2_ptr = new StiffenedGasEOS(2.35, 4e8);        // 液体
			break;
	}
	
	EquationOfState& eos1 = *eos1_ptr;
	EquationOfState& eos2 = *eos2_ptr;
	dx = (x_1 - x_0)/(nxPoints-1);
	dy = (y_1 - y_0)/(nyPoints-1);
	dz = (z_1 - z_0)/(nzPoints-1);

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells), u_con_n(nxCells, nyCells), Qe(nxCells, nyCells), Qs(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells),w_yL(nxCells, nyCells), w_yR(nxCells, nyCells);
	Array Bx_p0(nxCells,nyCells), By_0p(nxCells,nyCells);
	u_initialization(u_pri, Bx_p0, By_0p);
	upri2ucon(u_con, u_pri,eos1, eos2);
	ucon2upri(u_pri, u_con,eos1, eos2);
	std::vector<double> bRK(RKs);
	Array Ae(RKs, RKs), Ai(RKs, RKs);
	Tensor Hs(RKs, nxCells, nyCells);
	Array shock_eta(nxCells, nyCells);
	
	
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
	
	
	do{
		timestep++;
		bound_cond_u(u_pri);

		
		dt = fmin(caldt(u_pri,eos1,eos2), time_stop - t);
		t = t + dt;
		//std::cout << B_max << " " << U_max << endl;
		
		upri2ucon(u_con, u_pri, eos1, eos2);
		
		copy(u_con_n,u_con);
		
		for(int rks = 0; rks < RKs; rks++){
			copy(Qe,u_con_n);
			copy(Qs,u_con_n);
			
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

			adt =  Ai(rks,rks) * dt;

			// vanishing-phase 保护: 显式增量可能把 Qe/Qs 的消失相质量推成微负
			// (负相质量 -> 后续相压力/焓 NaN). clamp α 与相质量为正下限.
			auto clampQ=[&](Grid&Q,int i,int j){
				double a=Q(i,j,0); if(!std::isfinite(a)) a=0.5; Q(i,j,0)=fmax(ALPHA_MIN,fmin(a,1.0-ALPHA_MIN));
				double m1=Q(i,j,1),m2=Q(i,j,2);
				double rr=( std::isfinite(m1)?m1:0.0)+(std::isfinite(m2)?m2:0.0);
				double mf=ALPHA_MIN*fmax(rr,1e-30);
				if(!(std::isfinite(m1)&&m1>=mf)) Q(i,j,1)=mf;
				if(!(std::isfinite(m2)&&m2>=mf)) Q(i,j,2)=mf;
			};
			for(int i=0;i<nxCells+2*ghost_cells;i++) for(int j=0;j<nyCells+2*ghost_cells;j++){ clampQ(Qe,i,j); clampQ(Qs,i,j); }

			TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR, eos1, eos2);
			
			Path_RusanovC_Update(u_con, Qs, w_xL,w_xR,w_yL,w_yR,eos1,eos2);
			bound_cond_u(u_con);

			// vanishing-phase 保护: 刚性接触上显式通量可能把消失相质量推成微负,
			// 负相质量 -> 负相密度 -> getAlphaP/calcs 产 NaN. 在压力解之前 clamp
			// α 与相质量为正下限, 并把消失相压力重置为力学平衡(对主导相可忽略).
			for(int i=ghost_cells;i<nxCells+ghost_cells;i++) for(int j=ghost_cells;j<nyCells+ghost_cells;j++){
				double a=u_con(i,j,0); if(!std::isfinite(a)) a=0.5; u_con(i,j,0)=fmax(ALPHA_MIN,fmin(a,1.0-ALPHA_MIN));
				double m1=u_con(i,j,1),m2=u_con(i,j,2);
				double rr=(std::isfinite(m1)?m1:0.0)+(std::isfinite(m2)?m2:0.0);
				double mfloor=ALPHA_MIN*fmax(rr,1e-30);
				if(!(std::isfinite(m1)&&m1>=mfloor)) u_con(i,j,1)=mfloor;
				if(!(std::isfinite(m2)&&m2>=mfloor)) u_con(i,j,2)=mfloor;
			}
			bound_cond_u(u_con);
			reset_vanishing_phase(u_con, eos1, eos2);
			bound_cond_u(u_con);

			PressureSubsystem(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, Bx_p0, By_0p,adt,rks, eos1, eos2);
			clamp_pressure_floor(u_con, eos1, eos2);
			bound_cond_u(u_con);


			reset_vanishing_phase(u_con, eos1, eos2);
			PressureRelaxation_Boscheri(u_con, mu, eos1, eos2);
			clamp_pressure_floor(u_con, eos1, eos2);
			bound_cond_u(u_con);
			
			ucon2upri(u_pri, u_con,eos1, eos2);
			
			for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
				for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						Hs(rks,i,j,k) = (u_con(i,j,k) - Qs(i,j,k)) / adt;
					}
				}
			}
			bound_cond_u(u_con);
		}
		

		bound_cond_u(u_con);
		
		ucon2upri(u_pri, u_con,eos1, eos2);

		std::cout<<"[timestep, dt, time] = "  << timestep << " "<<dt << " " << t << endl;
}while(t < time_stop);
	
	ofstream outFileshock_de("/Users/chenzefeng/Desktop/shock_de.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileshock_de << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << shock_eta(i,j) << endl;
		}
	}
	outFileshock_de.close();
	
	ofstream outFilealpha("/Users/chenzefeng/Desktop/alpha.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilealpha << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,0) << endl;
		}
	}
	outFilealpha.close();
	
	ofstream outFilerho("/Users/chenzefeng/Desktop/6RK_rho.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerho << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,1) + u_con(i,j,2) << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx("/Users/chenzefeng/Desktop/6RK_vx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,3) << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy("/Users/chenzefeng/Desktop/6RK_vy.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,4) << endl;
		}
	}
	outFilevy.close();
	
	
	ofstream outFilep("/Users/chenzefeng/Desktop/6RK_p.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilep << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,5) + u_pri(i,j,6)) << endl;
		}
	}
	outFilep.close();
	
	ofstream outFilerhoe1("/Users/chenzefeng/Desktop/6RK_rhoe1.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerhoe1 << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,5) << endl;
		}
	}
	outFilerhoe1.close();
	
	ofstream outFilerhoe2("/Users/chenzefeng/Desktop/6RK_rhoe2.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerhoe2 << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,6) << endl;
		}
	}
	outFilerhoe2.close();
	
	
	
	ofstream outFiledata("/Users/chenzefeng/Desktop/IMEX_data.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFiledata << x_0 + dx * (i - ghost_cells + 0.5) << " "<< y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,0) << " " << u_pri(i,j,1)<< " " << u_pri(i,j,2)<< " " << u_pri(i,j,3)<< " " << u_pri(i,j,4)<< " " << u_pri(i,j,5)<< " " << u_pri(i,j,6)<< " " << (u_pri(i,j,5)+u_pri(i,j,6))<< endl;
		}
		//		}
	}
	outFiledata.close();
}