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



#define NUM_VARS 7
enum Formula{EnergyFormula, PressureFormula};
enum Test{testshockbubble};
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
Test test = testshockbubble;

const int nxCells = 1000;
const int nyCells = 1;
const int nzCells = 1;
const int ghost_cells = 3;

const int RKs = 1;
double CFL = 0.3;   // path-conservative 2nd-order MUSCL + IMEX 需要 CFL≲0.4，否则激波/接触处过冲
Reconstruction recon =First_Order;

bool add_Evisc = true;    // shock-activated chi-Laplacian on the mixture pressure (eq 53-54)
bool add_IEvisc = false;  // DIAGNOSTIC: IE Laplacian OFF

const double ALPHA_MIN = 1e-8;
const double mu = 1e30;
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
		virtual double getP(double rhoe) const = 0;
		virtual double getP_fromE(double E, double rho, double u, double v) const = 0;
		virtual double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,double u, double v, double alpha) const = 0;
		virtual	double getAlphaE_fromAlphaP(double alphaP, double alphaRho,double u, double v, double alpha) const = 0;
		virtual double getrhoe(double p) const = 0;
		virtual double getAlphaP(double alphaRhoe, double alpha) const = 0;
		virtual double internalEnergy(double rho, double p) const = 0;
		virtual double pressure(double rho, double e) const = 0;
		virtual double calcs(double rho, double p) const = 0;
		virtual double calcf_x(double rho, double p, double Bx, double By, double Bz) const = 0;
		virtual double calcf_y(double rho, double p, double Bx, double By, double Bz) const = 0;
		virtual void e_linear_coeff(double rho, double& c, double& d) const = 0;
		virtual EquationOfState* clone() const = 0;  // 可选：多态复制
};

//-----------------------------
// Stiffened Gas EOS
//-----------------------------
class StiffenedGasEOS : public EquationOfState {
private:
	const double gamma;
	const double pInf;
	
public:
	StiffenedGasEOS(double gamma_, double pInf_) : gamma(gamma_), pInf(pInf_) {}
	double getGamma() const override{ return gamma; }
	double getPInf() const override{ return pInf; }
	double getP(double rhoe) const override{ 
		return rhoe*(gamma-1) - gamma*pInf; 
	}
	// 传入总能密度 E,需要额外知道动能密度
	double getP_fromE(double E, double rho, double u, double v) const override{
		double ke = 0.5*rho*(u*u + v*v);     // 动能密度
		double rhoe = E - ke;                 // 内能密度
		return rhoe*(gamma-1) - gamma*pInf;
	}
	// 传入相总能密度 alphaRhoE = α_k ρ_k E_k,得到 α_k p_k
	// 需要: 相质量 alphaRho = α_k ρ_k, 速度 u,v, 体积分数 alpha
	double getAlphaP_fromAlphaE(double alphaRhoE, double alphaRho,
		double u, double v, double alpha) const override{
			double ke    = 0.5*alphaRho*(u*u + v*v);   // 相动能密度 ½ α_k ρ_k |u|²
			double alphaRhoe = alphaRhoE - ke;          // 相内能密度 α_k ρ_k e_k
			return alphaRhoe*(gamma-1) - alpha*gamma*pInf;   // α_k p_k
		}
	// 传入 α_k p_k,得到 α_k ρ_k E_k(相总能密度)
	// 需要: 相质量 alphaRho = α_k ρ_k, 速度 u,v, 体积分数 alpha
	double getAlphaE_fromAlphaP(double alphaP, double alphaRho,double u, double v, double alpha) const override{
			double alphaRhoe = (alphaP + alpha*gamma*pInf)/(gamma-1);  // 相内能密度
			double ke        = 0.5*alphaRho*(u*u + v*v);                // 相动能密度
			return alphaRhoe + ke;                                       // α_k ρ_k E_k
		}
	double getAlphaP(double alphaRhoe, double alpha) const override{
		return alphaRhoe*(gamma-1.0) - alpha*gamma*pInf;
	}
	double getrhoe(double p) const override{ 
		return (p) / (gamma - 1); 
	}
	double internalEnergy(double rho, double p) const override {
		return (p + gamma * pInf) / ((gamma - 1.0f) * rho);
	}
	
	double pressure(double rho, double e) const override {
		return (gamma - 1.0f) * rho * e - gamma * pInf;
	}
	
	double calcs(double rho, double p) const override {
		return std::sqrt(gamma * (p + pInf) / rho);
	}
	
	double calcf_x(double rho, double p, double Bx, double By, double Bz) const override {
		double B2 = Bx * Bx + By * By + Bz * Bz;
		return sqrt(0.5 * (((gamma * (p + pInf) + B2)/(rho)) + sqrt( ((gamma * (p + pInf) + B2) / rho) * ((gamma * (p + pInf) + B2) / rho) - 4 * ((gamma * (p + pInf) * Bx * Bx)/(rho * rho)))));
	}
	
	double calcf_y(double rho, double p, double Bx, double By, double Bz) const override {
		double B2 = Bx * Bx + By * By + Bz * Bz;
		return sqrt(0.5 * (((gamma * (p + pInf) + B2)/(rho)) + sqrt( ((gamma * (p + pInf) + B2) / rho) * ((gamma * (p + pInf) + B2) / rho) - 4 * ((gamma * (p + pInf) * By * By)/(rho * rho)))));
	}
	
	void e_linear_coeff(double rho, double& c, double& d) const override {
		c = 1.0f / ((gamma - 1.0f) * rho);
		d = (gamma * pInf) / ((gamma - 1.0f) * rho);
	}
	
	EquationOfState* clone() const override { return new StiffenedGasEOS(*this); }
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

	// ---- vanishing-phase guard: trace 相压力 slave 给在场相（只除在场相的大 α）----
	const double EPS = 1e-4;
	if(alpha1 < EPS){ double p2 = alpha2p2/alpha2; alpha1p1 = alpha1*p2; }  // 除 α2(≈1)
	if(alpha2 < EPS){ double p1 = alpha1p1/alpha1; alpha2p2 = alpha2*p1; }  // 除 α1(≈1)

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
	con[6] = eos2.getAlphaE_fromAlphaP(pri[6], con[2], vx, vy, alpha2);;
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
		for(int r = 0; r < NUM_VARS; r++){
			double Av = 0.0;
			for(int c = 0; c < NUM_VARS; c++) Av += Ah[r][c]*dQ[c];
			Dm[r] = 0.5*Av - 0.5*S*dQ[r];
			Dp[r] = 0.5*Av + 0.5*S*dQ[r];
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
		for(int r = 0; r < NUM_VARS; r++){
			double Av = 0.0;
			for(int c = 0; c < NUM_VARS; c++) Av += Ah[r][c]*dQ[c];
			Dm[r] = 0.5*Av - 0.5*S*dQ[r];
			Dp[r] = 0.5*Av + 0.5*S*dQ[r];
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
		
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
				for(int k = 0; k < NUM_VARS; k++){
					double slope_xR = u_pri(i + 1,j,k) - u_pri(i,j,k);
					double slope_xL = u_pri(i,j,k) - u_pri(i - 1,j,k);
					double slope_yR = u_pri(i,j + 1,k) - u_pri(i,j,k);
					double slope_yL = u_pri(i,j,k) - u_pri(i,j - 1,k);
					switch (recon) {
						case First_Order:
							wpri_xL(i,j,k) = u_pri(i,j,k);
							wpri_xR(i,j,k) = u_pri(i,j,k);
							wpri_yL(i,j,k) = u_pri(i,j,k);
							wpri_yR(i,j,k) = u_pri(i,j,k);
							break;
						case Second_Order:
							// x minmod
							if(slope_xL * slope_xR <= 0)            c = 0;
							else if(fabs(slope_xR) <= fabs(slope_xL)) c = slope_xR;
							else                                    c = slope_xL;
							wpri_xL(i,j,k) = u_pri(i,j,k) - 0.5 * c;
							wpri_xR(i,j,k) = u_pri(i,j,k) + 0.5 * c;
							// y minmod
							if(slope_yL * slope_yR <= 0)            c = 0;
							else if(fabs(slope_yR) <= fabs(slope_yL)) c = slope_yR;
							else                                    c = slope_yL;
							wpri_yL(i,j,k) = u_pri(i,j,k) - 0.5 * c;
							wpri_yR(i,j,k) = u_pri(i,j,k) + 0.5 * c;
							break;
					}
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
	const EquationOfState& eos2){
		double alpha1 = con[0];
		double alpha2 = 1.0 - alpha1;
		double rho1   = con[1] / alpha1;
		double rho2   = con[2] / alpha2;
		double rho    = con[1] + con[2];
		double u      = con[3] / rho;
		double v      = con[4] / rho;

		// phasic pressures from phasic TOTAL energies (con[5],con[6]) —
		// getAlphaP_fromAlphaE subtracts the kinetic part before applying SG EOS
		double p1 = eos1.getAlphaP_fromAlphaE(con[5], con[1], u, v, alpha1) / alpha1;
		double p2 = eos2.getAlphaP_fromAlphaE(con[6], con[2], u, v, alpha2) / alpha2;
		
		// phasic sound speeds:  c_k^2 = gamma_k (p_k + pinf_k) / rho_k
		double g1 = eos1.getGamma(), pinf1 = eos1.getPInf();
		double g2 = eos2.getGamma(), pinf2 = eos2.getPInf();
		double c1sq = g1 * (p1 + pinf1) / rho1;
		double c2sq = g2 * (p2 + pinf2) / rho2;
		
		// mass fractions
		double Y1 = con[1] / rho;   // = alpha1*rho1/rho
		double Y2 = con[2] / rho;
		
		double cf2 = Y1 * c1sq + Y2 * c2sq;
		return sqrt(fmax(cf2, 1e-30));
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
	const EquationOfState& eos1, const EquationOfState& eos2){
		double alpha1 = Q[0];
		double alpha2 = 1.0 - alpha1;
		double rho = Q[1] + Q[2];                       // = α1ρ1 + α2ρ2，O(1)
		double u = Q[3]/rho, v = Q[4]/rho;
		double a1p1 = eos1.getAlphaP_fromAlphaE(Q[5], Q[1], u, v, alpha1);  //
		double a2p2 = eos2.getAlphaP_fromAlphaE(Q[6], Q[2], u, v, alpha2);
		double g1 = eos1.getGamma(), pi1 = eos1.getPInf();
		double g2 = eos2.getGamma(), pi2 = eos2.getPInf();
		double rhoc2 = g1*(a1p1 + alpha1*pi1) + g2*(a2p2 + alpha2*pi2);
		return sqrt(fmax(rhoc2 / rho, 0.0));            // fmax 防 metastable
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
		
		// 相压力密度 α_k p_k(SG,不除 α)
		double ap1 = eos1.getAlphaP(arhoe1, alpha1);   // α₁p₁
		double ap2 = eos2.getAlphaP(arhoe2, alpha2);   // α₂p₂
		
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

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x,Grid& u_con, Grid& Qe, Grid& Qs, Array &chi,double h,const EquationOfState& eos1, const EquationOfState& eos2){
		const double gamma1 = eos1.getGamma();
		const double gamma2 = eos2.getGamma();
		const double pinf_1 = eos1.getPInf();
		const double pinf_2 = eos2.getPInf();
		
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
				
				// 算sigma
				double rho = u_con(i,j,1)+u_con(i,j,2);
				double Y1 = u_con(i,j,1)/rho;
				double Y2 = u_con(i,j,2)/rho;
				double alpha1 = u_con(i,j,0);
				double alpha2 = 1 - alpha1;
				double A1_n0 = (prip0_n[5] - prim0_n[5]) / (2 * dx);
				double A2_n0 = (pri0p_n[5] - pri0m_n[5]) / (2 * dy);
				double A3_n0 = (prip0_n[6] - prim0_n[6]) / (2 * dx);
				double A4_n0 = (pri0p_n[6] - pri0m_n[6]) / (2 * dy);
				double A1_n1 = (x(M_P(i + 1,j,0)) - x(M_P(i - 1,j,0))) / (2 * dx);
				double A2_n1 = (x(M_P(i,j + 1,0)) - x(M_P(i,j - 1,0))) / (2 * dy);
				double A3_n1 = (x(M_P(i + 1,j,1)) - x(M_P(i - 1,j,1))) / (2 * dx);
				double A4_n1 = (x(M_P(i,j + 1,1)) - x(M_P(i,j - 1,1))) / (2 * dy);
				double p_mix = x(M_P(i,j,0)) + x(M_P(i,j,1));
				double p_mix_p0 = x(M_P(i+1,j,0)) + x(M_P(i+1,j,1));
				double p_mix_m0 = x(M_P(i-1,j,0)) + x(M_P(i-1,j,1));
				double p_mix_0p = x(M_P(i,j+1,0)) + x(M_P(i,j+1,1));
				double p_mix_0m = x(M_P(i,j-1,0)) + x(M_P(i,j-1,1));
				double sigma_1 = -Y2*(u_con(i,j,3) / rho * A1_n1 - h / rho * 1 / (2*dx) * (p_mix_p0 - p_mix_m0) * A1_n0 + u_con(i,j,4) / rho * A2_n1 - h / rho * 1 / (2*dy) * (p_mix_0p - p_mix_0m) * A2_n0);
				double sigma_2 = Y1*(u_con(i,j,3) / rho * A3_n1 - h / rho * 1 / (2*dx) * (p_mix_p0 - p_mix_m0) * A3_n0 + u_con(i,j,4) / rho * A4_n1 - h / rho * 1 / (2*dy) * (p_mix_0p - p_mix_0m) * A4_n0);
				double sigma = sigma_1 + sigma_2;
				
				
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
				double rho_ratio1_00 = (u_con(i,j,1) / u_con(i,j,0)) / (u_con(i,j,1) + u_con(i,j,2));
				double rho_ratio1_p0 = (u_con(i + 1,j,1) / u_con(i + 1,j,0)) / (u_con(i + 1,j,1) + u_con(i + 1,j,2));
				double rho_ratio1_m0 = (u_con(i - 1,j,1) / u_con(i - 1,j,0)) / (u_con(i - 1,j,1) + u_con(i - 1,j,2));
				double rho_ratio1_0p = (u_con(i,j + 1,1) / u_con(i,j + 1,0)) / (u_con(i,j + 1,1) + u_con(i,j + 1,2));
				double rho_ratio1_0m = (u_con(i,j - 1,1) / u_con(i,j - 1,0)) / (u_con(i,j - 1,1) + u_con(i,j - 1,2));
				rho_ratio1_p0 = 0.5 * (rho_ratio1_00 + rho_ratio1_p0);
				rho_ratio1_m0 = 0.5 * (rho_ratio1_00 + rho_ratio1_m0);
				rho_ratio1_0p = 0.5 * (rho_ratio1_00 + rho_ratio1_0p);
				rho_ratio1_0m = 0.5 * (rho_ratio1_00 + rho_ratio1_0m);
				double rho_ratio2_00 = (u_con(i,j,2) / u_con(i,j,0)) / (u_con(i,j,1) + u_con(i,j,2));
				double rho_ratio2_p0 = (u_con(i + 1,j,2) / u_con(i + 1,j,0)) / (u_con(i + 1,j,1) + u_con(i + 1,j,2));
				double rho_ratio2_m0 = (u_con(i - 1,j,2) / u_con(i - 1,j,0)) / (u_con(i - 1,j,1) + u_con(i - 1,j,2));
				double rho_ratio2_0p = (u_con(i,j + 1,2) / u_con(i,j + 1,0)) / (u_con(i,j + 1,1) + u_con(i,j + 1,2));
				double rho_ratio2_0m = (u_con(i,j - 1,2) / u_con(i,j - 1,0)) / (u_con(i,j - 1,1) + u_con(i,j - 1,2));
				rho_ratio2_p0 = 0.5 * (rho_ratio2_00 + rho_ratio2_p0);
				rho_ratio2_m0 = 0.5 * (rho_ratio2_00 + rho_ratio2_m0);
				rho_ratio2_0p = 0.5 * (rho_ratio2_00 + rho_ratio2_0p);
				rho_ratio2_0m = 0.5 * (rho_ratio2_00 + rho_ratio2_0m);
				
				double alpha1_p0 = u_con(i + 1,j,0);
				double alpha1_m0 = u_con(i - 1,j,0);
				double alpha1_0p = u_con(i,j + 1,0);
				double alpha1_0m = u_con(i,j - 1,0);
				alpha1_p0 = 0.5 * (alpha1_p0 + alpha1);
				alpha1_m0 = 0.5 * (alpha1_m0 + alpha1);
				alpha1_0p = 0.5 * (alpha1_0p + alpha1);
				alpha1_0m = 0.5 * (alpha1_0m + alpha1);

				double alpha2_p0 = 1 - u_con(i + 1,j,0);
				double alpha2_m0 = 1 - u_con(i - 1,j,0);
				double alpha2_0p = 1 - u_con(i,j + 1,0);
				double alpha2_0m = 1 - u_con(i,j - 1,0);
				alpha2_p0 = 0.5 * (alpha2_p0 + alpha2);
				alpha2_m0 = 0.5 * (alpha2_m0 + alpha2);
				alpha2_0p = 0.5 * (alpha2_0p + alpha2);
				alpha2_0m = 0.5 * (alpha2_0m + alpha2);
				
				Ax(M_P(i,j,0)) = x(M_P(i,j,0)) / (gamma1 - 1) - h * h * ( 1 / (dx * dx) * (rho_ratio1_p0 * alpha1_p0 * h1_p0 * (p_mix_p0 - p_mix) - rho_ratio1_m0 * alpha1_m0 * h1_m0 * (p_mix - p_mix_m0)) + 1 / (dy * dy) * (rho_ratio1_0p * alpha1_0p * h1_0p * (p_mix_0p - p_mix) - rho_ratio1_0m * alpha1_0m * h1_0m * (p_mix - p_mix_0m)) ) + h * sigma ;
				Ax(M_P(i,j,1)) = x(M_P(i,j,1)) / (gamma2 - 1) - h * h * ( 1 / (dx * dx) * (rho_ratio2_p0 * alpha2_p0 * h2_p0 * (p_mix_p0 - p_mix) - rho_ratio2_m0 * alpha2_m0 * h2_m0 * (p_mix - p_mix_m0)) + 1 / (dy * dy) * (rho_ratio2_0p * alpha2_0p * h2_0p * (p_mix_0p - p_mix) - rho_ratio2_0m * alpha2_0m * h2_0m * (p_mix - p_mix_0m)) ) - h * sigma ;
				
				// ===== shock-detector Laplacian: -Δt/(γ_k-1)·∇·(κ ∇(α_k p_k)) =====
				// χ 已由 ShockDetector 预算 (Array chi, 含一层 ghost); κ 用对流速度
				if(add_Evisc){

				// 各格 |速度|（分母用各格混合密度）


				// κ at cells (x 用 |u|·dx, y 用 |v|·dy)
					double kx_00 = chi(i,  j  ) * vx_max * dx;
					double kx_p0 = chi(i+1,j  ) * vx_max * dx;
					double kx_m0 = chi(i-1,j  ) * vx_max * dx;
					double ky_00 = chi(i,  j  ) * vy_max * dy;
					double ky_0p = chi(i,  j+1) * vy_max * dy;
					double ky_0m = chi(i,  j-1) * vy_max * dy;
					
					// face-average κ
					double kxR = 0.5*(kx_00 + kx_p0), kxL = 0.5*(kx_m0 + kx_00);
					double kyT = 0.5*(ky_00 + ky_0p), kyB = 0.5*(ky_0m + ky_00);
					
					// q_k = α_k p_k 邻居值（未知量）
					double q1_00 = x(M_P(i,  j,  0)), q1_p0 = x(M_P(i+1,j,  0)), q1_m0 = x(M_P(i-1,j,  0));
					double q1_0p = x(M_P(i,  j+1,0)), q1_0m = x(M_P(i,  j-1,0));
					double q2_00 = x(M_P(i,  j,  1)), q2_p0 = x(M_P(i+1,j,  1)), q2_m0 = x(M_P(i-1,j,  1));
					double q2_0p = x(M_P(i,  j+1,1)), q2_0m = x(M_P(i,  j-1,1));
					
					// ∇·(κ∇q_k) 守恒型离散
					double lap_shock1 = 1.0/(dx*dx)*( kxR*(q1_p0-q1_00) - kxL*(q1_00-q1_m0) ) + 1.0/(dy*dy)*( kyT*(q1_0p-q1_00) - kyB*(q1_00-q1_0m) );
					double lap_shock2 = 1.0/(dx*dx)*( kxR*(q2_p0-q2_00) - kxL*(q2_00-q2_m0) ) + 1.0/(dy*dy)*( kyT*(q2_0p-q2_00) - kyB*(q2_00-q2_0m) );
					
					Ax(M_P(i,j,0)) += -h/(gamma1-1.0) * lap_shock1;   // 相1
					Ax(M_P(i,j,1)) += -h/(gamma2-1.0) * lap_shock2;   // 相2
				}   // if(add_Evisc)
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
		double g1 = eos1.getGamma(),  g2 = eos2.getGamma();
		double pinf1 = eos1.getPInf(), pinf2 = eos2.getPInf();
		
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
				
				// 相压力(SG: αp = (γ-1)αρe - αγp∞)
				double a1p1 = (g1-1.0)*arhoe1 - a1H*g1*pinf1;
				double a2p2 = (g2-1.0)*arhoe2 - a2H*g2*pinf2;
				double p1H  = a1p1/a1H;
				double p2H  = a2p2/a2H;
				
				// ---- 声速、阻抗、界面压力(sqrt 加 floor 防 p+p∞<0 出 NaN)----
				double c1 = sqrt(fmax(g1*(p1H + pinf1)/rho1H, 1e-12));
				double c2 = sqrt(fmax(g2*(p2H + pinf2)/rho2H, 1e-12));
				double Z1 = rho1H*c1, Z2 = rho2H*c2;
				double pI = (Z2*p1H + Z1*p2H)/(Z1 + Z2);
				
				// ---- 弛豫斜率 K_k = dp_k/dα_k(含 p_I 做功;SG: ρc²=γ(p+p∞))----
				double K1 = -(1.0/a1H)*( g1*(p1H+pinf1) + (g1-1.0)*(pI - p1H) );
				double K2 = -(1.0/a2H)*( g2*(p2H+pinf2) + (g2-1.0)*(pI - p2H) );
				
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
		
		Eigen::VectorXd Ax(dim), r(dim), w(dim);
		Eigen::MatrixXd V(dim, m + 1);
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m + 1, m);
		Eigen::VectorXd cs(m), sn(m), g(m + 1), y(m);
		
		for (int outer = 0; outer < GMRESiter; ++outer){
			matrix_vector_product(Ax, x, u_con, Qe, Qs, chi, h, eos1, eos2);
			r = b - Ax;
			double beta = r.norm();
			if (beta < GMREStol) return;
			
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
				if (std::abs(g(k + 1)) < GMREStol){ ++k; break; }
			}
			// back-substitution: H(0:k,0:k) y = g(0:k)
			for (int i = k - 1; i >= 0; --i){
				double s0 = g(i);
				for (int j = i + 1; j < k; ++j) s0 -= H(i, j) * y(j);
				y(i) = s0 / H(i, i);
			}
			for (int i = 0; i < k; ++i) x.noalias() += y(i) * V.col(i);
			
			if (std::abs(g(k)) < GMREStol) return;   // converged
		}
		
		matrix_vector_product(Ax, x, u_con, Qe, Qs, chi, h, eos1, eos2);
		std::cerr << "Warning: GMRES not converged, |r|/|b| = "
		<< ((b - Ax).norm() / std::max(norm_b, 1e-30)) << std::endl;
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
						
			double p_star_1 = -(alpha1 * gamma_1 * p_inf_1) / (gamma_1 - 1) + u_con(i,j,5);
			double p_star_2 = -(alpha2 * gamma_2 * p_inf_2) / (gamma_2 - 1) + u_con(i,j,6);
			
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
			
			
			double RHS_1 = p_star_1 - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j)) * u_con(i,j,1) / (rho * rho) - h * ( 1/dx * (0.5* (h1_p0 + h1_00) * 0.5 * ( rhovx_star(i+1,j)+rhovx_star(i,j)) * 0.5 * (u_con(i + 1,j,1) / rho_p0 + u_con(i,j,1) / rho)  - 0.5*(h1_m0+h1_00)*0.5*(rhovx_star(i-1,j)+rhovx_star(i,j))  * 0.5 * (u_con(i - 1,j,1) / rho_m0 + u_con(i,j,1) / rho)) + 1/dy * ( 0.5*(h1_0p+h1_00)*0.5*(rhovy_star(i,j+1)+rhovy_star(i,j))*0.5 * (u_con(i,j + 1,1) / rho_0p + u_con(i,j,1) / rho) - 0.5*(h1_0m+h1_00)*0.5*(rhovy_star(i,j-1)+rhovy_star(i,j)) *0.5 * (u_con(i,j - 1,1) / rho_0m + u_con(i,j,1) / rho) )) ;
			
			double RHS_2 = p_star_2 - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j)) * u_con(i,j,2) / (rho * rho) - h * ( 1/dx * (0.5* (h2_p0 + h2_00) * 0.5 * ( rhovx_star(i+1,j)+rhovx_star(i,j)) * 0.5 * (u_con(i + 1,j,2) / rho_p0 + u_con(i,j,2) / rho)  - 0.5*(h2_m0+h2_00)*0.5*(rhovx_star(i-1,j)+rhovx_star(i,j))  * 0.5 * (u_con(i - 1,j,2) / rho_m0 + u_con(i,j,2) / rho)) + 1/dy * ( 0.5*(h2_0p+h2_00)*0.5*(rhovy_star(i,j+1)+rhovy_star(i,j))*0.5 * (u_con(i,j + 1,2) / rho_0p + u_con(i,j,2) / rho) - 0.5*(h2_0m+h2_00)*0.5*(rhovy_star(i,j-1)+rhovy_star(i,j)) *0.5 * (u_con(i,j - 1,2) / rho_0m + u_con(i,j,2) / rho) )) ;
			
			
			b_p_vect(M_P(i,j,0)) = RHS_1;
			b_p_vect(M_P(i,j,1)) = RHS_2;
			
			p_vect(M_P(i,j,0)) = Qe_pri(i,j,5);
			p_vect(M_P(i,j,1)) = Qe_pri(i,j,6);
		}
	}

	GMRES_Solver(p_vect, b_p_vect, Pressure_Subsystem_AX, u_con, Qe, Qs, shock_eta,h, eos1,eos2,2 * nxCells*nyCells);
	
	
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
			
			double p_star_1 = -(alpha1 * gamma_1 * p_inf_1) / (gamma_1 - 1) + u_con(i,j,5);
			double p_star_2 = -(alpha2 * gamma_2 * p_inf_2) / (gamma_2 - 1) + u_con(i,j,6);
			
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
			
			
			double RHS_1 = p_star_1 - 0.5 * (u_con(i,j,3) * u_con(i,j,3) + u_con(i,j,4) * u_con(i,j,4)) * u_con(i,j,1) / (rho * rho) - h * ( 1/dx * (0.5* (h1_p0 + h1_00) * 0.5 * ( rhovx_star(i+1,j)+rhovx_star(i,j)) * 0.5 * (u_con(i + 1,j,1) / rho_p0 + u_con(i,j,1) / rho)  - 0.5*(h1_m0+h1_00)*0.5*(rhovx_star(i-1,j)+rhovx_star(i,j))  * 0.5 * (u_con(i - 1,j,1) / rho_m0 + u_con(i,j,1) / rho)) + 1/dy * ( 0.5*(h1_0p+h1_00)*0.5*(rhovy_star(i,j+1)+rhovy_star(i,j))*0.5 * (u_con(i,j + 1,1) / rho_0p + u_con(i,j,1) / rho) - 0.5*(h1_0m+h1_00)*0.5*(rhovy_star(i,j-1)+rhovy_star(i,j)) *0.5 * (u_con(i,j - 1,1) / rho_0m + u_con(i,j,1) / rho) )) ;
			
			double RHS_2 = p_star_2 - 0.5 * (u_con(i,j,3) * u_con(i,j,3) + u_con(i,j,4) * u_con(i,j,4)) * u_con(i,j,2) / (rho * rho) - h * ( 1/dx * (0.5* (h2_p0 + h2_00) * 0.5 * ( rhovx_star(i+1,j)+rhovx_star(i,j)) * 0.5 * (u_con(i + 1,j,2) / rho_p0 + u_con(i,j,2) / rho)  - 0.5*(h2_m0+h2_00)*0.5*(rhovx_star(i-1,j)+rhovx_star(i,j))  * 0.5 * (u_con(i - 1,j,2) / rho_m0 + u_con(i,j,2) / rho)) + 1/dy * ( 0.5*(h2_0p+h2_00)*0.5*(rhovy_star(i,j+1)+rhovy_star(i,j))*0.5 * (u_con(i,j + 1,2) / rho_0p + u_con(i,j,2) / rho) - 0.5*(h2_0m+h2_00)*0.5*(rhovy_star(i,j-1)+rhovy_star(i,j)) *0.5 * (u_con(i,j - 1,2) / rho_0m + u_con(i,j,2) / rho) )) ;
			
			
			b_p_vect(M_P(i,j,0)) = RHS_1;
			b_p_vect(M_P(i,j,1)) = RHS_2;
			
			p_vect(M_P(i,j,0)) = Qe_pri(i,j,5);
			p_vect(M_P(i,j,1)) = Qe_pri(i,j,6);
			
		}
	}
	GMRES_Solver(p_vect, b_p_vect, Pressure_Subsystem_AX, u_con, Qe, Qs, shock_eta,h, eos1,eos2,2 * nxCells*nyCells);
	
	for(int i=ghost_cells;i<ghost_cells+nxCells;i++){
		for(int j=ghost_cells;j<ghost_cells+nyCells;j++){
			double rho = u_con(i,j,1)+u_con(i,j,2);
			double uu = u_con(i,j,3)/rho, vv = u_con(i,j,4)/rho;
			double ke = 0.5*(uu*uu+vv*vv);
			double a1 = u_con(i,j,0), a2 = 1-a1;
			double a1p1 = p_vect(M_P(i,j,0)), a2p2 = p_vect(M_P(i,j,1));
			// α_kρ_ke_k = (α_kp_k + α_kγπ)/(γ-1);  α_kE_k = 内能 + 动能
			u_con(i,j,5) = (a1p1 + a1*gamma_1*p_inf_1)/(gamma_1-1) + u_con(i,j,1)*ke;
			u_con(i,j,6) = (a2p2 + a2*gamma_2*p_inf_2)/(gamma_2-1) + u_con(i,j,2)*ke;
		}
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
				if(a1 < 10.0*ALPHA_MIN || a2 < 10.0*ALPHA_MIN) continue;  // 近纯相跳过
				
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

// ============================================================================
//  vanishing-phase 清理：trace 相(α<EPS)的密度/压力 slave 给在场相
//  布局: con[5]=α1E1, con[6]=α2E2 (相总能, 含动能)
//  与 sixeq_boscheri_1d_ref.cpp 的 reset_vanishing_phase 等价(适配 2D 相总能)
// ============================================================================
void reset_vanishing_phase(Grid& u_con, const EquationOfState& eos1, const EquationOfState& eos2){
    const double g1 = eos1.getGamma(), pinf1 = eos1.getPInf();
    const double g2 = eos2.getGamma(), pinf2 = eos2.getPInf();
    const double EPS = 1e-4;
    for(int i = 0; i < nxCells + 2*ghost_cells; i++){
        for(int j = 0; j < nyCells + 2*ghost_cells; j++){
            double a1 = fmin(fmax(u_con(i,j,0), ALPHA_MIN), 1.0 - ALPHA_MIN);
            double a2 = 1.0 - a1;
            double m1 = u_con(i,j,1), m2 = u_con(i,j,2);
            double rho = m1 + m2;
            double u = u_con(i,j,3)/rho, v = u_con(i,j,4)/rho;
            double ke = 0.5*(u*u + v*v);
            double r1 = m1/a1, r2 = m2/a2;
            // 相内能密度 α_kρ_ke_k = 相总能 - 相动能
            double arhoe1 = u_con(i,j,5) - m1*ke;
            double arhoe2 = u_con(i,j,6) - m2*ke;
            // 相压力 (SG: αp = (γ-1)αρe - αγπ)
            double p1 = ((g1-1.0)*arhoe1 - a1*g1*pinf1)/a1;
            double p2 = ((g2-1.0)*arhoe2 - a2*g2*pinf2)/a2;

            if(a1 < EPS){            // phase 1 是 trace → 取齐 phase 2
                r1 = r2;  p1 = p2;
                m1 = a1*r1;  u_con(i,j,1) = m1;
                u_con(i,j,5) = a1*(p1 + g1*pinf1)/(g1-1.0) + m1*ke;
            }
            if(a2 < EPS){            // phase 2 是 trace → 取齐 phase 1
                r2 = r1;  p2 = p1;
                m2 = a2*r2;  u_con(i,j,2) = m2;
                u_con(i,j,6) = a2*(p2 + g2*pinf2)/(g2-1.0) + m2*ke;
            }
            u_con(i,j,0) = a1;
        }
    }
}

int main(int argc, char *argv[]){
	double gamma1,gamma2,pinf1,pinf2,a1,b1;
	switch (test) {
		case testshockbubble:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.0014;
			gamma1   = 1.4;
			gamma2   = 1.67;
			pinf1    = 0.0;
			pinf2    = 0.0;
			break;
	}
	
	StiffenedGasEOS eos1(gamma1, pinf1);
	StiffenedGasEOS eos2(gamma2, pinf2);
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

		reset_vanishing_phase(u_con, eos1, eos2);   // trace 相密度/压力 slave，防 vanishing-phase 噪声
		bound_cond_u(u_con);

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
			

			TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR, eos1, eos2);
			
			Path_RusanovC_Update(u_con, Qs, w_xL,w_xR,w_yL,w_yR,eos1,eos2);
			bound_cond_u(u_con);
			
			PressureSubsystem(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, Bx_p0, By_0p,adt,rks, eos1, eos2);
			bound_cond_u(u_con);
			
//			PressureRelaxation_Boscheri(u_con, mu, eos1, eos2);
			PressureRelaxation_Exact(u_con, eos1, eos2);
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

		// ---- 界面诊断:激波进氦时 dump x∈[0.40,0.56] 的 α1,ρ,u,p1,p2 ----
		static double next_dump = 0.0006;
		if(t >= next_dump){
			int j = ghost_cells;
			std::cerr << "### t="<<t<<" step="<<timestep<<"  x  a1  rho  u  p1  p2\n";
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				double x = x_0 + (i - ghost_cells + 0.5)*dx;
				if(x < 0.40 || x > 0.58) continue;
				double a1=u_pri(i,j,0), rho=a1*u_pri(i,j,1)+(1-a1)*u_pri(i,j,2);
				double p1=u_pri(i,j,5)/a1, p2=u_pri(i,j,6)/(1-a1);
				fprintf(stderr,"  %.4f  %.3e  %.4f  %8.2f  %9.1f  %9.1f\n",x,a1,rho,u_pri(i,j,3),p1,p2);
			}
			next_dump += 0.0003;
		}

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