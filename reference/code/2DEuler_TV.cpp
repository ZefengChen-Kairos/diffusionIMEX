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


#define NUM_VARS 4
enum Formula{EnergyFormula, PressureFormula};
enum Test{testToroModified,testToro1,testToro2,testToro3,testToro4,testToro5, test2d1, test2d2, test2d3, test2d4,testExplosion, testSmoothAcousticWave,  testDoubleShearLayer, testGreshoVortex, testWoodwardColellaBlastWave,testCollidingPulses,testRiemannLike, convergencetest, KH_instability, RPContact, RPLax, RPTwoShocks, RPTwoRare, RPFVVEM1, RPFVVEM2, RPFVVEM3, RPFVVEM4};
enum Reconstruction{NoRecon, Minmod, VanLeer, ENO3, WENO5};
enum shock_detector{Yes_add, No_add};
enum SolvePressure{no_correction, once_correction, all_correction};
enum dtCalculation{Fullydt, IMEXdt};
enum x_boundary_condition{x_Transmissive, x_Reflective, x_Periodic, x_Dirichlet};
enum y_boundary_condition{y_Transmissive, y_Reflective, y_Periodic, y_Dirichlet};
enum Add_Pvisc{MaterialSpeedp, Noaddpvis};
enum Energy_Update{ConservationofU, DefinitionofU};
enum Flux_Diss{Local_Diss, Global_Diss};


double vortex_Mach = 1e-1;

dtCalculation caltimestep = IMEXdt;
Test test = test2d1;
Formula formula = PressureFormula;
Flux_Diss flux_diss = Global_Diss;
x_boundary_condition x_bc;
y_boundary_condition y_bc;
SolvePressure PressureSolver = all_correction;
const int nxCells = 256;
const int nyCells = 256;
const int nzCells = 1;
const int ghost_cells = 5;

const int RKs = 2;

double CFL = 0.3;
int SpaceOrder;
Reconstruction recon = Minmod;
const int full_step = 10;
Add_Pvisc add_pvisc = MaterialSpeedp;
shock_detector add_shock_detector = Yes_add;
double shock_kappa1 = 1e-3;

double GMREStol = 1e-10;
double GMRESiter = 1e16;

// Energy Update
Energy_Update U_update = DefinitionofU;

int nxPoints = nxCells + 1, nyPoints = nyCells + 1, nzPoints = nzCells + 1;
double t, time_stop, dt, adt, adt2;
int timestep = 0;

double x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1, B_max_x = 0, B_max_y = 0, U_max_x = 0,U_max_y = 0,a_max_x = 0,a_max_y = 0, lambda = 0, lambda_x = 0, lambda_y = 0, fullydt = 0, imexdt = 0, p_ref = 1, p_epsilon = 1;
//double t, dt, time_stop = 0.2,x_0 = -1.0,y_0 = -1.0,x_1 = 1.0, y_1 = 1.0, CFL= 0.8;

double dx = (x_1 - x_0)/(nxPoints-1);
double dy = (y_1 - y_0)/(nyPoints-1);
double dz = (z_1 - z_0)/(nzPoints-1);
double gamma_plasma = 1.4;
double ch_x,ch_y, cp_x, cp_y, cA_x, cA_y, vx_max, vy_max;




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
void copystate(double *state_1, double *state_0);
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void ini_con2pri(double *pri, double *con);
void ini_pri2con(double *con, double *pri);
void upri2ucon(Grid &u_con, Grid &u_pri);
void ucon2upri(Grid &u_pri,Grid &u_con);
double calecon(double *con);
double calepri(double *pri);
void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p);
void bound_cond_u(Grid &u);
double caldt(Grid &u_pri);
void update(Grid &u_con,Grid &flux, char dir);
void store(Grid &u_con, Array& vxn,Array& vyn, Array& Bxn, Array& Byn);
void bound_cond_array(Array &u, x_boundary_condition xbc, y_boundary_condition ybc );
double calkinetic(double *con);
void ShockDetector(Array &shock_eta, Grid &u_con);
double calenthalpy(double *con);
void Roe_Flux(double *flux, double *uL, double *uR);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri, Array& eps_M, Array& eps_M_x,Array& eps_M_y,const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, Array& eps_M,Array& eps_M_x,Array& eps_M_y,int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime);



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
	double U = con[3];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[3];
	double e = calecon(con);
	double k = calkinetic(con);
	f[0] = rho * vx;
	f[1] = rho * vx * vx;
	f[2] = rho * vx * vy;
	f[3] = k * vx;
}

void f_pressure(double *f, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[3];
	double e = calecon(con);
	double h = calenthalpy(con);
	f[0] = 0;
	f[1] = p;
	f[2] = 0;
	f[3] = con[1] * h;
}


void f_full_con(double *f, double *con){
	double f_c[NUM_VARS], f_p[NUM_VARS];
	f_convective(f_c, con);
	f_pressure(f_p, con);
	for(int k = 0; k < NUM_VARS; k++){
		f[k] = f_c[k] + f_p[k];
	}
}

void g_convective(double *g, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double U = con[3];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[3];
	double e = calecon(con);
	double k = calkinetic(con);
	g[0] = rho * vy;
	g[1] = rho * vy * vx;
	g[2] = rho * vy * vy;
	g[3] = k * vy;
}

void g_pressure(double *g, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[3];
	double e = calecon(con);
	double h = calenthalpy(con);
	g[0] = 0;
	g[1] = 0;
	g[2] = p;
	g[3] = con[2] * h;
}


void g_full_con(double *g, double *con){
	double g_c[NUM_VARS], g_p[NUM_VARS];
	g_convective(g_c, con);
	g_pressure(g_p, con);
	for(int k = 0; k < NUM_VARS; k++){
		g[k] = g_c[k] + g_p[k];
	}
}



void ini_con2pri(double *pri, double *con){
	// convert conservative (rho rhov E) to primitive (rho v p)
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = (gamma_plasma - 1) * con[0] * e;
}

void ini_pri2con(double *con, double *pri){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	double e = calepri(pri);
	// rho rhov E
	double vx = pri[1];
	double vy = pri[2];
	double v_square = vx * vx + vy * vy;
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * pri[2];
	con[3] = pri[0] * e + 0.5 * pri[0] * v_square;
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

void RusanovC_x(double *flux, double *QL, double *QR){
	double fL_c[NUM_VARS], fR_c[NUM_VARS];
	double fL_p[NUM_VARS], fR_p[NUM_VARS];
	double fL[NUM_VARS],   fR[NUM_VARS];
	
	f_convective(fL_c, QL);
	f_convective(fR_c, QR);
	f_pressure  (fL_p, QL);
	f_pressure  (fR_p, QR);
	
	
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL_c[k] + fR_c[k]) - 0.5 * vx_max * (QR[k] - QL[k]);
	}
}

void RusanovC_y(double *flux, double *QL, double *QR){
	double gL_c[NUM_VARS], gR_c[NUM_VARS];
	double gL_p[NUM_VARS], gR_p[NUM_VARS];
	double gL[NUM_VARS],   gR[NUM_VARS];
	
	g_convective(gL_c, QL);
	g_convective(gR_c, QR);
	g_pressure  (gL_p, QL);
	g_pressure  (gR_p, QR);
	
		
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (gL_c[k] + gR_c[k]) - 0.5 * vy_max * (QR[k] - QL[k]);
	}
}


void Update(Grid &u_con, Grid &Qe, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_cellL[NUM_VARS], Q_cellR[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
				Q_cellL[k] = Qe(i,j,k); 
				Q_cellR[k] = Qe(i+1,j,k);
				//				std::cout <<  Q_xL[k] << " " << Q_xR[k] << endl;
			}
			
			RusanovC_x(fp,Q_xL, Q_xR);
								
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
				Q_cellL[k] = Qe(i - 1,j,k); 
				Q_cellR[k] = Qe(i,j,k);
				//				std::cout <<  Q_xL[k] << " " << Q_xR[k] << endl;
			}
			
			RusanovC_x(fm,Q_xL, Q_xR);
					
			
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				Q_cellL[k] = Qe(i,j,k); 
				Q_cellR[k] = Qe(i,j + 1,k);

				//				std::cout <<  Q_yL[k] << " " << Q_yR[k] << endl;
			}
			RusanovC_y(gp,Q_yL, Q_yR);
					
							
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
				Q_cellL[k] = Qe(i,j - 1,k); 
				Q_cellR[k] = Qe(i,j,k);
				//				std::cout <<  Q_yL[k] << " " << Q_yR[k] << endl;
			}
			
			
			RusanovC_y(gm,Q_yL, Q_yR);
					
			
			
			for(int k = 0; k < NUM_VARS; k++){
				u_con(i,j,k) = Qs(i,j,k) - adt / dx * (fp[k] - fm[k]) - adt / dy * (gp[k] - gm[k]);
			}
		}
	}
}



void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR){
	double c;
	double q0, q1, q2,q0L, q1L, q2L,q0R, q1R, q2R, wr0, wr1, wr2, wr0_bar, wr1_bar, wr2_bar, beta0, beta1, beta2, d0, d1, d2;
	double diff1, diff2, diff3, diff4;
	Grid u_pri(nxCells,nyCells);
	ucon2upri(u_pri, u_con);
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
			for(int k = 0; k < NUM_VARS; k++){
				double slope_xR = u_pri(i + 1,j,k) - u_pri(i,j,k);
				double slope_xL = u_pri(i,j,k) - u_pri(i - 1,j,k);
				
				double slope_yR = u_pri(i,j + 1,k) - u_pri(i,j,k);
				double slope_yL = u_pri(i,j,k) - u_pri(i,j - 1,k);
				switch (recon) {
					case NoRecon:
						SpaceOrder = 1;
						w_xL(i,j,k) = u_pri(i,j,k);
						w_xR(i,j,k) = u_pri(i,j,k);
						w_yL(i,j,k) = u_pri(i,j,k);
						w_yR(i,j,k) = u_pri(i,j,k);
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
						
						if(test == convergencetest){
							c = 0.5 * (slope_xL + slope_xR);
						}
						
						w_xL(i,j,k) = u_pri(i,j,k) - 0.5 * c;
						w_xR(i,j,k) = u_pri(i,j,k) + 0.5 * c;
						
						if(slope_yL * slope_yR <= 0){
							c = 0;
						}else if(fabs(slope_yR) <= fabs(slope_yL)){
							c = slope_yR;
						}else if(fabs(slope_yL) < fabs(slope_yR)){
							c = slope_yL;
						}
						
						if(test == convergencetest){
							c = 0.5 * (slope_yL + slope_yR);
						}
						w_yL(i,j,k) = u_pri(i,j,k) - 0.5 * c;
						w_yR(i,j,k) = u_pri(i,j,k) + 0.5 * c;
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
						w_xL(i,j,k) = u_pri(i,j,k) - 0.5 * c * slope_xL;
						w_xR(i,j,k) = u_pri(i,j,k) + 0.5 * c * slope_xR;
						
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
						w_yL(i,j,k) = u_pri(i,j,k) - 0.5 * c * slope_yL;
						w_yR(i,j,k) = u_pri(i,j,k) + 0.5 * c * slope_yR;
						break;
					case ENO3:
						SpaceOrder = 3;
						q0L = -1.0/6.0 * u_pri(i - 2, j, k) + 5.0/6.0 * u_pri(i - 1, j, k) + 1.0/3.0 * u_pri(i, j, k);
						q1L = 1.0/3.0 * u_pri(i - 1, j, k) + 5.0/6.0 * u_pri(i, j, k) - 1.0/6.0 * u_pri(i + 1, j, k);
						q2L = 11.0/6.0 * u_pri(i, j, k) - 7.0/6.0 * u_pri(i + 1, j, k) + 1.0/3.0 * u_pri(i + 2, j, k);
						
						
						q0R = 1.0/3.0 * u_pri(i - 2, j, k) - 7.0/6.0 * u_pri(i - 1, j, k) + 11.0/6.0 * u_pri(i, j, k);
						q1R = -1.0/6.0 * u_pri(i - 1, j, k) + 5.0/6.0 * u_pri(i, j, k) + 1.0/3.0 * u_pri(i + 1, j, k);
						q2R = 1.0/3.0 * u_pri(i, j, k) + 5.0/6.0 * u_pri(i + 1, j, k) - 1.0/6.0 * u_pri(i + 2, j, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u_pri(i, j, k) - u_pri(i - 1, j, k));
						diff2 = fabs(u_pri(i + 1, j, k) - u_pri(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u_pri(i + 2, j, k) - u_pri(i + 1, j, k)) - fabs(u_pri(i + 1, j, k) - u_pri(i, j, k)));
							diff4 = fabs(abs(u_pri(i + 1, j, k) - u_pri(i, j, k)) - fabs(u_pri(i, j, k) - u_pri(i - 1, j, k)));
							if (diff3 >= diff4){
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}else{
								w_xL(i, j, k) = q2L;
								w_xR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u_pri(i, j, k) - u_pri(i - 1, j, k)) - fabs(u_pri(i - 1, j, k) - u_pri(i - 2, j, k)));
							diff4 = fabs(abs(u_pri(i + 1, j, k) - u_pri(i, j, k)) - fabs(u_pri(i, j, k) - u_pri(i - 1, j, k)));
							if (diff4 >= diff3){
								w_xL(i, j, k) = q0L;
								w_xR(i, j, k) = q0R;
							}else{
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}
						}
						
						q0L = -1.0/6.0 * u_pri(i, j - 2, k) + 5.0/6.0 * u_pri(i, j - 1, k) + 1.0/3.0 * u_pri(i, j, k);
						q1L = 1.0/3.0 * u_pri(i, j - 1, k) + 5.0/6.0 * u_pri(i, j, k) - 1.0/6.0 * u_pri(i, j + 1, k);
						q2L = 11.0/6.0 * u_pri(i, j, k) - 7.0/6.0 * u_pri(i, j + 1, k) + 1.0/3.0 * u_pri(i, j + 2, k);
						
						
						q0R = 1.0/3.0 * u_pri(i, j - 2, k) - 7.0/6.0 * u_pri(i, j - 1, k) + 11.0/6.0 * u_pri(i, j, k);
						q1R = -1.0/6.0 * u_pri(i, j - 1, k) + 5.0/6.0 * u_pri(i, j, k) + 1.0/3.0 * u_pri(i, j + 1, k);
						q2R = 1.0/3.0 * u_pri(i, j, k) + 5.0/6.0 * u_pri(i, j + 1, k) - 1.0/6.0 * u_pri(i, j + 2, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u_pri(i, j, k) - u_pri(i, j - 1, k));
						diff2 = fabs(u_pri(i, j + 1, k) - u_pri(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u_pri(i, j + 2, k) - u_pri(i, j + 1, k)) - fabs(u_pri(i, j + 1, k) - u_pri(i, j, k)));
							diff4 = fabs(abs(u_pri(i, j + 1, k) - u_pri(i, j, k)) - fabs(u_pri(i, j, k) - u_pri(i, j - 1, k)));
							if (diff3 >= diff4){
								w_yL(i, j, k) = q1L;
								w_yR(i, j, k) = q1R;
							}else{
								w_yL(i, j, k) = q2L;
								w_yR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u_pri(i, j, k) - u_pri(i, j - 1, k)) - fabs(u_pri(i, j - 1, k) - u_pri(i, j - 2, k)));
							diff4 = fabs(abs(u_pri(i, j + 1, k) - u_pri(i, j, k)) - fabs(u_pri(i, j, k) - u_pri(i, j - 1, k)));
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
						q0 = -1.0/6.0 * u_pri(i - 2, j, k) + 5.0/6.0 * u_pri(i - 1, j, k) + 1.0/3.0 * u_pri(i, j, k);
						q1 = 1.0/3.0 * u_pri(i - 1, j, k) + 5.0/6.0 * u_pri(i, j, k) - 1.0/6.0 * u_pri(i + 1, j, k);
						q2 = 11.0/6.0 * u_pri(i, j, k) - 7.0/6.0 * u_pri(i + 1, j, k) + 1.0/3.0 * u_pri(i + 2, j, k);
						
						beta0 = 13.0/12.0 * pow((u_pri(i - 2, j, k) - 2 * u_pri(i - 1, j, k) + u_pri(i, j, k)),2.0) + 1.0/4.0 * pow((u_pri(i - 2, j, k) - 4 * u_pri(i - 1, j, k) + 3 * u_pri(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u_pri(i - 1, j, k) - 2 * u_pri(i, j, k) + u_pri(i + 1, j, k)),2.0) + 1.0/4.0 * pow((u_pri(i - 1, j, k) - u_pri(i + 1, j, k)),2.0);
						beta2 = 13.0/12.0 * pow((u_pri(i, j, k) - 2 * u_pri(i + 1, j, k) + u_pri(i + 2, j, k)),2.0) + 1.0/4.0 * pow((3 * u_pri(i, j, k) - 4 * u_pri(i + 1, j, k) + u_pri(i + 2, j, k)),2.0);
						
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
						
						q0 = 1.0/3.0 * u_pri(i - 2, j, k) - 7.0/6.0 * u_pri(i - 1, j, k) + 11.0/6.0 * u_pri(i, j, k);
						q1 = -1.0/6.0 * u_pri(i - 1, j, k) + 5.0/6.0 * u_pri(i, j, k) + 1.0/3.0 * u_pri(i + 1, j, k);
						q2 = 1.0/3.0 * u_pri(i, j, k) + 5.0/6.0 * u_pri(i + 1, j, k) - 1.0/6.0 * u_pri(i + 2, j, k);
						
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
						
						
						q0 = -1.0/6.0 * u_pri(i, j - 2, k) + 5.0/6.0 * u_pri(i, j - 1, k) + 1.0/3.0 * u_pri(i, j, k);
						q1 = 1.0/3.0 * u_pri(i, j - 1, k) + 5.0/6.0 * u_pri(i, j, k) - 1.0/6.0 * u_pri(i, j + 1, k);
						q2 = 11.0/6.0 * u_pri(i, j, k) - 7.0/6.0 * u_pri(i, j + 1, k) + 1.0/3.0 * u_pri(i, j + 2, k);
						
						beta0 = 13.0/12.0 * pow((u_pri(i, j - 2, k) - 2 * u_pri(i, j - 1, k) + u_pri(i, j, k)),2.0) + 1.0/4.0 * pow((u_pri(i, j - 2, k) - 4 * u_pri(i, j - 1, k) + 3 * u_pri(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u_pri(i, j - 1, k) - 2 * u_pri(i, j, k) + u_pri(i, j + 1, k)),2.0) + 1.0/4.0 * pow((u_pri(i, j - 1, k) - u_pri(i, j + 1, k)),2.0);
						beta2 = 13.0/12.0 * pow((u_pri(i, j, k) - 2 * u_pri(i, j + 1, k) + u_pri(i, j + 2, k)),2.0) + 1.0/4.0 * pow((3 * u_pri(i, j, k) - 4 * u_pri(i, j + 1, k) + u_pri(i, j + 2, k)),2.0);
						
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
						
						q0 = 1.0/3.0 * u_pri(i, j - 2, k) - 7.0/6.0 * u_pri(i, j - 1, k) + 11.0/6.0 * u_pri(i, j, k);
						q1 = -1.0/6.0 * u_pri(i, j - 1, k) + 5.0/6.0 * u_pri(i, j, k) + 1.0/3.0 * u_pri(i, j + 1, k);
						q2 = 1.0/3.0 * u_pri(i, j, k) + 5.0/6.0 * u_pri(i, j + 1, k) - 1.0/6.0 * u_pri(i, j + 2, k);
						
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
	
	upri2ucon(w_xL, w_xL);
	upri2ucon(w_xR, w_xR);
	upri2ucon(w_yL, w_yL);
	upri2ucon(w_yR, w_yR);
	//	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
	//		for(int k = 0; k < NUM_VARS; k++){
	//			wL(0, j, k) = wL(ghost_cells - 1, j, k);
	//			wR(0, j, k) = wR(ghost_cells - 1, j, k);
	//			wL(1, j, k) = wL(ghost_cells - 1, j, k);
	//			wR(1, j, k) = wR(ghost_cells - 1, j, k);
	//			
	//			wL(nxCells + ghost_cells + 1, j, k) = wL(nxCells + ghost_cells, j, k);
	//			wR(nxCells + ghost_cells + 1, j, k) = wR(nxCells + ghost_cells, j, k);
	//			wL(nxCells + ghost_cells + 2, j, k) = wL(nxCells + ghost_cells, j, k);
	//			wR(nxCells + ghost_cells + 2, j, k) = wR(nxCells + ghost_cells, j, k);
	//			
	//		}
	//	}
	bound_cond_u(w_xL);
	bound_cond_u(w_xR);
	bound_cond_u(w_yL);
	bound_cond_u(w_yR);
}


double calecon(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double U = con[3];
	double v_square = vx * vx + vy * vy;
	double e = (U - 0.5 * rho * v_square)/rho;
	return e;
}

double calepri(double *pri){
	double e = pri[3]/((gamma_plasma - 1) * pri[0]);
	return e;
}

vector<double> ConvergenceTestFunction(double x, double y, double t) {
	vector<double> f(NUM_VARS);
	const double g    = gamma_plasma;
	const double beta = 5.0;
	
	// Hu & Shu isentropic vortex (paper setup):
	//   background  (ρ∞, u∞, v∞, p∞) = (1, 1, 1, 1)
	//   domain      [0,10]² with periodic BCs
	const double u_inf = 1.0;
	const double v_inf = 1.0;
	
	// Vortex centre translates with the background flow,
	// wrapped periodically into [0,10]².
	const double xc = 5.0 + u_inf * t;
	const double yc = 5.0 + v_inf * t;
	double x_rel = x - xc;
	x_rel -= 10.0 * floor((x_rel + 5.0) / 10.0);
	double y_rel = y - yc;
	y_rel -= 10.0 * floor((y_rel + 5.0) / 10.0);
	const double r2 = x_rel*x_rel + y_rel*y_rel;
	
	// Perturbations (Eqs. 69–71)
	const double fac = (beta / (2.0 * M_PI)) * exp(0.5 * (1.0 - r2));
	const double du  = -y_rel * fac;
	const double dv  =  x_rel * fac;
	
	const double dT  = -(g - 1.0) * beta * beta
	/ (8.0 * g * M_PI * M_PI) * exp(1.0 - r2);
	
	// Isentropic relations (Eq. 72)
	f[0] = pow(1.0 + dT, 1.0 / (g - 1.0));   // ρ
	f[1] = u_inf + du;                       // u
	f[2] = v_inf + dv;                       // v
	f[3] = pow(1.0 + dT, g / (g - 1.0));     // p
	return f;
}

vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, 
	const double& yc, const double& dy,
	function<vector<double>(double, double,double)> func, int nGaussPoints, double t_final)
{
	vector<double> Weights(nGaussPoints);
	vector<double> Location(nGaussPoints);
	
	switch(nGaussPoints){
		case 1:
			Location = {0.0};
			Weights  = {1.0};
			break;
		case 2:
			Location = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
			Weights  = {1.0, 1.0};
			break;
		case 3:
			Location = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
			Weights  = {5.0/9.0, 8.0/9.0, 5.0/9.0};
			break;
		default:
			std::cerr << "Number of Gauss-Legendre nodes not implemented\n";
			std::exit(EXIT_FAILURE);
	}
	
	int nVars = func(xc, yc,t_final).size();
	vector<double> integral(nVars, 0.0);
	
	for(int i = 0; i < nGaussPoints; i++){
		double xG = xc + 0.5 * dx * Location[i];
		double wx = Weights[i];
		for(int j = 0; j < nGaussPoints; j++){
			double yG = yc + 0.5 * dy * Location[j];
			double wy = Weights[j];
			vector<double> f = func(xG, yG,t_final);
			for(int k = 0; k < nVars; k++){
				integral[k] += 0.25 * wx * wy * f[k];
			}
		}
	}
	
	return integral;
}

void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p){
	vector<double> result;
	double Mx;
	double alpha = 40;
	double M0 = 1.0/ 11.0;
	double L = 2.0 / M0;
	double rho0 = 0.955, rho1 = 2;
	double u0 = 2 * sqrt(gamma_plasma);
	double p0 = 1, p1 = 2 * gamma_plasma;
	double r;
	double velocity;
	double uL[NUM_VARS], uR[NUM_VARS], uIN[NUM_VARS], uOUT[NUM_VARS], Uavg[NUM_VARS];;
	double xG[2];double yG[2];double xg;double yg;
	
	
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
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			double r2 = (x - 2) * (x - 2) + (y - 2) * (y - 2);
			int sgn = x >= 0 ? 1 : -1;
			switch (test) {
				case testToroModified:
					if(x < 0.3){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.75;
						u(i,j,2) = 0;
						u(i,j,3) = 1;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.1;
					}
					break;

				case testToro1:
					if(x < 0.5){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.1;
					}
					break;
				case testToro2:
					if(x < 0.5){
						u(i,j,0) = 1;
						u(i,j,1) = -2;
						u(i,j,2) = 0;
						u(i,j,3) = 0.4;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 0;
						u(i,j,3) = 0.4;
					}
					break;
				case testToro3:
					if(x < 0.5){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1000;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.01;
					}
					break;
				case testToro4:
					if(x < 0.5){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.01;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 100;
					}
					break;
				case testToro5:
					if(x < 0.5){
						u(i,j,0) = 5.99924;
						u(i,j,1) = 19.5975;
						u(i,j,2) = 0;
						u(i,j,3) = 460.894;
					}else{
						u(i,j,0) = 5.99242;
						u(i,j,1) = -6.19633;
						u(i,j,2) = 0;
						u(i,j,3) = 46.0950;
					}
					break;
				case RPContact:
					if(x < 0.25){
						u(i,j,0) = 1000;
						u(i,j,1) = 1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 1e5;
					}else{
						u(i,j,0) = 0.01;
						u(i,j,1) = 1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 1e5;
					}
//					u(i,j,0) = 0.5 * (1000 + 0.01) + 0.5 * (1000 - 0.01) * erf((x - 0.25)/(dx));
					break;
				case RPLax:
					if(x < 0.5){
						u(i,j,0) = 0.445;
						u(i,j,1) = 0.698;
						u(i,j,2) = 0;
						u(i,j,3) = 3.528;
					}else{
						u(i,j,0) = 0.5;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.571;
					}
					break;
				case RPTwoShocks:
					if(x < 0.5){
						u(i,j,0) = 1.0;
						u(i,j,1) = 2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.1;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.1;
					}
					break;
				case RPTwoRare:
					if(x < 0.5){
						u(i,j,0) = 1.0;
						u(i,j,1) = -1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.4;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = 1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.4;
					}
					break;
				case RPFVVEM1:
					if(x < 0.5){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0;
						u(i,j,3) = 1.0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.1;
					}
					break;
				case RPFVVEM2:
					if(x < 0.5){
						u(i,j,0) = 1.0;
						u(i,j,1) = -1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.4;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = 1.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.4;
					}
					break;
				case RPFVVEM3:
					if(x < 0.5){
						u(i,j,0) = 0.445;
						u(i,j,1) = 0.698;
						u(i,j,2) = 0;
						u(i,j,3) = 3.528;
					}else{
						u(i,j,0) = 0.5;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.571;
					}
					break;
				case RPFVVEM4:{
					if(x < 0.8){
						u(i,j,0) = 1;
						u(i,j,1) = -19.59745;
						u(i,j,2) = 0;
						u(i,j,3) = 1000.0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -19.59745;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.01;
					}
					u(i,j,3) = 0.5 * (0.01 + 1000) + 0.5 * (0.01 - 1000) * erf((x - 0.8)/(dx));
					break;
				}
				case testExplosion:
					if(sqrt((x - 1.0) * (x - 1.0) + (y - 1.0) * (y - 1.0)) <= 0.4){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.1;
					}
					break;
				case testSmoothAcousticWave:
					u(i,j,0) = 1;
					u(i,j,1) = 0;
					u(i,j,2) = 0;
					u(i,j,3) = 1 + exp(-40 * r2);
					break;
				case testDoubleShearLayer:
					u(i,j,0) = 1;
					if(y < 0.5){
						u(i,j,1) = tanh(30 * (y - 0.25));
					}else{
						u(i,j,1) = tanh(30 * (0.75 - y));
					}
					u(i,j,2) = 0.05 * sin(2 * M_PI * x);
					u(i,j,3) = 1e5 / gamma_plasma;		
					break;
				case testWoodwardColellaBlastWave:
					if(x < 0.1){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1000;
					}else if(x < 0.9){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0.01;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 100;
					}
					break;
				case testGreshoVortex:
					u(i,j,0) = 1;
					r = sqrt(x * x + y * y);
					p0 = u(i,j,0) / (gamma_plasma * vortex_Mach * vortex_Mach);
					if(r < 0.2){
						velocity = 5 * r;
						u(i,j,3) = p0 + 25.0 / 2.0 * r * r;
					}else if(r < 0.4){
						velocity = 2 - 5 * r;
						u(i,j,3) = p0 + 25.0 / 2.0 * r * r + 4 * (1 - 5 * r - log(0.2) + log(r));
					}else{
						velocity = 0;
						u(i,j,3) = p0 - 2 + 4 * log(2);
					}
					
					if(r < 1e-14){
						u(i,j,1) = 0;
						u(i,j,2) = 0;
					} else {
						u(i,j,1) = -velocity * y / r;
						u(i,j,2) = velocity * x / r;
					}
					break;
				case testCollidingPulses:
					p_ref = (rho0 / gamma_plasma) * pow(u0 / M0, 2);
					u(i,j,0) = rho0 + M0 / 2.0 * rho1 * (1 - cos(2.0 * M_PI * x / L));
					u(i,j,1) = -u0 / 2.0 * sgn * (1 - cos(2.0 * M_PI * x / L));
					u(i,j,2) = 0;
					u(i,j,3) = (p0 + M0 / 2.0 * p1 * (1 - cos(2.0 * M_PI * x / L))) * p_ref;
					break;
				case testRiemannLike:
//					if(x < 0){
//						u(i,j,0) = 0.15;
//						u(i,j,1) = 21.55;
//						u(i,j,2) = 1.0;
//						u(i,j,3) = 0.28;
//					}else{
//						u(i,j,0) = 0.1;
//						u(i,j,1) = -26.45;
//						u(i,j,2) = 0;
//						u(i,j,3) =  0.1;
//					}
					
					p_epsilon = 0.99;
					
					
					if(x < 0.2 || x > 0.8){
						u(i,j,0) = 1;
						u(i,j,1) = 1 - p_epsilon / 2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 1 * 1 / p_epsilon;
					}else if(x > 0.25 && x < 0.75){
						u(i,j,0) = 1;
						u(i,j,1) = 1 + p_epsilon / 2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 1 * 1 / p_epsilon;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 1;
						u(i,j,2) = 0;
						u(i,j,3) = 1 * 1 / p_epsilon;
					}
					
					
//					uIN[0] = 0.15; 
//					uIN[1] = 21.55; 
//					uIN[2] = 1.0; 
//					uIN[3] = 0.28; 
//					
//					
//					uOUT[0] = 0.1; 
//					uOUT[1] = -26.45; 
//					uOUT[2] = 0; 
//					uOUT[3] = 0.1; 
//					
//					
//					xG[0] = x - 0.5 * dx * a;
//					xG[1] = x + 0.5 * dx * a;
//					yG[0] = y - 0.5 * dy * a;
//					yG[1] = y + 0.5 * dy * a;
//					
//					w = 0.25;
//					beta_s = 1 / (6 * dx);
//					for (int k = 0; k < NUM_VARS; ++k) Uavg[k] = 0.0;
//					
//					for (int jj = 0; jj < 2; ++jj) {
//						for (int ii = 0; ii < 2; ++ii) {
//							
//							double xg = xG[ii];
//							double yg = yG[jj]; // 这里不参与 s，但保留结构
//							
//							// 平滑阶跃：x<0 -> s≈1；x>0 -> s≈0
//							double s_x = 0.5 * (1.0 - std::tanh(beta_s * (xg - 0)));
//							
//							for (int k = 0; k < NUM_VARS; ++k) {
//								double U = uOUT[k] + (uIN[k] - uOUT[k]) * s_x;
//								Uavg[k] += w * U;
//							}
//						}
//					}
//					
//					for (int k = 0; k < NUM_VARS; ++k) {
//						u(i, j, k) = Uavg[k];
//					}

					break;
				case convergencetest:{
					// Background
					result = GaussLegendreNumericalIntegration(x, dx, y,dy,ConvergenceTestFunction, 3, 0);
					u(i,j,0) = result[0];
					u(i,j,1) = result[1];
					u(i,j,2) = result[2];
					u(i,j,3) = result[3];
					
					break;
				}
				case test2d1:
					if(x > 0.3 && y > 0.3){
						u(i,j,0) = 1.5;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1.5;
					}else if(x <= 0.3 && y > 0.3){
						u(i,j,0) = 0.5323;
						u(i,j,1) = 1.206;
						u(i,j,2) = 0;
						u(i,j,3) = 0.3;
					}else if(x <= 0.3 && y <= 0.3){
						u(i,j,0) = 0.138;
						u(i,j,1) = 1.206;
						u(i,j,2) = 1.206;
						u(i,j,3) = 0.029;
					}else{
						u(i,j,0) = 0.5323;
						u(i,j,1) = 0;
						u(i,j,2) = 1.206;
						u(i,j,3) = 0.3;
					}
					break;
				case test2d2:
					if(x > 0 && y > 0){
						u(i,j,0) = 1.1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 1.1;
					}else if(x <= 0 && y > 0){
						u(i,j,0) = 0.5065;
						u(i,j,1) = 0.8939;
						u(i,j,2) = 0;
						u(i,j,3) = 0.35;
					}else if(x <= 0 && y <= 0){
						u(i,j,0) = 1.1;
						u(i,j,1) = 0.8939;
						u(i,j,2) = 0.8939;
						u(i,j,3) = 1.1;
					}else{
						u(i,j,0) = 0.5065;
						u(i,j,1) = 0;
						u(i,j,2) = 0.8939;
						u(i,j,3) = 0.35;
					}
					break;
				case test2d3:
					if(x > 0 && y > 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = -0.75;
						u(i,j,2) = -0.5;
						u(i,j,3) = 1.0;
					}else if(x <= 0 && y > 0){
						u(i,j,0) = 2.0;
						u(i,j,1) = -0.75;
						u(i,j,2) = 0.5;
						u(i,j,3) = 1.0;
					}else if(x <= 0 && y <= 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.75;
						u(i,j,2) = 0.5;
						u(i,j,3) = 1.0;
					}else{
						u(i,j,0) = 3.0;
						u(i,j,1) = 0.75;
						u(i,j,2) = -0.5;
						u(i,j,3) = 1.0;
					}
					break;
				case test2d4:
					if(x > 0 && y > 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.75;
						u(i,j,2) = -0.5;
						u(i,j,3) = 1.0;
					}else if(x <= 0 && y > 0){
						u(i,j,0) = 2.0;
						u(i,j,1) = 0.75;
						u(i,j,2) = 0.5;
						u(i,j,3) = 1.0;
					}else if(x <= 0 && y <= 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = -0.75;
						u(i,j,2) = 0.5;
						u(i,j,3) = 1.0;
					}else{
						u(i,j,0) = 3.0;
						u(i,j,1) = -0.75;
						u(i,j,2) = -0.5;
						u(i,j,3) = 1.0;
					}
					break;
				case KH_instability:
					u(i, j, 0) = 1.0 - 0.25 * std::tanh(25.0 * (std::abs(y) - 0.5));
					u(i, j, 1) = -0.5 * std::tanh(25.0 * (std::abs(y) - 0.5));
					u(i, j, 2) = (1.0 / 100.0) * std::sin(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
					u(i, j, 3) = 1.0e4 / gamma_plasma;
					break;
			}
		}
	}
}


void bound_cond_u(Grid &u){
	double p0 = 1 / (gamma_plasma * vortex_Mach * vortex_Mach);
	double p = p0 - 2 + 4 * log(2);
	double E = p / (gamma_plasma - 1);
	
	
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
					
					
					u(ghost,j,3) = u(2 * ghost_cells - 1 - ghost,j,3);
					u(nxCells + ghost_cells + ghost,j,3) = u(nxCells + ghost_cells - 1 - ghost,j,3);
					
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
		case x_Dirichlet:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					if(test == testGreshoVortex){
						u(ghost,j,0) = 1;
						u(ghost,j,1) = 0;
						u(ghost,j,2) = 0;
						u(ghost,j,3) = E;
						
						u(nxCells + ghost_cells + ghost,j,0) = 1;
						u(nxCells + ghost_cells + ghost,j,1) = 0;
						u(nxCells + ghost_cells + ghost,j,2) = 0;
						u(nxCells + ghost_cells + ghost,j,3) = E;
					}else if(test == RPContact){
						double priL[NUM_VARS], priR[NUM_VARS],conL[NUM_VARS], conR[NUM_VARS];
						priL[0] = 1000;
						priL[1] = 0;
						priL[2] = 0;
						priL[3] = 10000;
						
						priR[0] = 0.01;
						priR[1] = 0;
						priR[2] = 0;
						priR[3] = 10000;
						
						ini_pri2con(conL, priL);
						ini_pri2con(conR, priR);
						
						u(ghost,j,0) = conL[0];
						u(ghost,j,1) = conL[1];
						u(ghost,j,2) = conL[2];
						u(ghost,j,3) = conL[3];
						
						u(nxCells + ghost_cells + ghost,j,0) = conR[0];
						u(nxCells + ghost_cells + ghost,j,1) = conR[1];
						u(nxCells + ghost_cells + ghost,j,2) = conR[2];
						u(nxCells + ghost_cells + ghost,j,3) = conR[3];
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
					
					u(i,ghost,1) = -u(i,2 * ghost_cells - 1 - ghost,1);
					u(i,nyCells + ghost_cells + ghost,1) = -u(i,nyCells + ghost_cells - 1 - ghost,1);
					
					u(i,ghost,2) = -u(i,2 * ghost_cells - 1 - ghost,2);
					u(i,nyCells + ghost_cells + ghost,2) = -u(i,nyCells + ghost_cells - 1 - ghost,2);
					
					u(i,ghost,3) = u(i,2 * ghost_cells - 1 - ghost,3);
					u(i,nyCells + ghost_cells + ghost,3) = u(i,nyCells + ghost_cells - 1 - ghost,3);
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
		case y_Dirichlet:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					if(test == testGreshoVortex){
						u(i,ghost,0) = 1;
						u(i,ghost,1) = 0;
						u(i,ghost,2) = 0;
						u(i,ghost,3) = E;
						
						u(i,nyCells + ghost_cells + ghost,0) = 1;
						u(i,nyCells + ghost_cells + ghost,1) = 0;
						u(i,nyCells + ghost_cells + ghost,2) = 0;
						u(i,nyCells + ghost_cells + ghost,3) = E;
					}else if(test == RPContact){
						double priL[NUM_VARS], priR[NUM_VARS],conL[NUM_VARS], conR[NUM_VARS];
						priL[0] = 1000;
						priL[1] = 0;
						priL[2] = 0;
						priL[3] = 10000;
						
						priR[0] = 0.01;
						priR[1] = 0;
						priR[2] = 0;
						priR[3] = 10000;
						
						ini_pri2con(conL, priL);
						ini_pri2con(conR, priR);
						
						u(i,ghost,0) = conL[0];
						u(i,ghost,1) = conL[1];
						u(i,ghost,2) = conL[2];
						u(i,ghost,3) = conL[3];
						
						u(i,nyCells + ghost_cells + ghost,0) = conR[0];
						u(i,nyCells + ghost_cells + ghost,1) = conR[1];
						u(i,nyCells + ghost_cells + ghost,2) = conR[2];
						u(i,nyCells + ghost_cells + ghost,3) = conR[3];
					}
					
				}
			}
			break;
	}
}

double caldt(Grid &u_pri){
	vx_max = 0;
	double dt;
	double lambda_11_x, lambda_22_x, lambda_33_x,lambda_44_x;
	double lambda_11_y, lambda_22_y, lambda_33_y,lambda_44_y;
	double lambda_1_x = 0, lambda_2_x = 0, lambda_3_x = 0, lambda_4_x = 0;
	double lambda_1_y = 0, lambda_2_y = 0, lambda_3_y = 0, lambda_4_y = 0;
	double vx_cs = 0, vy_cs = 0;
	lambda = 0, lambda_x = 0, lambda_y = 0;
	vx_max = 0;
	vy_max = 0;
	
	
	switch (caltimestep) {
		case Fullydt:
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
					double a = sqrt(gamma_plasma * u_pri(i,j,3) / u_pri(i,j,0));
					a_max_x = fmax(fabs(u_pri(i,j,1)) + a, fabs(lambda_x));
					a_max_y = fmax(fabs(u_pri(i,j,2)) + a, fabs(lambda_y));
					vx_max = fmax(fabs(vx_max), fabs(u_pri(i,j,1)));
					vy_max = fmax(fabs(vy_max), fabs(u_pri(i,j,2)));
				}
			}
			break;
		case IMEXdt:
			if(timestep < full_step + 1){
				for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
					for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
						double a = sqrt(gamma_plasma * u_pri(i,j,3) / u_pri(i,j,0));
						a_max_x = fmax(fabs(u_pri(i,j,1)) + a, fabs(lambda_x));
						a_max_y = fmax(fabs(u_pri(i,j,2)) + a, fabs(lambda_y));
						vx_max = fmax(fabs(vx_max), fabs(u_pri(i,j,1)));
						vy_max = fmax(fabs(vy_max), fabs(u_pri(i,j,2)));
					}
				}
			}else{
				for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
					for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
						vx_max = fmax(fabs(vx_max), fabs(u_pri(i,j,1)));
						vy_max = fmax(fabs(vy_max), fabs(u_pri(i,j,2)));
						
						a_max_x = vx_max;
						a_max_y = vy_max;
					}
				}
			}
			break;
	}

	dt = CFL * 1 / (a_max_x / dx + lambda_y / dy);
	return dt;
}

void bound_cond_array(Array &u, x_boundary_condition xbc = x_bc, y_boundary_condition ybc = y_bc){
	switch (xbc) {
		case x_Transmissive: case x_Dirichlet:
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
	
	switch(ybc){
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
	double h = e + pri[3] / pri[0];
	return h;
}

double calkinetic(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy);
	return kinetic_energy;
}


// 从原始变量计算声速 c = sqrt(γ p / ρ)
// pri[0] = ρ, pri[1] = u, pri[2] = v, pri[3] = p
void ShockDetector(Array &shock_eta, Grid &u_con) {
	// 转换守恒变量到原始变量(ρ, u, v, p)
	Grid u_pri(nxCells, nyCells);
	ucon2upri(u_pri, u_con);
	
	
	// 中间数组,先存 χ̃(公式 4.40),然后再做邻居扩展(公式 4.41)
	Array shock_eta_tilde(nxCells, nyCells);
	
	// ---- 第一步:计算 χ̃ (公式 4.40) ----
	for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			shock_eta_tilde(i, j) = 0.0;
			// 处理三种开关情况
			if (add_shock_detector == No_add) {
				shock_eta_tilde(i, j) = 1.0;
				continue;
			}
			
			
			// --- Yes_add 情况:正常计算 ---
			
			// 计算 3x3 邻域内的声速 c = sqrt(γ p / ρ),取最小值
			// u_pri(.,.,0) = ρ,  u_pri(.,.,3) = p
			double cf_min_nbr = 1.0e300;
			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					double rho_nbr = u_pri(i + di, j + dj, 0);
					double p_nbr   = u_pri(i + di, j + dj, 3);
					double cf      = std::sqrt(gamma_plasma * p_nbr / rho_nbr);
					if (cf < cf_min_nbr) cf_min_nbr = cf;
				}
			}
			
			// 计算速度散度 ∇·u(中心差分)
			// u_pri(.,.,1) = u, u_pri(.,.,2) = v
			double div_u = (u_pri(i + 1, j, 1) - u_pri(i - 1, j, 1)) / (2.0 * dx)
			+ (u_pri(i, j + 1, 2) - u_pri(i, j - 1, 2)) / (2.0 * dy);
			
			// 公式 (4.40):
			//   χ̃ = min[1, max(0, -(h·∇·u + k₁ a_min) / (k₁ a_min))]
			// 把网格尺度 h 吃进分子,让分子也是速度量纲,与 k₁·a_min 一致
			double h = std::min(dx, dy);
			double numerator   = h * div_u + shock_kappa1 * cf_min_nbr;
			double denominator = shock_kappa1 * cf_min_nbr;
			
			shock_eta_tilde(i, j) = std::min(1.0,
				std::max(0.0, -numerator / denominator));
		}
	}
	
	// 边界条件(让邻域扩展可以访问 ghost cell)
	bound_cond_array(shock_eta_tilde);
	
	// ---- 第二步:邻居扩展 (公式 4.41) ----
	//   χ_{i,j} = max over 3x3 neighborhood of χ̃
	// 这样激波"即将经过"的单元也提前被标记
	if (add_shock_detector == Yes_add) {
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
			}
		}
	} else {
		// No_add: 直接复制 χ̃ 的值
		for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
			for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
				shock_eta(i, j) = shock_eta_tilde(i, j);
			}
		}
	}
	bound_cond_array(shock_eta);
}

int M_P(int i, int j){
	int coord = 0;
	int xL, xR, yL, yR;
	int I, J;
	
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive:case x_Dirichlet: case x_Reflective:
				I = ghost_cells;
				break;
			case x_Periodic:
				I = nxCells + i;
				break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:case x_Reflective:
				I = nxCells + ghost_cells - 1;
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
			case y_Transmissive: case y_Dirichlet:case y_Reflective:
				J = ghost_cells;
				break;
			case y_Periodic:
				J = nyCells + j;
				break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive:case y_Dirichlet:case y_Reflective:
				J = nyCells + ghost_cells - 1;
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

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h){
	
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

			double p_bc = 1.0 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2.0);
			
			// 判断 (ii,jj) 是不是 Dirichlet ghost cell
			auto is_dir_ghost = [&](int ii, int jj) -> bool {
				return (x_bc == x_Dirichlet && (ii < ghost_cells || ii >= nxCells + ghost_cells))
				|| (y_bc == y_Dirichlet && (jj < ghost_cells || jj >= nyCells + ghost_cells));
			};
			
			// Ax 里只取未知数部分，ghost 返回 0
			auto get_x_unknown = [&](int ii, int jj) -> double {
				int II = ii, JJ = jj;
				double sign = 1.0;
				
				// x 方向 Dirichlet ghost → 镜像 + 取负
				if (x_bc == x_Dirichlet) {
					if (ii < ghost_cells) {
						II = 2*ghost_cells - 1 - ii;          // -1→0, -2→1 (左侧镜像)
						sign = -sign;
					} else if (ii >= nxCells + ghost_cells) {
						II = 2*(nxCells + ghost_cells) - 1 - ii;  // 右侧镜像
						sign = -sign;
					}
				}
				// y 方向 Dirichlet ghost → 镜像 + 取负
				if (y_bc == y_Dirichlet) {
					if (jj < ghost_cells) {
						JJ = 2*ghost_cells - 1 - jj;
						sign = -sign;
					} else if (jj >= nyCells + ghost_cells) {
						JJ = 2*(nyCells + ghost_cells) - 1 - jj;
						sign = -sign;
					}
				}
				return sign * x(M_P(II, JJ));
			};
			
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
			
			double p_i = x(M_P(i,j));
						
			
		
			// === composed 5-point operator ===
			
			// x方向 Step 1: pressure gradient at cells i-1, i, i+1
			double gradP_im = (1.0/dx) * (
				0.5*(p_i + get_x_unknown(i-1,j))- 0.5*(get_x_unknown(i-1,j) + get_x_unknown(i-2,j))
			);
			double gradP_i = (1.0/dx) * (
				0.5*(get_x_unknown(i+1,j) + p_i) - 0.5*(p_i + get_x_unknown(i-1,j))
			);
			double gradP_ip = (1.0/dx) * (
				0.5*(get_x_unknown(i+2,j) + get_x_unknown(i+1,j))- 0.5*(get_x_unknown(i+1,j) + p_i)
			);
			
			// x方向 Step 2: energy divergence
			double Lap_x = (1.0/dx) * (
				0.5*(hp0+h00) * 0.5*(gradP_i + gradP_ip)- 0.5*(hm0+h00) * 0.5*(gradP_im + gradP_i)
			);
			
			// y方向 Step 1
			double gradP_jm = (1.0/dy) * (
				0.5*(p_i + get_x_unknown(i,j-1))- 0.5*(get_x_unknown(i,j-1) + get_x_unknown(i,j-2))
			);
			double gradP_j = (1.0/dy) * (
				0.5*(get_x_unknown(i,j+1) + p_i)- 0.5*(p_i + get_x_unknown(i,j-1))
			);
			double gradP_jp = (1.0/dy) * (
				0.5*(get_x_unknown(i,j+2) + get_x_unknown(i,j+1))-  0.5*(get_x_unknown(i,j+1) + p_i)
			);
			
			// y方向 Step 2
			double Lap_y = (1.0/dy) * (
				0.5*(h0p+h00) * 0.5*(gradP_j + gradP_jp)- 0.5*(h0m+h00) * 0.5*(gradP_jm + gradP_j)
			);
			
			Ax(M_P(i,j)) = 1.0/(gamma_plasma - 1) * p_i - h*h * (Lap_x + Lap_y);
				
			double nu_x;
			double nu_y;
			
			switch (add_pvisc) {
				case Noaddpvis:
					break;
				case MaterialSpeedp:{
					switch (flux_diss) {
						case Global_Diss:
							nu_x = vx_max * dx;
							nu_y = vy_max * dy;
							break;
						case Local_Diss:
							nu_x = fabs(Qe(i,j,1) / Qe(i,j,0)) * dx;
							nu_y = fabs(Qe(i,j,2) / Qe(i,j,0)) * dy;
							break;
					}
					
					
					double lap_p = nu_x/(dx*dx) * ((get_x_unknown(i+1,j) - x(M_P(i,j)))
						- (x(M_P(i,j)) - get_x_unknown(i-1,j)))
					+ nu_y/(dy*dy) * ((get_x_unknown(i,j+1) - x(M_P(i,j)))
						- (x(M_P(i,j)) - get_x_unknown(i,j-1)));
					
					Ax(M_P(i,j)) -=  shock_eta(i,j) * 1.0/(gamma_plasma - 1) * h * lap_p;
					break;
				}
			}
		}
	}
}


void GMRES_Solver(Eigen::VectorXd& x, Eigen::VectorXd& b,
	std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, double)> Ax_func,
	Grid& u_con, Grid& Qe, Grid& Qs,double h, const int dim) {
		
		Eigen::VectorXd Ax(dim), r(dim), r_hat(dim), p(dim), s(dim), Ap(dim), As(dim);
		Ax.setZero();
		
		Ax_func(Ax, x, u_con, Qe, Qs, h);
		r = b - Ax;
		r_hat = r;
		p = r;
		
		
		
		double rho_old = r_hat.dot(r);
		double r0 = sqrt(r.dot(r));
		double Ax_norm = sqrt(Ax.dot(Ax));
		double b_norm = sqrt(b.dot(b));
		
		
		// 计算 uniform background
		double p_bg = 1.0 / (gamma_plasma * vortex_Mach * vortex_Mach) / (gamma_plasma - 1.0);
		double b_pert = 0;
		for(int i = 0; i < b.size(); i++) {
			double diff = b(i) - p_bg;
			b_pert += diff * diff;
		}
		//		std::cout << "[DIAG2] ||b-b_bg||=" << sqrt(b_pert) << std::endl;
		
		if (r0 < GMREStol) {
			//			std::cout << "[BiCGSTAB] skip (r0 < tol)" << std::endl;
			return;
		}
		
		int final_iter = 0;
		double final_r = r0;
		bool converged = false;
		bool breakdown = false;
		
		for (int iter = 0; iter < GMRESiter; iter++) {
			Ap.setZero();
			Ax_func(Ap, p, u_con, Qe, Qs, h);
			
			double rhat_dot_Ap = r_hat.dot(Ap);
			if (std::abs(rhat_dot_Ap) < 1e-300) {
				breakdown = true;
				break;
			}
			
			double alpha = rho_old / rhat_dot_Ap;
			s = r - alpha * Ap;
			
			if (sqrt(s.dot(s)) < GMREStol) {
				x += alpha * p;
				final_iter = iter;
				final_r = sqrt(s.dot(s));
				converged = true;
				break;
			}
			
			As.setZero();
			Ax_func(As, s, u_con, Qe, Qs, h);
			
			double As_dot_As = As.dot(As);
			if (As_dot_As < 1e-300) {
				breakdown = true;
				break;
			}
			
			double omega = As.dot(s) / As_dot_As;
			
			x += alpha * p + omega * s;
			r = s - omega * As;
			
			double rsnew = sqrt(r.dot(r));
			final_iter = iter;
			final_r = rsnew;
			
			if (rsnew < GMREStol) {
				converged = true;
				break;
			}
			
			double rho_new = r_hat.dot(r);
			double beta = (rho_new / rho_old) * (alpha / omega);
			p = r + beta * (p - omega * Ap);
			
			rho_old = rho_new;
		}
	}


void PressureSubsystem(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR,double h, const int rks){
//	Array shock_eta(nxCells,nyCells);
	double Bx_star, By_star, Bz_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_p_vect(nxCells*nyCells),p_vect(nxCells * nyCells);
	
	
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	
//	ShockDetector(shock_eta, Qs);
	Array rhovx_star(nxCells, nyCells);
	Array rhovy_star(nxCells, nyCells);
	Array p_solution(nxCells, nyCells);
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			rhovx_star(i,j) = u_con(i,j,1);
			rhovy_star(i,j) = u_con(i,j,2);
		}
	}
	bound_cond_array(rhovx_star);
	bound_cond_array(rhovy_star);
	
	if(x_bc == x_Dirichlet){
		for(int j = 0; j < nyCells + ghost_cells * 2; j++){
			for(int k = 0; k < ghost_cells; k++){
				rhovx_star(k,j) = 0;
				rhovx_star(nxCells + ghost_cells + k,j) = 0;
			}
		}
		
		for(int j = 0; j < nyCells + ghost_cells * 2; j++){
			for(int k = 0; k < ghost_cells; k++){
				rhovy_star(k,j) = 0;
				rhovy_star(nxCells + ghost_cells + k,j) = 0;
			}
		}
	}
	
	if(y_bc == y_Dirichlet){
		for(int i = 0; i < nxCells + ghost_cells * 2; i++){
			for(int ghost = 0; ghost < ghost_cells; ghost++){
				rhovx_star(i,ghost) = 0;
				rhovx_star(i,nyCells + ghost_cells + ghost) = 0;
			}
		}
		
		for(int i = 0; i < nxCells + ghost_cells * 2; i++){
			for(int ghost = 0; ghost < ghost_cells; ghost++){
				rhovy_star(i,ghost) = 0;
				rhovy_star(i,nyCells + ghost_cells + ghost) = 0;
			}
		}
	}	
	
//	std::cout << "First RHS: "<<endl;
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
			
//			double h00 = calenthalpy(con00_n);
//			double hp0 = calenthalpy(conp0_n);
//			double hm0 = calenthalpy(conm0_n);
//			double h0p = calenthalpy(con0p_n);
//			double h0m = calenthalpy(con0m_n);
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);

			
			double RHS =  u_con(i,j,3) - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j)) / u_con(i,j,0) - h * (1 / dx * ( 0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) -  0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * ( 0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));
			
			b_p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
			
			double bc_correction = 0.0;
			
			double p_bc = 1.0 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2.0);
			
			
			double nu_x = vx_max * dx;
			double nu_y = vy_max * dy;
			
			// 左邻居 Dirichlet
			auto get_x_const = [&](int ii, int jj) -> double {
				bool x_ghost = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
				bool y_ghost = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
				if (x_bc == x_Dirichlet && x_ghost) return p_bc;
				if (y_bc == y_Dirichlet && y_ghost) return p_bc;
				return 0.0;
			};
			
			double p_c = 0.0;   // interior 取 0
			
			// === x 方向 5 点 (copy 自 Ax，p_i→p_c, get_x_unknown→get_x_const) ===
			double gradP_im_c = (1.0/dx) * (
				0.5*(p_c + get_x_const(i-1,j))- 0.5*(get_x_const(i-1,j) + get_x_const(i-2,j))
			);
			double gradP_i_c = (1.0/dx) * (
				0.5*(get_x_const(i+1,j) + p_c)- 0.5*(p_c + get_x_const(i-1,j))
			);
			double gradP_ip_c = (1.0/dx) * (
				0.5*(get_x_const(i+2,j) + get_x_const(i+1,j))- 0.5*(get_x_const(i+1,j) + p_c)
			);
			double Lap_x_c = (1.0/dx) * (
				0.5*(hp0+h00) * 0.5*(gradP_i_c + gradP_ip_c)- 0.5*(hm0+h00) * 0.5*(gradP_im_c + gradP_i_c)
			);
			
			// === y 方向 5 点 ===
			double gradP_jm_c = (1.0/dy) * (
				0.5*(p_c + get_x_const(i,j-1))- 0.5*(get_x_const(i,j-1) + get_x_const(i,j-2))
			);
			double gradP_j_c = (1.0/dy) * (
				0.5*(get_x_const(i,j+1) + p_c)- 0.5*(p_c + get_x_const(i,j-1))
			);
			double gradP_jp_c = (1.0/dy) * (
				0.5*(get_x_const(i,j+2) + get_x_const(i,j+1))- 0.5*(get_x_const(i,j+1) + p_c)
			);
			double Lap_y_c = (1.0/dy) * (
				0.5*(h0p+h00) * 0.5*(gradP_j_c + gradP_jp_c)- 0.5*(h0m+h00) * 0.5*(gradP_jm_c + gradP_j_c)
			);
			
			double L_const = 1.0/(gamma_plasma - 1) * p_c - h*h * (Lap_x_c + Lap_y_c);
			
			// === 搬到 RHS：A_impl·p = b_true - L_const ===
			b_p_vect(M_P(i,j)) -= L_const;
			
			p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
		}
	}

	GMRES_Solver(p_vect, b_p_vect, Pressure_Subsystem_AX, u_con, Qe, Qs, h, nxCells*nyCells);
	
		
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			p_solution(i,j) = p_vect(M_P(i,j));
//			std::cout << p_solution(i,j) << endl;
		}
	}
	
	if(x_bc == x_Reflective){
		bound_cond_array(p_solution, x_Transmissive, y_Transmissive);
	}else{
		bound_cond_array(p_solution);
	}
	
	if(x_bc == x_Dirichlet){
		for(int j = 0; j < nyCells + ghost_cells * 2; j++){
			for(int k = 0; k < ghost_cells; k++){
				if(test == testGreshoVortex){
					p_solution(k,j) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
					p_solution(nxCells + ghost_cells + k,j) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
				}else if(test == RPContact){
					p_solution(k,j) = 1e5;
					p_solution(nxCells + ghost_cells + k,j) = 1e5;
				}
			}
		}
	}
	
	if(y_bc == y_Dirichlet){
		for(int i = 0; i < nxCells + ghost_cells * 2; i++){
			for(int ghost = 0; ghost < ghost_cells; ghost++){
				if(test == testGreshoVortex){
					p_solution(i,ghost) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
					p_solution(i,nyCells + ghost_cells + ghost) =  1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
				}else if(test == RPContact){
					p_solution(i,ghost) = 1e5;
					p_solution(i,nyCells + ghost_cells + ghost) = 1e5;
				}
			}
		}
	}	
	
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			double fL, fR,	flux_xm, flux_xp;
			fL =  p_solution(i - 1,j);
			fR =  p_solution(i,j);
			
			flux_xm = 0.5 * (fL + fR);
			
						
			fL =  p_solution(i,j);
			fR =  p_solution(i + 1,j);
			
			flux_xp = 0.5 * (fL + fR);
					
			u_con(i,j,1) = rhovx_star(i,j) - h / dx * (flux_xp - flux_xm);
			
			double gL, gR, flux_ym, flux_yp;
			gL = p_solution(i,j - 1);
			gR = p_solution(i,j);
			
			flux_ym = 0.5 * (gL + gR);
			
			gL = p_solution(i,j);
			gR = p_solution(i,j + 1);
			
			flux_yp = 0.5 * (gL + gR);
					
			
			u_con(i,j,2) = rhovy_star(i,j) - h / dy * (flux_yp - flux_ym);
		}
	}
	
	
	bound_cond_u(u_con);
	

	
	if(PressureSolver == all_correction || (PressureSolver == once_correction && rks == RKs - 1)){
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
				
//				double h00 = calenthalpy(con00_n);
//				double hp0 = calenthalpy(conp0_n);
//				double hm0 = calenthalpy(conm0_n);
//				double h0p = calenthalpy(con0p_n);
//				double h0m = calenthalpy(con0m_n);
				
				double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
				double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
				double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
				double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
				double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
	
				double RHS = u_con(i,j,3) - 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2)) / u_con(i,j,0) - h * (1 / dx * ( 0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) -  0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * (0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));

				b_p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
				
				double bc_correction = 0.0;
				
				double p_bc = 1.0 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2.0);
				
				double nu_x = vx_max * dx;
				double nu_y = vy_max * dy;
				
				// 左邻居 Dirichlet
				// 常数场: interior = 0, Dirichlet ghost = p_bc
				auto get_x_const = [&](int ii, int jj) -> double {
					bool x_ghost = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
					bool y_ghost = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
					if (x_bc == x_Dirichlet && x_ghost) return p_bc;
					if (y_bc == y_Dirichlet && y_ghost) return p_bc;
					return 0.0;
				};
				
				double p_c = 0.0;   // interior 处 p_const = 0
				
				// === x 方向 5 点 (与 Ax 完全镜像, 把 p_i→p_c, get_x_unknown→get_x_const) ===
				double gradP_im_c = (1.0/dx) * (
					0.5*(p_c + get_x_const(i-1,j))- 0.5*(get_x_const(i-1,j) + get_x_const(i-2,j))
				);
				double gradP_i_c = (1.0/dx) * (
					0.5*(get_x_const(i+1,j) + p_c)-  0.5*(p_c + get_x_const(i-1,j))
				);
				double gradP_ip_c = (1.0/dx) * (
					0.5*(get_x_const(i+2,j) + get_x_const(i+1,j)) - 0.5*(get_x_const(i+1,j) + p_c)
				);
				double Lap_x_c = (1.0/dx) * (
					0.5*(hp0+h00) * 0.5*(gradP_i_c  + gradP_ip_c)-0.5*(hm0+h00) * 0.5*(gradP_im_c + gradP_i_c)
				);
				
				// === y 方向 5 点 ===
				double gradP_jm_c = (1.0/dy) * (
					0.5*(p_c + get_x_const(i,j-1))-  0.5*(get_x_const(i,j-1) + get_x_const(i,j-2))
				);
				double gradP_j_c = (1.0/dy) * (
					0.5*(get_x_const(i,j+1) + p_c)- 0.5*(p_c + get_x_const(i,j-1))
				);
				double gradP_jp_c = (1.0/dy) * (
					0.5*(get_x_const(i,j+2) + get_x_const(i,j+1))- 0.5*(get_x_const(i,j+1) + p_c)
				);
				double Lap_y_c = (1.0/dy) * (
					0.5*(h0p+h00) * 0.5*(gradP_j_c  + gradP_jp_c)- 0.5*(h0m+h00) * 0.5*(gradP_jm_c + gradP_j_c)
				);
				
				double L_const = 1.0/(gamma_plasma - 1) * p_c - h*h * (Lap_x_c + Lap_y_c);
				
				b_p_vect(M_P(i,j)) -= L_const;
				
				p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
			}
		}
		GMRES_Solver(p_vect, b_p_vect, Pressure_Subsystem_AX, u_con, Qe, Qs,h, nxCells*nyCells);
	}
	
	
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			p_solution(i,j) = p_vect(M_P(i,j));
		}
	}
	
	if(x_bc == x_Reflective){
		bound_cond_array(p_solution, x_Transmissive, y_Transmissive);
	}else{
		bound_cond_array(p_solution);
	}
	
	if(x_bc == x_Dirichlet){
		for(int j = 0; j < nyCells + ghost_cells * 2; j++){
			for(int k = 0; k < ghost_cells; k++){
				if(test == testGreshoVortex){
					p_solution(k,j) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
					p_solution(nxCells + ghost_cells + k,j) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
				}else if(test == RPContact){
					p_solution(k,j) = 1e5;
					p_solution(nxCells + ghost_cells + k,j) = 1e5;
				}
			}
		}
	}
	
	if(y_bc == y_Dirichlet){
		for(int i = 0; i < nxCells + ghost_cells * 2; i++){
			for(int ghost = 0; ghost < ghost_cells; ghost++){
				if(test == testGreshoVortex){
					p_solution(i,ghost) = 1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
					p_solution(i,nyCells + ghost_cells + ghost) =  1 / (gamma_plasma * vortex_Mach * vortex_Mach) - 2 + 4 * log(2);
				}else if(test == RPContact){
					p_solution(i,ghost) = 1e5;
					p_solution(i,nyCells + ghost_cells + ghost) = 1e5;
				}
			}
		}
	}	

	
	switch(U_update){
		case DefinitionofU:
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					u_con(i,j,3) = 1 / (gamma_plasma - 1) * p_solution(i,j) + 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2)) / u_con(i,j,0);
				}
			}
			break;
		case ConservationofU:
//			std::cout << "Internal Energy: "<< endl;
			
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					//					double fL, fR,	flux_xm, flux_xp;
					//					double gL, gR, flux_ym, flux_yp;
					
					double e_00 = p_solution(i,j) / ((gamma_plasma - 1) * u_con(i,j,0));
					double e_p0 = p_solution(i + 1,j) / ((gamma_plasma - 1) * u_con(i + 1,j,0));
					double e_m0 = p_solution(i - 1,j) / ((gamma_plasma - 1) * u_con(i - 1,j,0));
					double e_0p = p_solution(i,j + 1) / ((gamma_plasma - 1) * u_con(i,j + 1,0));
					double e_0m = p_solution(i,j - 1) / ((gamma_plasma - 1) * u_con(i,j - 1,0));
					
					double h_00 = e_00 + p_solution(i,j) / u_con(i,j,0);
					double h_p0 = e_p0 + p_solution(i + 1,j) / u_con(i + 1,j,0);
					double h_m0 = e_m0 + p_solution(i - 1,j) / u_con(i - 1,j,0);
					double h_0p = e_0p + p_solution(i,j + 1) / u_con(i,j + 1,0);
					double h_0m = e_0m + p_solution(i,j - 1) / u_con(i,j - 1,0);
					
					
					double fL, fR, flux_xm, flux_xp;
					double gL, gR, flux_ym, flux_yp;
					fL = h_m0 * u_con(i-1,j,1);
					fR = h_00 * u_con(i,j,1);
					
					
					
					
					fL = h_00 * u_con(i,j,1);
					fR = h_p0 * u_con(i+1,j,1);
					
					
					
					
					gL = h_0m * u_con(i,j - 1,2);
					gR = h_00 * u_con(i,j,2);
					
					
					
					
					gL = h_00 * u_con(i,j,2);
					gR = h_0p * u_con(i,j + 1,2);
					
					
					
					u_con(i,j,3) = u_con(i,j,3) - h / (dx) * (flux_xp - flux_xm) - h / (dy) * (flux_yp - flux_ym);
					
					
					double internal_energy = u_con(i,j,3) - (0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2)) / u_con(i,j,0) );
					
//					if(internal_energy < 0){
//						std::cout << internal_energy << " ";
//					}
				}
			}
//			std::cout<<endl;
			break;
	}
	
	bound_cond_u(u_con);
}

double getsoundspeed(double *con){
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	return sqrt(gamma_plasma * pri[3] / pri[0]);
}

void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime){
	bool isPlotTime = false;
	
	std::string testname;
	switch(test){
		case testToroModified:
			testname = "testToroModified";
			break;
		case testToro1:
			testname = "testToro1";
			break;
		case testToro2:
			testname = "testToro2";
			break;
		case testToro3:
			testname = "testToro3";
			break;
		case testToro4:
			testname = "testToro4";
			break;
		case testToro5:
			testname = "testToro5";
			break;
		case RPContact:
			testname = "RPContact";
			break;
		case RPLax:
			testname = "RPLax";
			break;
		case RPTwoShocks:
			testname = "RPTwoShocks";
			break;
		case RPTwoRare:
			testname = "RPTwoRare";
			break;
		case RPFVVEM1:
			testname = "RPFVVEM1";
			break;
		case RPFVVEM2:
			testname = "RPFVVEM2";
			break;
		case RPFVVEM3:
			testname = "RPFVVEM3";
			break;
		case RPFVVEM4:
			testname = "RPFVVEM4";
			break;
		case testExplosion:
			testname = "testExplosion";
			break;
		case testSmoothAcousticWave:
			testname = "SmoothAcousticWave";
			break;
		case testDoubleShearLayer:
			testname = "DoubleShearLayer";
			break;
		case testWoodwardColellaBlastWave:
			testname = "WoodwardColellaBlastWave";
			break;
		case testGreshoVortex:
			testname = "GreshoVortex";
			break;
		case testCollidingPulses:
			testname = "CollidingPulses";
			break;
		case testRiemannLike:
			testname = "testRiemannLike";
			break;
		case convergencetest:
			testname = "convergencetest";
			break;
		case test2d1:
			testname = "test2d1";
			break;
		case test2d2:
			testname = "test2d2";
			break;
		case test2d3:
			testname = "test2d3";
			break;
		case test2d4:
			testname = "test2d4";
			break;
		case KH_instability:
			testname = "KH_instability";
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
		outFile << x_0 + (((double)i)-0.5)*dx << " ";
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
	
	for(int j = ghost_cells; j < ghost_cells+ nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			double varPrint;
			if(var < NUM_VARS){
				varPrint = u_con(i,j,var);
			}else if(var < NUM_VARS + 3){
				if(var - NUM_VARS + 1 == 4){
					varPrint = u_pri(i,j,var - NUM_VARS + 1) / p_ref * p_epsilon;
				}else{
					varPrint = u_pri(i,j,var - NUM_VARS + 1);
				}
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
	string flux_type;
	string update_type;
	string dt_type;
	string diss_type;
	string Laplacian_type;
	
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
	
	switch(U_update){
		case ConservationofU:
			update_type = "conservation";
			break;
		case DefinitionofU:
			update_type = "definition";
			break;
	}
	
	switch(add_pvisc){
		case MaterialSpeedp:
			switch(flux_diss){
				case Global_Diss:
					diss_type = "Global_Laplacian";
					break;
				case Local_Diss:
					diss_type = "Local_Laplacian";
					break;
			}
			switch(add_shock_detector){
				case Yes_add:
					Laplacian_type = "Yes_Detector";
					diss_type = diss_type + "/" + Laplacian_type;
					break;
				case No_add:
					Laplacian_type = "No_Detector";
					diss_type = diss_type + "/" + Laplacian_type;
					break;
			}
			break;
		case Noaddpvis:
			diss_type = "No_Laplacian";
			break;
	}
	
	switch(caltimestep){
		case Fullydt:
			dt_type = "Fullydt";
			break;
		case IMEXdt:
			dt_type = "IMEXdt" + std::to_string(full_step);
			break;
	}
	
	
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(0) << 100*CFL;
	std::string s = ss.str();  
	
	string vortex_MachNumber = "";
	
	if(test == testGreshoVortex){
		if(vortex_Mach == 1e-1){
			vortex_MachNumber = "1e-1";
		}else if(vortex_Mach == 1e-2){
			vortex_MachNumber = "1e-2";
		}else if(vortex_Mach == 1e-3){
			vortex_MachNumber = "1e-3";
		}
	}
	std::string Path = "/Users/chenzefeng/Desktop/Study/CamSC/phd/research/plot/"  + testname +vortex_MachNumber + "/" + diss_type + "/" + update_type + "/" + flux_type + "/"+ std::to_string(nxCells) + "x"+ std::to_string(nyCells) + "_"+ dt_type+"_RK" + std::to_string(RKs) + "_Order" + std::to_string(SpaceOrder) +"_CFL0" + s;
	
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
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(3,"U",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(4,"velX",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(5,"velY",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(6,"p",0,0,0));

	
	
	outFile << "FIELD FieldData " << variablesToPrint.size() << std::endl;
	
	for(unsigned int i = 0;i <variablesToPrint.size(); i++){
		printVTKVariable(outFile,u_con,u_pri,std::get<0>(variablesToPrint[i]),std::get<1>(variablesToPrint[i]));
	}
	
	
	outFile << "POINT_DATA " << (nxCells+1)*(nyCells+1) << std::endl;
	
	return;
}

int main(int argc, char *argv[]){
	
	double M0 = 1.0/ 11.0;
	double L = 2.0 / M0;
	
	switch (test) {
		case testToroModified:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testToro1:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testToro2:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testToro3:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.012;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testToro4:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.035;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testToro5:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.035;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPContact:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.5;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPLax:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.14;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPTwoShocks:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.8;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPTwoRare:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPFVVEM1:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPFVVEM2:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPFVVEM3:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.14;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RPFVVEM4:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.012;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testExplosion:
			x_0 = 0, x_1 = 2, y_0 = 0, y_1 = 2,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testSmoothAcousticWave:					
			x_0 = 0, x_1 = 4, y_0 = 0, y_1 = 4,z_0 = -1, z_1 = 1;
			time_stop = 1;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			// 256 * 256 0.45
			break;
		case testDoubleShearLayer:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 1.8;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			// 256 * 256 0.4
			break;
		case testWoodwardColellaBlastWave:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.038;
			gamma_plasma = 1.4;
			x_bc = x_Reflective;
			y_bc = y_Reflective;
			break;
		case testGreshoVortex:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.4 * M_PI;
			gamma_plasma = 1.4;
//			x_bc = x_Dirichlet;
//			y_bc = y_Dirichlet;
			
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			// 80 * 80, 0.4 
			break;
		case testCollidingPulses:
			x_0 = -L, x_1 = L, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
//			time_stop = 0.378;
			time_stop = 0.841;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testRiemannLike:
			x_0 = 0, x_1 = 1, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.05;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case convergencetest:
			x_0 = 0, x_1 = 10, y_0 = 0, y_1 = 10,z_0 = -1, z_1 = 1;
			time_stop = 1;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case test2d1:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.8;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test2d2:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test2d3:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.23;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test2d4:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.3;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case KH_instability:
			// 120 × 120
			x_0 = -1, x_1 = 1, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 5.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
	}
	
	
	dx = (x_1 - x_0)/(nxPoints-1);
	dy = (y_1 - y_0)/(nyPoints-1);
	dz = (z_1 - z_0)/(nzPoints-1);

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells), u_con_n(nxCells, nyCells), Qe(nxCells, nyCells), Qs(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells),w_yL(nxCells, nyCells), w_yR(nxCells, nyCells);
	Array Bx_p0(nxCells,nyCells), By_0p(nxCells,nyCells), shock_eta(nxCells,nyCells);
	u_initialization(u_pri, Bx_p0, By_0p);
	
	upri2ucon(u_con, u_pri);
	bound_cond_u(u_con);  
	ucon2upri(u_pri, u_con);
	std::vector<double> bRK(RKs);
	Array Ae(RKs, RKs), Ai(RKs, RKs);
	Tensor Hs(RKs, nxCells, nyCells);
	
	vector<double> dt_record;
	vector<double> imexdt_record;
	vector<double> fullydt_record;
	
	const int plotINTERVAL = 100;
	std::vector<double> IntermediatePlotTimes = {
		0.05 * time_stop,
		0.1 * time_stop,
		0.25 * time_stop,
		0.5 * time_stop,
		0.75 * time_stop,
		1.2,
		0.8,
		time_stop,
	};
	bool isFinalTime = false;
	
	IntermediatePlotIntermediatePlot(u_con, u_pri, timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);
	
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
		
		dt = fmin(caldt(u_pri), time_stop - t);
				
		t = t + dt;
		
		upri2ucon(u_con, u_pri);
		
	
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
			bound_cond_u(Qe);
			bound_cond_u(Qs);
			adt =  Ai(rks,rks) * dt;
			
			TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR);
			
			
			Update(u_con, Qe, Qs, w_xL,w_xR,w_yL,w_yR);
			bound_cond_u(u_con);
			
			
			PressureSubsystem(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, adt,rks);
			
			
			
			ucon2upri(u_pri, u_con);
			
			for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
				for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						Hs(rks,i,j,k) = (u_con(i,j,k) - Qs(i,j,k)) / adt;
					}
				}
			}
			bound_cond_u(u_con);
			
		}
		
		ucon2upri(u_pri, u_con);
		IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);
		std::cout<<"[timestep, dt, time] = "  << timestep << " "<< dt << " " << t<< endl;
		ShockDetector(shock_eta, u_con);
		
		
		dt_record.push_back(dt);
		fullydt_record.push_back(fullydt);
		imexdt_record.push_back(imexdt);
		
//		if(timestep > 0){
//			break;
//		}
	}while(t < time_stop);
	
	if(test == convergencetest){
		vector<double> result;
		double errorrho = 0;
		double errorvx = 0;
		double errorvy = 0;
		double errorp = 0;
		double errorrhovx = 0;
		double errorrhovy = 0;
		double errorE = 0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double x = x_0 + (i - ghost_cells + 0.5) * dx;
				double y = y_0 + (j - ghost_cells + 0.5) * dy;
				result = GaussLegendreNumericalIntegration(x, dx, y,dy,ConvergenceTestFunction, 3, time_stop);
				
				
				double con[NUM_VARS],pri[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					pri[k] = result[k];
				}
				ini_pri2con(con, pri);
				errorrho += dx * dy * pow((result[0] - u_pri(i,j,0)), 2.0);
				errorvx += dx * dy * pow((result[1] - u_pri(i,j,1)), 2.0);
				errorvy += dx * dy * pow((result[2] - u_pri(i,j,2)), 2.0);
				errorp += dx * dy * pow((result[3] - u_pri(i,j,3)), 2.0);
				
				errorrhovx += dx * dy * pow((con[1] - u_con(i,j,1)), 2.0);
				errorrhovy += dx * dy * pow((con[2] - u_con(i,j,2)), 2.0);
				errorE += dx * dy * pow((con[3] - u_con(i,j,3)), 2.0);
			}
		}
		errorrho = std::sqrt(errorrho);
		errorvx = std::sqrt(errorvx);
		errorvy = std::sqrt(errorvy);
		errorp = std::sqrt(errorp);
		
		errorrhovx = std::sqrt(errorrhovx);
		errorrhovy = std::sqrt(errorrhovy);
		errorE = std::sqrt(errorE);
		
		std::cout <<  "convergence " <<"ncells: " <<nxCells <<", Mach: "<<vortex_Mach<< endl;
		std::cout <<  "rho err: " <<errorrho <<", vx err: " <<errorvx<<", vy err: " << errorvy<<", p err: " << errorp<< endl;
		std::cout <<  "rhovx err: " <<errorrhovx<<", rhovy err: " << errorrhovy<<", E err: " << errorE<< endl;
	}
	
	ShockDetector(shock_eta, u_con);
	ofstream outFileshock_de("/Users/chenzefeng/Desktop/shock_de.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileshock_de << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << shock_eta(i,j) << endl;
		}
	}
	outFileshock_de.close();
	
	double totalt = 0;
	ofstream outFiledt_fullydt("/Users/chenzefeng/Desktop/dt_fullydt.dat");
	for (size_t i = 0; i < dt_record.size(); i++) {
		totalt += dt_record[i];
		outFiledt_fullydt << totalt << " " << dt_record[i] / fullydt_record[i] << endl;
	}
	outFiledt_fullydt.close();
	
	totalt = 0;
	ofstream outFiledt_imexdt("/Users/chenzefeng/Desktop/dt_imexdt.dat");
	for (size_t i = 0; i < dt_record.size(); i++) {
		totalt += dt_record[i];
		outFiledt_imexdt <<totalt << " " <<  dt_record[i] / imexdt_record[i]<< endl;
	}
	outFiledt_imexdt.close();
	
	
	ofstream outFilerho("/Users/chenzefeng/Desktop/RK_rho.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerho << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,0) << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx("/Users/chenzefeng/Desktop/RK_vx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,1) << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy("/Users/chenzefeng/Desktop/RK_vy.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,2) << endl;
		}
	}
	outFilevy.close();
	
	ofstream outFileE("/Users/chenzefeng/Desktop/RK_E.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileE << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,3) << endl;
		}
	}
	outFileE.close();
	
	
	ofstream outFilep("/Users/chenzefeng/Desktop/RK_p.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilep << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,3)) / p_ref * p_epsilon << endl;
		}
	}
	outFilep.close();
	
	ofstream outFilev("/Users/chenzefeng/Desktop/RK_v_mag.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilev << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2)) << endl;
		}
	}
	outFilev.close();
	
	ofstream outFilep_line("/Users/chenzefeng/Desktop/RK_p_line.dat");
	int j_mid = ghost_cells + nyCells / 2;
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			if(j == j_mid){
				outFilep_line << x_0 + dx * (i - ghost_cells + 0.5) << " "
				<< y_0 + dy * (j - ghost_cells + 0.5) << " "
				<< u_pri(i,j,3) << endl;
			}
		}
	}

	outFilep_line.close();
	
	ofstream outFilerho_line("/Users/chenzefeng/Desktop/RK_rho_line.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			if(j == j_mid){
				outFilerho_line << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,0) << endl;
			}
		}
	}
	outFilerho_line.close();
	
	
}