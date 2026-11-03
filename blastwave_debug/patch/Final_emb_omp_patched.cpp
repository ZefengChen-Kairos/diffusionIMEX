// ============================================================================
// Dual-parameter (eps_M, eps_B) IMEX finite volume scheme for ideal MHD, third order in space and time.
// Final_emb_omp.cpp: OpenMP version (compile with -fopenmp).  Final_emb_serial.cpp is the same code
// without the OpenMP pragmas.
//
// USAGE  Edit the RUN SETTINGS block below (test, mesh, CFL), compile, run in the output directory:
//          g++ -O3 -march=native -fopenmp -std=c++17 -I/usr/include/eigen3 Final_emb_omp.cpp -o mhd
//          OMP_NUM_THREADS=32 ./mhd
//        All output goes to ./plot/<test>/ (created if needed):
//          steplog_<test>.csv     one line per time step: step, t, dt, dt_MHD (fully explicit step from the
//                                 fast magnetosonic speed, same CFL), max|B|, max|B|/max|B_0|,
//                                 ||div B||_L2, ||div B||_inf
//          plt_<test>_<step>.vtk  at t = 0, at the output times of output_times() and at t = T
//          state_<test>_<step>.dat  interior cell averages (i j rho mx my mz E Bx By Bz) at the output times and at T
//        Output times: T/4, T/2, 3T/4 (Orszag-Tang: 0.5, 2, 3 / 2pi, with T = 5/2pi); the time step is shortened
//        to land on them exactly.
//        The divergence is the 4th-order central divergence of the cell-average field (the operator for which
//        div curl4 A = 0 holds exactly); for non-periodic boundaries the two cells next to the boundary are
//        excluded from the norms.
//
// TIME   SA-DIRK(3,4,3), stiffly accurate, partitioned evaluation of the implicit terms.  Per stage:
//          1. magnetic predictor (one GMRES solve), Newton-linearised Lorentz force and induction;
//          2. one pressure system with the constrained transport embedded in the operator, solved by defect
//             correction on the 2nd-order M-matrix operator with the nested Newton pressure floor;
//          3. momentum update; E = p/(g-1) + linearised kinetic and magnetic energy.
//        Time step from the eps-weighted explicit spectrum; first step fully explicit, growth <= 1.5 per step.
// SPACE  explicit: characteristic WENO3-Z + MP3 in the split basis, HLLEM flux; implicit: 4th-order FV
//        operators on cell averages; CT: A_z with HJ-WENO5-Z, B = curl4 A.  In 1D no CT (B_x constant).
// SHOCKS shock-indicator-gated Laplacian (1 - eps_B) |lambda_B| dx in the magnetic predictor.
// EPS    eps_M = 1 - s(M), eps_B = 1 - s(M_A), M0 = 0.15, M1 = 0.4, alpha = 14; face value = max of neighbours.
// ============================================================================
#include <iostream>
#include <climits>
#include <vector>
#include <map>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <chrono>
#define EIGEN_NO_DEBUG   // no index checks inside Eigen (the solver indexes its vectors itself)
#include <Eigen/Dense>
#include <Eigen/SparseLU>
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
enum MagField{UnStaggerCT, PotentialA, NoDivFree};
enum Test{test0, test1, test2, test3, test4, test5, test6, test7,test8,test9,test10, shocktest, convergencetest,test2d0, test1y, testOrszag_Tang, testblastwave, testblastwave1D,testRotor, testFieldLoop, testKHInstability, testdoubleshear, testSmooth2D, testVortex};
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

MagField magfield = UnStaggerCT;
// ============================== RUN SETTINGS ==================================
// The only place to edit: choose the test, the mesh and the CFL number, then compile and run.
Test test = testblastwave;
const int nxCells = 1000;
const int nyCells = 1000;
double CFL = 0.1;
const double KH_MX = 1e-3;             // testKHInstability: Mach number M_x of the stream
const double LOOP_A0 = 1e-3;           // testFieldLoop: amplitude A_0 of the field loop (1e-3, 1e-2, 1e-1, 1)
const double VORTEX_RHO0 = 1.0;        // testVortex: background density
const double VORTEX_BT   = 1.0;        // testVortex: magnetic amplitude
const double VORTEX_VT   = sqrt(2.0*M_PI);   // testVortex: velocity amplitude (1 for the stationary vortex)
const double VORTEX_U0   = 1.0;        // testVortex: uniform advection velocity (0 for the stationary vortex)
const double VORTEX_L    = 5.0;        // testVortex: domain [-L, L]^2 (10 for the stationary vortex)
// ==============================================================================
const int nzCells = 1;
const int ghost_cells = 4;   // o3: limited point<->average conversions need +-2 beyond the reconstruction stencil
const int RKs = 4;          // SA DIRK(3,4,3), Boscheri-Pareschi eq.(79)
int full_step = 1;     // 对齐 CASTOR:特征值 magnetosonic 只在 step 1,之后材料速度(硬编码 timestep==1)
                        // NOTE: was 'const int' but caldt()'s Fullydt branch does full_step=1e8 (pre-existing, unrelated to epsilon changes) -- removed const to fix compile

int SpaceOrder;
x_boundary_condition x_bc;
y_boundary_condition y_bc;

double eps_M_min = 0;    // physical floor (flux/envelope/time step). NOTE: the WENO3 characteristic basis has its OWN floor g_eps_basis_min (fix (3)); eps_M = M^14 underflows to ~0 for M < 0.07
double eps_M_max = 1;
double eps_B_min = 0;    // see eps_M_min note
double eps_B_max = 1;

int    g_gmres_outmax = 50;    // GMRES restart budget (x 25 Krylov vectors); lowered only for the unpreconditioned first try
bool   g_gmres_conv = false;   // set by GMRES_Solver: true if the tolerance was reached
int    g_gmres_iters = 0;      // set by GMRES_Solver: Krylov iterations of the last call
long   g_lu_builds = 0; double g_lu_time = 0.0;   // pressure preconditioner factorisations and their wall time
double GMRES_rtol = 1e-12;   // [o3s] relative tolerance: stop when |M^{-1}(b-Ax)| <= GMRES_rtol*|M^{-1}b|
double GMRESiter = 1e32;
int GMRESrest = 30;

double shock_kappa1 = 0.05;   // sensor threshold: chi > 0 only if h*div(u) < -k1*c_f (velocity jump > 5% of c_f per cell)
double shock_threshold = 0.9;
shock_threshold_Left threshold_Left = val_left;
shock_threshold_Right threshold_Right = zero;

int nxPoints = nxCells + 1, nyPoints = nyCells + 1, nzPoints = nzCells + 1;
double t, time_stop, dt, adt, adt2;
int timestep = 0;
double x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1, B_max_x = 0, B_max_y = 0, U_max_x = 0,U_max_y = 0,a_max_x = 0,a_max_y = 0, a_max = 0;

double dx = (x_1 - x_0)/(nxPoints-1);
double dy = (y_1 - y_0)/(nyPoints-1);
double dz = (z_1 - z_0)/(nzPoints-1);
double gamma_plasma = 5.0/3.0;
double ch_x,ch_y, cp_x, cp_y, cA_x, cA_y, cA, vx_max, vy_max, ch, v_max;

// --- Dirichlet ambient state (blast wave); E/A computed with gamma_plasma at runtime ---
double rho_amb = 1.0, p_amb = 0.1, Bx_amb = 70.0/sqrt(4.0*M_PI);
bool g_bc_prim = false;   // true while bound_cond_u is applied to a PRIMITIVE grid (slot 4 = p, not E)
// [removed: E_amb_val]
 // ambient total energy

// ===================== 压力保正:Casulli 嵌套 Newton(piecewise-linear V) =====================
// 思路:把 pressure subsystem 的对角项 V(p)=p/(gamma-1) 在 p<p_floor 处改成大斜率 Mbig 的直线,
//       于是 V 全局非线性 -> 用 Newton 外层 + GMRES 内层求解;大斜率把 p 弹回 p_floor 之上。
long   g_stage_stamp = 0;   // incremented once per fused-step stage; the AX operators cache
                             // ShockDetector/ucon2upri per stamp (Qe is constant within a stage --
                             // recomputing a full-domain shock detector on every matrix-vector
                             // product was a massive hidden cost)
bool   g_verbose_stages = false;  // silence the per-stage Internal/Solution array dumps (expensive I/O; re-enable for debugging)
bool   g_gmres_newton = false;   // true only during the Jacobian solves of the pressure Newton iteration (its own increment test decides convergence)
bool   floor_on   = true;           // DC restore: positivity ON by default. // 总开关:false 则完全退回你原来的线性单次 GMRES
double p_floor    = 1e-7;          // 压力地板
double Mbig_fac   = 1e9;           // lower bound of M; per solve M = max(Mbig_fac a, (V(p_floor) - b_min)_+ / ((1-theta) p_floor)) (Cor. 1 of 2608.15837)
double theta_floor = 0.9;          // Cor. 1: guarantees min p >= theta p_floor
double Mbig_cur   = 0.0;           // M of the current solve (set in Pressure_NestedNewton)           // 1e6 was QUANTITATIVELY insufficient: the piecewise-linear fixed point is p_floor - deficit/(Mbig*alpha); with deficit ~0.4 (magnetic-step Gibbs injection at low beta) that is NEGATIVE at 1e6. At 1e8 the fixed point sits 4e-9 below the floor: positive. Huge diagonal also decouples these rows -> faster, not slower, convergence.           // 大斜率倍率(相对 1/(gamma-1));GMRES 病态就调小
int    Pnewton_max = 60;   // Mbig corner nonlinearity needs more than 20 iters (observed hit-max on blastwave)           // Newton 最大迭代
double Pnewton_tol = 1e-10;        // Newton 残差容差
Eigen::VectorXd p_freeze;          // 当前 Newton 迭代值,供 Jacobian(=AX)读取斜率

// ---------------- eps-split CT: three face-flux groups ----------------
//   x-face (i,j) == (i+1/2,j) : By-row flux  -> E_z = -F
//   y-face (i,j) == (i,j+1/2) : Bx-row flux  -> E_z = +G
// _rhs  : B* share, (1-eps_B) per face, 4th order, differenced with stagger4
// _imp  : AX wave (1-eps_B)^2 + viscosity (1-eps_B), 2-point
// _expl : explicit eps_B share out of HLLEM_Update, 2-point (finite volume)
struct Array;
struct Grid;

// HLLEM anti-diffusion: LD-only (entropy + Alfven) with the classical entropy amplitude -- the
// complete 8-wave split decomposition was measured to be destroyed by the real delta_k (see paper).

// ---- WENO3 characteristic basis: no component-wise fallback ----
// true  : ALWAYS use the blocked split characteristic basis; only a genuinely
//         non-finite result (NaN/Inf) is guarded. The one unresolved degeneracy
//         (normal field ->0 with a LARGE transverse field) has a completeness
//         error of at most ~2e-6, far below the WENO truncation error at any
//         practical resolution, so falling back there was over-conservative.
// false : also fall back to component-wise when the completeness self-check
//         exceeds g_char_fallback_tol.
                                      // the completeness self-check is counted AND replaced
// eps_M floor used ONLY inside build_split_char_x (the characteristic basis).
// Measured completeness of the blocked basis vs eps_M (loop / KH / generic states):
//   1e-2: err ~1e-13, |L5| ~ 1e2..1e3    1e-4: err ~1e-11, |L5| ~ 1e3..6e4
//   1e-6: err ~1e-7,  |L5| ~ 1e7..6e8    <=1e-10: err O(1), |L5| ~ 1e15..1e19 (dead)
// The flux, the HLL envelope and the time step keep the physical eps_M, so the
// AP time step is unaffected; only the (arbitrary) reconstruction basis is floored.
double g_eps_basis_min = 1.0e-2;

// ---- FINITE-VOLUME HIGH-ORDER CORRECTIONS (both change the numerics) ----
// A finite volume scheme evolves CELL AVERAGES, but the reconstruction returns
// point (or face-average) values and the Riemann solver returns a point flux.
// Every conversion between the two is only second-order accurate unless it is
// corrected, so a third-order reconstruction is capped at second order in 2D.
// The corrections below are the standard convolution/deconvolution pair
//        qbar = q + (1/24) (d2x + d2y) q ,   q = qbar - (1/24)(d2x + d2y) qbar,
// with d2 the UNDIVIDED second difference; each is exact to O(h^4).
//
// g_ic_cell_average : convert the initial condition from the point values that
//        u_initialization writes into genuine cell averages.  MANDATORY for any
//        convergence study -- without it an O(h^2) error is injected at t=0 and
//        the measured order can never exceed 2.  Must be OFF for discontinuous
//        initial data, where the second difference is O(1).
// g_fv_face_quadrature : apply the transverse pair around the Riemann solve, so
//        that the differenced flux is the FACE AVERAGE rather than the flux at
//        the face midpoint.  Needed only for smooth 2D convergence studies; in
//        1D the face is a point and the correction vanishes identically.
//        Must be OFF at discontinuities, for the same reason.
// [o3 2026-09-26] The three switches are GONE: the corrections are always on,
// with the second difference LIMITED by minmod(d2_{i-1}, d2_i, d2_{i+1}).  In smooth
// regions the limited value equals d2 to O(h^3), which keeps third order; across a
// jump the three candidates have different signs and the correction vanishes, so
// no per-test switching is needed (at a discontinuity no scheme exceeds first
// order anyway).  The initial condition is a genuine 3x3 Gauss cell average.
// g_fv_variable_conversion : the reconstruction is carried out in primitive
//        variables, but the stored cell values are averages of the CONSERVATIVE
//        ones, and W(Qbar) is not the cell average of W -- the two differ by
//        O(h^2), which is the third of the conversions that cap the scheme at
//        second order.  With this flag the conversion is routed through the
//        point values,
//              Qbar -> (deconvolve) -> Q_pt -> W(Q_pt) -> (convolve) -> Wbar,
//        each step being exact to O(h^4), and the inverse pair is applied to the
//        reconstructed face values before they enter the Riemann solver.  As
//        for the other two, OFF at discontinuities.
double g_char_fallback_tol = 1.0e-6;

// counts how often the block characteristic basis had to fall back to component-wise
// (normal field ~0 with a large transverse field). Printed once per timestep; should be
// tiny or zero. The OLD cond>1e4 identity fallback fired over the whole Bperp~0 region.
// Quasi-periodic vector potential (fix (2)): with a net mean field the potential
// is A_z = <Bx> y - <By> x + (periodic part), so the periodic ghost fill must add
//   A_z(x, y + L_y) = A_z(x, y) + <Bx> L_y,   A_z(x + L_x, y) = A_z(x, y) - <By> L_x.
// Both jumps are set once from the initial mean field (exactly conserved).
double g_Az_yjump = 0.0, g_Az_xjump = 0.0;
double g_vortex_rho0 = VORTEX_RHO0;
double g_vortex_Bt = VORTEX_BT;
// [o3d] face-average PRODUCTS in the implicit fluxes: avg(f g) = avg(f) avg(g) + (dt^2/12) f_t g_t + O(h^4),
// t = transverse direction along the face.  Without the covariance term every product of two face
// averages (enthalpy x momentum in the energy RHS, enthalpy x pressure gradient in the elliptic
// operator) is only 2nd-order consistent in 2D -- measured: low-Mach vortex rho0 = 1e-2 gave EOC 1.95
// for p (u, B ~3) with the products uncorrected.
const bool g_use_ct = (nxCells != 1 && nyCells != 1);   // 1D runs: no constrained transport (B_n is constant, the transverse field obeys the induction equation, treated like B_z)
const int g_npicard = 3;            // [o3i] pressure solves per stage (momentum re-updated and CT after each): 3 = third order in time
Grid  *g_emb_Qe = nullptr;   // [emb] explicit stage state: fixes the upwind direction of the A transport
long   g_emb_iters = 0, g_emb_solves = 0, g_emb_floor = 0, g_emb_sweeps = 0, g_emb_nc = 0;
bool   g_emb_L2c = false;   // set while the embedded defect correction is running
Grid  *g_ct_A = nullptr, *g_ct_As = nullptr, *g_ct_Ae = nullptr;   // [o3m] main's A, A_s, A_e for the in-Picard CT

using Vec4 = Eigen::Vector4d;
using Vec8 = Eigen::Matrix<double, 8, 1>;
using Mat4 = Eigen::Matrix4d;
using Mat8 = Eigen::Matrix<double, 8, 8>;
inline double V_eos_slope(){ return 1.0/(gamma_plasma - 1.0); }            // 原线性斜率 alpha
inline double V_mod(double p){                                            // piecewise-linear V(p)
    double a = V_eos_slope();
    double M = Mbig_cur > 0.0 ? Mbig_cur : Mbig_fac * a;
    return (p >= p_floor) ? (a*p) : (M*(p - p_floor) + a*p_floor);
}

inline double V_slope(double p){                                          // dV/dp
    double a = V_eos_slope();
    double M = Mbig_cur > 0.0 ? Mbig_cur : Mbig_fac * a;
    return (p >= p_floor) ? a : M;
}

// [removed: V_slope]

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
		const long row = (long)NUM_VARS * (nxCells + 2 * ghost_cells), rows = (long)nRKS * (nyCells + 2 * ghost_cells);
#pragma omp parallel for schedule(static)
		for(long r = 0; r < rows; r++) for(long q = 0; q < row; q++) data[r*row + q] = 0.0;
	}
	
	double operator()(int r,int i, int j, int k)const
	{
		return data[r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int r,int i, int j, int k){
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
		const long row = (long)NUM_VARS * (nxCells + 2 * ghost_cells);
#pragma omp parallel for schedule(static)
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++) for(long q = 0; q < row; q++) data[j*row + q] = 0.0;   // first touch in row bands
	}
	
	double operator()(int i, int j, int k)const
	{
		return data[NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int i, int j, int k){
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
		const long row = (long)(nxCells + 2 * ghost_cells);
#pragma omp parallel for schedule(static)
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++) for(long q = 0; q < row; q++) data[j*row + q] = 0.0;
	}
	
	double operator()(int i, int j)const
	{
		return data[i + j * (nxCells + 2 * ghost_cells)];
	}
	double& operator()(int i, int j){
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
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void f_pressure(double *f, double *con);
void f_magnetic(double *f, double *con);
void g_convective(double *g, double *con);
void g_pressure(double *g, double *con);
void g_magnetic(double *g, double *con);
void ini_con2pri(double *pri, double *con);
void ini_pri2con(double *con, double *pri);
void upri2ucon(Grid &u_con, Grid &u_pri);
void ucon2upri(Grid &u_pri,Grid &u_con);
void half_time(double *halfL, double *halfR, double *L, double *R);
double shock_filter(double val);
void MagneticFlux_Pressure_x(double *Q, double *flux);
void MagneticFlux_Pressure_y(double *Q, double *flux);
// MagneticFlux_Update_Pressure declared near MagSubsystem_UpdateB (epsilon-aware signature)
void RusanovK_x(double *QL, double *QR, double *flux);
void RusanovK_y(double *QL, double *QR, double *flux);
void RusanovK_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
void RusanovK_Pressure_x(double *QL, double *QR, double *flux);
void RusanovK_Pressure_y(double *QL, double *QR, double *flux);
void RusanovK_Update_Pressure(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
void HLLEM_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,
                   Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y, double h);
double calecon(double *con);
double calepri(double *pri);
vector<double> ConvergenceTestFunction(double x);
vector<double> ConvergenceTestExact(double x);
vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<vector<double>(double)> func, int nGaussPoints);
void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p);
void bound_cond_u(Grid &u);
void bound_cond_A(Grid &A);
double caldt(Grid &u_pri, Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y);
double getsoundspeed(double *con);
double getTotalAlfvenSpeed(double *con);
double calMach(double *con);
double calAlfvenMach(double *con);
double f_epsM(double* con);
double f_epsB(double* con);
void bound_cond_epsilon(Array &u);
void getepsgrid(Grid &u_con, Array &eps_M, Array &eps_B,Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y);
void calch(Grid &u_pri);
void calca(Grid &u_pri);
double ini_calch(double *u_pri);
void calcp(Grid &u_pri);
void bound_cond_array(Array &u);
double calenthalpy(double *con);
double calkinetic(double *con);
void MagneticEigenvalues(Grid &u_pri);
void EnergyEigenvalues(Grid &u_pri);
void A_initial(Grid &A);
inline double minmod(double a, double b);
void Unstaggered_CT(Grid& u_con, Grid& A,Grid& A_s, Grid& A_e, double h);
// unknown index of cell (i,j) (ghost cells mapped by the boundary condition); interior cells take the inline fast path
int M_P_bc(int i, int j);
int M_B_bc(int i, int j, int k);
inline int M_P(int i, int j){
	if((unsigned)(i - ghost_cells) < (unsigned)nxCells && (unsigned)(j - ghost_cells) < (unsigned)nyCells)
		return (i - ghost_cells) + nxCells * (j - ghost_cells);
	return M_P_bc(i, j);
}
inline int M_B(int i, int j, int k){
	if((unsigned)(i - ghost_cells) < (unsigned)nxCells && (unsigned)(j - ghost_cells) < (unsigned)nyCells && k >= 5 && k <= 7)
		return (i - ghost_cells) + nxCells * (j - ghost_cells) + (k - 5) * (nxCells * nyCells);
	return M_B_bc(i, j, k);
}
void GMRES_Solver2(Eigen::VectorXd&x,Eigen::VectorXd& b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&,double)> matrix_vector_product, Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim);

// ---- OpenMP reduction helpers for the GMRES inner products (Eigen's dense
// vector ops are SIMD but single-threaded; the CGS dot chain is ~25% of a
// large GMRES solve and every dot is independent) ----
// Inner products: every thread sums a contiguous block in index order and the block sums are added in thread order,
// so the result is the same in every run with the same number of threads (an OpenMP reduction combines the
// partial sums in arrival order).  With one thread this is the plain sequential sum.
#ifdef _OPENMP
#include <omp.h>
#endif
static inline int omp_nthreads_now(){
#ifdef _OPENMP
	return omp_in_parallel() ? 1 : omp_get_max_threads();
#else
	return 1;
#endif
}
static inline int omp_tid(){
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}
// h[i] = v[i] . a for i < m in one pass over a
static void omp_multidot(const Eigen::VectorXd* const* v, const Eigen::VectorXd& a, int m, double* h){
	const int n = (int)a.size(), nt = omp_nthreads_now();
	std::vector<double> part((size_t)nt * m, 0.0);
#pragma omp parallel num_threads(nt)
	{
		const int t = omp_tid(), nth = nt;
		const int q0 = (int)((long long)n * t / nth), q1 = (int)((long long)n * (t + 1) / nth);
		double loc[32]; for(int i = 0; i < m; i++) loc[i] = 0.0;
		const double* ap = a.data();
		for(int q = q0; q < q1; q++){ const double aq = ap[q]; for(int i = 0; i < m; i++) loc[i] += (*v[i])[q] * aq; }
		for(int i = 0; i < m; i++) part[(size_t)t * m + i] = loc[i];
	}
	for(int i = 0; i < m; i++){ double s = 0.0; for(int t = 0; t < nt; t++) s += part[(size_t)t * m + i]; h[i] = s; }
}
inline double omp_dot(const Eigen::VectorXd& a, const Eigen::VectorXd& b){
	const Eigen::VectorXd* v[1] = {&a}; double h; omp_multidot(v, b, 1, &h); return h;
}
inline double omp_nrm(const Eigen::VectorXd& a){ return std::sqrt(omp_dot(a, a)); }
// element-wise vector kernels for GMRES (Eigen's vector expressions run on one thread)
inline void omp_sub(Eigen::VectorXd& out, const Eigen::VectorXd& a, const Eigen::VectorXd& b){
	const int n = (int)a.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++) out[q] = a[q] - b[q];
}
inline void omp_scale(Eigen::VectorXd& out, const Eigen::VectorXd& a, double c){
	const int n = (int)a.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++) out[q] = a[q]*c;
}
inline void omp_axpy(Eigen::VectorXd& y, double c, const Eigen::VectorXd& x){   // y += c x
	const int n = (int)y.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++) y[q] += c*x[q];
}
inline double omp_min(const Eigen::VectorXd& a){
	double m = 1e300; const int n = (int)a.size();
#pragma omp parallel for reduction(min:m) schedule(static)
	for(int q = 0; q < n; q++) m = fmin(m, a[q]);
	return m;
}
inline double omp_dist(const Eigen::VectorXd& a, const Eigen::VectorXd& b){   // |a - b|, summed like omp_multidot
	const int n = (int)a.size(), nt = omp_nthreads_now();
	std::vector<double> part(nt, 0.0);
#pragma omp parallel num_threads(nt)
	{
		const int t = omp_tid();
		const int q0 = (int)((long long)n * t / nt), q1 = (int)((long long)n * (t + 1) / nt);
		double loc = 0.0;
		for(int q = q0; q < q1; q++){ const double d = a[q] - b[q]; loc += d*d; }
		part[t] = loc;
	}
	double s = 0.0; for(int t = 0; t < nt; t++) s += part[t];
	return std::sqrt(s);
}
inline void omp_zero(Eigen::VectorXd& out){
	const int n = (int)out.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++) out[q] = 0.0;
}
inline void omp_copy(Eigen::VectorXd& out, const Eigen::VectorXd& a){
	const int n = (int)a.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++) out[q] = a[q];
}
inline void omp_axpy_sum(Eigen::VectorXd& x, const std::vector<double>& al, const std::vector<Eigen::VectorXd>& v, int m){   // x += sum_{i<m} al[i] v[i]
	const int n = (int)x.size();
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++){ double t = x[q]; for(int i = 0; i < m; i++) t += al[i]*v[i][q]; x[q] = t; }
}

void GMRES_Solver(Eigen::VectorXd &x, Eigen::VectorXd &b,
                   std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, Array&, Array&, Array&, Array&, double)> Ax_func,
                   Grid &u_con, Grid &Qe, Grid &Qs,
                   Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y,
                   double h,const int dim, std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond = {});
double htilde(double hL, double hR, double mL, double mR);
void divB_norms(Grid &u_con, double &L2, double &Linf);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, Array &eps_M, Array &eps_B,int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, Array &eps_M, Array &eps_B,std::string testname, int nstep, bool isFinalTime);

void copy(Grid &state_1, Grid &state_0){
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
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
	// F^p = (0, p, 0, 0, rho*h*vx, 0, 0, 0)^T  (paper eq.(2)/截图; h = e + p/rho)
	double rho = con[0];
	double vx = con[1] / con[0];
	double pri[8];
	ini_con2pri(pri, con);
	double p = pri[4];
	double e = calecon(con);
	double enthalpy = e + p/rho;
	f[0] = 0;
	f[1] = p;
	f[2] = 0;
	f[3] = 0;
	f[4] = rho * enthalpy * vx;
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
	// F^p (y方向) = (0, 0, p, 0, rho*h*vy, 0, 0, 0)^T
	double rho = con[0];
	double vy = con[2] / con[0];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double e = calecon(con);
	double enthalpy = e + p/rho;
	g[0] = 0;
	g[1] = 0;
	g[2] = p;
	g[3] = 0;
	g[4] = rho * enthalpy * vy;
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

void ini_con2pri(double *pri, double *con){
	// convert conservative (rho rhov E) to primitive (rho v p)
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = con[3] / con[0];
	pri[4] = (gamma_plasma - 1) * con[0] * e;
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
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
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
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
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

void ShockDetector(Array &shock_eta, Grid &u_con) {
	// 转换守恒变量到原始变量 (rho, u, v, w, p, Bx, By, Bz)
	static Grid u_pri(nxCells, nyCells);
	ucon2upri(u_pri, u_con);
	
	// 中间数组:先存 chi_tilde (公式 4.40),再做邻居扩展 (公式 4.41)
	static Array shock_eta_tilde(nxCells, nyCells);
	
	// ---- 第一步:计算 chi_tilde (公式 4.40) ----
#pragma omp parallel for collapse(2) schedule(static)
	for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
		for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
			
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
					
					double a2  = gamma_plasma * p_nbr / rho_nbr;
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
#pragma omp parallel for collapse(2) schedule(static)
	for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
		for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
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
	
	bound_cond_array(shock_eta);
}

// [removed: Energy_flux_x]

// [removed: Energy_flux_y]

// [removed: Energy_flux_Update]

// [removed: Pressure_flux_x]

// [removed: Pressure_flux_y]

// [removed: Pressure_flux_Update]

// [removed: RusanovC_x]

// [removed: RusanovC_y]

// [removed: RusanovC_Update]

// ============================================================
// 分裂特征结构的公共构件 (Document 1 §3 / Document 3 §2.2.1)
// HLLEM 和 compute_eigenvalues_x/y 共用,不再各写一份、容易走样。
// ============================================================

// Alfven 特征值对,闭式,只依赖 epsB。x/y 方向共用:
// 传 vn=法向速度(vx 或 vy)、bn2=法向 Alfven 速度平方(Bx²/rho 或 By²/rho)。
inline void SplitAlfvenEigenvalues(double vn, double epsB, double bn2, double &lamAm, double &lamAp){
	double dA = fmax((1.0-epsB)*(1.0-epsB)*vn*vn + 4.0*epsB*epsB*bn2, 0.0);
	double sA = sqrt(dA);
	lamAm = 0.5*((1.0+epsB)*vn - sA);
	lamAp = 0.5*((1.0+epsB)*vn + sA);
}

// 磁声波四次方程系数 (Document 1 eq.7 / Document 3 eq.18)。
// bn2=法向 Alfven 速度平方,bperp2=横向(垂直法向)Alfven 速度平方。
inline void SplitMagnetosonicCoeffs(double vn, double epsM, double epsB, double a2, double bn2, double bperp2,
	double &c3, double &c2, double &c1, double &c0){
		double s1 = (1.0+epsM)*vn, a1 = (1.0+epsB)*vn;
		double s0 = epsM*vn*vn - epsM*epsM*a2;
		double a0 = epsB*vn*vn - epsB*epsB*bn2;
		double q0 = epsM*vn*vn;
		double eB2bp2 = epsB*epsB*bperp2;
		c3 = -(s1+a1);
		c2 = s1*a1 + s0 + a0 - eB2bp2;
		c1 = -(s0*a1 + s1*a0) + eB2bp2*s1;
		c0 = s0*a0 - eB2bp2*q0;
	}

// 可靠的四次方程求根:伴随矩阵 + Eigen::EigenSolver,不用 Newton(理由见上一轮)。
inline Vec4 solve_quartic(double p3, double p2, double p1, double p0){
	Mat4 C = Mat4::Zero();
	C(1,0) = 1.0;  C(2,1) = 1.0;  C(3,2) = 1.0;
	C(0,3) = -p0;  C(1,3) = -p1;  C(2,3) = -p2;  C(3,3) = -p3;
	Eigen::EigenSolver<Mat4> solver(C, false);
	Vec4 r;
	for(int i = 0; i < 4; i++) r(i) = solver.eigenvalues()(i).real();
	std::sort(r.data(), r.data()+4);
	return r;
}

inline void compute_eigenvalues_x(const double* U, double eps_M, double eps_B, Vec8& lambda){
	double rho = U[0];
	double vx = U[1]/rho, vy = U[2]/rho, vz = U[3]/rho;
	double Bx = U[5], By = U[6], Bz = U[7];
	double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bp2 = By*By+Bz*Bz;
	double p  = (gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2);
	double a2 = gamma_plasma*p/rho, bx2 = Bx*Bx/rho, bp2 = Bp2/rho;
	
	lambda(0) = 0.0;
	lambda(1) = vx;
	
	double lamAm, lamAp;
	SplitAlfvenEigenvalues(vx, eps_B, bx2, lamAm, lamAp);
	lambda(2) = lamAm;
	lambda(3) = lamAp;
	
	double c3,c2,c1,c0;
	SplitMagnetosonicCoeffs(vx, eps_M, eps_B, a2, bx2, bp2, c3,c2,c1,c0);
	Vec4 ms = solve_quartic(c3,c2,c1,c0);
	lambda(4) = ms(0); lambda(5) = ms(1); lambda(6) = ms(2); lambda(7) = ms(3);
	
	std::sort(lambda.data(), lambda.data()+8);
}

inline void compute_eigenvalues_y(const double* U, double eps_M, double eps_B, Vec8& lambda){
	double rho = U[0];
	double vx = U[1]/rho, vy = U[2]/rho, vz = U[3]/rho;
	double Bx = U[5], By = U[6], Bz = U[7];
	double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bpt2 = Bx*Bx+Bz*Bz;
	double p  = (gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2);
	double a2 = gamma_plasma*p/rho, by2 = By*By/rho, bpt2 = Bpt2/rho;
	
	lambda(0) = 0.0;
	lambda(1) = vy;
	
	double lamAm, lamAp;
	SplitAlfvenEigenvalues(vy, eps_B, by2, lamAm, lamAp);
	lambda(2) = lamAm;
	lambda(3) = lamAp;
	
	double c3,c2,c1,c0;
	SplitMagnetosonicCoeffs(vy, eps_M, eps_B, a2, by2, bpt2, c3,c2,c1,c0);
	Vec4 ms = solve_quartic(c3,c2,c1,c0);
	lambda(4) = ms(0); lambda(5) = ms(1); lambda(6) = ms(2); lambda(7) = ms(3);
	
	std::sort(lambda.data(), lambda.data()+8);
}

// ================================================================================
// epsilon-HLLEM: delta_k-driven asymptotic-preserving Riemann solver
//   (dual-parameter epsM/epsB MHD -- Document 1 sec.8 / Document 3 sec.3.2 design)
//
// 设计(ROUTE A-analytic, 2026-07): 反耗散在真正的 SPLIT 基上做,但用 ANALYTIC
// 投影,不做 8x8 数值求逆 —— 因此没有条件数悬崖,B_perp->0 干净。
//
// 推导(已数值验证):
//  1) 令 beta_y=B_y/B_perp, beta_z=B_z/B_perp。原始变量 split 矩阵 A_p 的二维子空间
//     是 A_p 的不变子空间(数值:出 S 分量 3.7e-14),其上作用为
//         A_2 = [[u, -k],[-m, epsB*u]],  k = epsB*Bx/rho,  m = epsB*Bx
//     特征值即 lambda_A^{+-}(SplitAlfvenEigenvalues)。
//  2) 右/左特征对没有任何分母:  r_+- = (k, alpha_+-),  l_+- = (m, alpha_+-),
//     alpha_+- = u - lambda_A^{+-},  l.r = m*k + alpha^2 = epsB^2*Bx^2/rho + alpha^2。
//  3) 泛函 dv_perp = bz*dv - by*dw, dB*_perp = bz*dBy - by*dBz 湮灭所有非 Alfven 波
//     (数值:2.4e-16) —— 磁声波的 (v,w),(By,Bz) 都平行于 (by,bz),entropy 是纯 rho,
//     div 模是 Bx。所以投影是直接的,不需要知道另外六个波。
//  4) 右向量取 R_A = k*e_v + alpha*e_B,等于论文 eq.(22) 乘以 (-alpha/B_perp):
//     同一个特征向量,但 1/alpha 和 B_perp 两个奇异因子都被乘掉了。
//     分母 m*k+alpha^2 只在 epsB*Bx=0 且 alpha=0 时为零;由 alpha_+*alpha_- =
//     -epsB^2*c_a^2 知这仅发生在 Bx->0(Alfven 与 entropy 真简并),此时关掉该波。
//  5) entropy 右向量与 eps 无关(论文 eq.21),振幅用正则部分 dRho - dP/a^2。
//     严格的 split entropy 左向量含 (1-epsM)/epsM^2 项(论文自述 eta~1/epsM^2,
//     已数值确认),在 epsM->0 无良定极限(那时显式系统退化为弱双曲,严格双正交本身
//     失去意义),故只保留正则部分:epsM=1 时精确,且无奇异。
//
// 验证:与数值 L=R^{-1} 在良态态上的贡献 reldiff 6.1e-12(23392 态);eps=1 时约化
// 到 Roe-Balsara。magnetosonic 仍交给 HLL 包络(不参与反耗散),故不需要四次方程根。
// HLL 底通量用 SPLIT 显式通量 F~c;安全包络仍用全 MHD Roe fast speed(耗散上界)。
// ================================================================================

// dQ/dq: primitive q=(rho,u,v,w,p,Bx,By,Bz) -> conservative Q=(rho,rho u,...,U,B).
// Maps a split-basis primitive right eigenvector column into conservative space so
// the flux anti-diffusion F_HLL - coeff*sum delta_k alpha_k (M R_k^prim) is applied
// on conserved variables, exactly as eq.(32) requires.
inline void prim_eigvec_to_cons(const double rho, const double vx, const double vy, const double vz,
                                const double Bx, const double By, const double Bz,
                                const Vec8& rp, double* rc){
	const double q2 = vx*vx+vy*vy+vz*vz;
	const double dr=rp(0), du=rp(1), dv=rp(2), dw=rp(3), dp=rp(4), dBx=rp(5), dBy=rp(6), dBz=rp(7);
	rc[0] = dr;
	rc[1] = vx*dr + rho*du;
	rc[2] = vy*dr + rho*dv;
	rc[3] = vz*dr + rho*dw;
	// U = p/(g-1) + 1/2 rho|v|^2 + 1/2 |B|^2  ->  dU = dp/(g-1) + 1/2 q2 dr + rho(v.du) + (B.dB)
	rc[4] = dp/(gamma_plasma-1.0) + 0.5*q2*dr + rho*(vx*du+vy*dv+vz*dw) + (Bx*dBx+By*dBy+Bz*dBz);
	rc[5] = dBx;
	rc[6] = dBy;
	rc[7] = dBz;
}
// ================================================================================

// F~c = F^c + epsM*F^p + epsB*F^B (显式子系统通量,截图公式)
void SplitExplicitFlux_x(double *con, double epsM, double epsB, double *f){
	double fC[NUM_VARS], fP[NUM_VARS], fB[NUM_VARS];
	f_convective(fC, con);
	f_pressure(fP, con);   // 已修复: F^p 能量行是 rho*h*vx,不是 (U+p)*vx
	f_magnetic(fB, con);
	for(int k = 0; k < NUM_VARS; k++) f[k] = fC[k] + epsM * fP[k] + epsB * fB[k];
}
// [removed: SplitExplicitFlux_y]

// Alfven 特征值对,闭式,x/y 共用(HLLEM 和 compute_eigenvalues_x/y 共享同一个实现,定义见上文)

// Exact extreme magnetosonic speeds of the SPLIT system, by monotone Newton on the
// quartic (18). For a monic real-rooted quartic, p>0 and p''>0 to the right of the
// largest root, so Newton from a Cauchy bound descends monotonically to lambda_f+;
// symmetrically from the left for lambda_f-. ~10 scalar iterations (cap 200), no eigen-solver.
// Verified: matches the companion-matrix EigenSolver to 5.5e-13 over 60k states, and
// [lamFmin,lamFmax] envelopes the Alfven and entropy waves in 60000/60000 states
// (only the div mode lambda=0 can fall outside, and it carries zero flux: F~c[5]=0,
// so it is correctly excluded -- same convention as standard MHD wave-speed estimates).
inline void SplitFastSpeeds(double vn, double epsM, double epsB,
                            double a2, double bn2, double bperp2,
                            double &lamFmin, double &lamFmax){
	double c3,c2,c1,c0;
	SplitMagnetosonicCoeffs(vn, epsM, epsB, a2, bn2, bperp2, c3,c2,c1,c0);
	auto P  = [&](double x){ return (((x+c3)*x+c2)*x+c1)*x+c0; };
	auto dP = [&](double x){ return ((4.0*x+3.0*c3)*x+2.0*c2)*x+c1; };
	// Starting point: Fujiwara bound 2 max(|c3|, |c2|^(1/2), |c1|^(1/3), |c0/2|^(1/4)) >= every root.  It scales like the
	// wave speeds, so Newton needs ~10 steps; the Cauchy bound 1 + max|c_i| scales like speed^4 and left the iteration
	// unconverged at its cap for large speeds (blast-like states: 0.6% of the interfaces, envelope up to 20% too wide).
	const double R = 2.0*fmax(fmax(fabs(c3),sqrt(fabs(c2))), fmax(cbrt(fabs(c1)),sqrt(sqrt(0.5*fabs(c0))))) * (1.0 + 1e-12) + 1e-300;
	double x = R;
	for(int k=0;k<200;k++){ double d=dP(x); if(d<=0.0) break; double s=P(x)/d;
		if(!(s>0.0)) break; x-=s; if(fabs(s)<=1e-13*(fabs(x)+1.0)) break; }
	lamFmax = x;
	x = -R;
	for(int k=0;k<200;k++){ double d=dP(x); if(d>=0.0) break; double s=P(x)/d;
		if(!(s<0.0)) break; x-=s; if(fabs(s)<=1e-13*(fabs(x)+1.0)) break; }
	lamFmin = x;
}

// COMPLETE (all-8-wave) split anti-diffusion. Defined after build_split_char_x
// (it needs the SplitChar8 blocked basis, which lives further down); forward-declared
// here so HLLEM_x can call it. Returns false if the basis came out degenerate, in
// which case HLLEM_x falls back to plain HLL for that interface.

inline double HLLEM_lam_minus(double lam){ return 0.5 * (lam - fabs(lam)); }
inline double HLLEM_lam_plus (double lam){ return 0.5 * (lam + fabs(lam)); }
inline double HLLEM_clamp01(double x){ return fmax(0.0, fmin(1.0, x)); }

// [o3l] HLL base-flux weights of the last HLLEM_x call (w_L F(Q_L) + w_R F(Q_R) is the HLL average; the
// rest of the numerical flux is dissipation + anti-diffusion).  HLLEM_Update needs them to split the
// numerical flux into its F^c, F^p, F^B parts BEFORE the transverse convolution.
static thread_local double g_hll_wL = 0.5, g_hll_wR = 0.5;
int g_hllemflat = 2, g_hllemoff = 0;   // [fix C] HLLEM anti-diffusion: flattened by (1 - chi_face) (2: + pressure-jump sensor) / switched off (plain HLL)
static thread_local double g_face_chi = 0.0;   // shock indicator of the current face (set by HLLEM_Update)
Grid* g_hllem_Qe = nullptr;              // explicit stage state for the shock indicator
extern long g_rho_guard; long g_ent_lim = 0;   // counters: density point-value guard, entropy anti-diffusion limiter
double g_env_cap = 0.0;   // [W] time-step speed of the current sweep direction (a_max_x / a_max_y), set by HLLEM_Update
long g_env_wide = 0;      // [W] counter: faces whose envelope was widened
void HLLEM_x(double *QL, double *QR, double epsM, double epsB, double *flux){
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	double rhoL=priL[0], vxL=priL[1], vyL=priL[2], vzL=priL[3], pL=priL[4];
	double BxL=priL[5], ByL=priL[6], BzL=priL[7];
	double rhoR=priR[0], vxR=priR[1], vyR=priR[2], vzR=priR[3], pR=priR[4];
	double BxR=priR[5], ByR=priR[6], BzR=priR[7];
	
	// 显式分裂通量(构造 HLL 底通量,Document 3 eq.32)
	double FcL[NUM_VARS], FcR[NUM_VARS];
	SplitExplicitFlux_x(QL, epsM, epsB, FcL);
	SplitExplicitFlux_x(QR, epsM, epsB, FcR);
	
	// ---- Cargo-Gallice Roe 平均 (Document 3 eq.30) ----
	double sqrtL = sqrt(rhoL), sqrtR = sqrt(rhoR), sqsum = sqrtL + sqrtR;
	double rho_t = sqrtL * sqrtR;
	double vx_t = (sqrtL*vxL + sqrtR*vxR) / sqsum;
	double vy_t = (sqrtL*vyL + sqrtR*vyR) / sqsum;
	double vz_t = (sqrtL*vzL + sqrtR*vzR) / sqsum;
	double Bx_t = (sqrtR*BxL + sqrtL*BxR) / sqsum;   // B 用反权重
	double By_t = (sqrtR*ByL + sqrtL*ByR) / sqsum;
	double Bz_t = (sqrtR*BzL + sqrtL*BzR) / sqsum;
	
	double B2L = BxL*BxL+ByL*ByL+BzL*BzL, B2R = BxR*BxR+ByR*ByR+BzR*BzR;
	double H_L = (QL[4] + pL + 0.5*B2L) / rhoL;
	double H_R = (QR[4] + pR + 0.5*B2R) / rhoR;
	double H_t = (sqrtL*H_L + sqrtR*H_R) / sqsum;
	
	double q2_t = vx_t*vx_t + vy_t*vy_t + vz_t*vz_t;
	double B2_t = Bx_t*Bx_t + By_t*By_t + Bz_t*Bz_t;
	double a2_t = (gamma_plasma - 1.0) * (H_t - 0.5*q2_t - B2_t/rho_t);
	// UNPHYSICAL-STATE GUARD (root cause of the Brio-Wu blowup, found by flux
	// sentinel): a slight explicit overshoot at the compound wave leaves an
	// interface state with p ~ -0.013. The Roe-averaged a2_t then comes out
	// NEGATIVE, and the old fmax(a2_t,1e-12) floor turned alpha_C = dRho - dP/a2_t
	// into dP/1e-12 ~ 1e10 -- the "guard" itself was the amplifier that produced
	// the 1e10 mass flux and rho=-1e9. Correct handling (standard practice):
	// if either input pressure is non-positive or the Roe a2 is non-physical,
	// DISABLE the LD anti-diffusion entirely and fall back to plain HLL, which
	// is more dissipative but unconditionally robust. A degenerate state has no
	// business feeding a fine-grained characteristic decomposition anyway.
	bool ld_ok = (pL > 0.0) && (pR > 0.0) && (a2_t > 1e-10);
	a2_t = fmax(a2_t, 1e-12);   // keep positive for the wave-speed estimate below
	
	// ---- SPLIT Einfeldt envelope (was: full-MHD eps=1 Roe fast speed) ----
	// 论文 §3.2 要求 "as eps->0, delta_k->0 and the flux reverts to plain HLL with
	// O(|v|) dissipation" —— 这只有在包络本身用 SPLIT 波速时才成立。用 eps=1 的
	// c_f 会留下 O(c_f) 的耗散,并把 LLF 稳定性限制压成 CFL <~ eps_min。
	//
	// 必须对 L、Roe、R 三个态都取极值(Einfeldt 1988),不能只用 Roe 态:
	// 旧代码 min(uL,uR,~u-~cf) 只带 L/R 的物质速度,靠 eps=1 的 c_f 过大来兜底;
	// 换成紧的 split 速度后该缺陷立刻暴露 —— 实测 blast 界面(rho=1, p 1000 vs 0.1,
	// Bx=70/sqrt(4pi), eps=0.1)Roe 给 2.646 而热侧 L 需要 3.742,低估 1.414 倍,
	// HLL 的保正界被破坏 -> 负密度(ts=2 就 rho=-11)。
	double b2_t  = B2_t / rho_t;
	double bx2_t = Bx_t*Bx_t / rho_t;
	double bp2_t = fmax(b2_t - bx2_t, 0.0);

	double a2L = gamma_plasma*pL/rhoL, a2R = gamma_plasma*pR/rhoR;
	double bx2L = BxL*BxL/rhoL,        bx2R = BxR*BxR/rhoR;
	double bp2L = (ByL*ByL+BzL*BzL)/rhoL, bp2R = (ByR*ByR+BzR*BzR)/rhoR;

	double loT,hiT, loL,hiL, loR,hiR;
	double sL, sR;
	{	// SPLIT Einfeldt envelope: extreme eigenvalues of the split system over {L, Roe, R}
		SplitFastSpeeds(vx_t, epsM, epsB, a2_t,          bx2_t, bp2_t, loT, hiT);
		SplitFastSpeeds(vxL,  epsM, epsB, fmax(a2L,0.0), bx2L,  bp2L,  loL, hiL);
		SplitFastSpeeds(vxR,  epsM, epsB, fmax(a2R,0.0), bx2R,  bp2R,  loR, hiR);
		sL = fmin(fmin(loL, loT), loR);
		sR = fmax(fmax(hiL, hiT), hiR);
	}
	// [W] Strong density jump, J = |drho|/min(rho) > 1/2 (the J of the entropy anti-diffusion limiter below).  The split
	// envelope is O(|u| + eps c): at a stagnation point of a low-Mach region it vanishes, and so does the dissipation of the
	// plain HLL flux the limiter falls back to (measured: Orszag-Tang X-point, eps_M = 0, density holes grow until the point
	// values fail).  From J = 1/2 to J = 1 the envelope is widened linearly to the unsplit (eps = 1) Einfeldt envelope over
	// {L, Roe, R}, capped by the speed that sets the time step (the split spectral radius of this direction): the time step
	// is unchanged, and in a low-Mach flow the cap is O(|u|), so the dissipation keeps the asymptotic scaling.  Where
	// eps = 1 the split envelope already is the unsplit one and nothing changes.
	{	const double Jw = fabs(rhoR - rhoL)/fmin(rhoL, rhoR);
		if(Jw > 0.5 && g_env_cap > 0.0){
			const double th = fmin(1.0, (Jw - 0.5)/0.5);
			auto cfast = [](double a2, double bn2, double bp2){ const double s2 = a2 + bn2 + bp2; return sqrt(0.5*(s2 + sqrt(fmax(s2*s2 - 4.0*a2*bn2, 0.0)))); };
			const double cT = cfast(a2_t, bx2_t, bp2_t), cL = cfast(fmax(a2L,0.0), bx2L, bp2L), cR = cfast(fmax(a2R,0.0), bx2R, bp2R);
			const double sL1 = fmax(fmin(fmin(vxL - cL, vx_t - cT), vxR - cR), -g_env_cap);
			const double sR1 = fmin(fmax(fmax(vxL + cL, vx_t + cT), vxR + cR),  g_env_cap);
			if(th*sL1 < sL || th*sR1 > sR){
#pragma omp atomic
				g_env_wide++;
			}
			sL = fmin(sL, th*sL1); sR = fmax(sR, th*sR1);
		}
	}
	
	// 完全单侧的情况直接退化成分裂通量的迎风值(Document 3 eq.32)
	if(sL >= 0.0){ g_hll_wL = 1.0; g_hll_wR = 0.0; for(int k=0;k<NUM_VARS;k++) flux[k]=FcL[k]; return; }
		if(sR <= 0.0){ g_hll_wL = 0.0; g_hll_wR = 1.0; for(int k=0;k<NUM_VARS;k++) flux[k]=FcR[k]; return; }
			g_hll_wL = sR/(sR - sL); g_hll_wR = -sL/(sR - sL);
			
			double FHLL[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++)
				FHLL[k] = (sR*FcL[k] - sL*FcR[k] + sL*sR*(QR[k]-QL[k])) / (sR - sL);
			
			double coeff_hll = sL*sR/(sR-sL);
			
			// ==================== FULLY SPLIT: ALL 8 WAVES ====================
			// Complete split characteristic decomposition via the SAME blocked basis
			// WENO3 uses (SplitChar8). Because the sum is COMPLETE, the 1/epsM part of
			// the exact split entropy projection is cancelled by the magnetosonic waves
			// -- so no eps=1 amplitude is needed anywhere. delta_k auto-excludes the
			// extreme waves (lambda_k = s_L or s_R  =>  delta_k = 0).
			// ============ else: LD-only (entropy + Alfven), classical entropy ============
			
			// ---- LD 核心反耗散:entropy + Alfven+/- 在 SPLIT 基 R(epsM,epsB) 上 (ROUTE A) ----
			// ---- beta 旋转(B_perp=0 时取任意固定单位向量,公式保持有限)----
			double Bperp_t = sqrt(By_t*By_t + Bz_t*Bz_t);
			double beta_y, beta_z;
			const double tiny_bp = 1e-12;
			if(Bperp_t > tiny_bp){ beta_y = By_t/Bperp_t; beta_z = Bz_t/Bperp_t; }
			else { beta_y = 1.0; beta_z = 0.0; }

			// ---- entropy: R_ent 与 eps 无关(eq.21);振幅用 CLASSICAL dRho - dP/a^2 ----
			// 为什么不用精确的 split entropy 左向量(它带 1/(epsM a^2) 和
			// (1-epsM)/epsM^2 项,公式见 Prop.):因为 HLLEM 按设计只对 LD 波做反耗散
			// (Dumbser-Balsara 2016),而精确 split entropy 振幅的那些大项在 COMPLETE
			// 8 波分解里是被磁声波精确抵消的。实测(rho=1, a=26.46, Bx=19.75,
			// dRho=0, dP=-999.9):
			//     epsM      1.0     0.5     0.2     0.1
			//     alpha_ent 1.43    2.86    7.14   14.28   <- 精确 split,~1/epsM
			//     all-8 sum 0.000   0.000   0.000   0.000  <- 完整分解恒等于 dRho
			//     classical 1.43    1.43    1.43    1.43
			// LD-only 截断后 1/epsM 无人抵消,而 R_ent[0]=1 使它直接进密度通量 ->
			// blast 在 ts=2~3 出现 rho_min = -11 ... -30(与 14.28 同量级)。
			// Alfven 没有这个问题:它的二维子空间真不变、对偶泛函湮灭所有其他波,
			// 故 Alfven 的 split 投影截断安全,entropy 的不安全。
			double dRho = rhoR - rhoL, dP = pR - pL;
			double alpha_C = dRho - dP/a2_t;
			Vec8 Rent_p; Rent_p.setZero(); Rent_p(0) = 1.0;

			// ---- Alfven: 2x2 解析特征对,无分母 ----
			// The 2x2 block is A2 = [[u,-k],[-m, epsB*u]], k = epsB*Bx/rho,
			// m = epsB*Bx.  Generic case: r_pm = (k, alpha_pm), normalisation
			// l.r = m*k + alpha^2 = epsB^2 c_a^2 + alpha^2.
			//
			// Bx -> 0 IS NOT A PHYSICAL COALESCENCE, it is a degeneracy of this
			// PARAMETRISATION.  At Bx = 0 one has k = m = 0 and A2 becomes
			// DIAGONAL, diag(u, epsB*u): v_perp and B*_perp decouple exactly, with
			// the two speeds u and epsB*u, which are DISTINCT whenever epsB != 1.
			// The eigenvectors are simply e_v and e_B.  What fails is r = (k,alpha):
			// the branch with alpha = 0 (lambda = u, since alpha_+ = u - lambda_A+
			// -> 0) collapses to the zero vector, while the other branch,
			// alpha_- = (1-epsB)u, still correctly returns e_B.  The previous
			// threshold guard tested the two denominators separately, so it kept
			// the minus branch (den_- stays O(1)) and silently ZEROED the plus
			// branch -- discarding the lambda = u, pure-v_perp wave, which is
			// perfectly well defined there.  It was also INCONSISTENT with the
			// reconstruction, whose blocked basis already takes the diagonal limit
			// analytically (see build_split_char_x, S.decoupA).
			//
			// The switch is continuous, not a patch: as Bx -> 0 the generic branch
			// has alpha_+ = O(Bx^2) and den_+ ~ epsB^2 Bx^2/rho, so
			// alpha_A+ * R_A+ -> dv_perp * e_v, exactly the analytic branch below.
			//
			// The case is not exceptional: the blast wave has By == 0 everywhere at
			// t = 0, so the y-direction solve sits exactly on it, and Orszag-Tang
			// and Kelvin-Helmholtz spend extended regions there.  We therefore use
			// the same analytic branch and the same threshold as the reconstruction.
			double k_A = epsB*Bx_t/rho_t;    // k = epsB*Bx/rho
			double m_A = epsB*Bx_t;          // m = epsB*Bx
			double lambda_Am, lambda_Ap;
			SplitAlfvenEigenvalues(vx_t, epsB, bx2_t, lambda_Am, lambda_Ap);
			double alp_Ap = vx_t - lambda_Ap;   // alpha_+
			double alp_Am = vx_t - lambda_Am;   // alpha_-

			// 湮灭型泛函:自动只挑出 Alfven 分量(见函数头 3)
			double dv_perp = beta_z*(vyR-vyL) - beta_y*(vzR-vzL);
			double dB_perp = beta_z*(ByR-ByL) - beta_y*(BzR-BzL);

			double alpha_Ap, alpha_Am;
			Vec8 RAp_p, RAm_p;
			RAp_p.setZero(); RAm_p.setZero();

			const double spd_A   = fabs(vx_t) + sqrt(fmax(a2_t,0.0)) + sqrt(fmax(b2_t,0.0)) + 1e-300;
			const double ca_eff  = epsB*fabs(Bx_t)/sqrt(rho_t);
			const bool   decoupA = !(ca_eff > 1.0e-8*spd_A);

			if(decoupA){
				// A2 = diag(u, epsB*u).  lambda = u carries pure v_perp,
				// lambda = epsB*u carries pure B*_perp; both amplitudes are the
				// annihilating functionals themselves, with unit normalisation.
				lambda_Ap = vx_t;          // the v_perp mode
				lambda_Am = epsB*vx_t;     // the B*_perp mode
				alpha_Ap  = dv_perp;
				alpha_Am  = dB_perp;
				RAp_p(2) =  beta_z; RAp_p(3) = -beta_y;   // e_v
				RAm_p(6) =  beta_z; RAm_p(7) = -beta_y;   // e_B
			}else{
				// l.r = m*k + alpha^2 = epsB^2 c_a^2 + alpha^2 > 0 away from Bx = 0
				double den_Ap = m_A*k_A + alp_Ap*alp_Ap;
				double den_Am = m_A*k_A + alp_Am*alp_Am;
				alpha_Ap = (m_A*dv_perp + alp_Ap*dB_perp)/den_Ap;
				alpha_Am = (m_A*dv_perp + alp_Am*dB_perp)/den_Am;
				// R_A = k*e_v + alpha*e_B (primitive),无 1/alpha、无 1/B_perp
				RAp_p(2)= k_A*beta_z; RAp_p(3)= -k_A*beta_y;
				RAp_p(6)= alp_Ap*beta_z; RAp_p(7)= -alp_Ap*beta_y;
				RAm_p(2)= k_A*beta_z; RAm_p(3)= -k_A*beta_y;
				RAm_p(6)= alp_Am*beta_z; RAm_p(7)= -alp_Am*beta_y;
			}

			// ---- primitive -> conservative ----
			double Rc_ent[NUM_VARS], Rc_Am[NUM_VARS], Rc_Ap[NUM_VARS];
			prim_eigvec_to_cons(rho_t,vx_t,vy_t,vz_t,Bx_t,By_t,Bz_t, Rent_p, Rc_ent);
			prim_eigvec_to_cons(rho_t,vx_t,vy_t,vz_t,Bx_t,By_t,Bz_t, RAm_p,  Rc_Am);
			prim_eigvec_to_cons(rho_t,vx_t,vy_t,vz_t,Bx_t,By_t,Bz_t, RAp_p,  Rc_Ap);

			// ---- AP 缩放 delta_k:split 特征值 ----
			double lambda_C = vx_t;   // entropy 速度,split 与全 MHD 同

			double delta_C  = HLLEM_clamp01(1.0 - HLLEM_lam_minus(lambda_C)/sL  - HLLEM_lam_plus(lambda_C)/sR);
			double delta_Ap = HLLEM_clamp01(1.0 - HLLEM_lam_minus(lambda_Ap)/sL - HLLEM_lam_plus(lambda_Ap)/sR);
			double delta_Am = HLLEM_clamp01(1.0 - HLLEM_lam_minus(lambda_Am)/sL - HLLEM_lam_plus(lambda_Am)/sR);

			// Entropy-wave anti-diffusion across strong density jumps.  delta_C = 1 - lam-/sL - lam+/sR -> 1 where the
			// contact speed u -> 0, so at a stagnation point the density jump keeps no dissipation at all, and the split
			// envelope there is only O(|u|).  At the magnetic null of Orszag-Tang (v = B = 0, compressive stagnation flow)
			// a grid-scale density structure then grows from t ~ 0.3 (rho 0.9 next to 10 at 1024^2) until the point values
			// of rho fail.  The anti-diffusion is therefore reduced linearly from a relative jump |drho|/min(rho) of 1/2
			// and switched off from 1 on: never active in smooth flow (|drho|/rho = O(h)), plain HLL across strong contacts,
			// with the envelope widened by [W] above.
			{ const double J = fabs(rhoR - rhoL)/fmin(rhoL, rhoR); if(J > 0.5){ delta_C *= fmax(0.0, (1.0 - J)/0.5);
#pragma omp atomic
				g_ent_lim++; } }
			if(g_hllemoff){ delta_C = 0.0; delta_Ap = 0.0; delta_Am = 0.0; }
			else if(g_hllemflat){ const double fl = fmax(0.0, 1.0 - g_face_chi); delta_C *= fl; delta_Ap *= fl; delta_Am *= fl; }   // [fix C] flattener
			double coeff = sL*sR/(sR-sL);
			if(!ld_ok){ for(int k=0;k<NUM_VARS;k++) flux[k]=FHLL[k]; return; }   // degenerate state: plain HLL
			for(int k = 0; k < NUM_VARS; k++){
				double antidiff = delta_C*alpha_C*Rc_ent[k] + delta_Ap*alpha_Ap*Rc_Ap[k] + delta_Am*alpha_Am*Rc_Am[k];
				flux[k] = FHLL[k] - coeff*antidiff;
			}
		}
		
// y方向:循环置换 (x,y,z)->(y,z,x) 借 HLLEM_x 实现,不手推第二套符号
void HLLEM_y(double *QL, double *QR, double epsM, double epsB, double *flux){
	double QLr[NUM_VARS] = {QL[0], QL[2], QL[3], QL[1], QL[4], QL[6], QL[7], QL[5]};
	double QRr[NUM_VARS] = {QR[0], QR[2], QR[3], QR[1], QR[4], QR[6], QR[7], QR[5]};
	double fr[NUM_VARS];
	HLLEM_x(QLr, QRr, epsM, epsB, fr);
	// fr 排布: 0=质量,1=my通量,2=mz通量,3=mx通量,4=能量,5=By通量,6=Bz通量,7=Bx通量
	flux[0]=fr[0]; flux[1]=fr[3]; flux[2]=fr[1]; flux[3]=fr[2]; flux[4]=fr[4];
	flux[5]=fr[7]; flux[6]=fr[5]; flux[7]=fr[6];
}

// ================================================================================
// HLLEM_Update: replaces RusanovC_Update in the main loop. Same read pattern
// (w_xR(i,j)=left state of face i+1/2, w_xL(i+1,j)=right state -- matches
// WENO3Reconstruction's output convention, checked directly against the code
// before wiring this in), but uses HLLEM_x/y (SplitExplicitFlux-based HLL flux
// + AP-preserving delta_k anti-diffusion) instead of plain Rusanov, and passes
// the FULL split (epsM,epsB) at each face -- this intermediate momentum/energy
// state is what MagSubsystem_UpdateB's implicit B-solve reads as input.
// ================================================================================
// Whole-grid convolution / deconvolution between point values and cell averages.
//   sgn = +1 : point value  -> cell average
//   sgn = -1 : cell average -> point value
static inline double mm3(double a, double b, double c){
	if(a > 0.0 && b > 0.0 && c > 0.0) return std::min(a, std::min(b, c));
	if(a < 0.0 && b < 0.0 && c < 0.0) return std::max(a, std::max(b, c));
	return 0.0;
}
// S[k] = max |G(i,j,k)| over i0 <= i < i1, j0 <= j < j1, all k (parallel; max is order independent)
static void grid_absmax(const Grid& G, int i0, int i1, int j0, int j1, double* S){
	for(int k = 0; k < NUM_VARS; k++) S[k] = 0.0;
#pragma omp parallel
	{
		double loc[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) loc[k] = 0.0;
#pragma omp for schedule(static) nowait
		for(int j = j0; j < j1; j++){
			double m[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) m[k] = loc[k];
			for(int i = i0; i < i1; i++){ const double* c = G.data + NUM_VARS*(i + j*(G.nxCells + 2*ghost_cells)); for(int k = 0; k < NUM_VARS; k++) m[k] = std::max(m[k], std::fabs(c[k])); }   // (= fmax: m is never NaN)
			for(int k = 0; k < NUM_VARS; k++) loc[k] = m[k];
		}
#pragma omp critical
		for(int k = 0; k < NUM_VARS; k++) S[k] = fmax(S[k], loc[k]);
	}
}
// limited undivided second difference at the centre point of a 5-point line
// [o3b] The minmod gate is only allowed to fire when at least one candidate second
// difference is large on the scale of the variable (|d2| >= kappa*h*S): a genuine jump
// has |d2| ~ S, smooth data |d2| ~ h^2 S, grid-scale noise |d2| ~ delta << h S.  With the
// gate always armed (mm3 on everything) the on/off switching between neighbouring faces
// turned grid-scale density noise into a growing drift on the low-Mach KH (measured:
// 2.8e-3 -> 2.8e-4 at t=81, 64x32, once the conversions are linear there).
// [o3b] eps_M blend of the conversions: limited where the flow is explicit (shocks live at
// M >= 1, eps_M = 1), LINEAR in the implicit (low-Mach) share.  A minmod that fires on an
// unresolved low-Mach shear layer (2-cell KH layer: |d2 F| ~ S_F) switches on and off between
// neighbouring faces and injects grid-scale momentum; measured on the KH 64x32, t=81: limited
// flux convolution -> density drift 2.8e-3, linear -> 2.8e-4.
Array *g_epsM_cell = nullptr, *g_epsM_fx = nullptr, *g_epsM_fy = nullptr;
static inline double d2lin(double qm1, double q0, double qp1){ return qp1 - 2.0*q0 + qm1; }
static inline double d2lim(double qm2, double qm1, double q0, double qp1, double qp2, double thr){
	const double a = q0 - 2.0*qm1 + qm2, b = qp1 - 2.0*q0 + qm1, c = qp2 - 2.0*qp1 + q0;
	if(fabs(a) < thr && fabs(b) < thr && fabs(c) < thr) return b;
	return mm3(a, b, c);
}
// cell average <-> cell-centre point value, limited, applied where the 5-point stencil exists
// out = limited conversion of in (all cells written; out and in must be different grids)
static void fv_cell_lim_to(Grid &out, const Grid &in, double sgn){
	const Grid& tmp = in;
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	double S[NUM_VARS]; grid_absmax(tmp, ghost_cells, nxCells+ghost_cells, ghost_cells, nyCells+ghost_cells, S);
	const double kappa = 0.5;
	auto pres = [&](const Grid& G, int ii, int jj)->double{   // pressure of a conservative state (ideal gas)
		double r = G(ii,jj,0); if(!(r > 0.0)) return -1.0;
		double ke = 0.5*(G(ii,jj,1)*G(ii,jj,1)+G(ii,jj,2)*G(ii,jj,2)+G(ii,jj,3)*G(ii,jj,3))/r;
		double mb = 0.5*(G(ii,jj,5)*G(ii,jj,5)+G(ii,jj,6)*G(ii,jj,6)+G(ii,jj,7)*G(ii,jj,7));
		return (gamma_plasma - 1.0)*(G(ii,jj,4) - ke - mb);
	};
#pragma omp parallel for schedule(static)
	for(int j = 0; j < NY; j++)
		for(int i = 0; i < NX; i++){
			if(i < 2 || i >= NX-2 || j < 2 || j >= NY-2){ for(int k = 0; k < NUM_VARS; k++) out(i,j,k) = tmp(i,j,k); continue; }   // outer rings unchanged
			double w = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(std::max(ghost_cells,std::min(nxCells+ghost_cells-1,i)), std::max(ghost_cells,std::min(nyCells+ghost_cells-1,j))))) : 1.0;
			for(int pass = 0; pass < 2; pass++){
				for(int k = 0; k < NUM_VARS; k++){
					double d2 = w * d2lim(tmp(i-2,j,k),tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k),tmp(i+2,j,k), kappa*dx*S[k]) + (1.0 - w) * d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k));
					d2 += w * d2lim(tmp(i,j-2,k),tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k),tmp(i,j+2,k), kappa*dy*S[k]) + (1.0 - w) * d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k));
					out(i,j,k) = tmp(i,j,k) + sgn * d2 / 24.0;
				}
				// [o3b] positivity guard (eps_M is blind to a pressure jump in a fluid at rest, e.g. the
				// blast wave at t=0): if the linear branch takes more than half of rho or p away, redo
				// this cell fully limited.  Only meaningful for conservative states (sgn = -1).
				if(sgn < 0.0){
					double r0 = tmp(i,j,0), p0 = pres(tmp,i,j), r1 = out(i,j,0), p1 = pres(out,i,j);
					if(r1 < 0.5*r0 || p1 < 0.5*p0){
						if(pass == 0 && w < 1.0){ w = 1.0; continue; }
						// [o3q] still inadmissible fully limited (low-beta cell: p is a small difference of
						// large deconvolved terms): keep the cell average as the point value there (2nd-order
						// locally, admissible) -- a negative point p otherwise reached the Riemann solver.
						for(int k = 0; k < NUM_VARS; k++) out(i,j,k) = tmp(i,j,k);
					}
				}
				break;
			}
		}
}
// the limited conversion of a scalar field, sgn > 0 (fv_cell_lim_to applied to one slot)
static void fv_cell_lim_s(Array &out, const Array &in, double sgn){
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	double S = 0.0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:S)
	for(int j = ghost_cells; j < nyCells+ghost_cells; j++) for(int i = ghost_cells; i < nxCells+ghost_cells; i++) S = fmax(S, fabs(in(i,j)));
	const double kappa = 0.5;
#pragma omp parallel for schedule(static)
	for(int j = 0; j < NY; j++)
		for(int i = 0; i < NX; i++){
			if(i < 2 || i >= NX-2 || j < 2 || j >= NY-2){ out(i,j) = in(i,j); continue; }
			double w = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(std::max(ghost_cells,std::min(nxCells+ghost_cells-1,i)), std::max(ghost_cells,std::min(nyCells+ghost_cells-1,j))))) : 1.0;
			double d2 = w * d2lim(in(i-2,j),in(i-1,j),in(i,j),in(i+1,j),in(i+2,j), kappa*dx*S) + (1.0 - w) * d2lin(in(i-1,j),in(i,j),in(i+1,j));
			d2 += w * d2lim(in(i,j-2),in(i,j-1),in(i,j),in(i,j+1),in(i,j+2), kappa*dy*S) + (1.0 - w) * d2lin(in(i,j-1),in(i,j),in(i,j+1));
			out(i,j) = in(i,j) + sgn * d2 / 24.0;
		}
}
static void fv_cell_lim(Grid &Q, double sgn){
	static Grid out(nxCells, nyCells);   // computed out of place; becomes Q's storage (no copies)
	fv_cell_lim_to(out, Q, sgn);
	std::swap(Q.data, out.data);
}
// face midpoint value <-> face average along the face (transverse direction), limited
// [o3e] FINITE-VOLUME 4th-order face operators on CELL AVERAGES (the face4/stagger4 family above is the
// finite-DIFFERENCE family: exact for point values, O(h^2) off when fed cell averages -- measured on the
// low-Mach vortex as a 2nd-order pressure error while u, B were 3rd order).
//   facefv : face value (transverse average inherited) from cell averages, (-1,7,7,-1)/12
//   gradfv : face value of the normal derivative from cell averages, (1,-15,15,-1)/(12h)
//   the divergence of face-averaged fluxes is the plain two-point difference (F_{i+1/2}-F_{i-1/2})/h
template<typename Fn>
inline double facefvx(Fn get, int i, int j){ return (-get(i-1,j) + 7.0*get(i,j) + 7.0*get(i+1,j) - get(i+2,j)) / 12.0; }
template<typename Fn>
inline double facefvy(Fn get, int i, int j){ return (-get(i,j-1) + 7.0*get(i,j) + 7.0*get(i,j+1) - get(i,j+2)) / 12.0; }
template<typename Fn>
inline double gradfvx(Fn get, int i, int j, double invh){ return (get(i-1,j) - 15.0*get(i,j) + 15.0*get(i+1,j) - get(i+2,j)) * invh / 12.0; }
template<typename Fn>
inline double gradfvy(Fn get, int i, int j, double invh){ return (get(i,j-1) - 15.0*get(i,j) + 15.0*get(i,j+1) - get(i,j+2)) * invh / 12.0; }

// [o3n] derivative of the face splitting weight along the face normal, 4th-order central on the face array.
// The implicit momentum update is conservative, m = m** - h D[(1-eps_f) p_f]; the face value of that update is
// I_f D[(1-eps_f) p_f] = (1-eps_f) d_n p - p_f d_n eps_f + O(h^4).  The projected face momentum of the energy
// equation must carry the second term: without it the implicit enthalpy flux is off by h (1-eps) H p_f d_n eps,
// an O(dt) (FIRST-ORDER-IN-TIME, h-independent) error wherever eps_M varies in space -- measured on the vortex
// with rho0 = 0.1 (M in the switch band 0.25..0.45): rho, p EOC 0.97, error halving with dt; zero for constant eps.
inline double deps_fx(const Array& e, int i, int j){ return (-e(i+2,j) + 8.0*e(i+1,j) - 8.0*e(i-1,j) + e(i-2,j)) / (12.0*dx); }
inline double deps_fy(const Array& e, int i, int j){ return (-e(i,j+2) + 8.0*e(i,j+1) - 8.0*e(i,j-1) + e(i,j-2)) / (12.0*dy); }
static void fv_transverse_lim(Grid &F, char dir, double sgn, int off = 0, Grid* keep = nullptr){   // off = 1: F(i,j) belongs to face i-1/2 (the right state w_xL / w_yL); keep: receives the unconverted F
	static Grid scratch(nxCells, nyCells);
	Grid& out = keep ? *keep : scratch;   // result, computed out of place; becomes F's storage at the end (no copies)
	const Grid& tmp = F;
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	double S[NUM_VARS]; grid_absmax(tmp, ghost_cells, nxCells+ghost_cells, ghost_cells, nyCells+ghost_cells, S);
	const double kappa = 0.5;
#pragma omp parallel for schedule(static)
	for(int j = 0; j < NY; j++)
		for(int i = 0; i < NX; i++){
			if(i < 2 || i >= NX-2 || j < 2 || j >= NY-2){ for(int k = 0; k < NUM_VARS; k++) out(i,j,k) = tmp(i,j,k); continue; }   // outer rings unchanged
			Array* ef = (dir == 'x') ? g_epsM_fx : g_epsM_fy;
			const int fi = (dir == 'x') ? i - off : i, fj = (dir == 'y') ? j - off : j;   // [o3q] eps of the face this state belongs to
			const double w = ef ? fmin(1.0, fmax(0.0, (*ef)(std::max(ghost_cells-1,std::min(nxCells+ghost_cells-1,fi)), std::max(ghost_cells-1,std::min(nyCells+ghost_cells-1,fj))))) : 1.0;
			for(int k = 0; k < NUM_VARS; k++){
				const double thr = kappa * ((dir == 'x') ? dy : dx) * S[k];
				double d2 = (dir == 'x')
				          ? w * d2lim(tmp(i,j-2,k),tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k),tmp(i,j+2,k), thr) + (1.0 - w) * d2lin(tmp(i,j-1,k),tmp(i,j,k),tmp(i,j+1,k))
				          : w * d2lim(tmp(i-2,j,k),tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k),tmp(i+2,j,k), thr) + (1.0 - w) * d2lin(tmp(i-1,j,k),tmp(i,j,k),tmp(i+1,j,k));
				out(i,j,k) = tmp(i,j,k) + sgn * d2 / 24.0;
			}
		}
	std::swap(F.data, out.data);
}
// out = (in, converted): out = in + sgn (d2x + d2y) in / 24 on 1..N-2, out = in on the outermost ring.
// (out and in must be different grids.)
// rho_guard: in the same pass, the density point-value guard of rho_pt_guard (out, in) on every cell.
long g_rho_guard = 0;
static void fv_cell_to(Grid &out, const Grid &in, double sgn, bool rho_guard = false){
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	long ng = 0;
#pragma omp parallel for schedule(static) reduction(+:ng)
	for(int j = 0; j < NY; j++){
		if(j == 0 || j == NY-1){ for(int i = 0; i < NX; i++) for(int k = 0; k < NUM_VARS; k++) out(i,j,k) = in(i,j,k); }
		else{
			for(int k = 0; k < NUM_VARS; k++){ out(0,j,k) = in(0,j,k); out(NX-1,j,k) = in(NX-1,j,k); }
			for(int i = 1; i < NX-1; i++)
				for(int k = 0; k < NUM_VARS; k++){
					double d2 = in(i+1,j,k) - 2.0*in(i,j,k) + in(i-1,j,k)
					          + in(i,j+1,k) - 2.0*in(i,j,k) + in(i,j-1,k);
					out(i,j,k) = in(i,j,k) + sgn * d2 / 24.0;
				}
		}
		if(rho_guard) for(int i = 0; i < NX; i++) if(!(out(i,j,0) >= 0.5*in(i,j,0))){ out(i,j,0) = in(i,j,0); ng++; }
	}
	g_rho_guard += ng;
}
template<int N> struct FieldN {
	double* data; int NXg;
	FieldN(){
		NXg = nxCells + 2*ghost_cells; const int NYg = nyCells + 2*ghost_cells;
		data = new double[(size_t)NXg*NYg*N];
#pragma omp parallel for schedule(static)
		for(int j = 0; j < NYg; j++) for(int q = 0; q < N*NXg; q++) data[(size_t)N*NXg*j + q] = 0.0;   // zero, first touch by the owning thread
	}
	~FieldN(){ delete[] data; }
	double  operator()(int i, int j, int c) const { return data[N*(i + j*NXg) + c]; }
	double& operator()(int i, int j, int c)       { return data[N*(i + j*NXg) + c]; }
};
using Field3 = FieldN<3>;
// fv_cell_to for N components; In(i,j,c) is any accessor (a FieldN, or N slots of a Grid); rho_guard as in fv_cell_to (component 0)
template<int N, typename In> static void fvN_to(FieldN<N> &out, const In &in, double sgn, bool rho_guard = false){
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	long ng = 0;
#pragma omp parallel for schedule(static) reduction(+:ng)
	for(int j = 0; j < NY; j++){
		if(j == 0 || j == NY-1){ for(int i = 0; i < NX; i++) for(int c = 0; c < N; c++) out(i,j,c) = in(i,j,c); }
		else{
			for(int c = 0; c < N; c++){ out(0,j,c) = in(0,j,c); out(NX-1,j,c) = in(NX-1,j,c); }
			for(int i = 1; i < NX-1; i++)
				for(int c = 0; c < N; c++){
					double d2 = in(i+1,j,c) - 2.0*in(i,j,c) + in(i-1,j,c)
					          + in(i,j+1,c) - 2.0*in(i,j,c) + in(i,j-1,c);
					out(i,j,c) = in(i,j,c) + sgn * d2 / 24.0;
				}
		}
		if(rho_guard) for(int i = 0; i < NX; i++) if(!(out(i,j,0) >= 0.5*in(i,j,0))){ out(i,j,0) = in(i,j,0); ng++; }
	}
	if(rho_guard) g_rho_guard += ng;
}
template<int N> static void fvN(FieldN<N> &F, double sgn){
	static FieldN<N> tmp;
	fvN_to<N>(tmp, F, sgn);
	std::swap(F.data, tmp.data);
}
struct GridSlots { const Grid& G; int k0; double operator()(int i, int j, int c) const { return G(i,j,k0 + c); } };   // slots k0.. of a Grid
// the same conversion for a scalar field
static void fv_cell_s(Array &out, const Array &in, double sgn){
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
#pragma omp parallel for schedule(static)
	for(int j = 0; j < NY; j++){
		if(j == 0 || j == NY-1){ for(int i = 0; i < NX; i++) out(i,j) = in(i,j); continue; }
		out(0,j) = in(0,j); out(NX-1,j) = in(NX-1,j);
		for(int i = 1; i < NX-1; i++){
			double d2 = in(i+1,j) - 2.0*in(i,j) + in(i-1,j)
			          + in(i,j+1) - 2.0*in(i,j) + in(i,j-1);
			out(i,j) = in(i,j) + sgn * d2 / 24.0;
		}
	}
}
static void fv_cell(Grid &Q, double sgn){
	static Grid tmp(nxCells, nyCells);
	fv_cell_to(tmp, Q, sgn);
	std::swap(Q.data, tmp.data);   // the converted field becomes Q's storage (no copy back)
}

// Transverse convolution / deconvolution of a face-indexed quantity.
//   sgn = +1 : point value  -> face average      (convolution)
//   sgn = -1 : face average -> point value       (deconvolution)
// dir = 'x' : the face normal is x, so the transverse direction is y.

void HLLEM_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,
                   Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y, double h){
	// Face fluxes are stored before they are differenced, for two reasons: each
	// face is then solved ONCE instead of twice, and the transverse quadrature
	// correction below needs the flux as a function along the face.
	static Grid Fx(nxCells, nyCells), Fy(nxCells, nyCells);
	// [o3l] EPSILON-CONSISTENT TRANSVERSE CONVOLUTION.  The face flux is eps-weighted PER FACE,
	// F = P^c + eps_M,f P^p + eps_B,f P^B + R.  Convolving the weighted sum along the face mixes the
	// weights of the neighbouring faces into the second difference: where eps switches between
	// adjacent faces (M_A ~ 1 ring of the low-Alfven vortex) that adds (Delta eps) F^B / 24 = O(1)
	// to the flux of the switching faces -- measured as an O(h) local truncation error of the momentum
	// equation (rho0 = 1e-3: mx 1.2e-3 / 6.4e-4 / 3.2e-4 / 1.4e-4 at 32..256, each share alone 3.4+).
	// The parts are convolved separately and recombined with the weights of THIS face.
	static Grid FxC(nxCells, nyCells), FxP(nxCells, nyCells), FxB(nxCells, nyCells), FxR(nxCells, nyCells);
	static Grid FyC(nxCells, nyCells), FyP(nxCells, nyCells), FyB(nxCells, nyCells), FyR(nxCells, nyCells);
	static bool parts_init = false;
	if(!parts_init){ parts_init = true; }   // (zero from the Grid constructor)

	static Array hchi(nxCells, nyCells);
	if(g_hllemflat && g_hllem_Qe){
		ShockDetector(hchi, *g_hllem_Qe);   // [fix C] chi of the explicit stage state (ghosts filled by bound_cond_array)
		if(g_hllemflat >= 2){   // + pressure-jump sensor: s_p = min(1, max_nbr |dp| / min(p) / 0.5), cell := max(chi, s_p)
			static Array pc(nxCells, nyCells), sp(nxCells, nyCells);
			const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++){
				double q[NUM_VARS], pr[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) q[k] = (*g_hllem_Qe)(i,j,k); ini_con2pri(pr, q); pc(i,j) = fmax(pr[4], 1e-300);
			}
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
				const double p0 = pc(i,j); double d = 0.0;
				const double pn[4] = {pc(i+1,j), pc(i-1,j), pc(i,j+1), pc(i,j-1)};
				for(int q = 0; q < 4; q++) d = fmax(d, fabs(pn[q] - p0)/fmin(pn[q], p0));
				sp(i,j) = fmin(1.0, d/0.5);
			}
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++) hchi(i,j) = fmax(hchi(i,j), sp(i,j));
		}
	}
	g_env_cap = a_max_x;   // [W]
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells - 2; j < nyCells + ghost_cells + 2; j++){        // [o3] x-faces on +-2 extra rows
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
			double QL[NUM_VARS], QR[NUM_VARS], f[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){ QL[k] = w_xR(i,j,k); QR[k] = w_xL(i+1,j,k); }
			if(g_hllemflat){ const int jc = std::max(1, std::min(nyCells + 2*ghost_cells - 2, j)); g_face_chi = fmax(hchi(i,jc), hchi(i+1,jc)); }
			HLLEM_x(QL, QR, eps_M_grid_x(i,j), eps_B_grid_x(i,j), f);
			{	// [o3l] parts of the numerical flux
				double cL[NUM_VARS], cR[NUM_VARS], pL[NUM_VARS], pR[NUM_VARS], bL[NUM_VARS], bR[NUM_VARS];
				f_convective(cL, QL); f_convective(cR, QR); f_pressure(pL, QL); f_pressure(pR, QR); f_magnetic(bL, QL); f_magnetic(bR, QR);
				const double wL = g_hll_wL, wR = g_hll_wR, eM = eps_M_grid_x(i,j), eB = eps_B_grid_x(i,j);
				for(int k = 0; k < NUM_VARS; k++){
					const double Pc = wL*cL[k] + wR*cR[k], Pp = wL*pL[k] + wR*pR[k], Pb = wL*bL[k] + wR*bR[k];
					FxC(i,j,k) = Pc; FxP(i,j,k) = Pp; FxB(i,j,k) = Pb; FxR(i,j,k) = f[k] - (Pc + eM*Pp + eB*Pb);
				}
			}   // (Fx is assembled from the convolved parts below)
		}
	}
	g_env_cap = a_max_y;   // [W]
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 2; i < nxCells + ghost_cells + 2; i++){    // [o3] y-faces on +-2 extra columns
			double QL[NUM_VARS], QR[NUM_VARS], f[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){ QL[k] = w_yR(i,j,k); QR[k] = w_yL(i,j+1,k); }
			if(g_hllemflat){ const int ic = std::max(1, std::min(nxCells + 2*ghost_cells - 2, i)); g_face_chi = fmax(hchi(ic,j), hchi(ic,j+1)); }
			HLLEM_y(QL, QR, eps_M_grid_y(i,j), eps_B_grid_y(i,j), f);
			{	// [o3l] parts of the numerical flux
				double cL[NUM_VARS], cR[NUM_VARS], pL[NUM_VARS], pR[NUM_VARS], bL[NUM_VARS], bR[NUM_VARS];
				g_convective(cL, QL); g_convective(cR, QR); g_pressure(pL, QL); g_pressure(pR, QR); g_magnetic(bL, QL); g_magnetic(bR, QR);
				const double wL = g_hll_wL, wR = g_hll_wR, eM = eps_M_grid_y(i,j), eB = eps_B_grid_y(i,j);
				for(int k = 0; k < NUM_VARS; k++){
					const double Pc = wL*cL[k] + wR*cR[k], Pp = wL*pL[k] + wR*pR[k], Pb = wL*bL[k] + wR*bR[k];
					FyC(i,j,k) = Pc; FyP(i,j,k) = Pp; FyB(i,j,k) = Pb; FyR(i,j,k) = f[k] - (Pc + eM*Pp + eB*Pb);
				}
			}   // (Fy is assembled from the convolved parts below)
		}
	}
	// [o3] midpoint flux -> face-averaged flux (limited transverse second difference); always on
	{	// [o3l] convolve the parts, recombine with the weights of each face
		fv_transverse_lim(FxC, 'x', +1.0); fv_transverse_lim(FxP, 'x', +1.0); fv_transverse_lim(FxB, 'x', +1.0); fv_transverse_lim(FxR, 'x', +1.0);
		fv_transverse_lim(FyC, 'y', +1.0); fv_transverse_lim(FyP, 'y', +1.0); fv_transverse_lim(FyB, 'y', +1.0); fv_transverse_lim(FyR, 'y', +1.0);
		const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
#pragma omp parallel for collapse(2) schedule(static)
		for(int j = 2; j < NY-2; j++) for(int i = 2; i < NX-2; i++){
			int ic = std::max(ghost_cells-1, std::min(nxCells+ghost_cells-1, i)), jc = std::max(ghost_cells-1, std::min(nyCells+ghost_cells-1, j));
			const double eMx = eps_M_grid_x(ic,jc), eBx = eps_B_grid_x(ic,jc), eMy = eps_M_grid_y(ic,jc), eBy = eps_B_grid_y(ic,jc);
			for(int k = 0; k < NUM_VARS; k++){
				Fx(i,j,k) = FxC(i,j,k) + eMx*FxP(i,j,k) + eBx*FxB(i,j,k) + FxR(i,j,k);
				Fy(i,j,k) = FyC(i,j,k) + eMy*FyP(i,j,k) + eBy*FyB(i,j,k) + FyR(i,j,k);
			}
		}
	}

	// ALL 8 COMPONENTS: the induction rows k=5,6,7 of the eps_B-weighted
	// explicit flux must be applied here; the implicit magnetic sub-system only
	// carries the (1-eps_B) share.
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int k = 0; k < NUM_VARS; k++){
				u_con(i,j,k) = Qs(i,j,k)
				             - h / dx * (Fx(i,j,k) - Fx(i-1,j,k))
				             - h / dy * (Fy(i,j,k) - Fy(i,j-1,k));
			}
		}
	}
}

		// ================================================================================
		// WENO3-Z + MP3 特征重构 —— 完整组件清单
		// ================================================================================
		
		// ---------- 1. 可靠的四次方程求根:伴随矩阵法,x/y 共用(定义见上文,HLLEM/compute_eigenvalues_x 已引用同一个)----------
		
		// ---------- 2. Family-ordered 特征值: [f-, A-, s-, div, ent, s+, A+, f+] ----------
		// (顺序必须和下面 compute_R_prim_x/y 的列索引对齐,不能用全局排序)
		void compute_eigenvalues_family_x(const double* U, double eM, double eB, Vec8 &lam){
			double rho = U[0], vx = U[1]/rho, vy = U[2]/rho, vz = U[3]/rho;
			double Bx = U[5], By = U[6], Bz = U[7];
			double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bp2 = By*By+Bz*Bz;
			double p = std::max((gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2), 1e-14);
			double a2 = gamma_plasma*p/rho, bx2 = Bx*Bx/rho, bp2 = Bp2/rho;
			
			double dA = std::max((1.0-eB)*(1.0-eB)*vx*vx+4.0*eB*eB*bx2, 0.0);
			double sA = std::sqrt(dA);
			double lamAm = 0.5*((1.0+eB)*vx-sA), lamAp = 0.5*((1.0+eB)*vx+sA);
			
			double s1 = (1+eM)*vx, a1 = (1+eB)*vx;
			double s0 = eM*vx*vx-eM*eM*a2, a0 = eB*vx*vx-eB*eB*bx2, q0 = eM*vx*vx;
			double eB2bp2 = eB*eB*bp2;
			double c3 = -(s1+a1), c2 = s1*a1+s0+a0-eB2bp2, c1 = -(s0*a1+s1*a0)+eB2bp2*s1, c0 = s0*a0-eB2bp2*q0;
			Vec4 ms = solve_quartic(c3,c2,c1,c0);   // 升序: [f-,s-,s+,f+]
			
			lam(0)=ms(0); lam(1)=lamAm; lam(2)=ms(1); lam(3)=0.0;
			lam(4)=vx;    lam(5)=ms(2); lam(6)=lamAp; lam(7)=ms(3);
		}
		void compute_eigenvalues_family_y(const double* U, double eM, double eB, Vec8 &lam){
			double rho = U[0], vx = U[1]/rho, vy = U[2]/rho, vz = U[3]/rho;
			double Bx = U[5], By = U[6], Bz = U[7];
			double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bt2 = Bx*Bx+Bz*Bz;
			double p = std::max((gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2), 1e-14);
			double a2 = gamma_plasma*p/rho, by2 = By*By/rho, bt2 = Bt2/rho;
			
			double dA = std::max((1.0-eB)*(1.0-eB)*vy*vy+4.0*eB*eB*by2, 0.0);
			double sA = std::sqrt(dA);
			double lamAm = 0.5*((1.0+eB)*vy-sA), lamAp = 0.5*((1.0+eB)*vy+sA);
			
			double s1 = (1+eM)*vy, a1 = (1+eB)*vy;
			double s0 = eM*vy*vy-eM*eM*a2, a0 = eB*vy*vy-eB*eB*by2, q0 = eM*vy*vy;
			double eB2bt2 = eB*eB*bt2;
			double c3 = -(s1+a1), c2 = s1*a1+s0+a0-eB2bt2, c1 = -(s0*a1+s1*a0)+eB2bt2*s1, c0 = s0*a0-eB2bt2*q0;
			Vec4 ms = solve_quartic(c3,c2,c1,c0);
			
			lam(0)=ms(0); lam(1)=lamAm; lam(2)=ms(1); lam(3)=0.0;
			lam(4)=vy;    lam(5)=ms(2); lam(6)=lamAp; lam(7)=ms(3);
		}
		
		// ---------- 3. Primitive 空间 right-eigenvector 矩阵(含分母正则化)----------
		inline void compute_R_prim_x(const double* U, double eM, double eB, const Vec8& lam, Mat8& Rp){
			const double rho = U[0], inv_rho = 1.0/rho;
			const double vx = U[1]*inv_rho, vy = U[2]*inv_rho, vz = U[3]*inv_rho;
			const double Bx = U[5], By = U[6], Bz = U[7];
			const double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bp2 = By*By+Bz*Bz;
			const double p = std::max((gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2), 1.0e-14);
			const double a2 = gamma_plasma*p*inv_rho, bx2 = Bx*Bx*inv_rho, bp2 = Bp2*inv_rho;
			const double a = std::sqrt(a2), b_total = std::sqrt(bx2+bp2);
			const double scale = std::fabs(vx)+a+b_total+1.0;
			const double reg = 1.0e-12*scale, reg_DA = reg*scale;
			const double eB2 = eB*eB, eB2_Bx = eB2*Bx;
			
			Rp.setZero();
			Rp(5,3) = 1.0;   // div
			Rp(0,4) = 1.0;   // ent
			{
				const bool Bperp_zero = (Bp2 < 1.0e-24);
				if(Bperp_zero){
					Rp(6,1)=-1.0; Rp(7,1)=0.0; Rp(6,6)=-1.0; Rp(7,6)=0.0;
				} else {
					for(int which = 0; which < 2; which++){
						int col = which==0 ? 1 : 6;
						double aA = vx - lam(col);
						if(std::fabs(aA) < reg) aA = (aA >= 0.0 ? reg : -reg);
						double coef = eB*Bx*inv_rho/aA;
						Rp(2,col) = -coef*Bz;
						Rp(3,col) =  coef*By;
						Rp(6,col) = -Bz;
						Rp(7,col) =  By;
					}
				}
			}
			const double eM_rho_a2 = eM*rho*a2, eB2_BxBy = eB2_Bx*By, eB2_BxBz = eB2_Bx*Bz;
			const double eB_By = eB*By, eB_Bz = eB*Bz;
			int mscols[4] = {0,2,5,7};
			for(int idx = 0; idx < 4; idx++){
				int col = mscols[idx];
				double la = lam(col);
				double al = vx-la, be = eM*vx-la, ph = eB*vx-la;
				if(std::fabs(al) < reg) al = (al >= 0.0 ? reg : -reg);
				if(std::fabs(be) < reg) be = (be >= 0.0 ? reg : -reg);
				double DA = al*ph - eB2*bx2;
				if(std::fabs(DA) < reg_DA) DA = (DA >= 0.0 ? reg_DA : -reg_DA);
				double inv_al=1.0/al, inv_be=1.0/be, inv_DA=1.0/DA, inv_rho_DA=inv_rho*inv_DA;
				Rp(0,col) = -rho*inv_al;
				Rp(1,col) =  1.0;
				Rp(2,col) = -eB2_BxBy*inv_rho_DA;
				Rp(3,col) = -eB2_BxBz*inv_rho_DA;
				Rp(4,col) = -eM_rho_a2*inv_be;
				Rp(6,col) = -eB_By*al*inv_DA;
				Rp(7,col) = -eB_Bz*al*inv_DA;
			}
		}
		inline void compute_R_prim_y(const double* U, double eM, double eB, const Vec8& lam, Mat8& Rp){
			const double rho = U[0], inv_rho = 1.0/rho;
			const double vx = U[1]*inv_rho, vy = U[2]*inv_rho, vz = U[3]*inv_rho;
			const double Bx = U[5], By = U[6], Bz = U[7];
			const double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz, Bt2 = Bx*Bx+Bz*Bz;
			const double p = std::max((gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2), 1.0e-14);
			const double a2 = gamma_plasma*p*inv_rho, by2 = By*By*inv_rho, bt2 = Bt2*inv_rho;
			const double a = std::sqrt(a2), b_total = std::sqrt(by2+bt2);
			const double scale = std::fabs(vy)+a+b_total+1.0;
			const double reg = 1.0e-12*scale, reg_DA = reg*scale;
			const double eB2 = eB*eB, eB2_By = eB2*By;
			
			Rp.setZero();
			Rp(6,3) = 1.0;   // div
			Rp(0,4) = 1.0;   // ent
			{
				const bool Btrans_zero = (Bt2 < 1.0e-24);
				if(Btrans_zero){
					Rp(5,1)=-1.0; Rp(7,1)=0.0; Rp(5,6)=-1.0; Rp(7,6)=0.0;
				} else {
					for(int which = 0; which < 2; which++){
						int col = which==0 ? 1 : 6;
						double aA = vy - lam(col);
						if(std::fabs(aA) < reg) aA = (aA >= 0.0 ? reg : -reg);
						double coef = eB*By*inv_rho/aA;
						Rp(1,col) = -coef*Bz;
						Rp(3,col) =  coef*Bx;
						Rp(5,col) = -Bz;
						Rp(7,col) =  Bx;
					}
				}
			}
			const double eM_rho_a2 = eM*rho*a2, eB2_ByBx = eB2_By*Bx, eB2_ByBz = eB2_By*Bz;
			const double eB_Bx = eB*Bx, eB_Bz = eB*Bz;
			int mscols[4] = {0,2,5,7};
			for(int idx = 0; idx < 4; idx++){
				int col = mscols[idx];
				double la = lam(col);
				double al = vy-la, be = eM*vy-la, ph = eB*vy-la;
				if(std::fabs(al) < reg) al = (al >= 0.0 ? reg : -reg);
				if(std::fabs(be) < reg) be = (be >= 0.0 ? reg : -reg);
				double DA = al*ph - eB2*by2;
				if(std::fabs(DA) < reg_DA) DA = (DA >= 0.0 ? reg_DA : -reg_DA);
				double inv_al=1.0/al, inv_be=1.0/be, inv_DA=1.0/DA, inv_rho_DA=inv_rho*inv_DA;
				Rp(0,col) = -rho*inv_al;
				Rp(1,col) = -eB2_ByBx*inv_rho_DA;
				Rp(2,col) =  1.0;
				Rp(3,col) = -eB2_ByBz*inv_rho_DA;
				Rp(4,col) = -eM_rho_a2*inv_be;
				Rp(5,col) = -eB_Bx*al*inv_DA;
				Rp(7,col) = -eB_Bz*al*inv_DA;
			}
		}
		
		// ---------- 4. L,R 组装 + 条件数触发降级(验证过的关键加固)----------
		inline void compute_LR_prim_x(const double *U, double eM, double eB, const Vec8 &lam, Mat8 &Lp, Mat8 &Rp){
			compute_R_prim_x(U, eM, eB, lam, Rp);
			Eigen::JacobiSVD<Mat8> svd(Rp);
			double smax = svd.singularValues()(0);
			double smin = svd.singularValues()(7);
			double cond = (smin > 0.0) ? smax/smin : 1e300;
			const double cond_threshold = 1.0e4;   // 按测试态估的,建议在你实际状态范围内重新筛一遍
			if(cond > cond_threshold || !std::isfinite(cond)){
				Rp = Mat8::Identity();
				Lp = Mat8::Identity();
			} else {
				Lp = Rp.inverse();
			}
		}
		inline void compute_LR_prim_y(const double *U, double eM, double eB, const Vec8 &lam, Mat8 &Lp, Mat8 &Rp){
			compute_R_prim_y(U, eM, eB, lam, Rp);
			Eigen::JacobiSVD<Mat8> svd(Rp);
			double smax = svd.singularValues()(0);
			double smin = svd.singularValues()(7);
			double cond = (smin > 0.0) ? smax/smin : 1e300;
			const double cond_threshold = 1.0e4;
			if(cond > cond_threshold || !std::isfinite(cond)){
				Rp = Mat8::Identity();
				Lp = Mat8::Identity();
			} else {
				Lp = Rp.inverse();
			}
		}
		
		// ================================================================================
		// BLOCK-STRUCTURED SPLIT CHARACTERISTIC BASIS  (replaces compute_LR_prim_x/y and
		// its cond>1e4 -> IDENTITY fallback, which silently disabled characteristic
		// reconstruction wherever the 8x8 R was ill-conditioned -- i.e. over the whole
		// Bperp~0 region, e.g. the entire blast-wave centre.)
		//
		// A~_p decomposes into FOUR blocks, three of them analytic:
		//   (1) Bx        1D, trivial: row 5 AND column 5 of A~_p are identically zero,
		//                 so dBx is its own characteristic variable (lambda = 0).
		//   (2) Alfven    2D, span{e_v,e_B} with the beta rotation. Invariant (verified
		//                 out-of-S 3.7e-14) and its dual functionals dv_perp, dB*_perp
		//                 annihilate every other wave (2.4e-16) -> analytic projection,
		//                 no inverse. Eigenpairs have NO denominators:
		//                   r_+- = (k, alpha_+-),  l_+- = (m, alpha_+-),
		//                   k = epsB*Bx/rho,  m = epsB*Bx,  alpha_+- = u - lamA_+-,
		//                   l.r = epsB^2*Bx^2/rho + alpha^2.
		//   (3) entropy   1D, ANALYTIC left eigenvector. Solving l.A~_p = u*l in the
		//                 aligned sector with l_rho = 1 gives
		//                   l_u    = -rho*u*(1-epsM)/(epsM^2 a^2)
		//                   l_p    = -1/(epsM a^2)
		//                   l_vpar = l_u * Bperp/Bx
		//                 Verified DIRECTLY: ||l.A~_p - u*l|| = 5.3e-13 over 40000 states.
		//                 NOTE this is exactly the amplitude that must NOT be used in
		//                 HLLEM, where the LD-only truncation leaves its 1/epsM part
		//                 uncancelled and it floods the density flux. HERE the
		//                 decomposition is COMPLETE (all 8 waves), so the magnetosonic
		//                 waves cancel that part exactly -- verified: all-8 sum == dq.
		//   (4) magnetosonic 4D remainder inside the 5D aligned sector
		//                 (rho,u,v_par,p,B_par), itself invariant under A~_p (leak 2.2e-16).
		//                 Solved by a 5x4 least-squares.
		//
		// WHY THIS MATTERS (all measured):
		//   cond(R_8x8)         = O(epsM^-4) AND blows up as Bperp->0:
		//                         884 -> 2.1e5 -> 4.4e8 -> 2.2e9 at Bperp = 1,1e-2,1e-4,1e-8
		//   cond(R_5x5 aligned) = O(epsM^-2), COMPLETELY INSENSITIVE to Bperp (66,66,66,66)
		//   cond(R_ms 5x4, entropy analytically peeled) -> 1:
		//                         epsM 1 / 0.1 / 0.01 / 0.005 -> 1.41 / 2.47 / 1.12 / 1.06
		//   The whole Bperp degeneracy lives in the Alfven block (now analytic); the whole
		//   epsM ceiling is the entropy normalisation (now analytic). What remains is a
		//   4-wave solve with O(1) conditioning for every epsM and Bperp.
		//   => no fallback, no ceiling; eps_M_min becomes a physics choice, not a numerics one.
		//
		// The 5x5 aligned eigenproblem uses EigenSolver rather than the closed form of
		// eq.(23): the analytic magnetosonic columns carry 1/alpha, 1/beta, 1/D_A, and at
		// Bperp->0 the quartic factorises so D_A -> 0 for two roots (0/0). EigenSolver
		// handles that degeneracy, and a 5x5 is far cheaper than the 8x8 it replaces.
		// Optimisation path if needed: closed-form columns plus a separate Bperp->0 branch
		// (the aligned block then splits 3+2, both analytic).
		// ================================================================================
		struct SplitChar8 {
			Vec8   R[8];        // right eigenvectors, PRIMITIVE space
			double lam[8];      // matching eigenvalues (diagnostics)
			double by, bz;      // beta rotation
			double kA, mA;      // Alfven: k = epsB*Bx/rho, m = epsB*Bx
			double alA[2];      // alpha_+- = vn - lamA_+-   [0]=minus [1]=plus
			double denA[2];     // l.r = m*k + alpha^2
			bool   decoupA;     // true when the normal field -> 0 and the 2x2 is diagonal
			bool   degenerate;  // build-time completeness self-check failed -> component-wise
			Eigen::Matrix<double,5,5> L5;   // inverse of the aligned 5x5 right-eigenvector matrix

			// slot layout: 0..4 aligned (4 magnetosonic + entropy), 5 Alfven-, 6 Alfven+, 7 div(Bx)
			void project(const Vec8& dq, Vec8& c) const {
				if(degenerate){ c = dq; return; }               // component-wise fallback
				c(7) = dq(5);                                   // Bx, lambda = 0
				double dvp = bz*dq(2) - by*dq(3);               // annihilating functionals
				double dBp = bz*dq(6) - by*dq(7);
				if(decoupA){                                   // Bn -> 0: 2x2 is diagonal
					c(5) = dBp;                                // pure B*_perp mode, lambda = epsB*vn
					c(6) = dvp;                                // pure v_perp   mode, lambda = vn
				}else{
					c(5) = (mA*dvp + alA[0]*dBp)/denA[0];
					c(6) = (mA*dvp + alA[1]*dBp)/denA[1];
				}
				Eigen::Matrix<double,5,1> d5;
				d5(0) = dq(0);
				d5(1) = dq(1);
				d5(2) = by*dq(2) + bz*dq(3);
				d5(3) = dq(4);
				d5(4) = by*dq(6) + bz*dq(7);
				Eigen::Matrix<double,5,1> c5 = L5 * d5;
				for(int m = 0; m < 5; m++) c(m) = c5(m);
			}
			void reconstruct(const Vec8& c, Vec8& dq) const {
				dq.setZero();   // (degenerate: R = identity in the build frame, see the fallback)
				for(int m = 0; m < 8; m++) dq += c(m) * R[m];
			}
		};

		// embed an aligned-sector 5-vector (rho,u,v_par,p,B_par) into primitive 8-space
		static inline void embed_aligned(const Eigen::Matrix<double,5,1>& r5,
		                                 double by, double bz, Vec8& r8){
			r8.setZero();
			r8(0) = r5(0);
			r8(1) = r5(1);
			r8(2) = by*r5(2);  r8(3) = bz*r5(2);
			r8(4) = r5(3);
			r8(6) = by*r5(4);  r8(7) = bz*r5(4);
		}

		void build_split_char_x(const double* U, double eM_phys, double eB, SplitChar8& S){
			const double eM = std::max(eM_phys, g_eps_basis_min);   // fix (3): basis floor
			const double rho = U[0], inv_rho = 1.0/rho;
			const double vx = U[1]*inv_rho, vy = U[2]*inv_rho, vz = U[3]*inv_rho;
			const double Bx = U[5], By = U[6], Bz = U[7];
			const double q2 = vx*vx+vy*vy+vz*vz, B2 = Bx*Bx+By*By+Bz*Bz;
			const double p  = std::max((gamma_plasma-1.0)*(U[4]-0.5*rho*q2-0.5*B2), 1.0e-14);
			const double a2 = gamma_plasma*p*inv_rho, bx2 = Bx*Bx*inv_rho;
			const double Bperp = std::sqrt(By*By + Bz*Bz);
			const double spd = std::fabs(vx) + std::sqrt(a2) + std::sqrt(B2*inv_rho) + 1.0;

			if(Bperp > 1.0e-12*std::sqrt(B2 + 1.0)){ S.by = By/Bperp; S.bz = Bz/Bperp; }
			else { S.by = 1.0; S.bz = 0.0; }

			// ---------- (1) divergence mode, lambda = 0 ----------
			S.lam[7] = 0.0;
			S.R[7].setZero(); S.R[7](5) = 1.0;

			// ---------- (2) Alfven, analytic ----------
			// Generic case: coupled 2x2, eigenpairs with no denominators. As the NORMAL
			// field vanishes the pair degenerates -- alpha_A+ = vn - lamA+ -> 0 while
			// k,m -> 0 as well, so r_A+ -> 0 and the basis would lose a direction. In that
			// limit the 2x2 becomes DIAGONAL, [[vn,0],[0,epsB*vn]]: v_perp and B*_perp
			// decouple exactly, with eigenvalues vn and epsB*vn. Use that analytically.
			// (This is the same Alfven-entropy coalescence at Bn=0 as in full MHD, and it is
			//  not a corner case here: the blast wave has By == 0 everywhere at t=0, so the
			//  Y-direction reconstruction sits exactly on it.)
			const double ca_eff = eB*std::fabs(Bx)/std::sqrt(rho);
			S.decoupA = !(ca_eff > 1.0e-8*spd);
			if(!S.decoupA){
				double lamAm, lamAp;
				SplitAlfvenEigenvalues(vx, eB, bx2, lamAm, lamAp);
				S.lam[5] = lamAm;  S.lam[6] = lamAp;
				S.kA = eB*Bx*inv_rho;  S.mA = eB*Bx;
				S.alA[0] = vx - lamAm; S.alA[1] = vx - lamAp;
				for(int w = 0; w < 2; w++){
					S.denA[w] = S.mA*S.kA + S.alA[w]*S.alA[w];
					Vec8 &r = S.R[5+w]; r.setZero();
					r(2) =  S.kA*S.bz;      r(3) = -S.kA*S.by;
					r(6) =  S.alA[w]*S.bz;  r(7) = -S.alA[w]*S.by;
				}
			}else{
				S.kA = 0.0; S.mA = 0.0;
				S.lam[5] = eB*vx;  S.alA[0] = vx - S.lam[5];  S.denA[0] = 1.0;
				S.lam[6] = vx;     S.alA[1] = 0.0;            S.denA[1] = 1.0;
				S.R[5].setZero(); S.R[5](6) = S.bz;  S.R[5](7) = -S.by;   // pure B*_perp
				S.R[6].setZero(); S.R[6](2) = S.bz;  S.R[6](3) = -S.by;   // pure v_perp
			}

			// ---------- (3)+(4) aligned 5D sector: entropy + 4 magnetosonic together ----------
			Eigen::Matrix<double,5,5> A5; A5.setZero();
			A5(0,0)=vx; A5(0,1)=rho;
			A5(1,1)=vx; A5(1,3)=eM*inv_rho;      A5(1,4)=eB*Bperp*inv_rho;
			A5(2,2)=vx;                          A5(2,4)=-eB*Bx*inv_rho;
			A5(3,1)=eM*rho*a2; A5(3,3)=eM*vx;
			A5(4,1)=eB*Bperp;  A5(4,2)=-eB*Bx;   A5(4,4)=eB*vx;
			Eigen::EigenSolver<Eigen::Matrix<double,5,5>> es(A5);
			const Eigen::Matrix<std::complex<double>,5,5> EV = es.eigenvectors();   // once: every call rebuilds and normalises the complex matrix
			Eigen::Matrix<double,5,5> R5;
			for(int i = 0; i < 5; i++){
				Eigen::Matrix<double,5,1> r5;
				for(int k = 0; k < 5; k++) r5(k) = EV(k,i).real();
				double nn = r5.norm(); if(nn > 0.0) r5 /= nn;   // unit columns: best conditioning
				R5.col(i) = r5;
				S.lam[i] = es.eigenvalues()(i).real();
				embed_aligned(r5, S.by, S.bz, S.R[i]);
			}
			S.L5 = R5.inverse();

			// ---------- build-time completeness self-check ----------
			// The one degeneracy this construction does not resolve analytically is
			// NORMAL field -> 0 with a LARGE transverse field: the aligned block then
			// carries lambda = vn as a REPEATED eigenvalue (the decoupled v_par mode and
			// the entropy mode coincide) and EigenSolver returns a nearly parallel pair.
			// Verified reach of the construction (max rel error of
			// reconstruct(project(dq)) over 600 random states each):
			//     Bn = O(1)          : 6e-14 (epsM=1) ... 3e-12 (epsM=0.01)   <- exact
			//     Bn = O(1), Bperp->0: 3e-15 ... 3e-12                        <- exact
            //     Bn ~ 1e-6          : 5e-10 ... 2e-6                         <- degrades
            // so the fallback below is needed only when the normal field is ~6 orders
            // below the transverse one. Contrast the OLD cond>1e4 identity fallback,
            // which fired over the ENTIRE Bperp~0 region (all of the blast-wave centre)
            // and disabled characteristic reconstruction there completely.
			S.degenerate = false;
			{
				Vec8 probe; probe << 1.0,-0.7,0.4,-0.3,0.9,0.6,-0.5,0.2;
				probe /= probe.norm();
				Vec8 cc, bb;
				S.project(probe, cc);
				S.reconstruct(cc, bb);
				double err = (bb - probe).norm();
				// g_weno3_no_fallback: NEVER drop to component-wise on the completeness
				// tolerance -- only a genuinely non-finite basis is guarded. The single
				// unresolved degeneracy (normal field ->0, transverse field large) has a
				// completeness error <= ~2e-6, orders below the WENO truncation error, so
				// the old 1e-6 trip was over-conservative and threw away a valid basis.
				bool bad = !std::isfinite(err)
				        || (err > g_char_fallback_tol);
				if(bad){ S.degenerate = true;
					// [o3q] component-wise fallback: identity basis.  R must be set (not bypassed): the y-sweep
					// rotates R back to the (x,y,z) frame; bypassing it returned the ROTATED increments
					// (u got v's increment, Bx got By's, ...).
					for(int m = 0; m < 8; m++){ S.R[m].setZero(); S.R[m](m) = 1.0; }
				}
			}
		}

		// y-direction: cyclic permutation (x,y,z)->(y,z,x), same convention as HLLEM_y,
		// so the second set of signs is never hand-derived.
		static inline void rot_fwd(const double* q, double* r){
			r[0]=q[0]; r[1]=q[2]; r[2]=q[3]; r[3]=q[1]; r[4]=q[4]; r[5]=q[6]; r[6]=q[7]; r[7]=q[5];
		}
		static inline void rot_fwd8(const Vec8& q, Vec8& r){
			r(0)=q(0); r(1)=q(2); r(2)=q(3); r(3)=q(1); r(4)=q(4); r(5)=q(6); r(6)=q(7); r(7)=q(5);
		}
		static inline void rot_bwd8(const Vec8& r, Vec8& q){
			q(0)=r(0); q(1)=r(3); q(2)=r(1); q(3)=r(2); q(4)=r(4); q(5)=r(7); q(6)=r(5); q(7)=r(6);
		}
		void build_split_char_y(const double* U, double eM, double eB, SplitChar8& S){
			double Ur[NUM_VARS]; rot_fwd(U, Ur);
			build_split_char_x(Ur, eM, eB, S);
			for(int m = 0; m < 8; m++){ Vec8 t = S.R[m]; rot_bwd8(t, S.R[m]); }
		}

		// ================================================================================
		// COMPLETE (all-8-wave) split anti-diffusion for HLLEM. Forward-declared above
		// HLLEM_x. Uses the SAME blocked split basis as the reconstruction, so the
		// Riemann solver and WENO3 share one characteristic decomposition -- which is
		// what makes the whole scheme "split basis everywhere" with no eps=1 anywhere.
		//
		// Why ALL waves and not just the LD ones: the exact split entropy projection
		// carries a 1/epsM part that is cancelled ONLY by the magnetosonic waves (the
		// complete 8-wave sum reproduces dq identically). Truncating to LD leaves that
		// part uncancelled and it enters the density flux directly through R_ent[0]=1.
		// With the complete sum there is nothing to cancel and nothing to approximate.
		//
		// delta_k is the standard HLLEM AP coefficient. It vanishes identically for any
		// wave sitting on the envelope (lambda_k = s_L => lambda^- = s_L, lambda^+ = 0
		// => delta = 1 - 1 - 0 = 0), so the fast waves are auto-excluded and the HLL
		// robustness on the genuinely nonlinear extremes is retained.
		// ================================================================================

		static inline void project_y(const SplitChar8& S, const Vec8& dq, Vec8& c){
			Vec8 dqr; rot_fwd8(dq, dqr);
			S.project(dqr, c);
		}
		// ---------- 5. MP3 limiter ----------
		inline double minmod2(double a, double b){
			if(a*b <= 0.0) return 0.0;
			return (std::fabs(a) < std::fabs(b)) ? a : b;
		}
		inline double mp3_limiter(double qm1, double q0, double qp1, double qp2, double q_rec, bool from_left){
			constexpr double alpha = 2.0;
			constexpr double tol = 0;
			double c, cp, cm, cf;
			if(from_left){ c=q0; cp=qp1; cm=qm1; cf=qp2; }
			else{ c=qp1; cp=q0; cm=qp2; cf=qm1; }
			
			double q_MP = c + minmod2(cp-c, alpha*(c-cm));
			if((q_rec-c)*(q_rec-q_MP) <= tol) return q_rec;
			
			double d_c  = cm - 2.0*c  + cp;
			double d_cp = c  - 2.0*cp + cf;
			double d_M  = minmod2(d_c, d_cp);
			
			double q_UL = c + alpha*(c-cm);
			double q_AV = 0.5*(c+cp);
			double q_MD = q_AV - 0.5*d_M;
			double q_LC = c + 0.5*(c-cm) + (4.0/3.0)*d_M;
			
			double q_min = std::max(std::min({c,cp,q_MD}), std::min({c,q_UL,q_LC}));
			double q_max = std::min(std::max({c,cp,q_MD}), std::max({c,q_UL,q_LC}));
			
			return std::min(q_max, std::max(q_min, q_rec));
		}
		
		// ---------- 6. Positivity theta-scaling ----------
		inline double compute_theta_one_var(double q_c, double q_r, double floor){
			if(q_r >= floor) return 1.0;
			if(q_c <= floor) return 0.0;
			double theta = (floor-q_c)/(q_r-q_c);
			if(theta < 0.0) return 0.0;
			if(theta > 1.0) return 1.0;
			return theta;
		}
		inline double compute_theta(const double* u_c, const double* q_r){
			int RHO_IDX=0, P_IDX=4;
			double RHO_FLOOR=1e-13, P_FLOOR=1e-13;
			double theta_rho = compute_theta_one_var(u_c[RHO_IDX], q_r[RHO_IDX], RHO_FLOOR);
			double theta_p   = compute_theta_one_var(u_c[P_IDX],   q_r[P_IDX],   P_FLOOR);
			return std::min(theta_rho, theta_p);
		}
		inline double positivity_theta_scale(const double* u_c, double* q_r){
			double theta = compute_theta(u_c, q_r);
			if(theta < 1.0){
				for(int m = 0; m < NUM_VARS; m++) q_r[m] = u_c[m] + theta*(q_r[m]-u_c[m]);
			}
			return theta;
		}
		inline double positivity_theta_scale_vec(const Grid& u_pri, int i, int j, Vec8& q_recon){
			double u_c[NUM_VARS], q_r[NUM_VARS];
			for(int m = 0; m < NUM_VARS; m++){ u_c[m]=u_pri(i,j,m); q_r[m]=q_recon(m); }
			double theta = positivity_theta_scale(u_c, q_r);
			for(int m = 0; m < NUM_VARS; m++) q_recon(m) = q_r[m];
			return theta;
		}
		
		// ---------- 7. epsilon face->cell 插值 ----------
		inline double face2cell_x(const Array &ef, int i, int j){
			return (-ef(i-2,j) + 9.0*ef(i-1,j) + 9.0*ef(i,j) - ef(i+1,j)) / 16.0;
		}
		inline double face2cell_y(const Array &ef, int i, int j){
			return (-ef(i,j-2) + 9.0*ef(i,j-1) + 9.0*ef(i,j) - ef(i,j+1)) / 16.0;
		}
		
		// ---------- 8. 主驱动:WENO3Reconstruction ----------
		void WENO3Reconstruction(Grid &u_con,
			Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,
			Array &eps_M_x, Array &eps_M_y,
			Array &eps_B_x, Array &eps_B_y)
		{
			SpaceOrder = 3;
			static Grid u_pri(nxCells, nyCells);
			{	// [o3] cell average -> point value -> primitive -> cell average (limited).
				// W(Qbar) differs from the cell average of W at O(h^2).
				static Grid Qpt(nxCells, nyCells);
				fv_cell_lim_to(Qpt, u_con, -1.0);
				ucon2upri(u_pri, Qpt);
				fv_cell_lim(u_pri, +1.0);
				bound_cond_u(u_pri);
			}
			// [o3] WENO-Z epsilon = (h * S_m)^2 with S_m the scale of the m-th CHARACTERISTIC
			// variable (from the domain-wide scales of the primitives, mapped through the local
			// left eigenvectors).  Measured on the smooth 2D wave: eps=1e-40 gives L1 EOC 2.0
			// (the 3-point Z weights are O(1) off ideal in the h-wide bands around critical
			// points of each characteristic variable); eps=(hS)^2 keeps omega = d + O(h^2)
			// uniformly (tau/eps = O(h)) while a jump (tau ~ S^2 >> h^2 S^2) still gets ENO
			// weights.  Scale-invariant: eps scales with S^2.
			double S_scale[NUM_VARS];
			grid_absmax(u_pri, ghost_cells, nxCells + ghost_cells, ghost_cells, nyCells + ghost_cells, S_scale);
			auto eps_char = [&](SplitChar8& SCb, double hh, bool ydir, Vec8& eps_c){
				eps_c.setZero();
				for(int nn = 0; nn < NUM_VARS; nn++){
					Vec8 e; e.setZero(); e(nn) = S_scale[nn]; Vec8 ce;
					if(ydir){ Vec8 er; rot_fwd8(e, er); SCb.project(er, ce); }   // the y-basis expects the (x,y,z)->(y,z,x) rotated frame
					else SCb.project(e, ce);                                    // (bug: unrotated, the transverse-velocity channel got S_w = 0 -> eps ~ 0)
					for(int m = 0; m < NUM_VARS; m++) eps_c(m) += ce(m)*ce(m);
				}
				for(int m = 0; m < NUM_VARS; m++) eps_c(m) = hh*hh*eps_c(m) + 1e-300;
			};
			
			// ==================== X-direction ====================
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = ghost_cells - 2; j < nyCells + ghost_cells + 2; j++){   // [o3] +-2 rows: transverse stencils
				for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
					double Uavg[NUM_VARS];
					for(int m = 0; m < NUM_VARS; m++)
						Uavg[m] = 0.5*(u_con(i,j,m) + u_con(i+1,j,m));
					
					double eM = eps_M_x(i,j);
					double eB = eps_B_x(i,j);
					
					// block-structured split characteristic basis: no 8x8 inverse,
					// no cond>1e4 identity fallback (see SplitChar8 above)
					SplitChar8 SC;
					build_split_char_x(Uavg, eM, eB, SC);
					
					Vec8 dq_L, dq_C, dq_R;
					for(int m = 0; m < NUM_VARS; m++){
						dq_L(m) = u_pri(i,  j,m) - u_pri(i-1,j,m);
						dq_C(m) = u_pri(i+1,j,m) - u_pri(i,  j,m);
						dq_R(m) = u_pri(i+2,j,m) - u_pri(i+1,j,m);
					}
					Vec8 dc_L, dc_C, dc_R;
					SC.project(dq_L, dc_L); SC.project(dq_C, dc_C); SC.project(dq_R, dc_R);
					
					Vec8 incr_left, incr_right, eps_c; eps_char(SC, dx, false, eps_c);
					for(int m = 0; m < 8; m++){
						const double eps_w = eps_c(m);
						double b0L=dc_L(m)*dc_L(m), b1L=dc_C(m)*dc_C(m);
						double tauL=std::fabs(b0L-b1L), r;
						r = tauL/(b0L+eps_w); double a0L=(1.0/3.0)*(1.0+r*r);
						r = tauL/(b1L+eps_w); double a1L=(2.0/3.0)*(1.0+r*r);
						double sumL=a0L+a1L;
						incr_left(m) = 0.5*(a0L*dc_L(m)+a1L*dc_C(m))/sumL;
						
						double b0R=dc_C(m)*dc_C(m), b1R=dc_R(m)*dc_R(m);
						double tauR=std::fabs(b0R-b1R);
						r = tauR/(b0R+eps_w); double a0R=(2.0/3.0)*(1.0+r*r);
						r = tauR/(b1R+eps_w); double a1R=(1.0/3.0)*(1.0+r*r);
						double sumR=a0R+a1R;
						incr_right(m) = -0.5*(a0R*dc_C(m)+a1R*dc_R(m))/sumR;
					}
					
					Vec8 dq_L_recon, dq_R_recon;
					SC.reconstruct(incr_left,  dq_L_recon);
					SC.reconstruct(incr_right, dq_R_recon);
					
					Vec8 qL, qR;
					for(int m = 0; m < NUM_VARS; m++){
						qL(m) = u_pri(i,  j,m) + dq_L_recon(m);
						qR(m) = u_pri(i+1,j,m) + dq_R_recon(m);
					}
					
					for(int m = 0; m < NUM_VARS; m++){
						qL(m) = mp3_limiter(u_pri(i-1,j,m), u_pri(i,  j,m), u_pri(i+1,j,m), u_pri(i+2,j,m), qL(m), true);
						qR(m) = mp3_limiter(u_pri(i-1,j,m), u_pri(i,  j,m), u_pri(i+1,j,m), u_pri(i+2,j,m), qR(m), false);
					}

					positivity_theta_scale_vec(u_pri, i,   j, qL);
					positivity_theta_scale_vec(u_pri, i+1, j, qR);
					
					for(int m = 0; m < NUM_VARS; m++){
						w_xR(i,  j,m) = qL(m);
						w_xL(i+1,j,m) = qR(m);
					}
				}
			}
			
			// ==================== Y-direction ====================
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
				for(int i = ghost_cells - 2; i < nxCells + ghost_cells + 2; i++){   // [o3] +-2 columns
					double Uavg[NUM_VARS];
					for(int m = 0; m < NUM_VARS; m++)
						Uavg[m] = 0.5*(u_con(i,j,m) + u_con(i,j+1,m));
					
					double eM = eps_M_y(i,j);
					double eB = eps_B_y(i,j);
					
					SplitChar8 SC;
					build_split_char_y(Uavg, eM, eB, SC);
					
					Vec8 dq_L, dq_C, dq_R;
					for(int m = 0; m < NUM_VARS; m++){
						dq_L(m) = u_pri(i,j,  m) - u_pri(i,j-1,m);
						dq_C(m) = u_pri(i,j+1,m) - u_pri(i,j,  m);
						dq_R(m) = u_pri(i,j+2,m) - u_pri(i,j+1,m);
					}
					Vec8 dc_L, dc_C, dc_R;
					project_y(SC, dq_L, dc_L); project_y(SC, dq_C, dc_C); project_y(SC, dq_R, dc_R);
					
					Vec8 incr_left, incr_right, eps_c; eps_char(SC, dy, true, eps_c);
					for(int m = 0; m < 8; m++){
						const double eps_w = eps_c(m);
						double b0L=dc_L(m)*dc_L(m), b1L=dc_C(m)*dc_C(m);
						double tauL=std::fabs(b0L-b1L), r;
						r = tauL/(b0L+eps_w); double a0L=(1.0/3.0)*(1.0+r*r);
						r = tauL/(b1L+eps_w); double a1L=(2.0/3.0)*(1.0+r*r);
						double sumL=a0L+a1L;
						incr_left(m) = 0.5*(a0L*dc_L(m)+a1L*dc_C(m))/sumL;
						
						double b0R=dc_C(m)*dc_C(m), b1R=dc_R(m)*dc_R(m);
						double tauR=std::fabs(b0R-b1R);
						r = tauR/(b0R+eps_w); double a0R=(2.0/3.0)*(1.0+r*r);
						r = tauR/(b1R+eps_w); double a1R=(1.0/3.0)*(1.0+r*r);
						double sumR=a0R+a1R;
						incr_right(m) = -0.5*(a0R*dc_C(m)+a1R*dc_R(m))/sumR;
					}
					
					Vec8 dq_L_recon, dq_R_recon;
					SC.reconstruct(incr_left,  dq_L_recon);
					SC.reconstruct(incr_right, dq_R_recon);
					
					Vec8 qL, qR;
					for(int m = 0; m < NUM_VARS; m++){
						qL(m) = u_pri(i,j,  m) + dq_L_recon(m);
						qR(m) = u_pri(i,j+1,m) + dq_R_recon(m);
					}
					
					for(int m = 0; m < NUM_VARS; m++){
						qL(m) = mp3_limiter(u_pri(i,j-1,m), u_pri(i,j,  m), u_pri(i,j+1,m), u_pri(i,j+2,m), qL(m), true);
						qR(m) = mp3_limiter(u_pri(i,j-1,m), u_pri(i,j,  m), u_pri(i,j+1,m), u_pri(i,j+2,m), qR(m), false);
					}

					positivity_theta_scale_vec(u_pri, i, j,   qL);
					positivity_theta_scale_vec(u_pri, i, j+1, qR);
					
					for(int m = 0; m < NUM_VARS; m++){
						w_yR(i,j,  m) = qL(m);
						w_yL(i,j+1,m) = qR(m);
					}
				}
			}
			
			{	// [o3] The 1D reconstruction returns (x-point, y-average) face values of the
				// primitives; the Riemann solver needs the MIDPOINT value.  Deconvolve
				// transversally (limited), floor-check, convert pointwise -- and hand the
				// POINT states to HLLEM.  (The old flag path convolved the states back to
				// face averages AND convolved the flux again: double counting, O(h^2).)
				static Grid sx_L(nxCells, nyCells), sx_R(nxCells, nyCells), sy_L(nxCells, nyCells), sy_R(nxCells, nyCells);   // face averages before deconvolution
				fv_transverse_lim(w_xL, 'x', -1.0, 1, &sx_L);  fv_transverse_lim(w_xR, 'x', -1.0, 0, &sx_R);   // sx_*, sy_*: the unconverted states
				fv_transverse_lim(w_yL, 'y', -1.0, 1, &sy_L);  fv_transverse_lim(w_yR, 'y', -1.0, 0, &sy_R);
				const int NXc = nxCells + 2*ghost_cells, NYc = nyCells + 2*ghost_cells;
				Grid* ws[4] = {&w_xL, &w_xR, &w_yL, &w_yR}; Grid* sv[4] = {&sx_L, &sx_R, &sy_L, &sy_R};
				for(int q = 0; q < 4; q++)
#pragma omp parallel for schedule(static)
					for(int jj = 0; jj < NYc; jj++) for(int ii = 0; ii < NXc; ii++){
						Grid& W = *ws[q]; Grid& V = *sv[q];
						// [o3b] positivity guard: the deconvolution may not take more than half of rho or p
						// (blast wave at rest: eps_M = 0 makes the branch linear across a 1000:0.1 jump)
						if(!(W(ii,jj,0) > 0.5*V(ii,jj,0)) || !(W(ii,jj,4) > 0.5*V(ii,jj,4)) || !(W(ii,jj,0) > 1e-13) || !(W(ii,jj,4) > 1e-13)){
							for(int m = 0; m < NUM_VARS; m++) W(ii,jj,m) = V(ii,jj,m);   // keep the (MP3-limited) face average
						}
						// primitive -> conservative in the same pass (was upri2ucon of the whole grid)
						double pr[NUM_VARS], co[NUM_VARS];
						for(int m = 0; m < NUM_VARS; m++) pr[m] = W(ii,jj,m);
						ini_pri2con(co, pr);
						for(int m = 0; m < NUM_VARS; m++) W(ii,jj,m) = co[m];
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
	double e = pri[4]/((gamma_plasma - 1) * pri[0]);
	return e;
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
	// [o3] point-value initial condition (primitive) ...
	auto ic_point = [&](double x, double y, double* q){
		for(int k = 0; k < NUM_VARS; k++) q[k] = 0.0;
		switch (test) {
			case test0:
				if(x < 0){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 0.125;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}
				break;
			case test1:
				if(x < 0){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 0.75;
					q[6] = 1;
					q[7] = 0;
				}else{
					q[0] = 0.125;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 0.75;
					q[6] = -1;
					q[7] = 0;
				}
				break;
			case test1y:
				if(y < 0){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 1;
					q[6] = 0.75;
					q[7] = 0;
				}else{
					q[0] = 0.125;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = -1;
					q[6] = 0.75;
					q[7] = 0;
				}
				break;
			case test2:
				if(x < -0.1){
					q[0] = 1.08;
					q[1] = 1.2;
					q[2] = 0.01;
					q[3] = 0.5;
					q[4] = 0.95;
					q[5] = 0.564190;
					q[6] = 1.015541;
					q[7] = 0.564190;
				}else{
					q[0] = 0.9891;
					q[1] = -0.0131;
					q[2] = 0.0269;
					q[3] = 0.010037;
					q[4] =  0.97159;
					q[5] = 0.564190;
					q[6] = 1.135262;
					q[7] = 0.564923;
				}
				break;
			case test3:
				if(x < -0.1){
					q[0] = 1.7;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1.7;
					q[5] = 1.1;
					q[6] = 1;
					q[7] = 0;
				}else{
					q[0] = 0.2;
					q[1] = 0;
					q[2] = 0;
					q[3] = -1.496891;
					q[4] =  0.2;
					q[5] = 1.1;
					q[6] = 0.785887;
					q[7] = 0.618370;
				}
				break;
			case test4:
				if(x < 0){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 1.3;
					q[6] = 1;
					q[7] = 0;
				}else{
					q[0] = 0.4;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] =  0.4;
					q[5] = 1.3;
					q[6] = -1;
					q[7] = 0;
				}
				break;
			case test5:
				if(x < 0){
					q[0] = 0.15;
					q[1] = 21.55;
					q[2] = 1.0;
					q[3] = 1.0;
					q[4] = 0.28;
					q[5] = 0.05 / sqrt(4 * M_PI);
					q[6] = -2.0 / sqrt(4 * M_PI);
					q[7] = -1.0 / sqrt(4 * M_PI);
				}else{
					q[0] = 0.1;
					q[1] = -26.45;
					q[2] = 0;
					q[3] = 0;
					q[4] =  0.1;
					q[5] = 0.05 / sqrt(4 * M_PI);
					q[6] = 2.0 / sqrt(4 * M_PI);
					q[7] = 1.0 / sqrt(4 * M_PI);
				}
				break;
			case test6:
				if(x < 0){
					q[0] = 1.0;
					q[1] = 36.87;
					q[2] = -0.115;
					q[3] = -0.0386;
					q[4] = 1.0;
					q[5] = 4.0 / sqrt(4 * M_PI);
					q[6] = 4.0 / sqrt(4 * M_PI);
					q[7] = 1.0 / sqrt(4 * M_PI);
				}else{
					q[0] = 1;
					q[1] = -36.87;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 4 / sqrt(4 * M_PI);
					q[6] = 4 / sqrt(4 * M_PI);
					q[7] = 1 / sqrt(4 * M_PI);
				}
				break;
			case test7:
				if(x < 0){
					q[0] = 1.0 / (4 * M_PI);
					q[1] = -1;
					q[2] = 1;
					q[3] = -1;
					q[4] = 1.0;
					q[5] = 1.0 / sqrt(4 * M_PI);
					q[6] = -1.0 / sqrt(4 * M_PI);
					q[7] = 1.0 / sqrt(4 * M_PI);
				}else{
					q[0] = 1 / (4 * M_PI);
					q[1] = -1;
					q[2] = -1;
					q[3] = -1;
					q[4] = 1;
					q[5] = 1 / sqrt(4 * M_PI);
					q[6] = 1 / sqrt(4 * M_PI);
					q[7] = 1 / sqrt(4 * M_PI);
				}
				break;
			case test8:{
				if(x < 0){
					q[0] = 1.0;
					q[1] = -2.0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.45;
					q[5] = 0;
					q[6] = 0.5;
					q[7] = 0;
				}else{
					q[0] = 1.0;
					q[1] = 2.0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.45;
					q[5] = 0;
					q[6] = 0.5;
					q[7] = 0;
				}
				{
					double dlt = 4.0*dx;              // 抹宽约 4 个 cell(δ 越大越平滑)
					double vx  = 1 * sqrt(dx*dx);     // -2 -> +2 平滑过渡,去掉 x=0 速度间断
					q[0] = 1.0;
					q[1] = vx;
					q[2] = 0; q[3] = 0;
					q[4] = 0.45;
					q[5] = 0; q[6] = 0.5; q[7] = 0;
					break;
				}
			}
				break;
			case test9:
				if(x < 0){
					q[0] = 1.0;
					q[1] = 10;
					q[2] = 0;
					q[3] = 0;
					q[4] = 20.0;
					q[5] = 5.0 / sqrt(4 * M_PI);
					q[6] = 5.0 / sqrt(4 * M_PI);
					q[7] = 0;
				}else{
					q[0] = 1.0;
					q[1] = -10;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1.0;
					q[5] = 5 / sqrt(4 * M_PI);
					q[6] = 5 / sqrt(4 * M_PI);
					q[7] = 0;
				}
				break;
			case test10:
				if(x < 0){
					q[0] = 1.0;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1.0;
					q[5] = 0.7;
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 0.3;
					q[1] = 0.0;
					q[2] = 0.0;
					q[3] = 1.0;
					q[4] = 0.2;
					q[5] = 0.7;
					q[6] = 1.0;
					q[7] = 0;
				}
				break;
			case shocktest:
				if(x < 0){
					q[0] = 1;
					q[1] = 2;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 1;
					q[1] = -2;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}
				break;
			case convergencetest:
				
				result = GaussLegendreNumericalIntegration(x, dx, ConvergenceTestFunction, 2 );
				q[0] = result[0];
				q[1] = result[1];
				q[2] = result[2];
				q[3] = result[3];
				q[4] = result[4];
				q[5] = result[5];
				q[6] = result[6];
				q[7] = result[7];
				break;
			case test2d0:
				if(sqrt(x * x + y * y) <= 0.4){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 0.125;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}
				break;
			case testOrszag_Tang:
				q[0] = gamma_plasma * gamma_plasma;
				q[1] = -sin(2 * M_PI * y);
				q[2] = sin(2 * M_PI * x);
				q[3] = 0;
				q[4] = gamma_plasma;
				q[5] = -sin(2 * M_PI * y);
				q[6] = sin(4 * M_PI * x);
				q[7] = 0;
				break;
			case testblastwave: {
				double rad = sqrt(x*x + y*y);
				double rD = 0.1;
				double prsL = 1000.0, prsR = 0.1;
				double eps = 1 * sqrt(dx*dx + dy*dy);   // smoothing width
				
				q[0] = 1.0;
				q[1] = 0; q[2] = 0; q[3] = 0;
				q[4] = 0.5*(prsR+prsL) + 0.5*(prsR-prsL) * erf((rad - rD)/eps);
				q[5] = 70.0 / sqrt(4.0 * M_PI);   // B = 70 G (Boscheri-Thomann 2024), m = 70^2/(8 pi) = 195
				q[6] = 0.0;
				q[7] = 0.0;
				break;
			}
			case testblastwave1D:					
				if(fabs(x) <= 0.1){
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1000;
					q[5] = 1 / sqrt(4 * M_PI) * 100;
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 0.1;
					q[5] = 1 / sqrt(4 * M_PI) * 100;
					q[6] = 0;
					q[7] = 0;
				}
//					
				break;
			case testRotor:
				if(sqrt((x) * (x) + (y) * (y)) <= 0.1){
					q[0] = 10;
					q[1] = -10 * y;
					q[2] = 10 * x;
					q[3] = 0;
					q[4] = 1;
					q[5] = 2.5 / sqrt(4 * M_PI);
					q[6] = 0;
					q[7] = 0;
				}else{
					q[0] = 1;
					q[1] = 0;
					q[2] = 0;
					q[3] = 0;
					q[4] = 1;
					q[5] = 2.5 / sqrt(4 * M_PI);
					q[6] = 0;
					q[7] = 0;
				}
				
				break;
			case testFieldLoop:
				if(sqrt(x * x + y * y) < 0.3){
					q[0] = 1;
					q[1] = 2;
					q[2] = 1;
					q[3] = 0;
					q[4] = 1e5;
					q[5] = (-LOOP_A0 * y / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
					q[6] = (LOOP_A0 * x / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
					q[7] = 0;
				}else{
					q[0] = 1;
					q[1] = 2;
					q[2] = 1;
					q[3] = 0;
					q[4] = 1e5;
					q[5] = 0;
					q[6] = 0;
					q[7] = 0;
				}
				break;
			case testKHInstability:
				Mx = KH_MX;
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
				q[0] = 1.4;                         // rho = gamma = 1.4
				q[1] = Mx * (1.0 - 2.0 * eta);        // u = Mx(1 - 2 eta(y))
				q[2] = 0.1 * Mx * sin(2.0 * M_PI * x);// v = 0.1 Mx sin(2 pi x)
				q[3] = 0.0;                           // w = 0
				q[4] = 1.0;                           // p = 1
				q[5] = 0.1 * Mx;                      // Bx = 0.1 Mx (uniform)
				q[6] = 0.0;                           // By = 0
				q[7] = 0.0;   
				break;
			case testSmooth2D:   // [o3] entropy wave: rho advected by uniform (u,v)=(1,1), uniform p, B -> exact solution known
				q[0] = 1.0 + 0.2*sin(2.0*M_PI*(x + y)); q[1] = 1.0; q[2] = 1.0; q[3] = 0.0; q[4] = 0.1;
				q[5] = 0.1; q[6] = 0.1; q[7] = 0.0;
				break;
			case testVortex: {   // paper sec 4.1 analytic MHD vortex, u0=(1,1), p0=1, vt=sqrt(2pi), Bt=1, exact solution = translation
				const double r2 = x*x + y*y, ex = exp(0.5*(1.0 - r2)), vt = VORTEX_VT, Bt = g_vortex_Bt, rho0 = g_vortex_rho0;
				q[0] = rho0; q[1] = VORTEX_U0 + vt/(2.0*M_PI)*ex*(-y); q[2] = VORTEX_U0 + vt/(2.0*M_PI)*ex*x; q[3] = 0.0;
				q[4] = 1.0 + 0.5*exp(1.0 - r2)*(Bt*Bt/(4.0*M_PI*M_PI)*(1.0 - r2) - rho0*(vt/(2.0*M_PI))*(vt/(2.0*M_PI)));
				q[5] = Bt/(2.0*M_PI)*ex*(-y); q[6] = Bt/(2.0*M_PI)*ex*x; q[7] = 0.0;
				break; }
			case testdoubleshear:
				q[0] = 1;
				if(y < 0.5){
					q[1] = tanh(30 * (y - 0.25));
				}else{
					q[1] = tanh(30 * (0.75 - y));
				}
				q[2] = 0.05 * sin(2 * M_PI * x);
				q[4] = 1e5 / gamma_plasma;
				break;
			
		}
	};
	// ... averaged over each cell with 3x3 Gauss-Legendre in the CONSERVATIVE variables,
	// then converted back to primitives (upri2ucon in main reproduces the averages exactly).
	{
		const double gp[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
		const double gw[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double xc = x_0 + (i - ghost_cells + 0.5) * dx;
				double yc = y_0 + (j - ghost_cells + 0.5) * dy;
				double con[NUM_VARS] = {0,0,0,0,0,0,0,0}, q[NUM_VARS], cq[NUM_VARS], pri[NUM_VARS];
				if(test == convergencetest){          // 1D test: keeps its own 1D quadrature
					ic_point(xc, yc, q);
					for(int k = 0; k < NUM_VARS; k++) u(i,j,k) = q[k];
					continue;
				}
				for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++){
					ic_point(xc + 0.5*dx*gp[a], yc + 0.5*dy*gp[b], q);
					ini_pri2con(cq, q);
					for(int k = 0; k < NUM_VARS; k++) con[k] += gw[a]*gw[b]*cq[k];
				}
				ini_con2pri(pri, con);
				for(int k = 0; k < NUM_VARS; k++) u(i,j,k) = pri[k];
			}
		}
	}
	bound_cond_u(u);
	
	for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
			Bx_p0(i,j) = 0.5 * (u(i,j,5) + u(i + 1,j,5));
			By_0p(i,j) = 0.5 * (u(i,j,6) + u(i,j + 1,6));
		}
	}
	
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
			// Odd reflection.  At a reflecting wall normal to x one imposes
			// u -> -u, B_x -> -B_x, B_y -> +B_y.  With B_x = dA/dy and
			// B_y = -dA/dx this fixes the symmetry of A uniquely and in a
			// consistent way:  B_x odd requires dA/dy odd, hence A odd; and
			// B_y even requires dA/dx even, which under the mirror map (whose
			// Jacobian already contributes one sign change) again requires A
			// odd.  The additive constant in A is free, so the wall value may
			// be taken as zero and the reflection is A -> -A.
			//
			// A linear extrapolation here -- correct for a transmissive
			// boundary, where B is extrapolated at zeroth order -- would
			// produce a spurious normal field at the wall.
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < 3; k++){
					for(int g = 0; g < ghost_cells; g++){
						A(ghost_cells - 1 - g, j, k) = -A(ghost_cells + g, j, k);
					}
					int iR = nxCells + ghost_cells - 1;
					for(int g = 0; g < ghost_cells; g++){
						A(iR + 1 + g, j, k) = -A(iR - g, j, k);
					}
				}
			}
			break;
		}
		case x_Periodic: {
			// quasi-periodic (fix (2)): A_z carries -<By> x
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						const double jmp = (k == 2) ? g_Az_xjump : 0.0;   // A_z(x+Lx) = A_z(x) - <By> Lx
						A(ghost, j, k) = A(nxCells + ghost, j, k) + jmp;
						A(nxCells + ghost_cells + ghost, j, k) = A(ghost_cells + ghost, j, k) - jmp;
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
			// Odd reflection; see the x_Reflective branch for the derivation.
			// At a wall normal to y one imposes v -> -v, B_y -> -B_y,
			// B_x -> +B_x, which again forces A to be odd about the wall.
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int k = 0; k < 3; k++){
					for(int g = 0; g < ghost_cells; g++){
						A(i, ghost_cells - 1 - g, k) = -A(i, ghost_cells + g, k);
					}
					int jR = nyCells + ghost_cells - 1;
					for(int g = 0; g < ghost_cells; g++){
						A(i, jR + 1 + g, k) = -A(i, jR - g, k);
					}
				}
			}
			break;
		}
		case y_Periodic: {
			// quasi-periodic (fix (2)): A_z carries +<Bx> y
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						const double jmp = (k == 2) ? g_Az_yjump : 0.0;
						A(i, ghost, k) = A(i, nyCells + ghost, k) - jmp;
						A(i, nyCells + ghost_cells + ghost, k) = A(i, ghost_cells + ghost, k) + jmp;
					}
				}
			}
			break;
		}
	}
}

// ============================================================
// Dual-parameter epsilon computation (Document 3 sec.3.1), ported from
// your single-parameter epsilon-Euler version. Differences from the
// Euler original, and why:
//
//   1. calMach includes vz (con[3]): MHD carries a genuine 3rd (out-of-
//      plane) velocity component; M=||v||/a (Document 1 eq.6) uses the
//      full velocity, not just the in-plane vx,vy the Euler-only file has.
//   2. New calAlfvenMach / getTotalAlfvenSpeed: the existing getAlfvenSpeed()
//      in this file returns Bx/sqrt(rho) (the NORMAL/x-direction Alfven
//      speed -- used elsewhere for directional CFL), not the TOTAL
//      v_A=|B|/sqrt(rho) that M_A=||v||/v_A (Document 1 eq.6) needs. Wrote
//      a separate helper rather than repurpose getAlfvenSpeed, since other
//      code may depend on its existing (directional) meaning.
//   3. epsilon_min/epsilon_max didn't exist in this file -- added as
//      globals, defaulted to [0,1] (no artificial floor). Note: the WENO3
//      characteristic reconstruction's condition number was measured to
//      blow up for eps_min below ~0.05-0.1 (see earlier conversation) --
//      if that turns out to matter in practice, raising epsilon_min is a
//      one-line fix here, not a structural change.
//   4. bound_cond_array(Array&) only has the single-arg form here (unlike
//      epsilon-Euler's overload with xbc/ybc override params), and its
//      x_Reflective/y_Reflective branches NEGATE the value. Epsilon is a
//      positive weight, not a signed physical quantity -- negating it at
//      a reflective wall would be wrong. Added bound_cond_epsilon() as a
//      separate function (copy-based padding for Transmissive/Dirichlet/
//      Reflective, wrap-around for Periodic) rather than touch the shared
//      bound_cond_array used throughout the rest of the file.
// ============================================================

// EPSILON BOUNDS (per author instruction): epsilon_min MUST be > 0 when the
// explicit flux goes through HLLEM. At eps=0 the explicit split subsystem
// F^c + 0*F^p + 0*F^B degenerates to pure convection with an incomplete
// characteristic structure (no acoustic / magnetosonic / Alfven waves), so
// HLLEM's characteristic decomposition is not applicable. A strictly positive
// floor keeps a (scaled) full wave family alive everywhere.
// NOTE: kappa(R) ~ epsilon_min^{-4} for the split eigenvector matrix, so don't
// set this too small; 0.1 gives kappa ~ 1e4. TODO: sensitivity study.
// (eps==0 remains a valid configuration ONLY with the Rusanov flux path,
//  which needs no characteristic structure -- useful as a fully-implicit
//  baseline test.)
double epsilon_min = 0.1;
double epsilon_max = 1.0;

double getsoundspeed(double *con){
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	return sqrt(gamma_plasma * pri[4] / pri[0]);   // pri[4]=p, pri[0]=rho
}

double getTotalAlfvenSpeed(double *con){
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double B2 = pri[5]*pri[5] + pri[6]*pri[6] + pri[7]*pri[7];
	return sqrt(B2 / pri[0]);
}

double calMach(double *con){
	double soundspeed_0 = getsoundspeed(con);
	return sqrt(pow(con[1]/con[0],2) + pow(con[2]/con[0],2) + pow(con[3]/con[0],2)) / soundspeed_0;
}

double calAlfvenMach(double *con){
	double vA = fmax(getTotalAlfvenSpeed(con), 1e-12);   // guard B=0 (and v=0) -> MA=0 -> epsB floor, not NaN
	double v  = sqrt(pow(con[1]/con[0],2) + pow(con[2]/con[0],2) + pow(con[3]/con[0],2));
	return v / vA;
}

double g_sw_m0 = 0.15, g_sw_m1 = 0.40, g_sw_alpha = 14.0;   // switch parameters (paper: 0.15, 0.4, 14)
int g_epsmax3 = 0;   // 1: cell eps := max over the 3x3 neighbourhood before the face max
int g_pswitch = 0;   // 1: switching function driven by max(M, sqrt(|dp|/p)) and max(M_A, sqrt(|dm|/m))
inline double switchS(double M) {
	double SW_EPS0  = g_sw_m0;
	double SW_EPS1  = g_sw_m1;
	double SW_ALPHA = g_sw_alpha;
	
	if (M <= 0.0) return 1.0;
	if (M >= 1.0) return 0.0;
	
	if (M <= SW_EPS0) {
		return 1.0 - std::pow(M, SW_ALPHA);
	}
	if (M >= SW_EPS1) {
		return std::pow(1.0 - M, SW_ALPHA);
	}
	const double u    = (M - SW_EPS0) / (SW_EPS1 - SW_EPS0);
	const double bump = std::exp(1.0 - 1.0 / (1.0 - u * u));
	const double A    = (1.0 - std::pow(SW_EPS0, SW_ALPHA))
	-  std::pow(1.0 - SW_EPS1, SW_ALPHA);
	const double B    =  std::pow(1.0 - SW_EPS1, SW_ALPHA);
	return bump * A + B;
}

// M 越靠近/超过 1(声速量级),epsM 越大,越倾向显式
double f_epsM(double* con) {
	const double M   = calMach(con);
	const double eps = 1.0 - switchS(M);
	return std::fmin(eps_M_max, std::fmax(eps_M_min, eps));
}
// M_A 越靠近/超过 1(Alfven 速度量级),epsB 越大,越倾向显式
double f_epsB(double* con) {
	const double MA  = calAlfvenMach(con);
	const double eps = 1.0 - switchS(MA);
	return std::fmin(eps_B_max, std::fmax(eps_B_min, eps));
}

// epsilon 是正权重场,反射边界不应取负——单独一份,不碰共用的 bound_cond_array
void bound_cond_epsilon(Array &u){
	switch(x_bc){
		case x_Transmissive: case x_Dirichlet: case x_Reflective:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < ghost_cells; k++){
					u(k,j) = u(ghost_cells,j);
					u(nxCells + ghost_cells + k,j) = u(nxCells + ghost_cells - 1,j);
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
		case y_Transmissive: case y_Dirichlet: case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = u(i,ghost_cells);
					u(i,nyCells + ghost_cells + ghost) = u(i,nyCells + ghost_cells - 1);
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

// 双参数版本:同一次 cell 扫描里,epsM 用 calMach,epsB 用 calAlfvenMach,
// face 值都用 upwind-max(和你 Euler 版本、以及本文件其余地方的惯例一致)。
void getepsgrid(Grid &u_con, Array &eps_M, Array &eps_B,Array &eps_M_grid_x, Array &eps_M_grid_y,
                              Array &eps_B_grid_x, Array &eps_B_grid_y){
	static Array cM(nxCells, nyCells), cB(nxCells, nyCells);   // cell values on the interior plus one layer at +x, +y
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j <= ghost_cells + nyCells; j++)
		for(int i = ghost_cells; i <= ghost_cells + nxCells; i++){
			if(i == ghost_cells + nxCells && j == ghost_cells + nyCells) continue;
			double con[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++) con[k] = u_con(i,j,k);
			cM(i,j) = f_epsM(con); cB(i,j) = f_epsB(con);
			if(g_pswitch){   // [fix B] pressure-fluctuation-aware switch: M_eff = max(M, sqrt(max_nbr |dp|/p)), M_A,eff = max(M_A, sqrt(max_nbr |dm|/m))
				double pri[NUM_VARS]; ini_con2pri(pri, con);
				const double pc = fmax(pri[4], 1e-300), mc = fmax(0.5*(pri[5]*pri[5]+pri[6]*pri[6]+pri[7]*pri[7]), 1e-300);
				double dp = 0.0, dm = 0.0;
				const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
				for(int q = 0; q < 4; q++){
					const int ii = std::max(0, std::min(nxCells + 2*ghost_cells - 1, i + di[q])), jj = std::max(0, std::min(nyCells + 2*ghost_cells - 1, j + dj[q]));
					double cn[NUM_VARS], pn[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) cn[k] = u_con(ii,jj,k); ini_con2pri(pn, cn);
					const double pnb = fmax(pn[4], 1e-300), mnb = fmax(0.5*(pn[5]*pn[5]+pn[6]*pn[6]+pn[7]*pn[7]), 1e-300);
					dp = fmax(dp, fabs(pnb - pc)/fmin(pnb, pc)); dm = fmax(dm, fabs(mnb - mc)/fmin(mnb, mc));
				}
				const double Mp = sqrt(dp), Mm = sqrt(dm);
				const double eM = 1.0 - switchS(fmax(calMach(con), Mp)), eB = 1.0 - switchS(fmax(calAlfvenMach(con), Mm));
				cM(i,j) = fmin(eps_M_max, fmax(eps_M_min, eM)); cB(i,j) = fmin(eps_B_max, fmax(eps_B_min, eB));
			}
		}
	if(g_epsmax3){   // optional: widen eps by one cell in every direction (3x3 max), as the shock indicator chi
		static Array cM3(nxCells, nyCells), cB3(nxCells, nyCells);
		bound_cond_epsilon(cM); bound_cond_epsilon(cB);
#pragma omp parallel for schedule(static)
		for(int j = ghost_cells - 1; j <= ghost_cells + nyCells; j++)
			for(int i = ghost_cells - 1; i <= ghost_cells + nxCells; i++){
				double m = 0.0, b = 0.0;
				for(int dj = -1; dj <= 1; dj++) for(int di = -1; di <= 1; di++){ m = fmax(m, cM(i+di,j+dj)); b = fmax(b, cB(i+di,j+dj)); }
				cM3(i,j) = m; cB3(i,j) = b;
			}
#pragma omp parallel for schedule(static)
		for(int j = ghost_cells - 1; j <= ghost_cells + nyCells; j++)
			for(int i = ghost_cells - 1; i <= ghost_cells + nxCells; i++){ cM(i,j) = cM3(i,j); cB(i,j) = cB3(i,j); }
	}
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			eps_M_grid_x(i,j) = fmax(cM(i,j), cM(i+1,j));
			eps_M_grid_y(i,j) = fmax(cM(i,j), cM(i,j+1));
			
			eps_B_grid_x(i,j) = fmax(cB(i,j), cB(i+1,j));
			eps_B_grid_y(i,j) = fmax(cB(i,j), cB(i,j+1));
			
			eps_M(i,j) = cM(i,j);
			eps_B(i,j) = cB(i,j);
		}
	}
	
	bound_cond_epsilon(eps_M_grid_x);
	bound_cond_epsilon(eps_M_grid_y);
	bound_cond_epsilon(eps_B_grid_x);
	bound_cond_epsilon(eps_B_grid_y);
}

double caldt(Grid &u_pri, Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y){
	vx_max = 0, vy_max = 0;
	double dt;
	double pri[NUM_VARS], con[NUM_VARS];
	double smax_x = 0, smax_y = 0;
	
	{ double vxm = 0.0, vym = 0.0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:vxm,vym)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			vxm = fmax(vxm, fabs(u_pri(i,j,1)));   // (= the former |u| > |vx_max| test: vx_max >= 0)
			vym = fmax(vym, fabs(u_pri(i,j,2)));
		}
	}
	vx_max = vxm; vy_max = vym; }
	
	{	// time step from the eps-split explicit spectrum (paper eq. 40); the first full_step steps use the full MHD spectrum
			if(timestep < full_step + 1){ a_max_x = ch_x; a_max_y = ch_y; }
			else{
#pragma omp parallel for collapse(2) schedule(static) reduction(max:smax_x,smax_y)
				for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
					for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
						Vec8 lam_x, lam_y; double pri[NUM_VARS], con[NUM_VARS];
						for(int k = 0; k < NUM_VARS; k++) pri[k] = u_pri(i,j,k);
						ini_pri2con(con, pri);
						
						double eMx = fmax(eps_M_grid_x(i,j), eps_M_grid_x(i-1,j));
						double eBx = fmax(eps_B_grid_x(i,j), eps_B_grid_x(i-1,j));
						double eMy = fmax(eps_M_grid_y(i,j), eps_M_grid_y(i,j-1));
						double eBy = fmax(eps_B_grid_y(i,j), eps_B_grid_y(i,j-1));
						
						compute_eigenvalues_x(con, eMx, eBx, lam_x);
						compute_eigenvalues_y(con, eMy, eBy, lam_y);
						
						for(int k = 0; k < NUM_VARS; k++){
							smax_x = std::max(smax_x, std::fabs(lam_x(k)));
							smax_y = std::max(smax_y, std::fabs(lam_y(k)));
						}
					}
				}
				a_max_x = smax_x; a_max_y = smax_y;
			}
	}
	
	dt = CFL * 1 / (a_max_x / dx + a_max_y / dy);
	return dt;
}

void calch(Grid &u_pri){
	double ch_grid_x = 1e-8, ch_grid_y = 1e-8, ch_grid = 1e-8;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:ch_grid_x,ch_grid_y,ch_grid)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			double rho = u_pri(i,j,0);
			double p = u_pri(i,j,4);
			double Bx = u_pri(i,j,5);
			double By = u_pri(i,j,6);
			double Bz = u_pri(i,j,7);
			double B2 = Bx * Bx + By * By + Bz * Bz;
			double B = sqrt(B2);
			double v = sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			double cf_x = sqrt(0.5 * (((gamma_plasma * p + B2)/(rho)) + sqrt(fmax( ((gamma_plasma * p + B2) / rho) * ((gamma_plasma * p + B2) / rho) - 4 * ((gamma_plasma * p * Bx * Bx)/(rho * rho)), 0.0))));
			double cf_y = sqrt(0.5 * (((gamma_plasma * p + B2)/(rho)) + sqrt(fmax( ((gamma_plasma * p + B2) / rho) * ((gamma_plasma * p + B2) / rho) - 4 * ((gamma_plasma * p * By * By)/(rho * rho)), 0.0))));
			
			double cf = sqrt(0.5 * (((gamma_plasma * p + B2)/(rho)) + sqrt(fmax( ((gamma_plasma * p + B2) / rho) * ((gamma_plasma * p + B2) / rho) - 4 * ((gamma_plasma * p * B * B)/(rho * rho)), 0.0))));
			
			if(ch_grid_x < fabs(u_pri(i,j,1)) + cf_x){   // [o3q] |u| (was signed u: flows moving in -x/-y underestimated)
				ch_grid_x = fabs(u_pri(i,j,1)) + cf_x;
			}
			if(ch_grid_y < fabs(u_pri(i,j,2)) + cf_y){
				ch_grid_y = fabs(u_pri(i,j,2)) + cf_y;
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
#pragma omp parallel for collapse(2) schedule(static) reduction(max:ca_grid_x,ca_grid_y,ca_grid)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			double pri[NUM_VARS];
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

double calkinetic(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double kinetic_energy = 0.5 * rho*(vx * vx + vy * vy + vz * vz);
	return kinetic_energy;
}

void MagneticEigenvalues(Grid &u_pri){
	double bmax_x = 0;
	double bmax_y = 0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:bmax_x,bmax_y)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
			double b_x, b_y;
			b_x = 0.5 * (fabs(u_pri(i,j,1)) + sqrt(u_pri(i,j,1) * u_pri(i,j,1) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			
			b_y = 0.5 * (fabs(u_pri(i,j,2)) + sqrt(u_pri(i,j,2) * u_pri(i,j,2) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			bmax_x = fmax(bmax_x, b_x);   // (= the former |bmax| < |b| test: b >= 0)
			bmax_y = fmax(bmax_y, b_y);
		}
	}
	B_max_x = bmax_x;
	B_max_y = bmax_y;
}

void EnergyEigenvalues(Grid &u_pri){
	double umax_x = 0, umax_y = 0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:umax_x,umax_y)
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
			double u_x, u_y;
			double c = sqrt(gamma_plasma * u_pri(i,j,4) / u_pri(i,j,0));
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

// (face4/corner4/stagger4 templates are defined earlier, before Staggered_CT)

// Validated separately (standalone harness): face4/corner4 exact for bicubic
// polynomials + clean 4th-order convergence on smooth non-polynomial test functions;
// stagger4 exact for cubics + clean 4th-order convergence; the FULL nested
// construction used below (FluxBx4, see Magnetic_Subsystem_AX) gives observed
// self-convergence order -> 3.996, 3.999 on a smooth manufactured MHD-like field.
//
// IMPORTANT: naively differencing two ADJACENT 4th-order-accurate face/corner
// values (Fp12-Fm12)/h is only 2nd order overall (verified: order->2.00 exactly),
// regardless of how accurate Fp12/Fm12 themselves are -- that's why stagger4
// needs FOUR half-integer points, not two.
// ================================================================================
template<typename Fn>
inline double face4x(Fn get, int i, int j){ return (-get(i-1,j) + 9.0*get(i,j) + 9.0*get(i+1,j) - get(i+2,j)) / 16.0; }
template<typename Fn>
inline double face4y(Fn get, int i, int j){ return (-get(i,j-1) + 9.0*get(i,j) + 9.0*get(i,j+1) - get(i,j+2)) / 16.0; }
template<typename Fn>
inline double corner4(Fn get, int i, int j){   // corner at (i+1/2, j+1/2), tensor product of face4x
	double F0 = face4x(get, i, j-1), F1 = face4x(get, i, j), F2 = face4x(get, i, j+1), F3 = face4x(get, i, j+2);
	return (-F0 + 9.0*F1 + 9.0*F2 - F3) / 16.0;
}
inline double stagger4(double Fm32, double Fm12, double Fp12, double Fp32, double invh){
	return (Fm32 - 27.0*Fm12 + 27.0*Fp12 - Fp32) * invh / 24.0;
}

void A_initial(Grid &A){
	double C;
	double R;
	double r;
	// 注意: i 和 j 都要覆盖全部网格(包括ghost cells)
	auto A_point = [&](double x, double y, double& Ax_val, double& Ay_val, double& Az_val){
		Ax_val = 0.0; Ay_val = 0.0; Az_val = 0.0;
		switch(test) {
			case test0: 
			case shocktest:
			case convergencetest: 
			case test2d0:
			case testSmooth2D: Az_val = 0.1*y - 0.1*x; break;
			case testVortex: Az_val = g_vortex_Bt/(2.0*M_PI)*exp(0.5*(1.0 - x*x - y*y)); break;
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
				double C = LOOP_A0 / sqrt(4.0*M_PI);
				double R = 0.3;
				double r = sqrt(x*x + y*y);
				Az_val = (r < R) ? C*(R-r) : 0.0;
				break;
			}
				
			case testKHInstability:
				Az_val = 0.1 * KH_MX * y;
				break;
			
			default:
				break;
		}
	};
	// [o3] cell averages of the potential (3x3 Gauss): curl4 of cell averages is the
	// cell average of B to 4th order (linear operators commute with averaging)
	const double gp[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	const double gw[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
	for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
		for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
			double xc = x_0 + (i - ghost_cells + 0.5) * dx;
			double yc = y_0 + (j - ghost_cells + 0.5) * dy;
			double sx = 0.0, sy = 0.0, sz = 0.0, ax, ay, az;
			for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++){
				A_point(xc + 0.5*dx*gp[a], yc + 0.5*dy*gp[b], ax, ay, az);
				sx += gw[a]*gw[b]*ax; sy += gw[a]*gw[b]*ay; sz += gw[a]*gw[b]*az;
			}
			A(i, j, 0) = sx; A(i, j, 1) = sy; A(i, j, 2) = sz;
		}
	}
	bound_cond_A(A);
}

inline double minmod(double a, double b) {
	if (a * b <= 0.0) return 0.0;
	return (fabs(a) < fabs(b)) ? a : b;
}

// Density point value wherever a velocity m/rho is formed from point values (cellavg_velocity for the A transport,
// the magnetic operators via mfv_prepare, MagneticFlux_Update_FV, Magnetic_Sweep_FV, kemb_lin).
// The unlimited deconvolution rho - (d2x + d2y) rho / 24 can take a cell that is a narrow density minimum between
// much denser neighbours to rho_pt ~ 0 (measured on Orszag-Tang 1024^2 at t = 0.3234: rho_bar 0.94 with neighbours
// 2.9-10.5 gives rho_pt = 6e-3, |u_pt| = 38, and one stage later rho_pt <= 0, u ~ 1e300, E ~ 1e291).  Where the
// deconvolution removes more than half of rho, the cell average is used instead, the same admissibility threshold as
// in fv_cell_lim.  rho is fixed within a stage, so the functional stays affine in (m, B); in smooth flow the guard
// never acts and the scheme is unchanged.
static void rho_pt_guard(Grid& Qpt, Grid& Qbar){   // (also available inside fv_cell_to)
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	long n = 0;
#pragma omp parallel for collapse(2) schedule(static) reduction(+:n)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++)
		if(!(Qpt(i,j,0) >= 0.5*Qbar(i,j,0))){ Qpt(i,j,0) = Qbar(i,j,0); n++; }
	g_rho_guard += n;
}
// [o3j] CELL AVERAGES of the velocity u = m/rho via point values (ubar != mbar/rhobar at O(h^2) when rho varies):
// Qbar -> Q_pt (deconvolve, limited) -> u_pt = m_pt/rho_pt -> ubar (convolve).  Used by the A transport.
static void cellavg_velocity(Grid& u, Array& ub, Array& vb){
	static FieldN<3> Qpt; static Array Tu(nxCells, nyCells), Tv(nxCells, nyCells);   // point values of rho, m_x, m_y
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	fvN_to<3>(Qpt, GridSlots{u, 0}, -1.0, true);   // [emb] linear in the momentum; density point-value guard
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++){
		double r = fmax(Qpt(i,j,0), 1e-300);
		Tu(i,j) = Qpt(i,j,1)/r; Tv(i,j) = Qpt(i,j,2)/r;
	}
	fv_cell_s(ub, Tu, +1.0); fv_cell_s(vb, Tv, +1.0);
}
void Unstaggered_CT(Grid& u_con, Grid& A,Grid& A_s, Grid& A_e, double h) {
	// ================================================================
	// UNSTAGGERED (POTENTIAL-BASED) CT, following Rossmanith's MPACT as
	// adapted in Dematte et al.
	//
	// The z-component A of the vector potential obeys the pure ADVECTION
	// equation
	//        dA/dt + u dA/dx + v dA/dy = 0.
	// This is exactly dA/dt = -Ez rewritten: with Bx = dA/dy and By = -dA/dx,
	//        Ez = v Bx - u By = u dA/dx + v dA/dy.
	// Solving it directly in advective form -- rather than assembling a face
	// Ez and interpolating it down to the cell centre -- is what MPACT does,
	// and it is the right discretisation here.  The face->corner->cell chain
	// is central; its transfer function rolls off at Nyquist, so no upwind
	// dissipation reaches A and the grid-scale modes of A are undamped.  In
	// advective form the upwinding acts on A itself and the limiter applies
	// to it directly, which is why the earlier Ez-based version needed an
	// ad-hoc Rusanov term bolted onto A to stay clean.
	//
	// Note on the eps-splitting: A is advected with the velocity of the FINAL
	// state, i.e. after the explicit advective update and both implicit
	// sub-systems.  It therefore carries all three shares by which B is
	// advanced without any of them having to be weighted individually -- the
	// flux-dual bookkeeping needed by the staggered variant does not arise
	// here, because the advection equation is the physical one.
	//
	// RK PLACEMENT (fix (1), 2026-09-25).  The stage derivative stored in HsA is
	//        k_i = (A_i - A_s^i)/(a_ii dt) = -u_i . grad(A_e^i),
	// i.e. the transport term is evaluated on the EXPLICIT stage combination
	// A_e^i = A^n + dt sum_{j<i} ae_ij k_j, with the velocity of the stage
	// solution.  The update base stays the implicit combination A_s^i, so the
	// stiffly-accurate last stage still returns A^{n+1} = A^n + dt sum_i b_i k_i,
	// which is now the 3rd-order explicit method (Ae,b) of the IMEX pair.
	// (A is advected at the MATERIAL speed, so an explicit treatment costs
	//  nothing: the Alfven stiffness sits in the Lorentz force, not here.)
	// Previously the faces were built from A_s: with no implicit solve that is
	// forward Euler on the implicit base, whose effective tableau is the
	// strictly-lower part of Ai -- first order (sum b_i c_i = 0.064) and with
	// stability polynomial 1 + z + 0.064 z^2 - 0.079 z^3, unstable for
	// advection.  The field loop then GAINS magnetic energy (+19% at 256x128,
	// +75% at 512x256 by t=1) and develops the ring/striation pattern.
	//
	// B = curl A with 4th-order central differences, so that div(curl) is an
	// exact identity under the MATCHING 4th-order central divergence; a
	// 2nd-order divergence diagnostic reads O(dx^2) -- that is the measuring
	// stick, not the field.
	// ================================================================
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;

	auto As = [&](int ii,int jj)->double{           // faces from the EXPLICIT combination A_e (fix (1))
		int ic = std::max(0, std::min(NX-1, ii));
		int jc = std::max(0, std::min(NY-1, jj));
		return A_e(ic, jc, 2);
	};
	auto Uc = [&](int ii,int jj,int k)->double{
		int ic = std::max(0, std::min(NX-1, ii));
		int jc = std::max(0, std::min(NY-1, jj));
		return u_con(ic, jc, k);
	};

	// ---- advective update of A ----
	// Non-conservative form: A is not a density, so the update is
	// u dA/dx + v dA/dy and NOT d(uA)/dx + d(vA)/dy; the two differ by
	// A div(v), which does not vanish in compressible MHD.
	// [o3] (i) upwind derivatives of A by HJ-WENO5-Z (Jiang-Peng 2000 stencils, Borges Z-weights,
	//          eps = 1e-40): the 3-point WENO3 was anisotropic on the cone (interior lobes) and
	//          2nd order with strong dissipation at the kink (field loop rim ~5 cells);
	//      (ii) the cell average of the product u.grad(A) is NOT ubar.gradAbar: the
	//          covariance (h^2/12) grad(u).grad(G) is added with minmod-limited slopes,
	//          which is O(h^3)-exact in smooth regions and vanishes across a jump.
	auto wz5 = [&](double v1, double v2, double v3, double v4, double v5)->double{
		const double ew5 = 1e-40;
		double b1 = 13.0/12.0*(v1-2*v2+v3)*(v1-2*v2+v3) + 0.25*(v1-4*v2+3*v3)*(v1-4*v2+3*v3);
		double b2 = 13.0/12.0*(v2-2*v3+v4)*(v2-2*v3+v4) + 0.25*(v2-v4)*(v2-v4);
		double b3 = 13.0/12.0*(v3-2*v4+v5)*(v3-2*v4+v5) + 0.25*(3*v3-4*v4+v5)*(3*v3-4*v4+v5);
		double tau = fabs(b1-b3);
		double a1 = 0.1*(1+(tau/(b1+ew5))*(tau/(b1+ew5))), a2 = 0.6*(1+(tau/(b2+ew5))*(tau/(b2+ew5))), a3 = 0.3*(1+(tau/(b3+ew5))*(tau/(b3+ew5)));
		return (a1*(v1/3-7*v2/6+11*v3/6) + a2*(-v2/6+5*v3/6+v4/3) + a3*(v3/3+5*v4/6-v5/6))/(a1+a2+a3);
	};
	auto Dx = [&](int i,int j,int off)->double{ return (As(i+off+1,j) - As(i+off,j))/dx; };   // D_{i+off+1/2}
	auto Dy = [&](int i,int j,int off)->double{ return (As(i,j+off+1) - As(i,j+off))/dy; };
	auto Gx_at = [&](int i,int j)->double{
		double u = (g_emb_Qe) ? (*g_emb_Qe)(i,j,1)/(*g_emb_Qe)(i,j,0) : u_con(i,j,1)/u_con(i,j,0);   // [emb]
		return (u >= 0.0) ? wz5(Dx(i,j,-3),Dx(i,j,-2),Dx(i,j,-1),Dx(i,j,0),Dx(i,j,1))
		                  : wz5(Dx(i,j,2),Dx(i,j,1),Dx(i,j,0),Dx(i,j,-1),Dx(i,j,-2));
	};
	auto Gy_at = [&](int i,int j)->double{
		double v = (g_emb_Qe) ? (*g_emb_Qe)(i,j,2)/(*g_emb_Qe)(i,j,0) : u_con(i,j,2)/u_con(i,j,0);   // [emb]
		return (v >= 0.0) ? wz5(Dy(i,j,-3),Dy(i,j,-2),Dy(i,j,-1),Dy(i,j,0),Dy(i,j,1))
		                  : wz5(Dy(i,j,2),Dy(i,j,1),Dy(i,j,0),Dy(i,j,-1),Dy(i,j,-2));
	};
	auto mm2g = [&](double a, double b, double thr)->double{ if(fabs(a) < thr && fabs(b) < thr) return 0.5*(a+b); return (a*b > 0.0) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0; };
	// G = grad A_e (upwind WENO5) and its limited slopes depend only on A_e, on the upwind direction (from u_E when
	// g_emb_Qe is set) and on eps_M: all fixed within a stage, so they are computed once per stage and reused by every
	// CT of that stage (the embedded functional evaluates the CT once per defect-correction sweep).
	static Array Gx(nxCells, nyCells), Gy(nxCells, nyCells);
	static Array Sxx(nxCells, nyCells), Sxy(nxCells, nyCells), Syx(nxCells, nyCells), Syy(nxCells, nyCells);   // limited slopes of G
	static long G_stamp = -1; static const void* G_key[3] = {nullptr, nullptr, nullptr};
	const void* key[3] = {&A_e, g_emb_Qe, g_epsM_cell};
	const bool cacheable = (g_emb_Qe != nullptr);
	if(!(cacheable && G_stamp == g_stage_stamp && G_key[0] == key[0] && G_key[1] == key[1] && G_key[2] == key[2])){
#pragma omp parallel for collapse(2) schedule(static)
		for(int j = ghost_cells-1; j < nyCells + ghost_cells + 1; j++)
			for(int i = ghost_cells-1; i < nxCells + ghost_cells + 1; i++){ Gx(i,j) = Gx_at(i,j); Gy(i,j) = Gy_at(i,j); }
		double S_gx = 0.0, S_gy = 0.0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:S_gx,S_gy)
		for(int j = ghost_cells; j < nyCells+ghost_cells; j++) for(int i = ghost_cells; i < nxCells+ghost_cells; i++){ S_gx = fmax(S_gx, fabs(Gx(i,j))); S_gy = fmax(S_gy, fabs(Gy(i,j))); }
#pragma omp parallel for collapse(2) schedule(static)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				double wq = g_epsM_cell ? fmin(1.0, fmax(0.0, (*g_epsM_cell)(i,j))) : 1.0;
				auto sl = [&](double a, double b, double thr)->double{ return wq * mm2g(a, b, thr) + (1.0 - wq) * 0.5*(a+b); };
				// [o3t] the slopes of G = grad A use the SAME eps_M-blended, threshold-gated limiter as the slopes of u:
				// the unconditional minmod zeroed them at every smooth extremum of G (O(h^2) in h-wide bands, B residual
				// ~2.5th order on the stationary vortex); the gate still limits across a jump (|dG| > S_G/2).
				Sxx(i,j) = sl(Gx(i+1,j)-Gx(i,j), Gx(i,j)-Gx(i-1,j), 0.5*S_gx)/dx; Sxy(i,j) = sl(Gx(i,j+1)-Gx(i,j), Gx(i,j)-Gx(i,j-1), 0.5*S_gx)/dy;
				Syx(i,j) = sl(Gy(i+1,j)-Gy(i,j), Gy(i,j)-Gy(i-1,j), 0.5*S_gy)/dx; Syy(i,j) = sl(Gy(i,j+1)-Gy(i,j), Gy(i,j)-Gy(i,j-1), 0.5*S_gy)/dy;
			}
		if(cacheable){ G_stamp = g_stage_stamp; for(int q = 0; q < 3; q++) G_key[q] = key[q]; } else G_stamp = -1;
	}
	static Array Ubar(nxCells, nyCells), Vbar(nxCells, nyCells);
	cellavg_velocity(u_con, Ubar, Vbar);   // [o3j] cell-average velocity (was mbar/rhobar: O(h^2) off where rho varies)
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			auto U = [&](int ii,int jj)->double{ return Ubar(ii,jj); };
			auto V = [&](int ii,int jj)->double{ return Vbar(ii,jj); };
			double u = U(i,j), v = V(i,j);
			double T = u*Gx(i,j) + v*Gy(i,j);
			// the transport must be LINEAR in the velocity (it is embedded in the pressure operator): central slopes of u
			const double ux = 0.5*(U(i+1,j)-U(i-1,j))/dx, uy = 0.5*(U(i,j+1)-U(i,j-1))/dy, vx = 0.5*(V(i+1,j)-V(i-1,j))/dx, vy = 0.5*(V(i,j+1)-V(i,j-1))/dy;
			const double gxx = Sxx(i,j), gxy = Sxy(i,j), gyx = Syx(i,j), gyy = Syy(i,j);
			T += dx*dx/12.0*(ux*gxx + vx*gyx) + dy*dy/12.0*(uy*gxy + vy*gyy);
			A(i,j,2) = A_s(i,j,2) - h * T;
		}
	bound_cond_A(A);

	// ---- B = curl A, 4th-order central ----
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			u_con(i,j,5) =  ( -A(i,j+2,2) + 8.0*A(i,j+1,2) - 8.0*A(i,j-1,2) + A(i,j-2,2) ) / (12.0*dy);
			u_con(i,j,6) = -( -A(i+2,j,2) + 8.0*A(i+1,j,2) - 8.0*A(i-1,j,2) + A(i-2,j,2) ) / (12.0*dx);
		}
	bound_cond_u(u_con);
}

int M_P_bc(int i, int j){
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

int M_B_bc(int i, int j, int k){
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

// ================================================================================
// 4th-order interpolation/differentiation primitives for the magnetic operator.
// Validated separately (standalone harness): face4/corner4 exact for bicubic
// polynomials + clean 4th-order convergence on smooth non-polynomial test functions;
// stagger4 exact for cubics + clean 4th-order convergence; the FULL nested
// construction used below (FluxBx4, see Magnetic_Subsystem_AX) gives observed
// self-convergence order -> 3.996, 3.999 on a smooth manufactured MHD-like field.
//
// IMPORTANT: naively differencing two ADJACENT 4th-order-accurate face/corner
// values (Fp12-Fm12)/h is only 2nd order overall (verified: order->2.00 exactly),
// regardless of how accurate Fp12/Fm12 themselves are -- that's why stagger4
// needs FOUR half-integer points, not two.
// (face4/corner4/stagger4 templates defined earlier, before Staggered_CT)

// ================================================================================
// Shared field-accessor functors + per-component 4th-order Ax contribution
// functions for the magnetic operator. Refactored out of the original single
// ~300-line loop body; validated as a byte-identical (0.000e+00 diff) restructuring
// against the previous lambda-based implementation before being spliced in here.
//
// Design: face4x/face4y/corner4 (above) are templates needing any callable with
// operator(int,int)->double, so the field accessors are small functor structs
// (not lambdas) -- lets them be shared as ordinary values across all three
// component functions without std::function overhead or capture-list duplication.
// ================================================================================

// n_field: background ("n") state read from Qe, component comp (5=Bx,6=By,7=Bz)
struct NField4 {
	Grid &Qe; int comp;
	double operator()(int i, int j) const { return Qe(i,j,comp); }
};
// rho_field: background density from u_con
struct RhoField4 {
	Grid &u_con;
	double operator()(int i, int j) const { return u_con(i,j,0); }
};
// x_field: the UNKNOWN vector being solved for (GMRES trial vector), component comp.
// Needs its own reflection-sign handling (generalizes reflct_xL/xR/yL/yR to
// ARBITRARY offset -- 4th order reaches i+-3,j+-3, not just +-1). Qe/u_con don't
// need this: bound_cond_u already populated their ghosts with the correct sign,
// unlike the raw GMRES solve vector x.
struct XField4 {
	Eigen::VectorXd &x; int comp;
	static double xsign(int ii){
		if(x_bc == x_Reflective && (ii < ghost_cells || ii >= nxCells + ghost_cells)) return -1.0;
		return 1.0;
	}
	static double ysign(int jj){
		if(y_bc == y_Reflective && (jj < ghost_cells || jj >= nyCells + ghost_cells)) return -1.0;
		return 1.0;
	}
	double operator()(int i, int j) const { return xsign(i)*ysign(j) * x(M_B(i,j,comp)); }
};

// Bundles the six field accessors + the cross-component helper quantities that
// Bx/By/Bz's 4th-order construction share (Efun/Hfun used by 2 components each,
// the corner/cell product combos used by whichever components need them).
struct MagAX4Fields {
	NField4 nBx, nBy, nBz;
	RhoField4 nrho;
	XField4 xBx, xBy, xBz;

	MagAX4Fields(Grid &Qe, Grid &u_con, Eigen::VectorXd &x)
	  : nBx{Qe,5}, nBy{Qe,6}, nBz{Qe,7}, nrho{u_con}, xBx{x,5}, xBy{x,6}, xBz{x,7} {}

	// E(i,j) = 0.5*(Bxn.Bx+Byn.By+Bzn.Bz) - Byn*By, cell-centered (Bx, Bz share this)
	double Efun(int i,int j) const {
		double bxn=nBx(i,j), byn=nBy(i,j), bzn=nBz(i,j);
		double bx=xBx(i,j), by=xBy(i,j), bz=xBz(i,j);
		return 0.5*(bxn*bx+byn*by+bzn*bz) - byn*by;
	}
	// H(i,j) = 0.5*(...) - Bxn*Bx, cell-centered (By, Bz share this)
	double Hfun(int i,int j) const {
		double bxn=nBx(i,j), byn=nBy(i,j), bzn=nBz(i,j);
		double bx=xBx(i,j), by=xBy(i,j), bz=xBz(i,j);
		return 0.5*(bxn*bx+byn*by+bzn*bz) - bxn*bx;
	}
	// G at corner (i+1/2,j+1/2) = 0.5*(...) - Bxn*Bx (Bx only)
	double Gfun(int i,int j) const {
		double bxn=corner4(nBx,i,j), byn=corner4(nBy,i,j), bzn=corner4(nBz,i,j);
		double bx=corner4(xBx,i,j), by=corner4(xBy,i,j), bz=corner4(xBz,i,j);
		return 0.5*(bxn*bx+byn*by+bzn*bz) - bxn*bx;
	}
	// Gy at corner (i+1/2,j+1/2) = 0.5*(...) - Byn*By (By only)
	double Gyfun(int i,int j) const {
		double bxn=corner4(nBx,i,j), byn=corner4(nBy,i,j), bzn=corner4(nBz,i,j);
		double bx=corner4(xBx,i,j), by=corner4(xBy,i,j), bz=corner4(xBz,i,j);
		return 0.5*(bxn*bx+byn*by+bzn*bz) - byn*by;
	}
	double BxnBy(int i,int j)        const { return nBx(i,j)*xBy(i,j); }
	double BynBx(int i,int j)        const { return nBy(i,j)*xBx(i,j); }
	double BxnBz_cell(int i,int j)   const { return nBx(i,j)*xBz(i,j); }
	double BynBz_cell(int i,int j)   const { return nBy(i,j)*xBz(i,j); }
	double BynBx_corner(int i,int j) const { return corner4(nBy,i,j)*corner4(xBx,i,j); }
	double ByBz_corner(int i,int j)  const { return corner4(nBy,i,j)*corner4(xBz,i,j); }
	double BxnBz_corner(int i,int j) const { return corner4(nBx,i,j)*corner4(xBz,i,j); }
};

// ---- Bx: Ax(M_B(i,j,5)) contribution (outer y-derivative) ----

// ---- By: Ax(M_B(i,j,6)) contribution (outer x-derivative) ----
// NOT a direct i<->j swap of Bx's field roles -- checked the original By_1~By_4
// term by term (e.g. its coefficient is Byn but its corner product is still
// Bxn*By, not Byn*Bx), so this follows the actual formulas, not assumed symmetry.

// ---- Bz: Ax(M_B(i,j,7)) contribution (x-face flux + y-face flux, summed) ----
inline double FluxBzx4(const MagAX4Fields &F, int j, int ic){
	double cyd_BynBx = stagger4(F.BynBx_corner(ic,j-2),F.BynBx_corner(ic,j-1),F.BynBx_corner(ic,j),F.BynBx_corner(ic,j+1), 1.0/dy);
	double Hxd = stagger4(F.Hfun(ic-1,j),F.Hfun(ic,j),F.Hfun(ic+1,j),F.Hfun(ic+2,j), 1.0/dx);
	double rho_face = face4x(F.nrho,ic,j);
	double Bzn_face = face4x(F.nBz,ic,j);
	double Bxn_face = face4x(F.nBx,ic,j);
	double term1 = Bzn_face/rho_face * (cyd_BynBx - Hxd);

	double BxnBz_xd = stagger4(F.BxnBz_cell(ic-1,j),F.BxnBz_cell(ic,j),F.BxnBz_cell(ic+1,j),F.BxnBz_cell(ic+2,j), 1.0/dx);
	double ByBz_yd  = stagger4(F.ByBz_corner(ic,j-2),F.ByBz_corner(ic,j-1),F.ByBz_corner(ic,j),F.ByBz_corner(ic,j+1), 1.0/dy);
	double term2 = -Bxn_face/rho_face * (BxnBz_xd + ByBz_yd);
	return term1 + term2;
}
inline double FluxBzy4(const MagAX4Fields &F, int i, int jc){
	double BxnBy_xd = stagger4(F.BxnBy(i-2,jc),F.BxnBy(i-1,jc),F.BxnBy(i,jc),F.BxnBy(i+1,jc), 1.0/dx);
	double Eyd = stagger4(F.Efun(i,jc-1),F.Efun(i,jc),F.Efun(i,jc+1),F.Efun(i,jc+2), 1.0/dy);
	double rho_face = face4y(F.nrho,i,jc);
	double Bzn_face = face4y(F.nBz,i,jc);
	double Byn_face = face4y(F.nBy,i,jc);
	double term1 = Bzn_face/rho_face * (BxnBy_xd - Eyd);

	double BxnBz_xd2 = stagger4(F.BxnBz_corner(i-2,jc),F.BxnBz_corner(i-1,jc),F.BxnBz_corner(i,jc),F.BxnBz_corner(i+1,jc), 1.0/dx);
	double BynBz_yd = stagger4(F.BynBz_cell(i,jc-1),F.BynBz_cell(i,jc),F.BynBz_cell(i,jc+1),F.BynBz_cell(i,jc+2), 1.0/dy);
	double term2 = -Byn_face/rho_face * (BxnBz_xd2 + BynBz_yd);
	return term1 + term2;
}

// [removed: Magnetic_Subsystem_AX]

// [removed: Magnetic_Subsystem_potential_AX]

// [removed: Pressure_Subsystem_AX]

// ===================== 嵌套 Newton:保正压力求解 =====================
// 解非线性系统  g(p) = V_mod(p) + T p - b = 0,  其中 T = 焓 Laplacian + 黏性(线性,Qe 冻结)。
// 复用 Pressure_Subsystem_AX:
//   floor_on=false 时它返回线性算子 L p = (1/(gamma-1)) p + T p;
//   floor_on=true  时它返回 Newton Jacobian  J x = V_slope(p_freeze) x + T x。
// 故残差 g = L(p) - (1/(gamma-1)) p + V_mod(p) - b(仅内部 cell 做对角修正;Dirichlet 行保持 p-b)。
// [removed: Pressure_NestedNewton]

// ====================================================================

// [removed: Energy_Subsystem_AX]

void GMRES_Solver(
	Eigen::VectorXd &x,
	Eigen::VectorXd &b,
	std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, Array&, Array&, Array&, Array&, double)> Ax_func,
	Grid &u_con, Grid &Qe, Grid &Qs,
	Array &eps_M_grid_x, Array &eps_M_grid_y, Array &eps_B_grid_x, Array &eps_B_grid_y,
	double h,
	const int dim,
	std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond)
{
	const int JMAX = 25;
	const int outmax = g_gmres_outmax;   // restarts (x JMAX Krylov vectors); a warning is printed if the full budget is exhausted
	g_gmres_conv = true; g_gmres_iters = 0;
	// left preconditioner: applyM(out,in) = M^{-1} in  (identity if no precond)
	auto applyM = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		if(precond) precond(out, in); else omp_copy(out, in);
	};

	// Krylov workspace: allocated (and first touched) once per system size instead of 2 x 26 vectors per call
	// (at 1e6 unknowns that was ~420 MB of allocation and page faults per solve).  Every vector is written before
	// it is read.  GMRES is not re-entered from its operator or preconditioner; a nested call gets its own storage.
	struct GmresWork { std::vector<Eigen::VectorXd> v; Eigen::VectorXd Avj, w, Ax, r0, tmp; };   // A v_j is needed only in step j
	static std::map<int, GmresWork> pool; static int depth = 0;
	GmresWork local; GmresWork& Wk = (depth == 0) ? pool[dim] : local;
	struct DepthGuard { int& d; DepthGuard(int& x) : d(x) { d++; } ~DepthGuard(){ d--; } } depth_guard(depth);
	if((int)Wk.w.size() != dim || (int)Wk.v.size() != JMAX + 1){
		Wk.v.assign(JMAX + 1, Eigen::VectorXd());
		for(int q = 0; q <= JMAX; q++) Wk.v[q].resize(dim);
		Wk.Avj.resize(dim); Wk.w.resize(dim); Wk.Ax.resize(dim); Wk.r0.resize(dim); Wk.tmp.resize(dim);
	}
	std::vector<Eigen::VectorXd>& v = Wk.v; Eigen::VectorXd& Avj = Wk.Avj;
	const Eigen::VectorXd* vp[JMAX + 1]; for(int q = 0; q <= JMAX; q++) vp[q] = &v[q];
	Eigen::VectorXd& w = Wk.w; Eigen::VectorXd& Ax = Wk.Ax; Eigen::VectorXd& r0 = Wk.r0; Eigen::VectorXd& tmp = Wk.tmp;
	
	std::vector<std::vector<double>> H(JMAX + 2, std::vector<double>(JMAX + 2, 0.0));
	std::vector<double> alpha(JMAX + 2, 0.0);
	std::vector<double> gam(JMAX + 2, 0.0);
	std::vector<double> s(JMAX + 2, 0.0);
	std::vector<double> c(JMAX + 2, 0.0);
	
	double beta, pv;
	int j = 0;
	applyM(tmp, b);
	const double bnorm = omp_nrm(tmp);
	// Newton-Jacobian solves of the pressure floor: the Mbig = O(1e9) rows put the attainable relative residual near
	// 1e-8 (see the warning below), so the 1e-12 target is never reached and every such solve ran its full budget of
	// 1250 iterations (measured: blast wave 1000^2, ~4 min per step).  Newton needs only an accurate increment
	// relative to the current residual and tests convergence on its own increment, so these solves stop at 1e-7.
	const double rtol_eff = g_gmres_newton ? fmax(GMRES_rtol, 1e-7) : GMRES_rtol;
	const double GMREStol = bnorm > 0.0 ? rtol_eff*bnorm : 1e-300;   // relative (scale-invariant) stopping test
	
	for (int outer = 0; outer < outmax; outer++) {
		
		omp_zero(Ax);
		Ax_func(Ax, x, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
		omp_sub(tmp, b, Ax);
		applyM(r0, tmp);            // r0 = M^{-1}(b - A x)

		gam[0] = omp_nrm(r0);

		if (gam[0] <= GMREStol) {
//			std::cout << "GMRES converged: outer=" << outer
			return;
		}

		omp_scale(v[0], r0, 1.0/gam[0]);
		omp_zero(Ax);
		Ax_func(Ax, v[0], u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
		applyM(Avj, Ax);            // A v_0 = M^{-1} A v[0]
		
		// Reset H
		for (int ii = 0; ii < JMAX + 2; ii++)
			std::fill(H[ii].begin(), H[ii].end(), 0.0);
		
		for (j = 0; j < JMAX; j++) {
			g_gmres_iters++;
			
			// Arnoldi: classical Gram-Schmidt, all inner products in one pass
			{ double hv[JMAX + 1]; omp_multidot(vp, Avj, j + 1, hv); for (int i = 0; i <= j; i++) H[i][j] = hv[i]; }
			{   // w = A v_j - sum_i H[i][j] v[i]
				double c[JMAX + 1]; for(int i = 0; i <= j; i++) c[i] = H[i][j];
				const int n = (int)w.size();
#pragma omp parallel for schedule(static)
				for(int q = 0; q < n; q++){ double t = Avj[q]; for(int i = 0; i <= j; i++) t -= c[i]*v[i][q]; w[q] = t; }
			}
			H[j + 1][j] = omp_nrm(w);
			
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
				omp_scale(v[j + 1], w, 1.0/H[j + 1][j]);
				omp_zero(Ax);
				Ax_func(Ax, v[j + 1], u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
				applyM(Avj, Ax);         // A v_{j+1} = M^{-1} A v[j+1] (A v_j is no longer needed)
			} else {
				// Converged: back-substitute and update x
				for (int i = j; i >= 0; i--) {
					alpha[i] = gam[i];
					for (int k = i + 1; k <= j; k++) {
						alpha[i] -= H[i][k] * alpha[k];
					}
					alpha[i] /= H[i][i];
				}
				omp_axpy_sum(x, alpha, v, j + 1);
//				std::cout << "GMRES converged: outer=" << outer
//				<< " inner=" << j + 1
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
		omp_axpy_sum(x, alpha, v, JMAX);
		
		if (outer % 10 == 0 || outer < 3) {
//			std::cout << "GMRES restart " << outer
		}
	}
	
	omp_zero(Ax);
	Ax_func(Ax, x, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
	omp_sub(tmp, b, Ax); applyM(r0, tmp);
	// all callers sit inside an outer loop with its own residual test (P-Newton, defect correction, Picard), so a
	// 1e-12 target missed by roundoff (e.g. the Mbig-floored Newton rows, condition ~1e8) is harmless; warn only if large
	const double rr = omp_nrm(r0)/fmax(bnorm,1e-300);
	g_gmres_conv = false;
	if(rr > 1e-8 && !g_gmres_newton && outmax == 50)   // not for the P-Newton Jacobian solves (Mbig rows cap |r|/|b| near 1e-8; Newton tests its own increment); the magnetic solve IS reported
		std::cout << "WARNING: GMRES not converged after " << outmax << " restarts, |r|/|b| = " << rr << std::endl;
}

// [removed: MagSubsystem_UpdateA]

// [removed: PreSubsystem_Pressure]

// [removed: PreSubsystem_Energy]

// [removed: MagSubsystem_UpdateB]

// ============================================================
// Was forward-declared but never defined or called anywhere in the
// baseline file (Document 3 eq.47, step 4: explicit magnetic
// momentum/energy update using the newly-solved B). Called from
// the end of MagSubsystem_UpdateB above.
// ============================================================
// [removed: MagneticFlux_Update_Pressure]

// ||div B||_L2 and ||div B||_inf of the cell-average field with the 4th-order central divergence, the operator
// for which div curl4 A = 0 holds exactly.  For a non-periodic direction the two cells next to the boundary
// (where this stencil would use ghost values) are excluded.  In 1D B_x is constant and both norms are 0.
void divB_norms(Grid &u_con, double &L2, double &Linf){
	L2 = 0.0; Linf = 0.0;
	if(nxCells == 1 || nyCells == 1) return;
	const int i0 = (x_bc == x_Periodic) ? ghost_cells : ghost_cells + 2, i1 = (x_bc == x_Periodic) ? nxCells + ghost_cells : nxCells + ghost_cells - 2;
	const int j0 = (y_bc == y_Periodic) ? ghost_cells : ghost_cells + 2, j1 = (y_bc == y_Periodic) ? nyCells + ghost_cells : nyCells + ghost_cells - 2;
	for(int i = i0; i < i1; i++) for(int j = j0; j < j1; j++){
		const double d = ( -u_con(i+2,j,5) + 8.0*u_con(i+1,j,5) - 8.0*u_con(i-1,j,5) + u_con(i-2,j,5) ) / (12.0*dx)
		               + ( -u_con(i,j+2,6) + 8.0*u_con(i,j+1,6) - 8.0*u_con(i,j-1,6) + u_con(i,j-2,6) ) / (12.0*dy);
		L2 += dx*dy*d*d; Linf = fmax(Linf, fabs(d));
	}
	L2 = sqrt(L2);
}

std::string g_outdir = "";   // output directory, set in main to plot/<test>/
int g_plap = 1;      // [fix] 1 = shock-activated pressure Laplacian of paper 1 (eq. 45/85), (1-eps_M)-weighted, in L2 and L4
std::string g_tag = "";

// intermediate output times (VTK + state file); t = 0 and t = T are always written
std::vector<double> output_times(){
	if(test == testOrszag_Tang) return {0.5/(2.0*M_PI), 2.0/(2.0*M_PI), 3.0/(2.0*M_PI)};
	return {0.25*time_stop, 0.5*time_stop, 0.75*time_stop};
}

std::string test_name(){
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
			testname = "FieldLoop";
			break;
		case testSmooth2D:
			testname = "Smooth2D";
			break;
		case testVortex:
			testname = "Vortex";
			break;
		case testKHInstability:
			testname = "KHInstability";
			break;
		case testdoubleshear:
			testname = "doubleshear";
			break;
	}
	// parameters that distinguish runs of the same test go into the name, so that runs never overwrite each other
	std::ostringstream os; os << testname << "_" << nxCells << "x" << nyCells;
	if(test == testKHInstability) os << "_Mx" << KH_MX;
	if(test == testFieldLoop)     os << "_A" << LOOP_A0;
	if(test == testVortex)        os << "_rho" << VORTEX_RHO0 << "_Bt" << VORTEX_BT << "_U" << VORTEX_U0 << "_L" << VORTEX_L;
	os << "_CFL" << CFL; if(!g_tag.empty()) os << "_" << g_tag;
	testname = os.str();
	return testname;
}

void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri, Array&eps_M,Array&eps_B, const int iteration, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime){
	bool isPlotTime = false, isDumpTime = isFinalTime;
	const std::string testname = test_name();
	
	for(int k=0; k<(int)IntermediatePlotTimes.size(); k++){
		double Tdiff = std::abs(IntermediatePlotTimes[k]-t);
		if(Tdiff <= 1e-10*time_stop){isPlotTime = true; isDumpTime = true; break;}
	}
	if(iteration == 0) isPlotTime = true;
	if(isPlotTime || isFinalTime) printAllVariablesToVTK(u_con,u_pri, eps_M,eps_B,testname,iteration,isFinalTime);
	if(isDumpTime && iteration > 0){   // interior cell averages: i j rho mx my mz E Bx By Bz
		const std::string fn = g_outdir + "state_" + testname + "_" + std::to_string(iteration) + ".dat";
		FILE* fp = fopen(fn.c_str(), "w");
		if(fp){
			fprintf(fp, "# %s  %d x %d  t = %.10e  step = %d\n", testname.c_str(), nxCells, nyCells, t, iteration);
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++) for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				fprintf(fp, "%d %d", i - ghost_cells, j - ghost_cells);
				for(int k = 0; k < NUM_VARS; k++) fprintf(fp, " %.16e", u_con(i,j,k));
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
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
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, Array&eps_M, Array&eps_B,int var, std::string variableName){
	outFile << variableName << " 1 " << nxCells * nyCells << " double" << std::endl;
	
	int printWidth = 24;
	
	for(int jj = 0; jj < nyCells; jj++){
		for(int ii = 0; ii < nxCells; ii++){
			// [vtk] interior cells only: the old loop printed indices 0..N-1, i.e. the low-side GHOST layer
			// (ghost_cells = 4 cells: every plotted field was shifted by 4 cells, and eps_M/eps_B, which
			// are only filled in the interior, showed a ring of zeros at the left and bottom edges).
			const int i = ii + ghost_cells, j = jj + ghost_cells;
			double varPrint;
			if(var < NUM_VARS){
				varPrint = u_con(i,j,var);
			}else if(var < NUM_VARS + 4){
				varPrint = u_pri(i,j,var - NUM_VARS + 1);
			}else if(var < NUM_VARS + 5){
				varPrint = sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6)  + u_con(i,j,7) * u_con(i,j,7) );
			}else if(var < NUM_VARS + 6){
				varPrint = eps_M(i,j);
			}else if(var < NUM_VARS + 7){
				varPrint = eps_B(i,j);
			}
			outFile << std::setw(printWidth) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << varPrint << " ";
		}
		outFile << std::endl;
	}
	
	outFile << std::endl;
}

/**Prints as many variables as desired to a single VTK file.
*/
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, Array&eps_M, Array&eps_B,std::string testname, int nstep, bool isFinalTime){
	std::string fileString = g_outdir + "plt_" + testname + "_" + std::to_string(nstep) + ".vtk";
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
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(13,"eps_M",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(14,"eps_B",0,0,0));

	outFile << "FIELD FieldData " << variablesToPrint.size() << std::endl;
	
	for(unsigned int i = 0;i <variablesToPrint.size(); i++){
		printVTKVariable(outFile,u_con,u_pri,eps_M,eps_B,std::get<0>(variablesToPrint[i]),std::get<1>(variablesToPrint[i]));
	}

	outFile << "POINT_DATA " << (nxCells+1)*(nyCells+1) << std::endl;
	
	return;
}

// ============================================================
// PORTED FROM 2DMHD_EPSILON_HLLEM_WENO3_MP3.cpp (the half-finished file
// verified by the author to produce CORRECT Brio-Wu results): the FUSED
// implicit magnetic+pressure step and its exactly-paired operators.
// Ported as a unit because operator/RHS pairing consistency is precisely
// where the sequential-splitting pipeline kept going wrong.
// ============================================================

// ==================== RESTORED FROM ARCHIVE (positivity core) ====================
// 2nd-order compact M-MATRIX elliptic operator + nested-Newton positivity solver.
// Used as the SOLVER CORE of the defect-correction loop below: every inner solve
// is an M-matrix system (positivity mechanism applies), while the converged
// pressure satisfies the 4th-order composed discretization.
// =================================================================================
// [o3f] CELL-AVERAGE ENTHALPY for the implicit energy flux.  The naive h = g/(g-1) p(Qbar)/rho(Qbar)
// carries the O(h^2) nonlinear-map error p(Qbar) - pbar = -(g-1)(h^2/24)(rho|grad u|^2 + |grad B|^2 + ...),
// which at low Mach is DIVIDED BY rho in h and drives an O(h^2) error of the enthalpy flux
// (measured: low-Mach hydro vortex rho0 = 1e-2, EOC 2.0 for rho and 1.3-2.6 for p while the fully
// explicit path gives 2.9 for both).  Route through point values instead, as the explicit path does:
//   Qbar -> Q_pt (deconvolve) -> h_pt = g/(g-1) p(Q_pt)/rho(Q_pt) -> hbar (convolve),
// with p from Qe and rho from u_con (the operator-family convention of this file).
static Array g_Hbar(nxCells, nyCells);
static long  g_Hbar_stamp = -1;
static void build_Hbar(Grid& Qe, Grid& u_con){
	if(g_Hbar_stamp == g_stage_stamp) return;
	g_Hbar_stamp = g_stage_stamp;
	static Grid Qpt(nxCells, nyCells), Rpt(nxCells, nyCells); static Array Hs(nxCells, nyCells);
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	fv_cell_lim_to(Qpt, Qe, -1.0); fv_cell_lim_to(Rpt, u_con, -1.0);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++){
		double q[NUM_VARS], pr[NUM_VARS]; for(int k = 0; k < NUM_VARS; k++) q[k] = Qpt(i,j,k);
		ini_con2pri(pr, q);
		double hp = gamma_plasma/(gamma_plasma-1.0) * fmax(pr[4], p_floor) / fmax(Rpt(i,j,0), 1e-13);   // density floor = the theta-scaling floor (a 1e-3 floor broke rho0 = 1e-4)
		Hs(i,j) = hp;
	}
	fv_cell_lim_s(g_Hbar, Hs, +1.0);
	bound_cond_array(g_Hbar);
}

// [o3g] CELL AVERAGES of the internal energy and of KE + magnetic energy, via point values.
// the difference (h^2/24)(rho|grad u|^2 + |grad B|^2 + ...) is an O(h^2) pressure offset which, fed to the
// momentum update as grad p, produced an O(h^2) COMPRESSIVE velocity error (Helmholtz decomposition of the
// momentum error on the low-Mach vortex: solenoidal part 3rd order, curl-free part 2nd order).  The explicit
// path already routes W(Q) through point values; this does the same for the implicit energy balance.
// [o3s] avg(f(Q)) = f(Qbar) + [C(f(DQbar)) - f(C(DQbar))]: O(h^4) like C(f(DQbar)), but exact for the
// linear/quadratic bookkeeping at every wavenumber in the sense that the correction vanishes when C∘D = I
// is replaced by its actual discrete action. This removes the Nyquist mismatch (C∘D ≈ 0.89 at k=pi/h)
// between E and avg(KE+MB) that made the Picard sweep diverge at large Alfven CFL / low beta.
static inline double kem_of(const Grid& Q, int i, int j){
	double r = fmax(Q(i,j,0), 1e-300);
	return 0.5*(Q(i,j,1)*Q(i,j,1) + Q(i,j,2)*Q(i,j,2) + Q(i,j,3)*Q(i,j,3))/r
	     + 0.5*(Q(i,j,5)*Q(i,j,5) + Q(i,j,6)*Q(i,j,6) + Q(i,j,7)*Q(i,j,7));
}
static void cellavg_e_kemb(Grid& u, Array& ebar, Array& kembar){
	static Grid Qpt(nxCells, nyCells), Q2(nxCells, nyCells); static Array Ts(nxCells, nyCells), Tc(nxCells, nyCells);
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	fv_cell_lim_to(Qpt, u, -1.0);                 // DQbar
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++) Ts(i,j) = kem_of(Qpt, i, j);
	fv_cell_lim_s(Tc, Ts, +1.0);                  // C(f(DQbar))
	fv_cell_lim_to(Q2, Qpt, +1.0);                // C(DQbar)
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++){
		double kb = kem_of(u, i, j) + Tc(i,j) - kem_of(Q2, i, j);
		kembar(i,j) = kb; ebar(i,j) = u(i,j,4) - kb;
	}
}

// [emb] NEWTON-LINEARISED kinetic + magnetic energy about the explicit stage state Q_E:
//     kl(Q_E, Q) = (m_E.m - 1/2|m_E|^2)/rho + B_E.B - 1/2|B_E|^2     ( = KE + MB when Q_E = Q )
// as a CELL AVERAGE, with the same structure as cellavg_e_kemb: kl(Qbar) + [C kl(D Qbar) - kl(C D Qbar)],
// all conversions unlimited, so the result is affine in (m, B).
static inline double kl_of(const Grid& E, const Grid& Q, int i, int j){
	const double r = fmax(Q(i,j,0), 1e-300);
	const double m1 = E(i,j,1), m2 = E(i,j,2), m3 = E(i,j,3), b5 = E(i,j,5), b6 = E(i,j,6), b7 = E(i,j,7);
	return (m1*Q(i,j,1) + m2*Q(i,j,2) + m3*Q(i,j,3) - 0.5*(m1*m1 + m2*m2 + m3*m3))/r
	     + b5*Q(i,j,5) + b6*Q(i,j,6) + b7*Q(i,j,7) - 0.5*(b5*b5 + b6*b6 + b7*b7);
}
static void kemb_lin(Grid& Q, Grid& Qe, Array& out){
	static Grid Ep(nxCells, nyCells), Qp(nxCells, nyCells), E2(nxCells, nyCells), Q2(nxCells, nyCells); static Array Ts(nxCells, nyCells), Tc(nxCells, nyCells);
	static long ep_stamp = -1;
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	if(ep_stamp != g_stage_stamp){ fv_cell_to(Ep, Qe, -1.0); fv_cell_to(E2, Ep, +1.0); ep_stamp = g_stage_stamp; }
	fv_cell_to(Qp, Q, -1.0, true); fv_cell_to(Q2, Qp, +1.0);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++) Ts(i,j) = kl_of(Ep, Qp, i, j);
	fv_cell_s(Tc, Ts, +1.0);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++) out(i,j) = kl_of(Qe, Q, i, j) + Tc(i,j) - kl_of(E2, Q2, i, j);
}

static void Pressure_Subsystem_AX_ref(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,
                            Array& eps_M_grid_x, Array& eps_M_grid_y, Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	static Grid Qe_pri(nxCells, nyCells);
	static Array shock_eta(nxCells,nyCells);
	static long arch_stamp = -1;
	if(arch_stamp != g_stage_stamp){
		ucon2upri(Qe_pri, Qe);
		ShockDetector(shock_eta, Qe);
		arch_stamp = g_stage_stamp;
	}
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			// 4th-order composed gradient-then-divergence operator with (1-epsM)^2
			// weighting per Document 3's Remark 3.3 principle (confirmed for the
			// magnetic operator earlier): a compact single-formula stencil like
			// eq.82 is only consistent when epsM is spatially uniform. For a
			// genuinely varying epsM, need face-level gradients weighted BEFORE
			// the outer divergence, same structure as Magnetic_Subsystem_AX's
			// FluxBx4. Self-convergence + exact-analytic-target validated
			// separately: order -> 4.000, 4.006 (limited by double-precision
			// roundoff floor beyond that, not a discretization issue).
			// ENTHALPY FLOOR: the explicit HLLEM/WENO3 step can overshoot at
			// strong discontinuities and leave a slightly negative internal
			// energy in a cell (observed: -0.004 at the Brio-Wu compound wave,
			// ts=9 rks=1). A negative enthalpy entering this operator flips the
			// M-matrix sign condition and the elliptic solve amplifies it by
			// ~11 orders of magnitude (the actual Brio-Wu blowup mechanism,
			// pinpointed by stage probes). Clamping at 0 keeps every face
			// coefficient >= 0, preserving the M-matrix property; the tiny
			// negative-pressure cell simply decouples from the implicit
			// diffusion instead of poisoning the whole solve.
			auto henthalpy = [&](int ii,int jj)->double{
				double con_[NUM_VARS];
				for(int k=0;k<NUM_VARS;k++) con_[k]=Qe(ii,jj,k);
				return fmax(calenthalpy(con_) * Qe(ii,jj,0) / u_con(ii,jj,0), 0.0);
			};

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
			// V(p) 对角项:floor_on 时用冻结斜率 V_slope(p_freeze)(此函数即作 Newton Jacobian);
			// 否则退回原线性斜率 1/(gamma-1)。
			double Vdiag = floor_on ? V_slope( p_freeze(M_P(i,j)) ) : (1.0/(gamma_plasma - 1.0));

			// ============================================================
			// LHS elliptic operator: 2nd-order compact M-MATRIX form with
			// (1-epsM)^2 face weighting.
			//
			// WHY 2ND ORDER (deliberate, after the Brio-Wu blowup): the nested
			// Newton positivity proof (V_mod, Lemma 1/Thm 1) requires the
			// operator to be an M-matrix -- all off-diagonal entries <= 0. The
			// 3-point compact form satisfies this (each neighbor coefficient is
			// -h^2/dx^2 * w_face * h_face <= 0 since enthalpy>0, w>=0). ANY
			// 4th-order stencil (eq.82 or the composed double-stagger4) has
			// mixed-sign far coefficients (e.g. the -1/12 in (-1,16,-30,16,-1)/12)
			// and is NOT an M-matrix: positivity fails at strong shocks, which is
			// exactly what killed Brio-Wu at ts~9 (verified by ablation with
			// epsM forced to 0 -- the fully-implicit baseline with the 4th-order
			// operator blows up identically; with this 2nd-order form it doesn't).
			// floor_on "not helping" was the symptom: the theorem's premise was
			// gone, not the mechanism.
			//
			// ACCURACY: this does NOT compromise the 3rd-order target. The
			// elliptic term carries h^2, so h^2 * O(dx^2) = O(dx^4) local error,
			// one order BELOW the explicit WENO3's O(dx^3). The RHS enthalpy
			// divergence (h^1 coefficient) and the momentum gradient stay 4th
			// order where it actually matters.
			//
			// (1-epsM)^2 face weight: energy-row weight x momentum-gradient
			// weight, per the Remark 3.3 composed structure. At epsM=1 the
			// operator row degenerates to Vdiag*pii -- consistent with the
			// (1-epsM)-weighted RHS making the whole implicit step a no-op.
			double h00 = henthalpy(i,j);
			double hface_xp = 0.5*(h00 + henthalpy(i+1,j));
			double hface_xm = 0.5*(h00 + henthalpy(i-1,j));
			double hface_yp = 0.5*(h00 + henthalpy(i,j+1));
			double hface_ym = 0.5*(h00 + henthalpy(i,j-1));
			double wPx_p = (1.0-eps_M_grid_x(i,j))  *(1.0-eps_M_grid_x(i,j));
			double wPx_m = (1.0-eps_M_grid_x(i-1,j))*(1.0-eps_M_grid_x(i-1,j));
			double wPy_p = (1.0-eps_M_grid_y(i,j))  *(1.0-eps_M_grid_y(i,j));
			double wPy_m = (1.0-eps_M_grid_y(i,j-1))*(1.0-eps_M_grid_y(i,j-1));
			double xdiv = ( wPx_p*hface_xp*(getp(i+1,j) - pii) - wPx_m*hface_xm*(pii - getp(i-1,j)) ) / (dx*dx);
			double ydiv = ( wPy_p*hface_yp*(getp(i,j+1) - pii) - wPy_m*hface_ym*(pii - getp(i,j-1)) ) / (dy*dy);

			Ax(M_P(i,j)) = Vdiag * pii - h * h * (xdiv + ydiv);
			{	// [o3q] 2nd-order counterpart of L4's [o3n] term  -h^2 div[-w H p_f d_n eps_f]  (without it the
				// defect correction L2 p_{k+1} = b + (L2 - L4) p_k diverges where eps_M switches within ~1 cell at
				// large implicit acoustic CFL: rho(G) = 1.25 measured in a 1D model).  UPWIND face pressure:
				// p_f = p_R if c_f < 0 else p_L, c_f = w H d eps, keeps every off-diagonal <= 0 (M-matrix) and
				// only adds to the diagonal.
				const double cxp = (1.0-eps_M_grid_x(i,j))  *hface_xp*deps_fx(eps_M_grid_x, i,   j);
				const double cxm = (1.0-eps_M_grid_x(i-1,j))*hface_xm*deps_fx(eps_M_grid_x, i-1, j);
				const double cyp = (1.0-eps_M_grid_y(i,j))  *hface_yp*deps_fy(eps_M_grid_y, i, j  );
				const double cym = (1.0-eps_M_grid_y(i,j-1))*hface_ym*deps_fy(eps_M_grid_y, i, j-1);
				const double pxp = (cxp < 0.0) ? getp(i+1,j) : pii, pxm = (cxm < 0.0) ? pii : getp(i-1,j);
				const double pyp = (cyp < 0.0) ? getp(i,j+1) : pii, pym = (cym < 0.0) ? pii : getp(i,j-1);
				Ax(M_P(i,j)) += h*h*( (cxp*pxp - cxm*pxm)/dx + (cyp*pyp - cym*pym)/dy );
			}
			
			if(g_emb_L2c && g_emb_Qe){
				// [emb2] the linearised kinetic energy  kl = (m_E.m - |m_E|^2/2)/rho  with  m = m** - h G p  puts the
				// convective term  -h u_E . G p  into the embedded operator.  Its M-matrix discretisation: one-sided
				// differences on the side that makes every off-diagonal <= 0 and only adds to the diagonal
				//     -h u (p_{i+1} - p_i)/dx  (u > 0),     -h u (p_i - p_{i-1})/dx  (u < 0).
				const double ue = (*g_emb_Qe)(i,j,1)/u_con(i,j,0), ve = (*g_emb_Qe)(i,j,2)/u_con(i,j,0);
				const double wx = 0.5*((1.0-eps_M_grid_x(i,j)) + (1.0-eps_M_grid_x(i-1,j)));
				const double wy = 0.5*((1.0-eps_M_grid_y(i,j)) + (1.0-eps_M_grid_y(i,j-1)));
				Ax(M_P(i,j)) -= h*( wx*( fmax(ue,0.0)*(getp(i+1,j) - pii) + fmin(ue,0.0)*(pii - getp(i-1,j)) )/dx
				                  + wy*( fmax(ve,0.0)*(getp(i,j+1) - pii) + fmin(ve,0.0)*(pii - getp(i,j-1)) )/dy );
			}

			//Laplacian Operator
			double nu_x, nu_y;
			double nu_x_p, nu_y_p;
			double nu_x_m, nu_y_m;
			
			double lap_p;
			{	// shock-activated pressure Laplacian, chi-gated at the faces, kappa = |u|_max dx (paper: chi |u| dx)
				nu_x = vx_max * dx;   // [o3s] kappa = chi |u|_max dx (paper eq. 84; dimensionally L^2/T)
				nu_y = vy_max * dy;
				lap_p =  nu_x/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * (getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * (pii - getp(i - 1,j))) + nu_y/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * (pii - getp(i,j-1)));
			}
			if(g_plap){   // [fix] paper-1 eq.(45)/(85): V'(p) * dt * div(chi |u|_max dx grad p), (1-eps_M) per face, compact (M-matrix)
				const double a_ = 1.0/(gamma_plasma - 1.0);
				const double gxp = 0.5*(shock_eta(i+1,j)+shock_eta(i,j))*(1.0-eps_M_grid_x(i,j)),   gxm = 0.5*(shock_eta(i,j)+shock_eta(i-1,j))*(1.0-eps_M_grid_x(i-1,j));
				const double gyp = 0.5*(shock_eta(i,j+1)+shock_eta(i,j))*(1.0-eps_M_grid_y(i,j)),   gym = 0.5*(shock_eta(i,j)+shock_eta(i,j-1))*(1.0-eps_M_grid_y(i,j-1));
				const double lp = nu_x/(dx*dx)*(gxp*(getp(i+1,j)-pii) - gxm*(pii-getp(i-1,j))) + nu_y/(dy*dy)*(gyp*(getp(i,j+1)-pii) - gym*(pii-getp(i,j-1)));
				Ax(M_P(i,j)) -= a_*h*lp;
			}
		}
	}
}
// ---------------------------------------------------------------------------------------------------------
// The inner pressure operator above (L2 / L2c of the defect correction, plus the Newton diagonal) is a 5-point
// stencil whose coefficients are fixed within a stage: they depend on Q_E, the density, eps_M and h, and on whether
// the L2c convective term is on.  It is assembled once per stage and mode and applied as
//     (A x)_q = (V'_q + c0_q) x_q + cE_q x_E + cW_q x_W + cN_q x_N + cS_q x_S ,   V'_q = 1/(gamma-1) or V_slope(p_freeze_q),
// instead of recomputing the cell enthalpy (a conservative -> primitive conversion) five times per cell and product.
// Each assembly is checked against the matrix-free form on a random vector; on a mismatch (not expected) the
// matrix-free form is used for that stage.  The stencil also gives the diagonal for the Jacobi preconditioner and
// the Gershgorin radius that decides whether that preconditioner is enough.
// ---------------------------------------------------------------------------------------------------------
struct PStencil {
	long stamp = -1; double h = 0.0; const void* key[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
	bool ok = false;
	double gersh = 0.0;   // max_q sum_nb |c_nb,q| / |1/(gamma-1) + c0_q|, Gershgorin radius of the Jacobi-scaled operator
	std::vector<double> c0, cE, cW, cN, cS;
	std::vector<int> iE, iW, iN, iS;
	std::vector<char> id;   // identity row (strong Dirichlet boundary cell)
};
static PStencil g_pst[2];   // [g_emb_L2c]
static bool g_pst_warned = false;

// neighbour (ii, jj) -> unknown index and Dirichlet mirror sign, as getp in the matrix-free operator
static inline void pst_nb(int ii, int jj, int& idx, double& sg){
	int II = ii, JJ = jj; sg = 1.0;
	if(x_bc == x_Dirichlet){
		if(ii < ghost_cells){ II = 2*ghost_cells - 1 - ii; sg = -sg; }
		else if(ii >= nxCells + ghost_cells){ II = 2*(nxCells+ghost_cells) - 1 - ii; sg = -sg; }
	}
	if(y_bc == y_Dirichlet){
		if(jj < ghost_cells){ JJ = 2*ghost_cells - 1 - jj; sg = -sg; }
		else if(jj >= nyCells + ghost_cells){ JJ = 2*(nyCells+ghost_cells) - 1 - jj; sg = -sg; }
	}
	idx = M_P(II, JJ);
}

static void pst_apply(const PStencil& S, const Eigen::VectorXd& x, Eigen::VectorXd& Ax){
	const int n = nxCells*nyCells; const double a = 1.0/(gamma_plasma - 1.0);
	const bool fl = floor_on;
	const double* xp = x.data(); double* yp = Ax.data(); const double* pf = fl ? p_freeze.data() : nullptr;
#pragma omp parallel for schedule(static)
	for(int q = 0; q < n; q++){
		if(S.id[q]){ yp[q] = xp[q]; continue; }
		const double vd = fl ? V_slope(pf[q]) : a;
		yp[q] = (vd + S.c0[q])*xp[q] + S.cE[q]*xp[S.iE[q]] + S.cW[q]*xp[S.iW[q]] + S.cN[q]*xp[S.iN[q]] + S.cS[q]*xp[S.iS[q]];
	}
}

static void pst_assemble(PStencil& S, Grid& u_con, Grid& Qe, Array& eMx, Array& eMy, double h){
	const int n = nxCells*nyCells;
	if((int)S.c0.size() != n){
		S.c0.resize(n); S.cE.resize(n); S.cW.resize(n); S.cN.resize(n); S.cS.resize(n);
		S.iE.resize(n); S.iW.resize(n); S.iN.resize(n); S.iS.resize(n); S.id.resize(n);
	}
	static Array Hc(nxCells, nyCells);   // cell enthalpy (henthalpy of the matrix-free operator), interior + one ghost layer
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells - 1; j <= nyCells + ghost_cells; j++)
		for(int i = ghost_cells - 1; i <= nxCells + ghost_cells; i++){
			double con_[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++) con_[k] = Qe(i,j,k);
			Hc(i,j) = fmax(calenthalpy(con_) * Qe(i,j,0) / u_con(i,j,0), 0.0);
		}
	const bool emb = g_emb_L2c && g_emb_Qe;
	const double a = 1.0/(gamma_plasma - 1.0), hh = h*h;
	static Array pst_eta(nxCells, nyCells); static long pst_eta_stamp = -1;
	if(g_plap && pst_eta_stamp != g_stage_stamp){ ShockDetector(pst_eta, Qe); pst_eta_stamp = g_stage_stamp; }
	double gmax = 0.0;
#pragma omp parallel for schedule(static) reduction(max:gmax)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			const int q = M_P(i,j);
			if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ){
				S.id[q] = 1; S.c0[q] = 0.0; S.cE[q] = S.cW[q] = S.cN[q] = S.cS[q] = 0.0; S.iE[q] = S.iW[q] = S.iN[q] = S.iS[q] = q;
				continue;
			}
			S.id[q] = 0;
			double c0 = 0.0, cc[4] = {0.0, 0.0, 0.0, 0.0};   // E, W, N, S
			const double h00 = Hc(i,j);
			const double hfx_p = 0.5*(h00 + Hc(i+1,j)), hfx_m = 0.5*(h00 + Hc(i-1,j));
			const double hfy_p = 0.5*(h00 + Hc(i,j+1)), hfy_m = 0.5*(h00 + Hc(i,j-1));
			// - h^2 (xdiv + ydiv), (1-eps_M)^2 face weights
			const double wPx_p = (1.0-eMx(i,j))*(1.0-eMx(i,j)), wPx_m = (1.0-eMx(i-1,j))*(1.0-eMx(i-1,j));
			const double wPy_p = (1.0-eMy(i,j))*(1.0-eMy(i,j)), wPy_m = (1.0-eMy(i,j-1))*(1.0-eMy(i,j-1));
			const double kxp = hh*wPx_p*hfx_p/(dx*dx), kxm = hh*wPx_m*hfx_m/(dx*dx);
			const double kyp = hh*wPy_p*hfy_p/(dy*dy), kym = hh*wPy_m*hfy_m/(dy*dy);
			c0 += kxp + kxm + kyp + kym; cc[0] -= kxp; cc[1] -= kxm; cc[2] -= kyp; cc[3] -= kym;
			// [o3q] + h^2 div[c p_f] with the upwind face pressure p_f
			const double cxp = (1.0-eMx(i,j))  *hfx_p*deps_fx(eMx, i,   j), cxm = (1.0-eMx(i-1,j))*hfx_m*deps_fx(eMx, i-1, j);
			const double cyp = (1.0-eMy(i,j))  *hfy_p*deps_fy(eMy, i, j  ), cym = (1.0-eMy(i,j-1))*hfy_m*deps_fy(eMy, i, j-1);
			if(cxp < 0.0) cc[0] += hh*cxp/dx; else c0 += hh*cxp/dx;
			if(cxm < 0.0) c0 -= hh*cxm/dx;    else cc[1] -= hh*cxm/dx;
			if(cyp < 0.0) cc[2] += hh*cyp/dy; else c0 += hh*cyp/dy;
			if(cym < 0.0) c0 -= hh*cym/dy;    else cc[3] -= hh*cym/dy;
			if(emb){   // [emb2] - h u_E . G p, one-sided on the M-matrix side
				const double ue = (*g_emb_Qe)(i,j,1)/u_con(i,j,0), ve = (*g_emb_Qe)(i,j,2)/u_con(i,j,0);
				const double wx = 0.5*((1.0-eMx(i,j)) + (1.0-eMx(i-1,j))), wy = 0.5*((1.0-eMy(i,j)) + (1.0-eMy(i,j-1)));
				const double ex = h*wx/dx, ey = h*wy/dy;
				const double up = fmax(ue,0.0), um = fmin(ue,0.0), vp = fmax(ve,0.0), vm = fmin(ve,0.0);
				cc[0] -= ex*up; c0 += ex*(up - um); cc[1] += ex*um;
				cc[2] -= ey*vp; c0 += ey*(vp - vm); cc[3] += ey*vm;
			}
			if(g_plap){   // [fix] shock-activated pressure Laplacian (same as the matrix-free form)
				const double nu_x = vx_max*dx, nu_y = vy_max*dy;
				const double gxp = 0.5*(pst_eta(i+1,j)+pst_eta(i,j))*(1.0-eMx(i,j)),   gxm = 0.5*(pst_eta(i,j)+pst_eta(i-1,j))*(1.0-eMx(i-1,j));
				const double gyp = 0.5*(pst_eta(i,j+1)+pst_eta(i,j))*(1.0-eMy(i,j)),   gym = 0.5*(pst_eta(i,j)+pst_eta(i,j-1))*(1.0-eMy(i,j-1));
				const double kxp2 = a*h*nu_x/(dx*dx)*gxp, kxm2 = a*h*nu_x/(dx*dx)*gxm, kyp2 = a*h*nu_y/(dy*dy)*gyp, kym2 = a*h*nu_y/(dy*dy)*gym;
				c0 += kxp2 + kxm2 + kyp2 + kym2; cc[0] -= kxp2; cc[1] -= kxm2; cc[2] -= kyp2; cc[3] -= kym2;
			}
			int nb[4]; double sg[4];
			pst_nb(i+1, j, nb[0], sg[0]); pst_nb(i-1, j, nb[1], sg[1]); pst_nb(i, j+1, nb[2], sg[2]); pst_nb(i, j-1, nb[3], sg[3]);
			for(int m = 0; m < 4; m++){ cc[m] *= sg[m]; if(nb[m] == q){ c0 += cc[m]; cc[m] = 0.0; } }   // a neighbour mapped onto the cell itself
			S.c0[q] = c0; S.cE[q] = cc[0]; S.cW[q] = cc[1]; S.cN[q] = cc[2]; S.cS[q] = cc[3];
			S.iE[q] = nb[0]; S.iW[q] = nb[1]; S.iN[q] = nb[2]; S.iS[q] = nb[3];
			gmax = fmax(gmax, (fabs(cc[0]) + fabs(cc[1]) + fabs(cc[2]) + fabs(cc[3])) / fabs(a + c0));
		}
	S.gersh = gmax;
}

static PStencil& pst_get(Grid& u_con, Grid& Qe, Grid& Qs, Array& eMx, Array& eMy, Array& eBx, Array& eBy, double h){
	PStencil& S = g_pst[g_emb_L2c ? 1 : 0];
	const void* key[5] = {&u_con, &Qe, &eMx, &eMy, g_emb_L2c ? (const void*)g_emb_Qe : nullptr};
	bool same = (S.stamp == g_stage_stamp && S.h == h);
	for(int q = 0; q < 5; q++) same = same && (S.key[q] == key[q]);
	if(same) return S;
	pst_assemble(S, u_con, Qe, eMx, eMy, h);
	S.stamp = g_stage_stamp; S.h = h; for(int q = 0; q < 5; q++) S.key[q] = key[q];
	{	// check against the matrix-free operator (linear part)
		const int n = nxCells*nyCells;
		static Eigen::VectorXd z, y1, y2; z.resize(n); y1.resize(n); y2.resize(n);
#pragma omp parallel for schedule(static)
		for(int q = 0; q < n; q++) z[q] = (double)(((unsigned)q * 2654435761u) % 1000003u) / 500001.5 - 1.0;
		const bool sv = floor_on; floor_on = false;
		Pressure_Subsystem_AX_ref(y1, z, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
		pst_apply(S, z, y2);
		floor_on = sv;
		const double err = omp_dist(y1, y2), ref = omp_nrm(y1);
		S.ok = (err <= 1e-12*ref);
		if(!S.ok && !g_pst_warned){ std::cout << "WARNING: pressure stencil differs from the matrix-free operator (" << err/ref << "); using the matrix-free form" << std::endl; g_pst_warned = true; }
	}
	return S;
}

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,
                            Array& eps_M_grid_x, Array& eps_M_grid_y, Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	PStencil& S = pst_get(u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
	if(S.ok) pst_apply(S, x, Ax);
	else     Pressure_Subsystem_AX_ref(Ax, x, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
}
// ---------------------------------------------------------------------------------------------------------
// Preconditioner of the pressure GMRES.  The inner pressure operator (the compact second-order M-matrix of the
// defect correction, with the Newton diagonal) has a short stencil; it is assembled as a sparse matrix by probing
// (one matrix-free product per colour, colours separated by more than the stencil width) and factorised by sparse
// LU.  The assembled matrix is checked against the matrix-free operator on a random vector; if they disagree the
// preconditioner is switched off for that operator.  It is rebuilt once per stage and operator mode; within a stage
// it changes only through the Newton diagonal where the pressure floor is active, where it remains a preconditioner.
// The preconditioner changes only the convergence of GMRES, not the system that is solved.
// ---------------------------------------------------------------------------------------------------------
struct PressurePrecond {
	long stamp = -1; int mode = -1; bool ok = false;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
};
static PressurePrecond g_pprec[4];
static bool g_pprec_warned = false;

static void pprec_build(PressurePrecond& P, Grid& u_con, Grid& Qe, Grid& Qs, Array& eMx, Array& eMy, Array& eBx, Array& eBy, double h){
	{	// the current operator (Newton diagonal included) directly from the assembled stencil
		PStencil& S = pst_get(u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
		if(S.ok){
			const int n = nxCells*nyCells; const double a = 1.0/(gamma_plasma - 1.0); const bool fl = floor_on;
			std::vector<Eigen::Triplet<double>> T; T.reserve((size_t)n * 5);
			for(int q = 0; q < n; q++){
				if(S.id[q]){ T.emplace_back(q, q, 1.0); continue; }
				T.emplace_back(q, q, (fl ? V_slope(p_freeze[q]) : a) + S.c0[q]);
				if(S.cE[q] != 0.0) T.emplace_back(q, S.iE[q], S.cE[q]);
				if(S.cW[q] != 0.0) T.emplace_back(q, S.iW[q], S.cW[q]);
				if(S.cN[q] != 0.0) T.emplace_back(q, S.iN[q], S.cN[q]);
				if(S.cS[q] != 0.0) T.emplace_back(q, S.iS[q], S.cS[q]);
			}
			P.A.resize(n, n); P.A.setFromTriplets(T.begin(), T.end());
			P.lu.analyzePattern(P.A); P.lu.factorize(P.A); P.ok = (P.lu.info() == Eigen::Success);
			if(!P.ok && !g_pprec_warned){ std::cout << "WARNING: pressure preconditioner disabled (sparse LU failed)" << std::endl; g_pprec_warned = true; }
			return;
		}
	}
	// matrix-free operator only: assemble by probing
	const int nx = nxCells, ny = nyCells, n = nx*ny, r = 2;   // r: stencil radius covered by the probing
	auto period = [&](int N, bool per)->int{ if(N <= 2*r + 1) return N; if(!per) return 2*r + 1; for(int s = 2*r + 1; s <= N; s++) if(N % s == 0) return s; return N; };
	const bool px = (x_bc == x_Periodic), py = (y_bc == y_Periodic);
	const int sx = period(nx, px), sy = period(ny, py);
	std::vector<Eigen::Triplet<double>> T; T.reserve((size_t)n * 25);
	Eigen::VectorXd e(n), y(n);
	for(int cx = 0; cx < sx; cx++) for(int cy = 0; cy < sy; cy++){
		for(int J = 0; J < ny; J++) for(int I = 0; I < nx; I++) e(I + nx*J) = (I % sx == cx && J % sy == cy) ? 1.0 : 0.0;
		y.setZero();
		Pressure_Subsystem_AX(y, e, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
		for(int J = 0; J < ny; J++) for(int I = 0; I < nx; I++){
			const double v = y(I + nx*J); if(v == 0.0) continue;
			int ci = -1, cj = -1;   // the unique column of this colour within the stencil window of row (I,J)
			if(sx == nx) ci = cx; else for(int d = -r; d <= r; d++){ int c = I + d; if(px) c = (c % nx + nx) % nx; else if(c < 0 || c >= nx) continue; if(c % sx == cx){ ci = c; break; } }
			if(sy == ny) cj = cy; else for(int d = -r; d <= r; d++){ int c = J + d; if(py) c = (c % ny + ny) % ny; else if(c < 0 || c >= ny) continue; if(c % sy == cy){ cj = c; break; } }
			if(ci >= 0 && cj >= 0) T.emplace_back(I + nx*J, ci + nx*cj, v);
		}
	}
	P.A.resize(n, n); P.A.setFromTriplets(T.begin(), T.end());
	// check the assembled matrix against the matrix-free operator
	std::mt19937 gen(12345); std::uniform_real_distribution<double> U(-1.0, 1.0);
	Eigen::VectorXd z(n), Az(n); for(int q = 0; q < n; q++) z(q) = U(gen);
	Az.setZero(); Pressure_Subsystem_AX(Az, z, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
	const double err = (P.A * z - Az).norm() / (Az.norm() + 1e-300);
	P.ok = (err < 1e-10);
	if(P.ok){ P.lu.analyzePattern(P.A); P.lu.factorize(P.A); P.ok = (P.lu.info() == Eigen::Success); }
	if(!P.ok && !g_pprec_warned){ std::cout << "WARNING: pressure preconditioner disabled (assembly check " << err << ")" << std::endl; g_pprec_warned = true; }
}

// Pressure GMRES.
//  * diagonally dominant operator (Gershgorin radius of the Jacobi-scaled operator below 1/2, i.e. a small implicit
//    acoustic CFL or eps_M -> 1, e.g. the blast wave and Orszag-Tang): GMRES with the diagonal as preconditioner, which
//    converges at least like 2^-k and takes the pressure-floor rows (diagonal M ~ 1e9) exactly.  No factorisation;
//    each iteration is one parallel stencil product.  (Before, the floor rows defeated both the unpreconditioned
//    solve and any factorisation from an earlier iterate, and the blast wave refactorised every stage.)
//  * otherwise, the sparse-LU preconditioner built only when it is needed and kept while it works:
//  * no factorisation yet: solve with the diagonal preconditioner, budget 50 Krylov iterations; if that fails, factorise.
//  * a factorisation from an earlier stage or step exists: use it (the operator changes slowly from stage to
//    stage); if the preconditioned solve does not converge within 100 iterations, refactorise (a serial
//    factorisation at 1e6 unknowns costs ~30 s, more than 100 preconditioned iterations but far less than a stall)
//    and continue from the current iterate.
// Well-conditioned systems (Mach O(1)) never factorise; stiff ones (low Mach, low beta) factorise rarely instead
// of once per stage -- the sparse LU is serial and at 1e6 unknowns dominated the step.  Every path solves the same
// system to the same tolerance; only the iteration count changes.
static void pressure_gmres(Eigen::VectorXd& x, Eigen::VectorXd& b, Grid& u_con, Grid& Qe, Grid& Qs,
                           Array& eMx, Array& eMy, Array& eBx, Array& eBy, double h, const int dim){
	PStencil& S = pst_get(u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
	static Eigen::VectorXd dinv;
	auto jacobi = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		const int n = (int)in.size();
#pragma omp parallel for schedule(static)
		for(int q = 0; q < n; q++) out[q] = in[q]*dinv[q];
	};
	auto solve_jacobi = [&](){   // the diagonal of the current operator (Newton diagonal included)
		dinv.resize(dim); const double a = 1.0/(gamma_plasma - 1.0); const bool fl = floor_on;
#pragma omp parallel for schedule(static)
		for(int q = 0; q < dim; q++) dinv[q] = S.id[q] ? 1.0 : 1.0/((fl ? V_slope(p_freeze[q]) : a) + S.c0[q]);
		GMRES_Solver(x, b, Pressure_Subsystem_AX, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h, dim, jacobi);
	};
	const bool diag_dominant = S.ok && S.gersh < 0.5;
	if(diag_dominant){
		g_gmres_outmax = 4;
		solve_jacobi();
		g_gmres_outmax = 50;
		if(g_gmres_conv) return;
	}
	const int mode = (floor_on ? 1 : 0) + (g_emb_L2c ? 2 : 0);
	PressurePrecond& P = g_pprec[mode];
	auto build = [&](){ const auto t0 = std::chrono::steady_clock::now(); pprec_build(P, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h);
	                    P.stamp = g_stage_stamp; g_lu_builds++; g_lu_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
	auto solve = [&](){ if(P.ok) GMRES_Solver(x, b, Pressure_Subsystem_AX, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h, dim,
	                                          [&P](Eigen::VectorXd& out, const Eigen::VectorXd& in){ out = P.lu.solve(in); });
	                    else if(S.ok) solve_jacobi();
	                    else     GMRES_Solver(x, b, Pressure_Subsystem_AX, u_con, Qe, Qs, eMx, eMy, eBx, eBy, h, dim); };
	if(!P.ok){
		if(P.stamp == g_stage_stamp){ solve(); return; }   // factorisation failed this stage: no LU, full budget
		if(!diag_dominant){
			g_gmres_outmax = 2;
			solve();
			g_gmres_outmax = 50;
			if(g_gmres_conv) return;
		}
		build(); solve(); return;
	}
	// an existing factorisation (earlier stage, step or Newton iterate): try it with 100 iterations; if that fails,
	// refactorise at the current operator.  Within a Newton solve the Jacobian changes where the pressure floor
	// switches the diagonal between 1/(gamma-1) and Mbig, which can make a factorisation from the previous iterate
	// useless; restarted GMRES on that operator then stalls for its full budget, far longer than one refactorisation.
	g_gmres_outmax = 4;
	solve();
	g_gmres_outmax = 50;
	if(g_gmres_conv) return;
	build();
	solve();
}

void Pressure_NestedNewton(Eigen::VectorXd& p, Eigen::VectorXd& b,
                           Grid& u_con, Grid& Qe, Grid& Qs,
                           Array& eps_M_grid_x, Array& eps_M_grid_y, double h, const int dim){
	// eps_B isn't used by Pressure_Subsystem_AX (only pressure/epsM couples here) --
	// dummy zero-filled placeholder just to satisfy GMRES_Solver's uniform Ax_func
	// signature, shared with Magnetic_Subsystem_AX which does need eps_B.
	static Array dummy_eps_B_x(nxCells,nyCells), dummy_eps_B_y(nxCells,nyCells);

	if(!floor_on){                       // 关掉保正 -> 完全退回原来的单次线性 GMRES
		pressure_gmres(p, b, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, dummy_eps_B_x, dummy_eps_B_y, h, dim);
		return;
	}
	const double a = 1.0/(gamma_plasma - 1.0);
	{	// Cor. 1 of 2608.15837: M chosen from the right-hand side so that min p >= theta p_floor a priori
		const double bmin = omp_min(b), Vpe = p_floor * a;
		const double M_cor = std::max(0.0, Vpe - bmin) / ((1.0 - theta_floor) * p_floor);
		Mbig_cur = std::max(Mbig_fac * a, M_cor);
	}
	static Eigen::VectorXd g, dp, neg_g, Lp; g.resize(dim); dp.resize(dim); neg_g.resize(dim); Lp.resize(dim);
	p_freeze = p;                        // Newton 初值 = 当前压力
	int it = 0;
	for(it = 0; it < Pnewton_max; it++){
		// ---- 残差 g(p_freeze) ----
		bool save = floor_on; floor_on = false;                 // 让 AX 给线性算子 L p
		Pressure_Subsystem_AX(Lp, p_freeze, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, dummy_eps_B_x, dummy_eps_B_y, h);
		floor_on = save;
#pragma omp parallel for collapse(2) schedule(static)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				int idx = M_P(i,j);
				bool bnd = ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1))
				         || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1));
				if(bnd) g(idx) = Lp(idx) - b(idx);                          // 恒等行: p - p_amb
				else    g(idx) = Lp(idx) - a*p_freeze(idx) + V_mod(p_freeze(idx)) - b(idx);
			}
		}
		double res = omp_nrm(g);
		// RELATIVE criterion: with Mbig=1e8 the floored rows carry O(1e8*p) diagonal
		// entries, so an ABSOLUTE residual of 1e-10 is below double precision for
		// them -- Newton "hits max iters" while the solution has long since
		// converged to the floor. Scale by ||b|| (+1 to avoid the zero-b corner).
		double res_rel = res / (omp_nrm(b) + 1.0);
		if((res < Pnewton_tol || res_rel < 1e-12) && omp_min(p_freeze) >= 0.0){   // [o3q] the iterate, not the input guess
			break;
		}
		// ---- Jacobian 解:J dp = -g(floor_on=true,AX 即 Jacobian) ----
		omp_scale(neg_g, g, -1.0);
		omp_zero(dp);
		g_gmres_newton = true;
		pressure_gmres(dp, neg_g, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, dummy_eps_B_x, dummy_eps_B_y, h, dim);
		g_gmres_newton = false;
		omp_axpy(p_freeze, 1.0, dp);
		// INCREMENT criterion: the Mbig=1e8 floored rows carry a floating-point
		// residual FLOOR ~1e-7 (1e8 * eps_mach) that no residual test can reach.
		// A tiny Newton increment IS convergence, immune to that noise floor.
		if(omp_nrm(dp) <= 1e-12 * (omp_nrm(p_freeze) + 1.0) && omp_min(p_freeze) >= 0.0){ it = -1; break; }   // converged (flag)
	}
	if(it == Pnewton_max)   // only a genuine non-convergence (converged exits break with it < Pnewton_max or it = -1)
		std::cout << "  [P-Newton] WARNING: hit max iters (" << Pnewton_max << ")" << std::endl;
	p = p_freeze;
}

void Pressure_Subsystem_AX_HF(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,Array& eps_M_grid_x, Array& eps_M_grid_y,Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	build_Hbar(Qe, u_con);   // [o3f] (no-op if already built for this stage)
	static Array hf_eta(nxCells, nyCells); static long hf_eta_stamp = -1;
	if(g_plap && hf_eta_stamp != g_stage_stamp){ ShockDetector(hf_eta, Qe); hf_eta_stamp = g_stage_stamp; }
//			
//			
//			//Laplacian Operator
//			
//			
//			
//			switch (flux_diss) {
//				case Global_Diss:
//					switch(add_pvisc){
//						case MaterialSpeedp:
//							break;
//						case AlfvenSpeedp:
//							break;
//						case Noaddvisp:
//							break;
//					break;
//				case Local_Diss:
//					switch(add_pvisc){
//						case MaterialSpeedp:
//							
//							break;
//						case AlfvenSpeedp:
//							
//							break;
//						case Noaddvisp:
//							break;
//					break;
//					
	// ============================================================
	// 4TH-ORDER COMPOSED ELLIPTIC OPERATOR (user decision: positivity
	// M-matrix constraint dropped in favor of 4th order).
	//   face flux:  (1-epsM)_face^2 * face4(henthalpy) * G4(p)
	//     G4 = 4th-order staggered face gradient (stagger4 of cell p)
	//   cell value: stagger4 divergence of the four face fluxes.
	// PAIRING (the load-bearing property, learned the hard way): the
	// SAME face4/stagger4 family and the SAME enthalpy lambda are used
	// by (a) this LHS, (b) the RHS enthalpy divergence, and (c) the
	// momentum update below in PreMagSubsystem_Pressure_HF. epsM
	// weights are taken directly at faces (never interpolated).
	// CONSEQUENCE: this operator is NOT an M-matrix (mixed-sign far
	// coefficients) -> the nested-Newton positivity theorem premise is
	// gone. Positivity must come from elsewhere if needed later
	// (e.g. explicit-side limiting). See archived-nested-newton.cpp.
	// ============================================================
	auto henthA = [&](int ii,int jj)->double{
		// POSITIVITY GUARD on the operator coefficient: the RK stage combination Qe
		// carries no positivity guarantee (non-SSP DIRK weights), and a negative
		// p_Qe here makes the elliptic coefficient NEGATIVE -> anti-diffusion,
		// ill-posed operator, Newton/GMRES stall (observed: hit-max-iters every
		// solve on the 256^2 blastwave). Clamping to p_floor keeps the operator
		// well-posed; applied identically in LHS and RHS (operator-pairing law).
		return g_Hbar(ii,jj);   // [o3f] cell-average enthalpy (p_Qe / rho_ucon, both floored, via point values)
	};
	auto getpA = [&](int ii,int jj)->double{ return x(M_P(ii,jj)); };
	// [o3e] FV family: face-averaged enthalpy x face gradient (+ transverse covariance), two-point divergence
	auto Gx = [&](int ic,int jc)->double{ return gradfvx(getpA, ic, jc, 1.0/dx); };
	auto Gy = [&](int ic,int jc)->double{ return gradfvy(getpA, ic, jc, 1.0/dy); };
	auto FPx4 = [&](int ic, int j)->double{
		double w = 1.0 - eps_M_grid_x(ic,j);
		double prod = facefvx(henthA,ic,j) * Gx(ic,j);
		// [o3d] covariance of the two face averages along the face (y)
		prod += (facefvx(henthA,ic,j+1) - facefvx(henthA,ic,j-1)) * (Gx(ic,j+1) - Gx(ic,j-1)) / 48.0;
		// [o3n] face momentum of the conservative update: -(1-eps) G p + p_f d_x eps_f
		return w*w * prod - w * facefvx(henthA,ic,j) * facefvx(getpA,ic,j) * deps_fx(eps_M_grid_x, ic, j);
	};
	auto FPy4 = [&](int i, int jc)->double{
		double w = 1.0 - eps_M_grid_y(i,jc);
		double prod = facefvy(henthA,i,jc) * Gy(i,jc);
		prod += (facefvy(henthA,i+1,jc) - facefvy(henthA,i-1,jc)) * (Gy(i+1,jc) - Gy(i-1,jc)) / 48.0;
		return w*w * prod - w * facefvy(henthA,i,jc) * facefvy(getpA,i,jc) * deps_fy(eps_M_grid_y, i, jc);   // [o3n]
	};
	// face fluxes, each evaluated once: FXf(i,j) = face i+1/2 of row j, FYf(i,j) = face j+1/2 of column i
	static Array FXf(nxCells, nyCells), FYf(nxCells, nyCells);
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
			if(j >= ghost_cells) FXf(i,j) = FPx4(i, j);
			if(i >= ghost_cells) FYf(i,j) = FPy4(i, j);
		}
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double xdivA = (FXf(i,j) - FXf(i-1,j)) / dx;
			double ydivA = (FYf(i,j) - FYf(i,j-1)) / dy;
			Ax(M_P(i,j)) = 1/(gamma_plasma - 1) * x(M_P(i,j)) - h * h * (xdivA + ydivA);
			if(g_plap){   // [fix] shock-activated pressure Laplacian, identical to the inner operator's
				const double a_ = 1.0/(gamma_plasma - 1.0), nu_x = vx_max*dx, nu_y = vy_max*dy, pii = x(M_P(i,j));
				const double gxp = 0.5*(hf_eta(i+1,j)+hf_eta(i,j))*(1.0-eps_M_grid_x(i,j)),   gxm = 0.5*(hf_eta(i,j)+hf_eta(i-1,j))*(1.0-eps_M_grid_x(i-1,j));
				const double gyp = 0.5*(hf_eta(i,j+1)+hf_eta(i,j))*(1.0-eps_M_grid_y(i,j)),   gym = 0.5*(hf_eta(i,j)+hf_eta(i,j-1))*(1.0-eps_M_grid_y(i,j-1));
				const double lp = nu_x/(dx*dx)*(gxp*(x(M_P(i+1,j))-pii) - gxm*(pii-x(M_P(i-1,j)))) + nu_y/(dy*dy)*(gyp*(x(M_P(i,j+1))-pii) - gym*(pii-x(M_P(i,j-1))));
				Ax(M_P(i,j)) -= a_*h*lp;
			}
		}
	}
}

// ============================================================================
// [o3k] 4TH-ORDER FINITE-VOLUME IMPLICIT MAGNETIC PATH (third order in the low-Alfven regime).
// The compact (0.5 / 0.25 average) magnetic operator, its RHS and the face4/stagger4 momentum-energy
// update are the finite-DIFFERENCE / 2nd-order family on cell averages.  Where eps_B switches between
// neighbouring faces (M_A ~ 1 ring of the vortex) the implicit share (2nd order) and the explicit HLLEM
// share (WENO3 point states, transversally convolved) of the SAME flux differ by O(h^2), so the
// divergence at the switching cells carries an O(h) spike -- divided by rho in the momentum (rho0 = 1e-3:
// u error 0.34 / 0.050 / 0.025 at 32/64/128, dt-independent, grid-scale, localised in the ring).
// Everything here is routed through POINT VALUES (deconvolve -> pointwise product/flux -> convolve),
// face-averaged with facefv and differenced with the two-point divergence: every product is exact to
// O(h^4) without hand-written covariance terms, and the operator stays LINEAR in the unknown.
//   inner:  dm(x)  = -h div[(1-eps_B) F_lin(B_E, x)],   F_lin = 1/2 (B_E.x) I - B_E (x) x   (paper: |B|^2 ~ B^n.B^{n+1})
//           du     = dm / rho_I                          (point values)
//   outer:  A x    = x + h div[(1-eps_B) F_ind(du, B_E)] ,   F_ind(u,B) = u (x) B - B (x) u
//   RHS  :  b      = B*  - h div[(1-eps_B) F_ind(u*, B_E)],  u* = (m** - (1-eps_M) h grad p_E)/rho_I
// eps_B is the per-face splitting weight of THIS face at every nesting level (the old operator used
// (1-eps_B)^2 of the OUTER face for both levels).
// ============================================================================
static void mfv_zero(Grid& G){ const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++) for(int k = 0; k < NUM_VARS; k++) G(i,j,k) = 0.0; }
// out(i,j,k) = scale * div[(1-eps_B) facefv(F)]  on the interior, slots ks..ke; FX/FY hold CELL AVERAGES of the fluxes
static void mfv_div(Grid& FX, Grid& FY, Array& eBx, Array& eBy, int ks, int ke, Grid& out, double scale){
	// each weighted face average once (it was evaluated by both cells of the face); row buffers per thread
	const int g = ghost_cells, nk = ke - ks + 1;
#pragma omp parallel
	{
		std::vector<double> xf((size_t)(nxCells + 1) * nk), ym((size_t)nxCells * nk), yp((size_t)nxCells * nk);
		int jlast = -1000;
#pragma omp for schedule(static)
		for(int j = g; j < nyCells + g; j++){
			if(jlast == j - 1) ym.swap(yp);   // the upper faces of the previous row are this row's lower faces
			else for(int i = g; i < nxCells + g; i++) for(int k = ks; k <= ke; k++){
				auto fy = [&](int a,int b)->double{ return FY(a,b,k); };
				ym[(size_t)(i - g)*nk + (k - ks)] = (1.0 - eBy(i,j-1)) * facefvy(fy,i,j-1);
			}
			for(int i = g; i < nxCells + g; i++) for(int k = ks; k <= ke; k++){
				auto fy = [&](int a,int b)->double{ return FY(a,b,k); };
				yp[(size_t)(i - g)*nk + (k - ks)] = (1.0 - eBy(i,j)) * facefvy(fy,i,j);
			}
			for(int i = g - 1; i < nxCells + g; i++) for(int k = ks; k <= ke; k++){
				auto fx = [&](int a,int b)->double{ return FX(a,b,k); };
				xf[(size_t)(i - g + 1)*nk + (k - ks)] = (1.0 - eBx(i,j)) * facefvx(fx,i,j);
			}
			for(int i = g; i < nxCells + g; i++) for(int k = ks; k <= ke; k++){
				const double xpv = xf[(size_t)(i - g + 1)*nk + (k - ks)], xmv = xf[(size_t)(i - g)*nk + (k - ks)];
				const double ypv = yp[(size_t)(i - g)*nk + (k - ks)],     ymv = ym[(size_t)(i - g)*nk + (k - ks)];
				out(i,j,k) = scale * ((xpv - xmv)/dx + (ypv - ymv)/dy);
			}
			jlast = j;
		}
	}
}
// per-stage point-value caches
static Grid* g_mfv_Bn = nullptr;   // point values of Qe (B_E in slots 5..7)
static Grid* g_mfv_R  = nullptr;   // point values of u_con (rho_I in slot 0, m** in 1..3)
static long  g_mfv_stamp = -1;
static void mfv_prepare(Grid& u_con, Grid& Qe){
	if(!g_mfv_Bn){ g_mfv_Bn = new Grid(nxCells, nyCells); g_mfv_R = new Grid(nxCells, nyCells); }
	if(g_mfv_stamp == g_stage_stamp) return;
	g_mfv_stamp = g_stage_stamp;
	fv_cell_to(*g_mfv_Bn, Qe,    -1.0);
	fv_cell_to(*g_mfv_R,  u_con, -1.0, true);
}
// ---- N-component fields: the magnetic operator works on three components (B, or the momentum increment) and the
// energy / B_z sweep on four; an 8-slot Grid would move 8/N times the memory for the same arithmetic ----
static void fv_cell3(Field3 &F, double sgn){ fvN<3>(F, sgn); }
// bound_cond_u restricted to three slots starting at k0 (k0 = 1: momentum, k0 = 5: magnetic field)
static void bound_cond3(Field3 &u, int k0){
	switch(x_bc){
		case x_Transmissive:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++) for(int i = 0; i < ghost_cells; i++) for(int c = 0; c < 3; c++){
				u(i, j, c) = u(ghost_cells, j, c);
				u(nxCells + ghost_cells + i, j, c) = u(nxCells + ghost_cells - 1, j, c);
			}
			break;
		case x_Reflective:   // slots 1..3 and 5..7 are odd under the reflection
			for(int j = 0; j < nyCells + ghost_cells * 2; j++) for(int g = 0; g < ghost_cells; g++) for(int c = 0; c < 3; c++){
				u(g, j, c) = -u(2 * ghost_cells - 1 - g, j, c);
				u(nxCells + ghost_cells + g, j, c) = -u(nxCells + ghost_cells - 1 - g, j, c);
			}
			break;
		case x_Periodic:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++) for(int g = 0; g < ghost_cells; g++) for(int c = 0; c < 3; c++){
				u(g, j, c) = u(nxCells + g, j, c);
				u(nxCells + ghost_cells + g, j, c) = u(ghost_cells + g, j, c);
			}
			break;
		case x_Dirichlet:   // momentum 0, magnetic field transmissive
			for(int j = 0; j < nyCells + ghost_cells * 2; j++) for(int g = 0; g < ghost_cells; g++){
				const int L = g, R = nxCells + ghost_cells + g;
				for(int c = 0; c < 3; c++){
					if(k0 == 1){ u(L,j,c) = 0; u(R,j,c) = 0; }
					else{ u(L,j,c) = u(ghost_cells, j, c); u(R,j,c) = u(nxCells + ghost_cells - 1, j, c); }
				}
			}
			break;
	}
	switch(y_bc){
		case y_Transmissive:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++) for(int j = 0; j < ghost_cells; j++) for(int c = 0; c < 3; c++){
				u(i, j, c) = u(i, ghost_cells, c);
				u(i, nyCells + ghost_cells + j, c) = u(i, nyCells + ghost_cells - 1, c);
			}
			break;
		case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++) for(int g = 0; g < ghost_cells; g++) for(int c = 0; c < 3; c++){
				u(i, g, c) = -u(i, 2 * ghost_cells - 1 - g, c);
				u(i, nyCells + ghost_cells + g, c) = -u(i, nyCells + ghost_cells - 1 - g, c);
			}
			break;
		case y_Periodic:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++) for(int g = 0; g < ghost_cells; g++) for(int c = 0; c < 3; c++){
				u(i, g, c) = u(i, nyCells + g, c);
				u(i, nyCells + ghost_cells + g, c) = u(i, ghost_cells + g, c);
			}
			break;
		case y_Dirichlet:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++) for(int g = 0; g < ghost_cells; g++){
				const int B = g, T = nyCells + ghost_cells + g;
				for(int c = 0; c < 3; c++){
					if(k0 == 1){ u(i,B,c) = 0; u(i,T,c) = 0; }
					else{ u(i,B,c) = u(i, ghost_cells, c); u(i,T,c) = u(i, nyCells + ghost_cells - 1, c); }
				}
			}
			break;
	}
}
// mfv_div for an N-component field
template<int N> static void mfv_divN(const FieldN<N>& FX, const FieldN<N>& FY, Array& eBx, Array& eBy, FieldN<N>& out, double scale){
	const int g = ghost_cells, nk = N;
#pragma omp parallel
	{
		std::vector<double> xf((size_t)(nxCells + 1) * nk), ym((size_t)nxCells * nk), yp((size_t)nxCells * nk);
		int jlast = -1000;
#pragma omp for schedule(static)
		for(int j = g; j < nyCells + g; j++){
			if(jlast == j - 1) ym.swap(yp);
			else for(int i = g; i < nxCells + g; i++) for(int c = 0; c < nk; c++){
				auto fy = [&](int a,int b)->double{ return FY(a,b,c); };
				ym[(size_t)(i - g)*nk + c] = (1.0 - eBy(i,j-1)) * facefvy(fy,i,j-1);
			}
			for(int i = g; i < nxCells + g; i++) for(int c = 0; c < nk; c++){
				auto fy = [&](int a,int b)->double{ return FY(a,b,c); };
				yp[(size_t)(i - g)*nk + c] = (1.0 - eBy(i,j)) * facefvy(fy,i,j);
			}
			for(int i = g - 1; i < nxCells + g; i++) for(int c = 0; c < nk; c++){
				auto fx = [&](int a,int b)->double{ return FX(a,b,c); };
				xf[(size_t)(i - g + 1)*nk + c] = (1.0 - eBx(i,j)) * facefvx(fx,i,j);
			}
			for(int i = g; i < nxCells + g; i++) for(int c = 0; c < nk; c++){
				const double xpv = xf[(size_t)(i - g + 1)*nk + c], xmv = xf[(size_t)(i - g)*nk + c];
				const double ypv = yp[(size_t)(i - g)*nk + c],     ymv = ym[(size_t)(i - g)*nk + c];
				out(i,j,c) = scale * ((xpv - xmv)/dx + (ypv - ymv)/dy);
			}
			jlast = j;
		}
	}
}
static void mfv_div3(const Field3& FX, const Field3& FY, Array& eBx, Array& eBy, Field3& out, double scale){ mfv_divN<3>(FX, FY, eBx, eBy, out, scale); }
// per-stage compact copies of the point values prepared by mfv_prepare: B_E (slots 5..7 of Qe) and max(rho_pt, 1e-300)
static Field3* g_mfv_Bn3 = nullptr; static Array* g_mfv_Rr = nullptr; static long g_mfv3_stamp = -1;
static void mfv_compact(){
	if(!g_mfv_Bn3){ g_mfv_Bn3 = new Field3(); g_mfv_Rr = new Array(nxCells, nyCells); }
	if(g_mfv3_stamp == g_mfv_stamp) return;
	Grid& Bn = *g_mfv_Bn; Grid& R = *g_mfv_R; Field3& B3 = *g_mfv_Bn3; Array& Rr = *g_mfv_Rr;
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
#pragma omp parallel for schedule(static)
	for(int j = 0; j < NY; j++) for(int i = 0; i < NX; i++){ B3(i,j,0) = Bn(i,j,5); B3(i,j,1) = Bn(i,j,6); B3(i,j,2) = Bn(i,j,7); Rr(i,j) = fmax(R(i,j,0), 1e-300); }
	g_mfv3_stamp = g_mfv_stamp;
}
void Magnetic_Subsystem_AX_FV(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe, Grid& Qs,
                              Array& eps_M_grid_x, Array& eps_M_grid_y, Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	// three-component fields: X = B unknown, FX/FY = flux tensors, DM = momentum increment, OUT = induction divergence
	static Field3 X, FX, FY, DM, OUT;
	static Array shock_eta(nxCells,nyCells); static long eta_stamp = -1;
	if(eta_stamp != g_stage_stamp){ ShockDetector(shock_eta, Qe); eta_stamp = g_stage_stamp; }
	mfv_prepare(u_con, Qe); mfv_compact();
	const Field3& Bn3 = *g_mfv_Bn3; const Array& Rr = *g_mfv_Rr;   // B_E and max(rho_pt, 1e-300) at points
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	// 1. unknown as a field (B slots), boundary-conditioned, then point values
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){ X(i,j,0) = x(M_B(i,j,5)); X(i,j,1) = x(M_B(i,j,6)); X(i,j,2) = x(M_B(i,j,7)); }
	bound_cond3(X, 5);
	fv_cell3(X, -1.0);
	// 2. linearised Lorentz flux tensor at points:  F_ik = 1/2 (B_E . X) d_ik - B_E,i X_k   (components 1..3)
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double b5 = Bn3(i,j,0), b6 = Bn3(i,j,1), b7 = Bn3(i,j,2), x5 = X(i,j,0), x6 = X(i,j,1), x7 = X(i,j,2);
		const double bd = (b5*x5 + b6*x6 + b7*x7);   // [o3p] Newton: 1/2|B|^2 ~ B_E.B - 1/2|B_E|^2 (constant part on the RHS)
		FX(i,j,0) = bd - b5*x5;  FX(i,j,1) = -b5*x6;      FX(i,j,2) = -b5*x7;
		FY(i,j,0) = -b6*x5;      FY(i,j,1) = bd - b6*x6;  FY(i,j,2) = -b6*x7;
	}
	fv_cell3(FX, +1.0); fv_cell3(FY, +1.0);
	// 3. momentum increment (cell averages) with the inner faces' eps_B, then its point values / rho
	mfv_div3(FX, FY, eps_B_grid_x, eps_B_grid_y, DM, -h);
	bound_cond3(DM, 1);
	fv_cell3(DM, -1.0);
	// 4. induction fluxes of the velocity increment at points (components 5..7)
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double r = Rr(i,j);
		const double du = DM(i,j,0)/r, dv = DM(i,j,1)/r, dw = DM(i,j,2)/r;
		const double b5 = Bn3(i,j,0), b6 = Bn3(i,j,1), b7 = Bn3(i,j,2);
		FX(i,j,0) = 0.0;             FY(i,j,0) = dv*b5 - du*b6;
		FX(i,j,1) = du*b6 - dv*b5;   FY(i,j,1) = 0.0;
		FX(i,j,2) = du*b7 - dw*b5;   FY(i,j,2) = dv*b7 - dw*b6;
	}
	fv_cell3(FX, +1.0); fv_cell3(FY, +1.0);
	mfv_div3(FX, FY, eps_B_grid_x, eps_B_grid_y, OUT, +h);
	// 5. assemble, plus the chi-gated shock Laplacian (2nd-order compact, zero on smooth data)
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			const double eBx_p0 = 1.0 - eps_B_grid_x(i,j), eBx_m0 = 1.0 - eps_B_grid_x(i-1,j), eBy_0p = 1.0 - eps_B_grid_y(i,j), eBy_0m = 1.0 - eps_B_grid_y(i,j-1);
			const double gxp = 0.5*(shock_eta(i+1,j) + shock_eta(i,j)), gxm = 0.5*(shock_eta(i,j) + shock_eta(i-1,j));
			const double gyp = 0.5*(shock_eta(i,j+1) + shock_eta(i,j)), gym = 0.5*(shock_eta(i,j) + shock_eta(i,j-1));
			const double nu_x = B_max_x*dx, nu_y = B_max_y*dy;   // [o3s] kappa = chi |lambda_B| dx (paper eq. 84)
			for(int k = 5; k <= 7; k++){
				const double xc = x(M_B(i,j,k)), xe = x(M_B(i+1,j,k)), xw = x(M_B(i-1,j,k)), xn = x(M_B(i,j+1,k)), xs = x(M_B(i,j-1,k));
				const double lap = nu_x/(dx*dx)*(gxp*eBx_p0*(xe - xc) - gxm*eBx_m0*(xc - xw)) + nu_y/(dy*dy)*(gyp*eBy_0p*(xn - xc) - gym*eBy_0m*(xc - xs));
				const double diss = -h*lap;
				Ax(M_B(i,j,k)) = xc + OUT(i,j,k-5) + diss;
			}
		}
}
// RHS: b = B* - h div[(1-eps_B) F_ind(u*, B_E)],  u* = (m** - (1-eps_M) h grad p_E)/rho_I at points
static void Magnetic_RHS_FV(Grid& u_con, Grid& Qe, Grid& Qe_pri, Array& eps_M_grid_x, Array& eps_M_grid_y, Array& eps_B_grid_x, Array& eps_B_grid_y, double h, Eigen::VectorXd& b){
	static Grid U(nxCells,nyCells), FX(nxCells,nyCells), FY(nxCells,nyCells), OUT(nxCells,nyCells);
	static bool init = false; if(!init){ mfv_zero(U); mfv_zero(FX); mfv_zero(FY); mfv_zero(OUT); init = true; }
	mfv_prepare(u_con, Qe);
	Grid& Bn = *g_mfv_Bn; Grid& R = *g_mfv_R;
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	// [o3p] constant part of the Newton-linearised magnetic pressure, -1/2|B_E|^2 I: its momentum increment
	// -h div[(1-eps_B)(-1/2|B_E|^2 I)] belongs to the predictor velocity (the operator keeps the linear part).
	// With the paper's Picard linearisation |B|^2 ~ B^n.B^{n+1} the magnetic pressure in the predictor has half
	// the true stiffness while the momentum update applies the full force: the stage derivative of B~ is then
	// evaluated at a velocity O(dt) off the stage value -> second order in time and unstable for
	// dt k v_A(1-eps_B) >~ 3.5 (von Neumann of the full DIRK stage map); the Newton form is third order and
	// unconditionally stable (sym check), which matters at rho0 <= 1e-5 where v_A dt/dx ~ 2.
	static Grid DMC(nxCells,nyCells); static bool initc = false; if(!initc){ mfv_zero(DMC); initc = true; }
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double hb = -0.5*(Bn(i,j,5)*Bn(i,j,5) + Bn(i,j,6)*Bn(i,j,6) + Bn(i,j,7)*Bn(i,j,7));
		FX(i,j,1) = hb;  FX(i,j,2) = 0.0; FX(i,j,3) = 0.0;
		FY(i,j,1) = 0.0; FY(i,j,2) = hb;  FY(i,j,3) = 0.0;
	}
	fv_cell(FX, +1.0); fv_cell(FY, +1.0);
	mfv_div(FX, FY, eps_B_grid_x, eps_B_grid_y, 1, 3, DMC, -h);
	bound_cond_u(DMC);
	fv_cell(DMC, -1.0);
	// predictor velocity at points.  The (1-eps_M) h grad p_E term is an O(dt) predictor of the implicit
	// pressure correction; its spatial error enters B~ - B_CT multiplied by dt^2, so the compact
	// 2nd-order face form of the momentum update (cell quantity taken as a point value) is enough.
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double r = fmax(R(i,j,0), 1e-300);
		int ic = std::max(ghost_cells-1, std::min(nxCells+ghost_cells-1, i)), jc = std::max(ghost_cells-1, std::min(nyCells+ghost_cells-1, j));
		int icm = std::max(ghost_cells-1, std::min(nxCells+ghost_cells-1, i-1)), jcm = std::max(ghost_cells-1, std::min(nyCells+ghost_cells-1, j-1));
		const double gpx = ((1.0 - eps_M_grid_x(ic,jc))*0.5*(Qe_pri(i,j,4) + Qe_pri(i+1,j,4)) - (1.0 - eps_M_grid_x(icm,jc))*0.5*(Qe_pri(i-1,j,4) + Qe_pri(i,j,4)))/dx;
		const double gpy = ((1.0 - eps_M_grid_y(ic,jc))*0.5*(Qe_pri(i,j,4) + Qe_pri(i,j+1,4)) - (1.0 - eps_M_grid_y(ic,jcm))*0.5*(Qe_pri(i,j-1,4) + Qe_pri(i,j,4)))/dy;
		U(i,j,1) = (R(i,j,1) - h*gpx + DMC(i,j,1))/r; U(i,j,2) = (R(i,j,2) - h*gpy + DMC(i,j,2))/r; U(i,j,3) = (R(i,j,3) + DMC(i,j,3))/r;
	}
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double u = U(i,j,1), v = U(i,j,2), w = U(i,j,3), b5 = Bn(i,j,5), b6 = Bn(i,j,6), b7 = Bn(i,j,7);
		FX(i,j,5) = 0.0;           FY(i,j,5) = v*b5 - u*b6;
		FX(i,j,6) = u*b6 - v*b5;   FY(i,j,6) = 0.0;
		FX(i,j,7) = u*b7 - w*b5;   FY(i,j,7) = v*b7 - w*b6;
	}
	fv_cell(FX, +1.0); fv_cell(FY, +1.0);
	mfv_div(FX, FY, eps_B_grid_x, eps_B_grid_y, 5, 7, OUT, -h);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int k = 5; k <= 7; k++) b(M_B(i,j,k)) = u_con(i,j,k) + OUT(i,j,k);
}
// (1-eps_B)-weighted magnetic momentum / energy update with the solved B~ (u_test = {Qe hydro, B~})
static void MagneticFlux_Update_FV(Grid& u_con, Grid& u_test, Array& eps_B_grid_x, Array& eps_B_grid_y){
	static Grid Qp(nxCells,nyCells), FX(nxCells,nyCells), FY(nxCells,nyCells), OUT(nxCells,nyCells);
	static bool init = false; if(!init){ mfv_zero(Qp); mfv_zero(FX); mfv_zero(FY); mfv_zero(OUT); init = true; }
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	fv_cell_to(Qp, u_test, -1.0, true);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		double Q[NUM_VARS], F[NUM_VARS], G[NUM_VARS];
		for(int k = 0; k < NUM_VARS; k++) Q[k] = Qp(i,j,k);
		if(!(Q[0] > 0.0)) Q[0] = fmax(u_test(i,j,0), 1e-300);   // guard: point-value density of a limited state
		MagneticFlux_Pressure_x(Q, F); MagneticFlux_Pressure_y(Q, G);
		for(int k = 0; k < NUM_VARS; k++){ FX(i,j,k) = F[k]; FY(i,j,k) = G[k]; }
	}
	fv_cell(FX, +1.0); fv_cell(FY, +1.0);
	mfv_div(FX, FY, eps_B_grid_x, eps_B_grid_y, 1, 4, OUT, 1.0);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int k = 1; k <= 3; k++) u_con(i,j,k) -= adt * OUT(i,j,k);   // [o3r] momentum only; the energy row is applied in the Picard sweep
}
// [o3r] MAGNETIC BOOKKEEPING WITH ONE VELOCITY.  With e = E - KE - MB, the linearised pressure equation of a stage
// contains  -h B0^2 div(u_flux - u_ind): u_flux = velocity in the implicit magnetic energy flux, u_ind = velocity that
// evolved B.  Relative to the acoustic terms this is amplified by 1/beta.  The code used three velocities (u_E in the
// energy flux, the predictor u~ for Bz, the Picard velocity u_k for the CT'd in-plane B): at dt k v_A > ~1.3 and
// beta = 0.02 this was unstable (growth x4..x17 per step, cold magnetised plasma at rest).  Now, inside every Picard
// sweep and with the SAME current velocity u_k (and the explicit-stage field B_E, i.e. H(y_E, y_I) form):
//   in-plane B : CT (A advected with u_k)
//   Bz         : Bz = Bz* - h div[(1-eps_B) F_ind,z(u_k, B_E)]           (replaces the predictor B~_z)
//   energy     : E  = E*  - h div[(1-eps_B)(|B_E|^2 u_k - (u_k.B_E) B_E)]  (was applied once with u_E)
// The predictor B~ keeps only its role in the implicit Lorentz force.
static Array *g_Bx_st = nullptr, *g_By_st = nullptr;   // [1d] explicit-updated Bx, By of the current stage
static void Magnetic_Sweep_FV(Grid& u_con, Array& Bz_star, Array& E_star, Array& eps_B_grid_x, Array& eps_B_grid_y, double h){
	static FieldN<4> Qp, FX, FY, OUT;   // Qp: point values of slots 0..3 (rho, m); FX, FY, OUT: slots 4..7 (E, B)
	mfv_compact();
	const Field3& Bn = *g_mfv_Bn3;   // point values of B_E (prepared by the magnetic RHS of this stage)
	const int NX = nxCells + 2*ghost_cells, NY = nyCells + 2*ghost_cells;
	bound_cond_u(u_con);
	fvN_to<4>(Qp, GridSlots{u_con, 0}, -1.0, true);   // with the density point-value guard
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = 1; j < NY-1; j++) for(int i = 1; i < NX-1; i++){
		const double r = fmax(Qp(i,j,0), 1e-300);
		const double u = Qp(i,j,1)/r, v = Qp(i,j,2)/r, w = Qp(i,j,3)/r;
		const double b5 = Bn(i,j,0), b6 = Bn(i,j,1), b7 = Bn(i,j,2);
		const double b2 = b5*b5 + b6*b6 + b7*b7, ub = u*b5 + v*b6 + w*b7;
		for(int c = 0; c < 4; c++){ FX(i,j,c) = 0.0; FY(i,j,c) = 0.0; }
		FX(i,j,0) = b2*u - ub*b5;  FY(i,j,0) = b2*v - ub*b6;
		FX(i,j,3) = u*b7 - w*b5;   FY(i,j,3) = v*b7 - w*b6;
		if(!g_use_ct){ FY(i,j,1) = v*b5 - u*b6; FX(i,j,2) = u*b6 - v*b5; }   // [1d] in-plane induction with the same velocity and B_E
	}
	fvN<4>(FX, +1.0); fvN<4>(FY, +1.0);
	mfv_divN<4>(FX, FY, eps_B_grid_x, eps_B_grid_y, OUT, -h);
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		u_con(i,j,4) = E_star(i,j)  + OUT(i,j,0);
		u_con(i,j,7) = Bz_star(i,j) + OUT(i,j,3);
		if(!g_use_ct && g_Bx_st){ u_con(i,j,5) = (*g_Bx_st)(i,j) + OUT(i,j,1); u_con(i,j,6) = (*g_By_st)(i,j) + OUT(i,j,2); }
	}
	bound_cond_u(u_con);
}
void PreMagSubsystem_Pressure_HF(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, double h, Array& eps_M_grid_x, Array& eps_M_grid_y,Array& eps_B_grid_x, Array& eps_B_grid_y, const int rks){
	g_stage_stamp++;   // invalidate the per-stage AX caches
	static Array Bz_st(nxCells, nyCells), E_st(nxCells, nyCells);   // [o3r] explicit-updated Bz and E (before any implicit magnetic term)
	static Array Bx_st(nxCells, nyCells), By_st(nxCells, nyCells);   // [1d] explicit-updated in-plane field (used when there is no CT)
#pragma omp parallel for schedule(static)
	for(int j = 0; j < nyCells + 2*ghost_cells; j++) for(int i = 0; i < nxCells + 2*ghost_cells; i++){ Bz_st(i,j) = u_con(i,j,7); E_st(i,j) = u_con(i,j,4); Bx_st(i,j) = u_con(i,j,5); By_st(i,j) = u_con(i,j,6); }
	g_Bx_st = &Bx_st; g_By_st = &By_st;
	// ---- NEG-P DIAGNOSTIC (stage: explicit input) ----
	{
		static bool negp_reported_entry = false;
		if(!negp_reported_entry){   // first cell in (i, j) order with p <= 0, searched in parallel
			const int NYc = nyCells + 2*ghost_cells;
			long first = LONG_MAX;
#pragma omp parallel for collapse(2) schedule(static) reduction(min:first)
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				double rho=Qe(i,j,0), ke=0.5*(Qe(i,j,1)*Qe(i,j,1)+Qe(i,j,2)*Qe(i,j,2)+Qe(i,j,3)*Qe(i,j,3))/rho;
				double b2=Qe(i,j,5)*Qe(i,j,5)+Qe(i,j,6)*Qe(i,j,6)+Qe(i,j,7)*Qe(i,j,7);
				double pp=(gamma_plasma-1.0)*(Qe(i,j,4)-ke-0.5*b2);
				if(pp <= 0.0) first = std::min(first, (long)i*NYc + j);
			}
			if(first != LONG_MAX){
				const int i = (int)(first / NYc), j = (int)(first % NYc);
				double rho=Qe(i,j,0), ke=0.5*(Qe(i,j,1)*Qe(i,j,1)+Qe(i,j,2)*Qe(i,j,2)+Qe(i,j,3)*Qe(i,j,3))/rho;
				double b2=Qe(i,j,5)*Qe(i,j,5)+Qe(i,j,6)*Qe(i,j,6)+Qe(i,j,7)*Qe(i,j,7);
				double pp=(gamma_plasma-1.0)*(Qe(i,j,4)-ke-0.5*b2);
				std::cout << "NEGP stage=EXPLICIT-INPUT p=" << pp << " at ("<<i<<","<<j<<")" << std::endl; negp_reported_entry = true;
			}
		}
	}

	static Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	g_Hbar_stamp = -1;   // [o3f] force a rebuild of the cell-average enthalpy for this stage (u_con changed)
	static Eigen::VectorXd b_p_vect(nxCells*nyCells), b_Bt_vect(3 * nxCells * nyCells), B_vect(3 * nxCells * nyCells), p_vect(nxCells * nyCells);   // every entry is set before use

	typedef Eigen::Triplet<double> T;
	
	// [o3k] magnetic predictor B~: 4th-order FV right-hand side and operator, initial guess B_E
	Magnetic_RHS_FV(u_con, Qe, Qe_pri, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h, b_Bt_vect);
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){ B_vect(M_B(i,j,5)) = Qe(i,j,5); B_vect(M_B(i,j,6)) = Qe(i,j,6); B_vect(M_B(i,j,7)) = Qe(i,j,7); }
	
	GMRES_Solver(B_vect, b_Bt_vect, Magnetic_Subsystem_AX_FV, u_con, Qe, Qs,eps_M_grid_x, eps_M_grid_y,eps_B_grid_x, eps_B_grid_y, h, 3*nxCells*nyCells);

#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			u_con(i,j,5) = B_vect(M_B(i, j, 5));
			u_con(i,j,6) = B_vect(M_B(i, j, 6));
			u_con(i,j,7) = B_vect(M_B(i, j, 7));
		}
	}
	bound_cond_u(u_con);
	
	static Grid u_test(nxCells, nyCells);
#pragma omp parallel for schedule(static)
	for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
		for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
			u_test(i,j,0) = Qe(i,j,0);
			u_test(i,j,1) = Qe(i,j,1);
			u_test(i,j,2) = Qe(i,j,2);
			u_test(i,j,3) = Qe(i,j,3);
			u_test(i,j,4) = Qe(i,j,4);
			
			u_test(i,j,5) = u_con(i,j,5);
			u_test(i,j,6) = u_con(i,j,6);
			u_test(i,j,7) = u_con(i,j,7);
		}
	}
	
	MagneticFlux_Update_FV(u_con, u_test, eps_B_grid_x, eps_B_grid_y);   // [o3k] (1-eps_B) magnetic momentum/energy update with B~
		
	bound_cond_u(u_con);
	
	static Array rhovx_star(nxCells, nyCells);
	static Array rhovy_star(nxCells, nyCells);
	static Array rhovz_star(nxCells, nyCells);
	
#pragma omp parallel for schedule(static)
	for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			rhovx_star(i,j) = u_con(i,j,1);
			rhovy_star(i,j) = u_con(i,j,2);
			rhovz_star(i,j) = u_con(i,j,3);
			
		}
	}
	bound_cond_array(rhovx_star);
	bound_cond_array(rhovy_star);
	bound_cond_array(rhovz_star);
	build_Hbar(Qe, u_con);   // [o3f] cell-average enthalpy for the RHS assembly below (stamp was reset above)
	static Array g_ebar(nxCells, nyCells), g_kembar(nxCells, nyCells);
	// [o3m] Picard sweep = CT (A advected with the current stage velocity, B := curl4 A) -> pressure solve
	// -> momentum update.  Every pressure solve therefore balances the energy with a DIVERGENCE-FREE field
	// (the predictor B~ only enters the Lorentz force and the magnetic energy flux, which need the implicit
	// coupling); the CT is idempotent (A rebuilt from A_s), so the last CT after the last momentum update
	// leaves B consistent with the stored momentum.
	auto ct_now = [&](){   // [o3r] CT + Bz + magnetic energy flux, all with the current velocity
		if(g_use_ct && magfield == UnStaggerCT && g_ct_A) Unstaggered_CT(u_con, *g_ct_A, *g_ct_As, *g_ct_Ae, h);
		Magnetic_Sweep_FV(u_con, Bz_st, E_st, eps_B_grid_x, eps_B_grid_y, h);
	};
	// [emb] Phi(m) = E(m) - kl(Q_E; m, B(m)):  E(m) = E* - h div[(1-eps_B) F^B_U(u, B_E)],  B_perp(m) = curl4 A(m),
	// A(m) = A_s - h u.grad A_E (upwind direction of grad A_E from u_E),  B_z(m) = B_z* - h div[(1-eps_B) F_ind,z(u, B_E)].
	// Phi is AFFINE in m, and m = m** - h G p, so the stage energy equation
	//     L4 p - Phi(m** - h G p) = -h div[(1-eps_M) H m**_f]
	// is LINEAR in p: one GMRES solve of  A p := L4 p + Phi(m**) - Phi(m** - h G p) = b + Phi(m**).
	static Grid emb_W(nxCells, nyCells), emb_A(nxCells, nyCells);
	static Array emb_phi0(nxCells, nyCells), emb_phi(nxCells, nyCells), emb_kl(nxCells, nyCells);
	auto emb_eval = [&](Grid& S, Grid& Aout, Array& phi){   // S: state with the trial momentum; fills B, B_z, E and Phi
		if(g_use_ct) Unstaggered_CT(S, Aout, *g_ct_As, *g_ct_Ae, h);
		Magnetic_Sweep_FV(S, Bz_st, E_st, eps_B_grid_x, eps_B_grid_y, h);
		kemb_lin(S, Qe, emb_kl);
		const int NXa = nxCells + 2*ghost_cells, NYa = nyCells + 2*ghost_cells;
#pragma omp parallel for schedule(static)
		for(int j = 0; j < NYa; j++) for(int i = 0; i < NXa; i++) phi(i,j) = S(i,j,4) - emb_kl(i,j);
	};
	const bool emb_on = magfield == UnStaggerCT && g_ct_A;   // in 1D (g_use_ct = false) the same embedded solve runs without the CT map
	if(emb_on){
		copy(emb_W, u_con); copy(emb_A, *g_ct_A);
		emb_eval(emb_W, emb_A, emb_phi0);
#pragma omp parallel for schedule(static)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++) for(int i = 0; i < nxCells + 2*ghost_cells; i++) g_ebar(i,j) = emb_phi0(i,j);
	}else{
	ct_now();
	cellavg_e_kemb(u_con, g_ebar, g_kembar);   // [o3g] e_bar of the star state (E**, m**, B_CT(u**))
	}

	if(g_verbose_stages) std::cout << "rk = " <<rks + 1 <<", First Internal: "<<endl;
	vector<T> triplets_p;
	#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
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

			double h00 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j,4) / Qe(i,j,0) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i + 1,j,4) / Qe(i + 1,j,0) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i - 1,j,4) / Qe(i - 1,j,0) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j + 1,4) / Qe(i,j + 1,0) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j - 1,4) / Qe(i,j - 1,0) * Qe(i,j - 1,0) / u_con(i,j - 1,0);

			double p_star = u_con(i,j,4) - 0.5 * B2_00;
			double internal = p_star - 0.5*(rhovx_star(i,j)*rhovx_star(i,j) + rhovy_star(i,j)*rhovy_star(i,j) + rhovz_star(i,j)*rhovz_star(i,j))/u_con(i,j,0);
			double RHS;

			// ---- 4TH-ORDER enthalpy-flux divergence (upgrade of the 2nd-order
			// face-average form). SAME (1-epsM) single-power face weighting and
			// SAME enthalpy definition (gamma/(gamma-1) * p_Qe / rho_ucon) as the
			// verified 2nd-order fused step; only the spatial discretization is
			// raised: 4th-order face interpolation (face4x/y) for enthalpy and
			// momentum, 4-point staggered divergence (stagger4). This exact form
			// was validated to order 4.000 against analytic derivatives. The LHS
			// elliptic operator is ALSO 4th order: the momentum update carries h^1
			// and needs O(dx^3)+, and the elliptic operator is its composition with
			// the energy divergence, so it must live in the same face4/stagger4
			// family (pairing law).
			auto henth4 = [&](int ii,int jj)->double{
				// POSITIVITY GUARD on the operator coefficient: the RK stage combination Qe
				// carries no positivity guarantee (non-SSP DIRK weights), and a negative
				// p_Qe here makes the elliptic coefficient NEGATIVE -> anti-diffusion,
				// ill-posed operator, Newton/GMRES stall (observed: hit-max-iters every
				// solve on the 256^2 blastwave). Clamping to p_floor keeps the operator
				// well-posed; applied identically in LHS and RHS (operator-pairing law).
				return g_Hbar(ii,jj);   // [o3f] cell-average enthalpy (p_Qe / rho_ucon, both floored, via point values)
			};
			auto rvx4 = [&](int ii,int jj)->double{ return rhovx_star(ii,jj); };
			auto rvy4 = [&](int ii,int jj)->double{ return rhovy_star(ii,jj); };
			// [o3d] enthalpy flux = product of two face averages + transverse covariance (dt^2/12) h_t (rho u)_t
			auto HFx = [&](int ic)->double{   // [o3e] FV family
				double prod = facefvx(henth4,ic,j)*facefvx(rvx4,ic,j);
				prod += (facefvx(henth4,ic,j+1) - facefvx(henth4,ic,j-1)) * ((facefvx(rvx4,ic,j+1)) - (facefvx(rvx4,ic,j-1))) / 48.0;
				return (1.0-eps_M_grid_x(ic,j)) * prod;
			};
			auto HFy = [&](int jc)->double{
				double prod = facefvy(henth4,i,jc)*facefvy(rvy4,i,jc);
				prod += (facefvy(henth4,i+1,jc) - facefvy(henth4,i-1,jc)) * ((facefvy(rvy4,i+1,jc)) - (facefvy(rvy4,i-1,jc))) / 48.0;
				return (1.0-eps_M_grid_y(i,jc)) * prod;
			};
			double xdiv4 = (HFx(i) - HFx(i-1)) / dx;
			double ydiv4 = (HFy(j) - HFy(j-1)) / dy;
			RHS = g_ebar(i,j) - h * (xdiv4 + ydiv4);   // [o3g] cell-average internal energy

			if(internal < 1e-6){
				if(g_verbose_stages) std::cout << internal << " ";
			}
			
			b_p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
			
			p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = Qe_pri(i,j,4);
		}
	}
	
	std::cout << endl;

	if(emb_on){
		auto AX_emb = [&](Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& uc, Grid& qe, Grid& qs, Array& a1, Array& a2, Array& a3, Array& a4, double hh){
			const bool sv = floor_on; floor_on = false;
			Pressure_Subsystem_AX_HF(Ax, x, uc, qe, qs, a1, a2, a3, a4, hh);
			floor_on = sv;
			copy(emb_W, u_con);
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = ghost_cells; j < ghost_cells + nyCells; j++) for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				auto pm = [&](int ii,int jj)->double{ return x(M_P(ii,jj)); };
				auto FMx = [&](int ic)->double{ return (1.0 - eps_M_grid_x(ic,j)) * facefvx(pm,ic,j); };
				auto FMy = [&](int jc)->double{ return (1.0 - eps_M_grid_y(i,jc)) * facefvy(pm,i,jc); };
				emb_W(i,j,1) = rhovx_star(i,j) - h * (FMx(i) - FMx(i-1)) / dx;
				emb_W(i,j,2) = rhovy_star(i,j) - h * (FMy(j) - FMy(j-1)) / dy;
			}
			bound_cond_u(emb_W);
			emb_eval(emb_W, emb_A, emb_phi);
#pragma omp parallel for schedule(static)
			for(int j = ghost_cells; j < ghost_cells + nyCells; j++) for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				const bool bnd = ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1));
				if(!bnd) Ax(M_P(i,j)) += emb_phi0(i,j) - emb_phi(i,j);
			}
			g_emb_iters++;
		};
		g_emb_solves++;
		{
			// [emb2] UNIFIED POSITIVITY.  A_emb p = b is solved by defect correction on the M-matrix operator
			//     L2c p_{k+1} = b + (L2c - A_emb) p_k ,   L2c = L2 + upwind convective term,
			// every inner solve by the nested Newton with the floor: every iterate is admissible (p >= p_floor),
			// the converged one satisfies the embedded (4th-order, linearised) system wherever the floor is inactive.
			static Eigen::VectorXd dAx2(nxCells*nyCells), dAxE(nxCells*nyCells), dbeff(nxCells*nyCells), p_old(nxCells*nyCells);
#pragma omp parallel for schedule(static)
			for(int q = 0; q < nxCells*nyCells; q++) p_vect(q) = fmax(p_vect(q), p_floor);
			double dp_prev = 1e300, dp_now = 0.0; int idc = 0;
			g_emb_L2c = true;
			for(idc = 0; idc < 40; idc++){
				const bool sv = floor_on; floor_on = false;
				omp_zero(dAx2); omp_zero(dAxE);
				Pressure_Subsystem_AX(dAx2, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
				floor_on = sv;
				AX_emb(dAxE, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
				{ const int n = nxCells*nyCells;
#pragma omp parallel for schedule(static)
				  for(int q = 0; q < n; q++) dbeff[q] = b_p_vect[q] + dAx2[q] - dAxE[q]; }
				omp_copy(p_old, p_vect);
				Pressure_NestedNewton(p_vect, dbeff, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, h, nxCells*nyCells);
				g_emb_sweeps++;
				dp_now = omp_dist(p_vect, p_old) / (omp_nrm(p_vect) + 1e-300);
				if(dp_now < 1e-12) break;
				if(idc >= 3 && dp_now > dp_prev) break;
				dp_prev = dp_now;
			}
			g_emb_L2c = false;
			if(dp_now > 1e-10){ g_emb_nc++; std::cout << "  [EMB-DC] WARNING: not converged, |dp|/|p| = " << dp_now << " after " << idc + 1 << " sweeps" << std::endl; }
			{ const double pm = omp_min(p_vect); if(pm <= 0.0) std::cout << "NEGP stage=EMBEDDED pmin=" << pm << std::endl; }
		}
	}else
	{
		// ================= DEFECT-CORRECTION PRESSURE SOLVE =================
		// Goal: 4th-order accuracy AND the M-matrix positivity theorem, decoupled.
		// Solve L4 p = b by iterating   L2 p_{k+1} = b + (L2 - L4) p_k :
		// every inner solve runs Pressure_NestedNewton on the 2ND-ORDER M-MATRIX
		// operator (nested-Newton positivity machinery fully active, theorem
		// applies), and at convergence p satisfies the 4TH-ORDER discretization.
		// NOTE floor_on is forced OFF while forming the operator difference:
		// the archive AX's diagonal depends on floor_on (V_slope vs linear), and
		// letting the EOS modification leak into (L2 - L4) would corrupt the
		// correction; the EOS handling lives INSIDE Pressure_NestedNewton only.
		const int NDC = 3;
		Eigen::VectorXd dc_Ax2(nxCells*nyCells), dc_Ax4(nxCells*nyCells), dc_beff(nxCells*nyCells);
		// [o3q] residual-controlled: at least NDC sweeps, continue (up to NDCmax) while the update is not yet
		// at roundoff, and report non-convergence instead of silently accepting it.
		const int NDCmax = 12; double dp_prev = 1e300, dp_now = 0.0; int idc = 0;
		for(idc = 0; idc < NDCmax; idc++){
			bool dc_save = floor_on; floor_on = false;
			dc_Ax2.setZero(); dc_Ax4.setZero();
			Pressure_Subsystem_AX(dc_Ax2, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
			Pressure_Subsystem_AX_HF(dc_Ax4, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
			floor_on = dc_save;
			dc_beff = b_p_vect + dc_Ax2 - dc_Ax4;
			Eigen::VectorXd p_old = p_vect;
			Pressure_NestedNewton(p_vect, dc_beff, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, h, nxCells*nyCells);
			// ---- NEG-P DIAGNOSTIC (after inner positivity solve) ----
			if(p_vect.minCoeff() <= 0.0) std::cout << "NEGP stage=AFTER-INNER-SOLVE idc=" << idc << " pmin=" << p_vect.minCoeff() << std::endl;
			dp_now = (p_vect - p_old).norm() / (p_vect.norm() + 1e-300);
			if(idc + 1 >= NDC && dp_now < 1e-12) break;
			if(idc >= NDC && dp_now > dp_prev) break;   // not contracting: stop, reported below
			dp_prev = dp_now;
		}
		if(dp_now > 1e-10) std::cout << "  [DC] WARNING: defect correction not converged, |dp|/|p| = " << dp_now << " after " << idc + 1 << " sweeps" << std::endl;
	}

	if(g_verbose_stages) std::cout << "rk = " <<rks + 1<< ", First Solution: "<<std::endl;
	if(g_verbose_stages) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			if(g_verbose_stages && p_vect(M_P(i,j))  < 1e-6){
				if(g_verbose_stages) std::cout << p_vect(M_P(i,j))  << " ";
			}
		}
	}
	std::cout << endl;
	// [o3i] PICARD ITERATION OF THE (p, rho u) PAIR.  The stage equations are
	//     m = m** - a dt (1-eps) grad p ,   p/(g-1) + KE(m) = E** - a dt div[(1-eps) H m_f(p)] ,
	// coupled through KE(m).  The paper's steps 5-8 (solve with KE(m**) -> m1 = m** - a dt grad p~ ->
	// solve with KE(m1) -> E = p2/(g-1) + KE(m1)) keep the momentum of the FIRST solve: m1 - m_I =
	// -a dt grad(p~ - p_I) = O(dt^2) per stage (curl-free), i.e. FIRST order in time for the compressive
	// momentum, hence for rho (mass flux) and p (= c^2 drho at low Mach: measured on the rho0 = 1e-2 vortex,
	// p error halving with dt at 128^2 while u was h^3).  Symbolic check (sym/picard.py): re-updating the
	// momentum after every solve gains one power of dt per solve; 3 solves give m: O(dt^4), E: O(dt^5) per
	// stage, i.e. the DIRK's own order.  E is assembled from the LAST pressure and the KE of the PREVIOUS
	// momentum, which is exactly the conservative form E** - a dt div[(1-eps) H m_f(p_K)].
	auto momentum_update = [&](){
#pragma omp parallel for collapse(2) schedule(static)
	for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			// ============================================================
			// MOMENTUM UPDATE: 4TH-ORDER, SAME OPERATOR FAMILY AS THE LHS
			// (user decision: M-matrix positivity dropped for 4th order).
			//   rhov^{n+1} = rhov* - h * stagger4( (1-epsM)_face * face4(p) )
			// This is the 4th-order analogue of the verified 2nd-order
			// face-difference form -- NOT eq.(81)'s cell-point gradient, which
			// belongs to no operator family and produced the rho~1.23 plateau
			// (bisection-isolated). Pairing with the 4th-order composed elliptic
			// LHS in Pressure_Subsystem_AX_HF: same face4, same stagger4, epsM
			// weights taken directly at faces.
			// ============================================================
			auto pmom = [&](int ii,int jj)->double{ return p_vect(M_P(ii,jj)); };
			auto FMx4 = [&](int ic)->double{ return (1.0 - eps_M_grid_x(ic,j)) * facefvx(pmom,ic,j); };   // [o3e] FV family
			auto FMy4 = [&](int jc)->double{ return (1.0 - eps_M_grid_y(i,jc)) * facefvy(pmom,i,jc); };
			double gpx4 = (FMx4(i) - FMx4(i-1)) / dx;
			double gpy4 = (FMy4(j) - FMy4(j-1)) / dy;
			u_con(i,j,1) = rhovx_star(i,j) - h * gpx4;
			u_con(i,j,2) = rhovy_star(i,j) - h * gpy4;
		}
	}
	bound_cond_u(u_con);
	};   // momentum_update
	momentum_update();
	if(emb_on){   // [emb] B, B_z, E(mag) of the solved momentum; E_bar = p/(g-1) + kl is the conservative energy
		emb_eval(u_con, *g_ct_A, emb_phi);
#pragma omp parallel for schedule(static)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++) for(int i = 0; i < nxCells + 2*ghost_cells; i++) g_kembar(i,j) = emb_kl(i,j);
	}
//				//std::cout << Qe(i,j,k) << " ";
//				cout << setw(12)              // 每个数字宽度12字符
//				<< fixed                 	  // 固定小数点形式（非科学计数法）
//				<< setprecision(6)            // 保留6位小数

	bound_cond_u(u_con);
	
//	
	if(g_verbose_stages) std::cout << "rk = " <<rks + 1<< ", Second Internal: "<<endl;
	for(int ipic = 1; ipic < (emb_on ? 1 : g_npicard); ipic++)
	{
		ct_now();                                  // [o3m] B := curl A(m_{ipic})
		cellavg_e_kemb(u_con, g_ebar, g_kembar);   // [o3g] e_bar with the UPDATED momentum m_{ipic} (E still E**)
		vector<T> triplets_p_second;
		#pragma omp parallel for collapse(2) schedule(static)
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
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
	
//								
	
				double h00 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j,4) / Qe(i,j,0) * Qe(i,j,0) / u_con(i,j,0);
				double hp0 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i + 1,j,4) / Qe(i + 1,j,0) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
				double hm0 = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i - 1,j,4) / Qe(i - 1,j,0) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
				double h0p = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j + 1,4) / Qe(i,j + 1,0) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
				double h0m = gamma_plasma / (gamma_plasma - 1.0) * Qe_pri(i,j - 1,4) / Qe(i,j - 1,0) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
				
				double p_star = u_con(i,j,4) - 0.5 * B2_00;
				double internal = p_star - 0.5*(rhovx_star(i,j)*rhovx_star(i,j) + rhovy_star(i,j)*rhovy_star(i,j) + rhovz_star(i,j)*rhovz_star(i,j))/u_con(i,j,0);

				double RHS;
				
			// ---- 4TH-ORDER enthalpy-flux divergence (upgrade of the 2nd-order
				// face-average form). SAME (1-epsM) single-power face weighting and
				// SAME enthalpy definition (gamma/(gamma-1) * p_Qe / rho_ucon) as the
				// verified 2nd-order fused step; only the spatial discretization is
				// raised: 4th-order face interpolation (face4x/y) for enthalpy and
				// momentum, 4-point staggered divergence (stagger4). This exact form
				// was validated to order 4.000 against analytic derivatives. The LHS
				// elliptic operator deliberately STAYS 2nd-order compact (M-matrix /
				// positivity-proof requirement); it carries h^2 so h^2*O(dx^2)=O(dx^4)
				// -- the h^1 RHS below is what actually needed 4th order.
				auto henth4 = [&](int ii,int jj)->double{
					// POSITIVITY GUARD on the operator coefficient: the RK stage combination Qe
				// carries no positivity guarantee (non-SSP DIRK weights), and a negative
				// p_Qe here makes the elliptic coefficient NEGATIVE -> anti-diffusion,
				// ill-posed operator, Newton/GMRES stall (observed: hit-max-iters every
				// solve on the 256^2 blastwave). Clamping to p_floor keeps the operator
				// well-posed; applied identically in LHS and RHS (operator-pairing law).
				return g_Hbar(ii,jj);   // [o3f] cell-average enthalpy (p_Qe / rho_ucon, both floored, via point values)
				};
				auto rvx4 = [&](int ii,int jj)->double{ return rhovx_star(ii,jj); };
				auto rvy4 = [&](int ii,int jj)->double{ return rhovy_star(ii,jj); };
				// [o3d] enthalpy flux = product of two face averages + transverse covariance (dt^2/12) h_t (rho u)_t
				auto HFx = [&](int ic)->double{   // [o3e] FV family
					double prod = facefvx(henth4,ic,j)*facefvx(rvx4,ic,j);
					prod += (facefvx(henth4,ic,j+1) - facefvx(henth4,ic,j-1)) * ((facefvx(rvx4,ic,j+1)) - (facefvx(rvx4,ic,j-1))) / 48.0;
					return (1.0-eps_M_grid_x(ic,j)) * prod;
				};
				auto HFy = [&](int jc)->double{
					double prod = facefvy(henth4,i,jc)*facefvy(rvy4,i,jc);
					prod += (facefvy(henth4,i+1,jc) - facefvy(henth4,i-1,jc)) * ((facefvy(rvy4,i+1,jc)) - (facefvy(rvy4,i-1,jc))) / 48.0;
					return (1.0-eps_M_grid_y(i,jc)) * prod;
				};
				double xdiv4 = (HFx(i) - HFx(i-1)) / dx;
				double ydiv4 = (HFy(j) - HFy(j-1)) / dy;
				RHS = g_ebar(i,j) - h * (xdiv4 + ydiv4);   // [o3g] cell-average internal energy (with the updated momentum)

				if(internal < 1e-8){
					if(g_verbose_stages) std::cout << internal << " ";
				}
				b_p_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = RHS;
				
			}
		}
		std::cout<<endl;
		
		{
			// ================= DEFECT-CORRECTION PRESSURE SOLVE =================
			// Goal: 4th-order accuracy AND the M-matrix positivity theorem, decoupled.
			// Solve L4 p = b by iterating   L2 p_{k+1} = b + (L2 - L4) p_k :
			// every inner solve runs Pressure_NestedNewton on the 2ND-ORDER M-MATRIX
			// operator (nested-Newton positivity machinery fully active, theorem
			// applies), and at convergence p satisfies the 4TH-ORDER discretization.
			// NOTE floor_on is forced OFF while forming the operator difference:
			// the archive AX's diagonal depends on floor_on (V_slope vs linear), and
			// letting the EOS modification leak into (L2 - L4) would corrupt the
			// correction; the EOS handling lives INSIDE Pressure_NestedNewton only.
			const int NDC = 3;
			Eigen::VectorXd dc_Ax2(nxCells*nyCells), dc_Ax4(nxCells*nyCells), dc_beff(nxCells*nyCells);
			// [o3q] residual-controlled: at least NDC sweeps, continue (up to NDCmax) while the update is not yet
			// at roundoff, and report non-convergence instead of silently accepting it.
			const int NDCmax = 12; double dp_prev = 1e300, dp_now = 0.0; int idc = 0;
			for(idc = 0; idc < NDCmax; idc++){
				bool dc_save = floor_on; floor_on = false;
				dc_Ax2.setZero(); dc_Ax4.setZero();
				Pressure_Subsystem_AX(dc_Ax2, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
				Pressure_Subsystem_AX_HF(dc_Ax4, p_vect, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, h);
				floor_on = dc_save;
				dc_beff = b_p_vect + dc_Ax2 - dc_Ax4;
				Eigen::VectorXd p_old = p_vect;
				Pressure_NestedNewton(p_vect, dc_beff, u_con, Qe, Qs, eps_M_grid_x, eps_M_grid_y, h, nxCells*nyCells);
				// ---- NEG-P DIAGNOSTIC (after inner positivity solve) ----
				if(p_vect.minCoeff() <= 0.0) std::cout << "NEGP stage=AFTER-INNER-SOLVE idc=" << idc << " pmin=" << p_vect.minCoeff() << std::endl;
				dp_now = (p_vect - p_old).norm() / (p_vect.norm() + 1e-300);
				if(idc + 1 >= NDC && dp_now < 1e-12) break;
				if(idc >= NDC && dp_now > dp_prev) break;   // not contracting: stop, reported below
				dp_prev = dp_now;
			}
			if(dp_now > 1e-10) std::cout << "  [DC] WARNING: defect correction not converged, |dp|/|p| = " << dp_now << " after " << idc + 1 << " sweeps" << std::endl;
		}
		
		if(g_verbose_stages) std::cout << "rk = " <<rks + 1<< ", Second Solution: " << std::endl;
		if(g_verbose_stages) for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				if(g_verbose_stages && p_vect(M_P(i,j))  < 1e-6){
					if(g_verbose_stages) std::cout << p_vect(M_P(i,j))  << " ";
				}
			}
		}
		std::cout << endl;
		
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
//				
		momentum_update();   // [o3i] m_{ipic+1} = m** - a dt (1-eps) grad p_{ipic+1}
	}
	ct_now();   // [o3m] final CT: B consistent with the stored momentum m_K
	
	bound_cond_u(u_con);

	if(g_verbose_stages) std::cout << "rk = " <<rks + 1<< ", Internal Energy: "<< endl;
	{	// energy assembly: E_bar = p_K/(g-1) + avg(KE(m_{K-1}) + MB), the conservative form of the last solve
#pragma omp parallel for schedule(static)
			for(int j = ghost_cells; j < ghost_cells + nyCells; j++)
				for(int i = ghost_cells; i < ghost_cells + nxCells; i++)
					u_con(i,j,4) = 1.0/(gamma_plasma - 1.0) * p_vect(M_P(i,j)) + g_kembar(i,j);   // E = p/(g-1) + linearised KE + MB of the solved equation
			if(g_verbose_stages) for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					double internal_energy = u_con(i,j,4) - (0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7))) - (0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) );
					
					if(g_verbose_stages && internal_energy < 0){
						std::cout << "E_total=" << u_con(i,j,4) << " ";
						std::cout << "KE=" << (0.5 * (u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0)) << " ";
						std::cout << "ME=" << (0.5 * (u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))) << " ";
						std::cout << "IE=" << internal_energy;
					}
				}
			}
			std::cout<<endl;
	}
	
	bound_cond_u(u_con);
	
	if(nyCells != 1 && nxCells != 1){
		// TODO-2D: adapt to our 6-arg Staggered_CT(u_con,Bx_p0,By_0p,Bx_p0_Qs,By_0p_Qs,h) before running 2D
	}
	
	bound_cond_u(u_con);
}

// ============================================================================
// CollectMagFluxes_Imp -- AX wave flux + the EMF-representable part of the
// viscosity, evaluated on the CONVERGED B (u_con holds Btilde^{n+1}; the
// momenta are still the star values, since MagneticFlux_Update_Pressure_HF
// runs later). Signs read off the AX assembly:
//     Ax(...,6) = x(...,6) + (h*h/dx)*( eBx_p0^2*(By_1+By_2) - ... )
//     Ax(...,6) = Ax(...,6) - h*lap_By
// so with By^{n+1} = ... - (h/dx)*(F_p0 - F_m0):
//     F_p0 = +h*eBx_p0^2*(By_1+By_2) - nu_x/dx*eBx_p0*(By_{i+1}-By_i)
// The WAVE operator carries (1-eps_B)^2 -- the momentum update and the
// induction equation contribute one factor each -- while the VISCOSITY, being
// a single application, carries (1-eps_B)^1.
// By's y-diffusion and Bx's x-diffusion are not components of an EMF and are
// deliberately omitted: they are the self-check residual floor.
// ============================================================================

// 4th-order corner EMF: transverse (-1,9,9,-1)/16 interpolation of the face
// EMFs, both routes averaged. This is where the ACCURACY of the CT lives --
// the Stokes update below is an exact identity, so raising its difference
// operator buys nothing and destroys the divergence cancellation.
// No extra dissipation is added: it is already inside the collected fluxes.

// ============================================================================
// CT_EpsSplit -- two Stokes updates, one per difference operator.
//   group A: B* share      -> stagger4
//   group B: explicit+imp  -> 2-point
// Each is the curl of its own corner EMF with matched operators, so each
// preserves the staggered divergence separately; the sum does too.
// ============================================================================

int main(int argc, char *argv[]){
	for(int a = 1; a < argc; a++){ std::string k(argv[a]); auto eq = k.find('='); if(eq == std::string::npos) continue; std::string key = k.substr(0,eq); double v = atof(k.c_str()+eq+1);
		if(key=="cfl") CFL = v; else if(key=="epsMmin") eps_M_min = v; else if(key=="epsMmax") eps_M_max = v; else if(key=="epsBmin") eps_B_min = v; else if(key=="epsBmax") eps_B_max = v;
		else if(key=="plap") g_plap = (int)v; else if(key=="tag") g_tag = k.substr(eq+1); else if(key=="kappa1") shock_kappa1 = v;
		else if(key=="swm0") g_sw_m0 = v; else if(key=="swm1") g_sw_m1 = v; else if(key=="swalpha") g_sw_alpha = v; else if(key=="epsmax3") g_epsmax3 = (int)v; else if(key=="pswitch") g_pswitch = (int)v; else if(key=="hllemflat") g_hllemflat = (int)v; else if(key=="hllemoff") g_hllemoff = (int)v;
		else { std::cout << "unknown arg " << k << std::endl; return 1; } }
	std::cout << "ARGS cfl=" << CFL << " epsM=[" << eps_M_min << "," << eps_M_max << "] epsB=[" << eps_B_min << "," << eps_B_max << "] plap=" << g_plap << " kappa1=" << shock_kappa1 << " switch(m0,m1,alpha)=" << g_sw_m0 << "," << g_sw_m1 << "," << g_sw_alpha << " epsmax3=" << g_epsmax3 << " pswitch=" << g_pswitch << " hllemflat=" << g_hllemflat << " hllemoff=" << g_hllemoff << std::endl;
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
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test4:
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
			time_stop = 5.0/(2.0*M_PI);
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			gamma_plasma = 5.0/3.0;
			break;
		case testblastwave:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.01;
			gamma_plasma = 1.4;   // Boscheri-Thomann 2024
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
			gamma_plasma = 5.0/3.0;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testVortex:
			x_0 = -VORTEX_L, x_1 = VORTEX_L, y_0 = -VORTEX_L, y_1 = VORTEX_L, z_0 = -1, z_1 = 1;
			time_stop = 1;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testSmooth2D:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1, z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
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
			time_stop = 0.8 / KH_MX;   // t_f = t_max/6, t_max = 4.8/M_x
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
	Array eps_M(nxCells,nyCells), eps_B(nxCells,nyCells);
	g_epsM_cell = &eps_M;
	Array eps_M_x(nxCells,nyCells), eps_M_y(nxCells,nyCells);
	Array eps_B_x(nxCells,nyCells), eps_B_y(nxCells,nyCells);
	u_initialization(u_pri, Bx_p0, By_0p);
	upri2ucon(u_con, u_pri);
	// Point values -> cell averages.  u_initialization samples the initial
	// condition at the cell centre; a finite volume scheme evolves cell
	// averages, and the two differ by (h^2/24) * Laplacian.  Leaving that in
	// injects an O(h^2) error at t = 0 that no spatial order can recover, so a
	// convergence study without this conversion can never measure more than
	// second order.  Only for smooth initial data.
	{   // fix (2): jumps of the quasi-periodic potential from the initial mean field
		double bx0 = 0.0, by0 = 0.0; long n = 0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){ bx0 += u_con(i,j,5); by0 += u_con(i,j,6); n++; }
		bx0 /= n; by0 /= n;
		g_Az_yjump = (y_bc == y_Periodic) ? bx0 * (y_1 - y_0) : 0.0;
		g_Az_xjump = (x_bc == x_Periodic) ? by0 * (x_1 - x_0) : 0.0;
		if(g_Az_yjump != 0.0 || g_Az_xjump != 0.0)
			std::cout << "quasi-periodic A: <Bx>=" << bx0 << " <By>=" << by0
			          << " -> Az_yjump=" << g_Az_yjump << " Az_xjump=" << g_Az_xjump << std::endl;
	}
	A_initial(A);
	// Re-initialise the cell-centred magnetic field from curl A, so that the
	// discrete divergence-free condition holds from t = 0.  Without this the
	// initial B is the analytic one and does not satisfy the 4th-order central
	// divergence that curl A satisfies identically; the mismatch would appear
	// as a fixed O(dx^4) divergence floor for the whole run.
	if(g_use_ct && (magfield == UnStaggerCT || magfield == PotentialA)){   // 1D: keep the initial cell averages of B (no potential)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double mb_old = 0.5*(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6));
				u_con(i,j,5) =  ( -A(i,j+2,2) + 8.0*A(i,j+1,2) - 8.0*A(i,j-1,2) + A(i,j-2,2) ) / (12.0*dy);
				u_con(i,j,6) = -( -A(i+2,j,2) + 8.0*A(i+1,j,2) - 8.0*A(i-1,j,2) + A(i-2,j,2) ) / (12.0*dx);
				double mb_new = 0.5*(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6));
				u_con(i,j,4) += mb_new - mb_old;   // [o3] keep the internal energy: E was built with the quadrature B
			}
		bound_cond_u(u_con);
		ucon2upri(u_pri, u_con);
	}

	std::vector<double> bRK(RKs);
	Array Ae(RKs, RKs), Ai(RKs, RKs);
	Tensor Hs(RKs, nxCells, nyCells);
	Array shock_eta(nxCells, nyCells);
	// per-step log in the current directory
	auto maxB = [&](){ double m = 0.0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max:m)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++) for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		m = fmax(m, sqrt(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))); return m; };
	const double maxB0 = maxB();
	g_outdir = "plot/" + test_name() + "/";
	std::filesystem::create_directories(g_outdir);
	FILE* steplog = fopen((g_outdir + "steplog_" + test_name() + ".csv").c_str(), "w");
	fprintf(steplog, "step,t,dt,dt_MHD,maxB,maxB_over_maxB0,divB_L2,divB_Linf\n");
	{ double L2, Li; divB_norms(u_con, L2, Li); fprintf(steplog, "0,%.10e,0,0,%.10e,%.10e,%.6e,%.6e\n", t, maxB0, 1.0, L2, Li); }
	std::vector<double> IntermediatePlotTimes = output_times();
	bool isFinalTime = false;
	
	getepsgrid(u_con, eps_M, eps_B,eps_M_x, eps_M_y, eps_B_x, eps_B_y);
	IntermediatePlotIntermediatePlot(u_con, u_pri, eps_M,eps_B,timestep, IntermediatePlotTimes, isFinalTime);

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
	}else if(RKs == 4){
		// SA DIRK(3,4,3) -- Boscheri & Pareschi, JCP 434 (2021) eq.(79).
		// order-3 conditions verified; b == last row of Ai (stiffly accurate);
		// singly diagonal: a_ii = d for every stage.
		// [o3q] FULL-PRECISION coefficients.  The 6-digit values had sum(b) = 1 - 1e-6 and every order condition
		// off by 1e-7..1e-6: a dt- and h-independent error floor (~1e-7 in L2 on the vortex at T = 0.2).
		// d = root in (0,1) of x^3 - 3x^2 + 3/2 x - 1/6 (L-stable SDIRK3), c3 = (1+d)/2,
		// b2 = -(6d^2-16d+1)/4, b3 = (6d^2-20d+5)/4; explicit rows solve row sums c_E = c_I (stages 2..4) and
		// b.Ae.c = 1/6 with a42 = 1/2.  All implicit, explicit and coupling conditions of order 3 hold to 1e-17.
		const double d   = 0.43586652150845899942;
		const double a32 = 0.28206673924577050029;   // c3 - d
		const double b2  = 1.2084966491760100703;
		const double b3  = -0.64436317068446906975;
		Ae(0,0)=0.0;                  Ae(0,1)=0.0;                   Ae(0,2)=0.0;                   Ae(0,3)=0.0;
		Ae(1,0)=d;                    Ae(1,1)=0.0;                   Ae(1,2)=0.0;                   Ae(1,3)=0.0;
		Ae(2,0)=1.4377449737640151;   Ae(2,1)=-0.71981171300978564;  Ae(2,2)=0.0;                   Ae(2,3)=0.0;
		Ae(3,0)=0.91699329835202015;  Ae(3,1)=0.5;                   Ae(3,2)=-0.41699329835202015;  Ae(3,3)=0.0;

		Ai(0,0)=d;    Ai(0,1)=0.0;  Ai(0,2)=0.0;  Ai(0,3)=0.0;
		Ai(1,0)=0.0;  Ai(1,1)=d;    Ai(1,2)=0.0;  Ai(1,3)=0.0;
		Ai(2,0)=0.0;  Ai(2,1)=a32;  Ai(2,2)=d;    Ai(2,3)=0.0;
		Ai(3,0)=0.0;  Ai(3,1)=b2;   Ai(3,2)=b3;   Ai(3,3)=d;

		bRK[0]=0.0;  bRK[1]=b2;  bRK[2]=b3;  bRK[3]=d;
	}
	do{
		timestep++;
		g_bc_prim = true;  bound_cond_u(u_pri);  g_bc_prim = false;  // primitive ghosts (slot 4 = p)
		calch(u_pri);
		calcp(u_pri);
		calca(u_pri);
		
		Array dt_eps_M_x(nxCells,nyCells), dt_eps_M_y(nxCells,nyCells);
		Array dt_eps_B_x(nxCells,nyCells), dt_eps_B_y(nxCells,nyCells);
		getepsgrid(u_con, eps_M, eps_B,dt_eps_M_x, dt_eps_M_y, dt_eps_B_x, dt_eps_B_y);
		dt = caldt(u_pri, dt_eps_M_x, dt_eps_M_y, dt_eps_B_x, dt_eps_B_y);
		// DT GROWTH LIMITER: the material-speed CFL of the full_step logic is BLIND
		// to the "huge pressure jump, zero velocity" start-up of the blastwave --
		// at ts=2 the material a_max ~ 1 while the explicit sub-system still
		// carries eps~0.1 shares of the acoustic/Alfven fans (speeds ~5-7), and a
		// single 64x-larger step wrecked the magnetic linearization (B^n.B^{n+1}
		// assumes small dB): magE doubled in one step, IE went negative, Newton
		// saturated. Cap the growth at 1.5x per step: the transient is bridged
		// smoothly and the physical spectrum takes over within a few steps.
		{
			static double dt_prev = -1.0;
			if(dt_prev > 0.0 && dt > 1.5 * dt_prev) dt = 1.5 * dt_prev;
			dt_prev = dt;
		}
		dt = fmin(dt, time_stop - t);
		for(double tout : IntermediatePlotTimes) if(tout > t + 1e-12*time_stop && tout - t < dt) dt = tout - t;   // land on the output times
		if(fabs(dt - (time_stop - t)) < 1e-8){
			isFinalTime = true;
		}
		const double dt_mhd = CFL / (ch_x/dx + ch_y/dy);   // fully explicit step (fast magnetosonic speed), same CFL
		
		static Grid A_n_snapshot(nxCells, nyCells);
		copy(A_n_snapshot, A);

		MagneticEigenvalues(u_pri);
		EnergyEigenvalues(u_pri);
		t = t + dt;
		
		upri2ucon(u_con, u_pri);
		
		copy(u_con_n,u_con);
		
#pragma omp parallel for collapse(2) schedule(static)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++){
			for(int i = 0; i < nxCells + 2*ghost_cells; i++){
				Bx_p0_n(i,j) = Bx_p0(i,j); 
				By_0p_n(i,j) = By_0p(i,j);
			}
		}
		
		for(int rks = 0; rks < RKs; rks++){
			static Grid A_s(nxCells, nyCells), A_e(nxCells, nyCells);
			// stage combinations  Q_E = Q^n + dt sum_r a^E_r H_r,  Q_S = Q^n + dt sum_r a^I_r H_r  (and the same for A, B_x, B_y),
			// accumulated in the order r = 0, 1, .. per element; A starts from A_S
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
				for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
					for(int k = 0; k < NUM_VARS; k++){
						double qe = u_con_n(i,j,k), qs = qe;
						for(int r = 0; r < rks; r++){
							qe = qe + Ae(rks,r) * dt * Hs(r,i,j,k);
							qs = qs + Ai(rks,r) * dt * Hs(r,i,j,k);
						}
						Qe(i,j,k) = qe; Qs(i,j,k) = qs;
					}
					for(int k = 0; k < NUM_VARS; k++){
						double ae = A_n_snapshot(i,j,k), as = ae;
						if(k < 3) for(int r = 0; r < rks; r++){   // A 只有 3 分量
							ae += Ae(rks,r) * dt * HsA(r,i,j,k);
							as += Ai(rks,r) * dt * HsA(r,i,j,k);
						}
						A_e(i,j,k) = ae; A_s(i,j,k) = as; A(i,j,k) = as;
					}
					double bxe = Bx_p0_n(i,j), bye = By_0p_n(i,j), bxs = Bx_p0_n(i,j), bys = By_0p_n(i,j);
					for(int r = 0; r < rks; r++){
						bxe += Ae(rks,r) * dt * HsB(r,i,j,0);
						bye += Ae(rks,r) * dt * HsB(r,i,j,1);
						bxs += Ai(rks,r) * dt * HsB(r,i,j,0);
						bys += Ai(rks,r) * dt * HsB(r,i,j,1);
					}
					Bx_p0_Qe(i,j) = bxe; By_0p_Qe(i,j) = bye; Bx_p0_Qs(i,j) = bxs; By_0p_Qs(i,j) = bys;
				}
			}
			
			adt =  Ai(rks,rks) * dt;
			
			// epsM/epsB from the CURRENT explicit-stage state Qe (not the
			// previous full timestep's u_con) -- consistent with how
			// Magnetic_Subsystem_AX/MagSubsystem_UpdateB/TVDReconstruction
			// below all reference Qe as "the" explicit state for this stage.
			Array eps_M_grid_x(nxCells,nyCells), eps_M_grid_y(nxCells,nyCells);
			g_epsM_fx = &eps_M_grid_x; g_epsM_fy = &eps_M_grid_y;
			Array eps_B_grid_x(nxCells,nyCells), eps_B_grid_y(nxCells,nyCells);
			getepsgrid(Qe, eps_M, eps_B,eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y);

			WENO3Reconstruction(Qe, w_xL, w_xR, w_yL, w_yR, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y);

			// Convective Subsystem F^c + split explicit F^p,F^B (HLLEM: AP-preserving)
			g_hllem_Qe = &Qe;
			HLLEM_Update(u_con, Qs, w_xL,w_xR,w_yL,w_yR, eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, adt);
			bound_cond_u(u_con);

			// ============================================================
			// PressureFormula now goes through the PORTED fused implicit step
			// (PreMagSubsystem_Pressure_HF, verified-correct half-finished code):
			// solves B implicitly, applies the (1-epsB) magnetic momentum/energy
			// bookkeeping, solves p, updates momentum and energy -- all with
			// exactly-paired operators. EnergyFormula keeps the old pipeline.
			// ============================================================
			{
				g_ct_A = &A; g_ct_As = &A_s; g_ct_Ae = &A_e; g_emb_Qe = &Qe;   // [o3m] [emb]
				PreMagSubsystem_Pressure_HF(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, Bx_p0, By_0p, adt,
				                            eps_M_grid_x, eps_M_grid_y, eps_B_grid_x, eps_B_grid_y, rks);
				bound_cond_u(u_con);
				// keep Bx/By/Bz mirrors in sync for any downstream use
#pragma omp parallel for collapse(2) schedule(static)
				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
					for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
						Bx(i,j) = u_con(i,j,5);
						By(i,j) = u_con(i,j,6);
						Bz(i,j) = u_con(i,j,7);
					}
				}
			}
			// (legacy EnergyFormula pipeline removed; fused PressureFormula step above is the only maintained path)

			// (the constrained transport is part of the pressure solve of the stage: PreMagSubsystem_Pressure_HF leaves
			//  B = curl4 A(m^{n+1}) and E = p/(g-1) + kl; nothing is done to B or E here)
			bound_cond_u(u_con);
			
			// Modify Energy
			
			bound_cond_u(u_con);
			ucon2upri(u_pri, u_con);
			
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
				for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
					for(int k = 0; k < NUM_VARS; k++){
						Hs(rks,i,j,k) = (u_con(i,j,k) - Qs(i,j,k)) / adt;
					}
				}
			}
			g_epsM_fx = nullptr; g_epsM_fy = nullptr;   // the stage's face-eps arrays go out of scope below
			
#pragma omp parallel for collapse(2) schedule(static)
			for(int j = 0; j < nyCells + 2*ghost_cells; j++){
				for(int i = 0; i < nxCells + 2*ghost_cells; i++){
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
		
		{
			double L2, Li; divB_norms(u_con, L2, Li); const double mB = maxB();
			fprintf(steplog, "%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.6e,%.6e\n", timestep, t, dt, dt_mhd, mB, mB/(maxB0 > 0.0 ? maxB0 : 1.0), L2, Li);
			fflush(steplog);
			std::cout << "step " << timestep << "  t = " << t << "  dt = " << dt << "  dt/dt_MHD = " << dt/dt_mhd << "  |divB|_L2 = " << L2 << std::endl;
		}
		IntermediatePlotIntermediatePlot(u_con, u_pri, eps_M,eps_B, timestep, IntermediatePlotTimes, isFinalTime);
	}while(t < time_stop);
	fclose(steplog);
	std::cout << "pressure solves " << g_emb_solves << ", defect-correction sweeps " << g_emb_sweeps << ", not converged " << g_emb_nc
	          << ", density point-value guard applied " << g_rho_guard << " times, entropy anti-diffusion limited at " << g_ent_lim << " faces, envelope widened at " << g_env_wide << " faces, pressure LU builds " << g_lu_builds << " (" << g_lu_time << " s)" << std::endl;
	if(test == testVortex){   // paper Table 2: L2 errors of rho, u, p, Bx, Az against the exact (translated) cell averages
		const double gp[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)}; const double gw[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
		double l2[5] = {0,0,0,0,0}, li[5] = {0,0,0,0,0};
		// [o3i] cell averages of the NONLINEAR maps u = m/rho and p = (g-1)(E - KE - MB) via point values:
		// u(Qbar), p(Qbar) differ from ubar, pbar at O(h^2) -- the measurement itself would cap them at 2nd order.
		static Array u_pt(nxCells, nyCells), p_pt(nxCells, nyCells);
		bound_cond_u(u_con);
		for(int i = 1; i < nxCells + 2*ghost_cells - 1; i++) for(int j = 1; j < nyCells + 2*ghost_cells - 1; j++){
			double Qp[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++) Qp[k] = u_con(i,j,k) - (u_con(i+1,j,k) - 2.0*u_con(i,j,k) + u_con(i-1,j,k) + u_con(i,j+1,k) - 2.0*u_con(i,j,k) + u_con(i,j-1,k))/24.0;
			double ke = 0.5*(Qp[1]*Qp[1]+Qp[2]*Qp[2]+Qp[3]*Qp[3])/Qp[0], mb = 0.5*(Qp[5]*Qp[5]+Qp[6]*Qp[6]+Qp[7]*Qp[7]);
			u_pt(i,j) = Qp[1]/Qp[0]; p_pt(i,j) = (gamma_plasma-1.0)*(Qp[4] - ke - mb);
		}
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++) for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double xc = x_0 + (i - ghost_cells + 0.5) * dx, yc = y_0 + (j - ghost_cells + 0.5) * dy;
			double ex[5] = {0,0,0,0,0};
			for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++){
				double xg = xc + 0.5*dx*gp[a] - VORTEX_U0*t, yg = yc + 0.5*dy*gp[b] - VORTEX_U0*t;   // translated back, wrapped into the periodic box
				xg = xg - (x_1 - x_0)*floor((xg - x_0)/(x_1 - x_0)); yg = yg - (y_1 - y_0)*floor((yg - y_0)/(y_1 - y_0));
				const double r2 = xg*xg + yg*yg, e1 = exp(0.5*(1.0 - r2)), vt = VORTEX_VT, w = gw[a]*gw[b];
				ex[0] += w*g_vortex_rho0; ex[1] += w*(VORTEX_U0 + vt/(2.0*M_PI)*e1*(-yg));
				ex[2] += w*(1.0 + 0.5*exp(1.0 - r2)*(g_vortex_Bt*g_vortex_Bt/(4.0*M_PI*M_PI)*(1.0 - r2) - g_vortex_rho0*(vt/(2.0*M_PI))*(vt/(2.0*M_PI))));
				ex[3] += w*(g_vortex_Bt/(2.0*M_PI)*e1*(-yg)); ex[4] += w*(g_vortex_Bt/(2.0*M_PI)*e1);
			}
			double rho = u_con(i,j,0), uu = u_con(i,j,1)/rho, ke = 0.5*(u_con(i,j,1)*u_con(i,j,1)+u_con(i,j,2)*u_con(i,j,2)+u_con(i,j,3)*u_con(i,j,3))/rho;
			double mb = 0.5*(u_con(i,j,5)*u_con(i,j,5)+u_con(i,j,6)*u_con(i,j,6)+u_con(i,j,7)*u_con(i,j,7));
			double ubar = u_pt(i,j) + (u_pt(i+1,j) - 2.0*u_pt(i,j) + u_pt(i-1,j) + u_pt(i,j+1) - 2.0*u_pt(i,j) + u_pt(i,j-1))/24.0;   // [o3i]
			double pbar = p_pt(i,j) + (p_pt(i+1,j) - 2.0*p_pt(i,j) + p_pt(i-1,j) + p_pt(i,j+1) - 2.0*p_pt(i,j) + p_pt(i,j-1))/24.0;
			(void)uu; (void)ke; (void)mb;
			double num[5] = {rho, ubar, pbar, u_con(i,j,5), A(i,j,2)};
			for(int q = 0; q < 5; q++){ double e = num[q] - ex[q]; l2[q] += e*e*dx*dy; li[q] = fmax(li[q], fabs(e)); }
		}
		{	// POINT-VALUE errors: p(Qbar) and u(Qbar) differ from the cell averages of p, u at O(h^2)
			// (nonlinear maps), so the cell-average comparison above is capped at 2nd order for
			// them by the MEASUREMENT.  Deconvolve Qbar -> Q(x_c) (4th order, unlimited: smooth data)
			// and compare with the exact point values at the cell centre.
			double l2p[5] = {0,0,0,0,0};
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++) for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double Qp[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++)
					Qp[k] = u_con(i,j,k) - (u_con(i+1,j,k) - 2.0*u_con(i,j,k) + u_con(i-1,j,k) + u_con(i,j+1,k) - 2.0*u_con(i,j,k) + u_con(i,j-1,k))/24.0;
				double Ap = A(i,j,2) - (A(i+1,j,2) - 2.0*A(i,j,2) + A(i-1,j,2) + A(i,j+1,2) - 2.0*A(i,j,2) + A(i,j-1,2))/24.0;
				double xc = x_0 + (i - ghost_cells + 0.5) * dx - VORTEX_U0*t, yc = y_0 + (j - ghost_cells + 0.5) * dy - VORTEX_U0*t;
				xc = xc - (x_1 - x_0)*floor((xc - x_0)/(x_1 - x_0)); yc = yc - (y_1 - y_0)*floor((yc - y_0)/(y_1 - y_0));
				const double r2 = xc*xc + yc*yc, e1 = exp(0.5*(1.0 - r2)), vt = VORTEX_VT;
				double ex[5] = { g_vortex_rho0, VORTEX_U0 + vt/(2.0*M_PI)*e1*(-yc),
					1.0 + 0.5*exp(1.0 - r2)*(g_vortex_Bt*g_vortex_Bt/(4.0*M_PI*M_PI)*(1.0 - r2) - g_vortex_rho0*(vt/(2.0*M_PI))*(vt/(2.0*M_PI))),
					g_vortex_Bt/(2.0*M_PI)*e1*(-yc), g_vortex_Bt/(2.0*M_PI)*e1 };
				double rho = Qp[0], ke = 0.5*(Qp[1]*Qp[1]+Qp[2]*Qp[2]+Qp[3]*Qp[3])/rho, mb = 0.5*(Qp[5]*Qp[5]+Qp[6]*Qp[6]+Qp[7]*Qp[7]);
				double num[5] = {rho, Qp[1]/rho, (gamma_plasma-1.0)*(Qp[4]-ke-mb), Qp[5], Ap};
				for(int q = 0; q < 5; q++){ double e = num[q] - ex[q]; l2p[q] += e*e*dx*dy; }
			}
			{	// [diag] CONSERVED-variable cell-average errors: E and m_x are linear in the state, so the
			// cell-average comparison is exact (no nonlinear-map O(h^2) term): these measure the energy
			// and momentum equations directly.
			double l2E = 0.0, l2m = 0.0;
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++) for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double xc = x_0 + (i - ghost_cells + 0.5) * dx, yc = y_0 + (j - ghost_cells + 0.5) * dy, exE = 0.0, exm = 0.0;
				for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++){
					double xg = xc + 0.5*dx*gp[a] - VORTEX_U0*t, yg = yc + 0.5*dy*gp[b] - VORTEX_U0*t;
					xg = xg - (x_1 - x_0)*floor((xg - x_0)/(x_1 - x_0)); yg = yg - (y_1 - y_0)*floor((yg - y_0)/(y_1 - y_0));
					const double r2 = xg*xg + yg*yg, e1 = exp(0.5*(1.0 - r2)), vt = VORTEX_VT, w = gw[a]*gw[b];
					double rho = g_vortex_rho0, uu = VORTEX_U0 + vt/(2.0*M_PI)*e1*(-yg), vv = VORTEX_U0 + vt/(2.0*M_PI)*e1*xg;
					double pp = 1.0 + 0.5*exp(1.0 - r2)*(g_vortex_Bt*g_vortex_Bt/(4.0*M_PI*M_PI)*(1.0 - r2) - g_vortex_rho0*(vt/(2.0*M_PI))*(vt/(2.0*M_PI)));
					double bx = g_vortex_Bt/(2.0*M_PI)*e1*(-yg), by = g_vortex_Bt/(2.0*M_PI)*e1*xg;
					exE += w*(pp/(gamma_plasma - 1.0) + 0.5*rho*(uu*uu + vv*vv) + 0.5*(bx*bx + by*by)); exm += w*rho*uu;
				}
				double eE = u_con(i,j,4) - exE, em = u_con(i,j,1) - exm; l2E += eE*eE*dx*dy; l2m += em*em*dx*dy;
			}
			std::cout << "VORTEX-CONS N=" << nxCells << " t=" << t << " L2(cell avg): E=" << sqrt(l2E) << " mx=" << sqrt(l2m) << std::endl;
		}
		std::cout << "VORTEX-PT N=" << nxCells << " t=" << t << " L2(point): rho=" << sqrt(l2p[0]) << " u=" << sqrt(l2p[1]) << " p=" << sqrt(l2p[2]) << " Bx=" << sqrt(l2p[3]) << " Az=" << sqrt(l2p[4]) << std::endl;
		}
		std::cout << "VORTEX N=" << nxCells << " rho0=" << g_vortex_rho0 << " t=" << t
		          << " L2: rho=" << sqrt(l2[0]) << " u=" << sqrt(l2[1]) << " p=" << sqrt(l2[2]) << " Bx=" << sqrt(l2[3]) << " Az=" << sqrt(l2[4])
		          << " | Linf: rho=" << li[0] << " u=" << li[1] << " p=" << li[2] << " Bx=" << li[3] << " Az=" << li[4] << std::endl;
	}
	if(test == testSmooth2D){   // [o3] L1/Linf error of rho against the exact cell average (3x3 Gauss)
		const double gp[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)}; const double gw[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
		double l1 = 0.0, linf = 0.0, l2 = 0.0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++) for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double xc = x_0 + (i - ghost_cells + 0.5) * dx, yc = y_0 + (j - ghost_cells + 0.5) * dy, ex = 0.0;
			for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) ex += gw[a]*gw[b]*(1.0 + 0.2*sin(2.0*M_PI*(xc + 0.5*dx*gp[a] + yc + 0.5*dy*gp[b] - 2.0*t)));
			double e = fabs(u_con(i,j,0) - ex); l1 += e*dx*dy; l2 += e*e*dx*dy; linf = fmax(linf, e);
		}
		std::cout << "SMOOTH2D N=" << nxCells << " t=" << t << " L1(rho)=" << l1 << " L2(rho)=" << sqrt(l2) << " Linf(rho)=" << linf << std::endl;
	}
	
	ofstream outFilerho(g_outdir + "RK_rho.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerho << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,0) << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx(g_outdir + "RK_vx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,1) << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy(g_outdir + "RK_vy.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,2) << endl;
		}
	}
	outFilevy.close();
	
	ofstream outFileE(g_outdir + "RK_E.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileE << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,4) << endl;
		}
	}
	outFileE.close();

	ofstream outFilep(g_outdir + "RK_p.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilep << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,4)) << endl;
		}
	}
	outFilep.close();
	
	ofstream outFileBx(g_outdir + "RK_Bx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,5)) << endl;
		}
	}
	outFileBx.close();

	ofstream outFileBy(g_outdir + "RK_By.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,6)) << endl;
		}
	}
	outFileBy.close();
	
	ofstream outFileBz(g_outdir + "RK_Bz.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBz << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,7)) << endl;
		}
	}
	outFileBz.close();

	ofstream outFilev(g_outdir + "RK_v_mag.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilev << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2)) << endl;
		}
	}
	outFilev.close();
	
	ofstream outFilep_line(g_outdir + "RK_p_line.dat");
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
	
	ofstream outFilerho_line(g_outdir + "RK_rho_line.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			if(j == j_mid){
				outFilerho_line << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,0) << endl;
			}
		}
	}
	outFilerho_line.close();
	
	ofstream outFileepsM(g_outdir + "RK_eps_M.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileepsM << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << eps_M(i,j) << endl;
		}
	}
	outFileepsM.close();
	
	ofstream outFileepsB(g_outdir + "RK_eps_B.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileepsB << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << eps_B(i,j) << endl;
		}
	}
	outFileepsB.close();
}