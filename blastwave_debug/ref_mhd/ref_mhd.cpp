// Independent reference solver for the MHD blast wave (same initial data as the two IMEX codes):
// 2D ideal MHD, gamma = 1.4, rho = 1, u = 0, Bx = 70/sqrt(4 pi), p = 1000 (r <= 0.1) / 0.1,
// domain [-0.5,0.5]^2, transmissive boundaries, t_end = 0.01.
// Scheme: fully explicit, cell-centred, MUSCL (minmod on primitives) + HLL flux, SSP-RK2,
// Powell (8-wave) source terms for div B, CFL 0.4 on the fast magnetosonic speed.
// Purpose: a third, independent discretisation to judge which features of the two IMEX
// results are physical (converged) and which are numerical.
// Build: g++ -O3 -march=native -fopenmp -std=c++17 ref_mhd.cpp -o ref_mhd
// Run:   ./ref_mhd N [smooth=1] [cfl=0.4]  -> writes ref_N.dat (x y rho u v w p Bx By Bz)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
static const int NV = 8;
static const double GAM = 1.4;
struct Grid { int nx, ny, ng; std::vector<double> d; Grid(int nx_, int ny_, int ng_): nx(nx_), ny(ny_), ng(ng_), d((size_t)(nx_+2*ng_)*(ny_+2*ng_)*NV, 0.0) {}
	inline double* at(int i, int j){ return &d[((size_t)j*(nx+2*ng) + i)*NV]; } };
static inline void con2pri(const double* c, double* p){
	p[0] = c[0]; p[1] = c[1]/c[0]; p[2] = c[2]/c[0]; p[3] = c[3]/c[0];
	double ke = 0.5*(c[1]*c[1]+c[2]*c[2]+c[3]*c[3])/c[0], mb = 0.5*(c[5]*c[5]+c[6]*c[6]+c[7]*c[7]);
	p[4] = (GAM-1.0)*(c[4]-ke-mb); p[5] = c[5]; p[6] = c[6]; p[7] = c[7];
}
static inline void pri2con(const double* p, double* c){
	c[0] = p[0]; c[1] = p[0]*p[1]; c[2] = p[0]*p[2]; c[3] = p[0]*p[3];
	c[4] = p[4]/(GAM-1.0) + 0.5*p[0]*(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]) + 0.5*(p[5]*p[5]+p[6]*p[6]+p[7]*p[7]);
	c[5] = p[5]; c[6] = p[6]; c[7] = p[7];
}
// flux in direction d (0 = x, 1 = y) from primitive state
static inline void flux(const double* p, int d, double* f){
	double rho = p[0], u = p[1], v = p[2], w = p[3], pr = p[4], Bx = p[5], By = p[6], Bz = p[7];
	double B2 = Bx*Bx+By*By+Bz*Bz, pt = pr + 0.5*B2, E = pr/(GAM-1.0) + 0.5*rho*(u*u+v*v+w*w) + 0.5*B2, uB = u*Bx+v*By+w*Bz;
	if(d == 0){
		f[0] = rho*u; f[1] = rho*u*u + pt - Bx*Bx; f[2] = rho*u*v - Bx*By; f[3] = rho*u*w - Bx*Bz;
		f[4] = (E+pt)*u - Bx*uB; f[5] = 0.0; f[6] = u*By - v*Bx; f[7] = u*Bz - w*Bx;
	}else{
		f[0] = rho*v; f[1] = rho*v*u - By*Bx; f[2] = rho*v*v + pt - By*By; f[3] = rho*v*w - By*Bz;
		f[4] = (E+pt)*v - By*uB; f[5] = v*Bx - u*By; f[6] = 0.0; f[7] = v*Bz - w*By;
	}
}
static inline double cfast(const double* p, int d){
	double a2 = GAM*p[4]/p[0], B2 = (p[5]*p[5]+p[6]*p[6]+p[7]*p[7])/p[0], bn2 = (d==0 ? p[5]*p[5] : p[6]*p[6])/p[0];
	double s = a2 + B2; return sqrt(0.5*(s + sqrt(std::max(s*s - 4.0*a2*bn2, 0.0))));
}
static inline void hll(const double* pL, const double* pR, int d, double* F){
	double cL[NV], cR[NV], fL[NV], fR[NV]; pri2con(pL, cL); pri2con(pR, cR); flux(pL, d, fL); flux(pR, d, fR);
	double unL = pL[1+d], unR = pR[1+d], aL = cfast(pL, d), aR = cfast(pR, d);
	double sL = std::min(unL - aL, unR - aR), sR = std::max(unL + aL, unR + aR);
	if(sL >= 0.0){ for(int k=0;k<NV;k++) F[k] = fL[k]; return; }
	if(sR <= 0.0){ for(int k=0;k<NV;k++) F[k] = fR[k]; return; }
	for(int k=0;k<NV;k++) F[k] = (sR*fL[k] - sL*fR[k] + sL*sR*(cR[k]-cL[k]))/(sR - sL);
}
static inline double minmod(double a, double b){ return (a*b <= 0.0) ? 0.0 : (fabs(a) < fabs(b) ? a : b); }
int main(int argc, char** argv){
	int N = argc > 1 ? atoi(argv[1]) : 256; int smooth = argc > 2 ? atoi(argv[2]) : 1; double CFL = argc > 3 ? atof(argv[3]) : 0.4;
	const int ng = 2; const double L = 1.0, dx = L/N, dy = L/N, tend = 0.01, B0 = 70.0/sqrt(4.0*M_PI);
	Grid U(N, N, ng), U1(N, N, ng), U2(N, N, ng); Grid P(N, N, ng); Grid FX(N+1, N, ng), FY(N, N+1, ng);
	// initial condition: 3x3 Gauss cell average of the conservative variables (same as the IMEX code)
	const double gp[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)}, gw[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
	for(int j = 0; j < N; j++) for(int i = 0; i < N; i++){
		double xc = -0.5 + (i+0.5)*dx, yc = -0.5 + (j+0.5)*dy, c[NV] = {0,0,0,0,0,0,0,0};
		for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++){
			double x = xc + 0.5*dx*gp[a], y = yc + 0.5*dy*gp[b], r = sqrt(x*x+y*y), pr;
			if(smooth){ double eps = sqrt(dx*dx+dy*dy); pr = 0.5*(0.1+1000.0) + 0.5*(0.1-1000.0)*erf((r-0.1)/eps); }
			else pr = (r <= 0.1) ? 1000.0 : 0.1;
			double p[NV] = {1.0, 0, 0, 0, pr, B0, 0, 0}, cc[NV]; pri2con(p, cc);
			for(int k=0;k<NV;k++) c[k] += gw[a]*gw[b]*cc[k];
		}
		double* u = U.at(i+ng, j+ng); for(int k=0;k<NV;k++) u[k] = c[k];
	}
	auto bc = [&](Grid& G){   // transmissive
		for(int j = 0; j < N+2*ng; j++) for(int g = 0; g < ng; g++){ for(int k=0;k<NV;k++){ G.at(g,j)[k] = G.at(ng,j)[k]; G.at(N+ng+g,j)[k] = G.at(N+ng-1,j)[k]; } }
		for(int i = 0; i < N+2*ng; i++) for(int g = 0; g < ng; g++){ for(int k=0;k<NV;k++){ G.at(i,g)[k] = G.at(i,ng)[k]; G.at(i,N+ng+g)[k] = G.at(i,N+ng-1)[k]; } }
	};
	auto rhs = [&](Grid& Uc, Grid& Un, double dt){   // Un = Uc - dt * div F + dt * Powell source
		bc(Uc);
		#pragma omp parallel for
		for(int j = 0; j < N+2*ng; j++) for(int i = 0; i < N+2*ng; i++) con2pri(Uc.at(i,j), P.at(i,j));
		// x faces: FX(i,j) is the face between cell i-1 and i (i = ng .. N+ng)
		#pragma omp parallel for
		for(int j = ng; j < N+ng; j++) for(int i = ng; i <= N+ng; i++){
			double pL[NV], pR[NV];
			for(int k=0;k<NV;k++){ double *pm2 = P.at(i-2,j), *pm1 = P.at(i-1,j), *p0 = P.at(i,j), *pp1 = P.at(i+1,j);
				pL[k] = pm1[k] + 0.5*minmod(pm1[k]-pm2[k], p0[k]-pm1[k]); pR[k] = p0[k] - 0.5*minmod(p0[k]-pm1[k], pp1[k]-p0[k]); }
			if(pL[0] <= 0 || pL[4] <= 0){ for(int k=0;k<NV;k++) pL[k] = P.at(i-1,j)[k]; }
			if(pR[0] <= 0 || pR[4] <= 0){ for(int k=0;k<NV;k++) pR[k] = P.at(i,j)[k]; }
			hll(pL, pR, 0, FX.at(i,j));
		}
		#pragma omp parallel for
		for(int j = ng; j <= N+ng; j++) for(int i = ng; i < N+ng; i++){
			double pL[NV], pR[NV];
			for(int k=0;k<NV;k++){ double *pm2 = P.at(i,j-2), *pm1 = P.at(i,j-1), *p0 = P.at(i,j), *pp1 = P.at(i,j+1);
				pL[k] = pm1[k] + 0.5*minmod(pm1[k]-pm2[k], p0[k]-pm1[k]); pR[k] = p0[k] - 0.5*minmod(p0[k]-pm1[k], pp1[k]-p0[k]); }
			if(pL[0] <= 0 || pL[4] <= 0){ for(int k=0;k<NV;k++) pL[k] = P.at(i,j-1)[k]; }
			if(pR[0] <= 0 || pR[4] <= 0){ for(int k=0;k<NV;k++) pR[k] = P.at(i,j)[k]; }
			hll(pL, pR, 1, FY.at(i,j));
		}
		#pragma omp parallel for
		for(int j = ng; j < N+ng; j++) for(int i = ng; i < N+ng; i++){
			double* uc = Uc.at(i,j); double* un = Un.at(i,j); double* p = P.at(i,j);
			double divB = (P.at(i+1,j)[5] - P.at(i-1,j)[5])/(2*dx) + (P.at(i,j+1)[6] - P.at(i,j-1)[6])/(2*dy);
			double uB = p[1]*p[5]+p[2]*p[6]+p[3]*p[7];
			double S[NV] = {0, -p[5]*divB, -p[6]*divB, -p[7]*divB, -uB*divB, -p[1]*divB, -p[2]*divB, -p[3]*divB};   // Powell 8-wave
			for(int k=0;k<NV;k++) un[k] = uc[k] - dt/dx*(FX.at(i+1,j)[k] - FX.at(i,j)[k]) - dt/dy*(FY.at(i,j+1)[k] - FY.at(i,j)[k]) + dt*S[k];
		}
	};
	double t = 0.0; int step = 0;
	while(t < tend){
		bc(U); double smax = 0.0;
		#pragma omp parallel for reduction(max:smax)
		for(int j = ng; j < N+ng; j++) for(int i = ng; i < N+ng; i++){ double p[NV]; con2pri(U.at(i,j), p); smax = std::max(smax, (fabs(p[1]) + cfast(p,0))/dx + (fabs(p[2]) + cfast(p,1))/dy); }
		double dt = std::min(CFL/smax, tend - t);
		rhs(U, U1, dt); rhs(U1, U2, dt);
		#pragma omp parallel for
		for(int j = ng; j < N+ng; j++) for(int i = ng; i < N+ng; i++){ double* u = U.at(i,j); double* u2 = U2.at(i,j); for(int k=0;k<NV;k++) u[k] = 0.5*(u[k] + u2[k]); }
		t += dt; step++;
		if(step % 200 == 0) printf("step %d t=%g dt=%g\n", step, t, dt);
	}
	double pmin = 1e300, rmin = 1e300;
	std::string fn = "ref_" + std::to_string(N) + (smooth ? "" : "_sharp") + ".dat"; FILE* fp = fopen(fn.c_str(), "w");
	for(int i = ng; i < N+ng; i++) for(int j = ng; j < N+ng; j++){ double p[NV]; con2pri(U.at(i,j), p); pmin = std::min(pmin, p[4]); rmin = std::min(rmin, p[0]);
		fprintf(fp, "%.10e %.10e", -0.5+(i-ng+0.5)*dx, -0.5+(j-ng+0.5)*dy); for(int k=0;k<NV;k++) fprintf(fp, " %.10e", p[k]); fprintf(fp, "\n"); }
	fclose(fp); printf("done: %d steps, t=%g, pmin=%g rhomin=%g -> %s\n", step, t, pmin, rmin, fn.c_str());
	return 0;
}
