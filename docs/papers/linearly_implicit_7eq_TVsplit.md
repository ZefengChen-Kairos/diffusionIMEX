# Reading note: A linearly implicit shock capturing scheme for compressible two-phase flows at all Mach numbers

- **Authors:** Beatrice Battisti (LAMA, CNRS UMR 5127, Université Savoie Mont Blanc), Walter Boscheri (LAMA + Dept. Mathematics and Computer Science, University of Ferrara)
- **Venue:** *Journal of Computational Physics* 539 (2025) 114227, doi:10.1016/j.jcp.2025.114227 (open access, CC BY). Received 23 Feb 2025, accepted 6 Jul 2025.
- **Keywords:** all Mach solver, compressible two-phase flows, semi-implicit IMEX, asymptotic preserving, well-balanced, Baer–Nunziato.

> Scope of this note: equation-level transcription made for implementing a 6-equation diffuse-interface model with a transport/pressure-split IMEX scheme. Equation numbers are the paper's. Anything I could not read with certainty, or that looks like a typo in the paper, is marked **[uncertain]**. My own comments (not in the paper) are marked **[note]**.

**Big picture in one paragraph.** The paper solves the full 7-equation Baer–Nunziato (BN) model with a *semi-implicit, linearly implicit* scheme in the Boscheri–Pareschi (2021, "high order pressure-based semi-implicit IMEX schemes for the 3D Navier–Stokes equations at all Mach numbers") style, i.e. the splitting idea of Toro–Vázquez-Cendón is used for the **energy flux** (kinetic part explicit, enthalpy part implicit) and pressure gradients are implicit. The convective (nonlinear) fluxes are explicit, discretised by a finite-volume Rusanov / global Lax–Friedrichs flux whose dissipation uses **only material velocities**; the implicit pressure terms are discretised by **central finite differences** with no added numerical diffusion. The volume-fraction equation is solved **first and implicitly** (linear system, GMRES). Then a linear **wave equation for $(\phi_\alpha p_\alpha)^{n+1}$** per phase is solved by GMRES, followed by explicit-in-form updates of momentum and energy. Relaxation sources (pressure relaxation $\mu$, drag $\delta$) are handled afterwards by a Strang-type splitting with a semi-implicit cell-local linear solve. Second order in time: LSDIRK2(2,2,2) semi-implicit IMEX Runge–Kutta of Boscaroli–Russo style [41]. CFL is based on $\max(|u_f|,|u_s|)$ only. Proven AP (low-Mach mixture limit gives $\nabla\cdot(\phi_s u_s+\phi_f u_f)=0$) and well-balanced (preserves moving contact discontinuities with constant $p,u$).

**[note] Important clarification vs. the request's framing.** The paper does *not* split the flux into a "transport subsystem" and a "pressure subsystem" each with its own eigenstructure / Riemann solver in the classical Toro–Vázquez-Cendón (TV) sense. Instead the split is at the semi-discrete-in-time level: which terms are evaluated at $t^n$ (explicit) vs $t^{n+1}$ (implicit). The only eigenvalue information used is that the explicit convective part has eigenvalues $u_\alpha\cdot n$ (material velocity), used both for the CFL and for the Rusanov dissipation coefficient. The TV-type ingredient is the energy-flux split (9). Details below.

---

## 1. Governing equations (Section 2)

Domain $\Omega\subset\mathbb{R}^d$, $d\in\{1,2,3\}$, $\mathbf x=(x,y,z)$. Phases $\alpha\in\{f,s\}$ (f = fluid, s = solid). Saurel–Abgrall [15] form of BN:

$$
\begin{aligned}
&\partial_t(\rho_f\phi_f)+\nabla\cdot(\rho_f\phi_f\mathbf u_f)=0,\\
&\partial_t(\rho_f\phi_f\mathbf u_f)+\nabla\cdot(\rho_f\phi_f\mathbf u_f\otimes\mathbf u_f)+\nabla(p_f\phi_f)+p_I\nabla\phi_f=\delta(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_f\phi_f E_f)+\nabla\cdot(\phi_f\mathbf u_f(\rho_fE_f+p_f))+p_I\mathbf u_I\cdot\nabla\phi_f=\mu\,p_I(p_s-p_f)+\delta\,\mathbf u_I\cdot(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_s\phi_s)+\nabla\cdot(\rho_s\phi_s\mathbf u_s)=0,\\
&\partial_t(\rho_s\phi_s\mathbf u_s)+\nabla\cdot(\rho_s\phi_s\mathbf u_s\otimes\mathbf u_s)+\nabla(p_s\phi_s)-p_I\nabla\phi_s=-\delta(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_s\phi_s E_s)+\nabla\cdot(\phi_s\mathbf u_s(\rho_sE_s+p_s))-p_I\mathbf u_I\cdot\nabla\phi_s=-\mu\,p_I(p_s-p_f)-\delta\,\mathbf u_I\cdot(\mathbf u_s-\mathbf u_f),\\
&\partial_t\phi_s+\mathbf u_I\cdot\nabla\phi_s=\mu(p_s-p_f).
\end{aligned}\tag{1}
$$

**[uncertain / typo in paper]** The momentum equations are printed with a "$+p_I\nabla\phi_f$" on the LHS of the fluid momentum equation (written as $+\,p_I\nabla\phi_f$). Since $\nabla\phi_f=-\nabla\phi_s$ the fluid and solid interface forces are $-p_I\nabla\phi_s$ and $+p_I\nabla\phi_s$ on the RHS respectively... careful: with the printed signs the fluid momentum has LHS $+p_I\nabla\phi_f$, which is the *opposite* of the standard BN convention ($\partial_t(\alpha\rho u)+\dots+\nabla(\alpha p) = p_I\nabla\alpha$, i.e. LHS $-p_I\nabla\alpha$). Looking at the image, the fluid momentum line reads "$\nabla(p_f\phi_f)+p_I\nabla\phi_f$" and fluid energy "$+p_I\mathbf u_I\nabla\phi_f$", and the solid lines read "$-p_I\nabla\phi_s$", "$-p_I\mathbf u_I\nabla\phi_s$". Since $\nabla\phi_f=-\nabla\phi_s$, both phases then carry $-p_I\nabla\phi_s$ on the LHS for fluid (= $+p_I\nabla\phi_f$) and $-p_I\nabla\phi_s$ for solid — which means the interfacial forces do **not** cancel in the sum. The standard (momentum-conserving) BN is: fluid LHS $-p_I\nabla\phi_f$, solid LHS $-p_I\nabla\phi_s$, i.e. fluid LHS $= +p_I\nabla\phi_s$. Most likely the paper intends the standard BN and the fluid sign is a typo, or they intend "$\nabla\phi_f$" to be "$\nabla\phi_s$". Implement the standard, mixture-conservative form. Also the generic discrete form (12b) uses "$+\Delta t\,p_I^n\nabla\phi_\alpha^{n+1}$" for both phases, which is the standard form $\partial_t(\phi\rho u)+\dots=p_I\nabla\phi_\alpha$. So the generic scheme (12) is consistent with the standard BN; use that.

- Variables: phase density $\rho_\alpha$, pressure $p_\alpha$, velocity $\mathbf u_\alpha$, volume fraction $\phi_\alpha\in[0,1]$, $\phi_s+\phi_f=1$; specific total energy $E_\alpha=e_\alpha+k_\alpha$, $k_\alpha=\tfrac12|\mathbf u_\alpha|^2$.
- Sources: $\delta$ = drag (velocity relaxation), $\mu$ = pressure relaxation parameter; both possibly stiff ("frequently treated as instantaneous", $\mu,\delta\to\infty$).

### EOS (stiffened gas)

$$
e_\alpha=\frac{p_\alpha+\gamma_\alpha\pi_\alpha}{\rho_\alpha(\gamma_\alpha-1)},\tag{2}
$$
ideal gas for $\pi_\alpha=0$. Sound speed $c_\alpha=\sqrt{\gamma_\alpha(p_\alpha+\pi_\alpha)/\rho_\alpha}$.

**[note]** The key property used later: $\rho_\alpha\phi_\alpha e_\alpha=\frac{\phi_\alpha p_\alpha}{\gamma_\alpha-1}+\frac{\phi_\alpha\gamma_\alpha\pi_\alpha}{\gamma_\alpha-1}$ is *linear* in the unknown $(\phi_\alpha p_\alpha)$ once $\phi_\alpha^{n+1}$ is known. Hence the pressure system is linear.

### Interface closures (both allowed)

$$
p_I=p_f,\qquad \mathbf u_I=\mathbf u_s \tag{3a}
$$
(original BN [1], also [54,55]),
$$
p_I=\sum_\alpha\phi_\alpha p_\alpha,\qquad \mathbf u_I=\frac{\sum_\alpha\rho_\alpha\phi_\alpha\mathbf u_\alpha}{\sum_\alpha\rho_\alpha\phi_\alpha} \tag{3b}
$$
(mixture-weighted; [13,15,46,51]). Choice does not affect the numerical method.

### Eigenvalues (normal direction $\mathbf n$, from Andrianov–Warnecke [56])

$$
\lambda=\left[\mathbf u_f\cdot\mathbf n,\ \mathbf u_f\cdot\mathbf n\pm c_f,\ \mathbf u_s\cdot\mathbf n,\ \mathbf u_s\cdot\mathbf n\pm c_s,\ \mathbf u_I\cdot\mathbf n\right].\tag{4}
$$

### Compact form

$$
\partial_t\mathbf Q+\nabla\cdot\mathbf F(\mathbf Q)+\mathbf B(\mathbf Q)\cdot\nabla\mathbf Q=\mathbf S(\mathbf Q),\tag{5}
$$
$\mathbf Q=(\rho_f\phi_f,\rho_f\phi_f\mathbf u_f,\rho_f\phi_fE_f,\rho_s\phi_s,\rho_s\phi_s\mathbf u_s,\rho_s\phi_sE_s,\phi_s)^\top$.

Explicit CFL (what they want to avoid):
$$
\Delta t\le\mathrm{CFL}\min_\Omega\frac{\varrho}{\max(|\mathbf u_\alpha\cdot\mathbf n|\pm c_\alpha)},\tag{6}
$$
$\varrho$ = characteristic mesh size.

---

## 2. Non-dimensional form (Section 2.1)

Scalings with common reference values for both phases:
$\tilde{\mathbf x}=\mathbf x/L_0$, $\tilde{\mathbf u}=\mathbf u/\mathbf u_0$, $\tilde t=t/\tau_0=|\mathbf u_0|t/L_0$, $\tilde\rho=\rho/\rho_0$, $\tilde p=p/p_0$, $\tilde E=\rho_0E/p_0$ [as printed], $\tilde\delta=\tau_0\delta/\rho_0$, $\tilde\mu=\rho_0\mathbf u_0^2\tau_0\mu$. Strouhal $St=L_0/(\tau_0|\mathbf u_0|)\approx1$.

The dimensional-coefficient system is (7) (not transcribed term by term; it is (1) with each term multiplied by $\rho_0/\tau_0$, $\rho_0u_0/L_0$, $p_0/L_0$, etc.). With stiffness parameter
$$
\epsilon=\gamma_0M^2,\qquad M=|\mathbf u_0|/c_0,\quad c_0=\sqrt{\gamma_0p_0/\rho_0},
$$
dropping tildes:
$$
\begin{aligned}
&\partial_t(\rho_f\phi_f)+\nabla\cdot(\rho_f\phi_f\mathbf u_f)=0,\\
&\partial_t(\rho_f\phi_f\mathbf u_f)+\nabla\cdot(\rho_f\phi_f\mathbf u_f\otimes\mathbf u_f)+\tfrac1\epsilon\nabla(p_f\phi_f)+\tfrac1\epsilon p_I\nabla\phi_f=\delta(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_f\phi_fE_f)+\nabla\cdot(\phi_f\mathbf u_f(\rho_fE_f+p_f))+p_I\mathbf u_I\nabla\phi_f=\tfrac1\epsilon\mu p_I(p_s-p_f)+\epsilon\delta\,\mathbf u_I(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_s\phi_s)+\nabla\cdot(\rho_s\phi_s\mathbf u_s)=0,\\
&\partial_t(\rho_s\phi_s\mathbf u_s)+\nabla\cdot(\rho_s\phi_s\mathbf u_s\otimes\mathbf u_s)+\tfrac1\epsilon\nabla(p_s\phi_s)-\tfrac1\epsilon p_I\nabla\phi_s=-\delta(\mathbf u_s-\mathbf u_f),\\
&\partial_t(\rho_s\phi_sE_s)+\nabla\cdot(\phi_s\mathbf u_s(\rho_sE_s+p_s))-p_I\mathbf u_I\nabla\phi_s=-\tfrac1\epsilon\mu p_I(p_s-p_f)-\epsilon\delta\,\mathbf u_I(\mathbf u_s-\mathbf u_f),\\
&\partial_t\phi_s+\mathbf u_I\cdot\nabla\phi_s=\tfrac1\epsilon\mu(p_s-p_f).
\end{aligned}\tag{8}
$$
Here $E_\alpha$ is the non-dimensional total energy with $E_\alpha=e_\alpha+\epsilon k_\alpha$ implied by (9) **[uncertain: the paper does not write $E=e+\epsilon k$ explicitly, but (9) implies it]**.

### Energy-flux split (the Toro–Vázquez-Cendón / Boscheri–Pareschi ingredient)

Following [43]:
$$
\phi_\alpha\mathbf u_\alpha(\rho_\alpha E_\alpha+p_\alpha)=\epsilon\,\rho_\alpha\phi_\alpha k_\alpha\mathbf u_\alpha+\rho_\alpha\phi_\alpha h_\alpha\mathbf u_\alpha,\tag{9}
$$
$$
h_\alpha=e_\alpha+\frac{p_\alpha}{\rho_\alpha}.\tag{10}
$$
Kinetic part = convective (explicit); enthalpy part = "pressure" part (implicit, with $h$ frozen at $t^n$ and the mass flux $\rho\phi\mathbf u$ at $t^{n+1}$).

Stiffness: $1/\epsilon$ in momentum (acoustic waves) and in energy through the pressure work. Single global Mach number for the mixture (as in Re–Abgrall [46]).

---

## 3. Numerical scheme (Section 3)

Overall: solve homogeneous system (S=0) over $\Delta t$; then a relaxation operator (ODE, implicit). Time step from
$$
\Delta t\le\mathrm{CFL}\min_\Omega\frac{\varrho}{\max(\max(|\mathbf u_f\cdot\mathbf n|,|\mathbf u_s\cdot\mathbf n|))},\tag{11}
$$
**material speeds only**, independent of $\epsilon$.

### 3.1 First-order semi-discrete scheme in time (homogeneous part, dimensional form, $\epsilon$ absorbed)

$$
\begin{aligned}
(\rho_\alpha\phi_\alpha)^{n+1}&=(\rho_\alpha\phi_\alpha)^n-\Delta t\,\nabla\cdot\big((\rho_\alpha\mathbf u_\alpha)^n\phi_\alpha^{n+1}\big), &(12a)\\
(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1}&=(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n-\Delta t\,\nabla\cdot\big((\rho_\alpha\mathbf u_\alpha\otimes\mathbf u_\alpha)^n\phi_\alpha^{n+1}\big)-\Delta t\,\nabla(\phi_\alpha p_\alpha)^{n+1}+\Delta t\,p_I^n\nabla\phi_\alpha^{n+1}, &(12b)\\
(\rho_\alpha\phi_\alpha E_\alpha)^{n+1}&=(\rho_\alpha\phi_\alpha E_\alpha)^n-\Delta t\,\nabla\cdot\big((\rho_\alpha k_\alpha\mathbf u_\alpha)^n\phi_\alpha^{n+1}\big)-\Delta t\,\nabla\cdot\big(h_\alpha^n(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1}\big)+\Delta t\,p_I^n\mathbf u_I^n\cdot\nabla\phi_\alpha^{n+1}, &(12c)\\
\phi_s^{n+1}&=\phi_s^n-\Delta t\,\mathbf u_I^n\cdot\nabla\phi_s^{n+1}. &(12d)
\end{aligned}
$$

Explicit/implicit assignment (this is *the* split):

| Term | Time level | Comment |
|---|---|---|
| Mass flux $\rho_\alpha\mathbf u_\alpha\phi_\alpha$ | $(\rho_\alpha\mathbf u_\alpha)^n\,\phi_\alpha^{n+1}$ | explicit except $\phi$ (already known) |
| Momentum convective flux $\rho\mathbf u\otimes\mathbf u\,\phi$ | $(\rho\mathbf u\otimes\mathbf u)^n\phi^{n+1}$ | explicit |
| Pressure gradient $\nabla(\phi_\alpha p_\alpha)$ | $n+1$ | **implicit** |
| Interface force $p_I\nabla\phi_\alpha$ | $p_I^n\nabla\phi_\alpha^{n+1}$ | linearised (coefficient explicit) |
| Kinetic energy flux $\rho k\mathbf u\phi$ | $(\rho k\mathbf u)^n\phi^{n+1}$ | explicit |
| Enthalpy flux $\rho\phi h\mathbf u$ | $h^n(\rho\phi\mathbf u)^{n+1}$ | **implicit** in momentum, $h$ frozen |
| Interface work $p_I\mathbf u_I\cdot\nabla\phi$ | $p_I^n\mathbf u_I^n\cdot\nabla\phi^{n+1}$ | linearised |
| Volume fraction transport $\mathbf u_I\cdot\nabla\phi_s$ | $\mathbf u_I^n\cdot\nabla\phi_s^{n+1}$ | **implicit** (needed for well-balancing, Sec. 3.7) |

Interface quantities and the enthalpy are always explicit → scheme is **linearly implicit** (no nonlinearity in time). Volume fraction (12d) is solved first, so $\phi^{n+1}$ is then available to all other equations.

Explicit abbreviations:
$$
\rho_\alpha^\star=(\rho_\alpha\phi_\alpha)^n-\Delta t\nabla\cdot((\rho_\alpha\mathbf u_\alpha)^n\phi_\alpha^{n+1}),\tag{13a}
$$
$$
\mathbf q_\alpha^\star=(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n-\Delta t\nabla\cdot((\rho_\alpha\mathbf u_\alpha\otimes\mathbf u_\alpha)^n\phi_\alpha^{n+1}),\tag{13b}
$$
$$
E_\alpha^\star=(\rho_\alpha\phi_\alpha E_\alpha)^n-\Delta t\nabla\cdot((\rho_\alpha k_\alpha\mathbf u_\alpha)^n\phi_\alpha^{n+1}).\tag{13c}
$$
Then
$$
(\rho_\alpha\phi_\alpha)^{n+1}=\rho_\alpha^\star,\tag{14a}
$$
$$
(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1}=\mathbf q_\alpha^\star-\Delta t\nabla(\phi_\alpha p_\alpha)^{n+1}+\Delta t\,p_I^n\nabla\phi_\alpha^{n+1},\tag{14b}
$$
$$
(\rho_\alpha\phi_\alpha e_\alpha)^{n+1}+(\rho_\alpha\phi_\alpha k_\alpha)^{n+1}=E_\alpha^\star-\Delta t\nabla\cdot(h_\alpha^n(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1})+\Delta t\,p_I^n\mathbf u_I^n\nabla\phi_\alpha^{n+1}.\tag{14c}
$$

Semi-implicit linearisation of kinetic energy (Boscheri–Pareschi [43]):
$$
(\rho_\alpha\phi_\alpha k_\alpha)^{n+1}:=\frac12\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}\cdot(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1}.\tag{15}
$$

### Pressure (wave) equation

Insert (14b) into (14c), use EOS (2) for $\rho\phi e$:
$$
\frac{(\phi_\alpha p_\alpha)^{n+1}}{\gamma_\alpha-1}+\frac{\phi_\alpha^{n+1}\gamma_\alpha\pi_\alpha}{\gamma_\alpha-1}+\frac12\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}\cdot\Big(\mathbf q_\alpha^\star-\Delta t\nabla(\phi_\alpha p_\alpha)^{n+1}+\Delta t\,p_I^n\nabla\phi_\alpha^{n+1}\Big)
$$
$$
=E_\alpha^\star-\Delta t\nabla\cdot(h_\alpha^n\mathbf q_\alpha^\star)+\Delta t^2\nabla\cdot\big(h_\alpha^n\nabla(\phi_\alpha p_\alpha)^{n+1}\big)-\Delta t^2\nabla\cdot\big(h_\alpha^np_I^n\nabla\phi_\alpha^{n+1}\big)+\Delta t\,p_I^n\mathbf u_I^n\cdot\nabla\phi_\alpha^{n+1}.\tag{16}
$$
Unknowns to LHS:
$$
\frac{(\phi_\alpha p_\alpha)^{n+1}}{\gamma_\alpha-1}+\frac{\phi_\alpha^{n+1}\gamma_\alpha\pi_\alpha}{\gamma_\alpha-1}-\frac{\Delta t}{2}\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}\cdot\nabla(\phi_\alpha p_\alpha)^{n+1}-\Delta t^2\nabla\cdot\big(h_\alpha^n\nabla(\phi_\alpha p_\alpha)^{n+1}\big)
$$
$$
=E_\alpha^\star-\frac12\mathbf q_\alpha^\star\cdot\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}-\frac{\Delta t}{2}p_I^n\nabla\phi_\alpha^{n+1}\cdot\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}-\Delta t\nabla\cdot(h_\alpha^n\mathbf q_\alpha^\star)-\Delta t^2\nabla\cdot\big(h_\alpha^np_I^n\nabla\phi_\alpha^{n+1}\big)+\Delta t\,p_I^n\mathbf u_I^n\cdot\nabla\phi_\alpha^{n+1}.\tag{17}
$$

Remarks in the paper:
- Unknown is $(\phi_\alpha p_\alpha)^{n+1}$ (not $p_\alpha$), one decoupled linear system **per phase** (the phases couple only via $\phi^{n+1}$, $p_I^n$, $\mathbf u_I^n$ and later via relaxation).
- Stiffened gas → linear relation $\rho e\leftrightarrow p$ → linear system. For general EOS, nonlinearity is only on the diagonal; use nested Newton (Brugnano–Casulli [58]).
- The pressure from (17) is used only as a **stable numerical flux** for the *conservative* updates of momentum and total energy (14b), (14c), unlike non-conservative pressure updates (Re–Abgrall [46]).
- Enthalpy must be discretised as
$$
h_\alpha^n:=\frac{\rho_\alpha^n h_\alpha^n}{\rho_\alpha^{n+1}}\tag{18}
$$
to preserve constant velocity and pressure (well-balancing), as in [43]. **[note]** I.e. use $(\rho h)^n/\rho^{n+1}$ where $\rho^{n+1}$ is known from (14a); with $\phi$ factored it should be $(\rho_\alpha\phi_\alpha h_\alpha)^n/(\rho_\alpha\phi_\alpha)^{n+1}$ **[uncertain which; printed as $\rho^n_\alpha h^n_\alpha/\rho^{n+1}_\alpha$]**.

**Remark 1 (thermodynamic consistency).** After (17), total energy is updated by (14c); the resulting total energy equals
$$
(\rho_\alpha\phi_\alpha E_\alpha)^{n+1}=\frac{(\phi_\alpha p_\alpha)^{n+1}}{\gamma_\alpha-1}+\frac{\phi_\alpha^{n+1}\gamma_\alpha\pi_\alpha}{\gamma_\alpha-1}+\frac12\frac{(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^n}{(\rho_\alpha\phi_\alpha)^n}\cdot(\rho_\alpha\phi_\alpha\mathbf u_\alpha)^{n+1},
$$
i.e. thermodynamic inconsistency of $\mathcal O(\Delta t^q)$ ($q$ = order). Fixes: (i) drop the linearisation (15) and do Picard/fixed-point iterations [42,59]; (ii) Boscheri et al. [60]: solve pressure system twice — first with linearised kinetic energy to get a *pressure flux*, update momentum by (14b), then solve the pressure system again for the *pressure state* using the now-known $k^{n+1}$. Two linear solves per stage → strict consistency.

### 3.2 Second order in time: semi-implicit IMEX-RK (Boscarino–Filbet–Russo [41])

Autonomous partitioned form:
$$
\partial_t\mathbf Q=\mathcal H(\mathbf Q_E(t),\mathbf Q_I(t)),\qquad\mathbf Q_0=\mathbf Q(t=t_0),\tag{19}
$$
$\mathbf Q_E$ = arguments treated explicitly, $\mathbf Q_I$ = implicitly. The first-order scheme (12) is an instance. $s$-stage Butcher tableaux
$$
\begin{array}{c|c}\hat c&\hat A\\\hline&\hat b^T\end{array}\qquad\begin{array}{c|c}c&A\\\hline&b^T\end{array}\tag{20}
$$
($\hat A,A\in\mathbb R^{s\times s}$), left explicit, right implicit. Used: stiffly accurate **LSDIRK2(2,2,2)**:
$$
\begin{array}{c|cc}0&0&0\\\hat c&\hat c&0\\\hline&1-\beta&\beta\end{array}\qquad
\begin{array}{c|cc}\beta&\beta&0\\1&1-\beta&\beta\\\hline&1-\beta&\beta\end{array}\tag{21}
$$
$\beta=1-1/\sqrt2$, $\hat c=1/(2\beta)$.

**[note] How the semi-implicit IMEX of [41] works (not re-derived in this paper; summary from Boscarino–Filbet–Russo / Boscheri–Pareschi):** at stage $i$, explicit arguments are $\mathbf Q_E^{(i)}=\mathbf Q^n+\Delta t\sum_{j<i}\hat a_{ij}\mathcal H(\mathbf Q_E^{(j)},\mathbf Q_I^{(j)})$ and the implicit stage solves $\mathbf Q_I^{(i)}=\mathbf Q^n+\Delta t\sum_{j<i}a_{ij}\mathcal H(\mathbf Q_E^{(j)},\mathbf Q_I^{(j)})+\Delta t\,a_{ii}\mathcal H(\mathbf Q_E^{(i)},\mathbf Q_I^{(i)})$, which is linear in $\mathbf Q_I^{(i)}$ — each stage is one "first-order-like" step (12)–(17) with effective time step $a_{ii}\Delta t$ and modified explicit RHS. Stiffly accurate: $\mathbf Q^{n+1}=\mathbf Q_I^{(s)}$. Note $\hat c=1/(2\beta)\approx1.707>1$ — the second explicit stage is evaluated outside $[t^n,t^{n+1}]$ in the explicit tableau. The paper does not spell out the stage-level formulas; see [41,43]. **[uncertain: exact stage algebra in this paper]**.

### 3.3 Fully discrete scheme (Cartesian grid)

Domain $[x_{\min},x_{\max}]\times[y_{\min},y_{\max}]\times[z_{\min},z_{\max}]$,
$$
\Delta x=\frac{x_{\max}-x_{\min}}{N_x},\ \Delta y=\frac{y_{\max}-y_{\min}}{N_y},\ \Delta z=\frac{z_{\max}-z_{\min}}{N_z},\tag{22}
$$
$N_e=N_xN_yN_z$, cell $C_{i,j,k}$, faces at $i\pm1/2$ etc. Presented in 1D ($x$), extension dimension-by-dimension.

#### Volume fraction (solved first, implicit, fluctuation form)
$$
\phi_{s,i}^{n+1}+\frac{\Delta t}{2\Delta x}\Big(u^n_{I,i+1/2}(\phi^{n+1}_{s,i+1}-\phi^{n+1}_{s,i})-u^n_{I,i-1/2}(\phi^{n+1}_{s,i-1}-\phi^{n+1}_{s,i})\Big)=\phi^n_{s,i},\tag{23}
$$
**[uncertain sign]** as printed the second fluctuation is $(\phi_{s,i-1}-\phi_{s,i})$ with a minus in front, giving $u_{I,i-1/2}(\phi_i-\phi_{i-1})$, i.e. the central approximation $\frac{1}{2\Delta x}[u_{i+1/2}(\phi_{i+1}-\phi_i)+u_{i-1/2}(\phi_i-\phi_{i-1})]\approx u\,\partial_x\phi$. Consistent.
$$
u^n_{I,i+1/2}=\tfrac12(u^n_{I,i+1}+u^n_{I,i}),\qquad u^n_{I,i-1/2}=\tfrac12(u^n_{I,i}+u^n_{I,i-1}),\tag{24}
$$
with (3a) or (3b). Preserves constant $p,u$ across a traveling contact [46]; equivalent to a path-conservative scheme (Parés [61]) with straight-line path and midpoint rule. Pure central → no upwinding (dissipation added optionally, (36)). Linear system solved by **GMRES**, tolerance $10^{-12}$.

#### Explicit convective fluxes (finite volume, Rusanov / global LF)
$$
m_i^\star=m_i^n-\frac{\Delta t}{\Delta x}\big(f^m_{i+1/2}-f^m_{i-1/2}\big),\tag{25}
$$
$$
f^m_{i+1/2}=\tfrac12\big(f(m^n_{i+1})+f(m^n_i)\big)-\tfrac12a^n_\Omega(m^n_{i+1}-m^n_i),\qquad
f^m_{i-1/2}=\tfrac12\big(f(m^n_i)+f(m^n_{i-1})\big)-\tfrac12a^n_\Omega(m^n_i-m^n_{i-1}),\tag{26}
$$
$$
a^n_\Omega=\max_\Omega\big(\max(|u^n_f|,|u^n_s|)\big),\tag{27}
$$
**global** LF with **only convective eigenvalues** (no sound speed). Energy dissipation acts on the jump of the **kinetic** energy $\rho\phi k$, not total energy (so internal energy is not diffused).

Note: the physical flux is $f(m)$ with $\phi^{n+1}$ inside, e.g. $f^{\rho_\alpha\phi_\alpha}=(\rho_\alpha u_\alpha)^n\phi^{n+1}_\alpha$ **[uncertain: whether the dissipation jump uses $(\rho\phi)^n$ — it is printed $m^n$, i.e. yes the old conserved variable]**.

$$
\rho^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha)^n_i-\frac{\Delta t}{\Delta x}\big(f^{\rho_\alpha\phi_\alpha}_{i+1/2}-f^{\rho_\alpha\phi_\alpha}_{i-1/2}\big)=(\rho_\alpha\phi_\alpha)^{n+1}_i,\tag{28a}
$$
$$
q^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha u_\alpha)^n_i-\frac{\Delta t}{\Delta x}\big(f^{\rho_\alpha\phi_\alpha u_\alpha}_{i+1/2}-f^{\rho_\alpha\phi_\alpha u_\alpha}_{i-1/2}\big),\tag{28b}
$$
$$
E^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha E_\alpha)^n_i-\frac{\Delta t}{\Delta x}\big(f^{\rho_\alpha\phi_\alpha k_\alpha}_{i+1/2}-f^{\rho_\alpha\phi_\alpha k_\alpha}_{i-1/2}\big).\tag{28c}
$$

#### Implicit terms: central finite differences — discrete pressure equation
$$
\frac{(\phi_\alpha p_\alpha)^{n+1}_i}{\gamma_\alpha-1}+\frac{\phi^{n+1}_{\alpha,i}\gamma_\alpha\pi_\alpha}{\gamma_\alpha-1}-\frac{\Delta t}{4\Delta x}\frac{(\rho_\alpha\phi_\alpha u_\alpha)^n_i}{(\rho_\alpha\phi_\alpha)^n_i}\Big((\phi_\alpha p_\alpha)^{n+1}_{i+1}-(\phi_\alpha p_\alpha)^{n+1}_{i-1}\Big)-\Delta t^2\nabla\cdot\big(h^n_\alpha\nabla(\phi_\alpha p_\alpha)^{n+1}\big)
$$
$$
=E^\star_{\alpha,i}-\frac12q^\star_{\alpha,i}\frac{(\rho_\alpha\phi_\alpha u_\alpha)^n_i}{(\rho_\alpha\phi_\alpha)^n_i}-\frac{\Delta t}{4\Delta x}\Big(p^n_{I,i+1/2}(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i})-p^n_{I,i-1/2}(\phi^{n+1}_{\alpha,i-1}-\phi^{n+1}_{\alpha,i})\Big)\frac{(\rho_\alpha\phi_\alpha u_\alpha)^n_i}{(\rho_\alpha\phi_\alpha)^n_i}
$$
$$
-\frac{\Delta t}{2\Delta x}\Big((h^n_\alpha q^\star_\alpha)_{i+1}-(h^n_\alpha q^\star_\alpha)_{i-1}\Big)-\Delta t^2\nabla\cdot\big(h^n_\alpha p^n_I\nabla\phi^{n+1}_\alpha\big)
+\frac{\Delta t}{2\Delta x}\Big(p^n_{I,i+1/2}u^n_{I,i+1/2}(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i})-p^n_{I,i-1/2}u^n_{I,i-1/2}(\phi^{n+1}_{\alpha,i-1}-\phi^{n+1}_{\alpha,i})\Big).\tag{29}
$$
- Non-conservative $\nabla\phi^{n+1}$ terms use the same fluctuation form as (23) (with interface values $p_{I,i\pm1/2}$, $(p_Iu_I)_{i\pm1/2}$ — presumably arithmetic averages as in (24) **[uncertain: averaging of $p_I$ not stated explicitly]**). Note the factor $\frac{\Delta t}{4\Delta x}$ in the kinetic term = $\frac12\cdot\frac{\Delta t}{2\Delta x}$.
- **[note]** Sign oddity: for the fluctuation in (29) the printed form is $p_{I,i+1/2}(\phi_{i+1}-\phi_i)-p_{I,i-1/2}(\phi_{i-1}-\phi_i)$, consistent with (23).

The Laplacian-type term. Naive consistent discretisation
$$
\Delta t^2\nabla\cdot(h^n_\alpha\nabla(\phi_\alpha p_\alpha)^{n+1})=\frac{\Delta t^2}{2\Delta x}\left(h^n_{\alpha,i+1}\frac{(\phi_\alpha p_\alpha)^{n+1}_{i+2}-(\phi_\alpha p_\alpha)^{n+1}_i}{2\Delta x}-h^n_{\alpha,i-1}\frac{(\phi_\alpha p_\alpha)^{n+1}_i-(\phi_\alpha p_\alpha)^{n+1}_{i-2}}{2\Delta x}\right)\tag{30}
$$
has wide stencil $\{i-2,\dots,i+2\}$ (and odd-even decoupling **[note]**). Instead, compact form from second-degree Lagrange interpolation on $\{i-1,i,i+1\}$ ([43]):
$$
\Delta t^2\nabla\cdot(h^n_\alpha\nabla(\phi_\alpha p_\alpha)^{n+1})=\frac{\Delta t^2}{\Delta x^2}\begin{bmatrix}h^n_{\alpha,i-1}&h^n_{\alpha,i}&h^n_{\alpha,i+1}\end{bmatrix}
\begin{bmatrix}3/4&-1&1/4\\0&0&0\\1/4&-1&3/4\end{bmatrix}
\begin{bmatrix}(\phi_\alpha p_\alpha)^{n+1}_{i-1}\\(\phi_\alpha p_\alpha)^{n+1}_i\\(\phi_\alpha p_\alpha)^{n+1}_{i+1}\end{bmatrix}+\mathcal O(\Delta x^2).\tag{31}
$$
**[note]** Expanding: $\frac{\Delta t^2}{\Delta x^2}\big[h_{i-1}(\tfrac34P_{i-1}-P_i+\tfrac14P_{i+1})+h_{i+1}(\tfrac14P_{i-1}-P_i+\tfrac34P_{i+1})\big]$ with $P=\phi_\alpha p_\alpha$; for constant $h$ this is $h(P_{i-1}-2P_i+P_{i+1})/\Delta x^2$; $h_i$ drops out (middle row zero). The matrix is tridiagonal in 1D (2d+1 point stencil in d dimensions, dimension by dimension), **non-symmetric** (convective term $\frac{\Delta t}{4\Delta x}u(\cdot)$ and the weights), hence GMRES. Sanity check of row: coefficients of $P_{i-1}$: $\tfrac34h_{i-1}+\tfrac14h_{i+1}$; $P_i$: $-(h_{i-1}+h_{i+1})$; $P_{i+1}$: $\tfrac14h_{i-1}+\tfrac34h_{i+1}$. Row sum zero → constant $P$ is in the kernel (well-balanced).

The term $\Delta t^2\nabla\cdot(h^n_\alpha p^n_I\nabla\phi^{n+1}_\alpha)$ uses the same operator (31) (with $h\to h p_I$ and $P\to\phi$) **[uncertain: whether $h p_I$ is the nodal product in the weight vector]**; proven to preserve $p,u$ across contacts (Sec. 3.7).

The system (29) is solved directly for $(\phi_\alpha p_\alpha)^{n+1}$ via **GMRES**. Then:
$$
(\rho_\alpha\phi_\alpha u_\alpha)^{n+1}_i=q^\star_{\alpha,i}-\frac{\Delta t}{2\Delta x}\big((\phi_\alpha p_\alpha)^{n+1}_{i+1}-(\phi_\alpha p_\alpha)^{n+1}_{i-1}\big)+\frac{\Delta t}{2\Delta x}\Big(p^n_{I,i+1/2}(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i})-p^n_{I,i-1/2}(\phi^{n+1}_{\alpha,i-1}-\phi^{n+1}_{\alpha,i})\Big),\tag{32a}
$$
$$
(\rho_\alpha\phi_\alpha E_\alpha)^{n+1}_i=E^\star_{\alpha,i}-\frac{\Delta t}{2\Delta x}\Big((h^n_\alpha(\rho_\alpha\phi_\alpha u_\alpha)^{n+1})_{i+1}-(h^n_\alpha(\rho_\alpha\phi_\alpha u_\alpha)^{n+1})_{i-1}\Big)+\frac{\Delta t}{2\Delta x}\Big(p^n_{I,i+1/2}u^n_{I,i+1/2}(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i})-p^n_{I,i-1/2}u^n_{I,i-1/2}(\phi^{n+1}_{\alpha,i-1}-\phi^{n+1}_{\alpha,i})\Big).\tag{32b}
$$
**[note]** (32a)'s pressure gradient is a central difference $\pm1$ (span 2Δx), but the Laplacian in (29) uses compact (31) — so (29) is *not* literally "(32a) substituted into (14c)"; the discrete energy update (32b) uses the central difference of $h(\rho\phi u)^{n+1}$, which composes to the wide stencil (30). Hence the energy computed by (32b) is not exactly the one implied by (29). This mismatch is part of the $\mathcal O(\Delta t^q)$ thermodynamic inconsistency (Remark 1) / an $\mathcal O(\Delta x^2)$ one. Worth noting for implementation.

### 3.4 Second-order numerical dissipation (MUSCL-type reconstruction for the explicit fluxes)

Central FD operators are already 2nd order. Only the dissipation in (26) must be made $\mathcal O(\Delta x^2)$:
$$
w^n_i(x):=m^n_i+(x-x_i)\alpha_{x,i},\tag{33}
$$
minmod:
$$
\alpha_{x,i}=\begin{cases}0&\text{if }\alpha_{x,\ell}\alpha_{x,r}\le0\\\alpha_{x,\ell}&\text{if }|\alpha_{x,\ell}|\le|\alpha_{x,r}|\\\alpha_{x,r}&\text{if }|\alpha_{x,\ell}|>|\alpha_{x,r}|\end{cases},\qquad\alpha_{x,\ell}=\frac{m^n_i-m^n_{i-1}}{\Delta x},\quad\alpha_{x,r}=\frac{m^n_{i+1}-m^n_i}{\Delta x}.\tag{34}
$$
Dissipation jumps use extrapolated boundary values:
$$
D^m_{i+1/2}=\tfrac12a^n_\Omega\big(w^n_{i+1}(x_{i+1/2})-w^n_i(x_{i+1/2})\big),\qquad D^m_{i-1/2}=\tfrac12a^n_\Omega\big(w^n_i(x_{i-1/2})-w^n_{i-1}(x_{i-1/2})\big),\tag{35}
$$
**central part of the flux untouched** (still $\frac12(f(m_{i+1})+f(m_i))$ with cell averages). Same applies to any artificial viscosity term.

#### Numerical dissipation in the volume-fraction equation (optional, for strong $\phi$ gradients)
Explicit artificial viscosity added on the RHS of (23):
$$
\phi^{n+1}_{s,i}+\frac{\Delta t}{2\Delta x}\Big(u^n_{I,i+1/2}(\phi^{n+1}_{s,i+1}-\phi^{n+1}_{s,i})-u^n_{I,i-1/2}(\phi^{n+1}_{s,i-1}-\phi^{n+1}_{s,i})\Big)=\phi^n_{s,i}+\frac{\Delta t}{\Delta x^2}\Big(\nu^+(\phi^n_{s,i+1}-\phi^n_{s,i})-\nu^-(\phi^n_{s,i}-\phi^n_{s,i-1})\Big),\tag{36}
$$
$$
\nu^\pm(\sigma)=\frac{\sigma^n_i+\sigma^n_{i\pm1}}{2},\qquad\sigma^n_i=\tfrac12\Delta x\,a^n_\Omega.\tag{37}
$$
(Since $a_\Omega$ is global, $\nu^\pm=\tfrac12\Delta x a_\Omega$; equivalent to the first-order Rusanov dissipation on $\phi$; the MUSCL version via (35) presumably applies too **[uncertain]**.)

**Crucial consistency requirement:** this dissipation must also be added in the pressure system (29) and in the conservative total-energy update (32b), with rescaled coefficients
$$
\sigma^n_{i,\star}=\sigma^n_i\,\frac{\partial(\rho_\alpha\phi_\alpha e_\alpha)^n_i}{\partial\phi^n_{\alpha,i}}.\tag{38}
$$
(For stiffened gas at fixed $p$: $\partial(\rho\phi e)/\partial\phi=(p_\alpha+\gamma_\alpha\pi_\alpha)/(\gamma_\alpha-1)$ — **[note]**, paper's (71) shows $p_0/(\gamma-1)$ for ideal gas.) I.e. the diffusion of $\phi$ must be mirrored by a diffusion of internal energy $\propto\nabla^2\phi$, otherwise pressure oscillations appear at contacts (Abgrall's pressure-oscillation issue). Details in Sec. 3.7.

#### Numerical dissipation in the energy equation (optional, high Mach)
Kinetic-energy convective flux augmented by a pressure-jump viscosity:
$$
f^{\rho_\alpha\phi_\alpha k_\alpha}_{i+1/2}=\tfrac12\Big(((\rho_\alpha k_\alpha u_\alpha)^n\phi^{n+1}_\alpha)_{i+1}+((\rho_\alpha k_\alpha u_\alpha)^n\phi^{n+1}_\alpha)_i\Big)-\tfrac12a^n_\Omega\Big((\rho_\alpha\phi_\alpha k_\alpha)^n_{i+1}-(\rho_\alpha\phi_\alpha k_\alpha)^n_i\Big)-\frac{1}{\Delta x}\nu^+(\sigma_{\star\star})\Big(p^n_{\alpha,i+1}-p^n_{\alpha,i}\Big),\tag{39}
$$
$$
\sigma^n_{i,\star\star}=\sigma^n_i\,\frac{\partial(\rho_\alpha\phi_\alpha e_\alpha)^n_i}{\partial p^n_{\alpha,i}}.\tag{40}
$$
(For stiffened gas $\partial(\rho\phi e)/\partial p=\phi/(\gamma-1)$ — **[note]**.) Only active for very high Mach flows (it is not scaled with $\epsilon$-independent velocity only... the paper notes in Sec. 3.6 that it is only activated at very high Mach).

### 3.5 Relaxation solver (sources)

Homogeneous-step result labelled $\mathbf Q^H$. Strang splitting [63] **[note: as written it is a first-order Lie step per stage; the text says "Strang"]**; the ODE step:
$$
\frac{\mathbf Q^{n+1}_i-\mathbf Q^H_i}{\Delta t}=\mathbf S(\mathbf Q^{n+1}_i),\qquad i=1,\dots,N_x,\tag{41}
$$
discretised semi-implicitly (linear, no nonlinearity), in $x$:
$$
\mathbf S(\mathbf Q^{n+1}_i)=\begin{bmatrix}0\\ \delta(u^{n+1}_{s,i}-u^{n+1}_{f,i})\\ \mu\,p^H_{I,i}(p^H_{s,i}-p^H_{f,i})+\delta\,u^H_{I,i}(u^{n+1}_{s,i}-u^{n+1}_{f,i})\\ 0\\ -\delta(u^{n+1}_{s,i}-u^{n+1}_{f,i})\\ -\mu\,p^H_{I,i}(p^H_{s,i}-p^H_{f,i})-\delta\,u^H_{I,i}(u^{n+1}_{s,i}-u^{n+1}_{f,i})\\ \mu(p^H_{s,i}-p^H_{f,i})\end{bmatrix}.\tag{42}
$$
Key points:
- **Pressures are frozen at the homogeneous solution $p^H$** (taken from the wave equation (29)) → pressure relaxation is effectively **explicit** in this paper: $\phi_s^{n+1}=\phi_s^H+\Delta t\,\mu(p^H_s-p^H_f)$, energy exchange $\pm\Delta t\,\mu p^H_I(p^H_s-p^H_f)$. **[note]** This is *not* an instantaneous/stiff pressure relaxation; for $\mu\Delta t\gg1$ it would overshoot. RP7 uses $\mu=10^2$, $\delta=10^3$ (dimensionless) and works, but with stiff $\mu\to\infty$ you need a proper pressure-relaxation procedure (e.g. Saurel–Petitpas–Berry / Pelanti–Shyue) — the paper does not provide one.
- Velocities implicit: cell-local 2×2 linear system (per direction) for $(u^{n+1}_s,u^{n+1}_f)$, inspired by Re–Abgrall [46]: with $m_\alpha=(\rho_\alpha\phi_\alpha)^H$ (unchanged by sources),
  $m_f u_f^{n+1}=m_fu_f^H+\Delta t\,\delta(u_s^{n+1}-u_f^{n+1})$, $m_su_s^{n+1}=m_su_s^H-\Delta t\,\delta(u_s^{n+1}-u_f^{n+1})$ **[note: closed form: mixture momentum conserved, $\Delta u^{n+1}=\Delta u^H/(1+\Delta t\,\delta(1/m_s+1/m_f))$]**.
- Then energy sources computed with $u^{n+1}$ (drag work uses $u_I^H$).
- 2nd order in space automatic (cell-centred = midpoint rule); time accuracy via the IMEX scheme (21) **[uncertain: how the relaxation step is embedded in the IMEX stages is not detailed — the text suggests the sources are part of the implicit operator within (21)]**.

### Summary of the linearly implicit scheme (per stage / per step)
1. Solve implicitly (23) for $\phi^{n+1}_s$ (GMRES).
2. Explicit convective fluxes with FV (25)–(26): obtain $\rho^\star_\alpha$, $q^\star_\alpha$, $E^\star_\alpha$; note $\rho^\star_\alpha=(\rho_\alpha\phi_\alpha)^{n+1}$.
3. Solve linearly implicit wave equation (29) for $(\phi_\alpha p_\alpha)^{n+1}$ (GMRES), each phase.
4. Update momentum (32a) → $(\rho_\alpha\phi_\alpha u_\alpha)^{n+1}$.
5. Update energy (32b) → $(\rho_\alpha\phi_\alpha E_\alpha)^{n+1}$.
6. Relax with sources, semi-implicit discretisation (42).

**[note]** Order matters: $\phi$ first (it enters every flux as $\phi^{n+1}$), then explicit fluxes (which need $\phi^{n+1}$), then pressure, then momentum (needs $p^{n+1}$), then energy (needs momentum$^{n+1}$, and $h^n$ via (18) needs $\rho^{n+1}$ from step 2).

### 3.6 Asymptotic-preserving property

AP: for $\epsilon\to0$ with $\Delta t$ independent of $\epsilon$, the scheme gives a consistent, stable discretisation of the limit model.

**Continuous limit.** Bounded $\Omega$, impermeable walls
$$
\mathbf u_\alpha\cdot\mathbf n=0\ \text{on }\partial\Omega.\tag{43}
$$
Expansion $m=m_{(0)}+\epsilon m_{(1)}+\epsilon^2m_{(2)}+\dots+\mathcal O(\epsilon^k)$ (44), well-prepared data. Volume fraction eq. with finite $\mu$:
$$
\mathcal O(\epsilon^{-1}):\ \mu(p_{s(0)}-p_{f(0)})=0\Rightarrow p_{s(0)}=p_{f(0)}=p_{I(0)}=p_{(0)},\tag{45}
$$
$$
\mathcal O(\epsilon^0):\ \partial_t\phi_{\alpha(0)}+\mathbf u_{I(0)}\cdot\nabla\phi_{\alpha(0)}=\pm\mu(p_{s(1)}-p_{f(1)}).\tag{46}
$$
Continuity: $\partial_t(\rho_{\alpha(0)}\phi_{\alpha(0)})+\nabla\cdot(\rho_{\alpha(0)}\phi_{\alpha(0)}\mathbf u_{\alpha(0)})=0$ (47). Momentum:
$$
\mathcal O(\epsilon^{-1}):\ \nabla(p_{\alpha(0)}\phi_{\alpha(0)})-p_{I(0)}\nabla\phi_{\alpha(0)}=0,\tag{48a}
$$
$$
\mathcal O(\epsilon^0):\ \partial_t(\rho_{\alpha(0)}\phi_{\alpha(0)}\mathbf u_{\alpha(0)})+\nabla\cdot(\rho_{\alpha(0)}\phi_{\alpha(0)}\mathbf u_{\alpha(0)}\otimes\mathbf u_{\alpha(0)})+\nabla(p_{\alpha(1)}\phi_{\alpha(1)})-p_{I(1)}\nabla\phi_{\alpha(1)}=\mp\delta(\mathbf u_{s(0)}-\mathbf u_{f(0)}).\tag{48b}
$$
Summing (48a) over phases with (45) and $\sum\phi=1$: $\nabla p_{(0)}=0$, $p_{(0)}=p_{(0)}(t)$ (49). Energy:
$$
\mathcal O(\epsilon^0):\ \partial_t(\rho_{\alpha(0)}\phi_{\alpha(0)}e_{\alpha(0)})+\nabla\cdot(\rho_{\alpha(0)}\phi_{\alpha(0)}h_{\alpha(0)}\mathbf u_{\alpha(0)})-p_{I(0)}\mathbf u_{I(0)}\cdot\nabla\phi_{\alpha(0)}=\mp\mu p_{I(1)}(p_{s(1)}-p_{f(1)}),\tag{50a}
$$
$$
\mathcal O(\epsilon^1):\ \partial_t(\rho_{\alpha(0)}\phi_{\alpha(0)}k_{\alpha(0)})+\nabla\cdot(\rho_{\alpha(0)}\phi_{\alpha(0)}k_{\alpha(0)}\mathbf u_{\alpha(0)})=\mp\delta\,\mathbf u_{I(0)}(\mathbf u_{s(0)}-\mathbf u_{f(0)}).\tag{50b}
$$
Ideal gas ($\pi=0$), integrate (50a) over $\Omega$:
$$
\frac{1}{\gamma_\alpha-1}\int_\Omega\frac{d(\phi_{\alpha(0)}p_{(0)})}{dt}dx+\frac{\gamma_\alpha}{\gamma_\alpha-1}p_{(0)}\int_{\partial\Omega}\phi_{\alpha(0)}\mathbf u_{\alpha(0)}\cdot\mathbf n\,dS-p_{(0)}\int_\Omega\mathbf u_{I(0)}\cdot\nabla\phi_{\alpha(0)}dx=\mp\mu\int_\Omega p_{I(1)}(p_{s(1)}-p_{f(1)})dx.\tag{51}
$$
Sum over phases + walls + (46):
$$
\frac{\gamma_f+\gamma_s-2}{(\gamma_s-1)(\gamma_f-1)}|\Omega|\frac{dp_{(0)}}{dt}-p_{(0)}\int_\Omega\mathbf u_{I(0)}\nabla(\phi_{s(0)}+\phi_{f(0)})dx=-\int_\Omega p_{I(1)}\mathbf u_{I(0)}\nabla(\phi_{s(0)}+\phi_{f(0)})dx,\tag{52}
$$
$$
\frac{dp_{(0)}}{dt}=0.\tag{53}
$$
**[uncertain]** (52) coefficient as printed; it looks like it should be a weighted sum of $\phi_\alpha/(\gamma_\alpha-1)$; I transcribe as printed.
Then (50a) for the mixture gives
$$
p_{(0)}\nabla\cdot(\phi_{s(0)}\mathbf u_{s(0)}+\phi_{f(0)}\mathbf u_{f(0)})=0,\tag{54}
$$
and with instantaneous relaxation $\mathbf u_{s(0)}=\mathbf u_{f(0)}=\mathbf u_{I(0)}=\mathbf u_{(0)}$:
$$
\nabla\cdot\mathbf u_{(0)}=0.\tag{55}
$$
Limit = incompressible Euler (for a single phase).

**Discrete AP.** Insert expansion into (12a)–(12d) plus sources (42):
$$
\mathcal O(\epsilon^{-1}):\ \mu(p^H_{s(0)}-p^H_{f(0)})=0\Rightarrow p^H_{s(0)}=p^H_{f(0)}=p^H_{I(0)}=p^H_{(0)}.\tag{56}
$$
From (12b): $\nabla(p^H_{s(0)}\phi^{n+1}_{s(0)}+p^H_{f(0)}\phi^{n+1}_{f(0)})=p^n_{I(0)}\nabla(\phi_{s(0)}+\phi_{f(0)})^{n+1}$, $\nabla(p^H_{(0)}(\phi_{s(0)}^{n+1}+\phi_{f(0)}^{n+1}))=p^n_{I(0)}\cdot0\Rightarrow\nabla p^H_{(0)}=0$ (57). Leading order discrete system:
$$
\begin{aligned}
(\rho_{\alpha(0)}\phi_{\alpha(0)})^{n+1}&=(\rho_{\alpha(0)}\phi_{\alpha(0)})^n-\Delta t\nabla\cdot((\rho_{\alpha(0)}\mathbf u_{\alpha(0)})^n\phi^{n+1}_{\alpha(0)}),&(58a)\\
(\rho_{\alpha(0)}\phi_{\alpha(0)}\mathbf u_{\alpha(0)})^{n+1}&=(\dots)^n-\Delta t\nabla\cdot((\rho\mathbf u\otimes\mathbf u)^n_{\alpha(0)}\phi^{n+1}_{\alpha(0)})-\Delta t\nabla(\phi^{n+1}_{\alpha(1)}p^H_{\alpha(1)})+\Delta t\,p^n_{I(1)}\nabla\phi^{n+1}_{\alpha(1)}\mp\Delta t\,\delta(\mathbf u^{n+1}_{s(0)}-\mathbf u^{n+1}_{f(0)}),&(58b)\\
(\rho_{\alpha(0)}\phi_{\alpha(0)}e_{\alpha(0)})^{n+1}&=(\dots)^n-\Delta t\nabla\cdot(h^n_{\alpha(0)}(\rho_{\alpha(0)}\phi_{\alpha(0)}\mathbf u_{\alpha(0)})^{n+1})+\Delta t\,p^n_{I(0)}\mathbf u^n_{I(0)}\nabla\phi^{n+1}_{\alpha(0)}\mp\Delta t\,\mu p^H_{I(1)}(p^H_{s(1)}-p^H_{f(1)}),&(58c)\\
\phi^{n+1}_{\alpha(0)}&=\phi^n_{\alpha(0)}-\Delta t\,\mathbf u^n_{I(0)}\nabla\phi^{n+1}_{\alpha(0)}\pm\Delta t\,\mu(p^H_{s(1)}-p^H_{f(1)}),&(58d)
\end{aligned}
$$
kinetic terms at $\mathcal O(\epsilon)$:
$(\rho\phi k)^{n+1}_{\alpha(0)}=(\rho\phi k)^n_{\alpha(0)}-\Delta t\nabla\cdot((\rho\phi k\mathbf u)^n_{\alpha(0)}\phi^{n+1}_{\alpha(0)})\mp\delta\mathbf u^{n+1}_{I(0)}(\mathbf u^{n+1}_{s(0)}-\mathbf u^{n+1}_{f(0)})$ (59) [as printed].
Integrating (58c) with ideal gas, enthalpy (18):
$$
p^{n+1}_{(0)}\int_\Omega\phi^{n+1}_{\alpha(0)}dx=p^n_{(0)}\int_\Omega\phi^n_{\alpha(0)}dx-\Delta t\,p^n_{(0)}\int_{\partial\Omega}\phi^{n+1}_{\alpha(0)}\mathbf u^{n+1}_{\alpha(0)}\cdot\mathbf n\,dS+\Delta t\,p^n_{(0)}\int_\Omega\mathbf u^n_{I(0)}\nabla\phi^{n+1}_{\alpha(0)}dx\mp\Delta t\,\mu\int_\Omega p^H_{I(1)}(p^H_{s(1)}-p^H_{f(1)})dx\tag{60}
$$
[the $\gamma$ factors appear to have been absorbed — as printed]. Summing phases + walls + (56): $p^{n+1}_{(0)}=p^n_{(0)}$ (61). Then
$$
p_{(0)}\phi^{n+1}_{\alpha(0)}=p_{(0)}\phi^n_{\alpha(0)}-\Delta t\,p_{(0)}\nabla\cdot(\phi_{\alpha(0)}\mathbf u_{\alpha(0)})^{n+1}+\Delta t\,p^n_{(0)}\mathbf u^n_{I(0)}\nabla\phi^{n+1}_{\alpha(0)}\mp\Delta t\,\mu p^H_{I(1)}(p^H_{s(1)}-p^H_{f(1)}),\tag{62}
$$
which for the mixture gives $\nabla\cdot(\phi_{s(0)}\mathbf u_{s(0)}+\phi_{f(0)}\mathbf u_{f(0)})^{n+1}=0$, the discrete analogue of (54). → AP.

**Remark 2.** If the phases are treated independently (not as a mixture), relaxation terms vanish and this analysis does not apply; separate Mach numbers per phase needed ([51]). Single-phase case: recovers the AP scheme of Boscheri–Pareschi [43].

Numerical dissipation does not dominate for small $\epsilon$ (cf. Dellacherie [22]): the viscosity in Sec. 3.4 is proportional only to the material velocity, so independent of $\epsilon$ by construction. Extra pressure-jump dissipation (39) only activated for very high Mach.

### 3.7 Well-balanced property (moving contact discontinuity)

1D, $\Omega=[x_L,x_R]$, discontinuity at $x_D$:
$$
\rho_{\alpha,i}(t_0)=\begin{cases}\rho_{\alpha_L}&x\le x_D\\\rho_{\alpha_R}&x>x_D\end{cases},\quad
\phi_{\alpha,i}(t_0)=\begin{cases}\phi_{\alpha_L}&x\le x_D\\\phi_{\alpha_R}&x>x_D\end{cases},\quad u_{\alpha,i}(t_0)=u_0,\quad p_{\alpha,i}(t_0)=p_0.\tag{63}
$$
Sources vanish. Convective subsystem (28) gives
$$
\rho^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha)^{n+1}_i,\tag{64a}
$$
$$
q^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha u_\alpha)^n_i-\frac{\Delta t}{\Delta x}(f^{\rho\phi u}_{i+1/2}-f^{\rho\phi u}_{i-1/2})^n=u_0\Big((\rho_\alpha\phi_\alpha)^n_i-\frac{\Delta t}{\Delta x}(f^{\rho\phi}_{i+1/2}-f^{\rho\phi}_{i-1/2})^n\Big)=u_0(\rho_\alpha\phi_\alpha)^{n+1}_i,\tag{64b}
$$
$$
E^\star_{\alpha,i}=(\rho_\alpha\phi_\alpha e_\alpha)^n_i+\frac{u_0^2}{2}(\rho_\alpha\phi_\alpha)^n_i-\frac{u_0^2}{2}\frac{\Delta t}{\Delta x}(f^{\rho\phi}_{i+1/2}-f^{\rho\phi}_{i-1/2})^n=(\rho_\alpha\phi_\alpha e_\alpha)^n_i+\frac{u_0^2}{2}(\rho_\alpha\phi_\alpha)^{n+1}_i,\tag{64c}
$$
$$
\phi^\star_{\alpha,i}=\phi^n_{\alpha,i}-\frac{\Delta t}{2\Delta x}u_0(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i-1}):=\phi^{n+1}_{\alpha,i}.\tag{64d}
$$
Wave equation (17) with ideal gas, discrete (65): LHS
$p_0\frac{\phi^{n+1}_{\alpha,i}}{\gamma-1}-\frac{\Delta t}{4\Delta x}u_0p_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})-\frac{\Delta t^2}{\Delta x^2}p_0\Big((\tfrac34h_{i-1}+\tfrac14h_{i+1})\phi^{n+1}_{i-1}-(h_{i-1}+h_{i+1})\phi^{n+1}_i+(\tfrac14h_{i-1}+\tfrac34h_{i+1})\phi^{n+1}_{i+1}\Big)$;
RHS
$p_0\frac{\phi^n_{\alpha,i}}{\gamma-1}+\frac{u_0^2}{2}(\rho\phi)^{n+1}_i-\frac{u_0^2}{2}(\rho\phi)^{n+1}_i-\frac{\Delta t}{4\Delta x}p_0u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})-\frac{\Delta t}{2\Delta x}u_0\big((\rho\phi)^{n+1}_{i+1}h_{i+1}-(\rho\phi)^{n+1}_{i-1}h_{i-1}\big)+\frac{\Delta t}{2\Delta x}p_0u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})-\frac{\Delta t^2}{\Delta x^2}p_0\big(\text{same compact operator on }\phi^{n+1}\big)$.

The two $\Delta t^2$ terms cancel **because $\nabla\cdot(h p_I\nabla\phi)$ uses the same compact operator (31) as $\nabla\cdot(h\nabla(\phi p))$** — this is *the* design requirement. Using (18), $h^n_\alpha=\frac{\rho^n_\alpha}{\rho^{n+1}_\alpha}(e^n_\alpha+p^n_\alpha/\rho^n_\alpha)$ so $(\rho\phi)^{n+1}_{i\pm1}h^n_{i\pm1}=(\rho_\alpha e_\alpha)^n_{i\pm1}\phi^{n+1}_{i\pm1}+\dots$; printed result: $-\frac{\Delta t}{2\Delta x}u_0\big((\rho_\alpha e_\alpha)^n_{i+1}\phi^{n+1}_{\alpha,i+1}-(\rho_\alpha e_\alpha)^n_{i-1}\phi^{n+1}_{\alpha,i-1}\big)$ **[uncertain: the $p/\rho$ part is lumped into $\rho e$ notation; for ideal gas $\rho h=\gamma p/(\gamma-1)$]**. Left with
$$
p_0\frac{\phi^{n+1}_{\alpha,i}}{\gamma-1}=p_0\frac{\phi^n_{\alpha,i}}{\gamma-1}-\frac{\Delta t}{2\Delta x}u_0\big((\rho_\alpha e_\alpha)^n_{i+1}\phi^{n+1}_{\alpha,i+1}-(\rho_\alpha e_\alpha)^n_{i-1}\phi^{n+1}_{\alpha,i-1}\big)-\frac{\Delta t}{2\Delta x}u_0p_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}p_0u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1}),\tag{66}
$$
last two cancel; with $\rho e=p_0/(\gamma-1)$:
$$
p_0\frac{\phi^{n+1}_{\alpha,i}}{\gamma-1}=p_0\frac{1}{\gamma-1}\Big(\phi^n_{\alpha,i}-\frac{\Delta t}{2\Delta x}u_0(\phi^{n+1}_{\alpha,i+1}-\phi^{n+1}_{\alpha,i-1})\Big)=p_0\frac{\phi^{n+1}_{\alpha,i}}{\gamma-1},\tag{67}
$$
by (64d) ⇒ $p_0$ solves the pressure system. Momentum and energy:
$$
(\rho_\alpha\phi_\alpha u_\alpha)^{n+1}_i=q^\star_{\alpha,i}-\frac{\Delta t}{2\Delta x}p_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}p_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})=u_0(\rho_\alpha\phi_\alpha)^{n+1}_i,\tag{68a}
$$
$$
(\rho_\alpha\phi_\alpha E_\alpha)^{n+1}_i=E^\star_{\alpha,i}-\frac{\Delta t}{2\Delta x}u_0\Big(\frac{p_0}{\gamma-1}+p_0\Big)(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}p_0u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})=\frac{p_0}{\gamma-1}\phi^{n+1}_{\alpha,i}+\frac{u_0^2}{2}(\rho_\alpha\phi_\alpha)^{n+1}_i.\tag{68b}
$$

**[note] Implementation checklist for well-balancing:** (i) $\phi$ eq. and all $p_I\nabla\phi$, $p_Iu_I\nabla\phi$ terms use the *same* central fluctuation form; (ii) mass, momentum and kinetic-energy convective fluxes all carry the same $\phi^{n+1}$ and the same Rusanov dissipation coefficient (so $q^\star=u_0\rho^\star$, $E^\star$ kinetic part $=\tfrac{u_0^2}{2}\rho^\star$); (iii) energy dissipation acts on $\rho\phi k$ only; (iv) enthalpy via (18); (v) the two $\Delta t^2$ second-derivative terms use the identical compact operator (31).

**With numerical dissipation.** $a_\Omega=u_0$, $\nu^\pm=\tfrac12\Delta x u_0$:
$$
\phi^{n+1}_{s,i}=\phi^n_{s,i}-\frac{\Delta t}{2\Delta x}u_0(\phi^{n+1}_{s,i+1}-\phi^{n+1}_{s,i-1})+\frac{\Delta t}{2\Delta x}u_0(\phi^n_{s,i+1}-2\phi^n_{s,i}+\phi^n_{s,i-1}).\tag{69}
$$
Must be added to RHS of (65), rescaled by
$$
\frac{\partial(\rho_\alpha\phi_\alpha e_\alpha)^n_i}{\partial\phi^n_{\alpha,i}}=\frac{\partial}{\partial\phi^n_{\alpha,i}}\Big(\frac{(\phi_\alpha p_\alpha)^n_i}{\gamma_\alpha-1}\Big)=\frac{p_0}{\gamma_\alpha-1},\tag{71}
$$
giving (70): $p_0\phi^{n+1}/(\gamma-1)=\frac{p_0}{\gamma-1}(\phi^n-\frac{\Delta t}{2\Delta x}u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}u_0(\phi^n_{s,i+1}-2\phi^n_{s,i}+\phi^n_{s,i-1}))$. Same added term to total energy (68b) → (72):
$(\rho\phi E)^{n+1}_i=E^\star_i-\frac{\Delta t}{2\Delta x}u_0(\frac{p_0}{\gamma-1}+p_0)(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}p_0u_0(\phi^{n+1}_{i+1}-\phi^{n+1}_{i-1})+\frac{\Delta t}{2\Delta x}\frac{p_0}{\gamma-1}u_0(\phi^n_{s,i+1}-2\phi^n_{s,i}+\phi^n_{s,i-1})=\frac{p_0}{\gamma-1}\phi^{n+1}_{\alpha,i}+\frac{u_0^2}{2}(\rho\phi)^{n+1}_i$.
**[note]** In (69)–(72) the dissipation for phase $\alpha$ is written with $\phi_s$; for the fluid phase it must be applied with $\phi_f$ (same magnitude since $\phi_f=1-\phi_s$, Laplacian flips sign).

---

## 4. Numerical results (Section 4)

Scheme labelled **SIFV-BN**. Mach numbers in $[10^{-3},10^1]$. $\Delta t$ from (11) with **CFL = 0.9**. If initial velocity is zero, the **first** time step uses the acoustic CFL (6). Dimensionless form used unless stated.

### 4.1 Convergence: smooth (isentropic) vortex
- Modified Hu–Shu vortex [66], setup as in Dumbser et al. [67] (details deferred there).
- $\Omega=[-5,5]^2$, periodic. Ideal gases $\gamma_f=1.35$, $\gamma_s=1.4$. Mean advection $\mathbf v_c=(2,2)$; $t_f=5$ (one full period back to origin).
- Background pressure $p_0\in[1,10^4]$ to vary Mach. Fig. 1: solid pressure for $p_0=1,10^2,10^4$ on $\varrho=1/16$ — no spurious oscillations, same accuracy.
- $L_2$ errors vs. exact solution, Table 1:

| $N_x$ | $\varrho$ | $L_2(\rho_s\phi_s)$ | ord | $L_2(\rho_s\phi_su_s)$ | ord | $L_2(\phi_s)$ | ord |
|---|---|---|---|---|---|---|---|
| **$p_0=1$** |||||||
| 20 | 0.5 | 4.2250E−01 | – | 1.0690E+00 | – | 1.6308E−01 | – |
| 40 | 0.25 | 1.6787E−01 | 1.33 | 3.4438E−01 | 1.63 | 6.3713E−02 | 1.36 |
| 80 | 0.125 | 4.6463E−02 | 1.85 | 8.7884E−02 | 1.97 | 1.6970E−02 | 1.91 |
| 160 | 0.0625 | 1.1957E−02 | 1.96 | 2.5092E−02 | 1.81 | 4.2119E−03 | 2.01 |
| **$p_0=10$** |||||||
| 20 | 0.5 | 3.5478E−02 | – | 1.0726E+00 | – | 1.6402E−01 | – |
| 40 | 0.25 | 1.3797E−02 | 1.36 | 3.5107E−01 | 1.61 | 6.3683E−02 | 1.36 |
| 80 | 0.125 | 3.2493E−03 | 2.09 | 8.8656E−02 | 1.99 | 1.6530E−02 | 1.95 |
| 160 | 0.0625 | 7.6859E−04 | 2.08 | 2.5126E−02 | 1.82 | 4.0399E−03 | 2.03 |
| **$p_0=10^2$** |||||||
| 20 | 0.5 | 4.9215E−03 | – | 1.0725E+00 | – | 1.6412E−01 | – |
| 40 | 0.25 | 1.6694E−03 | 1.56 | 3.5247E−01 | 1.61 | 6.3743E−02 | 1.36 |
| 80 | 0.125 | 4.2384E−04 | 1.98 | 8.9128E−02 | 1.98 | 1.6527E−02 | 1.95 |
| 160 | 0.0625 | 1.0384E−04 | 2.03 | 2.5228E−02 | 1.82 | 4.0357E−03 | 2.03 |
| **$p_0=10^3$** |||||||
| 20 | 0.5 | 2.8493E−03 | – | 1.0675E+00 | – | 1.6386E−01 | – |
| 40 | 0.25 | 1.2504E−03 | 1.19 | 3.4870E−01 | 1.61 | 6.3264E−02 | 1.37 |
| 80 | 0.125 | 3.3387E−04 | 1.90 | 8.7773E−02 | 1.99 | 1.6292E−02 | 1.96 |
| 160 | 0.0625 | 8.2761E−05 | 2.01 | 2.4899E−02 | 1.82 | 3.9670E−03 | 2.04 |
| **$p_0=10^4$** |||||||
| 20 | 0.5 | 9.4160E−03 | – | 1.0305E+00 | – | 1.6343E−01 | – |
| 40 | 0.25 | 4.2722E−03 | 1.14 | 3.2406E−01 | 1.67 | 6.2452E−02 | 1.39 |
| 80 | 0.125 | 9.1133E−04 | 2.23 | 7.9132E−02 | 2.03 | 1.5855E−02 | 1.98 |
| 160 | 0.0625 | 1.4772E−04 | 2.63 | 2.2924E−02 | 1.79 | 3.8419E−03 | 2.05 |

($N_x$ here is cells per direction; grid $N_x\times N_x$. Momentum order saturates ≈1.8.)

### 4.2 1D Riemann problems
$\Omega=[0,1]$, $N_x=1000$, Dirichlet BCs, $x_D=0.5$. Reference: explicit 2nd-order TVD FV on 10 000 cells (RP7: exact Euler Riemann solver with different $\gamma$ either side of contact). Table 2:

| RP | EOS f | EOS s | $(\delta,\mu)$ | side | $\rho_f$ | $u_f$ | $p_f$ | $\rho_s$ | $u_s$ | $p_s$ | $\phi_s$ | $t_f$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RP0 | 1.4, 0 | 1.4, 0 | 0, 0 | L | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 0.8 | 0.2 |
| | | | | R | 0.5 | 1.0 | 1.0 | 0.5 | 1.0 | 1.0 | 0.4 | |
| RP1 [69] | 1.4, 0 | 1.4, 0 | 0, 0 | L | 0.5 | 0.0 | 1.0 | 1.0 | 0.0 | 1.0 | 0.4 | 0.1 |
| | | | | R | 1.5 | 0.0 | 2.0 | 2.0 | 0.0 | 2.0 | 0.8 | |
| RP2 [69] | 1.4, 0 | 3.0, 100 | 0, 0 | L | 1.5 | 0.0 | 2.0 | 800.0 | 0.0 | 500.0 | 0.4 | 0.1 |
| | | | | R | 1.0 | 0.0 | 1.0 | 1000.0 | 0.0 | 600.0 | 0.3 | |
| RP3 [69] | 1.4, 0 | 1.4, 0 | 0, 0 | L | 1.0 | 0.0 | 1.0 | 1.0 | 0.9 | 2.5 | 0.9 | 0.1 |
| | | | | R | 1.2 | 1.0 | 2.0 | 1.0 | 0.0 | 1.0 | 0.2 | |
| RP4 [70] | 1.35, 0 | 3.0, 3400 | 0, 0 | L | 2.0 | 0.0 | 3.0 | 1900.0 | 0.0 | 10.0 | 0.2 | 0.15 |
| | | | | R | 1.0 | 0.0 | 1.0 | 1950.0 | 0.0 | 1000.0 | 0.9 | |
| RP5 [70] | 1.4, 0 | 1.4, 0 | 0, 0 | L | 0.2 | 0.0 | 0.3 | 1.0 | 0.0 | 1.0 | 0.8 | 0.2 |
| | | | | R | 1.0 | 0.0 | 1.0 | 1.0 | 0.0 | 1.0 | 0.3 | |
| RP6 [68] | 1.4, 0 | 1.4, 0 | 0, 0 | L | 0.5806 | 1.5833 | 1.375 | 0.2068 | 1.4166 | 0.0416 | 0.1 | 0.1 |
| | | | | R | 0.4890 | −0.70138 | 0.986 | 2.2263 | 0.9366 | 6.0 | 0.2 | |
| RP7 [54] | 1.67, 0 | 1.4, 0 | $10^3$, $10^2$ | L | 1.0 | 0.0 | 1.0 | 1.0 | 0.0 | 1.0 | 0.99 | 0.2 |
| | | | | R | 0.125 | 0.0 | 0.1 | 0.125 | 0.0 | 0.1 | 0.01 | |

(EOS columns: $\gamma,\pi$.) Results:
- **RP0** (moving contact, Fig. 2): constant $p=1$ and $u=1$ of both phases preserved to plotting precision ($\pm10^{-6}$ axis range shows flat lines); density and $\phi_s$ jumps captured (smeared over ~5–6 cells). Numerical confirmation of well-balancing.
- **RP1–RP4** (Figs. 3–6): shock tubes; rarefactions and shocks captured well, slightly more difficulty for RP3 and RP4 (RP4 shows mild smearing / small deviation in $\rho_f$ around the solid contact).
- **RP5, RP6** (Figs. 7, 8): small oscillations, but agreement very good (RP6 shows small wiggles near $x\approx0.4$–0.45 in $\rho_s,\rho_f,u_f$).
- **RP7**: classical Sod with stiff relaxation (pressure + velocity), compared with exact Euler solution with different $\gamma$ either side of the contact (Fig. 9) — satisfying.
- Authors' caveat: method is not designed specifically for sharp discontinuities but handles high-Mach regime.

### Explosion problem with/without drag relaxation (2D)
- $\Omega=[-1,1]^2$, $N_x\times N_y=1000\times1000$. Inside circle $R=0.4$: L state of RP2, outside: R state of RP2. $t_f=0.15$. Cases $\delta=10^5$ and $\delta=0$ ($\mu$ presumably 0 **[uncertain]**).
- Reference: explicit 2nd-order TVD FV on 10 000 cells (1D radial — presumably a 1D radial reference with geometric source **[uncertain]**). Fig. 10: excellent agreement; for $\delta=10^5$ phase velocities relax to the mixture velocity. Cylindrical symmetry preserved.
- Efficiency vs fully explicit: Table 5 for $\delta=10^5$, speed-up factor **43** (due to the stiff drag handled implicitly).

### 4.3 2D Riemann problems (start)
MPI-parallel code. Two 2D Riemann problems from Dumbser–Boscheri [54] (inspired by Zhang–Zheng [71], Kurganov–Tadmor [72]). Four constant states in the quadrants of $\Omega=[-0.5,0.5]^2$, reflective walls everywhere. Parameters in Table 3 (next pages).
