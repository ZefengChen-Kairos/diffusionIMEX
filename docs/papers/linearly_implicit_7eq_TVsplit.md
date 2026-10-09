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
