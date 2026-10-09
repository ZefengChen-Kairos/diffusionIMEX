# Reading note: Saurel, Petitpas & Berry (2009) — the "6+1" (6-equation + total energy) relaxation method

**Title:** *Simple and efficient relaxation methods for interfaces separating compressible fluids, cavitating flows and shocks in multiphase mixtures*
**Authors:** Richard Saurel (Polytech'Marseille / IUSTI / Univ. Institute of France), Fabien Petitpas (Polytech'Marseille / IUSTI), Ray A. Berry (Idaho National Laboratory)
**Venue:** *Journal of Computational Physics* **228** (2009) 1678–1712. doi:10.1016/j.jcp.2008.11.002
(Received 16 Apr 2008, revised 13 Oct 2008, accepted 3 Nov 2008, online 13 Nov 2008.)

File in repo: `6Plus1TotalEnergy.pdf` (35 pages). Equation numbers below are the paper's (Roman section numeral + index, e.g. (II.6), (III.4)).

> Note on scope vs. our project. This paper is the origin of the **single-velocity, two-pressure 6-equation model + redundant mixture total energy equation** ("6+1"). It treats **only stiff mechanical (pressure) relaxation** in the main text; heat/mass transfer (thermal/chemical relaxation) are **not** part of this paper's relaxation step [they appear in later works, e.g. Saurel–Petitpas–Abgrall 2008 / Zein et al. 2010 / Pelanti–Shyue 2014]. Shock treatment in mixtures uses Petitpas et al. (2007) "artificial heat exchange" (Section 5). The time integration is an explicit Godunov hyperbolic step followed by an instantaneous relaxation step (Strang/Lie-type operator splitting), **not** an IMEX or TV split.

---

## 1. Motivation (Sections 1–2.1)

### 1.1 The 5-equation (Kapila et al.) model (II.1)

Two fluids, single velocity, single pressure:

$$
\begin{aligned}
&\frac{\partial \alpha_1}{\partial t} + u\frac{\partial \alpha_1}{\partial x}
 = \frac{\rho_2 c_2^2-\rho_1 c_1^2}{\dfrac{\rho_2 c_2^2}{\alpha_2}+\dfrac{\rho_1 c_1^2}{\alpha_1}}\,\frac{\partial u}{\partial x},\\
&\partial_t(\alpha\rho)_1 + \partial_x\big((\alpha\rho)_1 u\big)=0,\qquad
\partial_t(\alpha\rho)_2 + \partial_x\big((\alpha\rho)_2 u\big)=0,\\
&\partial_t(\rho u) + \partial_x(\rho u^2+p)=0,\qquad
\partial_t(\rho E) + \partial_x\big((\rho E+p)u\big)=0,
\end{aligned}\tag{II.1}
$$

$E=e+\tfrac12u^2$; $\alpha$ volume fraction, $\rho$ mixture density, $u$ velocity, $p$ mixture pressure, $e$ mixture internal energy.

Mixture internal energy:
$$ e = Y_1 e_1(\rho_1,p) + Y_2 e_2(\rho_2,p),\qquad Y_k=\frac{(\alpha\rho)_k}{\rho},\quad \rho=(\alpha\rho)_1+(\alpha\rho)_2. \tag{II.2}$$

Each fluid has a convex EOS $e_k=e_k(\rho_k,p)$, sound speed $c_k=c_k(\rho_k,p)$. Mixture pressure from solving (II.2).

**Stiffened gas EOS:**
$$ p_k = (\gamma_k-1)\rho_k e_k - \gamma_k p_{\infty k}. \tag{II.3}$$

Resulting mixture EOS (used later for the reinitialisation step):
$$ p(\rho,e,\alpha_1,\alpha_2)=\frac{\rho e-\left(\dfrac{\alpha_1\gamma_1p_{\infty1}}{\gamma_1-1}+\dfrac{\alpha_2\gamma_2p_{\infty2}}{\gamma_2-1}\right)}{\dfrac{\alpha_1}{\gamma_1-1}+\dfrac{\alpha_2}{\gamma_2-1}}. \tag{II.4}$$

Phase entropies satisfy $ds_k/dt=0$, $k=1,2$ (smooth flow).

Shock relations for the 5-eq. model (from ref. [44], Saurel–Le Métayer–Massoni–Gavrilyuk 2007):
$$
Y_k=Y_k^0,\quad \rho(u-\sigma)=\rho^0(u^0-\sigma)=m,\quad p-p^0+m^2(v-v^0)=0,\quad
e_k-e_k^0+\frac{p+p^0}{2}(v_k-v_k^0)=0, \tag{II.5}
$$
$\sigma$ shock speed, superscript 0 = unshocked state.

### 1.2 Difficulties of the 5-eq. model (listed by authors)
- Non-conservative; numerical convergence to the exact shock solution is very hard; cell averages of non-conservative variables lack meaning; refs [43,37] needed a relaxation procedure + artificial heat exchange in shock layers.
- Volume-fraction positivity under shocks/strong rarefactions (liquid–gas: liquid pressure can go negative, then gas sound speed fails) — frequent in cavitation.
- Wood sound speed $\frac{1}{\rho c_{eq}^2}=\frac{\alpha_1}{\rho_1c_1^2}+\frac{\alpha_2}{\rho_2c_2^2}$ is **non-monotonic** in $\alpha$ (Fig. 1: water–air at atmospheric conditions, $c_{eq}$ drops to ~tens of m/s in the mixture).

Illustration (Fig. 2): 1 m shock tube, interface at $x=0.8$ m; water $\rho=1000$, $\gamma=4.4$, $p_\infty=6\times10^8$ Pa; air $\rho=10$, $\gamma=1.4$, $p_\infty=0$; left $\alpha_{air}=10^{-6}$, $p=10^9$ Pa; right reversed fractions, $p=10^5$ Pa; $u=0$; $t=220\,\mu$s; 1000 cells; Lagrange–relaxation method of Petitpas et al. [37] vs exact. Mach number oscillation in the diffuse interface (because of non-monotone $c_{eq}$).

Fig. 3: advection of a water–air interface at 50 m/s, first-order: equilibrium $c_{eq}$ dips well below both pure-fluid values in the diffusion zone, frozen $c_f^2=Y_{water}c_{water}^2+Y_{air}c_{air}^2$ is monotone. Fig. 4: schematic of the acoustic delay $\tau$ produced by low $c_{eq}$ in the diffuse zone.

Goals: volume-fraction positivity; cavitating flows (Mach 0.01–100); Riemann problem easy; correct wave transmission through diffuse interfaces; unstructured meshes.

---

## 2. The 6-equation model (Section 2.2)

Obtained from Baer–Nunziato in the stiff **velocity** relaxation limit only (first reduced model of Kapila et al.); pressure non-equilibrium kept. Authors stress: a *step-model* to solve the 5-eq. model, not a physical model per se.

### 2.1 Flow model (II.6)
$$
\begin{aligned}
&\frac{\partial\alpha_1}{\partial t}+u\frac{\partial\alpha_1}{\partial x}=\mu(p_1-p_2),\\
&\frac{\partial\alpha_1\rho_1}{\partial t}+\frac{\partial\alpha_1\rho_1u}{\partial x}=0,\qquad
\frac{\partial\alpha_2\rho_2}{\partial t}+\frac{\partial\alpha_2\rho_2u}{\partial x}=0,\\
&\frac{\partial\rho u}{\partial t}+\frac{\partial\big(\rho u^2+(\alpha_1p_1+\alpha_2p_2)\big)}{\partial x}=0,\\
&\frac{\partial\alpha_1\rho_1e_1}{\partial t}+\frac{\partial\alpha_1\rho_1e_1u}{\partial x}+\alpha_1p_1\frac{\partial u}{\partial x}=-p_I\,\mu(p_1-p_2),\\
&\frac{\partial\alpha_2\rho_2e_2}{\partial t}+\frac{\partial\alpha_2\rho_2e_2u}{\partial x}+\alpha_2p_2\frac{\partial u}{\partial x}=+p_I\,\mu(p_1-p_2).
\end{aligned}\tag{II.6}
$$

- $\mu>0$: pressure relaxation coefficient (taken $\to+\infty$ numerically).
- **Interfacial pressure** (asymptotic limit of the symmetric 7-eq. model of Saurel et al. [42] at equal velocities):
$$ p_I=\frac{Z_2p_1+Z_1p_2}{Z_1+Z_2},\qquad Z_k=\rho_kc_k \text{ (acoustic impedance)}. $$
- Interfacial velocity: single velocity $u$ (the $u\,\partial_x\alpha_1$ transport term); there is no separate $u_I$.
- Mixture pressure: $p=\alpha_1p_1+\alpha_2p_2$. Mixture density $\rho=\alpha_1\rho_1+\alpha_2\rho_2$.

### 2.2 Redundant mixture total energy (II.7)
Summing the internal-energy eqs. with mass and momentum:
$$
\frac{\partial\rho\left(Y_1e_1+Y_2e_2+\tfrac12u^2\right)}{\partial t}
+\frac{\partial u\left(\rho\left(Y_1e_1+Y_2e_2+\tfrac12u^2\right)+(\alpha_1p_1+\alpha_2p_2)\right)}{\partial x}=0. \tag{II.7}
$$
Important numerically: used to correct the inaccurate non-conservative internal-energy equations in the presence of shocks.

### 2.3 Entropy and sound speed
$$
\alpha_1\rho_1T_1\frac{ds_1}{dt}=\mu(p_1-p_2)^2\frac{Z_1}{Z_1+Z_2},\qquad
\alpha_2\rho_2T_2\frac{ds_2}{dt}=\mu(p_1-p_2)^2\frac{Z_2}{Z_1+Z_2},
$$
so mixture entropy $s=Y_1s_1+Y_2s_2$ is non-decreasing (this is why $p_I$ is impedance-weighted).

Frozen mixture sound speed (monotone in $Y$, $\alpha$):
$$ c_f^2 = Y_1c_1^2+Y_2c_2^2 .$$
Strictly hyperbolic with waves $u-c_f,\ u,\ u+c_f$.

### 2.4 Shock relations (Section 2.2.2)
Authors' stance: even with known shock relations a non-conservative system's numerical solution hardly converges to them, so precise relations aren't needed (the 6-eq. model only serves to approximate the 5-eq. one). But admissibility requirements (from [44]): mixture energy conservation; tangency of mixture Hugoniot and isentrope; correct single-phase limit; symmetry; entropy production.

Mass: $\alpha_1\rho_1(u-\sigma)=\alpha_1^0\rho_1^0(u^0-\sigma)=m_1$, $\alpha_2\rho_2(u-\sigma)=\alpha_2^0\rho_2^0(u^0-\sigma)=m_2$, $m=m_1+m_2$.
Momentum: $p-p^0+m^2(v-v^0)=0$ with $p=\alpha_1p_1+\alpha_2p_2$.
Mixture energy: $e-e^0+\frac{p+p^0}{2}(v-v^0)=0$, $e=Y_1e_1+Y_2e_2$, $v=Y_1v_1+Y_2v_2$, $v_k=1/\rho_k$.
Volume fraction (no relaxation): $\alpha_1=\alpha_1^0$.
**Proposed phasic jump conditions:**
$$ e_k-e_k^0+\frac{p_k+p_k^0}{2}(v_k-v_k^0)=0. \tag{II.8}$$
Check: summing $Y_k\times$(II.8) and using $Y_k=\alpha_kv/v_k$ and $\alpha_k=\alpha_k^0$ gives $e-e^0+\frac{\alpha_1p_1+\alpha_2p_2+\alpha_1p_1^0+\alpha_2p_2^0}{2}(v-v^0)=0$, i.e. mixture energy jump. Tangency holds because $\alpha$ is constant across acoustic waves (no relaxation) and phasic Hugoniots are tangent to phasic isentropes. Entropy: each phase on its own Hugoniot ⇒ mixture entropy increases.

### 2.5 Asymptotic limit (Section 2.2.3)
$\mu\to\infty$ recovers the 5-eq. model: proof in Appendix B (see below).

---

## 3. Numerical method (Section 3)

### 3.0 The 7-equation system actually solved (III.1)
$$
\begin{aligned}
&\partial_t\alpha_1+u\,\partial_x\alpha_1=\mu(p_1-p_2),\\
&\partial_t(\alpha_1\rho_1e_1)+\partial_x(\alpha_1\rho_1e_1u)+\alpha_1p_1\,\partial_xu=-p_I\mu(p_1-p_2),\\
&\partial_t(\alpha_2\rho_2e_2)+\partial_x(\alpha_2\rho_2e_2u)+\alpha_2p_2\,\partial_xu=p_I\mu(p_1-p_2),\\
&\partial_t(\alpha_1\rho_1)+\partial_x(\alpha_1\rho_1u)=0,\quad \partial_t(\alpha_2\rho_2)+\partial_x(\alpha_2\rho_2u)=0,\\
&\partial_t(\rho u)+\partial_x\big(\rho u^2+(\alpha_1p_1+\alpha_2p_2)\big)=0,\\
&\partial_t\Big(\rho(Y_1e_1+Y_2e_2+\tfrac12u^2)\Big)+\partial_x\Big(u\big(\rho(Y_1e_1+Y_2e_2+\tfrac12u^2)+(\alpha_1p_1+\alpha_2p_2)\big)\Big)=0,
\end{aligned}\tag{III.1}
$$
with $p_I=\frac{Z_2p_1+Z_1p_2}{Z_1+Z_2}$, EOS $e_k=e_k(\rho_k,p_k)$, and approximate shock relations (II.8).

### 3.1 Acoustic (linearised) Riemann solver (Section 3.1.1)
Without relaxation and without the redundant energy eq., primitive form $W_t+A(W)W_x=0$, $W=(\alpha_1,s_1,s_2,u,p_1,p_2)^T$:
$$
A(W)=\begin{pmatrix}
u&0&0&0&0&0\\
0&u&0&0&0&0\\
0&0&u&0&0&0\\
\frac{p_1-p_2}{\rho}&0&0&u&\frac{\alpha_1}{\rho}&\frac{\alpha_2}{\rho}\\
0&0&0&\rho_1c_1^2&u&0\\
0&0&0&\rho_2c_2^2&0&u
\end{pmatrix}.
$$
Eigenvalues: $\lambda_0=u$ (×4), $\lambda_{1,2}=u\mp c$ with
$$ c^2=Y_1c_1^2+Y_2c_2^2. \tag{III.2}$$

Characteristic relations:
- along $dx/dt=u$: $d\alpha_1/dt=0,\ ds_1/dt=0,\ ds_2/dt=0$;
- along $dx/dt=u-c$: $-\frac{p_1-p_2}{\rho c}\frac{d\alpha_1}{dt}+\frac{du}{dt}-\frac{\alpha_1}{\rho c}\frac{dp_1}{dt}-\frac{\alpha_2}{\rho c}\frac{dp_2}{dt}=0$;
- along $dx/dt=u+c$: $\frac{p_1-p_2}{\rho c}\frac{d\alpha_1}{dt}+\frac{du}{dt}+\frac{\alpha_1}{\rho c}\frac{dp_1}{dt}+\frac{\alpha_2}{\rho c}\frac{dp_2}{dt}=0$.

With $Z=\rho c$ frozen across each acoustic wave: right-facing wave: $\alpha_{1R}^*=\alpha_{1R}$, $s_{kR}^*=s_{kR}$, $(\alpha_1p_1+\alpha_2p_2)_R^*-Z_Ru_R^*=(\alpha_1p_1+\alpha_2p_2)_R-Z_Ru_R$, $Z_R=\rho_Rc_R$; left-facing: analogous with $+Z_Lu$. Interface conditions $(\alpha_1p_1+\alpha_2p_2)^*_L=(\cdot)^*_R=p^*$, $u^*_L=u^*_R=u^*$ give
$$
u^*=\frac{Z_Lu_L+Z_Ru_R+p_L-p_R}{Z_L+Z_R},\qquad
p^*=\frac{Z_Rp_L+Z_Lp_R+Z_RZ_L(u_L-u_R)}{Z_L+Z_R}, \tag{III.3}
$$
with $p=\alpha_1p_1+\alpha_2p_2$, $Z=\rho c$, $\rho=\alpha_1\rho_1+\alpha_2\rho_2$, $c^2=Y_1c_1^2+Y_2c_2^2$. Star phase densities from the entropy jumps (isentropes). Good for subsonic/weak-shock flows and boundary conditions; HLLC preferred (arbitrary shocks, positive, arbitrary convex EOS).

### 3.2 HLLC-type solver (Section 3.1.2)
Davis wave speeds:
$$ S_R=\max(u_L+c_L,\,u_R+c_R),\qquad S_L=\min(u_L-c_L,\,u_R-c_R), $$
with $c$ the frozen mixture sound speed (III.2). Contact speed (HLL estimate):
$$ S_M=\frac{(\rho u^2+p)_L-(\rho u^2+p)_R-S_L(\rho u)_L+S_R(\rho u)_R}{(\rho u)_L-(\rho u)_R-S_L\rho_L+S_R\rho_R}, $$
mixture $\rho$, mixture $p=\alpha_1p_1+\alpha_2p_2$.
Star states:
$$
(\alpha_k\rho_k)^*_R=(\alpha_k\rho_k)_R\frac{S_R-u_R}{S_R-S_M},\qquad
(\alpha_k\rho_k)^*_L=(\alpha_k\rho_k)_L\frac{S_L-u_L}{S_L-S_M},
$$
$$
p^*=p_R+\rho_Ru_R(u_R-S_R)-\rho_R^*S_M(S_M-S_R),\qquad \rho_R^*=\sum_k(\alpha_k\rho_k)_R^*,
$$
$$
E_R^*=\frac{\rho_RE_R(u_R-S_R)+p_Ru_R-p^*S_M}{\rho_R^*(S_M-S_R)},\qquad
E_L^*=\frac{\rho_LE_L(u_L-S_L)+p_Lu_L-p^*S_M}{\rho_L^*(S_M-S_L)},\qquad E=Y_1e_1+Y_2e_2+\tfrac12u^2.
$$
[Transcription note: the paper writes $E_R^*$ with numerator "$\rho_RE_R(u_R-S_R)+p_Ru_R-p^*S_M$"; this is the standard HLLC energy relation.]
Volume fraction unchanged across acoustic waves: $\alpha_{kR}^*=\alpha_{kR}$, $\alpha_{kL}^*=\alpha_{kL}$. Hence phase densities
$$ \rho^*_{kR}=\rho_{kR}\frac{u_R-S_R}{S_M-S_R} $$
(and similarly for L). **Phase internal energies via phasic Hugoniot (II.8)**; for SG:
$$
p_k^*(\rho_k^*)=(p_k+p_{\infty k})\frac{(\gamma_k-1)\rho_k-(\gamma_k+1)\rho_k^*}{(\gamma_k-1)\rho_k^*-(\gamma_k+1)\rho_k}-p_{\infty k},
$$
then $e^*_k=e_k(p^*_k,\rho^*_k)$ from EOS. (These star $e_k$ are used for the fluxes $(\alpha\rho e)_k^*$ of the non-conservative energy equations.)

### 3.3 First-order Godunov-type scheme (Section 3.2.1)
Conservative part ($U=((\alpha\rho)_1,(\alpha\rho)_2,\rho u,\rho E)^T$, $F=((\alpha\rho)_1u,(\alpha\rho)_2u,\rho u^2+p,(\rho E+p)u)^T$, $E=Y_1e_1+Y_2e_2+\frac12u^2$, $p=\alpha_1p_1+\alpha_2p_2$):
$$ U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}\big(F^*(U_i^n,U_{i+1}^n)-F^*(U_{i-1}^n,U_i^n)\big). $$
Volume fraction (Godunov for advection equations):
$$ \alpha_{1i}^{n+1}=\alpha_{1i}^n-\frac{\Delta t}{\Delta x}\Big((u\alpha_1)^*_{i+1/2}-(u\alpha_1)^*_{i-1/2}-\alpha^n_{1i}\big(u^*_{i+1/2}-u^*_{i-1/2}\big)\Big). $$
Guarantees volume-fraction positivity during the hyperbolic step (star values are upwinded by $S_M$ sign).
Non-conservative internal energies — simplest approx., product $(\alpha p)_{ki}^n$ frozen over the step:
$$ (\alpha\rho e)_{ki}^{n+1}=(\alpha\rho e)_{ki}^n-\frac{\Delta t}{\Delta x}\Big((\alpha\rho eu)^*_{k,i+1/2}-(\alpha\rho eu)^*_{k,i-1/2}+(\alpha p)^n_{ki}\big(u^*_{i+1/2}-u^*_{i-1/2}\big)\Big). $$
Inaccuracy in internal energies "not so crucial": they only estimate phase pressures at the end of the hyperbolic step before relaxation; relaxation gives a first correction; a second correction uses the mixture total energy.

### 3.4 Uniform-flow (Abrgall) test, analytic (Section 3.2.2)
For uniform $u$, $p_1=p_2=p$: mass ⇒ $\rho_i^{n+1}=\rho_i^n-\frac{\Delta t}{\Delta x}u(\rho^*_{i+1/2}-\rho^*_{i-1/2})$; momentum ⇒ $(\rho u)^{n+1}_i=u\rho^{n+1}_i$ so $u$ stays uniform. Internal energy with SG, $\rho_ke_k=\frac{p_k+\gamma_kp_{\infty k}}{\gamma_k-1}$:
$$ (\alpha(p+\gamma p_\infty))^{n+1}_{ki}=(\alpha(p+\gamma p_\infty))^n_{ki}-\frac{\Delta t}{\Delta x}u\big((\alpha(p+\gamma p_\infty))^*_{k,i+1/2}-(\cdot)^*_{k,i-1/2}\big), $$
and combined with the discrete $\alpha$ update $\alpha^{n+1}_{ki}=\alpha^n_{ki}-\frac{\Delta t}{\Delta x}u(\alpha^*_{k,i+1/2}-\alpha^*_{k,i-1/2})$ gives $p^{n+1}_{ki}=p$. **Interface conditions are preserved exactly** (needs the volume-fraction and energy discretisations to be consistent — same upwind star states). For Mie–Grüneisen $\rho_ke_k=\frac{p_k+\gamma_kp_{\infty k}(\rho_k)}{\gamma_k-1}$ the same is observed in practice since Godunov mass/volume fraction updates give a prolonged density field at the interface (locally reduces to SG).

### 3.5 Relaxation step (Section 3.3) — stiff pressure relaxation
Solve the ODE system with $\mu\to+\infty$:
$$
\partial_t\alpha_1=\mu(p_1-p_2),\quad \partial_t(\alpha_1\rho_1e_1)=-p_I\mu(p_1-p_2),\quad \partial_t(\alpha_2\rho_2e_2)=p_I\mu(p_1-p_2),
$$
$$ \partial_t(\alpha_1\rho_1)=\partial_t(\alpha_2\rho_2)=\partial_t(\rho u)=\partial_t(\rho E)=0. $$
Manipulation gives
$$ \frac{\partial e_k}{\partial t}+p_I\frac{\partial v_k}{\partial t}=0,\quad k=1,2. $$
Integral form:
$$ e_k-e_k^0+\hat p_{Ik}(v_k-v_k^0)=0,\qquad \hat p_{Ik}=\frac{1}{v_k-v_k^0}\int_0^{\Delta t}p_I\frac{\partial v_k}{\partial t}\,dt. $$
Summing: $Y_1e_1-Y_1e_1^0+Y_2e_2-Y_2e_2^0+\hat p_{I1}(Y_1v_1-Y_1v_1^0)+\hat p_{I2}(Y_2v_2-Y_2v_2^0)=0$; mixture mass gives $(Y_1v_1-Y_1v_1^0)+(Y_2v_2-Y_2v_2^0)=0$; so
$e-e^0+(\hat p_{I1}-\hat p_{I2})(Y_1v_1-Y_1v_1^0)=0$. Mixture energy conservation ⇒ **$\hat p_{I1}=\hat p_{I2}=\hat p_I$**. Possible estimates: $\hat p_I=p_I^0$ (initial interfacial pressure) or $\hat p_I=p$ (relaxed pressure); both compatible with entropy inequality [43]; difference negligible in practice (illustrated in results).

Algebraic system (3 unknowns $v_1,v_2,p$):
$$ e_k(p,v_k)-e_k^0(p_k^0,v_k^0)+\hat p_I(v_k-v_k^0)=0,\quad k=1,2, $$
closed by saturation $\sum_k\alpha_k=1$, i.e. $\sum_k(\alpha\rho)_kv_k=1$ ($(\alpha\rho)_k$ constant during relaxation).
For SG this gives explicitly
$$ v_k(p)=v_k^0\,\frac{p_k^0+\gamma_kp_{\infty k}+(\gamma_k-1)\hat p_I}{p+\gamma_kp_{\infty k}+(\gamma_k-1)\hat p_I}, $$
and a single scalar equation for $p$:
$$ \sum_k(\alpha\rho)_k\,v_k(p)=1. \tag{III.4}$$
Solve with Newton (recommended in summary). Then $v_k$, $\alpha_k=(\alpha\rho)_kv_k$.
[Implementation remark (mine): with $\hat p_I=p_I^0$ fixed, (III.4) is monotone decreasing in $p$ for $p>-\gamma_kp_{\infty k}-(\gamma_k-1)\hat p_I$; for two SG fluids it can be multiplied out to a quadratic in $p$ — the paper does not write that quadratic here [check later pages/appendix]. With $\hat p_I=p$ the equation is still scalar but nonlinear.]

In Miller–Puckett [32] the relaxed pressure is used directly; no guarantee it agrees with mixture EOS/mixture energy ⇒ next correction.

### 3.6 Reinitialisation step (Section 3.4) — the "+1" total energy correction
With the relaxed volume fractions, compute the mixture pressure from the **mixture total energy** (conservative, accurate at shocks) through the mixture EOS (II.4):
$$ p(\rho,e,\alpha_1,\alpha_2)=\frac{\rho e-\left(\frac{\alpha_1\gamma_1p_{\infty1}}{\gamma_1-1}+\frac{\alpha_2\gamma_2p_{\infty2}}{\gamma_2-1}\right)}{\frac{\alpha_1}{\gamma_1-1}+\frac{\alpha_2}{\gamma_2-1}},\qquad e=E-\tfrac12u^2 \text{ from the conservative } \rho E. $$
Valid in pure fluids and in the diffuse zone; guarantees correct/conservative wave dynamics on both sides of the interface. Then reset phase internal energies:
$$ e_k=e_k(p,\alpha_k\rho_k,\alpha_k)\quad[\text{i.e. } e_k(p,\rho_k),\ \rho_k=(\alpha\rho)_k/\alpha_k]. \tag{III.5}$$

### 3.7 Algorithm summary (Section 3.5)
1. At each cell boundary solve the Riemann problem of (III.1) (HLLC recommended).
2. Evolve all variables with the Godunov-type scheme of 3.2 (incl. $\alpha_1$, $(\alpha\rho e)_k$, and $\rho E$).
3. Relaxation: solve (III.4) for the relaxed pressure and especially volume fraction (Newton).
4. Compute mixture pressure with (II.4) from total energy.
5. Reset internal energies with their EOS (III.5).
6. Next time step.

---

## 4. Tests and validations (Section 4)

### 4.1 Advection of an interface in uniform $p$, $u$ (Fig. 5)
1 m tube, discontinuity at $x=0.5$ m, $u=100$ m/s, $p=10^5$ Pa. Left: water $\rho=1000$, $\gamma=4.4$, $p_\infty=6\times10^8$ Pa, $\alpha_{water}=1-\varepsilon$; right: air $\rho=10$, $\gamma=1.4$, $p_\infty=0$, $\alpha_{water}=\varepsilon$, $\varepsilon=10^{-8}$. $t=2.79$ ms, 200 cells, second order (Superbee limiter, App. C). Excellent agreement, oscillation free; Mach number with Wood sound speed shows a spike in diffuse zone, frozen-c Mach monotone. Relaxation plays no role here (flow is in equilibrium) — the interface-condition preservation is due to the consistent Godunov discretisation (Section 3.2.2).

### 4.2 Shock tube with Mie–Grüneisen (Cochran–Chan) EOS (Figs. 6–7)
Single fluid (liquid nitromethane) treated as two-phase to cure $p_\infty(\rho)$ nonlinearity oscillations of Godunov [43]. High-pressure chamber 20 GPa, low-pressure 0.2 MPa; densities 1134 and 1200 kg/m³; $\alpha_1=1-\varepsilon$ left, $\alpha_1=\varepsilon$ right, $\varepsilon=10^{-8}$; same EOS for both "phases":
$$ p(\rho,e)=\rho\Gamma(e-e_k(\rho))+p_k(\rho),\quad
e_k(\rho)=\frac{A_1}{\rho_{ref}(E_1-1)}\Big(\frac{\rho}{\rho_{ref}}\Big)^{E_1-1}-\frac{A_2}{\rho_{ref}(E_2-1)}\Big(\frac{\rho}{\rho_{ref}}\Big)^{E_2-1},\quad
p_k(\rho)=A_1\Big(\frac{\rho}{\rho_{ref}}\Big)^{E_1}-A_2\Big(\frac{\rho}{\rho_{ref}}\Big)^{E_2}.
$$
$\Gamma=1.19$, $\rho_{ref}=1134$, $A_1=0.819181\times10^9$ Pa, $A_2=1.50835\times10^9$ Pa, $E_1=4.52969$, $E_2=1.42144$. $t=67\,\mu$s, 500 cells (Fig. 6), exact Euler solution matched; Fig. 7 magnified $p$ (~9.82 GPa) and $u$ (~1770 m/s) around contact: oscillation-free. [Interface location not stated explicitly in the text I read; figure suggests 0.5 m — uncertain.]

### 4.3 Water–air shock tubes
#### 4.3.1 Moderate pressure ratio, high density ratio (Figs. 8–9)
1 m tube, interface at $x=0.75$ m. Water $\rho=1000$, $\gamma=4.4$, $p_\infty=6\times10^8$; air $\rho=1$ kg/m³, $\gamma=1.4$, $p_\infty=0$. Left: $\alpha_{air}=10^{-6}$, $p=1$ GPa; right reversed fractions, $p=0.1$ MPa; $u=0$. $t=240\,\mu$s; 1000 cells (Fig. 8, 2nd order van Leer limiter) and 100 cells (Fig. 9). Density ratio 1000, pressure ratio 10 000. Excellent agreement with exact Euler solution.
Text after Fig. 9: in this and all subsequent tests strong pressure waves propagate, so relaxation terms in the $\alpha$ and energy equations and the pressure relaxation step become important; robustness and convergence (to the weak solution in the unsteady build-up) are improved by pressure relaxation. 100-cell run (Fig. 9) shows more diffusion but "good agreement".

#### 4.3.2 Extreme conditions (Fig. 10)
Same as 4.3.1 but left pressure $10^{12}$ Pa (1 TPa = 10 Mbar) and air density 10 kg/m³ (caption: density ratio 100, pressure ratio $10^7$). $t=8.3\,\mu$s, 1000 cells. Velocity reaches ~32 000 m/s, Mach ~2.7. Shows robustness/convergence vs exact Euler solution.

### 4.4 Influence of $\hat p_I$ (Fig. 11)
Fig. 8 test re-run with $\hat p_I=p_I^0$ (left column) and $\hat p_I=p$ (right column): no visible difference in $p$, $u$, $\rho$.

### 4.5 Cavitation (expansion) test (Fig. 12)
1 m tube, liquid water at atmospheric pressure, $\rho=1000$ kg/m³, small air fraction $\alpha_{air}=10^{-2}$ everywhere; velocity discontinuity at $x=0.5$ m: $u=-100$ m/s left, $+100$ m/s right. $t=1.85$ ms, 1000 uniform cells. EOS: [not restated; presumably the water/air SG parameters above — uncertain]. Strong rarefactions; liquid pressure would go negative; gas fraction increases → cavitation pocket, two interfaces appear dynamically (air fraction ≈ 1 between ~0.3–0.67 m, pressure ≈ 0). Excellent agreement with the exact solution of the 5-eq. model [37]. Cavitation here is a purely **mechanical relaxation** effect (no mass transfer).

### 4.6 2D shock–bubble interaction (Figs. 13–14)
Layes & Le Métayer [27] experiment: Mach 1.5 shock in air (atmospheric, $P=101325$ Pa) hitting a helium bubble. $\rho_{air}=1.29$, $\rho_{He}=0.167$ kg/m³, ideal gases $\gamma_{air}=1.4$, $\gamma_{He}=1.67$. Domain 0.3 m × 0.08 m, bubble diameter 0.04 m, high-pressure air chamber at right behind a diaphragm (shock travels leftward). Numerical schlieren vs experiment at $t\approx134/137$, $274/291$, $344/350$, $554/577\,\mu$s: good qualitative agreement. Mesh not stated in text [uncertain].

### 4.7 Cavitating Richtmyer–Meshkov instability (Figs. 15–17)
Domain 3 m × 1 m; 900 × 400 cells. Left nearly pure water, right nearly pure gas, separated by a circular-arc interface of radius 0.6 m centred at $(x,y)=(1.2,0.5)$ m. Both phases initial velocity $-200$ m/s; top, bottom, left boundaries are solid walls. Water $\rho=1000$, $\gamma=4.4$, $p_\infty=6\times10^8$ Pa; gas $\rho=100$ kg/m³, $\gamma=1.8$, $p_\infty=0$. Left chamber $\alpha_{gas}=10^{-6}$; right chamber $\alpha_{water}=10^{-6}$. $P_0=1$ atm (Fig. 15). Density contours at $t=0,2,3.1,4.8,6.4,8.6$ ms (Fig. 16); gas fraction at 8.6 ms (Fig. 17). Wall impact → shock → RMI, jet elongation creates expansion waves; cavitation pockets appear near the wall and in the jet core, modifying jet/spike shape.

---

## 5. Shocks in multiphase mixtures: artificial heat exchange (Section 5)

Needed **only** for shocks propagating through genuine mixtures (not for interfaces between nearly pure fluids or cavitation). Reason: in a smeared numerical shock the thermodynamic path differs from the Hugoniot (Fig. 18: epoxy, single phase — end state still correct thanks to conservation). In mixtures the extra degree of freedom ($\alpha$) means the EOS changes with each weak numerical shock; the volume fraction in the shock layer is wrong, hence the mixture EOS (II.4) is wrong, and errors accumulate: the end state is **not** on the mixture Hugoniot.

### 5.1 Epoxy–spinel, moderate pressure ratio (Fig. 19)
1 m tube, interface at $x=0.6$ m, same mixture both sides. Epoxy: $\rho=1185$ kg/m³, $\gamma=2.43$, $p_\infty=5.3\times10^9$ Pa. Spinel: $\rho=3622$ kg/m³, $\gamma=1.62$, $p_\infty=141\times10^9$ Pa. $\alpha_{epoxy}=0.5954$ ($\alpha_{spinel}=1-\alpha_{epoxy}$) both chambers. $p_L=10^{10}$ Pa, right atmospheric (pressure ratio $10^5$), all at rest. 500 cells, $t=80\,\mu$s. Converges to the exact 5-eq. solution (with (II.5) shock relations) **without** artificial heat exchange.

### 5.2 Extreme pressure ratio $2\times10^6$ (Fig. 20)
Same setup, $p_L=2\times10^{11}$ Pa [inferred from ratio $2\times10^6$ and the figure's $2\times10^6$ bar]. Numerical solution does **not** converge to the exact 5-eq. solution: wrong shocked mixture density, volume fraction, spinel density, Mach — due to incorrect partition of internal energy in the shock layer.

### 5.3 Artificial heat exchange in the 6-equation model
Extra pressure $\pm q$ in the internal energy equations:
$$
\frac{\partial\alpha_1\rho_1e_1}{\partial t}+\frac{\partial\alpha_1\rho_1e_1u}{\partial x}+(\alpha_1p_1+q)\frac{\partial u}{\partial x}=0,\qquad
\frac{\partial\alpha_2\rho_2e_2}{\partial t}+\frac{\partial\alpha_2\rho_2e_2u}{\partial x}+(\alpha_2p_2-q)\frac{\partial u}{\partial x}=0,
$$
$$ q=\eta\!\left(\frac{\partial u}{\partial x}\right)\nu(v),\qquad \eta=\begin{cases}1 & \partial u/\partial x<0\\ 0&\text{otherwise}\end{cases}, $$
$\nu(v)$ the heat-exchange function (active only in compressive shock layers). More convenient (and mesh-independent) rewriting:
$$
\frac{\partial\alpha_1\rho_1e_1}{\partial t}+\frac{\partial\big(\alpha_1\rho_1e_1+\eta(\tfrac{\partial u}{\partial x})\mu(v)\big)u}{\partial x}+\alpha_1p_1\frac{\partial u}{\partial x}=0,\qquad
\frac{\partial\alpha_2\rho_2e_2}{\partial t}+\frac{\partial\big(\alpha_2\rho_2e_2-\eta(\tfrac{\partial u}{\partial x})\mu(v)\big)u}{\partial x}+\alpha_2p_2\frac{\partial u}{\partial x}=0.
$$
($\mu(v)$ here is a *different* function from the relaxation coefficient — notation clash in paper.) $\mu(v)$ must be predetermined for a given mixture (method in Petitpas et al. [37]); depends on the initial mixture state and on the scheme's numerical shock smearing. For epoxy–spinel (Fig. 21): piecewise-linear table and fit
$$ \tilde\mu(v)=\exp\!\big(-5.64\times10^{7}\,v^2+5.34\times10^{3}\,v+25.6\big)\ \ [\text{J/m}^3] $$
for $v\in[2.65\times10^{-4},4.61\times10^{-4}]$ m³/kg, corresponding to piston velocity 0–4200 m/s and pressure 1–880 000 atm. [Coefficients read from caption "$-5.6410^7$, $5.3410^3$" — interpreted as $5.64\times10^7$, $5.34\times10^3$; plotted values ~0–12 GJ/m³, units GJ/m³ on axis — sign/units slightly uncertain.]

### 5.4 Rerun with artificial heat exchange (Fig. 22)
Fig. 20 problem, 500 cells, artificial heat exchange in the shock layer only: converges to the exact 5-eq. solution. The $\mu(v)$ of Fig. 21 gives converged results for any shock strength in 1–880 000 atm; mesh independence "guaranteed".

## 6. Conclusion (Section 6)
6-eq. relaxation model simplifies numerics of the Kapila model; algorithm simple/robust; extensions to solid–fluid coupling (Gavrilyuk et al. elastic model [17]) planned.

---

## Appendix A — Link with barotropic cavitation models

Barotropic simplification of the 5-eq. model (shocks absent/weak):
$$ \partial_ts_k+u\partial_xs_k=0\ (k=1,2),\quad \partial_t(\alpha\rho)_k+\partial_x((\alpha\rho)_ku)=0,\quad \partial_t(\rho u)+\partial_x(\rho u^2+p)=0, \tag{A.1}$$
closed by $p_1(\rho_1,s_1)=p_2(\rho_2,s_2)=p$ (A.2). With uniform entropies $s_k=s_k^0$ → three conservation laws (A.3) (two masses + momentum). SG isentrope = Tait EOS: $\frac{p_k+p_{\infty k}}{\rho_k^{\gamma_k}}=\frac{p_0+p_{\infty k}}{\rho_{0k}^{\gamma_k}}$,
$$ p_k=p_k^S(\rho_k)=\Big(\frac{\rho_k}{\rho_{0k}}\Big)^{\gamma_k}(p_0+p_{\infty k})-p_{\infty k}, \tag{A.4}$$
closure $p_1^S(\rho_1)=p_2^S(\rho_2)$ (A.5), i.e. scalar equation in $\alpha_1$:
$$ f(\alpha_1)=\Big(\frac{(\alpha\rho)_1}{\alpha_1\rho_{01}}\Big)^{\gamma_1}(p_0+p_{\infty1})-p_{\infty1}-\Big(\frac{(\alpha\rho)_2}{(1-\alpha_1)\rho_{02}}\Big)^{\gamma_2}(p_0+p_{\infty2})+p_{\infty2}=0. \tag{A.6}$$
Cavitation = mechanical relaxation (gas fraction growth), not mass transfer. Wood sound speed again non-monotone.

**Relaxation version (A.7):**
$$
\partial_t\alpha_1+u\partial_x\alpha_1=\mu\big(p_1^S(\rho_1)-p_2^S(\rho_2)\big),\quad
\partial_t(\alpha\rho)_k+\partial_x((\alpha\rho)_ku)=0,\quad
\partial_t(\rho u)+\partial_x(\rho u^2+\alpha_1p_1+\alpha_2p_2)=0, \tag{A.7}
$$
frozen sound speed $c_f^2=Y_1c_1^2+Y_2c_2^2$ (monotone); $\mu\to\infty$ recovers (A.3)+(A.5). Algorithm: HLLC Riemann solver for (A.7) without relaxation; Godunov update of $W=(\alpha_1,(\alpha\rho)_1,(\alpha\rho)_2,\rho u)$; relax by solving (A.6) with Newton; repeat.

## Appendix B — Asymptotic limit $\mu\to\infty$
Expansion $f=f^0+\varepsilon f^1$, $\mu=\mu_0/\varepsilon$. Lagrangian form ($d/dt=\partial_t+u\partial_x$):
$$ \frac{d\alpha_1}{dt}=\mu(p_1-p_2),\quad \alpha_1\rho_1\frac{de_1}{dt}+\alpha_1p_1\frac{\partial u}{\partial x}=-p_I\mu(p_1-p_2),\quad \alpha_2\rho_2\frac{de_2}{dt}+\alpha_2p_2\frac{\partial u}{\partial x}=p_I\mu(p_1-p_2). $$
With $e_1=e_1(\rho_1,p_1)$ and the phase mass eq. $\frac{d\rho_1}{dt}=-\frac{\rho_1}{\alpha_1}\mu(p_1-p_2)-\rho_1\frac{\partial u}{\partial x}$:
$$ \frac{dp_1}{dt}+\rho_1\,\frac{\frac{p_1}{\rho_1^2}-\frac{\partial e_1}{\partial\rho_1}\big|_{p_1}}{\frac{\partial e_1}{\partial p_1}\big|_{\rho_1}}\frac{\partial u}{\partial x}=-\frac{\rho_1}{\alpha_1}\,\frac{\frac{p_I}{\rho_1^2}-\frac{\partial e_1}{\partial\rho_1}\big|_{p_1}}{\frac{\partial e_1}{\partial p_1}\big|_{\rho_1}}\mu(p_1-p_2). $$
Defining
$$ c_1^2=\frac{\frac{p_1}{\rho_1^2}-\frac{\partial e_1}{\partial\rho_1}\big|_{p_1}}{\frac{\partial e_1}{\partial p_1}\big|_{\rho_1}},\qquad
c_{I1}^2=\frac{\frac{p_I}{\rho_1^2}-\frac{\partial e_1}{\partial\rho_1}\big|_{p_1}}{\frac{\partial e_1}{\partial p_1}\big|_{\rho_1}}, $$
$$ \frac{dp_1}{dt}+\rho_1c_1^2\frac{\partial u}{\partial x}=-\frac{\rho_1c_{I1}^2}{\alpha_1}\mu(p_1-p_2),\qquad
\frac{dp_2}{dt}+\rho_2c_2^2\frac{\partial u}{\partial x}=\frac{\rho_2c_{I2}^2}{\alpha_2}\mu(p_1-p_2). $$
(These pressure-evolution equations with relaxation source are the basis for any implicit/IMEX pressure-relaxation treatment.)
Order $1/\varepsilon$: $p_1^0=p_2^0=p^0$ ⇒ $p_I^0=p^0$, $c_{I1}^{0\,2}=c_1^{0\,2}$, $c_{I2}^{0\,2}=c_2^{0\,2}$. Order 0:
$$ \frac{dp^0}{dt}+\rho_1^0c_1^{0\,2}\frac{\partial u^0}{\partial x}=-\frac{\rho_1^0c_1^{0\,2}}{\alpha_1^0}\mu_0(p_1^1-p_2^1),\qquad
\frac{dp^0}{dt}+\rho_2^0c_2^{0\,2}\frac{\partial u^0}{\partial x}=\frac{\rho_2^0c_2^{0\,2}}{\alpha_2^0}\mu_0(p_1^1-p_2^1), $$
difference ⇒ $\mu_0(p_1^1-p_2^1)=\dfrac{\rho_2^0c_2^{0\,2}-\rho_1^0c_1^{0\,2}}{\frac{\rho_1^0c_1^{0\,2}}{\alpha_1^0}+\frac{\rho_2^0c_2^{0\,2}}{\alpha_2^0}}\dfrac{\partial u^0}{\partial x}$ [the paper writes $p_1^1-p_2^1=\ldots$ omitting $\mu_0$ — typo/absorbed], so
$$ \frac{d\alpha_1^0}{dt}=\frac{\rho_2^0c_2^{0\,2}-\rho_1^0c_1^{0\,2}}{\frac{\rho_1^0c_1^{0\,2}}{\alpha_1^0}+\frac{\rho_2^0c_2^{0\,2}}{\alpha_2^0}}\frac{\partial u^0}{\partial x}, $$
i.e. the Kapila 5-eq. model.

## Appendix C — Second order (MUSCL–Hancock in primitive variables)
Hyperbolic part (C.1) = (II.6) without relaxation. Primitive form, $W=(\alpha_1,\rho_1,\rho_2,u,p_1,p_2)^T$:
$$
\partial_t\alpha_1+u\partial_x\alpha_1=0,\ \ \partial_t\rho_k+u\partial_x\rho_k+\rho_k\partial_xu=0,\ \
\partial_tu+u\partial_xu+\frac{p_1-p_2}{\rho}\partial_x\alpha_1+\frac{\alpha_1}{\rho}\partial_xp_1+\frac{1-\alpha_1}{\rho}\partial_xp_2=0,\ \
\partial_tp_k+u\partial_xp_k+\rho_kc_k^2\partial_xu=0,
$$
$W_t+A(W)W_x=0$ (C.2) with
$$ A(W)=\begin{pmatrix}u&0&0&0&0&0\\0&u&0&\rho_1&0&0\\0&0&u&\rho_2&0&0\\ \frac{p_1-p_2}{\rho}&0&0&u&\frac{\alpha_1}{\rho}&\frac{1-\alpha_1}{\rho}\\0&0&0&\rho_1c_1^2&u&0\\0&0&0&\rho_2c_2^2&0&u\end{pmatrix}. $$
Steps: (C.1) slopes $\Delta_i^-=(W_i-W_{i-1})/\Delta x$, $\Delta_i^+=(W_{i+1}-W_i)/\Delta x$, limited $\Delta_i=\xi(\Delta_i^-,\Delta_i^+)$ (minmod, van Leer, Superbee). (C.2) $W_{i,L}=W_i-\frac{\Delta x}{2}\Delta_i$, $W_{i,R}=W_i+\frac{\Delta x}{2}\Delta_i$; half-step predictor
$$ W^{n+1/2}_{i,L,R}=W^n_{i,L,R}+\frac12\frac{\Delta t}{\Delta x}A(W_i^n)\big[W^n_{i,L}-W^n_{i,R}\big]. $$
(C.3) Riemann problems $F^*_{i-1/2}=F^*(W^{n+1/2}_{i-1,R},W^{n+1/2}_{i,L})$, $F^*_{i+1/2}=F^*(W^{n+1/2}_{i,R},W^{n+1/2}_{i+1,L})$, also giving $\alpha^*_{k,i\pm1/2}$, $u^*_{i\pm1/2}$, $(\alpha\rho e)^*_{k,i\pm1/2}$ (C.4). (C.4) Evolution — same formulas as first order with these star values, **the non-conservative products $\alpha_{1i}^n$ and $(\alpha p)^n_{ki}$ still taken at time $n$ cell-centred values**:
$$ \alpha_{1i}^{n+1}=\alpha_{1i}^n-\frac{\Delta t}{\Delta x}\Big((u\alpha_1)^*_{i+1/2}-(u\alpha_1)^*_{i-1/2}-\alpha^n_{1i}(u^*_{i+1/2}-u^*_{i-1/2})\Big), $$
$$ (\alpha\rho e)_{ki}^{n+1}=(\alpha\rho e)^n_{ki}-\frac{\Delta t}{\Delta x}\Big((\alpha\rho eu)^*_{k,i+1/2}-(\alpha\rho eu)^*_{k,i-1/2}+(\alpha p)^n_{ki}(u^*_{i+1/2}-u^*_{i-1/2})\Big). $$
[No CFL numbers are given anywhere in the paper — uncertain/unspecified.]

## Appendix D — Multi-D finite volume (structured/unstructured)
Conservative part $\partial_t\int_{V_i}U+\int_A\mathbf H\cdot\mathbf n\,dA=0$ (D.1), $U=((\alpha\rho)_1,(\alpha\rho)_2,\rho u,\rho v,\rho E)^T$, $F=((\alpha\rho)_1u,(\alpha\rho)_2u,\rho u^2+p,\rho uv,(\rho E+p)u)^T$, $G=((\alpha\rho)_1v,(\alpha\rho)_2v,\rho uv,\rho v^2+p,(\rho E+p)v)^T$, $E=Y_1e_1+Y_2e_2+\frac12u^2+\frac12v^2$, $p=\alpha_1p_1+\alpha_2p_2$. Polygon with segments $[A_s,A_{s+1}]$, length $L_s$, normal $\mathbf n_s=(\cos\theta_s,\sin\theta_s)$:
$$ U_i^{n+1}=U_i^n-\frac{\Delta t}{V_i}\sum_{s=1}^NL_s\big(F_s^*\cos\theta_s+G_s^*\sin\theta_s\big), $$
$$ \alpha_{k,i}^{n+1}=\alpha^n_{k,i}-\frac{\Delta t}{V_i}\sum_sL_s\Big[(u\alpha_k)^*_s\cos\theta_s+(v\alpha_k)^*_s\sin\theta_s-\alpha^n_{k,i}(u^*_s\cos\theta_s+v^*_s\sin\theta_s)\Big], $$
$$ (\alpha\rho e)^{n+1}_{k,i}=(\alpha\rho e)^n_{k,i}-\frac{\Delta t}{V_i}\sum_sL_s\Big[((\alpha\rho e)_ku)^*_s\cos\theta_s+((\alpha\rho e)_kv)^*_s\sin\theta_s+(\alpha p)^n_{k,i}(u^*_s\cos\theta_s+v^*_s\sin\theta_s)\Big]. $$
Riemann problems solved in the face-normal direction.

---

## 7. Implementation remarks (for the diffusionIMEX 6-eq. + TV-split project)

1. **What this paper gives vs. what we need.** It provides the 6-eq. model (II.6), $p_I$ closure, redundant total energy (II.7), HLLC with phasic Hugoniot star energies, the positivity-preserving Godunov discretisation of $u\partial_x\alpha$ and $\alpha_kp_k\partial_xu$, and **instantaneous pressure relaxation** (III.4) + **total-energy reinitialisation** (II.4)/(III.5). It contains **no thermal or chemical (Gibbs) relaxation**, no IMEX, no TV splitting, no low-Mach treatment — take those from Pelanti–Shyue (2014), Saurel et al. (2008, metastable liquids [45]), Zein et al. (2010), Toro–Vázquez-Cendón (2012).
2. **Consistency of the $\alpha$ and $(\alpha\rho e)_k$ discretisations is what makes the uniform-$u$,$p$ interface test exact** (Section 3.2.2): both use the same Riemann star states $\alpha^*$, $u^*$ and the cell-centred, time-$n$ values $\alpha^n_{1i}$, $(\alpha p)^n_{ki}$ multiplying $(u^*_{i+1/2}-u^*_{i-1/2})$. Keep this structure in any TV-split variant (the pressure subsystem must update $\alpha$, $(\alpha\rho e)_k$ with the *same* face velocities).
3. **Total-energy correction is essential at shocks**: after relaxation, recompute $p$ from conservative $\rho E$ via the mixture EOS (II.4) with relaxed $\alpha_k$, then overwrite $e_k=e_k(p,\rho_k)$. Phases end in exact pressure equilibrium $p_1=p_2=p$ by construction. This is the "+1" in "6+1". Note: the relaxed pressure from (III.4) and the reinitialised $p$ from (II.4) differ in general; the paper keeps the (II.4) one.
4. **Relaxation solve:** with SG and fixed $\hat p_I=p_I^0$, $v_k(p)$ is explicit and (III.4) is a monotone scalar equation; Newton is used. Multiplying out gives a quadratic for two phases (not written in the paper; derive yourself): with $A_k=(\alpha\rho)_kv_k^0(p_k^0+\gamma_kp_{\infty k}+(\gamma_k-1)\hat p_I)=\alpha_k^0(p_k^0+\gamma_kp_{\infty k}+(\gamma_k-1)\hat p_I)$, $B_k=\gamma_kp_{\infty k}+(\gamma_k-1)\hat p_I$, solve $\frac{A_1}{p+B_1}+\frac{A_2}{p+B_2}=1$, i.e. $p^2+(B_1+B_2-A_1-A_2)p+B_1B_2-A_1B_2-A_2B_1=0$, take the root with $p+B_k>0$. [My derivation, not in paper.] Choice $\hat p_I=p_I^0$ vs $\hat p_I=p$ is immaterial in practice (Fig. 11).
5. **Positivity:** the hyperbolic $\alpha$ update is a Godunov advection scheme → stays in $[0,1]$ under CFL; after relaxation $\alpha_k=(\alpha\rho)_kv_k(p)>0$ if $p+B_k>0$. Pure phases are approximated with $\varepsilon=10^{-8}$ (or $10^{-6}$) of the other phase; no special abnormal-case handling is described. Non-monotone Wood sound speed is avoided because waves are computed with the frozen $c_f$.
6. **Sound speed in Riemann solver** = frozen $c^2=Y_1c_1^2+Y_2c_2^2$ (III.2), not Wood.
7. **HLLC star phase energies** come from each phase's **Hugoniot (II.8)** (closed form for SG), not from isentropes — needed for admissible single-phase shock limits.
8. **Mixture shocks** with large pressure ratio (≳10^6) need the artificial heat-exchange function $\mu(v)$ calibrated per mixture and per scheme — a known weakness; irrelevant for interfaces between nearly pure fluids and cavitation.
9. **Mie–Grüneisen EOS:** the method removes the $p_\infty(\rho)$-induced oscillations by treating a single fluid as two "phases" across the contact (Section 4.2); mixture EOS (II.4) must then be replaced by the generalisation for that EOS [the paper does not write it out — uncertain how $p$ is reinitialised for CC EOS; presumably solving $\rho e=\sum\alpha_k\rho_ke_k(p,\rho_k)$ for $p$].
10. CFL numbers, boundary condition details, and 2D mesh for the shock–bubble test are not given in the paper.
