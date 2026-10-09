# Paper 6 — Reading note

**Title:** A mixture-energy-consistent six-equation two-phase numerical model for fluids with interfaces, cavitation and evaporation waves
**Authors:** Marica Pelanti (ENSTA ParisTech, Palaiseau, France), Keh-Ming Shyue (National Taiwan University, Taipei)
**Venue:** Journal of Computational Physics 259 (2014) 331–357. doi:10.1016/j.jcp.2013.12.003
(Received 31 Mar 2013, revised 30 Oct 2013, accepted 1 Dec 2013, online 10 Dec 2013.) The PDF ends with an "update to this article" (a corrigendum) — see end of note.

**Keywords:** multiphase compressible flow, mechanical relaxation, thermo-chemical relaxation, cavitation, phase transition, finite volume, wave propagation algorithms, Riemann solvers.

## 0. One-paragraph summary

Single-velocity 6-equation model of Saurel–Petitpas–Berry (SPB 2009, ref [9]) with stiff pressure relaxation, but written with **phasic total energies** $E_k=\mathcal E_k+\tfrac12\rho_k u^2$ instead of phasic internal energies. The sum of the two phasic total-energy equations is exactly the conservative mixture-total-energy equation, so a conservative discretisation of the *conservative part* of each phasic energy equation plus a *symmetric* ($\pm\Sigma$) discretisation of the non-conservative term automatically yields a conservative discrete mixture energy, and the relaxed pressure is then consistent with the mixture EOS ("mixture-energy-consistent"). No extra mixture-energy equation is needed (unlike SPB [9] and Zein et al. [17]). Heat and mass transfer are added as relaxation terms $\mathcal Q=\theta(T_2-T_1)$, $\dot m=\nu(g_2-g_1)$, with $\theta,\nu\to\infty$ only at interfaces (and, for $\nu$, only under metastable conditions $T_{liq}>T_{sat}$). Scheme: fractional steps — (1) homogeneous hyperbolic system by second-order wave-propagation (LeVeque, CLAWPACK) with HLLC (preferred) or Roe solver; (2) instantaneous pressure relaxation (algebraic, closed form quadratic for SG EOS); (3) if active, p–T relaxation or p–T–g relaxation, solved algebraically from equilibrium conditions (saturation curve, quadratic in $T$ etc.). Stiffened gas EOS for each phase.

---

## 1. Model (Section 2)

### 1.1 Phasic-internal-energy formulation (SPB [9]) — Eq. (1)

$$
\begin{aligned}
&\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mu(p_1-p_2), &(1a)\\
&\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=0, &(1b)\\
&\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=0, &(1c)\\
&\partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u)+\nabla(\alpha_1p_1+\alpha_2p_2)=0, &(1d)\\
&\partial_t(\alpha_1\mathcal E_1)+\nabla\cdot(\alpha_1\mathcal E_1\vec u)+\alpha_1p_1\nabla\cdot\vec u=-\mu p_I(p_1-p_2), &(1e)\\
&\partial_t(\alpha_2\mathcal E_2)+\nabla\cdot(\alpha_2\mathcal E_2\vec u)+\alpha_2p_2\nabla\cdot\vec u=\mu p_I(p_1-p_2). &(1f)
\end{aligned}
$$

- $\alpha_k$: volume fraction of phase $k=1,2$, $\alpha_1+\alpha_2=1$; $\rho_k$ phasic density; $p_k$ phasic pressure; $\mathcal E_k=\rho_k\varepsilon_k$ phasic internal energy per unit volume, $\varepsilon_k$ specific internal energy.
- $\rho=\alpha_1\rho_1+\alpha_2\rho_2$ mixture density, $\vec u$ single flow velocity.
- $\mu>0$ pressure relaxation parameter; interface pressure
  $$p_I=\frac{Z_2p_1+Z_1p_2}{Z_1+Z_2},\qquad Z_k=\rho_kc_k \ \text{(acoustic impedance)}.$$
  NB: since relaxation is taken with $\mu\to\infty$, the exact form of $p_I$ does not matter in the end (the relaxation step uses an algebraic closure; see §3.3 below).
- Instantaneous relaxation $\mu\to\infty$ → Kapila et al. 5-equation model [18]. Authors argue 6-eq + stiff relaxation is numerically easier than 5-eq: the 5-eq volume-fraction equation has a $K\nabla\cdot u$ non-conservative term depending on phasic impedances, hard for Riemann solvers, positivity of $\alpha$ hard to preserve in strong shocks/rarefactions.

**Closure.** EOS per phase $p_k=p_k(\mathcal E_k,\rho_k)$. Sound speed
$c_k=\sqrt{\kappa_k h_k+\chi_k}$, $h_k=(\mathcal E_k+p_k)/\rho_k$, $\kappa_k=\partial p_k(\mathcal E_k,\rho_k)/\partial\mathcal E_k$, $\chi_k=\partial p_k(\mathcal E_k,\rho_k)/\partial\rho_k$.

**Mixture sound speed of the 6-eq model (frozen):**
$$c=\sqrt{Y_1c_1^2+Y_2c_2^2},\qquad Y_k=\alpha_k\rho_k/\rho \tag{2}$$
Monotone in $\alpha$ (unlike Wood's speed of the 5-eq model) — no double sonic points in the diffuse zone; this is an advantage for approximate Riemann solvers.

A caloric law $T_k=T_k(p_k,\rho_k)$ is also needed when heat/mass transfer is included (incomplete vs complete EOS).

### 1.2 Stiffened gas EOS — Eq. (3)

$$
\begin{aligned}
p_k(\mathcal E_k,\rho_k)&=(\gamma_k-1)\mathcal E_k-\gamma_k\pi_k-(\gamma_k-1)\eta_k\rho_k, &(3a)\\
T_k(p_k,\rho_k)&=\frac{p_k+\pi_k}{C_{vk}\rho_k(\gamma_k-1)}, &(3b)\\
s_k(p_k,T_k)&=C_{vk}\log\frac{T_k^{\gamma_k}}{(p_k+\pi_k)^{\gamma_k-1}}+\eta'_k, &(3c)\\
g_k(p_k,T_k)&=(\gamma_kC_{vk}-\eta'_k)T_k-C_{vk}T_k\log\frac{T_k^{\gamma_k}}{(p_k+\pi_k)^{\gamma_k-1}}+\eta_k, &(3d)
\end{aligned}
$$
with $g_k=h_k-T_ks_k$ (Gibbs free energy / chemical potential), parameters $\gamma_k,\pi_k,\eta_k,C_{vk}$, constant $\eta'_k$.
For SG: $\kappa_k=\gamma_k-1$, $\chi_k=-(\gamma_k-1)\eta_k$, $c_k=\sqrt{\gamma_k\frac{p_k+\pi_k}{\rho_k}}$.

### 1.3 Mixture EOS — Eqs. (4)–(5)

Mixture specific internal energy $\varepsilon=Y_1\varepsilon_1+Y_2\varepsilon_2$; per unit volume $\mathcal E=\rho\varepsilon=\alpha_1\mathcal E_1+\alpha_2\mathcal E_2$. Under the isobaric assumption $p_1=p_2=p$:
$$\mathcal E=\alpha_1\mathcal E_1(p,\rho_1)+\alpha_2\mathcal E_2(p,\rho_2) \tag{4}$$
which defines implicitly $p=p(\mathcal E,\rho_1,\rho_2,\alpha_1)$. For SG EOS (explicit, ref [28]):
$$
p(\mathcal E,\rho_1,\rho_2,\alpha_1)=\frac{\mathcal E-(\alpha_1\rho_1\eta_1+\alpha_2\rho_2\eta_2)-\left(\frac{\alpha_1\gamma_1\pi_1}{\gamma_1-1}+\frac{\alpha_2\gamma_2\pi_2}{\gamma_2-1}\right)}{\frac{\alpha_1}{\gamma_1-1}+\frac{\alpha_2}{\gamma_2-1}}. \tag{5}
$$
Hyperbolicity: in the thermodynamically stable region the model (1) is hyperbolic (real eigenvalues, complete eigenvectors; App. A). (1) is the asymptotic limit of the 7-eq two-velocity Saurel–Abgrall model for instantaneous kinetic relaxation. Eigenvalue order is known a priori ($u-c,u,u,u,u,u+c$ in 1D) — easy wave decomposition.

### 1.4 Phasic-total-energy formulation (the paper's model) — Eq. (6)

$E_k=\mathcal E_k+\tfrac12\rho_k\vec u\cdot\vec u$ (phasic total energy per unit volume, note uses $\rho_k$ not $\alpha_k\rho_k$; the conserved-like variable is $\alpha_kE_k$).

$$
\begin{aligned}
&\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mu(p_1-p_2), &(6a)\\
&\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=0, &(6b)\\
&\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=0, &(6c)\\
&\partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u)+\nabla(\alpha_1p_1+\alpha_2p_2)=0, &(6d)\\
&\partial_t(\alpha_1E_1)+\nabla\cdot(\alpha_1E_1\vec u+\alpha_1p_1\vec u)+\Sigma(q,\nabla q)=-\mu p_I(p_1-p_2), &(6e)\\
&\partial_t(\alpha_2E_2)+\nabla\cdot(\alpha_2E_2\vec u+\alpha_2p_2\vec u)-\Sigma(q,\nabla q)=\mu p_I(p_1-p_2), &(6f)
\end{aligned}
$$
with non-conservative term
$$
\Sigma(q,\nabla q)=-\vec u\cdot\big(Y_2\nabla(\alpha_1p_1)-Y_1\nabla(\alpha_2p_2)\big)
=-\vec u\cdot\big((Y_2p_1+Y_1p_2)\nabla\alpha_1+\alpha_1Y_2\nabla p_1-\alpha_2Y_1\nabla p_2\big). \tag{6g}
$$
($q$ = vector of unknowns.) Sum of (6e)+(6f) gives the conservative mixture total energy equation, $E=\mathcal E+\tfrac12\rho u^2=\alpha_1E_1+\alpha_2E_2$:
$$\partial_tE+\nabla\cdot(E\vec u+\alpha_1p_1\vec u+\alpha_2p_2\vec u)=0. \tag{7}$$

*Check (my derivation, not in paper):* $\Sigma$ arises because phasic kinetic energy $\tfrac12\alpha_k\rho_k u^2$ evolves with $u\cdot\nabla(\alpha_1p_1+\alpha_2p_2)$ weighted by $Y_k$; $\alpha_1\rho_1 u\cdot D u/Dt=-Y_1 u\cdot\nabla P$ with $P=\alpha_1p_1+\alpha_2p_2$, and the phase-1 internal energy eq. contributes $\alpha_1p_1\nabla\cdot u$; combining: $-Y_1u\cdot\nabla P+\nabla\cdot(\alpha_1p_1u)-\alpha_1p_1\nabla\cdot u... $ gives $u\cdot\nabla(\alpha_1p_1)-Y_1u\cdot\nabla P = u\cdot(Y_2\nabla(\alpha_1p_1)-Y_1\nabla(\alpha_2p_2))=-\Sigma$. Consistent with (6g).

### 1.5 Model with heat and mass transfer — Eqs. (8)–(10)

$$
\begin{aligned}
&\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mu(p_1-p_2)+\frac{\dot m}{\rho_I}, &(8a)\\
&\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\dot m, &(8b)\\
&\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\dot m, &(8c)\\
&\partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u)+\nabla(\alpha_1p_1+\alpha_2p_2)=0, &(8d)\\
&\partial_t(\alpha_1E_1)+\nabla\cdot(\alpha_1E_1\vec u+\alpha_1p_1\vec u)+\Sigma(q,\nabla q)=-\mu p_I(p_1-p_2)+\mathcal Q+e_I\dot m, &(8e)\\
&\partial_t(\alpha_2E_2)+\nabla\cdot(\alpha_2E_2\vec u+\alpha_2p_2\vec u)-\Sigma(q,\nabla q)=\mu p_I(p_1-p_2)-\mathcal Q-e_I\dot m, &(8f)
\end{aligned}
$$
$$\mathcal Q=\theta(T_2-T_1),\qquad \dot m=\nu(g_2-g_1). \tag{9}$$
$\theta,\nu$: thermal and chemical relaxation parameters, infinite at selected locations, zero elsewhere:
$$
\theta=\begin{cases}\infty & \text{if }\epsilon_I\le\alpha_1\le1-\epsilon_I,\\ 0&\text{otherwise},\end{cases}\quad(10a)
\qquad
\nu=\begin{cases}\infty & \text{if }\epsilon_I\le\alpha_1\le1-\epsilon_I\text{ and }T_{liq}>T_{sat},\\ 0&\text{otherwise}.\end{cases}\quad(10b)
$$
$\epsilon_I$ identifies interface locations (e.g. $\epsilon_I=10^{-4}$). Phase 1 = liquid in the convention of (8a)/(8b) [the sign: $\dot m=\nu(g_2-g_1)>0$ when $g_2>g_1$, mass goes to phase 1; note $\dot m>0$ means phase 2→1]. Hmm: in later sections phase 1 is liquid and phase 2 vapour [confirm with tests].

Subcharacteristic condition proved for the 3 relaxation levels (p; p,T; p,T,g) in refs [33,34] (Flåtten & Lund; Lund).

**Saturation curve** from $g_1=g_2$ with SG EOS:
$$A_S+\frac{B_S}{T}+C_S\log T+D_S\log(p+\pi_1)-\log(p+\pi_2)=0, \tag{11a}$$
$$A_S=\frac{C_{p1}-C_{p2}+\eta'_2-\eta'_1}{C_{p2}-C_{v2}},\quad B_S=\frac{\eta_1-\eta_2}{C_{p2}-C_{v2}},\quad C_S=\frac{C_{p2}-C_{p1}}{C_{p2}-C_{v2}},\quad D_S=\frac{C_{p1}-C_{v1}}{C_{p2}-C_{v2}}. \tag{11b}$$
$C_{pk}=C_{vk}\gamma_k$. SG parameters fitted to the experimental saturation curve (Le Métayer et al. [36]). The theoretical $p$–$T$ curve is used in the algorithm (§4.2).

Interfacial density $\rho_I$ and interfacial specific total energy $e_I$: no need to specify them, thanks to infinite relaxation + algebraic treatment (Section 4). Symmetry of sources ⇒ mixture mass & total energy conserved as in (6).

### 1.6 Compact form — Eq. (12)

$$\partial_tq+\nabla\cdot f(q)+\sigma(q,\nabla q)=\psi_\mu(q)+\psi_\theta(q)+\psi_\nu(q) \tag{12a}$$
$$
q=\begin{pmatrix}\alpha_1\\ \alpha_1\rho_1\\ \alpha_2\rho_2\\ \rho\vec u\\ \alpha_1E_1\\ \alpha_2E_2\end{pmatrix},\quad
f(q)=\begin{pmatrix}0\\ \alpha_1\rho_1\vec u\\ \alpha_2\rho_2\vec u\\ \rho\vec u\otimes\vec u+(\alpha_1p_1+\alpha_2p_2)\mathbb I\\ \alpha_1(E_1+p_1)\vec u\\ \alpha_2(E_2+p_2)\vec u\end{pmatrix},\quad
\sigma(q,\nabla q)=\begin{pmatrix}\vec u\cdot\nabla\alpha_1\\0\\0\\0\\ \Sigma(q,\nabla q)\\ -\Sigma(q,\nabla q)\end{pmatrix} \tag{12b}
$$
$$
\psi_\mu=\begin{pmatrix}\mu(p_1-p_2)\\0\\0\\0\\-\mu p_I(p_1-p_2)\\ \mu p_I(p_1-p_2)\end{pmatrix},\quad
\psi_\theta=\begin{pmatrix}0\\0\\0\\0\\ \theta(T_2-T_1)\\ -\theta(T_2-T_1)\end{pmatrix},\quad
\psi_\nu=\begin{pmatrix}\nu\frac{g_2-g_1}{\rho_I}\\ \nu(g_2-g_1)\\ -\nu(g_2-g_1)\\0\\ \nu e_I(g_2-g_1)\\ -\nu e_I(g_2-g_1)\end{pmatrix}. \tag{12c}
$$

---

## 2. Numerical method with mechanical relaxation only (Section 3)

Model (13): $\partial_tq+\nabla\cdot f(q)+\sigma(q,\nabla q)=\psi_\mu(q)$. Fractional steps:

1. **Homogeneous hyperbolic system** over $\Delta t$:
   $\partial_tq+\nabla\cdot f(q)+\sigma(q,\nabla q)=0$. (14)
2. **Stiff mechanical relaxation**, $\mu\to\infty$: $\partial_tq=\psi_\mu(q)$ (15). Drives to $p_1=p_2=p$. Partial densities, mixture momentum, mixture total energy, mixture internal energy are constant during this step. $\alpha_1$, $p$, and $\alpha_k\mathcal E_k$ with $\mathcal E_k=\mathcal E_k(p,(\alpha_k\rho_k)/\alpha_k)$ are updated, then go back to step 1.

### 2.1 Mixture-energy consistency (Sec. 3.1, Definition 3.1)

Superscript 0: after step 1; superscript $*$: after relaxation. $\rho_k^*$... precisely $(\alpha_k\rho_k)^*=(\alpha_k\rho_k)^0$ [paper writes $\rho^*_k=\rho^0_k$ — as printed; I believe it means partial densities], $(\rho\vec u)^*=(\rho\vec u)^0$, $\mathcal E^*=\mathcal E^0$, $E^*=E^0$. $E^{0,C}$: discrete mixture total energy from a conservative approximation of (7).

**Definition 3.1.** Scheme is *mixture-energy-consistent* if
(i) Mixture total energy conservation consistency:
$$E^0=E^{0,C},\qquad E^0=(\alpha_1E_1)^0+(\alpha_2E_2)^0. \tag{16a}$$
(ii) Relaxed pressure consistency:
$$\mathcal E^{0,C}=\alpha_1^*\mathcal E_1\!\left(p^*,\frac{(\alpha_1\rho_1)^0}{\alpha_1^*}\right)+\alpha_2^*\mathcal E_2\!\left(p^*,\frac{(\alpha_2\rho_2)^0}{\alpha_2^*}\right),\qquad \mathcal E^{0,C}=E^{0,C}-\frac{(\rho\vec u)^0\cdot(\rho\vec u)^0}{2\rho^0}. \tag{16b}$$
For SG EOS this reads
$$p^*=\frac{\mathcal E^{0,C}-\big((\alpha_1\rho_1)^0\eta_1+(\alpha_2\rho_2)^0\eta_2\big)-\left(\frac{\alpha_1^*\gamma_1\pi_1}{\gamma_1-1}+\frac{\alpha_2^*\gamma_2\pi_2}{\gamma_2-1}\right)}{\frac{\alpha_1^*}{\gamma_1-1}+\frac{\alpha_2^*}{\gamma_2-1}}. \tag{17}$$

How to get it: in step 1 use a standard conservative scheme for the conservative part of each phasic energy eq., $\partial_t(\alpha_kE_k)+\nabla\cdot(\alpha_kE_k\vec u+\alpha_kp_k\vec u)$, and discretise $\Sigma$ **symmetrically** ($+\Sigma$ in eq. 5, $-\Sigma$ in eq. 6) so they cancel in the sum. Property (ii) then follows from the relaxation procedure of §3.3. With the internal-energy formulation (1) this is hard; refs [9,17] add an extra conservation law for $E$ and use $E^{0,C}$ to correct the thermodynamic state via the mixture EOS — but this does not in general guarantee (ii), i.e. $p^*=p(\mathcal E^{0,C},\alpha_1^*,\rho_1^*,\rho_2^*)$.

### 2.2 Wave propagation scheme, 1D (Sec. 3.2.1)

1D system $\partial_tq+\partial_xf(q)+\sigma(q,\partial_xq)=0$. Uniform grid $\Delta x$, $Q_i^n$ cell average.
$$Q_i^{n+1}=Q_i^n-\frac{\Delta t}{\Delta x}\big(\mathcal A^+\Delta Q_{i-1/2}+\mathcal A^-\Delta Q_{i+1/2}\big)-\frac{\Delta t}{\Delta x}\big(F^h_{i+1/2}-F^h_{i-1/2}\big) \tag{18}$$
Fluctuations $\mathcal A^\mp\Delta Q_{i+1/2}$ from local Riemann problems; $F^h$ second-order correction fluxes.
Riemann solution: $M$ waves $\mathcal W^l$, speeds $s^l$, $l=1..M$ ($M\le6$; HLLC $M=3$, Roe $M=6$). $f$-waves $\mathcal Z^l=s^l\mathcal W^l$.
$$\Delta q\equiv q_r-q_\ell=\sum_{l=1}^M\mathcal W^l \tag{19}$$
Conservation for conserved components $\xi=2,3,4$ ($\alpha_k\rho_k$, $\rho u$):
$$\Delta f^{(\xi)}\equiv f^{(\xi)}(q_r)-f^{(\xi)}(q_\ell)=\sum_l s^l\mathcal W^l_\xi=\sum_l\mathcal Z^l_\xi \tag{20}$$
Mixture total energy:
$$\Delta f_E\equiv f_E(q_r)-f_E(q_\ell)=\sum_l s^l(\mathcal W^l_5+\mathcal W^l_6)=\sum_l(\mathcal Z^l_5+\mathcal Z^l_6),\qquad f_E=u(E+\alpha_1p_1+\alpha_2p_2). \tag{21}$$
$$\mathcal A^\pm\Delta Q_{i+1/2}=\sum_{l=1}^M(s^l_{i+1/2})^\pm\mathcal W^l_{i+1/2}, \quad s^+=\max(s,0),\ s^-=\min(s,0) \tag{22}$$
$$F^h_{i+1/2}=\frac12\sum_{l=1}^M|s^l_{i+1/2}|\left(1-\frac{\Delta t}{\Delta x}|s^l_{i+1/2}|\right)\widetilde{\mathcal W}^l_{i+1/2} \tag{23}$$
($\widetilde{\mathcal W}$ = limited wave, written $\mathcal W^{h,l}$ in paper.)

### 2.3 HLLC-type solver (Sec. 3.2.2)

Key observation: one only needs $(\alpha_1E_1)^0+(\alpha_2E_2)^0=E^{0,C}$; individual phasic energies are **reset** at the start of each step 1 using relaxed $p^*$, $\alpha_1^*$ from step 2 (which by construction sum to the correct mixture energy). The individual phasic energy values from step 1 only serve as initial condition for the relaxation. Hence simplest approach: **neglect $\Sigma$** in the Riemann solver, use standard HLLC on the conservative portion + advection of $\alpha_1$.

Three waves, speeds
$$s^1=S_\ell,\quad s^2=S^\star,\quad s^3=S_r \tag{24}$$
separating $q_\ell,q_{\star\ell},q_{\star r},q_r$. Davis estimates:
$$S_\ell=\min(u_\ell-c_\ell,u_r-c_r),\qquad S_r=\max(u_\ell+c_\ell,u_r+c_r) \tag{25}$$
($c$ = frozen mixture sound speed (2)).
$$S^\star=\frac{p_r-p_\ell+\rho_\ell u_\ell(S_\ell-u_\ell)-\rho_ru_r(S_r-u_r)}{\rho_\ell(S_\ell-u_\ell)-\rho_r(S_r-u_r)} \tag{26}$$
($p=\alpha_1p_1+\alpha_2p_2$ presumably [uncertain: paper writes $p_r,p_\ell$; mixture pressure $\alpha_1p_1+\alpha_2p_2$ is natural since relaxed state has $p_1=p_2$].)
Rankine–Hugoniot for $\xi=2,\dots,6$:
$$f^{(\xi)}(q_r)-f^{(\xi)}(q_{\star r})=S_r(q_r-q_{\star r}),\ (27a)\quad f^{(\xi)}(q_{\star\ell})-f^{(\xi)}(q_\ell)=S_\ell(q_{\star\ell}-q_\ell),\ (27b)\quad f^{(\xi)}(q_{\star r})-f^{(\xi)}(q_{\star\ell})=S^\star(q_{\star r}-q_{\star\ell}).\ (27c)$$
Middle states, $\iota=\ell,r$:
$$
q_{\star\iota}=\begin{pmatrix}
\alpha_{1,\iota}\\
(\alpha_1\rho_1)_\iota\frac{S_\iota-u_\iota}{S_\iota-S^\star}\\
(\alpha_2\rho_2)_\iota\frac{S_\iota-u_\iota}{S_\iota-S^\star}\\
\rho_\iota\frac{S_\iota-u_\iota}{S_\iota-S^\star}S^\star\\
(\alpha_1\rho_1)_\iota\frac{S_\iota-u_\iota}{S_\iota-S^\star}\left(\frac{E_{1,\iota}}{\rho_{1,\iota}}+(S^\star-u_\iota)\left(S^\star+\frac{p_\iota}{\rho_{1,\iota}(S_\iota-u_\iota)}\right)\right)\\
(\alpha_2\rho_2)_\iota\frac{S_\iota-u_\iota}{S_\iota-S^\star}\left(\frac{E_{2,\iota}}{\rho_{2,\iota}}+(S^\star-u_\iota)\left(S^\star+\frac{p_\iota}{\rho_{2,\iota}(S_\iota-u_\iota)}\right)\right)
\end{pmatrix} \tag{28}
$$
[uncertain: in the energy rows the paper prints $p_\iota$ without phase index; it may be $p_{1,\iota}$ / $p_{2,\iota}$ — but the RH condition (27) with flux $\alpha_k(E_k+p_k)u$ gives $p_{k,\iota}$; check against corrigendum. With relaxed initial data $p_1=p_2=p$ anyway.]
Waves:
$$\mathcal W^1=q_{\star\ell}-q_\ell,\quad\mathcal W^2=q_{\star r}-q_{\star\ell},\quad\mathcal W^3=q_r-q_{\star r}. \tag{29}$$
(27) for $\xi=2,3,4$ ⇒ (20) with $M=3$. For $\xi=5,6$:
$$\Delta f^{(\xi)}=\sum_{l=1}^3s^l\mathcal W^l_\xi=\sum_l\mathcal Z^l_\xi,\quad \Delta f^{(4+k)}=u\alpha_k(E_k+p_k) \tag{30}$$
Since $f_E=f^{(5)}+f^{(6)}$, (21) holds.

More sophisticated variants: add a non-conservative correction $\Phi^l$ to the f-wave energy components: $\mathcal Z^{\Sigma,l}_{4+k}=\mathcal Z^l_{4+k}+(-1)^k\Phi^l$, $k=1,2$ — any $\Phi^l$ preserves (21) by symmetry. Numerically no remarkable difference was observed, so **the simplest method neglecting $\pm\Sigma$ is adopted**.

### 2.4 Roe-type solver (Sec. 3.2.3)

Roe matrix $\hat A(q_\ell,q_r)\in\mathbb R^{6\times6}$, averaged version of $A(q)$ in quasi-linear form $\partial_tq+\partial_xf+\sigma=\partial_tq+A(q)\partial_xq=0$ — here $\pm\Sigma$ **is** included (implicitly defines meaning of non-conservative terms). Required:
$$\sum_{l=1}^6\hat A_{\xi l}(q_r-q_\ell)^{(l)}=f^{(\xi)}(q_r)-f^{(\xi)}(q_\ell),\ \xi=2,3,4 \tag{31a}$$
$$\sum_{l=1}^6(\hat A_{5l}+\hat A_{6l})(q_r-q_\ell)^{(l)}=f_E(q_r)-f_E(q_\ell) \tag{31b}$$
Waves $\mathcal W^\xi=\hat\beta_\xi\hat r_\xi$, $s^\xi=\hat\lambda_\xi$, $\xi=1..6$ (32), $\hat\beta$ = projection coefficients of $q_r-q_\ell$ on Roe eigenvectors. Roe matrix and eigenstructure for SG in Appendix B. HLLC more robust; preferred for 1D and normal waves in 2D; Roe eigenstructure used in 2D to compute transverse fluctuations.

### 2.5 2D wave propagation (Sec. 3.2.4)

$$
Q^{n+1}_{ij}=Q^n_{ij}-\frac{\Delta t}{\Delta x}\big(\mathcal A^+\Delta Q_{i-1/2,j}+\mathcal A^-\Delta Q_{i+1/2,j}\big)-\frac{\Delta t}{\Delta y}\big(\mathcal B^+\Delta Q_{i,j-1/2}+\mathcal B^-\Delta Q_{i,j+1/2}\big)
-\frac{\Delta t}{\Delta x}(F^h_{i+1/2,j}-F^h_{i-1/2,j})-\frac{\Delta t}{\Delta y}(G^h_{i,j+1/2}-G^h_{i,j-1/2}) \tag{33}
$$
Normal fluctuations by HLLC; transverse fluctuations by projecting on Roe eigenvectors of orthogonal direction (App. B). Generalises to logically rectangular curvilinear grids.

### 2.6 Stiff mechanical relaxation step (Sec. 3.3)

ODE system (limit $\mu\to\infty$):
$$
\partial_t\alpha_1=\mu(p_1-p_2),\ (34a)\quad \partial_t(\alpha_k\rho_k)=0,\ (34b,c)\quad \partial_t(\rho\vec u)=0,\ (34d)\quad
\partial_t(\alpha_1E_1)=-\mu p_I(p_1-p_2),\ (34e)\quad \partial_t(\alpha_2E_2)=\mu p_I(p_1-p_2).\ (34f)
$$
Hence
$$(\alpha_k\rho_k)^*=(\alpha_k\rho_k)^0,\ k=1,2,\qquad (\rho\vec u)^*=(\rho\vec u)^0 \tag{35}$$
so $\rho^*=\rho^0$, $\vec u^*=\vec u^0$. Combining (34a),(34e),(34f),(35):
$$\partial_t(\alpha_1E_1)=\partial_t(\alpha_1\mathcal E_1)=-p_I\partial_t\alpha_1,\qquad \partial_t(\alpha_2E_2)=\partial_t(\alpha_2\mathcal E_2)=p_I\partial_t\alpha_1. \tag{36}$$
Sum: $\partial_tE=\partial_t\mathcal E=0$ ⇒ $E^*=E^0$, $\mathcal E^*=\mathcal E^0$.

**Approximation:** linear variation of interface pressure with $\alpha_1$:
$$p_I=p_I^0+\frac{p_I^*-p_I^0}{\alpha_1^*-\alpha_1^0}(\alpha_1-\alpha_1^0). \tag{37}$$
Integrating (36):
$$
\begin{aligned}
(\alpha_1E_1)^*-(\alpha_1E_1)^0&=(\alpha_1\mathcal E_1)^*-(\alpha_1\mathcal E_1)^0=-\frac{p_I^0+p_I^*}{2}(\alpha_1^*-\alpha_1^0), &(38a)\\
(\alpha_2E_2)^*-(\alpha_2E_2)^0&=(\alpha_2\mathcal E_2)^*-(\alpha_2\mathcal E_2)^0=\frac{p_I^0+p_I^*}{2}(\alpha_1^*-\alpha_1^0). &(38b)
\end{aligned}
$$
(37) ⇔ approximating $p_I\approx\bar p_I=\frac{p_I^0+p_I^*}{2}$ in (36) (also in [17]).
Impose $p_1^*=p_2^*=p_I^*=p^*$, i.e. $\mathcal E_k^*=\mathcal E_k(p^*,(\alpha_k\rho_k)^0/\alpha_k^*)$. With (35) and $\alpha_1+\alpha_2=1$, (38) gives two equations for the two unknowns $\alpha_1^*$, $p^*$. For SG EOS: a **quadratic for $p^*$** (always has a physically admissible root), then $\alpha_1^*$ explicitly — Appendix C.1. $p_I^0$ computed from (1)-closure $p_I=(Z_2p_1+Z_1p_2)/(Z_1+Z_2)$ with step-1 values [the paper uses $p_I^0$ from the state after step 1; formula of $p_I$ given in §2.1].

**Why consistent:** $p^*,\alpha_1^*$ satisfy the mixture internal energy equation (4) by construction (sum of (38a,b) gives $\mathcal E^*=\mathcal E^0$), hence $p^*=p(\mathcal E^0,\alpha_1^*,(\alpha_1\rho_1)^0/\alpha_1^*,(\alpha_2\rho_2)^0/\alpha_2^*)$. Since the wave propagation scheme gives $E^0=(\alpha_1E_1)^0+(\alpha_2E_2)^0=E^{0,C}$, $\mathcal E^0=E^0-\frac{(\rho u)^0\cdot(\rho u)^0}{2\rho^0}$ agrees with mixture total energy conservation ⇒ properties (i),(ii) hold.

After relaxation the phasic energies are reset: $(\alpha_kE_k)^*=\alpha_k^*\mathcal E_k(p^*,\rho_k^*)+\tfrac12(\alpha_k\rho_k)^0|\vec u^0|^2$ [this reset is implied by the text in §3.2.2: "the individual phasic energies are re-set at the beginning of Step 1 by using the relaxed variables $p^*$ and $\alpha_1^*$"].

---

## 3. Thermal and chemical relaxation (Section 4)

Assumption (as in [17]): mechanical relaxation time $1/\mu\ll 1/\theta,1/\nu$; thermal and chemical relaxation occur under pressure equilibrium. Heat/mass transfer handled **after** the mechanical relaxation step (15). **Step 3**, activated at liquid–vapour interfaces only, two sub-steps:

- **3a. Thermal relaxation**: solve $\partial_tq=\psi_\mu(q)+\psi_\theta(q)$ (39) with $\mu,\theta\to\infty$. Relaxes to $T_1=T_2=T$ under mechanical equilibrium. The updated $T,p$ are used to check the metastability condition $T_{liq}=T>T_{sat}(p)$. If metastable, go to 3b; otherwise update volume fraction, pressure, phasic energies and return to Step 1.
- **3b. Thermo-chemical relaxation**: if metastable, solve $\partial_tq=\psi_\mu+\psi_\theta+\psi_\nu$ (40) with $\mu,\theta,\nu\to\infty$: $T_1=T_2=T$, $g_1=g_2$, $p_1=p_2$. Mixture density, momentum, total and internal energy constant; partial densities vary. Outputs $\alpha_1$, $p$, $\alpha_k\rho_k$, $\alpha_k\mathcal E_k$ with $\mathcal E_k=\mathcal E_k(p,(\alpha_k\rho_k)/\alpha_k)$.

Nomenclature in §5: "$p$–$pT$–$pTg$" = Steps 2, 3a, 3b.

### 3.1 Thermal relaxation step (Sec. 4.1)

ODEs (41a–f): $\partial_t\alpha_1=\mu(p_1-p_2)$; $\partial_t(\alpha_k\rho_k)=0$; $\partial_t(\rho\vec u)=0$; $\partial_t(\alpha_1E_1)=-\mu p_I(p_1-p_2)+\theta(T_2-T_1)$; $\partial_t(\alpha_2E_2)=\mu p_I(p_1-p_2)-\theta(T_2-T_1)$.
Initial state: superscript $*$ (after pressure relaxation); final: $**$.
$$(\alpha_k\rho_k)^{**}=(\alpha_k\rho_k)^0,\quad(\rho\vec u)^{**}=(\rho\vec u)^0,\quad\mathcal E^{**}=\mathcal E^0,\quad E^{**}=E^0 \tag{42}$$
Equilibrium conditions:
$$p_1^{**}=p_2^{**}\equiv p^{**},\qquad T_1^{**}=T_2^{**}\equiv T^{**}. \tag{43}$$
(42)+(43) → algebraic system for $\alpha_1^{**},\rho_k^{**},p^{**},T^{**}$; for SG EOS reduces to a **quadratic in $p^{**}$** (Appendix C.2).

### 3.2 Thermo-chemical relaxation step (Sec. 4.2)

Activated if metastable: $T^{**}>T_{sat}(p^{**})$, $T_{sat}$ computed by solving the nonlinear $p$–$T$ saturation equation (11) by **Newton's method**. ODEs, $\mu,\theta,\nu\to\infty$:
$$
\begin{aligned}
&\partial_t\alpha_1=\mu(p_1-p_2)+\nu\frac{g_2-g_1}{\rho_I}, &(44a)\\
&\partial_t(\alpha_1\rho_1)=\nu(g_2-g_1),\quad \partial_t(\alpha_2\rho_2)=-\nu(g_2-g_1), &(44b,c)\\
&\partial_t(\rho\vec u)=0, &(44d)\\
&\partial_t(\alpha_1E_1)=-\mu p_I(p_1-p_2)+\theta(T_2-T_1)+\nu e_I(g_2-g_1), &(44e)\\
&\partial_t(\alpha_2E_2)=\mu p_I(p_1-p_2)-\theta(T_2-T_1)-\nu e_I(g_2-g_1). &(44f)
\end{aligned}
$$
Initial: $**$; final: $\circledast$ (full equilibrium).
$$\rho^\circledast=\rho^0,\quad(\rho\vec u)^\circledast=(\rho\vec u)^0,\quad\mathcal E^\circledast=\mathcal E^0,\quad E^\circledast=E^0. \tag{45}$$
Mass & energy conservation:
$$\rho^0=\alpha_1^\circledast\rho_1^\circledast+\alpha_2^\circledast\rho_2^\circledast,\ (46a)\qquad \mathcal E^0=\alpha_1^\circledast\mathcal E_1^\circledast+\alpha_2^\circledast\mathcal E_2^\circledast.\ (46b)$$
Equilibrium:
$$p_1^\circledast=p_2^\circledast\equiv p^\circledast,\ (47a)\quad T_1^\circledast=T_2^\circledast\equiv T^\circledast,\ (47b)\quad g_1^\circledast=g_2^\circledast.\ (47c)$$
(46),(47), $\alpha_1+\alpha_2=1$ and $\mathcal E_k=\mathcal E_k(p_k,\rho_k)$ → algebraic system for $\alpha_1^\circledast,p^\circledast,T^\circledast,\rho_1^\circledast,\rho_2^\circledast$; reduced to a single nonlinear equation in $p^\circledast$, solved by Newton (Appendix C.3). **No need to specify $\rho_I$ or $e_I$.** Same approach used in [23]. Consistency: relaxed thermodynamic state satisfies mixture energy (4) and $E^\circledast=E^{**}=E^*=E^0$ conservation-consistent.

<!-- CONTINUED BELOW -->
