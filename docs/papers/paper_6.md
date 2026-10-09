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

*Check (my derivation, not in paper):* with $P=\alpha_1p_1+\alpha_2p_2$, phasic kinetic energy obeys $\partial_t(\tfrac12\alpha_k\rho_ku^2)+\nabla\cdot(\tfrac12\alpha_k\rho_ku^2\vec u)=-Y_k\vec u\cdot\nabla P$. Adding to (1e): $\partial_t(\alpha_1E_1)+\nabla\cdot(\alpha_1E_1\vec u+\alpha_1p_1\vec u)=\vec u\cdot\nabla(\alpha_1p_1)-Y_1\vec u\cdot\nabla P=\vec u\cdot(Y_2\nabla(\alpha_1p_1)-Y_1\nabla(\alpha_2p_2))=-\Sigma$. Consistent with (6g).

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
$\epsilon_I$ identifies interface locations (e.g. $\epsilon_I=10^{-4}$; see [8] for choice of tolerances). Sign convention: $\dot m=\nu(g_2-g_1)>0$ transfers mass from phase 2 to phase 1 (mass flows toward the phase of lower Gibbs energy). The paper does not fix which phase is liquid in Sec. 2.3; in the tests phases are labelled $l$/$v$. Since the relaxation is instantaneous and solved algebraically, the labelling is immaterial for the algorithm.

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
Impose $p_1^*=p_2^*=p_I^*=p^*$, i.e. $\mathcal E_k^*=\mathcal E_k(p^*,(\alpha_k\rho_k)^0/\alpha_k^*)$. With (35) and $\alpha_1+\alpha_2=1$, (38) gives two equations for the two unknowns $\alpha_1^*$, $p^*$. For SG EOS: a **quadratic for $p^*$** (always has a physically admissible root), then $\alpha_1^*$ explicitly — Appendix C.1. $p_I^0$ = interface pressure $p_I=(Z_2p_1+Z_1p_2)/(Z_1+Z_2)$ evaluated with the step-1 (superscript 0) state — confirmed by re-deriving (C.2), see §4 below.

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

---

## 4. Appendix C — closed-form relaxed states for SG EOS (implementation core)

### C.1 Mechanical relaxation (Step 2)

Quadratic for $p^*$:
$$a(p^*)^2+bp^*+d=0,$$
$$
a=1+\gamma_2\alpha_1^0+\gamma_1\alpha_2^0,\qquad
b=C_1\alpha_2^0+C_2\alpha_1^0-(1+\gamma_2)\alpha_1^0p_1^0-(1+\gamma_1)\alpha_2^0p_2^0,\qquad
d=-\big(C_2\alpha_1^0p_1^0+C_1\alpha_2^0p_2^0\big),
$$
$$C_1=2\gamma_1\pi_1+(\gamma_1-1)p_I^0,\qquad C_2=2\gamma_2\pi_2+(\gamma_2-1)p_I^0.$$
(The subscript in the printed $C_k$ looks like "$p_1^0$"/"$p_I^0$"; I re-derived C.2 from (38a) with SG EOS and it requires **$p_I^0$**, the interface pressure $p_I=(Z_2p_1+Z_1p_2)/(Z_1+Z_2)$ evaluated with the post-hyperbolic-step state. Derivation check: (38a) with $\alpha_k\mathcal E_k=\alpha_k\frac{p_k+\gamma_k\pi_k}{\gamma_k-1}+\alpha_k\rho_k\eta_k$ gives exactly (C.2) with $C_1$ containing $p_I^0$, and $\sum\alpha_k^*=1$ gives the leading coefficient $-2(1+\gamma_2\alpha_1^0+\gamma_1\alpha_2^0)$, matching $a$.)
$$p^*=\frac{-b+\sqrt{b^2-4ad}}{2a} \tag{C.1}$$
Root retained = physically meaningful positive pressure; note $C_1,C_2>0$, $a>0$, $d<0$, $b^2-4ad>0$. [Positivity of $C_k$, $d<0$ requires $p_k^0,p_I^0$ such that these hold — e.g. $p_k^0>0$ or at least $p_I^0>-2\gamma_k\pi_k/(\gamma_k-1)$; with negative pressures in liquid (tension) check sign of $d$. My remark.]
$$\alpha_1^*=\frac{(\gamma_1-1)p^*+2p_1^0+C_1}{(\gamma_1+1)p^*+C_1}\,\alpha_1^0. \tag{C.2}$$
Then $\rho_k^*=(\alpha_k\rho_k)^0/\alpha_k^*$, $\alpha_k\mathcal E_k^*=\alpha_k^*\mathcal E_k(p^*,\rho_k^*)$. [Alternatively, $p^*$ can be recomputed from the mixture EOS (17) with $\alpha_1^*$ — identical analytically.]

### C.2 Thermal relaxation (Step 3a)

From (42),(43):
$$a'(p^{**})^2+b'p^{**}+d'=0 \tag{C.3}$$
$$a'=C_{v1}(\alpha_1\rho_1)^0+C_{v2}(\alpha_2\rho_2)^0,$$
$$
\begin{aligned}
b'=\;&\eta_1C_{v1}(\gamma_1-1)\big((\alpha_1\rho_1)^0\big)^2+\eta_2C_{v2}(\gamma_2-1)\big((\alpha_2\rho_2)^0\big)^2+(\alpha_1\rho_1)^0C_{v1}(\gamma_1\pi_1+\pi_2)\\
&+(\alpha_2\rho_2)^0C_{v2}(\gamma_2\pi_2+\pi_1)+(\alpha_1\rho_1)^0(\alpha_2\rho_2)^0\big(\eta_1C_{v2}(\gamma_2-1)+\eta_2C_{v1}(\gamma_1-1)\big)\\
&-\mathcal E^0\big(C_{v1}(\gamma_1-1)(\alpha_1\rho_1)^0+C_{v2}(\gamma_2-1)(\alpha_2\rho_2)^0\big),
\end{aligned}
$$
$$
\begin{aligned}
d'=\;&\eta_1C_{v1}(\gamma_1-1)\pi_2\big((\alpha_1\rho_1)^0\big)^2+\eta_2C_{v2}(\gamma_2-1)\pi_1\big((\alpha_2\rho_2)^0\big)^2+\big((\alpha_1\rho_1)^0C_{v1}\gamma_1\\
&+(\alpha_2\rho_2)^0C_{v2}\gamma_2\big)\pi_2\pi_1+(\alpha_1\rho_1)^0(\alpha_2\rho_2)^0\big(\eta_1C_{v2}(\gamma_2-1)\pi_1+\eta_2C_{v1}(\gamma_1-1)\pi_2\big)\\
&-\mathcal E^0\big(C_{v1}(\gamma_1-1)\pi_2(\alpha_1\rho_1)^0+C_{v2}(\gamma_2-1)\pi_1(\alpha_2\rho_2)^0\big).
\end{aligned}
$$
[Transcribed as printed. I am somewhat unsure about the exact pairing in the cross terms of $b'$/$d'$ — **[uncertain]**; recommended to re-derive symbolically: unknowns $\alpha_1,p,T$ with $\rho_k=(p+\pi_k)/((\gamma_k-1)C_{vk}T)$, $\alpha_k=(\alpha_k\rho_k)^0/\rho_k$, $\alpha_1+\alpha_2=1$, $\mathcal E^0=\sum\alpha_k\frac{p+\gamma_k\pi_k}{\gamma_k-1}+\sum(\alpha_k\rho_k)^0\eta_k$. Eliminating $T$ via $\alpha_1+\alpha_2=1$: $T=\frac{(\gamma_1-1)C_{v1}(\alpha_1\rho_1)^0}{p+\pi_1}+\frac{(\gamma_2-1)C_{v2}(\alpha_2\rho_2)^0}{p+\pi_2}$, then substitution into the energy equation gives a quadratic in $p$.]
$$\alpha_1^{**}=\frac{C_{v1}(\gamma_1-1)(p^{**}+\pi_2)(\alpha_1\rho_1)^0}{C_{v1}(\gamma_1-1)(p^{**}+\pi_2)(\alpha_1\rho_1)^0+C_{v2}(\gamma_2-1)(p^{**}+\pi_1)(\alpha_2\rho_2)^0} \tag{C.4}$$
$$T^{**}=\frac{(p^{**}+\pi_1)\alpha_1^{**}}{(\gamma_1-1)C_{v1}(\alpha_1\rho_1)^0}=\frac{(p^{**}+\pi_2)(1-\alpha_1^{**})}{(\gamma_2-1)C_{v2}(\alpha_2\rho_2)^0} \tag{C.5}$$
(Root choice for C.3: not stated explicitly in the paper; presumably the positive/physical root, analogous to C.1 [uncertain].)

### C.3 Thermo-chemical relaxation (Step 3b)

Using (46),(47) and $\mathcal E_k=\mathcal E_k(p_k,\rho_k)$: first a quadratic for $T^\circledast$ as a function of $p^\circledast$:
$$a_p(p^\circledast)(T^\circledast)^2+b_p(p^\circledast)T^\circledast+d_p(p^\circledast)=0, \tag{C.6}$$
$$
\begin{aligned}
a_p(p)&=\rho^0C_{v1}C_{v2}\big((\gamma_2-1)(p+\gamma_1\pi_1)-(\gamma_1-1)(p+\gamma_2\pi_2)\big),\\
b_p(p)&=\mathcal E^0\big((\gamma_1-1)C_{v1}(p+\pi_2)-(\gamma_2-1)C_{v2}(p+\pi_1)\big)\\
&\quad+\rho^0\big((\gamma_2-1)C_{v2}\eta_1(p+\pi_1)-(\gamma_1-1)C_{v1}\eta_2(p+\pi_2)\big)\\
&\quad+C_{v2}(p+\pi_1)(p+\gamma_2\pi_2)-C_{v1}(p+\pi_2)(p+\gamma_1\pi_1),\\
d_p(p)&=(\eta_2-\eta_1)(p+\pi_1)(p+\pi_2).
\end{aligned}
$$
$$T^\circledast(p^\circledast)=\frac{-b_p+\sqrt{b_p^2-4a_pd_p}}{2a_p} \tag{C.7}$$
Substitute into the Gibbs equilibrium (47c)/(11):
$$A_S+\frac{B_S}{T^\circledast(p^\circledast)}+C_S\log T^\circledast(p^\circledast)+D_S\log(p^\circledast+\pi_1)-\log(p^\circledast+\pi_2)=0 \tag{C.8}$$
— single scalar nonlinear equation in $p^\circledast$, solved by Newton (Sec. 4.2). Then $T^\circledast$ from (C.7) and
$$\rho_k^\circledast=\frac{p^\circledast+\pi_k}{(\gamma_k-1)C_{vk}T^\circledast},\ k=1,2,\qquad \alpha_1^\circledast=\frac{\rho^*-\rho_2^\circledast}{\rho_1^\circledast-\rho_2^\circledast}. \tag{C.9}$$
($\rho^*=\rho^0$.) Then $(\alpha_k\rho_k)^\circledast=\alpha_k^\circledast\rho_k^\circledast$ and phasic energies from EOS.

(Derivation sketch, mine: from (46a) $\alpha_1=(\rho-\rho_2)/(\rho_1-\rho_2)$; (46b) with $\mathcal E_k=\frac{p+\gamma_k\pi_k}{\gamma_k-1}+\rho_k\eta_k$ and $\rho_k=\frac{p+\pi_k}{(\gamma_k-1)C_{vk}T}$, multiply through by $T$ → quadratic in $T$.)

**Abnormal case (end of Appendix C, p.356):** the $pTg$ equilibrium solution of (46)–(47) **may not be physically admissible** (e.g. $\alpha_1^\circledast\notin[0,1]$). Then assume the mixture is composed nearly of the species $k$ having the highest entropy $s_k^{**}$ (idea from [14], Faccanoni et al.). Fix $\alpha_k^\circledast=\bar\alpha_k^\circledast=1-\epsilon$ (e.g. $\epsilon=10^{-8}$), and solve for $p^\circledast,T^\circledast,\rho_1^\circledast,\rho_2^\circledast$ from conservation (46) and $p,T$ equilibrium (47a,b) only (no Gibbs equality). For SG EOS this reduces to a single quadratic in $p^\circledast$ (or $T^\circledast$); select the physically admissible root that **maximises the total mixture entropy** $s^\circledast=Y_1^\circledast s_1^\circledast+Y_2^\circledast s_2^\circledast$.

---

## 5. Appendix A — eigenstructure (1D homogeneous system)

$$\lambda_1=u-c,\quad\lambda_{2..5}=u,\quad\lambda_6=u+c, \tag{A.1}$$
$c$ = mixture sound speed (2). Right eigenvectors (columns), with $H_k=h_k+u^2/2$, $\Pi_k=-\rho_kc_k^2+p_k(1+\kappa_k)$, $\kappa_k=\partial p_k/\partial\mathcal E_k$; component order $(\alpha_1,\alpha_1\rho_1,\alpha_2\rho_2,\rho u,\alpha_1E_1,\alpha_2E_2)$:
$$
R=\begin{pmatrix}
0&0&0&0&1&0\\
Y_1&0&0&1&0&Y_1\\
Y_2&0&1&0&0&Y_2\\
u-c&0&u&u&0&u+c\\
Y_1(H_1-uc)&-\frac{\kappa_2}{\kappa_1}&\frac{\kappa_2}{\kappa_1}H_2-\frac{c_2^2}{\kappa_1}&H_1-\frac{c_1^2}{\kappa_1}&\frac{\Pi_1-\Pi_2}{\kappa_1}&Y_1(H_1+uc)\\
Y_2(H_2-uc)&1&0&0&0&Y_2(H_2+uc)
\end{pmatrix} \tag{A.2}
$$
(For SG: $\Pi_k=\gamma_k\pi_k$ as stated in App. B — check: $-\rho c^2+\gamma p=-\gamma(p+\pi)+\gamma p=-\gamma\pi$; so the sign convention gives $\Pi_k=-\gamma_k\pi_k$ from the A definition but App. B says "$\Pi_k=\gamma_k\pi_k$" [uncertain sign; matters only through $\Pi_1-\Pi_2$ and $-\Pi_1+\Pi_2$ in $\hat A$ — re-derive before use].)

## 6. Appendix B — Roe matrix (SG EOS)

$\hat A=A(\hat u,\hat Y_1,\widetilde{uY_1},\widehat{Y_1H_1},\widehat{Y_2H_2})$:
$$
\hat A=\begin{pmatrix}
\hat u&0&0&0&0&0\\
0&\hat u\hat Y_2&-\hat u\hat Y_1&\hat Y_1&0&0\\
0&-\hat u\hat Y_2&\hat u\hat Y_1&\hat Y_2&0&0\\
-\Pi_1+\Pi_2&\hat A_{42}&\hat A_{43}&2\hat u-(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2})&\kappa_1&\kappa_2\\
(-\Pi_1+\Pi_2)\widetilde{uY_1}&\hat A_{52}&\hat A_{53}&\widehat{Y_1H_1}-\widetilde{uY_1}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2})&\hat u+\kappa_1\widetilde{uY_1}&\widetilde{uY_1}\kappa_2\\
(-\Pi_1+\Pi_2)\widetilde{uY_2}&\hat A_{62}&\hat A_{63}&\widehat{Y_2H_2}-\widetilde{uY_2}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2})&\widetilde{uY_2}\kappa_1&\hat u+\kappa_2\widetilde{uY_2}
\end{pmatrix} \tag{B.1a}
$$
(SG: $\Pi_k=\gamma_k\pi_k$, $\kappa_k=\gamma_k-1$, $\chi_k=-(\gamma_k-1)\eta_k$.)
$$
\begin{aligned}
\hat A_{42}&=-\hat u^2+\chi_1-\kappa_1\tfrac{\hat u^2}{2}+\hat u(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}), &(B.1b)\\
\hat A_{43}&=-\hat u^2+\chi_2-\kappa_2\tfrac{\hat u^2}{2}+\hat u(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}), &(B.1c)\\
\hat A_{52}&=-\widehat{Y_1H_1}\hat u+\widetilde{uY_1}\big(\chi_1-\kappa_1\tfrac{\hat u^2}{2}\big)+\hat u\widetilde{uY_1}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}), &(B.1d)\\
\hat A_{53}&=-\widehat{Y_1H_1}\hat u+\widetilde{uY_1}\big(\chi_2-\kappa_2\tfrac{\hat u^2}{2}\big)+\hat u\widetilde{uY_1}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}), &(B.1e)\\
\hat A_{62}&=-\widehat{Y_2H_2}\hat u+\widetilde{uY_2}\big(\chi_1-\kappa_1\tfrac{\hat u^2}{2}\big)+\hat u\widetilde{uY_2}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}), &(B.1f)\\
\hat A_{63}&=-\widehat{Y_2H_2}\hat u+\widetilde{uY_2}\big(\chi_2-\kappa_2\tfrac{\hat u^2}{2}\big)+\hat u\widetilde{uY_2}(\kappa_1\widetilde{uY_1}+\kappa_2\widetilde{uY_2}). &(B.1g)
\end{aligned}
$$
Averages (Roe, $\sqrt\rho$ weights):
$$\hat u=\frac{u_\ell\sqrt{\rho_\ell}+u_r\sqrt{\rho_r}}{\sqrt{\rho_\ell}+\sqrt{\rho_r}},\ (B.2a)\quad \hat Y_k=\frac{Y_{k\ell}\sqrt{\rho_\ell}+Y_{kr}\sqrt{\rho_r}}{\sqrt{\rho_\ell}+\sqrt{\rho_r}},\ (B.2b)\quad \overline{uY_k}=\frac{(uY_k)_\ell\sqrt{\rho_\ell}+(uY_k)_r\sqrt{\rho_r}}{\sqrt{\rho_\ell}+\sqrt{\rho_r}},\ (B.2c)$$
$$\widetilde{uY_k}=\tfrac12(\hat u\hat Y_k+\overline{uY_k}),\ (B.2d)\qquad \widehat{Y_kH_k}=\frac{(Y_kH_k)_\ell\sqrt{\rho_\ell}+(Y_kH_k)_r\sqrt{\rho_r}}{\sqrt{\rho_\ell}+\sqrt{\rho_r}}.\ (B.2e)$$
(In print, (B.2c) defines a $\sqrt\rho$-average of $uY_k$ — notation distinguishes it from the $\widetilde{uY_k}$ of (B.2d) [notation partially uncertain].) $\hat Y_1+\hat Y_2=1$, $\widetilde{uY_1}+\widetilde{uY_2}=\hat u$.
$$\hat c=\sqrt{\widetilde{Y_1c_1^2}+\widetilde{Y_2c_2^2}},\qquad \widetilde{Y_kc_k^2}=\kappa_k\big(\widehat{Y_kH_k}-\tfrac{\hat u^2}{2}\hat Y_k\big)+\chi_k\hat Y_k. \tag{B.2f}$$
Eigenvalues $\hat\lambda_1=\hat u-\hat c$, $\hat\lambda_{2..5}=\hat u$, $\hat\lambda_6=\hat u+\hat c$ (B.3). Eigenvectors:
$$
\hat R=\begin{pmatrix}
0&0&0&0&1&0\\
\hat Y_1&0&0&1&0&\hat Y_1\\
\hat Y_2&0&1&0&0&\hat Y_2\\
\hat u-\hat c&0&\hat u&\hat u&0&\hat u+\hat c\\
\widehat{Y_1H_1}-\widetilde{uY_1}\hat c&-\frac{\kappa_2}{\kappa_1}&-\frac{\chi_2}{\kappa_1}+\frac{\kappa_2}{\kappa_1}\frac{\hat u^2}{2}&-\frac{\chi_1}{\kappa_1}+\frac{\hat u^2}{2}&\frac{\Pi_1-\Pi_2}{\kappa_1}&\widehat{Y_1H_1}+\widetilde{uY_1}\hat c\\
\widehat{Y_2H_2}-\widetilde{uY_2}\hat c&1&0&0&0&\widehat{Y_2H_2}+\widetilde{uY_2}\hat c
\end{pmatrix} \tag{B.4}
$$
Wave strengths ($\Delta(\cdot)=(\cdot)_r-(\cdot)_\ell$):
$$
\begin{aligned}
\hat\beta_{1,6}&=\frac{\Delta(\alpha_1p_1+\alpha_2p_2)\mp\hat c\sqrt{\rho_\ell\rho_r}\,\Delta u}{2\hat c^2}, &(B.5a)\\
\hat\beta_2&=-\frac{\Delta(\alpha_1p_1+\alpha_2p_2)}{\hat c^2}\widehat{Y_2H_2}+\frac{\hat u^2}{2}\Delta(\alpha_2\rho_2)+\Delta(\alpha_2\mathcal E_2), &(B.5b)\\
\hat\beta_3&=\Delta(\alpha_2\rho_2)-\hat Y_2\frac{\Delta(\alpha_1p_1+\alpha_2p_2)}{\hat c^2}, &(B.5c)\\
\hat\beta_4&=\Delta(\alpha_1\rho_1)-\hat Y_1\frac{\Delta(\alpha_1p_1+\alpha_2p_2)}{\hat c^2}, &(B.5d)\\
\hat\beta_5&=\Delta\alpha_1. &(B.5e)
\end{aligned}
$$
[(B.5b): printed as "$\Delta(\alpha_2\mathcal E_2)$" — might be $\Delta(\alpha_2E_2)$; transcribed as seen, **[uncertain]**.]

**2D (B.6–B.8):** plane-wave Riemann problem in $x$ with $\vec u=(u,v)$, 7 unknowns $A^x\in\mathbb R^{7\times7}$; kinetic energy terms use $(u^2+v^2)/2$. Eigenvalues $\hat u-\hat c$, $\hat u$ (×5), $\hat u+\hat c$.
$$
\hat R^x=\begin{pmatrix}
0&0&0&0&1&0&0\\
\hat Y_1&0&0&1&0&0&\hat Y_1\\
\hat Y_2&0&1&0&0&0&\hat Y_2\\
\hat u-\hat c&0&\hat u&\hat u&0&0&\hat u+\hat c\\
\hat v&0&0&0&0&1&\hat v\\
\widehat{Y_1H_1}-\widetilde{uY_1}\hat c&-\frac{\kappa_2}{\kappa_1}&-\frac{\chi_2}{\kappa_1}+\frac{\kappa_2}{\kappa_1}\hat{\mathcal K}-\hat v\frac{\hat V}{\kappa_1}&-\frac{\chi_1}{\kappa_1}+\hat{\mathcal K}-\hat v\frac{\hat V}{\kappa_1}&\frac{\Pi_1-\Pi_2}{\kappa_1}&\frac{\hat V}{\kappa_1}&\widehat{Y_1H_1}+\widetilde{uY_1}\hat c\\
\widehat{Y_2H_2}-\widetilde{uY_2}\hat c&1&0&0&0&0&\widehat{Y_2H_2}+\widetilde{uY_2}\hat c
\end{pmatrix} \tag{B.7}
$$
$\hat v$ Roe average (B.8a); $\overline{vY_k}$, $\widetilde{vY_k}=\tfrac12(\hat v\hat Y_k+\overline{vY_k})$ (B.8b,c); $\hat c=\sqrt{\widetilde{Y_1c_1^2}+\widetilde{Y_2c_2^2}}$, $\widetilde{Y_kc_k^2}=\kappa_k(\widehat{Y_kH_k}-\hat{\mathcal K}\hat Y_k)+\chi_k\hat Y_k$, $\hat{\mathcal K}=\frac{\hat u^2+\hat v^2}{2}$ (B.8d); $\hat V=\kappa_1\widetilde{vY_1}+\kappa_2\widetilde{vY_2}$ (B.8e).

---

## 7. Test cases (Section 5)

General: CFL called Courant number $\varkappa$. Relaxation levels compared: "$p$ relaxation" (mechanical only), "$p$–$pT$" (mechanical + thermal at interfaces, no mass transfer), "$p$–$pT$–$pTg$" (paper's physical model: Steps 2, 3a, 3b with activation (10)), "$p$–$pTg$" (alternative: skip 3a; go from Step 2 directly to 3b when $\epsilon_I\le\alpha_1\le1-\epsilon_I$ and $T_k>T_{sat}(p_{eq})$ for either $k=1$ or $k=2$, refs [8,22]).

**Table 1 — SG parameters, water:**

| Phase | $\gamma$ | $\pi$ [Pa] | $\eta$ [J/kg] | $\eta'$ [J/(kg K)] | $C_v$ [J/(kg K)] |
|---|---|---|---|---|---|
| liquid | 2.35 | $10^9$ | $-1167\times10^3$ | 0 | 1816 |
| vapor | 1.43 | 0 | $2030\times10^3$ | $-23.4\times10^3$ | 1040 |

**Table 2 — SG parameters, dodecane:**

| Phase | $\gamma$ | $\pi$ [Pa] | $\eta$ [J/kg] | $\eta'$ [J/(kg K)] | $C_v$ [J/(kg K)] |
|---|---|---|---|---|---|
| liquid | 2.35 | $4\times10^8$ | $-775.269\times10^3$ | 0 | 1077.7 |
| vapor | 1.025 | 0 | $-237.547\times10^3$ | $-24.4\times10^3$ | 1956.45 |

### 5.1 Water cavitation tube (1D, from Saurel et al. [8])
- Unit-length tube, liquid $\rho_l=1150$ kg/m³, $p=10^5$ Pa, uniform small vapour $\alpha_v=10^{-2}$ everywhere. Water SG (Table 1). From (3): $T_l=354.728$ K; vapour density $\rho_v=0.63$ kg/m³ assuming $T_v=T_l$.
- Velocity discontinuity at $x=0.5$ m.
- **Case a:** $u=-2$ m/s left, $u=2$ m/s right. 5000 cells (as in [17]), $\varkappa=0.5$ ([17] needed $\varkappa\approx0.15$). Results at $t=3.2$ ms (Fig. 1): $p$, $u$, $\alpha_v$, $Y_v$. Two rarefactions → pressure drop. $p$–$pT$ and $p$–$pT$–$pTg$ leading rarefaction edges slower than $p$ and $p$–$pTg$ (thermal transfer activated → lower equilibrium sound speed, consistent with Flåtten–Lund [33]). In $p$–$pT$–$pTg$ and $p$–$pTg$: two additional evaporation waves (slower) — liquid expanded until saturation pressure reached (pressure plateau ~0.5×10^5 Pa in Fig.1), small $Y_v$ generated (~1.3×10⁻⁴ peak). Pure $p$ relaxation: pressure dips to ~0.05×10^5.
- **Case b:** $u=\mp500$ m/s. 5000 cells, $\varkappa=0.1$ (higher Courant deteriorates pressure accuracy at evaporation fronts due to stiffness). Fig. 2 at $t=0.58$ ms: $\alpha_v\to1$ in centre, $Y_v$ up to ~0.145; evaporation fronts clearly visible. Comparison: Zein et al. [17] needed $\varkappa=0.03$. Agrees qualitatively with [8,17].

### 5.2 Dodecane liquid–vapour shock tube (1D, from [8])
- Unit length, discontinuity at $x=0.75$ m; left liquid dodecane, right vapour dodecane; each contains small amount of other phase, $\epsilon=10^{-8}$.
- $(\rho_v,\rho_l,u,p,\alpha_v)_L=(2\ \mathrm{kg/m^3},500\ \mathrm{kg/m^3},0,10^8\ \mathrm{Pa},10^{-8})$, $(\rho_v,\rho_l,u,p,\alpha_v)_R=(2,500,0,10^5\ \mathrm{Pa},1-10^{-8})$.
- Table 2 SG parameters. $t=473\ \mu$s (Fig. 3: $\rho$ (log), $u$, $\alpha_v$, $Y_v$, $p$ (log, bar), $T_v-T_l$). 2000 cells. Grid/CFL: [the text says "computation performed by using a grid with 2000 cells"; CFL not stated [uncertain]].
- Levels: $p$, $p$–$pT$–$pTg$, $p$–$pTg$. Structure: left rarefaction, right-going contact, shock; with phase change an extra evaporation wave between rarefaction and contact, producing a liquid–vapour saturation region at higher speed (velocity ~310 m/s vs ~140 m/s). No visible difference between $p$–$pT$–$pTg$ and $p$–$pTg$. In the rarefaction (almost pure liquid) only mechanical relaxation acts → flow out of thermal equilibrium.
- **Corrigendum (JCP 288 (2015) 196–197, doi:10.1016/j.jcp.2015.01.036):** the $T_v-T_l$ panel of Fig. 3 was wrong (post-processing visualisation settings, not the method). Corrected plot shows $T_v-T_l$ ≈ 10 (×10^5? axis label hard to read) on the left plateau, ~9 in the middle, dropping to 0 near $x\approx0.76$–0.8; the $p$ relaxation curve slightly lower (~8.8) and extends to ~0.8 [values/multiplier **uncertain** from image]. Thanks to T. Flåtten and S. L. Bore.

### 5.3 Cavitating Richtmyer–Meshkov instability (2D, from [9]; no heat/mass transfer)
- Domain $[0,3]\times[0,1]$ m; curved interface; left nearly pure gas ($\alpha_w=10^{-6}$), right nearly pure water ($\alpha_g=10^{-6}$) [left-zone/right-zone assignment as in text: "left zone filled with nearly pure water ($\alpha_g=10^{-6}$) and a right zone filled with nearly pure gas ($\alpha_w=10^{-6}$)" — **as printed**; but Fig. 4 shows gas (yellow, $\alpha_g\approx1$) on the left and water on the right with the interface bulging into water [inconsistency, uncertain]].
- $p=10^5$ Pa; $\rho_g=100$ kg/m³, $\rho_w=1000$ kg/m³; all $u=-200$ m/s.
- SG: water $\gamma_w=4.4$, $\pi_w=6\times10^8$ Pa, $\eta_w=\eta'_w=0$; gas $\gamma_g=1.4$, $\pi_g=0$, $\eta_g=\eta'_g=0$.
- Top, bottom, left: solid walls; right: outflow. Grid 300×100. Outputs at $t=2,3.1,6.4,8.6$ ms (Fig. 4: gas volume fraction, mixture pressure contours). Flow impinging on left wall → right-going shock through curved interface → RM instability, expansion waves, elongating jet; expansions near left wall produce cavitation pockets. Agrees with [9].

### 5.4 Vapour bubble compression (2D, setup like Caro et al. [53], see [14])
- Unit square, stationary vapour bubble radius 0.2 m at centre, liquid-filled; wall-mounted piston at left boundary, fixed walls elsewhere. Dodecane (Table 2). $p_l=p_v=10^5$ Pa, $T_l=T_v=600$ K → $\rho_l=458.338$, $\rho_v=3.408$ kg/m³. $\alpha_l=10^{-6}$ in bubble, $\alpha_v=10^{-6}$ in liquid. Piston velocity $u_p=100$ m/s (moving boundary treatment per [19]).
- Full thermo-chemical model; 200×200 grid; Fig. 5 at $t=0.4,0.8,1.2$ ms: vapour mass fraction and pressure contours; bubble deforms to kidney shape and condenses; no spurious pressure oscillations at interfaces.

### 5.5 High-velocity underwater projectile (2D mapped grid, from [8], also [10])
- Water at $u=600$ m/s flows left→right over a hexagonal-section obstacle; $p=10^5$ Pa, $\rho_l=1500$ kg/m³ [as printed], $\alpha_v=10^{-3}$. Water SG (Table 1).
- Symmetry: compute upper half, symmetry BC at $y=0$. 600×200 cells over $[0,1.5]\times[0,0.5]$ mapped to curvilinear grid (Fig. 6).
- With/without heat and mass transfer. Steady state at $t=10$ ms (Fig. 7: $\alpha_v$, $Y_v$, $\rho$, $p$). Without thermo-chemical effects: only mechanical cavitation, $Y_v$ constant, pressure in wake ~$10^{-3}$ bar; with it, evaporation fronts, pressure decreases until saturation ~0.5 bar. Consistent with [8].

### 5.6 High-pressure fuel injector (2D mapped grid, from [8])
- Nozzle length 10 cm, height 4 cm, throat height 1.2 cm; converging/diverging angles 45° and 10°. Initial discontinuity at $\bar x=0.8$ cm: left liquid dodecane at 1000 bar and 640 K (density 570 kg/m³); right dodecane vapour at atmospheric conditions (1 bar). $\alpha_v=10^{-4}$ in the liquid initially. Table 2 SG parameters.
- 400×160 grid; $t=600\ \mu$s (Fig. 8: $\alpha_v$, $Y_v$, $\rho$, vapour temperature, $p$). Cavitation pocket size similar with/without phase transition, but $Y_v$ and vapour temperature differ markedly; without thermo-chemical effects, sharp vapour-temperature variation across interfaces (vapour T up to ~1000+ K vs 500–580 K with transfer).

---

## 8. Conclusions (Section 6)

Phasic-total-energy formulation ⇒ conservation-consistent mixture energy without an extra equation; relaxation by simple algebraic solvers; consistent equilibrium pressure; robust, allows larger Courant numbers than [17]. Future: powder compaction [54], low-Mach strategies (large impedance contrast liquid/vapour).

---

## 9. Implementation remarks (for the 6-eq IMEX / TV-split project)

1. **Use $\alpha_kE_k$ (phasic total energy) as unknowns.** Discretise the conservative part $\nabla\cdot(\alpha_k(E_k+p_k)\vec u)$ conservatively and the $\pm\Sigma$ term symmetrically (or simply drop it in the hyperbolic step, as the authors do with HLLC). Then $\sum_k\alpha_kE_k$ is exactly conserved and the relaxed pressure is EOS-consistent. This directly maps onto an IMEX split: explicit hyperbolic part on $(\alpha_1,\alpha_k\rho_k,\rho u,\alpha_kE_k)$; "implicit" part = algebraic instantaneous relaxation.
2. **Dropping $\Sigma$ is only OK because relaxation follows every step** (the phasic energies are reset from $p^*,\alpha_1^*$). Without stiff pressure relaxation everywhere, $\Sigma$ would be needed. If one uses finite $\mu$ (non-instantaneous) relaxation in an IMEX scheme, this argument weakens [my remark].
3. **Pressure relaxation is applied everywhere** (all cells, every step), with the closed-form quadratic (C.1)–(C.2). $p_I^0$ uses impedance-weighted interface pressure from the post-hyperbolic state. Check $b^2-4ad>0$ and $\alpha_1^*\in(0,1)$; negative pressures (liquid under tension, $\pi_l\gg|p|$) can make $p_k^0<0$ — the paper asserts the root is always admissible, but guard anyway.
4. **Thermal & chemical relaxation only at interfaces** $\epsilon_I\le\alpha_1\le1-\epsilon_I$ ($\epsilon_I=10^{-4}$), and chemical only when metastable $T_{liq}>T_{sat}(p)$; $T_{sat}(p)$ from Newton on (11a). In the "$p$–$pTg$" variant thermal step is skipped and $pTg$ activated if $T_k>T_{sat}(p_{eq})$ for either phase.
5. **pTg solver:** nested: for given $p$, $T(p)$ from quadratic (C.6)–(C.7); outer Newton on $G(p)=$ (C.8). Need good initial guess (e.g. $p^{**}$) and safeguarding (bracketing; $p+\pi_2>0$, $T>0$, $a_p\ne0$). Then $\alpha_1$ from (C.9); if $\alpha_1\notin[0,1]$ ⇒ abnormal case: fix $\alpha_k=1-\epsilon$ ($\epsilon=10^{-8}$) for the higher-entropy species $k$ and solve p–T equilibrium only, choosing the root maximising mixture entropy.
6. **Relaxations conserve** $\alpha_k\rho_k$ (pressure, thermal) or $\rho$ (chemical), $\rho\vec u$, $E$, $\mathcal E$. Kinetic energy is untouched (single velocity).
7. HLLC (Davis wave speeds with frozen 6-eq sound speed $c=\sqrt{Y_1c_1^2+Y_2c_2^2}$, monotone — no Wood-speed issues) is more robust than Roe. Roe eigenstructure only needed for transverse fluctuations in 2D wave propagation.
8. Stiffness at evaporation fronts: they used $\varkappa=0.5$ (weak expansion) but $\varkappa=0.1$ for strong expansion ($u=\pm500$) for accuracy, not stability per se. Still much larger than Zein et al. (0.15 / 0.03).
9. Residual phase volume fractions used: $10^{-2}$ (water tube, physical vapour content), $10^{-8}$ (dodecane), $10^{-6}$ (RM, bubble), $10^{-3}$ (projectile), $10^{-4}$ (injector).
10. Possible sign/notation issues to double-check (marked [uncertain] above): $p_\iota$ in HLLC energy rows of (28); $\Pi_k$ sign in App. A vs B; $\Delta(\alpha_2\mathcal E_2)$ vs $\Delta(\alpha_2E_2)$ in (B.5b); cross terms in $b',d'$ of (C.3) (re-derive symbolically).
