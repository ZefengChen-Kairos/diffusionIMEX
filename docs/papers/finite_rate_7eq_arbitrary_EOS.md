# Efficient solution of mechanical and thermo-chemical finite-rate relaxation for compressible liquid-vapor flows under arbitrary equations of state

- **Authors:** Giuseppe Sirianni, Barbara Re (corresponding author, barbara.re@polimi.it)
- **Affiliation:** Department of Aerospace Science and Technology, Politecnico di Milano, via La Masa 34, 20156 Milano, Italy
- **Venue:** *Journal of Computational Physics*, Vol. 545 (2026), article 114467
- **DOI:** https://doi.org/10.1016/j.jcp.2025.114467
- **History:** received 18 March 2025; revised 22 August 2025; accepted 15 October 2025; available online 21 October 2025. Open access (CC BY).
- **Keywords:** Baer-Nunziato model, mass transfer, arbitrary-rate relaxation, stiff ODE solver, compressible liquid-vapor flows, nonlinear relaxation parameters.

> Reading note for the diffusionIMEX project (6-equation Saurel/Pelanti-type model with TV-split IMEX and finite-rate relaxation $\mu,\theta,\nu$). Equation numbers $(n)$ are the paper's. Things I could not read with certainty from the PDF are marked **[uncertain]**.

---

## 0. One-paragraph summary

The paper solves the **full 7-equation Baer-Nunziato (BN) model** with **finite-rate** velocity ($\lambda$), pressure ($\mu$), temperature ($\theta$) and Gibbs-free-energy ($\nu$) relaxation, for **any EOS** (tested with stiffened gas and Span-Wagner multiparameter EOS through CoolProp). Time integration is a Strang split: hyperbolic step (HLLC + MUSCL) then a relaxation ODE step. The relaxation step is the contribution: an **adaptive sub-stepping solver** (extending Chiocchetti & Müller [30], who did only mechanical relaxation) in which, in each sub-step, all relaxation coefficients are frozen at a midpoint state $\bar{\mathbf w}=(\mathbf w^m+\mathbf w^{m+1})/2$ (found by fixed-point iteration), the ODEs are rewritten in terms of the **interphase disequilibria** $\Delta_1P,\Delta_1T,\Delta_1g$ (which then obey a *linear* system with constant coefficients), and that linear system is solved **analytically** (closed form generated with a symbolic toolbox, shipped as a `.cpp` file). An extra physically-motivated assumption — pressure relaxes faster than thermo-chemical relaxation — decouples $(\Delta_1T,\Delta_1g)$ from $\Delta_1P$, making the linear system sequential (2x2 then scalar). The sub-step size is controlled by a weighted linearisation-error measure $\delta$. In 0-D tests it matches LSODA to similar error with only 3-10% of the time steps; in full BN computations it reproduces literature results (dodecane and water Riemann problems, water depressurization with Span-Wagner, two-phase nozzle with MDM droplets using nonlinear relaxation-rate closures).

---

## 1. Context / motivation (Sec. 1)

- Hierarchy of BN-derived models: 7-eq full disequilibrium [2]; 6-eq single-velocity (Saurel et al. [4]); 5-eq Kapila [5] (pressure + velocity equilibrium); 4-eq [6-8] (+ thermal equilibrium; only chemical disequilibrium); HRM (Bilicki et al. [9]), 3-eq HEM [10,11].
- BN usually solved by operator splitting [2,16,22]; relaxation ODE system is stiff due to multiple separate time scales. RK45 / LSODA [24] "prohibitively expensive" for CFD.
- Most earlier works assume infinitely fast relaxation; the instantaneous assumption cannot predict metastable states, so finite-rate is wanted. Arbitrary-rate mechanical relaxation already in [25,26].
- **Pelanti and coworkers [21,27,28]**: efficient arbitrary-rate *chemical* relaxation in the single-velocity 6-eq model, general EOS, using **semi-analytical exponential solutions** for the chemical part while assuming mechanical and thermal equilibrium. Authors' critique: the exponential solutions are assumed a priori or derived from reduced models (Kapila/4-eq), so the relaxation ODEs are *not* solved time-accurately; cannot be extended directly to full BN; kinetic/potential energy exchange in mass transfer neglected by construction (velocity equilibrium) [27].
- **Chiocchetti & Müller [30]**: ODE integrator for the mechanical relaxation operator of BN based on the analytical solution of a linearised version of the relaxation ODEs + adaptive sub-stepping. This paper extends the idea to concurrent mechanical + thermo-chemical relaxation, claimed to be the first such operator. The linearisation needs extra assumptions because thermo-chemical relaxation tightly couples the ODEs. They replicate one of Pelanti's tests for comparison.
- Focus: the solver and verification; relaxation-parameter closures are not studied in depth, but the last test (nozzle) uses non-constant (nonlinear) parameters.

---

## 2. Governing equations (Sec. 2)

### 2.1 BN model with relaxation, Eq. (1)

$$
\frac{\partial \mathbf q}{\partial t} + \nabla\cdot\mathbf F + \mathbf G\cdot\nabla\alpha_1 = \mathbf S, \tag{1}
$$

$$
\mathbf q=\begin{bmatrix}\alpha_1\\ \alpha\rho_k\\ \alpha\rho\mathbf u_k\\ \alpha E^t_k\end{bmatrix},\quad
\mathbf F=\begin{bmatrix}\mathbf 0\\ \alpha\rho\mathbf u_k\\ \alpha_k(\rho_k\mathbf u_k\otimes\mathbf u_k+P_k\mathbf I)\\ \alpha_k(E^t_k+P_k)\mathbf u_k\end{bmatrix},\quad
\mathbf G=\begin{bmatrix}\mathbf u_I\\ \mathbf 0\\ (-1)^kP_I\mathbf I\\ (-1)^kP_I\mathbf u_I\end{bmatrix},
$$

$$
\mathbf S=\begin{bmatrix}
\mu\Delta_1P-\dfrac{\nu}{\rho_I}\Delta_1 g\\[4pt]
-\nu\Delta_k g\\[2pt]
-\lambda\Delta_k\mathbf u-\nu\mathbf u_I\Delta_k g\\[2pt]
-\mu P_I\Delta_kP-\lambda\mathbf u_I\cdot\Delta_k\mathbf u-\nu e^t_I\Delta_k g-\theta\Delta_kT
\end{bmatrix}.
$$

Variables / notation:
- $k\in\{1,2\}$; **$k=1$ vapour, $k=2$ liquid**. Subscript $(\cdot)_k$ at the end of a group applies to the whole group: $\alpha\rho_k\equiv\alpha_k\rho_k$, $\alpha\rho\mathbf u_k\equiv\alpha_k\rho_k\mathbf u_k$.
- $\alpha_k$ volume fractions ($\alpha_2=1-\alpha_1$), $\alpha\rho_k$ partial densities, $\alpha\rho\mathbf u_k$ partial momenta, $\alpha E^t_k=\alpha\rho_k(e_k+|\mathbf u_k|^2/2)$ partial total energies.
- $P_k$ phasic pressures, $T_k$ temperatures, $g_k=e_k+P_k/\rho_k-T_ks_k$ Gibbs free energy per unit mass, $s_k$ specific entropy.
- **Jump notation:** $\Delta_1 y=-\Delta_2 y=(y_1-y_2)$, i.e. $\Delta_k y = y_k - y_{k'}$.
- Relaxation parameters (positive, finite): $\mu$ pressure, $\lambda$ velocity, $\nu$ Gibbs free energy, $\theta$ temperature.
- Mass transfer $\dot m=-\nu\Delta_1 g$; heat transfer $Q=-\theta\Delta_1T$.
- **Interface closures:** original BN: $\mathbf u_I=\mathbf u_2$ (liquid velocity), $P_I=P_1$ (vapour pressure). Each phase has its own EOS.
- Sign convention: with $\mathbf G$ on the LHS and $(-1)^k$, phase 1 momentum gets $+P_I\nabla\alpha_1$ on the RHS (standard BN).

Shorthand for partial derivatives of phasic internal energy, Eq. (2):
$$
(e_{\rho,\varphi})_k=\left(\frac{\partial e_k}{\partial\rho_k}\right)_{\varphi_k},\qquad
(e_{\varphi,\rho})_k=\left(\frac{\partial e_k}{\partial\varphi_k}\right)_{\rho_k}, \tag{2}
$$
with $\varphi_k$ a generic phasic thermodynamic variable (temperature, pressure or entropy...).

### 2.2 Interfacial density and energy closures, Eq. (3)

Thermo-chemical sources need $\rho_I$ and $e_I = e^t_I-|\mathbf u_I|^2/2$. The BN model doesn't assume mechanical/thermal equilibrium, so there is "no straightforward way" to close them; the authors borrow closures derived for reduced models:
- $\rho_I$ from Saurel et al. [15] (mechanical equilibrium; satisfies 2nd law for the mixture);
- $e_I$ from Zein et al. [16] (thermo-mechanical equilibrium).

$$
\rho_I=\frac{\sum_k\big(\rho_kc_k^2/\alpha_k\big)}{\sum_k\big(c_k^2/\alpha_k\big)},\qquad
e_I=\frac{1}{\sum_k(\alpha\rho_kc_{v,k})^{-1}}\sum_k\left[\frac{e_k+\rho_k(e_{\rho,T})_k(1-\rho_k/\rho_I)}{\alpha\rho_kc_{v,k}}\right], \tag{3}
$$
$c_{v,k}$ phasic isochoric specific heat, $c_k$ phasic sound speed. The scheme does not depend on these closures; better closures flagged as future work.

### 2.3 Operator splitting, Eq. (4)

Strang splitting between convection and relaxation over $[t^n,t^{n+1}]$:
- (4a) hyperbolic operator $\mathbf q^{*,n+1}=\mathcal H(\mathbf q^n)$: $\partial_t\mathbf q+\nabla\cdot\mathbf F+\mathbf G\cdot\nabla\alpha_1=\mathbf 0$;
- (4b) relaxation operator $\mathbf q^{n+1}=\mathcal R(\mathbf q^{*,n+1})$: $d\mathbf q/dt=\mathbf S$.

(The text says "Strang splitting" but the summary (4a)-(4b) is written as a single H then R sequence. **[uncertain]** whether a symmetric half-step R/H/R arrangement is actually used.)

Hyperbolic step: **HLLC-type Riemann solver of Lochon et al. [31]** (for BN), second order via **Schwendeman [32] MUSCL** reconstruction with **Barth-Jespersen [33] limiter**; time update with the conservative explicit time-integration scheme presented for Euler in [34], recently extended to BN [35].

### 2.4 Thermodynamic models (Sec. 2.1)

The relaxation operator is valid for **any EOS**; two are used.

**Stiffened gas (SG)**, Eq. (5)-(6):
$$
P_k(\rho_k,e_k)=(\gamma_k-1)\rho_k(e_k-q_{0,k})-\gamma_kP_{\infty,k},\qquad
T_k(\rho_k,P_k)=\frac{P_k+P_{\infty,k}}{(\gamma_k-1)c_{v,k}\rho_k}, \tag{5}
$$
$$
g_k(T_k,P_k)=(\gamma_kc_{v,k}-q'_{0,k})T_k-c_{v,k}T_k\log\!\left(\frac{T_k^{\gamma_k}}{(P_k+P_{\infty,k})^{\gamma_k-1}}\right)+q_{0,k}, \tag{6}
$$
five parameters per fluid ($\gamma,P_\infty,q_0,q'_0,c_v$) from experimental data [15].

**Span-Wagner (SW)** multiparameter EOS [36], via **CoolProp** [37]. Dimensionless Helmholtz free energy, Eq. (7):
$$
\begin{aligned}
\Phi_k(\delta_k,\tau_k)=\;&\log\delta_k+a^{(1)}_k\log\tau_k+a^{(2)}_k+a^{(3)}_k\tau_k+a^{(4)}_k\tau_k^{-1}+a^{(5)}_k\tau_k^{-2}+a^{(6)}_k\tau_k^{-3}+a^{(7)}_k\log\!\big(1-e^{-a^{(8)}_k\tau_k}\big)\\
&+\sum_{j=1}^{6}N^{(j)}_k\delta_k^{b^{(j)}_k}\tau_k^{c^{(j)}_k}
+\sum_{j=7}^{32}N^{(j)}_k\delta_k^{b^{(j)}_k}\tau_k^{c^{(j)}_k}e^{-\delta_k^{d^{(j)}_k}}
+\sum_{j=33}^{36}N^{(j)}_k\delta_k^{b^{(j)}_k}\tau_k^{c^{(j)}_k}e^{-\phi^{(j)}_k(\delta_k-1)^2-\beta^{(j)}_k(\tau_k-\kappa^{(j)}_k)^2},
\end{aligned}\tag{7}
$$
$\delta_k=\rho_k/\rho^c_k$, $\tau_k=T^c_k/T_k$ (scaled with critical values). All thermodynamic properties are analytical derivatives of (7); $g_k=e_k+P_k/\rho_k-T_ks_k$. (Superscript limits/indices on the exponential terms read from a small figure: **[uncertain]** exact symbol names.) Note: SW has **$(\rho,T)$ as natural variables**, which drives the choice $\varphi_k=T_k$ (Sec. 3).

---

## 3. Relaxation ODE system (Sec. 2.2)

### 3.1 Conservative form, Eq. (8)

$$
\begin{cases}
\dfrac{d\alpha_1}{dt}=\mu\Delta_1P-\dfrac{\nu}{\rho_I}\Delta_1g\\[6pt]
\dfrac{d\alpha\rho_k}{dt}=-\nu\Delta_kg\\[6pt]
\dfrac{d\alpha\rho\mathbf u_k}{dt}=-\lambda\Delta_k\mathbf u-\nu\mathbf u_I\Delta_kg\\[6pt]
\dfrac{d\alpha E^t_k}{dt}=-\mu P_I\Delta_kP-\lambda\mathbf u_I\cdot\Delta_k\mathbf u-\nu e^t_I\Delta_kg-\theta\Delta_kT.
\end{cases}\tag{8}
$$

**Mixture conservation is built in**: since $\Delta_1y=-\Delta_2y$, the sources cancel when summed over $k$. Hence $\rho_m=\alpha\rho_1+\alpha\rho_2$, $\rho\mathbf u_m=\alpha\rho\mathbf u_1+\alpha\rho\mathbf u_2$ and $E^t_m=\alpha E^t_1+\alpha E^t_2$ are constant during relaxation. (Mass: the solver stores only $\alpha\rho_1$ and recovers $\alpha\rho_2=\rho_m-\alpha\rho_1$, so mixture mass is exact. **Energy/momentum: not addressed explicitly in the paper.** The solver evolves $\mathbf u_1,\mathbf u_2,T_1,T_2$ independently from linearised ODEs with coefficients frozen at $\bar{\mathbf w}$; with frozen, phase-specific $\bar K$ the discrete updates of $\alpha E^t_1$ and $\alpha E^t_2$ need not cancel exactly, so $E^t_m$ (and $\rho\mathbf u_m$, since $K_{\mathbf u_I}^{\mathbf u_k}$ contains $1/\alpha\rho_k$ frozen) is conserved only up to the linearisation error controlled by $\delta^{\max}$. Nowhere in the paper is an a-posteriori re-projection onto $E^t_m$ described — **[my inference; the paper is silent]**. Ref. [35] (Sirianni, Re, Abgrall, "Mixture-conservative temperature-based BN solver", C&F 2025) is the authors' hyperbolic solver that handles conservation with temperature as variable; the relaxation step here might rely on something similar but this is not stated.)

### 3.2 Primitive-variable form, Eq. (9) (derivation in Appendix A)

Assume each phasic EOS is written $e_k=e_k(\rho_k,\varphi_k)$, $\varphi_k\in\{P_k,T_k,g_k,\dots\}$:

$$
\begin{cases}
\dfrac{d\alpha_1}{dt}=\mu\Delta_1P-\dfrac{\nu}{\rho_I}\Delta_1g\\[6pt]
\dfrac{d\alpha\rho_k}{dt}=-\nu\Delta_kg\\[6pt]
\dfrac{d\mathbf u_k}{dt}=-\dfrac{\lambda}{\alpha\rho_k}\Delta_k\mathbf u+\dfrac{\nu\Delta_kg}{\alpha\rho_k}(\mathbf u_k-\mathbf u_I)\\[6pt]
\dfrac{d\varphi_k}{dt}=\left(\dfrac{\rho_k^2(e_{\rho,\varphi})_k-P_I}{\alpha\rho_k(e_{\varphi,\rho})_k}\mu\right)\Delta_kP
-\left(\dfrac{\theta}{\alpha\rho_k(e_{\varphi,\rho})_k}\right)\Delta_kT\\[6pt]
\qquad+\left\{\dfrac{\nu}{\alpha\rho_k(e_{\varphi,\rho})_k}\left[(e^t_k-e^t_I)+\left(\rho_k-\dfrac{\rho_k^2}{\rho_I}\right)(e_{\rho,\varphi})_k+\mathbf u_k\cdot(\mathbf u_I-\mathbf u_k)\right]\right\}\Delta_kg\\[6pt]
\qquad+\left(\dfrac{\lambda}{\alpha\rho_k(e_{\varphi,\rho})_k}\right)(\mathbf u_k-\mathbf u_I)\cdot\Delta_k\mathbf u.
\end{cases}\tag{9}
$$

(Sanity check I did: with fixed phasic mass, $d\rho_k=-(\rho_k/\alpha_k)d\alpha_k$ and $\alpha\rho_k\,de_k=-P_Id\alpha_k$ give $\alpha\rho_k(e_{\varphi,\rho})_kd\varphi_k=(\rho_k^2(e_{\rho,\varphi})_k-P_I)d\alpha_k$ — consistent with the $\mu$ term, since $d\alpha_k/dt=\mu\Delta_kP$.)

### 3.3 Compact form with coefficients $K$, Eq. (10)

$$
\begin{aligned}
\frac{d\alpha_1}{dt}&=K_P^{\alpha_1}\Delta_1P+K_g^{\alpha_1}\Delta_1g &(10a)\\
\frac{d\alpha\rho_k}{dt}&=K_g^{\alpha\rho_k}\Delta_kg &(10b)\\
\frac{d\mathbf u_k}{dt}&=K_{\mathbf u}^{\mathbf u_k}\Delta_k\mathbf u+K_{\mathbf u_I}^{\mathbf u_k}(\mathbf u_k-\mathbf u_I) &(10c)\\
\frac{d\varphi_k}{dt}&=K_P^{\varphi_k}\Delta_kP+K_T^{\varphi_k}\Delta_kT+K_g^{\varphi_k}\Delta_kg+K_{\mathbf u}^{\varphi_k}(\mathbf u_k-\mathbf u_I)\cdot\Delta_k\mathbf u &(10d)
\end{aligned}
$$

$K^{eq}_y$ = coefficient multiplying jump $\Delta_k y$ in equation $eq$. E.g. $K_P^{\alpha_1}=\mu$, $K_g^{\alpha_1}=-\nu/\rho_I$, $K_g^{\alpha\rho_k}=-\nu$, $K_{\mathbf u}^{\mathbf u_k}=-\lambda/\alpha\rho_k$, $K_{\mathbf u_I}^{\mathbf u_k}=\nu\Delta_kg/\alpha\rho_k$ (note: contains $\Delta_k g$ itself — a nonlinear coefficient), $K_T^{\varphi_k}=-\theta/(\alpha\rho_k(e_{\varphi,\rho})_k)$. All listed in Table A.8 (Appendix A). Every coefficient is in principle time-dependent.

Remarks from the paper: three thermodynamic variables per phase relax (P, T, g) and are linked through the EOS; Gibbs free energy is highly nonlinear (especially for cubic / multiparameter EOS); Gibbs relaxation also produces mass and momentum transfer ($K_g^{\alpha\rho_k}\Delta_kg$, $K_{\mathbf u_I}^{\mathbf u_k}(\mathbf u_k-\mathbf u_I)$), so system (10) is much more coupled than the mechanical-only case.

### 3.4 ODEs for interphase disequilibria, Eqs. (11)-(13) (idea from Pelanti [28])

$$
\frac{d(\Delta_1\varphi)}{dt}=\frac{d\varphi_1}{dt}-\frac{d\varphi_2}{dt}
=(K_P^{\varphi_1}+K_P^{\varphi_2})\Delta_1P+(K_T^{\varphi_1}+K_T^{\varphi_2})\Delta_1T+(K_g^{\varphi_1}+K_g^{\varphi_2})\Delta_1g
+\big[K_{\mathbf u}^{\varphi_1}\mathbf u_1+K_{\mathbf u}^{\varphi_2}\mathbf u_2-(K_{\mathbf u}^{\varphi_1}+K_{\mathbf u}^{\varphi_2})\mathbf u_I\big]\cdot\Delta_1\mathbf u \tag{11}
$$
(sums appear because $\Delta_2=-\Delta_1$). Define $K_y^{\Delta_1\varphi}=K_y^{\varphi_1}+K_y^{\varphi_2}$ and
$S_{\mathbf u}^{\Delta_1\varphi}(\mathbf u_1,\mathbf u_2)=\big[K_{\mathbf u}^{\varphi_1}\mathbf u_1+K_{\mathbf u}^{\varphi_2}\mathbf u_2-(K_{\mathbf u}^{\varphi_1}+K_{\mathbf u}^{\varphi_2})\mathbf u_I\big]\cdot\Delta_1\mathbf u$:
$$
\frac{d(\Delta_1\varphi)}{dt}=K_P^{\Delta_1\varphi}\Delta_1P+K_T^{\Delta_1\varphi}\Delta_1T+K_g^{\Delta_1\varphi}\Delta_1g+S_{\mathbf u}^{\Delta_1\varphi}(\mathbf u_1,\mathbf u_2). \tag{12}
$$
$S_{\mathbf u}$ is nonlinear in both velocities. Writing (12) for $\varphi=P,T,g$ (i.e. using each as the generic variable in turn):
$$
\begin{cases}
\dfrac{d(\Delta_1P)}{dt}=K_P^{\Delta_1P}\Delta_1P+K_T^{\Delta_1P}\Delta_1T+K_g^{\Delta_1P}\Delta_1g+S_{\mathbf u}^{\Delta_1P}(\mathbf u_1,\mathbf u_2)\\[6pt]
\dfrac{d(\Delta_1T)}{dt}=K_P^{\Delta_1T}\Delta_1P+K_T^{\Delta_1T}\Delta_1T+K_g^{\Delta_1T}\Delta_1g+S_{\mathbf u}^{\Delta_1T}(\mathbf u_1,\mathbf u_2)\\[6pt]
\dfrac{d(\Delta_1g)}{dt}=K_P^{\Delta_1g}\Delta_1P+K_T^{\Delta_1g}\Delta_1T+K_g^{\Delta_1g}\Delta_1g+S_{\mathbf u}^{\Delta_1g}(\mathbf u_1,\mathbf u_2).
\end{cases}\tag{13}
$$
Key point: all interphase disequilibria $\Delta_1\varphi$ that appear in (10) are unknowns of (13) and appear **linearly** on its RHS.

---

## 4. Relaxation solver (Sec. 3)

### 4.1 Sub-time stepping framework (Fig. 1, Algorithm 1)

Inside $[t^n,t^{n+1}]$ (hyperbolic $\Delta t$), sub-steps $t^m\to t^{m+1}$. At each sub-step a linearised version of the ODEs is solved analytically; a fixed-point iteration determines the linearisation state; the sub-step size adapts to an estimate of the linearisation error.

- State vector: $\mathbf w=[\alpha_1,\alpha\rho_1,\mathbf u_1,\varphi_1,\mathbf u_2,\varphi_2]$. Together with $\alpha_2=1-\alpha_1$ and **mass conservation $\alpha\rho_2=\rho_m-\alpha\rho_1$**, it determines the state uniquely.
- Linearisation state: $\bar{\mathbf w}=(\mathbf w^{m+1}+\mathbf w^m)/2$ (midpoint, implicit -> fixed-point iterations).
- $\tilde{\mathbf w}(t;t^m,\mathbf w^m,\bar{\mathbf w})$: analytical solution at time $t$ of the IVP with data $\mathbf w^m$ at $t^m$, linearised around $\bar{\mathbf w}$.

**Algorithm 1 (pseudocode, transcribed):**
```
w^m = w^{*,n+1}
t^m = t^n ; first guess Δt = 1/|K_P^{Δ1P}(w^m)|          # init sub-stepping
while t^m < t^{n+1}:
    j = 0
    w^{m+1}(j) = w^m                                      # init fixed point
    while r(j) > r_max and j <= j_max:
        j = j+1
        compute w̄                                          # = (w^m + w^{m+1}(j-1))/2
        w^{m+1}(j) = w̃(t^m+Δt; t^m, w^m, w̄)              # analytical solution
        compute r(j)
    if r(j) <= r_max and w^{m+1}(j) physical:
        compute δ
        if δ <= δ_max:                                     # linearisation error OK
            w^m = w^{m+1}(j);  t^m = t^m + Δt
            Δt = min{ t^{n+1}-t^m , 0.8 Δt (δ_max/δ) }
    if r(j) > r_max or δ > δ_max or j > j_max or w^{m+1}(j) unphysical:
        w^{m+1}(0) = w^m ;  Δt = Δt/3                      # reject, shrink
w^{n+1} = w^m
```
Initial sub-step guess $\Delta t=1/|K_P^{\Delta_1P}(\mathbf w^m)|$ = the pressure-relaxation time scale (the fastest).

**Fixed-point convergence**, Eq. (14):
$$
r(j)=\max_i\left(\frac{|w_i^{m+1}(j)-w_i^{m+1}(j-1)|}{|w_i^{m+1}(j)|+|w_i^{m+1}(j-1)|+\varepsilon}\right)\le r^{\max}=10^{-12}. \tag{14}
$$

**Linearisation-error measure** (from [30]), $\delta=\max(\bar\delta,\delta^{m+1})$, Eq. (15):
$$
\bar\delta=\max_i\left(\mathscr W_i\frac{|D_i^m-\bar D_i|}{|D_i^m|+|\bar D_i|+\varepsilon}\right),\qquad
\delta^{m+1}=\max_i\left(\mathscr W_i\frac{|D_i^m-D_i^{m+1}|}{|D_i^m|+|D_i^{m+1}|+\varepsilon}\right), \tag{15}
$$
$D_i^m=D_i(\mathbf w^m)$, $\bar D_i=D_i(\bar{\mathbf w})$, $D_i^{m+1}=D_i(\mathbf w^{m+1}(j))$; $D$ = vector of all linearised quantities on which a linearisation bound is applied; $\mathscr W$ = weights. Choosing these is much harder than in the mechanical-only case — too strict over-constrains (tiny sub-steps). Definitions in Appendix B.

**Admissibility checks**: $0\le\alpha\le1$, densities positive, phasic state not beyond the **spinodal**, i.e. $c_k^2>0$. If inadmissible, reject and recompute from $t^m$ with sub-step reduced to **one third** (user parameter; halving = more conservative/cheaper trade-off; 1/3 "reasonable compromise").

Next step: $\Delta t=0.8\,\Delta t(\delta^{\max}/\delta)$, safety factor 0.8. $j^{\max}=10$, convergence typically in **3-5 iterations**.

### 4.2 Analytical solution of the linearised ODEs (Sec. 3.1)

**Linearisation assumption:** every coefficient is constant within a sub-step, evaluated at $\bar{\mathbf w}$: $K_y^{eq}(\mathbf w)\simeq K_y^{eq}(\bar{\mathbf w})=\bar K_y^{eq}$. Notation: $\tilde{\mathbf w}(t)$ for $\tilde{\mathbf w}(t;t^m,\mathbf w^m,\bar{\mathbf w})$.

**Step 1 — velocities, Eq. (16).** With $K_{\mathbf u_I}^{\mathbf u_k}(\mathbf w)\simeq\bar K_{\mathbf u_I}^{\mathbf u_k}$ (this freezes the $\nu\Delta_kg/\alpha\rho_k$ factor), velocity ODEs decouple from everything else:
$$
\begin{cases}
\dfrac{d\mathbf u_1}{dt}=\bar K_{\mathbf u}^{\mathbf u_1}\Delta_1\mathbf u+\bar K_{\mathbf u_I}^{\mathbf u_1}(\mathbf u_1-\mathbf u_I)\\[6pt]
\dfrac{d\mathbf u_2}{dt}=\bar K_{\mathbf u}^{\mathbf u_2}\Delta_2\mathbf u+\bar K_{\mathbf u_I}^{\mathbf u_2}(\mathbf u_2-\mathbf u_I)
\end{cases}\tag{16}
$$
a 2x2 linear constant-coefficient system (per component) solved exactly -> $\tilde{\mathbf u}_k(t)$. (With $\mathbf u_I=\mathbf u_2$ the last term of the second eq vanishes.)

**Why not linearise (10d) directly:** even with constant $\bar K$, $d\varphi_k/dt$ in (10d) remains highly nonlinear because $\Delta P,\Delta T,\Delta g$ are nonlinear EOS functions of $(\rho_k,\varphi_k)$, Eq. (17):
$$
\begin{cases}
\dfrac{d\varphi_1}{dt}=\bar K_P^{\varphi_1}\big(P_1(\rho_1,\varphi_1)-P_2(\rho_2,\varphi_2)\big)+\bar K_T^{\varphi_1}\big(T_1(\rho_1,\varphi_1)-T_2(\rho_2,\varphi_2)\big)\\
\qquad+\bar K_g^{\varphi_1}\big(g_1(\rho_1,\varphi_1)-g_2(\rho_2,\varphi_2)\big)+\bar K_{\mathbf u}^{\varphi_1}(\mathbf u_1-\mathbf u_I)\cdot\Delta_1\mathbf u\\
\dfrac{d\varphi_2}{dt}=\bar K_P^{\varphi_2}\big(P_2(\rho_2,\varphi_2)-P_1(\rho_1,\varphi_1)\big)+\bar K_T^{\varphi_2}\big(T_2(\rho_2,\varphi_2)-T_1(\rho_1,\varphi_1)\big)\\
\qquad+\bar K_g^{\varphi_2}\big(g_2(\rho_2,\varphi_2)-g_1(\rho_1,\varphi_1)\big)+\bar K_{\mathbf u}^{\varphi_2}(\mathbf u_2-\mathbf u_I)\cdot\Delta_2\mathbf u.
\end{cases}\tag{17}
$$
E.g. for $\varphi_k=P_k$, $\Delta_1T$ and $\Delta_1g$ remain nonlinear in $(\rho,P)$; for $\varphi_k=T_k$, $\Delta_1P$ and $\Delta_1g$ remain nonlinear. Appendix A.1 writes (17) for $\varphi_k=P_k$ with SG to show this even for the simplest EOS.

**Step 2 — disequilibrium system.** Under $K\simeq\bar K$, (13) becomes a linear non-homogeneous 3x3 ODE system in $(\Delta_1P,\Delta_1T,\Delta_1g)$. Increases the number of equations, but avoids further linearisation.

**Time-scale separation assumption** (physical grounds): pressure relaxation is faster than thermo-chemical relaxation, while thermal and chemical relaxation evolve on comparable time scales. This allows a *mild* decoupling -> sequential treatment. *This does not imply instantaneous pressure relaxation in the final solution.* Resulting linearised system, Eq. (18):
$$
\begin{cases}
\begin{cases}
\dfrac{d(\Delta_1T)}{dt}=\cancel{\bar K_P^{\Delta_1T}\Delta_1P}^{\,0}+\bar K_T^{\Delta_1T}\Delta_1T+\bar K_g^{\Delta_1T}\Delta_1g+\bar S_{\mathbf u}^{\Delta_1T}(\tilde{\mathbf u}_1(t),\tilde{\mathbf u}_2(t))\\[6pt]
\dfrac{d(\Delta_1g)}{dt}=\cancel{\bar K_P^{\Delta_1g}\Delta_1P}^{\,0}+\bar K_T^{\Delta_1g}\Delta_1T+\bar K_g^{\Delta_1g}\Delta_1g+\bar S_{\mathbf u}^{\Delta_1g}(\tilde{\mathbf u}_1(t),\tilde{\mathbf u}_2(t))
\end{cases}\\[14pt]
\dfrac{d(\Delta_1P)}{dt}=\bar K_P^{\Delta_1P}\Delta_1P+\bar K_T^{\Delta_1P}\Delta_1T+\bar K_g^{\Delta_1P}\Delta_1g+\bar S_{\mathbf u}^{\Delta_1P}(\tilde{\mathbf u}_1(t),\tilde{\mathbf u}_2(t))
\end{cases}\tag{18}
$$
i.e. the pressure-disequilibrium coupling into the $\Delta T$ and $\Delta g$ equations is dropped (one-way coupling P <- T,g). Thanks to $\bar K_{\mathbf u}^{\varphi_1},\bar K_{\mathbf u}^{\varphi_2}$ constant and known $\tilde{\mathbf u}_k(t)$, the forcing $\bar S^{\Delta_1\varphi}_{\mathbf u}(\tilde{\mathbf u}_1(t),\tilde{\mathbf u}_2(t))$ is a completely known (nonlinear-in-time) function.

**Step 2a — thermo-chemical 2x2 subsystem**, Eq. (19):
$$
\frac{d(\Delta_1T)}{dt}=\bar K_T^{\Delta_1T}\Delta_1T+\bar K_g^{\Delta_1T}\Delta_1g+\bar S_{\mathbf u}^{\Delta_1T}(\tilde{\mathbf u}_1,\tilde{\mathbf u}_2),\qquad
\frac{d(\Delta_1g)}{dt}=\bar K_T^{\Delta_1g}\Delta_1T+\bar K_g^{\Delta_1g}\Delta_1g+\bar S_{\mathbf u}^{\Delta_1g}(\tilde{\mathbf u}_1,\tilde{\mathbf u}_2) \tag{19a,b}
$$
solved exactly -> $\Delta_1\tilde T(t)$, $\Delta_1\tilde g(t)$.

**Step 2b — pressure scalar ODE**, Eq. (20):
$$
\frac{d(\Delta_1P)}{dt}=\bar K_P^{\Delta_1P}\Delta_1P+\bar K_T^{\Delta_1P}\Delta_1\tilde T(t)+\bar K_g^{\Delta_1P}\Delta_1\tilde g(t)+\bar S_{\mathbf u}^{\Delta_1P}(\tilde{\mathbf u}_1,\tilde{\mathbf u}_2) \tag{20}
$$
exact solution -> $\Delta_1\tilde P(t)$.

**Step 3 — update the state by quadrature**, Eq. (21) (the RHS of linearised (10) is now a known function of time):
$$
\begin{aligned}
\mathbf u_k^{m+1}&=\tilde{\mathbf u}_k(t^{m+1}) &(21a)\\
\alpha_1^{m+1}&=\alpha_1^m+\bar K_P^{\alpha_1}\int_{t^m}^{t^{m+1}}\Delta_1\tilde P(t)\,dt+\bar K_g^{\alpha_1}\int_{t^m}^{t^{m+1}}\Delta_1\tilde g(t)\,dt &(21b)\\
\alpha\rho_1^{m+1}&=\alpha\rho_1^m+\bar K_g^{\alpha\rho_1}\int_{t^m}^{t^{m+1}}\Delta_1\tilde g(t)\,dt &(21c)\\
\varphi_k^{m+1}&=\varphi_k^m+\bar K_P^{\varphi_k}\int\Delta_k\tilde P\,dt+\bar K_T^{\varphi_k}\int\Delta_k\tilde T\,dt+\bar K_g^{\varphi_k}\int\Delta_k\tilde g\,dt+\bar K_{\mathbf u}^{\varphi_k}\int(\tilde{\mathbf u}_k-\tilde{\mathbf u}_I)\cdot\Delta_k\tilde{\mathbf u}\,dt &(21d)
\end{aligned}
$$
All integrals computed **analytically with a symbolic toolbox**; result $\tilde{\mathbf w}(t;t^m,\mathbf w^m,\bar{\mathbf w})$ evaluated directly in Algorithm 1 — no numerical quadrature and no nonlinear solve (other than the fixed-point loop). Expressions too long for the paper: supplied as a `.cpp` supplementary file.

**Choice of $\varphi_k$:** phasic **temperature** $T_k$, so $\mathbf w=[\alpha_1,\alpha\rho_1,\mathbf u_1,T_1,\mathbf u_2,T_2]$. Reason: $(\rho,T)$ are the natural variables of the (expensive) SW Helmholtz EOS -> no EOS inversion. Using $P_k$ or $s_k$ would work with minimal change but cost EOS inversions.

**Coefficients vary between sub-steps** (expected); likewise $\mu,\lambda,\theta,\nu$ are never assumed constant and may depend on $t$ and $\mathbf w$ (frozen at $\bar{\mathbf w}$ within a sub-step).

#### Implementation remarks for diffusionIMEX (my notes)
- Limits: as rates $\to0$ all $K\to0$, $\tilde{\mathbf w}=\mathbf w^m$ (frozen), and the first-guess $\Delta t=1/|K_P^{\Delta_1P}|\to\infty$ is clipped by $t^{n+1}-t^m$. As rates $\to\infty$ the exponential solutions of the linear systems tend to their equilibrium (particular) solutions, but the *coefficients* are frozen at $\bar{\mathbf w}$, so accuracy depends on the fixed-point/linearisation-error control; the paper mimics instantaneous relaxation with large finite parameters (e.g. $\mu=10^4$, $\theta=10^{12}$ ...). The analytical exponential solution itself is unconditionally stable for negative eigenvalues, which is what makes large rates affordable.
- For a 6-eq single-velocity model the velocity step (16) and $S_{\mathbf u}$ terms disappear; the core reduces to (19)-(21) with $\varphi_k=T_k$.

---


## 5. Appendix A — derivation for a generic variable $\varphi$, coefficient table, SG example

### A.1-A.7 Generic derivation
Write (8) generically as sums of source contributions $S_j$:
$$
\frac{d\alpha_k}{dt}=\sum_jS_j^{\alpha_k},\quad \frac{d\alpha\rho_k}{dt}=\sum_jS_j^{\alpha\rho_k},\quad \frac{d\alpha\rho\mathbf u_k}{dt}=\sum_j\mathbf S_j^{\alpha\rho\mathbf u_k},\quad \frac{d\alpha E^t_k}{dt}=\sum_jS_j^{\alpha E^t_k}. \tag{A.1}
$$
$$
\frac{d\rho_k}{dt}=\frac1{\alpha_k}\Big(\frac{d\alpha\rho_k}{dt}-\rho_k\frac{d\alpha_k}{dt}\Big)=\frac1{\alpha_k}\Big(\sum_jS_j^{\alpha\rho_k}-\rho_k\sum_jS_j^{\alpha_k}\Big), \tag{A.2}
$$
$$
\frac{d\mathbf u_k}{dt}=\frac1{\alpha\rho_k}\Big(\frac{d\alpha\rho\mathbf u_k}{dt}-\mathbf u_k\frac{d\alpha\rho_k}{dt}\Big)=\frac1{\alpha\rho_k}\Big(\sum_j\mathbf S_j^{\alpha\rho\mathbf u_k}-\mathbf u_k\sum_jS_j^{\alpha\rho_k}\Big), \tag{A.3}
$$
$$
\frac{d\alpha E^t_k}{dt}=\alpha\rho_k\frac{d}{dt}\Big(e_k+\tfrac12\mathbf u_k^2\Big)+e^t_k\frac{d\alpha\rho_k}{dt}
=\alpha\rho_k\Big((e_{\varphi,\rho})_k\frac{d\varphi_k}{dt}+(e_{\rho,\varphi})_k\frac{d\rho_k}{dt}\Big)+\alpha\rho\mathbf u_k\cdot\frac{d\mathbf u_k}{dt}+e^t_k\frac{d\alpha\rho_k}{dt}, \tag{A.4}
$$
($e^t_k=e_k+|\mathbf u_k|^2/2$ specific total energy)
$$
\frac{d\varphi_k}{dt}=-\frac{(e_{\rho,\varphi})_k}{(e_{\varphi,\rho})_k}\frac{d\rho_k}{dt}-\frac{1}{(e_{\varphi,\rho})_k}\mathbf u_k\cdot\frac{d\mathbf u_k}{dt}-\frac{e^t_k}{\alpha\rho_k(e_{\varphi,\rho})_k}\frac{d\alpha\rho_k}{dt}+\frac{1}{\alpha\rho_k(e_{\varphi,\rho})_k}\frac{d\alpha E^t_k}{dt}. \tag{A.5}
$$
Substituting (A.1)-(A.3):
$$
\frac{d\varphi_k}{dt}=-\frac{(e_{\rho,\varphi})_k}{\alpha_k(e_{\varphi,\rho})_k}\Big[\sum_jS_j^{\alpha\rho_k}-\rho_k\sum_jS_j^{\alpha_k}\Big]
-\frac{1}{\alpha\rho_k(e_{\varphi,\rho})_k}\mathbf u_k\cdot\Big[\sum_j\mathbf S_j^{\alpha\rho\mathbf u_k}-\mathbf u_k\sum_jS_j^{\alpha\rho_k}\Big]
-\frac{e^t_k}{\alpha\rho_k(e_{\varphi,\rho})_k}\sum_jS_j^{\alpha\rho_k}+\frac{1}{\alpha\rho_k(e_{\varphi,\rho})_k}\sum_jS_j^{\alpha E^t_k}. \tag{A.6}
$$
Valid for any relaxation model of the form (A.1) — **directly reusable for a 6-eq model** (drop the velocity-relaxation sources, single $\mathbf u$). Inserting the sources of (8) gives (A.7) = (9), with the velocity equation written as $\frac{d\mathbf u_k}{dt}=\frac{1}{\alpha\rho_k}\big[-\lambda\Delta_k\mathbf u+(\mathbf u_k-\mathbf u_I)\nu\Delta_kg\big]$.

### Table A.8 — coefficients $K^{eq}_y$
| $K^{eq}_y$ | definition | $K^{eq}_y$ | definition |
|---|---|---|---|
| $K_P^{\alpha_1}$ | $\mu$ | $K_g^{\alpha_1}$ | $-\nu/\rho_I$ |
| $K_g^{\alpha\rho_k}$ | $-\nu$ | $K_{\mathbf u}^{\mathbf u_k}$ | $-\lambda/\alpha\rho_k$ |
| $K_{\mathbf u_I}^{\mathbf u_k}$ | $\nu\Delta_kg/\alpha\rho_k$ | $K_P^{\varphi_k}$ | $\dfrac{\rho_k^2(e_{\rho,\varphi})_k-P_I}{\alpha\rho_k(e_{\varphi,\rho})_k}\,\mu$ |
| $K_{\mathbf u}^{\varphi_k}$ | $\dfrac{\lambda}{\alpha\rho_k(e_{\varphi,\rho})_k}$ | $K_T^{\varphi_k}$ | $-\dfrac{\theta}{\alpha\rho_k(e_{\varphi,\rho})_k}$ |
| $K_g^{\varphi_k}$ | $\dfrac{\nu}{\alpha\rho_k(e_{\varphi,\rho})_k}\Big[(e^t_k-e^t_I)+\Big(\rho_k-\dfrac{\rho_k^2}{\rho_I}\Big)(e_{\rho,\varphi})_k+\mathbf u_k\cdot(\mathbf u_I-\mathbf u_k)\Big]$ | | |

Then $K_y^{\Delta_1\varphi}=K_y^{\varphi_1}+K_y^{\varphi_2}$ where **$\varphi$ in (13) is taken as $P$, $T$ or $g$ respectively** — i.e. to build (13) with $\varphi_k=T_k$ as stored variable one still needs $(e_{\rho,P})_k,(e_{P,\rho})_k$, $(e_{\rho,g})_k,(e_{g,\rho})_k$ etc. for each EOS (the paper does not spell out how these derivative sets are computed for SW; presumably through CoolProp derivatives / Jacobian transformations — **[uncertain]**).

### A.1 Stiffened-gas example, Eq. (A.8) ($\varphi_k=P_k$)
With $T_k=\dfrac{P_k+P_{\infty,k}}{\rho_kc_{v,k}(\gamma_k-1)}$ and
$$
g_k(\rho_k,P_k)=\frac{\gamma_k}{\gamma_k-1}\frac{P_k+P_{\infty,k}}{\rho_k}+q_{0,k}-\frac{P_k+P_{\infty,k}}{\rho_k(\gamma_k-1)}\left[\log\!\left(\frac{P_k+P_{\infty,k}}{(\rho_kc_{v,k}(\gamma_k-1))^{\gamma_k}}\right)+\frac{q'_{0,k}}{c_{v,k}}\right],
$$
the linearised pressure ODEs read
$$
\frac{dP_1}{dt}=\bar K_P^{P_1}(P_1-P_2)+\bar K_T^{P_1}\big(T_1(\rho_1,P_1)-T_2(\rho_2,P_2)\big)+\bar K_g^{P_1}\big(g_1(\rho_1,P_1)-g_2(\rho_2,P_2)\big)+\bar K_{\mathbf u}^{P_1}(\mathbf u_1-\mathbf u_I)\cdot\Delta_1\mathbf u, \tag{A.8}
$$
and symmetrically for $P_2$ — still strongly nonlinear even for SG. (I verified the $g_k(\rho,P)$ expression from (5)-(6).)

## 6. Appendix B — linearisation-quality measure (Table B.9)

$D_i$ = linearised coefficients; weights $\mathscr W_i$ chosen by **trial and error** to harmonise error vs number of sub-steps across all tests; depend on **condensation ($g_1>g_2$) vs evaporation ($g_1<g_2$)**. Not claimed optimal.

| $D_i$ | $\mathscr W_i$ condensation | $\mathscr W_i$ evaporation | appears in |
|---|---|---|---|
| $K_P^{\Delta_1P}$ | 1 | $5\cdot10^{-2}$ | (20) |
| $K_T^{\Delta_1P}$ | 1 | $5\cdot10^{-2}$ | (20) |
| $K_g^{\Delta_1P}$ | 4 | $20\cdot10^{-2}$ | (20) |
| $K_T^{\Delta_1T}$ | 1 | $5\cdot10^{-2}$ | (19a) |
| $K_g^{\Delta_1T}$ | 0 | 0 | (19a) |
| $K_T^{\Delta_1g}$ | 1 | $5\cdot10^{-2}$ | (19b) |
| $K_g^{\Delta_1g}$ | 1 | $5\cdot10^{-2}$ | (19b) |
| $K_{\mathbf u}^{P_k}$ | 1 | $5\cdot10^{-2}$ | (20) |
| $K_{\mathbf u}^{T_k}$ | 1 | $5\cdot10^{-2}$ | (19a) |
| $K_{\mathbf u}^{g_k}$ | 1 | $5\cdot10^{-2}$ | (19b) |
| $K_{\mathbf u}^{\mathbf u_k}$ | 10 | 1 | (16), (21d) |
| $K_g^{\alpha_1}$ | 1 | $5\cdot10^{-2}$ | (21b) |
| $K_P^{\varphi_k}$ | 1 | $5\cdot10^{-2}$ | (21d) |
| $K_T^{\varphi_k}$ | 1 | $5\cdot10^{-2}$ | (21d) |
| $K_g^{\varphi_k}$ | 1 | $5\cdot10^{-2}$ | (21d) |

(Row labels read from small print; the $K_{\mathbf u}^{P_k},K_{\mathbf u}^{T_k},K_{\mathbf u}^{g_k}$ rows are the velocity coefficients of the $\varphi=P,T,g$ equations entering $S_{\mathbf u}^{\Delta_1\varphi}$ — **[uncertain]** on exact sub/superscripts.) Note evaporation weights are ~20x looser than condensation.

---

## 7. Verification of the relaxation solver alone (0-D, Sec. 4)

Setup common to all three: relaxation operator (4b) only, high initial disequilibrium, both phases evolve to full mechanical + thermo-chemical equilibrium. **SW EOS via CoolProp** for all. Reference ("exact") = **LSODA** with rel/abs tolerance $2.221\cdot10^{-14}$. Convergence study: $\delta^{\max}_i=2^{-i}$, $i=0,1,2,\dots$; errors ($L_1$) vs number of sub-steps $N_t$. **Constant** relaxation parameters here. Units as printed: $\mu$ [1/(Pa s)], $\lambda$ [kg/(m$^3$ s)], $\theta$ [Pa/(K s)], $\nu$ printed as "kg$^3$/J·m$^3$·s" (**[uncertain]**; dimensional analysis of $\dot m=-\nu\Delta_1g$ gives kg$^2$/(J m$^3$ s)).

### Test 4.1 — R134a pure condensation (Table 1, Figs. 2-3)
| | Phase 1 (vapour, SW [38]) | Phase 2 (liquid, SW [38]) |
|---|---|---|
| $\alpha$ | 0.5 | 0.5 |
| $\rho$ | 75.595 kg/m$^3$ | 1159.329 kg/m$^3$ |
| $u$ | 10 m/s | −5 m/s |
| $T$ | 329.81 K | 310.26 K |

$t_{\rm final}=10^{-4}$ s; $\mu=10^4$, $\lambda=10^9$, $\theta=10^{12}$, $\nu=10^2$. Plots with $\delta^{\max}=2^{-3}$. Initial $P_1\approx15$ bar, $P_2\approx10$ bar (from Fig. 2). Three separate relaxation time scales visible (log time axis $10^{-9}$-$10^{1}\,\mu$s). Adaptive solver time-accurate with **28 sub-steps**; final state on saturation line (~10 bar, ~313 K; $\alpha_1$ 50% -> ~48.3%, $Y_1$ ~6.1% -> ~3.9%). Errors (Fig. 3): adaptive much lower absolute error than LSODA at few sub-steps; adaptive convergence **first order** (velocity almost 2nd in first refinements); LSODA higher order.

### Test 4.2 — MM (hexamethyldisiloxane) pure evaporation (Table 2, Figs. 4-5)
| | Phase 1 (MM vapour [39]) | Phase 2 (MM liquid [39]) |
|---|---|---|
| $\alpha$ | 0.1 | 0.9 |
| $\rho$ | 5.384 kg/m$^3$ | 642.123 kg/m$^3$ |
| $u$ | 0 m/s | 30 m/s |
| $T$ | 384.26 K | 398.99 K |

$t_{\rm final}=10^{-4}$ s; $\mu=10^4$, $\lambda=10^9$, $\theta=10^{12}$, $\nu=10^3$. Time-accurate with **11 sub-steps**; convergence more erratic: after ~20 sub-steps higher than 2nd order, then approaches 1st order. Up to "1000% less error" than LSODA at 100 sub-steps; LSODA higher rate. Final ~2 bar, ~399 K, $\alpha_1$ 10% -> ~10.19%.

### Test 4.3 — D6 mixed evaporation-condensation (Table 3, Figs. 6-8)
| | Phase 1 (D6 vapour [40]) | Phase 2 (D6 liquid [40]) |
|---|---|---|
| $\alpha$ | 0.1 | 0.9 |
| $\rho$ | 11.459 kg/m$^3$ | 707.314 kg/m$^3$ |
| $u$ | 20 m/s | 21 m/s |
| $T$ | 525.01 K | 516.51 K |

$t_{\rm final}=10^{-4}$ s; $\mu=10^4$, $\lambda=10^6$, $\theta=10^{12}$, $\nu=10^5$. **13 sub-steps**; convergence > 2nd order initially then between 1st and 2nd. Fig. 8: $\delta^{\max}=2^0\ldots2^{-7}$ -> $N_t=6,8,11,14,20,30,53,89$; even 6 sub-steps give acceptable (stable, physical) results. $\alpha_1$ goes up (evaporation, peak ~10.02%) then down (condensation) to ~9.99%.

**Observations across 0-D tests:** sub-step adaptation refines mostly in the **chemical-relaxation time-scale range** => Gibbs relaxation is the most nonlinear process. **LSODA could not give a physically acceptable solution with < 100 sub-steps (positivity issues)**; adaptive solver acceptable with as few as 8 (conclusions) and lower absolute error than LSODA below ~300 sub-steps.

---

## 8. Full BN solver tests (Sec. 5)

All with $\delta^{\max}=0.1$, $r^{\max}=10^{-12}$. Hyperbolic: HLLC [31] + MUSCL [32] + Barth-Jespersen [33], explicit [34,35]. Strang splitting with relaxation (4b) by the adaptive solver.

### 5.1a Dodecane shock tube, SG EOS (Table 4, Fig. 9) — from Saurel et al. [15]
SG parameters:
| | vapour dodecane | liquid dodecane |
|---|---|---|
| $\gamma$ | 1.025 | 2.35 |
| $c_v$ [J/(kg K)] | 1956.45 | 1077.7 |
| $P_\infty$ [Pa] | 0 | $4\cdot10^8$ |
| $q_0$ [J/kg] | $-237.547\cdot10^3$ | $-775.269\cdot10^3$ |
| $q'_0$ [J/(kg K)] | $-24.4\cdot10^3$ | 0 |

Initial data: left $x\in[-0.5,0.25)$ m: $\alpha_1=10^{-6}$, $\alpha_2=1-10^{-6}$, $u=0$, $P_1=P_2=10^8$ Pa, $T_1=475$ K, $T_2=687.335$ K; right $x\in(0.25,0.5]$ m: $\alpha_1=1-10^{-6}$, $\alpha_2=10^{-6}$, $u=0$, $P=10^5$ Pa, $T_1=1022.26$ K, $T_2=486.023$ K. (Left liquid at high pressure, right vapour at 1 bar; non-flowing phase temperatures differ.)
Relaxation "sufficiently high to mimic instantaneous": $\mu=10^4$, $\lambda=10^{11}$, $\theta=10^{12}$, $\nu=10^2$. $\Delta x=0.001$ m (1000 cells), CFL = 1 (as printed), $t_{\rm final}=473\,\mu$s.
Result: $Y_1$, $P_m$ (log scale $10^3\to10^1\to$ few bar $\to1$), $\rho_m$, $u_m$ agree well with Pelanti's 6-eq exponential semi-analytical solver [28] ($\theta,\nu\to\infty$). (Evaporation front; mixture velocity up to ~310 m/s.)

### 5.1b Water cavitation (double expansion), SG EOS (Table 5, Fig. 10) — from Pelanti & Shyue [22]
| | vapour water | liquid water |
|---|---|---|
| $\gamma$ | 1.43 | 2.35 |
| $c_v$ | 1040 | 1816 |
| $P_\infty$ | 0 | $10^9$ |
| $q_0$ | $2030\cdot10^3$ | $-1167\cdot10^3$ |
| $q'_0$ | $-23.4\cdot10^3$ | 0 |

Initial: $\alpha_1=0.01$, $\alpha_2=0.99$, $P=10^5$ Pa, $T=354.728$ K everywhere; $u=-2$ m/s for $x\in[-0.5,0)$, $u=+2$ m/s for $x\in(0,0.5]$. $\mu=10^4$, $\lambda=10^{11}$, $\theta=10^{12}$, $\nu=10^4$; $\Delta x=0.001$ m, CFL = 1, $t_{\rm final}=3.2$ ms.
Three models: mechanical only ($\theta=0,\nu=0$), mechanical + thermal ($\theta\to\infty,\nu=0$), mechanical + thermo-chemical ($\theta\to\infty,\nu\to\infty$). **Chemical relaxation activated only if the liquid is metastable** (as in [22]) — i.e. a switch, not just a rate. Great agreement with Pelanti et al. [22]. Without mass transfer, **spurious oscillations in $Y_1$** appear: due to the complete BN model with two velocities and the complete HLLC solver [31] — the two contact discontinuities of liquid and vapour propagate at different velocities -> local composition variations even without mass transfer; bounded, confined near the initial discontinuity, two orders of magnitude smaller than cavitation-induced variation. Without thermo-chemical relaxation the pressure drops below the liquid saturation pressure (~0.05-0.1 bar at centre); with it, the pressure stabilises at saturation (~0.5 bar) while evaporating ($Y_1$ peak ~$2\times10^{-4}$, $\alpha_1$ peak ~0.4).

### 5.2 Water rapid depressurization — Bartak [41] experiment (Table 6, Figs. 11-12)
Long thin high-pressure vessel with liquid water, closed except on the left where an outlet to ambient opens at $t=0$; sudden expansion -> metastable liquid, long-lasting chemical disequilibrium. Goal: show the solver handles **time-dependent closures for $\nu$** and reproduces Pelanti [28].
- Domain: ambient chamber $x\in[-1.7,0)$ m, high-pressure vessel $x\in(0,1.7]$ m, diaphragm at $x=0$, pressure probe at 48 mm right of the diaphragm (vessel length to probe annotated 1.652 m in Fig. 11). **2D mesh $1000\times2$** control volumes (to model the right wall of the vessel rather than 1D).
- SW EOS for water vapour/liquid [44] (IAPWS-95 via CoolProp).
- Initial: ambient: $\alpha_1=1-10^{-3}$, $\alpha_2=10^{-3}$, $u=0$, $P=10^5$ Pa, $T_1=T_2=T^{\rm sat}(P)$; vessel: $\alpha_1=10^{-3}$, $\alpha_2=1-10^{-3}$, $u=0$, $P=1.25\cdot10^7$ Pa, $T_1=T^{\rm sat}(P)$, $T_2=563.13$ K (liquid at 125 bar, 563.13 K).
- $\mu=10^4$, $\lambda=10^{11}$, $\theta=10^{12}$ (the $\mu$ exponent is faint in the PDF — **[uncertain] $10^4$**); $\nu$ varied: $10^2,10^{-1},10^{-2}$ and Pelanti's closure
$$\nu=\alpha_1^{0.6}\Big[\big(P^{\rm sat}-P_2\big)\big/\big(P^{c}-P^{\rm sat}\big)\Big]^{1.76}$$
(as printed; $P^c$ critical pressure; units/prefactor as in [28] — **[uncertain]** whether a dimensional constant is absorbed). $\Delta x=0.0034$ m, CFL = 1, $t_{\rm final}=0.0175$ s.
- Results: with the $\nu(P,P^{\rm sat},\alpha)$ closure + SW, agrees well with experiment and with Pelanti [28] (who used NASG EOS + same closure): pressure undershoot to ~30-35 bar then recovery to ~40-43 bar. $\nu=10^2$ (fastest) almost overlaps HEM/IF97 results of De Lorenzo et al. [43] — equilibrium maintained, pressure follows saturation (~73 bar plateau), no metastability. $\nu=10^{-1},10^{-2}$ don't reach equilibrium nor match pressure drop/recovery. Message: the **closure of $\nu$ is crucial**, more than the EOS accuracy.

### 5.3 2D converging-diverging supersonic nozzle with evaporating MDM droplets (Table 7, Figs. 13-14)
Geometry/conditions from single-phase MDM experiments of Spinelli et al. [45]; inject 5% volume of MDM droplets with diameter $D_p\in\{1,10,50,100,200\}\,\mu$m. Not validation (no data), showcase of fully finite-rate nonlinear closures. All quantities in closures evaluated at $\mathbf w^{*,n+1}$ (post-hyperbolic state, i.e. frozen over the relaxation step).
- Fluids (SW [50]): MDM vapour $k_1=0.0353$ W/(m K), $\eta_1=1.5148\cdot10^{-5}$ Pa s; MDM droplets $k_2=0.06475$ W/(m K), $\eta_2=8.7675\cdot10^{-5}$ Pa s.
- Inlet: $\alpha_1=0.95$, $P_1^{\rm tot}=904388$ Pa, $T_1^{\rm tot}=542.13$ K. Outlet: $P=200000$ Pa. $\Delta x\simeq0.0012$ m, CFL = 1, $t_{\rm final}=15\,\mu$s (as printed; Fig. 13 caption says $t=15$ ms — **[uncertain]** which).

Closures:
$$
\mu=\frac{6\alpha_2}{\pi D_p\rho_2c_2}\quad\text{(Saurel et al. [46], acoustic impedance of liquid)}, \tag{22}
$$
$$
\lambda_{\rm drag}=3\pi\eta_1D_pf_{\rm drag}n_p,\qquad n_p=\frac{6\alpha_2}{\pi D_p^3}\quad\text{(McGrath et al. [47], dense particle cloud)}, \tag{23}
$$
$$
f_{\rm drag}=\begin{cases}8.33\dfrac{\alpha_2}{\alpha_1}+0.0972\,{\rm Re}_R & \alpha_1<0.8\\ f_{\rm base}\,\alpha_1^{-\xi}&\alpha_1\ge0.8\end{cases},\quad
f_{\rm base}=\begin{cases}1+\dfrac{{\rm Re}_R^{0.687}}{6}&{\rm Re}_R<1000\\ 0.0183\,{\rm Re}_R&{\rm Re}_R\ge1000\end{cases},
$$
$$
\xi=3.7-0.65\exp\!\Big(-\tfrac12\big(1.5-\log_{10}{\rm Re}_R\big)^2\Big),\qquad {\rm Re}_R=\frac{\rho_1D_p|\mathbf u_1-\mathbf u_2|}{\eta_1}.
$$
Heat transfer droplet/vapour [48]:
$$
Q=\xi A_I\big(T^{\rm sat}(P_2)-T_1\big)\;\rightarrow\;\theta\simeq\frac{Q}{T_2-T_1}, \tag{24}
$$
$A_I=6\alpha_2/D_p$; heat transfer coefficient (symbol $\xi$ reused) $=Nu\,k_1/D_p$ (printed "$Nu/k_1D_p$", read as $Nu\,k_1/D_p$ — **[uncertain]**), $Nu$ from Ranz-Marshall [49]. Mass transfer from interface heat balance:
$$
\dot m=\frac{Q}{h_2-h_1^{\rm sat}(P_2)}\;\rightarrow\;\nu\simeq\frac{\dot m}{g_2-g_1}, \tag{25}
$$
$h$ specific enthalpy. In this evaporation, $T_2<T^{\rm sat}(P_2)<T_1$ -> $Q<0$ (vapour heats liquid), $\dot m>0$ (liquid -> vapour). (Note the rate-from-flux definitions $\theta=Q/(T_2-T_1)$, $\nu=\dot m/(g_2-g_1)$ can be singular as $\Delta T,\Delta g\to0$; the paper does not discuss regularisation — **[my remark]**.)
Instantaneous comparison: $\mu=10^5$, $\lambda=10^{11}$, $\theta=10^{12}$, $\nu=10^4$.

Results (Fig. 13, liquid volume fraction): pressure and thermal relaxation times from the closures are almost instantaneous for all $D_p$ -> associated disequilibrium negligible. Gibbs relaxation time decreases with the droplet diameter (paper's wording; i.e. it is shorter for smaller droplets): small droplets evaporate before the throat; droplets > 10 µm follow nearly straight trajectories to the outlet (inertia, velocity disequilibrium) and do not follow the expanding vapour in the divergent; droplets > 1 µm tend to cluster on the convergent walls.

**CPU cost (Fig. 14, averaged over 100 time steps, normalised by hyperbolic step):** LSODA run at least restrictive tolerances still robust. Larger droplets/less stiff: LSODA relaxation ~200% of hyperbolic CPU; instantaneous: LSODA ~28400%. Adaptive: ~60% of hyperbolic in less stiff cases (~330% speed-up vs LSODA), roughly constant across cases except instantaneous ~200% (speed-up 14200%).

---

## 9. Conclusions (Sec. 6)
- Efficient, time-accurate relaxation solver for 7-eq BN with simultaneous mechanical + thermo-chemical relaxation, generic EOS; extends Chiocchetti & Müller [30].
- Augmented system of linearised ODEs including three extra ODEs for interphase disequilibria in $P,T,g$.
- 0-D: LSODA needs >100 sub-steps for physical solution; proposed solver acceptable with 8; lower error than LSODA below ~300 sub-steps; convergence between 1st and 2nd order (lower than LSODA) but competitive esp. for stiff problems.
- Speed-up over LSODA 330% (mild stiffness) to 14200% (very stiff).
- Claimed first full non-equilibrium BN with finite-rate thermo-chemical relaxation + simultaneous relaxation under generic EOS. Extension to three-phase (inert third phase) would need further approximations.
- Supplementary `.cpp` with analytical solution of (21); data on Zenodo [51] (zenodo.org/records/12178930, doi 10.5281/zenodo.12178930).
- Funding: EU Horizon Europe MSCA-2021-PF 101066019 (NI2PhORC) (B. Re); NextGenerationEU (G. Sirianni).

---

## 10. Implementation pitfalls, positivity, abnormal cases — and relevance to a 6-eq TV-split IMEX code

From the paper:
1. **Admissibility**: after every sub-step check $0\le\alpha\le1$, $\rho_k>0$, $c_k^2>0$ (not beyond spinodal — important for SW/cubic EOS in metastable regions). Reject and retry with $\Delta t/3$.
2. **Fixed-point divergence** (> $j^{\max}=10$) also triggers rejection. Typically 3-5 iterations.
3. **Weights $\mathscr W_i$ are hand-tuned** and condensation/evaporation dependent; too-strict linearisation bounds over-constrain sub-steps.
4. **Gibbs relaxation is the hardest** (most nonlinear) — adaptation concentrates there.
5. **Time-scale separation assumption** ($\Delta P$ eq. one-way coupled) is a modelling approximation inside the linearised sub-step: if $\theta$ or $\nu$ become faster than $\mu$, the dropped $\bar K_P^{\Delta_1T}\Delta_1P$, $\bar K_P^{\Delta_1g}\Delta_1P$ terms could matter — not analysed in the paper **[my remark]**.
6. **Chemical relaxation switch**: in the cavitation test it is turned on only when the liquid is metastable (following [22]).
7. **Spurious $Y_1$ oscillations** in 7-eq with complete HLLC (two contacts) even without mass transfer — inherent to the two-velocity model, not the relaxation solver.
8. LSODA (generic stiff solver) fails positivity with few steps — motivates structure-exploiting analytical solutions.
9. Interface closures $\rho_I$, $e_I$ (3) are borrowed from equilibrium-reduced models; acknowledged as provisional.
10. Analytical expressions (exponentials of 2x2 systems, time integrals of products with velocity exponentials) are very long; generated symbolically (`.cpp`). For the 2x2 $(\Delta T,\Delta g)$ system, eigenvalues of $\begin{bmatrix}\bar K_T^{\Delta_1T}&\bar K_g^{\Delta_1T}\\ \bar K_T^{\Delta_1g}&\bar K_g^{\Delta_1g}\end{bmatrix}$ may be complex or (nearly) repeated; the paper does not discuss these degenerate branches — an implementer must handle repeated eigenvalues / $K\to0$ limits robustly (e.g. via $\varphi$-functions $(e^{z}-1)/z$) **[my remark]**.

Mapping to diffusionIMEX (6-eq, single velocity, rates $\mu,\theta,\nu$):
- Drop (16) and all $S_{\mathbf u}$, $K_{\mathbf u}$, $\lambda$ terms; $\mathbf u_I=\mathbf u$. $K_g^{\varphi_k}$ loses the kinetic term. 6-eq usually evolves phasic energies rather than a mixture $E$ with two velocities; derivation (A.1)-(A.6) applies verbatim.
- **Limits:** rate -> 0: all $K\to0$, $\tilde{\mathbf w}\equiv\mathbf w^m$, exact frozen scheme (sub-step = whole step). Rate -> ∞: linear exponential solutions tend to the (frozen-coefficient) equilibrium of the linear system; nothing blows up since the solutions are exponentials with negative real-part eigenvalues — but the result is the equilibrium **of the linearisation at $\bar{\mathbf w}$**, so to get the true instantaneous equilibrium you rely on the fixed point making $\bar{\mathbf w}$ consistent and on the linearisation-error control; the initial guess $\Delta t=1/|K_P^{\Delta_1P}|\to0$ and step growth by $\delta^{\max}/\delta$ governs cost. The paper only emulates ∞ with large finite rates (e.g. $\mu=10^4$-$10^5$, $\theta=10^{12}$, $\nu=10^2$-$10^4$); no asymptotic-preserving proof is given.
- Mixture energy conservation is **not** guaranteed exactly by the linearised update (see Sec. 3.1 note); for an IMEX scheme one would likely recompute one phase's energy from $E_m$ after the relaxation step.
- Choice of stored thermodynamic variable: $T_k$ chosen for SW ($(\rho,T)$-explicit); with SG/NASG, $P_k$ or $T_k$ are both cheap.
