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

