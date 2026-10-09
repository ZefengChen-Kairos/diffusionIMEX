# Reading note: Semi-implicit method for two-phase flows in pressure disequilibrium

- **Title:** *Semi-implicit method for two-phase flows in pressure disequilibrium*
- **Authors:** Solène Schropff (corresponding), Fabien Petitpas, Eric Daniel — Aix Marseille Univ, CNRS, IUSTI, Marseille, France
- **Venue:** *Journal of Computational Physics* **547** (2026) 114545. DOI: [10.1016/j.jcp.2025.114545](https://doi.org/10.1016/j.jcp.2025.114545)
- **Dates:** received 17 Mar 2025, revised 15 Nov 2025, accepted 22 Nov 2025, available online 28 Nov 2025. Open access (CC BY).
- **Keywords:** two-phase flows, compressible flows, low-Mach flows, implicit-explicit, IMEX scheme, pressure-based.

> Notation in this note: $k$ = phase index ($k=1,2$; general $\Pi$ phases), $i$ = cell index, $n$ = time level,
> $(n+\star)$ = intermediate level after the implicit (pressure) operator, $pred$ = level after the explicit
> (advection) operator, $r$ = Picard iteration index. Items I could not read with certainty are marked **[uncertain]**.

## 0. One-paragraph orientation (how this relates to our TV-split plan)

The paper builds a **1D collocated, first-order-in-time, operator-split IMEX scheme** for the
6-equation (velocity-equilibrium, pressure-disequilibrium) model of Saurel–Petitpas–Berry (2009) [their ref 19],
with the Pelanti-style/Saurel energy correction of [20]. The flux splitting is a **Toro–Vázquez-Cendón (TV) /
Zha–Bilgen (ZB)–type** splitting extended to the 6-equation model: a *pressure subsystem* (acoustic, eigenvalues
$0,\pm A$, independent of $u$) solved implicitly, and an *advection subsystem* (all eigenvalues $u$) solved
explicitly. So despite the user's guess, **it *is* a TV/ZB-style flux-split scheme** — the splitting is applied as a
Lie (sequential) operator splitting $[C][\mathcal E][\mathcal I]$, not as a Runge–Kutta IMEX with a Butcher tableau.
The implicit part is linearised by **Picard iteration**, giving per Picard step a **heptadiagonal linear system of size
$2N$** for the two phase pressures, solved by banded Gauss–Seidel. Pressure **relaxation is done afterwards,
independently** (as in [19,20]), not inside the IMEX step.

## 1. Introduction / motivation (pp. 1–3)

- Explicit FV schemes for the 6-eq. model are constrained by the acoustic CFL with the *frozen* mixture sound
  speed $c_f$; a small amount of liquid raises $c_f$ towards the liquid value, so $\Delta t$ collapses in
  weakly-compressible/low-Mach flows; also the explicit compressible discretisation is not asymptotic-preserving
  in the low-Mach limit [21–24].
- Fully implicit: very non-linear, costly, and dissipative for slow waves. Semi-implicit methods [27–39] replace the
  acoustic restriction by a convective one.
- Pressure-based methods (Harlow–Amsden [27,28], SIMPLE extensions [36], all-speed schemes [44,50,51,52]).
  Conservative semi-implicit = split flux vector hyperbolically into an implicit part and an explicit part;
  flux-splitting literature for Euler [53–56], Baer–Nunziato [57], Kapila [58]; semi-implicit conservative schemes [39,59–67].
- Contribution: reformulate the **Zha–Bilgen splitting [54] as proposed by Toro & Vázquez [55]** for the velocity-equilibrium
  (6-eq) model — "no flux-splitting has been proposed for this model so far". Pressure subsystem implicit (Picard), semi-discrete
  collocated 1D formulation à la [63,64] for a compact stencil; advection subsystem explicit (material CFL);
  internal energies corrected as in [20] for total-energy conservation. General EOS formulation.

## 2. Two-phase flow model and mathematical properties

### 2.1 Flow equations (single velocity, two pressures) — Eq. (1)

$$
\begin{cases}
\partial_t \alpha_1 + \mathbf u\cdot\nabla\alpha_1 = \mu(p_1-p_2)\\
\partial_t \alpha_2 + \mathbf u\cdot\nabla\alpha_2 = -\mu(p_1-p_2)\\
\partial_t(\alpha_1\rho_1) + \nabla\cdot(\alpha_1\rho_1\mathbf u) = 0\\
\partial_t(\alpha_2\rho_2) + \nabla\cdot(\alpha_2\rho_2\mathbf u) = 0\\
\partial_t(\alpha_1\rho_1 e_1) + \nabla\cdot(\alpha_1\rho_1 e_1\mathbf u) + \alpha_1 p_1\nabla\cdot\mathbf u = -\mu p_I(p_1-p_2)\\
\partial_t(\alpha_2\rho_2 e_2) + \nabla\cdot(\alpha_2\rho_2 e_2\mathbf u) + \alpha_2 p_2\nabla\cdot\mathbf u = \mu p_I(p_1-p_2)\\
\partial_t(\rho\mathbf u) + \nabla\cdot(\rho\mathbf u\otimes\mathbf u) + \nabla p = \mathbf 0
\end{cases}\tag{1}
$$

- $\alpha_k$ volume fraction, $\rho_k$ phase density, $p_k$ phase pressure, $e_k$ phase specific internal energy;
  saturation $\sum_k\alpha_k=1$.
- Mixture: $\rho=\sum_k\alpha_k\rho_k$, $p=\sum_k\alpha_k p_k$, $e=\sum_k Y_k e_k$, $Y_k=\alpha_k\rho_k/\rho$,
  $E=e+\tfrac12\mathbf u^2$. Common velocity $\mathbf u$.
- $\mu$ = pressure relaxation coefficient, $0\le\mu\le+\infty$ (chosen per problem, [20]). $p_I$ interfacial pressure.
- Note: only mechanical (pressure) relaxation appears; **no thermal or chemical relaxation** in this paper.

Total energy (redundant, conservative), obtained by combining the phase internal energies and momentum:
$$\partial_t(\rho E) + \nabla\cdot\big((\rho E + p)\mathbf u\big) = 0. \tag{2}$$

**EOS.** Each phase has its own convex EOS $e_k=e_k(p_k,\rho_k)$. General sound speed
$$c_k^2 = \frac{p_k}{\rho_k^2 e_{p,k}} - \frac{e_{\rho,k}}{e_{p,k}},\tag{3}$$
$$e_{p,k} = \left.\frac{\partial e_k}{\partial p_k}\right)_{\rho_k},\qquad e_{\rho,k} = \left.\frac{\partial e_k}{\partial \rho_k}\right)_{p_k}.\tag{4}$$
Stiffened gas (SG):
$$e_k(p_k,\rho_k) = \frac{p_k + \gamma_k p_{k,\infty}}{(\gamma_k-1)\rho_k} + e_{k,ref},\tag{5}$$
$$c_k^2 = \frac{\gamma_k(p_k+p_{k,\infty})}{\rho_k}.\tag{6}$$
(Paper mentions extension to NASG [68] / Mie–Grüneisen [69,70]; parameters from reference curves [71].)

Entropy equations:
$$
\begin{cases}
\alpha_1\rho_1 T_1 \dfrac{d s_1}{dt} = (p_1 - p_I)\mu(p_1-p_2)\\[4pt]
\alpha_2\rho_2 T_2 \dfrac{d s_2}{dt} = -(p_2 - p_I)\mu(p_1-p_2)
\end{cases}\tag{7}
$$
Interfacial pressure: asymptotic limit of the 7-eq. model's [72]; to ensure entropy production, a possible choice is
$$p_I = \frac{Z_2 p_1 + Z_1 p_2}{Z_1+Z_2},\qquad Z_k=\rho_k c_k$$ [73,74].

Hyperbolicity: waves $u,\ u-c_f,\ u+c_f$ ([75] for analysis), with **frozen** sound speed
$$c_f^2 = \sum_k Y_k c_k^2,\tag{8}$$
monotone in volume fraction.

### 2.2 Splitting into convective and acoustic parts (TV/ZB-style)

Relaxation terms omitted. Primitive vector
$$\mathbf W = (\alpha_1,\alpha_2,\rho_1,\rho_2,p_1,p_2,\mathbf u)^T.\tag{9}$$

**Pressure subsystem** (10):
$$
\begin{cases}
\partial_t\alpha_1 = 0,\quad \partial_t\alpha_2 = 0\\
\partial_t(\alpha_1\rho_1)=0,\quad \partial_t(\alpha_2\rho_2)=0\\
\partial_t(\alpha_1\rho_1e_1) + \alpha_1p_1\nabla\cdot\mathbf u = 0\\
\partial_t(\alpha_2\rho_2e_2) + \alpha_2p_2\nabla\cdot\mathbf u = 0\\
\partial_t(\rho\mathbf u) + \nabla p = \mathbf 0
\end{cases}\tag{10}
$$
Primitive form $\partial_t\mathbf W + \mathbf Q_p(\mathbf W)\nabla\cdot\mathbf W=0$ with (rows ordered as $\mathbf W$)
$$
\mathbf Q_p(\mathbf W)=\begin{pmatrix}
0&0&0&0&0&0&0\\
0&0&0&0&0&0&0\\
0&0&0&0&0&0&0\\
0&0&0&0&0&0&0\\
0&0&0&0&0&0&p_1\Gamma_1\\
0&0&0&0&0&0&p_2\Gamma_2\\
\frac{p_1}{\rho}&\frac{p_2}{\rho}&0&0&\frac{\alpha_1}{\rho}&\frac{\alpha_2}{\rho}&0
\end{pmatrix}.\tag{11}
$$
Eigenvalues $\lambda=(0,0,0,0,0,-A,A)^T$ (12) with
$$A^2 = \frac1\rho\sum_k\alpha_kp_k\Gamma_k,\tag{13}$$
$\Gamma_k=(\rho_k e_{p,k})^{-1}$ the Grüneisen parameter [76]. Hyperbolic if $p_k\Gamma_k>0$. For SG:
$$A^2 = \sum_k\alpha_kp_k(\gamma_k-1)\;=\;c_f^2 - \frac{p+\sum_k\alpha_k\gamma_kp_{k,\infty}}{\rho}\tag{14}$$
**[uncertain]**: the middle expression as printed is $\sum_k\alpha_kp_k(\gamma_k-1)$ without a $1/\rho$; from (13) with SG,
$\Gamma_k=\gamma_k-1$, so it should be $\frac1\rho\sum_k\alpha_kp_k(\gamma_k-1)$. Treat (13) as authoritative. Ideal gas: $A^2=c_f^2-p/\rho$.
The waves $\pm A$ do not depend on $u$ ⇒ the pressure subsystem wave pattern is always subsonic.
*Pitfall:* $A^2$ can become small or negative if $p_k\le 0$ (SG liquids under tension) — hyperbolicity of (10) needs $p_k\Gamma_k>0$.

Mixture total energy of the pressure subsystem:
$$\partial_t(\rho E) + \nabla\cdot(p\mathbf u) = 0.\tag{15}$$

**Advection subsystem** (16):
$$
\begin{cases}
\partial_t\alpha_1 + \nabla\cdot(\alpha_1\mathbf u) - \alpha_1\nabla\cdot\mathbf u = 0\\
\partial_t\alpha_2 + \nabla\cdot(\alpha_2\mathbf u) - \alpha_2\nabla\cdot\mathbf u = 0\\
\partial_t(\alpha_1\rho_1) + \nabla\cdot(\alpha_1\rho_1\mathbf u)=0\\
\partial_t(\alpha_2\rho_2) + \nabla\cdot(\alpha_2\rho_2\mathbf u)=0\\
\partial_t(\alpha_1\rho_1e_1) + \nabla\cdot(\alpha_1\rho_1e_1\mathbf u)=0\\
\partial_t(\alpha_2\rho_2e_2) + \nabla\cdot(\alpha_2\rho_2e_2\mathbf u)=0\\
\partial_t(\rho\mathbf u) + \nabla\cdot(\rho\mathbf u\otimes\mathbf u)=\mathbf 0
\end{cases}\tag{16}
$$
(10)+(16) = (1) without relaxation; $\alpha_k$ equation written in the equivalent conservative-plus-correction form.
Primitive matrix
$$
\mathbf Q_a(\mathbf W)=\begin{pmatrix}
u&0&0&0&0&0&0\\
0&u&0&0&0&0&0\\
0&0&u&0&0&0&\rho_1\\
0&0&0&u&0&0&\rho_2\\
0&0&0&0&u&0&-e_{\rho,1}\rho_1^2\Gamma_1\\
0&0&0&0&0&u&-e_{\rho,2}\rho_2^2\Gamma_2\\
0&0&0&0&0&0&u
\end{pmatrix},\tag{17}
$$
eigenvalues $\lambda=(u,u,u,u,u,u,u)^T$ (18) — weakly hyperbolic (single repeated eigenvalue; not
diagonalisable in general) **[note: my remark, not stated in the paper]**.
Mixture total energy of the advection subsystem:
$$\partial_t(\rho E) + \nabla\cdot(\rho E\mathbf u)=0.\tag{19}$$
(15)+(19) = (2).

## 3. Classical explicit FV scheme (reference; method of [20])

Hyperbolic part of (1) for $k$ phases, augmented with total energy:
$$
\begin{cases}
\partial_t\alpha_k + \nabla\cdot(\alpha_k\mathbf u) - \alpha_k\nabla\cdot\mathbf u = 0\\
\partial_t(\alpha_k\rho_k) + \nabla\cdot(\alpha_k\rho_k\mathbf u)=0\\
\partial_t(\alpha_k\rho_ke_k) + \nabla\cdot(\alpha_k\rho_ke_k\mathbf u) + \alpha_kp_k\nabla\cdot\mathbf u = 0\\
\partial_t(\rho\mathbf u) + \nabla\cdot(\rho\mathbf u\otimes\mathbf u) + \nabla p = \mathbf 0\\
\partial_t(\rho E) + \nabla\cdot((\rho E+p)\mathbf u) = 0
\end{cases}\tag{20}
$$
Compact: $\partial_t\mathbf U + \nabla\cdot\mathbf F(\mathbf U) + \mathbf H(\mathbf U)\nabla\cdot\mathbf u = \mathbf 0$ (21) with
$$
\mathbf U=\begin{pmatrix}\alpha_k\\\alpha_k\rho_k\\\alpha_k\rho_ke_k\\\rho\mathbf u\\\rho E\end{pmatrix},\quad
\mathbf F=\begin{pmatrix}\alpha_k\mathbf u\\\alpha_k\rho_k\mathbf u\\\alpha_k\rho_ke_k\mathbf u\\\rho\mathbf u\otimes\mathbf u+p\mathbf I\\(\rho E+p)\mathbf u\end{pmatrix},\quad
\mathbf H=\begin{pmatrix}-\alpha_k\\0\\\alpha_kp_k\\\mathbf 0\\0\end{pmatrix}.\tag{22}
$$
**Prediction:**
$$\mathbf U_i^{pred} = \mathbf U_i^n - \frac{\Delta t}{\Delta x}\left(\mathbf F^*_{i+\frac12}-\mathbf F^*_{i-\frac12}\right) - \frac{\Delta t}{\Delta x}\mathbf H_i^n\left(u^*_{i+\frac12}-u^*_{i-\frac12}\right),\tag{23}$$
fluxes and $u^*$ from HLLC [19,77]. CFL: $\Delta t^* = t^{n+1}-t^n = \Delta x/|\sigma|$ (24), $\sigma=\max|u\pm c_f|$
(the CFL number multiplies this in practice — see test cases).

**Energy correction** (of [20]): assume $\alpha_k,\ \alpha_k\rho_k,\ \rho u,\ \rho E$ at $n+1$ equal their $pred$ values; then
$$(\alpha_k\rho_ke_k)^{n+1} = (\alpha_k\rho_ke_k)^{pred} + \alpha_k^n\,\delta_\varepsilon,\tag{25}$$
$$\delta_\varepsilon = \delta_{\rho E} - \delta_{\frac12\rho u^2} - \sum_k\delta^{pred}_{\alpha_k\rho_ke_k},\tag{26a}$$
$$\delta_{\rho E} = (\rho E)^{n+1}-(\rho E)^n,\tag{26b}$$
$$\delta_{\frac12\rho u^2} = \frac12\left(\frac{((\rho u)^{n+1})^2}{\rho^{n+1}} - \frac{((\rho u)^n)^2}{\rho^n}\right),\tag{26c}$$
$$\delta^{pred}_{\alpha_k\rho_ke_k} = (\alpha_k\rho_ke_k)^{pred}-(\alpha_k\rho_ke_k)^n.\tag{26d}$$
So the total internal-energy defect is redistributed with weight $\alpha_k^n$ (sum of weights = 1), making
$\sum_k(\alpha_k\rho_ke_k)^{n+1} = (\rho E)^{n+1} - \tfrac12((\rho u)^{n+1})^2/\rho^{n+1}$ exactly
(given $\sum_k(\alpha_k\rho_ke_k)^n = (\rho E)^n-\tfrac12(\rho u)^{n,2}/\rho^n$ holds at $n$).

## 4. Semi-implicit (IMEX) scheme

Split system with total energy:
$$\partial_t\mathbf U + \nabla\cdot(\mathbf P(\mathbf U)+\mathbf A(\mathbf U)) + (\mathbf H_p(\mathbf U)+\mathbf H_a(\mathbf U))\nabla\cdot\mathbf u=\mathbf 0,\tag{27}$$
$$
\mathbf U=\begin{pmatrix}\alpha_k\\\alpha_k\rho_k\\\alpha_k\rho_ke_k\\\rho\mathbf u\\\rho E\end{pmatrix},\
\mathbf P=\begin{pmatrix}0\\0\\0\\p\mathbf I\\p\mathbf u\end{pmatrix},\
\mathbf A=\begin{pmatrix}\alpha_k\mathbf u\\\alpha_k\rho_k\mathbf u\\\alpha_k\rho_ke_k\mathbf u\\\rho\mathbf u\otimes\mathbf u\\\rho E\mathbf u\end{pmatrix},\
\mathbf H_p=\begin{pmatrix}0\\0\\\alpha_kp_k\\\mathbf 0\\0\end{pmatrix},\
\mathbf H_a=\begin{pmatrix}-\alpha_k\\0\\0\\\mathbf 0\\0\end{pmatrix}.\tag{28}
$$
An explicit flux-split (non-IMEX) variant is in Appendix B: its wave speed $|A|$ (14) does not depend on $u$, but still acoustic.

**Operator splitting (Lie/Godunov, first order in time):**
$$\mathbf U^{n+1} = [C][\mathcal E][\mathcal I](\mathbf U^n),\tag{29}$$
applied right-to-left over one $\Delta t$:
- implicit (pressure) operator $\mathbf U^{n+\star} = [\mathcal I](\mathbf U^n)$, Eq. (37);
- explicit (advection) operator $\mathbf U^{pred} = [\mathcal E](\mathbf U^{n+\star})$, Eq. (46);
- energy correction $\mathbf U^{n+1}=[C](\mathbf U^{pred})$, Eq. (49).
No Butcher tableau: it is not a multi-stage RK-IMEX; time accuracy is first order (Lie splitting, backward-Euler-like
pressure step, forward-Euler advection step — but advection uses MUSCL-Hancock, which is second order on its own).

### 4.1 Pressure subsystem: implicit operator $[\mathcal I]$

Compact form: $\partial_t\mathbf U + \nabla\cdot\mathbf P(\mathbf U) + \mathbf H_p(\mathbf U)\nabla\cdot\mathbf u = \mathbf 0$ (30).
Idea: $\alpha_k$, $\alpha_k\rho_k$ trivial; combine energy + momentum + EOS into a pressure evolution equation per phase, solve implicitly; iterate.

#### 4.1.1 Pressure equations
Primitive form of (10):
$$
\begin{cases}
\partial_t\alpha_k = 0\\
\partial_t\rho_k = 0\\
\partial_t p_k + p_k\Gamma_k\nabla\cdot\mathbf u = 0\\
\partial_t\mathbf u + \frac1\rho\sum_k(p_k\nabla\alpha_k + \alpha_k\nabla p_k) = \mathbf 0
\end{cases}\tag{31}
$$
(The pressure equation coefficient is $p_k\Gamma_k$ — **not** $\rho_kc_k^2$: in the pressure subsystem the $e_k$
equation lacks the convective part, so $D_tp_k = -p_k\Gamma_k\nabla\cdot u$; for SG $p_k\Gamma_k=(\gamma_k-1)p_k$. Note:
for SG this coefficient does not include $p_\infty$ — consistent with (14).)

Time-implicit discretisation between $n$ and $n+\star$:
$$
\begin{cases}
\alpha^{n+\star}_{k,i} = \alpha^n_{k,i}\\
\rho^{n+\star}_{k,i} = \rho^n_{k,i}\\
p^{n+\star}_{k,i} = p^n_{k,i} - \Delta t\,(p_{k,i}\Gamma_{k,i})^{n+\star}\,\nabla\cdot\mathbf u_i^{n+\star}\\
\mathbf u_i^{n+\star} = \mathbf u_i^n - \Delta t\,\dfrac{1}{\rho_i^{n+\star}}\sum_k\left(p_{k,i}\nabla\alpha_{k,i}+\alpha_{k,i}\nabla p_{k,i}\right)^{n+\star}
\end{cases}\tag{32}
$$
Non-linear because $\Gamma^{n+\star}_k=\Gamma(e^{n+\star}_{p,k})=\Gamma(p^{n+\star}_k,\rho^n_k)$ and $u^{n+\star}=u(p^{n+\star})$.
**Picard iteration** [61,63,79]: superscript $r$ = known at current Picard iterate, $r+1$ = unknown. Using $\alpha^{n+\star}=\alpha^n$, $\rho^{n+\star}=\rho^n$:
$$
\begin{cases}
p^{n+\star,r+1}_{k,i} = p^n_{k,i} - \Delta t\,(p_{k,i}\Gamma_{k,i})^{n+\star,r}\,\nabla\cdot\mathbf u_i^{n+\star,r+1}\\
\mathbf u_i^{n+\star,r+1} = \mathbf u_i^n - \Delta t\,\dfrac{1}{\rho_i^n}\sum_k\left(p^{n+\star,r+1}_{k,i}\nabla\alpha^n_{k,i}+\alpha^n_{k,i}\nabla p^{n+\star,r+1}_{k,i}\right)
\end{cases}\tag{33}
$$
i.e. the coefficient $p_k\Gamma_k$ is lagged (Picard), everything else linear in $p^{r+1}$.

Substitute velocity into the pressure equation **before** spatially discretising $\nabla\cdot u$ (compact stencil, [64]):
$$p^{n+\star,r+1}_{k,i} - (\Delta t)^2\nu^{n+\star,r}_{k,i}\sum_{l=1}^{\Pi}\left(\eta^n_{l,i}p^{n+\star,r+1}_{l,i-1} + \zeta^n_{l,i}p^{n+\star,r+1}_{l,i} + \beta^n_{l,i}p^{n+\star,r+1}_{l,i+1}\right) = p^n_{k,i} - \Delta t\,\nu^{n+\star,r}_{k,i}\nabla\cdot\mathbf u^n_i.\tag{34}$$
Coefficients $\nu,\eta,\zeta,\beta$ in Appendix C (presumably $\nu_k = p_k\Gamma_k$ — to verify against App. C).

Two-phase ($\Pi=2$):
$$
(k=1):\ a^{n+\star,r}_{1,i}p^{n+\star,r+1}_{1,i-1} + b^{n+\star,r}_{1,i}p^{n+\star,r+1}_{2,i-1} + a^{n+\star,r}_{2,i}p^{n+\star,r+1}_{1,i}
+ b^{n+\star,r}_{2,i}p^{n+\star,r+1}_{2,i} + a^{n+\star,r}_{3,i}p^{n+\star,r+1}_{1,i+1} + b^{n+\star,r}_{3,i}p^{n+\star,r+1}_{2,i+1} = c^{n+\star,r}_i,\tag{35a}
$$
$$
(k=2):\ m^{n+\star,r}_{1,i}p^{n+\star,r+1}_{1,i-1} + q^{n+\star,r}_{1,i}p^{n+\star,r+1}_{2,i-1} + m^{n+\star,r}_{2,i}p^{n+\star,r+1}_{1,i}
+ q^{n+\star,r}_{2,i}p^{n+\star,r+1}_{2,i} + m^{n+\star,r}_{3,i}p^{n+\star,r+1}_{1,i+1} + q^{n+\star,r}_{3,i}p^{n+\star,r+1}_{2,i+1} = g^{n+\star,r}_i.\tag{35b}
$$
Linear system (heptadiagonal, size $2N$, unknowns interleaved $(p_{1,i},p_{2,i})$):
$$\mathbf A^{n+\star,r}\,\mathbf p^{n+\star,r+1} = \mathbf B^{n+\star,r},\tag{36}$$
solved with **Gauss–Seidel adapted to banded matrices**.

#### 4.1.2 Implicit numerical scheme (FV form of the pressure subsystem)
$$\mathbf U_i^{n+\star} = \mathbf U_i^n - \frac{\Delta t}{\Delta x}\left(\mathbf P^{*,n+\star}_{i+\frac12}-\mathbf P^{*,n+\star}_{i-\frac12}\right) - \frac{\Delta t}{\Delta x}\mathbf H^{n+\star}_{p,i}\left(u^{*,n+\star}_{i+\frac12}-u^{*,n+\star}_{i-\frac12}\right).\tag{37}$$

**Picard loop:**
1. Init ($r=1$): $p^{n+\star,1}_k=p^n_k$, $e^{n+\star,1}_{p,k}=e^n_{p,k}$. If $\mathbf p^n$ already solves (36), no update needed.
2. Solve (36) for $p^{n+\star,r+1}_k$. Update momentum:
$$(\rho u)^{n+\star,r+1}_i = (\rho u)^n_i - \frac{\Delta t}{\Delta x}\,\frac{\sum_k\alpha^n_{k,i+1}p^{n+\star,r+1}_{k,i+1} - \sum_k\alpha^n_{k,i-1}p^{n+\star,r+1}_{k,i-1}}{2}.\tag{38}$$
   (Centered collocated pressure gradient — i.e. interface pressure $p_{i+1/2}=(p_i+p_{i+1})/2$ with mixture $p=\sum\alpha_kp_k$.)
3. Update $e^{n+\star,r+1}_{p,k}=e_{p,k}(p^{n+\star,r+1}_k,\rho^n_k)$ via EOS (5). Next Picard iterate.
4. Stop when the standard deviation between pressures of two consecutive Picard iterates is below a tolerance, or max iterations reached.
   Then $p^{n+\star} = \sum_k\alpha^n_kp^{n+\star}_k$.

**After convergence — energies and conservative updates:**
$$(\alpha_k\rho_ke_k)^{n+\star}_i = (\alpha_k\rho_ke_k)^n_i - \frac{\Delta t}{\Delta x}\alpha^n_{k,i}p^{n+\star}_{k,i}\left(u^{*,n+\star}_{i+\frac12}-u^{*,n+\star}_{i-\frac12}\right),\tag{39}$$
$$u^{*,n+\star}_{i+\frac12} = \frac{(\rho u)^{n+\star}_i + (\rho u)^{n+\star}_{i+1}}{\rho^n_i+\rho^n_{i+1}},\tag{40}$$
$$(\rho E)^{n+\star}_i = (\rho E)^n_i - \frac{\Delta t}{\Delta x}\left((pu)^{n+\star}_{i+\frac12}-(pu)^{n+\star}_{i-\frac12}\right),\tag{41}$$
$$p^{n+\star}_i = \sum_k\alpha^n_{k,i}p^{n+\star}_{k,i},\tag{42}$$
$$(pu)^{n+\star}_{i+\frac12} = p^{n+\star}_{i+\frac12}u^{n+\star}_{i+\frac12},\quad p^{n+\star}_{i+\frac12}=\frac{p^{n+\star}_i+p^{n+\star}_{i+1}}2,\tag{43}$$
with $u^{n+\star}_{i+1/2}$ as in (40). This preserves $p$ and $u$ across a contact.

**Checkerboard:** collocated $p,u$ — not observed in their 1D tests; authors note it may appear in multi-D [22,78,83]
and that pressure relaxation may mitigate it [15,84]. *(Pitfall for a 2D extension: centred $\pm1$ pressure gradient
(38) with face-averaged velocity (40) is the classic odd–even decoupling configuration; Rhie–Chow-type stabilisation may be needed.)*

**Summary of $[\mathcal I]$:**
$$
\mathbf U^{n+\star}_i=\begin{cases}
\alpha^{n+\star}_{k,i} = \alpha^n_{k,i}\\
(\alpha_k\rho_k)^{n+\star}_i = (\alpha_k\rho_k)^n_i\\
(\alpha_k\rho_ke_k)^{n+\star}_i = (\alpha_k\rho_ke_k)^n_i - \frac{\Delta t}{\Delta x}\alpha^n_{k,i}p^{n+\star}_{k,i}\left(u^{*,n+\star}_{i+\frac12}-u^{*,n+\star}_{i-\frac12}\right)\\
(\rho u)^{n+\star}_i = (\rho u)^n_i - \frac{\Delta t}{\Delta x}\frac{\sum_k\alpha^n_{k,i+1}p^{n+\star}_{k,i+1}-\sum_k\alpha^n_{k,i-1}p^{n+\star}_{k,i-1}}{2}\\
(\rho E)^{n+\star}_i = (\rho E)^n_i - \frac{\Delta t}{\Delta x}\left((pu)^{n+\star}_{i+\frac12}-(pu)^{n+\star}_{i-\frac12}\right)
\end{cases}\tag{44}
$$
Note: phase internal energy (39) uses the **implicit** phase pressure $p^{n+\star}_k$ times the *new* face-velocity divergence — this is a
backward-Euler work term. The phase pressure that one would recompute from $(\alpha_k\rho_ke_k)^{n+\star}$ via the EOS is not
exactly the one solved for in (36) (the Picard pressure equation uses $\nabla\cdot u$ discretised from (34)/App. C); any mismatch is
absorbed by the later energy correction $[C]$ **[my interpretation]**.

### 4.2 Advection subsystem: explicit operator $[\mathcal E]$
Compact: $\partial_t\mathbf U + \nabla\cdot\mathbf A(\mathbf U) + \mathbf H_a(\mathbf U)\nabla\cdot\mathbf u = \mathbf 0$ (45).
$$\mathbf U^{pred}_i = \mathbf U^{n+\star}_i - \frac{\Delta t}{\Delta x}\left(\mathbf A^{*,n+\star}_{i+\frac12}-\mathbf A^{*,n+\star}_{i-\frac12}\right) - \frac{\Delta t}{\Delta x}\mathbf H^{n+\star}_{a,i}\left(u^{*,n+\star}_{i+\frac12}-u^{*,n+\star}_{i-\frac12}\right).\tag{46}$$
**Rusanov flux:**
$$\mathbf A^{*,n+\star}_{i+\frac12} = \frac12\left(\mathbf A(\mathbf U^{n+\star,+}_{i+\frac12})+\mathbf A(\mathbf U^{n+\star,-}_{i+\frac12})\right) - \frac12|\sigma|\left(\mathbf U^{n+\star,+}_{i+\frac12}-\mathbf U^{n+\star,-}_{i+\frac12}\right),\tag{47}$$
$\mathbf U^{\mp}_{i+1/2}$ = left/right extrapolated states at $x_{i+1/2}$ (first order: $\mathbf U^-=\mathbf U_i$, $\mathbf U^+=\mathbf U_{i+1}$).
Because Rusanov is dissipative, they use **second-order MUSCL–Hancock TVD** [39,78] (App. D).
Intercell velocity for the non-conservative term: $u^{*,n+\star}_{i+1/2} = (u^{n+\star,-}_{i+1/2}+u^{n+\star,+}_{i+1/2})/2$,
$u^{\pm}=(\rho u)^\pm/\rho^\pm$.
Wave-speed magnitude:
$$\sigma = D\times\max\left(|u^{n+\star,-}_{i+\frac12}|,\ |u^{n+\star,+}_{i+\frac12}|\right),\tag{48}$$
$D\ge1$ an extra-dissipation coefficient to avoid spurious oscillations / keep stability across operators.
**Time step:** the same $\Delta t$ is used in $[\mathcal I]$ and $[\mathcal E]$; the next $\Delta t$ is computed from (24) with $\sigma$ from (48)
(material speed only), i.e. $\Delta t = \text{CFL}\cdot\Delta x/\max|\sigma|$ **[the CFL factor is applied in the test cases]**.
*Pitfall:* if $u\equiv0$ initially, $\sigma=0$ ⇒ $\Delta t=\infty$; need a fallback (see test cases for how they set $\Delta t$).

### 4.3 Energy correction operator $[C]$
$$
\mathbf U^{n+1}_i=\begin{cases}
(\alpha_k)^{n+1}_i = (\alpha_k)^{pred}_i\\
(\alpha_k\rho_k)^{n+1}_i = (\alpha_k\rho_k)^{pred}_i\\
(\rho u)^{n+1}_i = (\rho u)^{pred}_i\\
(\rho E)^{n+1}_i = (\rho E)^{pred}_i\\
(\alpha_k\rho_ke_k)^{n+1}_i = (\alpha_k\rho_ke_k)^{pred}_i + \alpha^n_{k,i}\delta_\varepsilon
\end{cases}\tag{49}
$$
$\delta_\varepsilon$ from (26), with all "pred" quantities from the semi-implicit prediction (46) (i.e. after both $[\mathcal I]$ and $[\mathcal E]$; $n$ quantities from start of step).

### 4.4 Summary
Fig. 1: functional diagram of the IMEX step (29). **Pressure relaxation is performed after the IMEX scheme, fully independently**,
with the methods of [19,20] (i.e. stiff instantaneous relaxation $\mu\to\infty$ solved as an algebraic problem per cell).


**Fig. 1 (functional diagram), read literally:**
1. Start from $\mathbf U^n$. Test "is $p^n$ solution?" (of (36)). If yes: $\mathbf U^{n+\star}\leftarrow\mathbf U^n$ (skip implicit stage).
2. Else $\mathbf U^r\leftarrow\mathbf U^n$; loop: build $A(\mathbf U^r)$, $B(\mathbf U^r)$; solve $A p^{r+1}=B$; compute $\mathbf U^{r+1}(p^{r+1})$;
   test $\|p^{r+1}-p^r\|<\varepsilon$ with **$\varepsilon=10^{-6}$**; if no, $\mathbf U^r\leftarrow\mathbf U^{r+1}$ and repeat; if yes $\mathbf U^{n+\star}\leftarrow\mathbf U^{r+1}$.
   (Whether the norm is absolute or relative is not specified; text says "standard deviation between the pressure of two consecutive Picard iterations" **[uncertain: likely relative/normalised; with $p\sim10^5$–$10^{10}$ Pa an absolute $10^{-6}$ would be very strict]**.)
3. $[\mathcal E]$: compute $\mathbf A^*_{i+1/2}(\mathbf U^{n+\star})$, $u^*_{i+1/2}(\mathbf U^{n+\star})$, $\sigma(\mathbf U^{n+\star})$, and the **next** $\Delta t^*(\sigma)$ → $\mathbf U^{pred}$.
4. $[C]$: $\delta_\varepsilon(\mathbf U^{pred})$, correct $(\alpha_k\rho_ke_k)^{n+1}$ → $\mathbf U^{n+1}$; then $\mathbf U^n\leftarrow\mathbf U^{n+1}$, $\Delta t\leftarrow\Delta t^*$.

## 5. Numerical results

### Common setup
- SG EOS parameters (**Table 1**):

| Fluid | $\gamma$ | $p_\infty$ (Pa) | $e_{ref}$ (J/kg) |
|---|---|---|---|
| Air | 1.4 | 0 | 0 |
| Water | 4.4 | $6\cdot10^8$ | 0 |
| Epoxy | 2.43 | $5.3\cdot10^9$ | 0 |
| Spinel | 1.62 | $141\cdot10^9$ | 0 |

- Domain: pipe $L=1$ m, constant section; $N_x=L/\Delta x$; meshes **200 and 1000** cells (one case also 5000).
- **IMEX time step:** $\Delta t_{IMEX} = 0.5\times\Delta x/|\sigma|$, $\sigma$ from (48) (i.e. CFL 0.5 on $D\max|u|$).
  Dissipation coefficient $D$ fixed per test case (independent of mesh), given in tables.
- Riemann problem initial data, interface at $x_{inter}$.
- Reference explicit scheme: Godunov-type FV (Sect. 3) with HLLC + second-order MUSCL–Hancock; $\Delta t_{HLLC} = 0.9\,\Delta x/\max|u\pm c_f|$.
- Performance metric: $R=\Delta t_{IMEX}/\Delta t_{HLLC}$ (at last iteration). Full perf comparison in App. E.
- *Note for implementation:* for a quiescent initial state ($u=0$ everywhere, e.g. shock tubes) $\sigma=0$ at $t=0$; the paper
  does not say how the first $\Delta t$ is chosen **[uncertain — not specified in the main text; check App. E / D]**.

### 5.1 Single-phase tests (**Table 2**)

| RP | $t_f$ (s) | $\rho_L$ | $u_L$ | $p_L$ (Pa) | $x_{inter}$ | $\rho_R$ | $u_R$ | $p_R$ (Pa) | $D$ |
|---|---|---|---|---|---|---|---|---|---|
| Low-Mach air advection | 0.5 | 1000 | 1 | $10^5$ | 0.25 | 0.01 | 1 | $10^5$ | 1.1 |
| Low-Mach water shock tube | $10^{-4}$ | 1000 | 0 | $10^8$ | 0.5 | 1000 | 15 | $0.98\cdot10^8$ | 40 |
| Lax shock tube (air) | 0.14 | 0.445 | 0.698 | 3.528 | 0.5 | 0.5 | 0 | 0.571 | 1.5 |
| Strong air shock tube | $2.5\cdot10^{-6}$ | 1 | 0 | $10^{10}$ | 0.5 | 0.125 | 0 | 0.1 | 1.2 |
| Two rarefaction fans (air) | 0.15 | 1 | $-1$ | 0.4 | 0.5 | 1 | 1 | 0.4 | 1 |

(Units kg/m³, m/s, Pa, m. Small-print table; values read from a low-resolution render — the water shock tube $u_R=15$ and $p_R=0.98\cdot10^8$
are consistent with Fig. 3 (velocity 0→8→15 m/s, pressure ≈ $10^8$/$8.6\cdot10^7$/$9.8\cdot10^7$). **[uncertain on exact digits of low-Mach water $p_R$]**.)

**Results:**
- *Low-Mach air advection* (Fig. 2, $t_f=0.5$ s): density step of five orders of magnitude advected at 1 m/s, uniform $p$. Mach $2.7\cdot10^{-4}$ to $8.5\cdot10^{-2}$.
  $u$, $p$ stay exactly uniform (no spurious oscillations); differences only in interface diffusion. IMEX (1000) is *less diffused* than HLLC (1000) because it takes far fewer steps.
  For 1000 cells $R=1890$; wall time $t_{HLLC}\approx 513$ s vs $t_{IMEX}=0.45$ s.
- *Low-Mach water shock tube* (Fig. 3, from [85], $t_f=10^{-4}$ s): water, tiny $p$ and $u$ jumps; Mach max $8.6\cdot10^{-3}$, interface at 8.07 m/s.
  Discontinuities correctly located, plateaus correct; IMEX shows **diffused acoustic discontinuities** (because the explicit dissipation is scaled with material speed and $D=40$ is large), converges with refinement. Demonstrates handling of $p_\infty\ne0$.
- *Lax shock tube* (Fig. 4, $t_f=0.14$): right shock, left rarefaction, contact at 1.53 m/s; Mach 0→0.95. Correct positions/plateaus; converges.
  **Effect of $D$** (Fig. 5, NX=200): $D=1$ → spurious oscillations at the shock/contact (density overshoot ≈1.35 near $x\approx0.84$); $D=1.5$ mitigates; $D=2$ removes them, with slight extra loss of precision at waves (e.g. contact at 0.73 m) and smaller $\Delta t$ (since $\sigma\propto D$). Compromise accepted.
- *Strong shock tube* (Fig. 6, from [86], $t_f=2.5\cdot10^{-6}$ s): pressure ratio 11 orders of magnitude; Mach 0→1.9; material wave $118\cdot10^3$ m/s. Slight density undershoot across the contact, improved by refinement; supersonic regime handled.
- *Two rarefaction fans* (Fig. 7, $t_f=0.15$): $u=\mp1$; Mach 0→1.34. Almost indistinguishable from exact; **non-physical oscillations in the central (near-vacuum, $u\approx0$) zone** for 200 cells (visible in $\rho$, $u$, $p$, Mach), mitigated by refinement.
  *(Interpretation: at the centre $u\approx0$ so the Rusanov dissipation $\propto D|u|$ vanishes — the material-speed-scaled numerical viscosity offers no damping where $u\to0$.)*

#### 5.1.1 Comparison with another low-Mach scheme ([8], fully implicit + preconditioning, pressure-disequilibrium model)
Liquid–liquid shock tube: $L=1$ m, $x_{inter}=0.5$, pure water both sides, $\rho_L=\rho_R=1000$ kg/m³, $p_L=10^6$ Pa, $p_R=10^5$ Pa, $t=0.15$ ms, NX=200.
Fig. 8: both get correct jumps/mean positions; implicit+preconditioning smears shocks/expansions much more; IMEX smears the contact more (due to $D$ in the explicit advection). IMEX shows small overshoots at the shock fronts in $u$.
(D value for this case not stated in the text **[uncertain]**.)

### 5.2 Two-phase tests
All with **instantaneous pressure relaxation** as in [19] ⇒ recovers the 5-eq. Kapila model, whose exact solutions are references.
Mach number uses the **Wood** (mechanical-equilibrium) sound speed:
$$\frac{1}{\rho c^2} = \sum_k\frac{\alpha_k}{\rho_kc_k^2}.\tag{50}$$

**Table 3** ($\varepsilon_\alpha = 10^{-7}$; phase 1 listed first):

| RP | $t_f$ (s) | $\alpha_{1,L}$ | $\rho_{1,L}$ | $\rho_{2,L}$ | $u_L$ | $p_L$ | $x_{inter}$ | $\alpha_{1,R}$ | $\rho_{1,R}$ | $\rho_{2,R}$ | $u_R$ | $p_R$ | $D$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Water(1)–air(2) advection | 0.25 | $1-\varepsilon_\alpha$ | 1000 | 0.01 | 1 | $10^5$ | 0.25 | $\varepsilon_\alpha$ | 1000 | 0.01 | 1 | $10^5$ | 3.2 |
| Water(1)–air(2) shock tube | $241\cdot10^{-6}$ | $1-\varepsilon_\alpha$ | 1000 | 50 | 0 | $10^9$ | 0.7 | $\varepsilon_\alpha$ | 1000 | 50 | 0 | $10^5$ | 2 |
| Epoxy(1)–spinel(2) shock tube | $30\cdot10^{-6}$ | 0.5954 | 1185 | 3622 | 0 | $10^{10}$ | 0.6 | 0.5954 | 1185 | 3622 | 0 | $10^5$ | 6 |

(Pressures in Pa, both phases at the same initial pressure.)

**Results:**
- *Low-Mach water–air advection* (Fig. 9, $t_f=0.25$ s): uniform $p=1$ bar, $u=1$ m/s; Mach ≈ $2.7\cdot10^{-3}$. $u$, $p$, $\rho_{water}$ stay uniform;
  volume-fraction discontinuity advected correctly with numerical diffusion. Mach spike up to 0.042 at the interface — not an oscillation, a consequence of the non-monotone Wood speed (50) in mixture cells (motivates the pressure-disequilibrium model, see [20]).
  $R=517.72$; wall time IMEX 0.66 s vs explicit ≈257.71 s (1000 cells). IMEX(1000) slightly sharper than HLLC(1000).
- *Water–air shock tube* (Fig. 10, from [20], adapted from [19], $t_f=241\,\mu$s; 200, 1000, **5000** cells): Mach 0 → 0.26 → max 1.84; interface moves at 482 m/s.
  Shock and interface nearly merge on density due to diffusion. **Pressure undershoot at the tail of the rarefaction** (log scale exaggerates), also seen in the literature [20,75]; its spatial extent shrinks with refinement (5000 cells shows a narrow dip).
- *Epoxy–spinel shock tube* (Fig. 11, from [19], $t_f=30\,\mu$s): same mixture both sides, different pressures; interface at 550 m/s, Mach max 0.136. Discontinuities located correctly, plateaus correct, convergence with refinement (see [19,70] for convergence of non-conservative multiphase shocks).

**Convergence study** (Fig. 12): relative $L_1$ error of mixture density at $t_f$, normalised as in [66]:
$$L_{1,\rho}(\text{NX}) = \frac{\sum_{i=1}^{\text{NX}}|\rho^{num}_i(t_f)-\rho^{ana}(x_i,t_f)|\Delta x}{\sum_{i=1}^{\text{NX}}|\rho^{ana}(x_i,t_0)|\Delta x}.\tag{51}$$
Regression $L_{1,\rho}(\text{NX}) = a\,\text{NX}^b$ with $b\approx-0.7$ (typical sub-first-order rate for discontinuous solutions).

Fig. 11 (epoxy–spinel, 200/1000 cells): mixture density, velocity (plateau ≈ 550 m/s), pressure, Mach (≈0.127/0.136 plateaus), epoxy volume fraction (0.5954 → ≈0.632 / ≈0.538), epoxy density (1185 → ≈1000 / ≈1525). Good agreement.

Fig. 12 fitted regressions (log–log, NX = 200, 1000):
- Water–air advection: $L_{1,\rho} = 1.22\,\text{NX}^{-0.65}$
- Water–air shock tube: $L_{1,\rho} = 0.84\,\text{NX}^{-0.72}$
- Epoxy–spinel shock tube: $L_{1,\rho} = 0.19\,\text{NX}^{-0.72}$

## 6. Conclusion (paraphrase)
- Semi-implicit scheme built from hyperbolic flux-splitting removes the acoustic restriction for single- and two-phase single-velocity flows; extension to $\Pi$ phases is straightforward (pressure equations already derived for $\Pi$ phases); general EOS.
- Semi-discrete spatial discretisation ⇒ compact stencil ⇒ banded matrix with small width. Explicit part second order (Rusanov + MUSCL–Hancock).
- Multi-D extension "presents no theoretical difficulties" apart from managing wider-band matrices.
- Big gains for low-Mach cases ($R\gg1$); for shock tubes $R\gtrsim1$ or $R\lesssim1$ (little or no gain).
- Future work: verify asymptotic-preserving property, semi-implicit boundary conditions, multi-D, steady-state cases.

## Appendix A. Low-Mach asymptotic analysis (of the hyperbolic part of (20))
Scaling (reference values $(\cdot)_r$, $t_r = x_r/u_r$, $c_{f,r}$ frozen speed (8)):
$$\tilde\alpha_k=\alpha_k,\ \tilde\rho_k=\rho_k/\rho_r,\ \tilde p_k = p_k/(\rho_rc_{f,r}^2),\ \tilde e_k = e_k/c_{f,r}^2,\ \tilde{\mathbf u}=\mathbf u/u_r,\ \tilde{\mathbf x}=\mathbf x/x_r,\ \tilde t=t/t_r\tag{A.1}$$
$$\tilde p=\sum_k\tilde\alpha_k\tilde p_k = p/(\rho_rc_{f,r}^2),\quad \tilde\rho=\sum_k\tilde\alpha_k\tilde\rho_k=\rho/\rho_r,\quad \tilde e=\sum_k(\tilde\alpha_k\tilde\rho_k\tilde e_k)/\tilde\rho = e/c_{f,r}^2.\tag{A.2}$$
With $M_r=u_r/c_{f,r}$:
$$
\begin{cases}
\partial_{\tilde t}\tilde\alpha_k + \tilde{\mathbf u}\cdot\nabla_{\tilde x}\tilde\alpha_k = 0\\
\partial_{\tilde t}(\tilde\alpha_k\tilde\rho_k) + \nabla_{\tilde x}\cdot(\tilde\alpha_k\tilde\rho_k\tilde{\mathbf u}) = 0\\
\partial_{\tilde t}(\tilde\alpha_k\tilde\rho_k\tilde e_k) + \nabla_{\tilde x}\cdot(\tilde\alpha_k\tilde\rho_k\tilde e_k\tilde{\mathbf u}) + \tilde\alpha_k\tilde p_k\nabla_{\tilde x}\cdot\tilde{\mathbf u} = 0\\
\partial_{\tilde t}(\tilde\rho\tilde{\mathbf u}) + \nabla_{\tilde x}\cdot(\tilde\rho\tilde{\mathbf u}\otimes\tilde{\mathbf u}) + \frac{1}{M_r^2}\nabla_{\tilde x}\tilde p = \mathbf 0
\end{cases}\tag{A.3}
$$
Expand every variable as $\phi = \phi_0 + M_r\phi_1 + M_r^2\phi_2 + \mathcal O(M_r^3)$ (A.4); leading orders give (A.5).
Cancelling $M_r^{-2}$, $M_r^{-1}$ terms: $\nabla\tilde p_0 = \mathbf 0$, $\nabla\tilde p_1=\mathbf 0$ (A.6); with well-prepared data/BCs [4,6,47], $\tilde p_1\equiv0$:
$\tilde p = \tilde p_0(\tilde t) + M_r^2\tilde p_2(\tilde x,\tilde t) + \mathcal O(M_r^3)$ (A.7).
Total energy: $\tilde E = \tilde e + \tfrac12M_r^2\tilde{\mathbf u}^2$, $\tilde E_0=\tilde e_0=\sum_k\tilde\alpha_{k,0}\tilde\rho_{k,0}\tilde e_{k,0}/\tilde\rho_0$ (A.8);
$\partial_{\tilde t}(\tilde\rho\tilde E) + \nabla\cdot((\tilde\rho\tilde E+\tilde p)\tilde{\mathbf u})=0$, $\tilde E = E/c_{f,r}^2$ (A.9);
at leading order $\partial_{\tilde t}(\tilde\rho_0\tilde e_0) + \nabla\cdot(\tilde\rho_0\tilde e_0\tilde{\mathbf u}_0) + \tilde p_0\nabla\cdot\tilde{\mathbf u}_0 + \mathcal O(M_r)=0$ (A.10).
EOS expansion $\tilde e_k = \tilde e_{k,0} + M_r\tilde e_{k,1} + \mathcal O(M_r^2)$, $\tilde e_{k,1} = \tilde e_{p,k,0}\tilde p_{k,1} + \tilde e_{\rho,k,0}\tilde\rho_{k,1}$ (A.11). Phase-$k$ pressure equation:
$$\tilde\alpha_{k,0}\partial_{\tilde t}\tilde p_{k,0} + \tilde\alpha_{k,0}\tilde{\mathbf u}_0\cdot\nabla\tilde p_{k,0} + \tilde\alpha_{k,0}\tilde\rho_{k,0}\tilde c_{k,0}^2\nabla\cdot\tilde{\mathbf u}_0 = 0,\tag{A.12}$$
$\tilde c_{k,0}^2 = (\tilde p_{k,0}/\tilde\rho_{k,0}^2 - \tilde e_{\rho,k,0})/\tilde e_{p,k,0}$. Summing with $\nabla\tilde p_0=0$:
$$\partial_{\tilde t}\tilde p_0 + \tilde\rho_0\tilde c_{f,0}^2\nabla\cdot\tilde{\mathbf u}_0 = 0,\quad \tilde c_{f,0}^2=\sum_k\tilde\alpha_{k,0}\tilde\rho_{k,0}\tilde c_{k,0}^2/\tilde\rho_0.\tag{A.13}$$
Open boundaries ⇒ $\tilde p_0$ constant ⇒ $\nabla\cdot\tilde{\mathbf u}_0 = 0$ (A.14), and incompressible limit
$$\nabla\cdot\tilde{\mathbf u}_0 = 0,\qquad \tilde\rho_0\partial_{\tilde t}\tilde{\mathbf u}_0 + \tilde\rho_0(\tilde{\mathbf u}_0\cdot\nabla)\tilde{\mathbf u}_0 + \nabla\tilde p_2 = \mathbf 0,\tag{A.15}$$
with pure transport of $\tilde\alpha_{k,0}$, $\tilde\rho_{k,0}$, $\tilde p_{k,0}$ (A.16).
(The paper does **not** prove the discrete AP property — explicitly left to future work.)

## Appendix B. Explicit flux-splitting scheme (TV-type, for reference; not the IMEX scheme)
$$\mathbf U^{pred}_i = \mathbf U^n_i - \frac{\Delta t}{\Delta x}\left(\left(\mathbf P^*_{i+\frac12}+\mathbf A^*_{i+\frac12}\right)-\left(\mathbf P^*_{i-\frac12}+\mathbf A^*_{i-\frac12}\right)\right) - \frac{\Delta t}{\Delta x}\mathbf H^n_i\left(u^*_{i+\frac12}-u^*_{i-\frac12}\right),\tag{B.1}$$
followed by the Sect. 3 energy correction.
**Pressure flux** $\mathbf P^*_{i+1/2} = (0,0,0,p^*,p^*u^*)^T$ (B.2). Riemann problem of the pressure subsystem (Fig. B.1): two acoustic waves $\sigma_L=-A$, $\sigma_R=+A$ and a stationary contact on the $t$-axis; always subsonic. $\alpha_k$, $\rho_k$, $e_k$ invariant across $\pm A$; Rankine–Hugoniot $\mathbf P^* = \mathbf P + \sigma(\mathbf U^*-\mathbf U)$ gives
$$u^* = \frac{(p_L - m_Lu_L) - (p_R - m_Ru_R)}{m_R - m_L},\tag{B.3a}$$
$$p^* = \frac{m_R(p_L-m_Lu_L) - m_L(p_R-m_Ru_R)}{m_R-m_L},\tag{B.3b}$$
$m_L=\rho_L\sigma_L$, $m_R=\rho_R\sigma_R$, Davis estimates [89]: $\sigma_L = \min(-A_L,-A_R)$, $\sigma_R=\max(A_L,A_R)$, $A$ from (13)/(14).
**Advection flux** (TV-style upwinding with $u^*$ from (B.3a)):
$$\mathbf A^*_{i+\frac12} = u^*_{i+\frac12}\begin{cases}\mathbf U^n_i & \text{if } u^*_{i+\frac12}>0\\ \mathbf U^n_{i+1}&\text{else}\end{cases}\tag{B.4}$$
(i.e. $\mathbf A^*$ = $u^*$ times the upwind state of the advected quantities $(\alpha_k,\alpha_k\rho_k,\alpha_k\rho_ke_k,\rho u,\rho E)$.)
Algorithm: (1) $p^*,u^*$ from (B.3) → $\mathbf P^*$; (2) $\mathbf A^*$ from (B.4); (3) $\mathbf U^{pred}$ from (B.1) with $u^*$ from (B.3a); (4) energy correction.
Stability: $|A|$-based CFL (still acoustic). **This is the closest piece to a fully explicit TV splitting of the 6-eq. model** and is directly reusable as the explicit baseline.

## Appendix C. Coefficients of the implicit subsystem
Insert velocity into pressure equation:
$$p^{n+\star,r+1}_{k,i} = p^n_{k,i} - \Delta t(p_{k,i}\Gamma_{k,i})^{n+\star,r}\left[\nabla\cdot\mathbf u^n_i - \Delta t\,\nabla\cdot\sum_k\left(\frac{\nabla\alpha^n_{k,i}}{\rho^n_i}p^{n+\star,r+1}_{k,i} + \frac{\alpha^n_{k,i}}{\rho^n_i}\nabla p^{n+\star,r+1}_{k,i}\right)\right].$$
$$\nu^{n+\star,r}_{k,i} = \frac{1}{\rho^n_{k,i}}\left(\frac{p_{k,i}}{e_{p,k,i}}\right)^{n+\star,r} = (p_{k,i}\Gamma_{k,i})^{n+\star,r}.\tag{C.1}$$
(For SG: $e_{p,k} = 1/((\gamma_k-1)\rho_k)$ ⇒ $\nu_k = (\gamma_k-1)p_k$. Note it is **linear in $p_k$** for SG, so the Picard lag only concerns $\nu_k$ ∝ $p_k^{r}$.)
$$p^{r+1}_{k,i} - (\Delta t)^2\nu^{r}_{k,i}\sum_l\nabla\cdot\left(\frac{\nabla\alpha^n_{l,i}}{\rho^n_i}p^{r+1}_{l,i} + \frac{\alpha^n_{l,i}}{\rho^n_i}\nabla p^{r+1}_{l,i}\right) = p^n_{k,i} - \Delta t\,\nu^r_{k,i}\nabla\cdot\mathbf u^n_i,$$
expanded:
$$p^{r+1}_{k,i} - (\Delta t)^2\nu^r_{k,i}\sum_l\left[\frac{\alpha^n_{l,i}}{\rho^n_i}\nabla\cdot\nabla p^{r+1}_{l,i} + \left(\nabla\Big(\frac{\alpha^n_{l,i}}{\rho^n_i}\Big)+\frac{\nabla\alpha^n_{l,i}}{\rho^n_i}\right)\cdot\nabla p^{r+1}_{l,i} + \nabla\cdot\Big(\frac{\nabla\alpha^n_{l,i}}{\rho^n_i}\Big)p^{r+1}_{l,i}\right] = p^n_{k,i} - \Delta t\,\nu^r_{k,i}\nabla\cdot\mathbf u^n_i.$$
1D operators on the unknown: $\nabla p_{k,i} = \frac{p_{k,i+1}-p_{k,i-1}}{2\Delta x}$, $\nabla\cdot\nabla p_{k,i} = \frac{p_{k,i+1}-2p_{k,i}+p_{k,i-1}}{(\Delta x)^2}$
(**compact 3-point Laplacian**, not the 5-point composition of two centred gradients — this is the point of substituting before discretising).
Resulting (34) coefficients (written with phase index; in (34) they carry the summed index $l$):
$$\eta^n_{k,i} = \frac{1}{\rho^n_i}\left[\frac{\alpha^n_{k,i}}{(\Delta x)^2} - \frac{1}{2\Delta x}\left(2\nabla\alpha^n_{k,i} - \frac{\alpha^n_{k,i}}{\rho^n_i}\nabla\rho^n_i\right)\right],\tag{C.2a}$$
$$\zeta^n_{k,i} = \frac{1}{\rho^n_i}\left[\nabla\cdot\nabla\alpha^n_{k,i} - \frac{\nabla\alpha^n_{k,i}\cdot\nabla\rho^n_i}{\rho^n_i} - 2\frac{\alpha^n_{k,i}}{(\Delta x)^2}\right],\tag{C.2b}$$
$$\beta^n_{k,i} = \frac{1}{\rho^n_i}\left[\frac{\alpha^n_{k,i}}{(\Delta x)^2} + \frac{1}{2\Delta x}\left(2\nabla\alpha^n_{k,i} - \frac{\alpha^n_{k,i}}{\rho^n_i}\nabla\rho^n_i\right)\right].\tag{C.2c}$$
(I re-derived these: $\nabla(\alpha/\rho)+\nabla\alpha/\rho = (2\nabla\alpha - \alpha\nabla\rho/\rho)/\rho$ and $\nabla\cdot(\nabla\alpha/\rho) = \Delta\alpha/\rho - \nabla\alpha\cdot\nabla\rho/\rho^2$ — consistent.)

Two-phase coefficients of (35) (all $\nu$ at $(n+\star,r)$, all $\eta,\zeta,\beta$ at $n$):
$$
\begin{aligned}
&a_{1,i} = -(\Delta t)^2\nu_{1,i}\eta_{1,i}, && b_{1,i} = -(\Delta t)^2\nu_{1,i}\eta_{2,i},\\
&a_{2,i} = 1-(\Delta t)^2\nu_{1,i}\zeta_{1,i}, && b_{2,i} = -(\Delta t)^2\nu_{1,i}\zeta_{2,i},\\
&a_{3,i} = -(\Delta t)^2\nu_{1,i}\beta_{1,i}, && b_{3,i} = -(\Delta t)^2\nu_{1,i}\beta_{2,i},\\
&c_i = p^n_{1,i} - \Delta t\,\nu_{1,i}\nabla\cdot\mathbf u^n_i,\\[3pt]
&m_{1,i} = -(\Delta t)^2\nu_{2,i}\eta_{1,i}, && q_{1,i} = -(\Delta t)^2\nu_{2,i}\eta_{2,i},\\
&m_{2,i} = -(\Delta t)^2\nu_{2,i}\zeta_{1,i}, && q_{2,i} = 1-(\Delta t)^2\nu_{2,i}\zeta_{2,i},\\
&m_{3,i} = -(\Delta t)^2\nu_{2,i}\beta_{1,i}, && q_{3,i} = -(\Delta t)^2\nu_{2,i}\beta_{2,i},\\
&g_i = p^n_{2,i} - \Delta t\,\nu_{2,i}\nabla\cdot\mathbf u^n_i.
\end{aligned}
$$
1D discretisation of the explicit/known operators (C.3):
$$\nabla\alpha_{k,i} = \frac{\alpha_{k,i+1}-\alpha_{k,i-1}}{2\Delta x},\quad \nabla\cdot\nabla\alpha_{k,i} = \frac{\alpha_{k,i+1}-2\alpha_{k,i}+\alpha_{k,i-1}}{(\Delta x)^2},\quad \nabla\rho_i = \frac{\rho_{i+1}-\rho_{i-1}}{2\Delta x},\quad \nabla\cdot\mathbf u_i = \frac{u_{i+1}-u_{i-1}}{2\Delta x}.\tag{C.3}$$
Note $\nabla\cdot u^n$ in the RHS is a **centred, collocated** divergence (wide stencil on $u$), while the Laplacian on $p$ is compact.
No boundary treatment of the matrix rows is given in the paper **[not specified]** (future work mentions "semi-implicit boundary conditions").

## Appendix D. MUSCL–Hancock for the advection subsystem (Toro [78])
1. **Reconstruction** in primitive variables $\mathbf W^{n+\star} = (\alpha_1,\alpha_2,\rho_1,\rho_2,p_1,p_2,u)^{n+\star}$:
   $\mathbf W^{n+\star}_{i,L} = \mathbf W^{n+\star}_i - \tfrac12\Delta_i$, $\mathbf W^{n+\star}_{i,R} = \mathbf W^{n+\star}_i + \tfrac12\Delta_i$ (D.1), slopes $\Delta_i$ with **minmod** limiter [90].
2. **Half-time prediction** (local, no Riemann problem):
$$\overline{\mathbf U}^{n+\star}_i = \mathbf U^{n+\star}_i - \frac{\Delta t/2}{\Delta x}\left(\mathbf A(\mathbf U(\mathbf W^{n+\star}_{i,R})) - \mathbf A(\mathbf U(\mathbf W^{n+\star}_{i,L}))\right) - \frac{\Delta t/2}{\Delta x}\mathbf H^{n+\star}_{a,i}\left(u^{n+\star}_{i,R}-u^{n+\star}_{i,L}\right).\tag{D.2}$$
   Convert $\overline{\mathbf U}_i\to\overline{\mathbf W}_i$ and re-extrapolate with the same slopes: $\overline{\mathbf W}_{i,L}$, $\overline{\mathbf W}_{i,R}$ (as in (D.1)).
   **[uncertain]**: whether "as in (D.1)" re-uses $\Delta_i$ or recomputes slopes — standard MUSCL–Hancock reuses $\Delta_i$.
3. **Riemann problem:** $\mathbf U^{n+\star,+}_{i-1/2} \equiv \overline{\mathbf U}(\overline{\mathbf W}_{i,L})$, $\mathbf U^{n+\star,-}_{i+1/2}\equiv\overline{\mathbf U}(\overline{\mathbf W}_{i,R})$; Rusanov flux (47) with $\sigma$ from (48) on these states.
   (Note: the half-step in (D.2) uses only the advection fluxes; the pressure contribution was already applied in $[\mathcal I]$.)

## Appendix E. Performance (Tables E.1, E.2)

**Table E.1** — time steps (s):

| Case | $\Delta t^{200}_{IMEX}$ | $\Delta t^{200}_{HLLC}$ | $\Delta t^{1000}_{IMEX}$ | $\Delta t^{1000}_{HLLC}$ | $R$ |
|---|---|---|---|---|---|
| Low-Mach air advection | $2.27\cdot10^{-3}$ | $1.20\cdot10^{-6}$ | $4.55\cdot10^{-4}$ | $2.40\cdot10^{-7}$ | 1890 |
| Low-Mach water shock tube | $4.17\cdot10^{-6}$ | $2.55\cdot10^{-6}$ | $8.33\cdot10^{-7}$ | $5.09\cdot10^{-7}$ | 1.64 |
| Lax shock tube | $1.09\cdot10^{-3}$ | $9.56\cdot10^{-4}$ | $2.17\cdot10^{-4}$ | $1.91\cdot10^{-4}$ | 1.14 |
| Strong air shock tube | $1.75\cdot10^{-8}$ | $2.10\cdot10^{-8}$ | $3.52\cdot10^{-9}$ | $4.21\cdot10^{-9}$ | 0.83 |
| Two rarefaction fans | $2.50\cdot10^{-3}$ | $2.57\cdot10^{-3}$ | $5.00\cdot10^{-4}$ | $5.15\cdot10^{-4}$ | 0.97 |
| Water–air advection | $6.25\cdot10^{-4}$ | $1.21\cdot10^{-6}$ | $1.25\cdot10^{-4}$ | $2.41\cdot10^{-7}$ | 517.7 |
| Water–air shock tube | $2.57\cdot10^{-6}$ | $1.70\cdot10^{-6}$ | $5.13\cdot10^{-7}$ | $3.39\cdot10^{-7}$ | 1.51 |
| Epoxy–spinel shock tube | $7.55\cdot10^{-7}$ | $5.75\cdot10^{-7}$ | $1.51\cdot10^{-7}$ | **[uncertain: rendered like $3.39\cdot10^{-7}$; $R=1.32$ implies ≈ $1.14\cdot10^{-7}$]** | 1.32 |

Sanity check with $\Delta t_{IMEX} = 0.5\Delta x/(D\max|u|)$, $\Delta x = 5\cdot10^{-3}$ (200 cells): air advection $0.0025/(1.1\cdot1)=2.27\cdot10^{-3}$ ✓; water shock tube $0.0025/(40\cdot15)=4.17\cdot10^{-6}$ ✓;
Lax $0.0025/(1.5\cdot1.53)=1.09\cdot10^{-3}$ ✓; water–air shock $0.0025/(2\cdot\approx486)$ ✓; epoxy–spinel $0.0025/(6\cdot552)$ ✓;
**water–air advection: $0.0025/6.25\cdot10^{-4} = 4 = D\cdot1$ ⇒ $D=4$, not 3.2 as I read in Table 3 [uncertain: either the table value or my reading is off]**.
This confirms $\sigma = D\max|u|$ and CFL = 0.5, and that $D$ directly shrinks the time step.

**Table E.2** — wall-clock times (s):

| Case | $t^{200}_{IMEX}$ | $t^{200}_{HLLC}$ | $t^{1000}_{IMEX}$ | $t^{1000}_{HLLC}$ |
|---|---|---|---|---|
| Low-Mach air advection | 0.03 | 20.22 | 0.45 | 513.50 |
| Low-Mach water shock tube | 0.01 | 0.01 | 0.11 | 0.06 |
| Lax shock tube | 0.04 | 0.02 | 0.50 | 0.20 |
| Strong air shock tube | 0.06 | 0.01 | 0.72 | 0.18 |
| Two rarefaction fans | 0.02 | 0.01 | 0.20 | 0.10 |
| Water–air advection | 0.04 | 10.21 | 0.66 | 257.71 |
| Water–air shock tube | 0.06 | 0.02 | 0.92 | 0.20 |
| Epoxy–spinel shock tube | 0.02 | 0.01 | 0.22 | 0.09 |

Text: semi-implicit generally $R\gtrsim1$, hugely better at very low Mach ($R\gg1$); for higher-Mach cases $R\lesssim1$ and IMEX wall time is 2–5× **larger** (Picard + Gauss–Seidel cost per step). All runs < 1 s except explicit low-Mach advection.

---

## 7. Implementation remarks, pitfalls, positivity & conservation (my synthesis)

**What the scheme actually is.**
- First-order Lie splitting $[C]\circ[\mathcal E]\circ[\mathcal I]$; implicit = backward Euler on the pressure subsystem with Picard-lagged $\nu_k=p_k\Gamma_k$; explicit = forward-Euler/MUSCL–Hancock on advection. **No IMEX-RK, no Butcher tableau, no Strang splitting**, no stiff relaxation inside the IMEX step (pressure relaxation is a separate instantaneous step after $[C]$, as in [19]; no thermal/chemical relaxation at all).
- Unlike the classical TV splitting (Appendix B) where pressure fluxes come from a Riemann solver, the implicit stage uses **centred collocated** differences (no upwinding/Riemann solver for $p$), dissipation comes only from the explicit Rusanov with $\sigma=D\max|u|$ and from the backward-Euler acoustic step.

**Conservation.**
- Conserved: $\alpha_k\rho_k$ (untouched in $[\mathcal I]$, conservative in $[\mathcal E]$), $\rho u$ (centred flux $\tfrac12(p_i+p_{i+1})$ with $p=\sum\alpha_kp_k$, since (38) = difference of face averages), $\rho E$ (flux $(pu)_{i+1/2}$ in $[\mathcal I]$ and Rusanov in $[\mathcal E]$). Phase internal energies are non-conservative and fixed by $[C]$ (weights $\alpha^n_k$).
- Subtle: the momentum update (38) uses **$\alpha^n_{k,i\pm1}p_{k,i\pm1}$** in a centred difference, while the matrix (34)/App. C uses the expanded product rule ($\nabla\alpha\,p + \alpha\nabla p$, compact Laplacian). The two are consistent to truncation error but not identical discretely; the velocity implied by the linear system ≠ the velocity from (38)/(40). Energy consistency is restored only via $[C]$.
- Energy (41) uses face velocity (40) $u_{i+1/2} = (\rho u_i+\rho u_{i+1})/(\rho_i+\rho_{i+1})$ at $n+\star$ — need $(\rho u)^{n+\star}$ from (38) *before* (39)–(41).

**Positivity / robustness.**
- No positivity proof for $\alpha_k$, $\rho_k$, $p_k$ or $e_k$. $\alpha_k$, $\alpha_k\rho_k$ positivity in $[\mathcal E]$ follows from the Rusanov CFL only if $\sigma\ge\max|u|$ ($D\ge1$) and CFL ≤ 0.5 — and the non-conservative $-\alpha_k\nabla\cdot u$ term can still break it **[my remark]**.
- Phase pressures from the implicit step can go negative for large $\Delta t$ (expansion), especially with SG liquids where $p_k$ may legitimately be < 0 but $p_k\Gamma_k>0$ is needed for hyperbolicity of (10) and for $\nu_k>0$ (else the operator $1-\Delta t^2\nu\,\mathcal L$ loses its M-matrix character).
- Gauss–Seidel convergence: the $2\times2$ cell blocks couple phases strongly ($m_{2,i}=-(\Delta t)^2\nu_2\zeta_1$ can be large when $\alpha_1\approx1$ and $\nu_2$ large), so the matrix is **not diagonally dominant** in general for large $\Delta t/\Delta x$; also $\eta,\beta$ can change sign where $|\nabla\alpha|$ or $|\nabla\rho|/\rho$ is large (sharp interfaces: when $\Delta x\,|2\nabla\alpha - \alpha\nabla\rho/\rho| > 2\alpha$). Expect convergence trouble at sharp density jumps; a direct banded/block-tridiagonal solver (block Thomas on $2\times2$ blocks — the system is exactly block-tridiagonal in 1D) is more robust than GS **[recommendation]**.
- Large $D$ (up to 40 for low-Mach water, 6 for epoxy–spinel) needed to suppress oscillations — $D$ is a tuning parameter, multiplies dissipation and divides $\Delta t$. Oscillations where $u\approx0$ (two-rarefaction case) because dissipation $\propto|u|$ vanishes.
- Acoustic waves are smeared (backward Euler + large acoustic CFL), visible in the water shock tube; shock tubes give little speed-up.
- Checkerboard risk (collocated $p,u$, centred gradient) acknowledged for multi-D.

**Things not specified in the paper** (decide yourself): boundary conditions for the implicit system; first $\Delta t$ when $u\equiv0$; exact Picard norm (absolute vs relative; $\varepsilon=10^{-6}$); max Picard iterations; Gauss–Seidel tolerance; which $\alpha$ weights $\alpha^n$ vs $\alpha^{pred}$ in $[C]$ (paper writes $\alpha^n_{k,i}$); pressure relaxation algorithm details (refer to [19,20]); the reading $A^2$ formula (14) lacking $1/\rho$.

**Relevance to a TV-split IMEX for the 6-equation model with relaxations:**
- Their pressure/advection split (10)/(16) and pressure-subsystem eigenstructure ($0^5,\pm A$, $A^2=\rho^{-1}\sum\alpha_kp_k\Gamma_k$) and the advection-subsystem matrix (17) are directly reusable; App. B gives explicit TV fluxes ($p^*,u^*$ via Davis-bounded acoustic solver, upwind $u^*$ advection).
- Their implicit pressure equation per phase uses $\nu_k=p_k\Gamma_k$ (not $\rho_kc_k^2$), because $p_k$ evolution in the pressure subsystem lacks the $\rho$-variation part — note this is the TV/ZB-specific "pressure" wave speed $A<c_f$.
- They keep relaxation outside the IMEX step (Lie splitting); a stiff-relaxation-coupled IMEX (e.g. ARS/SSP IMEX-RK with relaxation implicit) would be a genuine extension beyond this paper.
- Energy correction $[C]$ with weights $\alpha_k$ is the mechanism to retain total-energy conservation with non-conservative phase energies.

