# Reading note: Pelanti (2022), "Arbitrary-rate relaxation techniques for the numerical modeling of compressible two-phase flows with heat and mass transfer"

- **Author:** Marica Pelanti (IMSIA, CEA, CNRS, EDF, ENSTA Paris, Institut Polytechnique de Paris, 828 Boulevard des Maréchaux, 91120 Palaiseau, France). Email: marica.pelanti@ensta-paris.fr
- **Venue:** *International Journal of Multiphase Flow* **153** (2022) 104097
- **DOI:** https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
- Received 13 June 2021; revised 21 Feb 2022; accepted 9 Apr 2022; online 20 Apr 2022.
- MSC: 65M08, 76T10. Keywords: multiphase compressible flows, relaxation processes, liquid–vapor phase transition, finite volume schemes, Riemann solvers.

Notation conventions in this note: equation numbers in parentheses, e.g. (1a), are the paper's. Anything I could not read with certainty from the PDF is marked **[uncertain]**. My own comments (not in the paper) are marked *Note:*.

---

## 0. Abstract / one-paragraph summary

Single-velocity six-equation model (2 phasic masses, 2 phasic total energies, 1 volume fraction, mixture momentum) with relaxation source terms for volume (pressure), heat (temperature) and mass (Gibbs free energy / chemical potential) transfer. Fractional-step: homogeneous hyperbolic part by a second-order wave-propagation scheme with an HLLC-type Riemann solver; then a sequence of three ODE systems for the relaxation sources. Novelty: relaxation techniques (i) capable of **arbitrary relaxation time** (finite-rate or instantaneous heat and mass transfer) and (ii) applicable to a **general EOS**. The key idea is to write the relaxation ODEs using the equations of the *p-relaxed* (5-eq.) and *pT-relaxed* (4-eq.) models, freeze some coefficients over the step, and obtain semi-exact **exponential** solutions. The equilibrium pressure is always recomputed from the mixture energy relation (12) to keep mixture-energy-consistency.

---

## 1. Introduction (p. 1–2) — key points

- Hierarchy: Baer–Nunziato 7-eq. (two velocities) → 6-eq. single-velocity model (instantaneous kinetic equilibrium) → 5-eq. *p*-relaxed (Kapila et al. 2001) → 4-eq. *pT*-relaxed → 3-eq. *pTg*-relaxed = HEM. Also mentions HRM (Bilicki–Kestin; Downar-Zapolski) which is outside this hierarchy (mechanical equilibrium + vapour at saturation, thermal disequilibrium).
- The numerical method for the 6-eq. model must also approximate the relaxed models when instantaneous relaxations are activated.
- Literature relaxation procedures usually assume (a) instantaneous processes (algebraic equilibrium systems) and (b) stiffened gas (SG) EOS. Instantaneous thermo-chemical relaxation is inadequate for e.g. fast depressurisations, delayed vaporisation, metastable states (Downar-Zapolski 1996, Faucher 2000, Lund–Aursand 2012, De Lorenzo 2017).
- Previous own work (De Lorenzo et al. 2019, 2021): finite-rate heat/mass transfer with general EOS (IAPWS-97) solved with **explicit Runge–Kutta with adaptive step** — not suitable for stiff problems; implicit solvers too expensive. Appendix of De Lorenzo 2019 briefly proposed exponential solutions but *a priori* assumed exponential decay of $p,T,g$ differences; here the exponential systems are **derived** from the equations of the parent/relaxed two-phase models assuming some quantities constant during the relaxation step.
- Related exponential-solution approaches: Downar-Zapolski 1996, Helzel 2000, Pelanti–LeVeque 2006, Aursand 2014, De Lorenzo 2017, 2019.
- Procedures are simple, robust, usable for other models of the BN hierarchy, and guarantee consistency of relaxed states with the mixture pressure law (mixture-energy-consistent, Pelanti & Shyue 2014).
- HLLC solver can also be reinterpreted as a Suliciu-type solver (De Lorenzo et al. 2018).

---

## 2. Six-equation model (Section 2, p. 2–3)

### 2.1 Variables

Phases $k=1,2$, common velocity $\vec u$. $\alpha_k,\rho_k,p_k,\varepsilon_k$: volume fraction, density, pressure, specific internal energy. $\alpha_1+\alpha_2=1$.

- $\mathcal E_k=\rho_k\varepsilon_k$ (phasic internal energy per unit volume); $E_k=\mathcal E_k+\rho_k\frac{|\vec u|^2}{2}$ (phasic total energy per unit volume).
- $\rho=\sum_k\alpha_k\rho_k$; $\mathcal E=\sum_k\alpha_k\mathcal E_k=\rho\varepsilon$; $E=\sum_k\alpha_kE_k=\mathcal E+\frac12\rho|\vec u|^2$.
- $Y_k=\alpha_k\rho_k/\rho$ (mass fraction), $\varepsilon=Y_1\varepsilon_1+Y_2\varepsilon_2$.
- $p_m=\alpha_1p_1+\alpha_2p_2$ (mixture pressure); $p$ = mixture equilibrium pressure; $p_I$ interface pressure.
- $h_k=(\mathcal E_k+p_k)/\rho_k$, $h=Y_1h_1+Y_2h_2$.
- $c_k=\sqrt{(\partial p_k/\partial\rho_k)_{s_k}}$; frozen mixture sound speed $c_f=\sqrt{Y_1c_1^2+Y_2c_2^2}$.
- $T_k$ phasic temperature, $T$ equilibrium temperature; $s_k$, $s=Y_1s_1+Y_2s_2$; $g_k$ phasic chemical potential (specific Gibbs energy), $g_I$ interface chemical potential.
- $\Gamma_k=(\partial p_k/\partial\mathcal E_k)_{\rho_k}$ (Grüneisen); $\chi_k=(\partial p_k/\partial\rho_k)_{\mathcal E_k}$.
- $\phi_k=(\partial\rho_k/\partial T_k)_{p_k}=-\rho_k\beta_k$ ($\beta_k$ thermal expansion coefficient); $\zeta_k=(\partial\rho_k/\partial p_k)_{T_k}=\rho_k\mathcal K_{Tk}$ ($\mathcal K_{Tk}$ isothermal compressibility).
- $\mathcal K_{Sk}=1/(\rho_kc_k^2)$ isentropic compressibility.
- $\kappa_{pk}=T_k(\partial s_k/\partial T_k)_{p_k}=(\partial h_k/\partial T_k)_{p_k}$; $\kappa_{vk}=T_k(\partial s_k/\partial T_k)_{\rho_k}=(\partial\varepsilon_k/\partial T_k)_{\rho_k}$.
- $C_{pk}=\alpha_k\rho_k\kappa_{pk}$: phasic *extensive* heat capacity at constant pressure.

### 2.2 Equations (1a)–(1g)

$$\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mathcal P \tag{1a}$$
$$\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\mathcal M \tag{1b}$$
$$\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\mathcal M \tag{1c}$$
$$\partial_t(\rho\vec u)+\nabla\cdot\big(\rho\vec u\otimes\vec u+(\alpha_1p_1+\alpha_2p_2)\mathbb I\big)=0 \tag{1d}$$
$$\partial_t(\alpha_1E_1)+\nabla\cdot\big(\alpha_1(E_1+p_1)\vec u\big)+\Sigma=-p_I\mathcal P+\mathcal Q+\Big(g_I+\frac{|\vec u|^2}{2}\Big)\mathcal M \tag{1e}$$
$$\partial_t(\alpha_2E_2)+\nabla\cdot\big(\alpha_2(E_2+p_2)\vec u\big)-\Sigma=p_I\mathcal P-\mathcal Q-\Big(g_I+\frac{|\vec u|^2}{2}\Big)\mathcal M \tag{1f}$$
Non-conservative term:
$$\Sigma=-\vec u\cdot\vec\Xi,\qquad \vec\Xi=Y_2\nabla(\alpha_1p_1)-Y_1\nabla(\alpha_2p_2). \tag{1g}$$

$\mathcal P,\mathcal Q,\mathcal M$ = volume, heat, mass transfer:
$$\mathcal P=\mu(p_1-p_2),\qquad \mathcal Q=\vartheta(T_2-T_1),\qquad \mathcal M=\nu(g_2-g_1). \tag{2}$$
$\mu,\vartheta,\nu\ge0$ are parameters (or more generally functions) giving the rates of mechanical, thermal, chemical relaxation. **In this paper mechanical relaxation is always instantaneous ($\mu\to+\infty$)**; the parent model with instantaneous pressure relaxation approximates the pressure-equilibrium model (9). Thermal and chemical relaxation: no assumption, finite or infinite $\vartheta,\nu$.

Second-law (entropy production) sufficient conditions (Flåtten & Lund 2011): $\mu,\vartheta,\nu\ge0$ and
$$p_I\in\{\min(p_1,p_2),\max(p_1,p_2)\},\qquad g_I\in\{\min(g_1,g_2),\max(g_1,g_2)\} \tag{3}$$
*Note:* the paper's braces mean the interval $[\min,\max]$, i.e. convex combinations. Choice used in tests: $p_I=\dfrac{Z_2^ap_1+Z_1^ap_2}{Z_1^a+Z_2^a}$ with $Z_k^a=\rho_kc_k$ acoustic impedance (Saurel et al. 2003, 2009). Alternative: $p_I=\alpha_1p_2+\alpha_2p_1$ (Saurel–Abgrall 1999). $g_I$ is **not needed** in the numerical scheme (see later).

EOS: each phase needs $p_k(\mathcal E_k,\rho_k)$ and $T_k(p_k,\rho_k)$ (complete EOS when heat/mass transfer present; incomplete suffices otherwise).

### 2.3 Hyperbolicity and sound speeds

Eigenvalues in direction $\vec n$: $\lambda_{1,5+d}=\vec u\cdot\vec n\mp c_f$, $\lambda_l=\vec u\cdot\vec n$ for $l=2,\dots,4+d$ (multiplicity $3+d$).
$$c_f=\sqrt{Y_1c_1^2+Y_2c_2^2} \tag{4}$$
$$c_k=\sqrt{\Gamma_kh_k+\chi_k},\quad h_k=(\mathcal E_k+p_k)/\rho_k \tag{5}$$
$$\Gamma_k=\Big(\frac{\partial p_k}{\partial\mathcal E_k}\Big)_{\rho_k},\qquad \chi_k=\Big(\frac{\partial p_k}{\partial\rho_k}\Big)_{\mathcal E_k}. \tag{6}$$

Mixture total energy conservation (sum of (1e)+(1f)):
$$\partial_tE+\nabla\cdot\big(E\vec u+(\alpha_1p_1+\alpha_2p_2)\vec u\big)=0. \tag{7}$$

Phasic pressure equations:
$$\partial_tp_k+\vec u\cdot\nabla p_k+\rho_kc_k^2\nabla\cdot\vec u=\frac{(-1)^{k-1}}{\alpha_k}\Big\{-[\Gamma_k(\mathcal E_k+p_I)+\chi_k\rho_k]\mathcal P+\Gamma_k\mathcal Q+(\Gamma_kg_I+\chi_k)\mathcal M\Big\},\quad k=1,2. \tag{8}$$

---

## 3. Hierarchy of relaxed models (Section 3, p. 3–5)

### 3.1 Five-equation *p*-relaxed model ($\mu\to+\infty$, $p_1=p_2=p$)

Kapila et al. (2001) model ($4+d$ equations):
$$\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1-\frac{\alpha_1\alpha_2}{D}(\rho_2c_2^2-\rho_1c_1^2)\nabla\cdot\vec u=\frac1D(\alpha_2\Gamma_1+\alpha_1\Gamma_2)\mathcal Q+\frac1D(\alpha_2\omega_{I1}+\alpha_1\omega_{I2})\mathcal M \tag{9a}$$
$$\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\mathcal M \tag{9b}$$
$$\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\mathcal M \tag{9c}$$
$$\partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u+p\mathbb I)=0 \tag{9d}$$
$$\partial_tE+\nabla\cdot((E+p)\vec u)=0 \tag{9e}$$
with
$$D=\alpha_1\rho_2c_2^2+\alpha_2\rho_1c_1^2 \tag{10}$$
$$\omega_{Ik}=\Gamma_k(g_I-h_k)+c_k^2=\Gamma_kg_I+\chi_k,\quad k=1,2. \tag{11}$$
Derivation in Appendix A. Mixture pressure law $p=p(\mathcal E,\rho_1,\rho_2,\alpha_1)$ determined (given phasic energy laws $\mathcal E_k(p_k,\rho_k)$) by the mixture energy relation
$$\mathcal E=\alpha_1\mathcal E_1(p,\rho_1)+\alpha_2\mathcal E_2(p,\rho_2) \tag{12}$$
(isobaric condition). Wood speed of sound:
$$c_p=\Big(\rho\sum_{k=1}^2\frac{\alpha_k}{\rho_kc_k^2}\Big)^{-1/2}. \tag{13}$$
Note $\frac1D=\frac{\rho c_p^2}{\rho_1c_1^2\rho_2c_2^2}$.

Pressure equation:
$$\partial_tp+\vec u\cdot\nabla p+\rho c_p^2\nabla\cdot\vec u=\frac1D\Big[(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\mathcal Q+(\rho_2c_2^2\omega_{I1}-\rho_1c_1^2\omega_{I2})\mathcal M\Big]. \tag{14}$$

Phasic temperature equations ($k=1,2$):
$$\partial_tT_k+\vec u\cdot\nabla T_k+\frac{\rho c_p^2}{\phi_k}\Big(-\zeta_k+\frac{1}{c_k^2}\Big)\nabla\cdot\vec u
=\frac{1}{\phi_kD}\Big[(-1)^k\frac{\rho_k}{\alpha_k}(\alpha_2\Gamma_1+\alpha_1\Gamma_2)-\zeta_k(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big]\mathcal Q$$
$$\qquad+\frac{1}{\phi_k}\Big\{\frac{(-1)^{k-1}}{\alpha_k}+\frac1D\Big[(-1)^k\frac{\rho_k}{\alpha_k}(\alpha_2\omega_{I1}+\alpha_1\omega_{I2})-\zeta_k(\rho_2c_2^2\omega_{I1}-\rho_1c_1^2\omega_{I2})\Big]\Big\}\mathcal M \tag{15}$$
with
$$\phi_k=\Big(\frac{\partial\rho_k}{\partial T_k}\Big)_{p_k}=-\rho_k\beta_k,\qquad \zeta_k=\Big(\frac{\partial\rho_k}{\partial p_k}\Big)_{T_k}=\rho_k\mathcal K_{Tk}, \tag{16}$$
$$\beta_k=\frac{\Gamma_k\kappa_{pk}}{c_k^2}=\frac{\Gamma_kC_{pk}}{c_k^2\alpha_k\rho_k},\qquad \mathcal K_{Tk}=\mathcal K_{Sk}+\frac{\beta_k^2T_k}{\rho_k\kappa_{pk}}=\frac{1}{\rho_kc_k^2}+\frac{\beta_k^2T_k\alpha_k}{C_{pk}}, \tag{17}$$
$\mathcal K_{Sk}=1/(\rho_kc_k^2)$, $C_{pk}=\alpha_k\rho_k\kappa_{pk}$, $\kappa_{pk}=\partial h_k/\partial T_k|_{p_k}=T_k\,\partial s_k/\partial T_k|_{p_k}$.

### 3.2 Four-equation *pT*-relaxed model ($\mu,\vartheta\to+\infty$; $p_1=p_2=p$, $T_1=T_2=T$)

$$\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\mathcal M,\ \ \partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\mathcal M,\ \ \partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u+p\mathbb I)=0,\ \ \partial_tE+\nabla\cdot((E+p)\vec u)=0. \tag{18a–d}$$
Mixture pressure from (12) with isothermal condition $T_1(p,\rho_1)=T_2(p,\rho_2)$. Speed of sound:
$$\frac{1}{c_{pT}^2}=\frac{1}{c_p^2}+\frac{\rho TC_{p1}C_{p2}}{C_{p1}+C_{p2}}\Big(\frac{\Gamma_2}{\rho_2c_2^2}-\frac{\Gamma_1}{\rho_1c_1^2}\Big)^2. \tag{19}$$
Volume fraction, temperature and pressure equations:
$$\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1+\rho c_{pT}^2\Big[\alpha_1\alpha_2\Big(\frac{1}{\rho_2c_2^2}-\frac{1}{\rho_1c_1^2}\Big)+\frac{TC_{p1}C_{p2}}{C_{p1}+C_{p2}}\Big(\frac{\Gamma_2}{\rho_2c_2^2}-\frac{\Gamma_1}{\rho_1c_1^2}\Big)\Big(\frac{\alpha_1\Gamma_2}{\rho_2c_2^2}+\frac{\alpha_2\Gamma_1}{\rho_1c_1^2}\Big)\Big]\nabla\cdot\vec u=\mathcal MS_\alpha \tag{20}$$
(verified on a 200-dpi render)
$$\partial_tT+\vec u\cdot\nabla T+\frac{\rho c_{pT}^2T}{C_{p1}+C_{p2}}\Big(\frac{C_{p1}\Gamma_1}{\rho_1c_1^2}+\frac{C_{p2}\Gamma_2}{\rho_2c_2^2}\Big)\nabla\cdot\vec u=\mathcal MS_T \tag{21}$$
$$\partial_tp+\vec u\cdot\nabla p+\rho c_{pT}^2\nabla\cdot\vec u=\mathcal MS_p \tag{22}$$
with
$$S_\alpha=\frac{1}{D_T}\Big[\Big(\frac{\alpha_1}{\Gamma_1}+\frac{\alpha_2}{\Gamma_2}\Big)(\alpha_1\phi_1+\alpha_2\phi_2)+\alpha_1\alpha_2\Big(\frac{\chi_1}{\Gamma_1}-\frac{\chi_2}{\Gamma_2}\Big)(\phi_1\zeta_2-\phi_2\zeta_1)\Big] \tag{23a}$$
$$S_T=\frac{1}{D_T}\Big[\Big(\frac{\chi_2}{\Gamma_2}-\frac{\chi_1}{\Gamma_1}\Big)(\alpha_1\zeta_1\rho_2+\alpha_2\zeta_2\rho_1)+\Big(\frac{\rho_1c_1^2}{\Gamma_1}-\frac{\rho_2c_2^2}{\Gamma_2}\Big)(\alpha_1\zeta_1+\alpha_2\zeta_2)+\Big(\frac{\alpha_1}{\Gamma_1}+\frac{\alpha_2}{\Gamma_2}\Big)(\rho_2-\rho_1)\Big] \tag{23b}$$
$$S_p=\frac{1}{D_T}\Big[\Big(\frac{\chi_1}{\Gamma_1}-\frac{\chi_2}{\Gamma_2}\Big)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1)+\Big(\frac{\rho_2c_2^2}{\Gamma_2}-\frac{\rho_1c_1^2}{\Gamma_1}\Big)(\alpha_1\phi_1+\alpha_2\phi_2)\Big] \tag{23c}$$
$$D_T=\alpha_1\alpha_2\Big(\frac{\rho_1c_1^2}{\Gamma_1}-\frac{\rho_2c_2^2}{\Gamma_2}\Big)(\phi_1\zeta_2-\phi_2\zeta_1)+\Big(\frac{\alpha_1}{\Gamma_1}+\frac{\alpha_2}{\Gamma_2}\Big)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1). \tag{23d}$$
Derivation in Appendix B. (23a)–(23c) verified on a 200-dpi render.

### 3.3 Three-equation *pTg*-relaxed model (HEM; $\mu,\vartheta,\nu\to+\infty$)

$$\partial_t\rho+\nabla\cdot(\rho\vec u)=0,\quad \partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u+p\mathbb I)=0,\quad \partial_tE+\nabla\cdot((E+p)\vec u)=0. \tag{24a–c}$$
$p=p(\mathcal E,\rho)$ from (12) + $T_1(p,\rho_1)=T_2(p,\rho_2)$ + $g_1(p,T)=g_2(p,T)$. Sound speed:
$$\frac{1}{c_{pTg}^2}=\frac{1}{c_{pT}^2}+\frac{\rho T}{C_{p1}+C_{p2}}\Big[\frac{\Gamma_1C_{p1}}{\rho_1c_1^2}+\frac{\Gamma_2C_{p2}}{\rho_2c_2^2}-\frac1T\Big(\frac{dT}{dp}\Big)_{\rm sat}(C_{p1}+C_{p2})\Big]^2. \tag{25}$$
Subcharacteristic ordering:
$$c_{pTg}\le c_{pT}\le c_p\le c_f. \tag{26}$$

---

## 4. Equations of state (Section 4, p. 5–6)

Both phases are always present; pure phase $k$ approximated by $\alpha_k=1-\epsilon$, $0<\epsilon\ll1$. Each EOS must be valid over the whole thermodynamic domain of the problem (admissible $\rho_k,\mathcal E_k>0$ [as written], $\kappa_{vk}=T_k(\partial s_k/\partial T_k)_{\rho_k}>0$, $(\partial p/\partial v)_{T_k}<0$, $v=1/\rho$). Single liquid–vapour EOS (Peng–Robinson, GERG-2008, IAPWS-IF97) have a spinodal instability region; usable only if domain covered by the two-phase mixture domain of each EOS. The kinetic approach used here avoids unstable regions.

### 4.1 Mie–Grüneisen
$$p(\mathcal E,\rho)=\Gamma(\rho)(\mathcal E-\rho\varepsilon_r(\rho))+p_r(\rho). \tag{27}$$
For pressure-equilibrium mixtures with MG phases, explicit mixture pressure from (12):
$$p(\mathcal E,\rho_1,\rho_2,\alpha_1)=\frac{\mathcal E-\big(\alpha_1\rho_1\varepsilon_{r1}(\rho_1)+\alpha_2\rho_2\varepsilon_{r2}(\rho_2)\big)+\Big(\alpha_1\frac{p_{r1}(\rho_1)}{\Gamma_1(\rho_1)}+\alpha_2\frac{p_{r2}(\rho_2)}{\Gamma_2(\rho_2)}\Big)}{\frac{\alpha_1}{\Gamma_1(\rho_1)}+\frac{\alpha_2}{\Gamma_2(\rho_2)}}. \tag{28}$$
Important: no implicit solve for the pressure.

### 4.2 JWL
$\Gamma(\rho)=\Gamma_0$ (29a);
$\varepsilon_r(\rho)=\frac{a}{r_1\rho_0}e^{-r_1\rho_0/\rho}+\frac{b}{r_2\rho_0}e^{-r_2\rho_0/\rho}-\varepsilon_0$ (29b);
$p_r(\rho)=a\,e^{-r_1\rho_0/\rho}+b\,e^{-r_2\rho_0/\rho}$ (29c). Parameters $\Gamma_0,\rho_0,r_1,r_2,\varepsilon_0$ (and $a,b$).

### 4.3 Noble–Abel stiffened gas (NASG) (Le Métayer & Saurel 2016)
$$p(\mathcal E,\rho)=\frac{\gamma-1}{1-\rho b}(\mathcal E-\eta\rho)-\gamma\varpi, \tag{30a}$$
$$T(p,\rho)=\frac{(1-\rho b)(p+\varpi)}{\kappa_v\rho(\gamma-1)}. \tag{30b}$$
$\gamma,\varpi,\eta,b,\kappa_v$ constants; $b$ = covolume ($b=0$: stiffened gas). MG form with
$$\Gamma(\rho)=\frac{\gamma-1}{1-\rho b},\quad \varepsilon_r=\eta,\quad p_r=-\gamma\varpi. \tag{31}$$
$$s(p,T)=\kappa_v\log\frac{T^\gamma}{(p+\varpi)^{\gamma-1}}+\bar\eta, \tag{32a}$$
$$h(T,p)=\kappa_pT+bp+\eta, \tag{32b}$$
$$g(p,T)=(\gamma\kappa_v-\bar\eta)T-\kappa_vT\log\frac{T^\gamma}{(p+\varpi)^{\gamma-1}}+\eta+bp, \tag{32c}$$
with $\kappa_p=\gamma\kappa_v$; $\bar\eta$ constant.
$$c=\sqrt{\gamma\frac{p+\varpi}{\rho(1-\rho b)}} \tag{33}$$
$$\phi=-\frac{\rho}{T}(1-\rho b),\qquad \zeta=\frac{1}{T(\gamma-1)\kappa_v+b(p+\varpi)}(1-\rho b). \tag{34}$$
*Note:* (34) follows from $\rho=(p+\varpi)/(\kappa_v(\gamma-1)T+b(p+\varpi))$.

#### 4.3.1 Saturation curve
Equilibrium $p_1=p_2=p$, $T_1=T_2=T$, $g_1(p,T)=g_2(p,T)$ gives
$$A_s+\frac{B_s}{T}+C_s\log T+D_s\log(p+\varpi_1)-\log(p+\varpi_2)+\frac{pE_s}{T}=0, \tag{35}$$
$$A_s=\frac{\kappa_{p1}-\kappa_{p2}-\bar\eta_1+\bar\eta_2}{\kappa_{p2}-\kappa_{v2}},\ B_s=\frac{\eta_1-\eta_2}{\kappa_{p2}-\kappa_{v2}},\ C_s=\frac{\kappa_{p2}-\kappa_{p1}}{\kappa_{p2}-\kappa_{v2}},\ D_s=\frac{\kappa_{p1}-\kappa_{v1}}{\kappa_{p2}-\kappa_{v2}},\ E_s=\frac{b_1-b_2}{\kappa_{p2}-\kappa_{v2}}. \tag{36}$$
NASG parameters fitted so the theoretical saturation curves match experiments in a temperature range (Le Métayer & Saurel 2016). Tables 3–5 (Section 8): water parameters from LMS2016; Table 4: slightly modified dodecane set.

---

## 5. Numerical method (Section 5, p. 6–7)

Compact form, $q\in\mathbb R^{5+d}$:
$$\partial_tq+\nabla\cdot f(q)+\sigma(q,\nabla q)=\psi_\mu(q)+\psi_\vartheta(q)+\psi_\nu(q), \tag{37a}$$
$$q=\begin{bmatrix}\alpha_1\\\alpha_1\rho_1\\\alpha_2\rho_2\\\rho\vec u\\\alpha_1E_1\\\alpha_2E_2\end{bmatrix},\
f(q)=\begin{bmatrix}0\\\alpha_1\rho_1\vec u\\\alpha_2\rho_2\vec u\\\rho\vec u\otimes\vec u+(\alpha_1p_1+\alpha_2p_2)\mathbb I\\\alpha_1(E_1+p_1)\vec u\\\alpha_2(E_2+p_2)\vec u\end{bmatrix},\
\sigma(q,\nabla q)=\begin{bmatrix}\vec u\cdot\nabla\alpha_1\\0\\0\\0\\\Sigma\\-\Sigma\end{bmatrix} \tag{37b}$$
$$\psi_\mu=\begin{bmatrix}\mathcal P\\0\\0\\0\\-p_I\mathcal P\\p_I\mathcal P\end{bmatrix},\
\psi_\vartheta=\begin{bmatrix}0\\0\\0\\0\\\mathcal Q\\-\mathcal Q\end{bmatrix},\
\psi_\nu=\begin{bmatrix}0\\\mathcal M\\-\mathcal M\\0\\(g_I+\frac{|\vec u|^2}{2})\mathcal M\\-(g_I+\frac{|\vec u|^2}{2})\mathcal M\end{bmatrix}. \tag{37c}$$

Characteristic times $\tau_\mu,\tau_\vartheta,\tau_\nu$ with underlying assumption $\tau_\mu\ll\tau_\vartheta\ll\tau_\nu$ (Kapila et al. 2001). Algorithm per time step:

1. **Homogeneous hyperbolic system** $\partial_tq+\nabla\cdot f(q)+\sigma(q,\nabla q)=0$ (38). Output: superscript $^0$.
2. **Relaxation steps:**
   - **2(a)** Instantaneous mechanical relaxation: solve $\partial_tq=\psi_\mu(q)$ (39) in the limit $\mu\to+\infty$. Output: $^*$.
   - **2(b)** Thermal relaxation: $\partial_tq=\psi_\mu(q)+\psi_\vartheta(q)$ (40) with $\mu\to+\infty$ — drives toward thermal equilibrium while maintaining pressure equilibrium. Output: $^{**}$.
   - **2(c)** Chemical relaxation: $\partial_tq=\psi_\mu(q)+\psi_\vartheta(q)+\psi_\nu(q)$ (41) with $\mu\to+\infty$ — toward full thermodynamic equilibrium, maintaining pressure equilibrium. Output: $^\otimes$.

Step 2(a) always activated (mechanical equilibrium flows). Steps 2(b), 2(c) activated or not depending on the problem, possibly only at selected locations, typically interfaces identified by $\min_k\alpha_k>\epsilon$ (e.g. $\epsilon=10^{-6}$). If thermal and chemical relaxation are unconditionally instantaneous, the scheme approximates the *pTg*-relaxed model (24).

*Note (for TV-split IMEX users):* the paper's splitting is first-order Godunov (Lie) splitting: hyperbolic → p → (p,T) → (p,T,g), applied once per step with $\Delta t$; no Strang splitting is mentioned [to be confirmed in later pages].

### 5.1 Mixture-energy-consistency (Pelanti & Shyue 2014)
With $^\#\in\{*,**,\otimes\}$ and $E^{0,C}$ the discrete mixture total energy from a conservative discretisation of (7), the scheme is mixture-energy-consistent if:

(i) mixture total energy conservation consistency:
$$E^0=E^{0,C},\qquad E^0=(\alpha_1E_1)^0+(\alpha_2E_2)^0; \tag{42}$$
(ii) relaxed pressure consistency (consistency with mixture pressure law (12)):
$$\mathcal E^0=\alpha_1^\#\mathcal E_1(p^\#,\rho_1^\#)+\alpha_2^\#\mathcal E_2(p^\#,\rho_2^\#),\qquad \mathcal E^0=E^0-\frac{(\rho\vec u)^0\cdot(\rho\vec u)^0}{2\rho^0}. \tag{43}$$
(42) is guaranteed by the HLLC method.

---

## 6. Homogeneous system solver (Section 6, p. 7)

1D ($d=1$), $q\in\mathbb R^6$ (setting $\vec u=u$, $\nabla=\partial_x$). Second-order wave-propagation scheme (LeVeque 1997, 2002):
$$Q_i^{n+1}=Q_i^n-\frac{\Delta t}{\Delta x}\big(\mathcal A^+\Delta Q_{i-1/2}+\mathcal A^-\Delta Q_{i+1/2}\big)-\frac{\Delta t}{\Delta x}\big(F^h_{i+1/2}-F^h_{i-1/2}\big). \tag{44}$$
Riemann solution: $\mathcal M$ waves $\mathcal W^l$ with speeds $s^l$ ($\mathcal M\gtrless3N$ in general [the paper's symbol; for a two-phase HLLC, 3 waves]).
$$\mathcal A^\pm\Delta Q_{i+1/2}=\sum_{l=1}^{\mathcal M}(s^l_{i+1/2})^\pm\mathcal W^l_{i+1/2},\quad s^+=\max(s,0),\ s^-=\min(s,0) \tag{45}$$
$$F^h_{i+1/2}=\frac12\sum_{l=1}^{\mathcal M}|s^l_{i+1/2}|\Big(1-\frac{\Delta t}{\Delta x}|s^l_{i+1/2}|\Big)\widetilde{\mathcal W}^l_{i+1/2} \tag{46}$$
$\widetilde{\mathcal W}$ = limited waves (limiter function, LeVeque 2002). Riemann solver: simple HLLC-type solver of Pelanti & Shyue (2014), Appendix C.

---

## 7. Relaxation processes (Section 7)

For any relaxation process:
$$\partial_t\rho=0,\quad\partial_t(\rho\vec u)=0,\quad\partial_tE=0 \tag{47a–c}$$
$$\Rightarrow\ \rho,\ \vec u,\ E,\ \mathcal E=\text{const.} \tag{48}$$
If chemical relaxation is not active, also $\partial_t(\alpha_k\rho_k)=0$:
$$\alpha_k\rho_k=\text{const.} \tag{49}$$
Mechanical and thermal steps: two independent unknowns ($\alpha_1$ and equilibrium $p$). Chemical step: three (partial densities also vary). For consistency with (12), the equilibrium pressure $p$ in **all** relaxation steps must satisfy (43).

### 7.1 Instantaneous relaxation (recap of Pelanti & Shyue 2014)
Use invariants (48) (+ (49) for mechanical/thermal) and equilibrium conditions to obtain an algebraic system:

- **2(a)** mechanical: (48), (49) + $p_1=p_2=p$. Also integrate the phasic energy equations $\partial_t\mathcal E_k=(-1)^kp_I\partial_t\alpha_1$, $k=1,2$, between states $^0$ and $^*$, with assumption $p_I=\beta_pp^0+(1-\beta_p)p^*$, $\beta_p\in[0,1]$ (P&S 2014: $\beta_p=\tfrac12$; Pelanti 2017: $\beta_p=0$). Algebraic system for e.g. $\alpha_1^*$, $p^*$. For SG EOS reduces to a quadratic in $p^*$. [*Note:* $p^0$ here presumably a phase-dependent $p_k^0$ — the paper writes "$p_I=\beta_pp^0+(1-\beta_p)p^*$" — uncertain whether $p^0$ means $p_I^0$.]
- **2(b)** mechanical+thermal: (48), (49) + $p_1=p_2=p$ + $T_1=T_2=T$. Unknowns e.g. $\alpha_1^{**}$, $p^{**}$. SG: quadratic for $p^{**}$.
- **2(c)** mechanical+thermal+chemical: (48) + $p_1=p_2=p$ + $T_1=T_2=T$ + $g_1=g_2$. Unknowns e.g. $\alpha_1^\otimes,p^\otimes,T^\otimes$. In general iterative.

All these satisfy (43) by construction → mixture-energy-consistent.

### 7.2 Arbitrary-rate relaxation processes

Mechanical relaxation is always instantaneous, but a new semi-exact exponential procedure is given for it too (advantageous for complex EOS; potentially usable for finite-rate pressure relaxation, e.g. Pandare et al. 2020 — not pursued).

**Core idea.** Thermal and chemical relaxation must preserve pressure equilibrium ($\psi_\mu$ with $\mu\to\infty$ present in (40), (41)), and chemical relaxation with $\vartheta\to\infty$ must keep thermal equilibrium. So a simple fractional step integrating each source separately is not possible, and solving all relaxation terms together (very different time scales) is complicated. Instead: replace (40) by $\partial_tq=\bar\psi_\vartheta$ where $\bar\psi_\vartheta$ models thermal relaxation **under pressure equilibrium**, and (41) by $\partial_tq=\bar\psi_\nu$ modelling chemical relaxation under pressure equilibrium and a thermal constraint (in particular thermal equilibrium). These are built from the equations of the *p*-relaxed (9) and *pT*-relaxed (18) models. With simplifying assumptions, analytical semi-exact exponential solutions follow.

Efficiency: each step uses explicit expressions for the independent variables of relaxed states **except for the equilibrium pressure**, which is defined implicitly by (12) (iterative in general; explicit for MG-type EOS via (28)). Generally cheaper than the instantaneous procedures of 7.1; even for SG, the 7.1 chemical relaxation needs iteration, while the new technique uses only explicit formulas for SG.

#### 7.2.1 Mechanical relaxation (new procedure)
ODEs from (1a) and (8) (only $\mathcal P$ active):
$$\partial_t\alpha_1=\mu(p_1-p_2), \tag{50a}$$
$$\partial_tp_1=\mu\frac{1}{\alpha_1}\big[\Gamma_1(\mathcal E_1+p_I)+\chi_1\rho_1\big](p_2-p_1), \tag{50b}$$
$$\partial_tp_2=-\mu\frac{1}{\alpha_2}\big[\Gamma_2(\mathcal E_2+p_I)+\chi_2\rho_2\big](p_2-p_1), \tag{50c}$$
$\mu>0$ constant; initial condition = homogeneous solution (superscript $^0$). Define
$$\frac{1}{\xi_k^p}=\frac{1}{\alpha_k}\big[\Gamma_k(\mathcal E_k+p_I)+\chi_k\rho_k\big]=\frac{1}{\alpha_k}\big[\Gamma_k(p_I-p_k)+\rho_kc_k^2\big],\quad k=1,2. \tag{51}$$
Approximation: $\xi_k^p$ constant in time, $\xi_k^p=\xi_k^{p0}$:
$$\partial_tp_1=\mu\frac{1}{\xi_1^{p0}}(p_2-p_1),\qquad \partial_tp_2=-\mu\frac{1}{\xi_2^{p0}}(p_2-p_1). \tag{52a,b}$$
For $\Delta p=p_2-p_1$:
$$\partial_t\Delta p=-\mu\Big(\frac{1}{\xi_1^{p0}}+\frac{1}{\xi_2^{p0}}\Big)\Delta p, \tag{53}$$
$$\Delta p^*=\Delta p^0e^{-K_p\Delta t}, \tag{54}$$
$$K_p=\mu\Big(\frac{1}{\xi_1^{p0}}+\frac{1}{\xi_2^{p0}}\Big). \tag{55}$$
$K_p>0$ assumed (physical: compressed phase with lowest pressure, expanded phase with highest pressure, $\mu>0$). Integrate (50a) with (54):
$$\alpha_1^*=\alpha_1^0-\frac{\Delta p^0}{\frac{1}{\xi_1^{p0}}+\frac{1}{\xi_2^{p0}}}\big(1-e^{-K_p\Delta t}\big). \tag{56}$$
Limit $\mu\to+\infty$:
$$\alpha^*_{1,\mu\to\infty}=\alpha_1^0-\frac{\Delta p^0}{\frac{1}{\xi_1^{p0}}+\frac{1}{\xi_2^{p0}}}, \tag{57}$$
$$p^*_{\mu\to\infty}=\frac{\xi_1^{p0}p_1^0+\xi_2^{p0}p_2^0}{\xi_1^{p0}+\xi_2^{p0}}. \tag{58}$$
Since the pair (57)–(58) does **not** in general satisfy (43) (because of the approximations), the procedure: update $\alpha_1^*$ by (57) [the exponential solution], and **recompute the equilibrium pressure $p^*$ from the energy relation**
$$\mathcal E^0=\alpha_1^*\mathcal E_1\Big(p^*,\frac{(\alpha_1\rho_1)^0}{\alpha_1^*}\Big)+\alpha_2^*\mathcal E_2\Big(p^*,\frac{(\alpha_2\rho_2)^0}{\alpha_2^*}\Big). \tag{59}$$
Convenient with complex EOS (P&S 2014 procedure may need an iterative method). Appendix D shows that this procedure ensures velocity and pressure invariance at material interfaces at least with SG EOS [continued next page].

*Note:* $p_I$ in (51) is evaluated at state $^0$ (with impedance-weighted $p_I$ in tests) [uncertain whether the paper uses $p_I^0$ or another value — see later pages].
used [i.e. the SG EOS]. A similar mechanical relaxation procedure was proposed in De Lorenzo et al. (2019), but there the exponential solution for $\alpha_1$ was assumed *a priori*, whereas here it is deduced from the model equations by assuming $\xi_k^p$ constant during relaxation. Also, in De Lorenzo 2019 the instantaneous case was an exponential with very small relaxation time; here the analytical limit (57) is used.

**Summary of the mechanical step (as I read it):**
1. From state $^0$ compute $\xi_k^{p0}$ via (51) (needs $\Gamma_k,\rho_kc_k^2,p_k,p_I$ at $^0$).
2. $\alpha_1^*$ from (57) (or (56) for finite $\mu$).
3. $\rho_k^*=(\alpha_k\rho_k)^0/\alpha_k^*$; solve (59) for $p^*$ (explicit via (28) for MG EOS).
4. Phasic energies: $\mathcal E_k^*=\mathcal E_k(p^*,\rho_k^*)$, so $(\alpha_kE_k)^*=\alpha_k^*\mathcal E_k^*+\tfrac12(\alpha_k\rho_k)^0|\vec u|^2$. [*Note:* the paper does not spell out step 4 explicitly in Section 7.2.1; it is implied by the relaxed state being defined by $\alpha_1^*,p^*$ and the invariants.]
5. If $\alpha_1^*\notin(0,1)$ (unphysical), set $\alpha_1=\bar\alpha_1$ (see "complete evaporation" below).

#### 7.2.2 Thermal relaxation (finite or infinite $\vartheta$, under $p_1=p_2$)

Built from the 5-eq. *p*-relaxed model (9). From (9a) and (15) with only $\mathcal Q$:
$$\partial_t\alpha_1=\vartheta\frac{Z}{D}(T_2-T_1), \tag{60a}$$
$$\partial_tT_1=\vartheta\frac{1}{\phi_1D}\Big[-\frac{\rho_1}{\alpha_1}Z-\zeta_1(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big](T_2-T_1), \tag{60b}$$
$$\partial_tT_2=\vartheta\frac{1}{\phi_2D}\Big[\frac{\rho_2}{\alpha_2}Z-\zeta_2(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big](T_2-T_1), \tag{60c}$$
$D$ from (10) and
$$Z=\alpha_2\Gamma_1+\alpha_1\Gamma_2. \tag{61}$$
Initial condition: pressure-equilibrium state from the mechanical step ($^*$). Define
$$\frac{1}{\xi_k^T}=-\frac{\frac{\rho_k}{\alpha_k}Z+(-1)^{k-1}\zeta_k(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)}{\phi_kD},\quad k=1,2. \tag{62}$$
(Checked: consistent with (60b,c) given $\partial_tT_1=\vartheta\xi_1^{-1}\Delta T$, $\partial_tT_2=-\vartheta\xi_2^{-1}\Delta T$.)
Freeze $\xi_k^T=\xi_k^{T*}$:
$$\partial_tT_1=\vartheta\frac{1}{\xi_1^{T*}}(T_2-T_1),\qquad \partial_tT_2=-\vartheta\frac{1}{\xi_2^{T*}}(T_2-T_1). \tag{63a,b}$$
$\vartheta>0$ constant; **if $\vartheta$ is a function of variables, freeze it at state $^*$.** For $\Delta T=T_2-T_1$:
$$\partial_t\Delta T=-\vartheta\Big(\frac{1}{\xi_1^{T*}}+\frac{1}{\xi_2^{T*}}\Big)\Delta T, \tag{64}$$
$$\Delta T^{**}=\Delta T^*e^{-K_T\Delta t}, \tag{65}$$
$$K_T=\vartheta\Big(\frac{1}{\xi_1^{T*}}+\frac{1}{\xi_2^{T*}}\Big). \tag{66}$$
$K_T>0$ assumed (heat flows from hotter to colder phase if $\vartheta>0$). Assuming $Z/D$ constant, integrating (60a):
$$\alpha_1^{**}=\alpha_1^*+\Big(\frac ZD\Big)^*\frac{\Delta T^*}{\frac{1}{\xi_1^{T*}}+\frac{1}{\xi_2^{T*}}}\big(1-e^{-K_T\Delta t}\big). \tag{67}$$
Instantaneous limit $\vartheta\to+\infty$:
$$\alpha_{1,\vartheta\to\infty}=\alpha_1^*+\Big(\frac ZD\Big)^*\frac{\Delta T^*}{\frac{1}{\xi_1^{T*}}+\frac{1}{\xi_2^{T*}}}, \tag{68}$$
$$T_{\vartheta\to\infty}=\frac{\xi_1^{T*}T_1^*+\xi_2^{T*}T_2^*}{\xi_1^{T*}+\xi_2^{T*}}. \tag{69}$$
Update: use $\alpha_1^{**}$ from (67) (or (68)), and $p^{**}$ from the energy relation
$$\mathcal E^0=\alpha_1^{**}\mathcal E_1\Big(p^{**},\frac{(\alpha_1\rho_1)^0}{\alpha_1^{**}}\Big)+\alpha_2^{**}\mathcal E_2\Big(p^{**},\frac{(\alpha_2\rho_2)^0}{\alpha_2^{**}}\Big) \tag{70}$$
to ensure mixture-energy-consistency. *Note:* (69) is not used to set the state; the temperatures after the step follow from $T_k(p^{**},\rho_k^{**})$, so in the instantaneous limit $T_1^{**}=T_2^{**}$ only approximately (the paper's figures show the temperatures overlap to plotting accuracy).

*Note:* with $\phi_k<0$ (normal fluids, $\beta_k>0$), the sign of $\xi_k^T$ is the sign of $[\rho_kZ/\alpha_k\pm\dots]$; $K_T>0$ is an *assumption*. Implementation should guard $K_T\le0$ (not discussed explicitly in the paper; only $\alpha_1$ outside $(0,1)$ is treated, see below).

#### 7.2.3 Chemical relaxation (finite or infinite $\nu$, under $p_1=p_2$ and $T_1=T_2$)

Chemical relaxation is coupled with a thermal relaxation assumed much faster, so it is modelled under both mechanical and thermal equilibrium, using the ODEs of the 4-eq. *pT*-relaxed model (18). Unknowns: $\alpha_1$, $\alpha_1\rho_1$, $g_1$, $g_2$. From (18) and (23): $\partial_t\alpha_1=S_\alpha\mathcal M$, $\partial_t(\alpha_1\rho_1)=\mathcal M$; and $\partial_tg_k=\frac{1}{\rho_k}\partial_tp_k-s_k\partial_tT_k$ (Gibbs relation $dg=v\,dp-s\,dT$) with $\partial_tT=S_T\mathcal M$, $\partial_tp=S_p\mathcal M$:
$$\partial_t\alpha_1=\nu S_\alpha(g_2-g_1), \tag{71a}$$
$$\partial_tg_1=\nu\Big(\frac{S_p}{\rho_1}-s_1S_T\Big)(g_2-g_1), \tag{71b}$$
$$\partial_tg_2=\nu\Big(\frac{S_p}{\rho_2}-s_2S_T\Big)(g_2-g_1), \tag{71c}$$
$$\partial_t(\alpha_1\rho_1)=\nu(g_2-g_1). \tag{71d}$$
**$g_I$ does not appear** (as in the *pT*-relaxed model), so no definition is needed. Initial condition: state $^{**}$. Define
$$\frac{1}{\xi_k^g}=\frac{(-1)^{k-1}}{\rho_k}S_p-s_kS_T,\quad k=1,2. \tag{72 as printed}$$
**[uncertain / probable typo]**: consistency with (71b,c) and with (73) ($\partial_tg_2=-\nu\xi_2^{-1}\Delta g$) requires $\frac{1}{\xi_k^g}=(-1)^{k-1}\Big(\frac{S_p}{\rho_k}-s_kS_T\Big)$, i.e. $\frac1{\xi_2^g}=-\frac{S_p}{\rho_2}+s_2S_T$. The printed (72) has $-s_2S_T$ for $k=2$. I recommend the version derived from (71): then
$$K_g=\nu\Big[S_p\Big(\frac1{\rho_1}-\frac1{\rho_2}\Big)-(s_1-s_2)S_T\Big].$$
Freeze $\xi_k^g=\xi_k^{g**}$:
$$\partial_tg_1=\nu\frac{1}{\xi_1^{g**}}(g_2-g_1),\qquad \partial_tg_2=-\nu\frac{1}{\xi_2^{g**}}(g_2-g_1). \tag{73a,b}$$
$\nu>0$ constant; **if $\nu$ is a function, freeze at state $^{**}$.** For $\Delta g=g_2-g_1$:
$$\partial_t\Delta g=-\nu\Big(\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}\Big)\Delta g, \tag{74}$$
$$\Delta g^\otimes=\Delta g^{**}e^{-K_g\Delta t}, \tag{75}$$
$$K_g=\nu\Big(\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}\Big). \tag{76}$$
$K_g>0$ assumed (mass goes from phase with highest chemical potential to lowest if $\nu>0$). Integrating (71d) and (71a) (with $S_\alpha$ constant $=S_\alpha^{**}$):
$$(\alpha_1\rho_1)^\otimes=(\alpha_1\rho_1)^{**}+\frac{\Delta g^{**}}{\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}}\big(1-e^{-K_g\Delta t}\big), \tag{77}$$
$$\alpha_1^\otimes=\alpha_1^{**}+S_\alpha^{**}\frac{\Delta g^{**}}{\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}}\big(1-e^{-K_g\Delta t}\big). \tag{78}$$
Instantaneous limit $\nu\to+\infty$:
$$(\alpha_1\rho_1)_{\nu\to\infty}=(\alpha_1\rho_1)^{**}+\frac{\Delta g^{*}}{\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}}, \tag{79}$$
$$\alpha_{1,\nu\to\infty}=\alpha_1^{*}+S_\alpha^{**}\frac{\Delta g^{*}}{\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}}, \tag{80}$$
**[probable typos in the paper: (79)–(80) print $\Delta g^*$ and $\alpha_1^*$; by consistency with (77)–(78) these should be $\Delta g^{**}$ and $\alpha_1^{**}$.]**
$$g_{\nu\to\infty}=\frac{\xi_1^{g**}g_1^{**}+\xi_2^{g**}g_2^{**}}{\xi_1^{g**}+\xi_2^{g**}}. \tag{81}$$
Update: $(\alpha_1\rho_1)^\otimes$, $\alpha_1^\otimes$ from (77),(78) (or (79),(80)); $(\alpha_2\rho_2)^\otimes=\rho^0-(\alpha_1\rho_1)^\otimes$ [implied by (48)]; $p^\otimes$ from
$$\mathcal E^0=\alpha_1^\otimes\mathcal E_1\Big(p^\otimes,\frac{(\alpha_1\rho_1)^\otimes}{\alpha_1^\otimes}\Big)+\alpha_2^\otimes\mathcal E_2\Big(p^\otimes,\frac{(\alpha_2\rho_2)^\otimes}{\alpha_2^\otimes}\Big). \tag{82}$$
Remark: (71) was derived assuming $T_1=T_2$, but holds more generally for chemical relaxation at **constant temperature difference** $T_2-T_1=\Delta T$ (since then $\partial T_1=\partial T_2$, and the derivation of (23) in Appendix B extends with phase-$k$ variables functions of $p$ and $T_k$, $T_2=T_1+\Delta T$). So the procedure can model chemical relaxation with thermal disequilibrium, provided consistent relaxation times are chosen (thermal relaxation physically faster than chemical).

*Note:* the coefficients $S_\alpha,S_T,S_p$ (23) at state $^{**}$ need $\Gamma_k,\chi_k,\phi_k,\zeta_k,\rho_kc_k^2,\rho_k,\alpha_k$; $s_k$ also needed for (72). For finite-rate $\vartheta$ followed by chemical relaxation, the state $^{**}$ is not in thermal equilibrium; the paper's remark above justifies using (23) with each phase's own $T_k$ [my interpretation].

#### Complete evaporation or condensation
Mass transfer may make one phase disappear; then (78) (or (80)) gives $\alpha_1^\otimes\notin(0,1)$. Both phases must always be present: a pure phase $k$ is modelled by $\alpha_k=1-\epsilon_\alpha$, $0<\epsilon_\alpha\ll1$ (e.g. $\epsilon_\alpha=10^{-8}$). If $\alpha_1^\otimes\notin(0,1)$: one phase vanishes; set $\alpha_1^\otimes=\bar\alpha_1$ with $\bar\alpha_1=\epsilon_\alpha$ if $s_2^{**}>s_1^{**}$ (vanishing of phase 1), or $\bar\alpha_1=1-\epsilon_\alpha$ if $s_1^{**}>s_2^{**}$ (vanishing of phase 2). [*Note:* $s_k$ here is the phasic **entropy**, used to identify the vapour (higher entropy) as the phase that survives on evaporation — as printed; the logic presumably: if $\alpha_1$ overshoots, the phase being destroyed is identified from entropy ordering. Slightly ambiguous; an implementation can equivalently choose based on sign of $\alpha_1^\otimes$: $<0\Rightarrow\epsilon_\alpha$, $>1\Rightarrow1-\epsilon_\alpha$.] Then determine the relaxation parameter giving $\bar\alpha_1$ through (78):
$$\bar K_g=-\frac{1}{\Delta t}\log\left(1-\frac{\bar\alpha_1-\alpha_1^{**}}{S_\alpha^{**}\dfrac{\Delta g^{**}}{\frac{1}{\xi_1^{g**}}+\frac{1}{\xi_2^{g**}}}}\right). \tag{83}$$
Then set $(\alpha_1\rho_1)^\otimes$ via (77) with $K_g=\bar K_g$ (keeps mass and volume updates consistent). If an unphysical $\alpha_1$ is computed in the pressure or temperature relaxation procedures, simply set $\alpha_1=\bar\alpha_1$.

*Note:* (83) requires the argument of log in $(0,1]$; it is when $\bar\alpha_1$ lies between $\alpha_1^{**}$ and the overshooting $\alpha_{1,\nu\to\infty}$. In the instantaneous case this replaces the infinite rate with a finite effective rate $\bar K_g$.

---

## 8. Numerical experiments (Section 8)

Implemented with the basic routines of **CLAWPACK** (LeVeque). All computations: second-order wave-propagation scheme with **minmod limiter**. For some 1D Riemann problems, exact solutions of the *p*-relaxed model (9) (Petitpas et al. 2007) and of the *pT*-relaxed model (18) are plotted, computed by extending Kamm (2015) and Israel et al. (2014) (single-phase Euler, general EOS).

### EOS parameter tables

**Table 2 — JWL (explosive / detonation gases):** $\rho_0=1840$ kg/m³, $\Gamma_0=0.25$, $a=854.5\times10^9$ Pa, $b=20.5\times10^9$ Pa, $r_1=4.6$, $r_2=1.35$, $\varepsilon_0=0$ J/kg. (Taken from Saurel & Abgrall 1999.)

**Table 3 — NASG water, 300–500 K:**

| Phase | $\gamma$ | $\varpi$ [Pa] | $\eta$ [J/kg] | $\bar\eta$ [J/(kg K)] | $\kappa_v$ [J/(kg K)] | $b$ [m³/kg] |
|---|---|---|---|---|---|---|
| Liquid | 1.187 | $7028\times10^5$ | $-1\,177\,788$ | 0 | 3610 | $6.61\times10^{-4}$ |
| Vapor | 1.467 | 0 | $2\,077\,616$ | 14 317 | 955 | 0 |

**Table 4 — NASG dodecane, 400–600 K (slightly modified from LMS2016):**

| Phase | $\gamma$ | $\varpi$ [Pa] | $\eta$ [J/kg] | $\bar\eta$ [J/(kg K)] | $\kappa_v$ [J/(kg K)] | $b$ [m³/kg] |
|---|---|---|---|---|---|---|
| Liquid | 1.206 | $1681\times10^5$ | $-996\,054$ | 0 | 2532 | $7.51\times10^{-4}$ |
| Vapor | 1.021 | 0 | $-384\,592$ | $-4301$ | 2274 | 0 |

**Table 5 — NASG water, 350–550 K:**

| Phase | $\gamma$ | $\varpi$ [Pa] | $\eta$ [J/kg] | $\bar\eta$ [J/(kg K)] | $\kappa_v$ [J/(kg K)] | $b$ [m³/kg] |
|---|---|---|---|---|---|---|
| Liquid | 1.387 | $8899\times10^5$ | $-1\,244\,191$ | 0 | 3202 | $4.78\times10^{-4}$ |
| Vapor | 1.954 | 0 | $2\,287\,484$ | 6417 | 462 | 0 |

### 8.1 Tests with only (instantaneous) mechanical relaxation

#### 8.1.1 Detonation gas–water shock tube (1D)
- Domain $[-10,10]$ m, discontinuity at $x=0$.
- Left: detonation gases, $\rho_{\rm gas}=2000$ kg/m³, $p_{\rm gas}=4.6406\cdot10^{10}$ Pa; JWL (Table 2).
- Right: liquid water, $\rho_{\rm liq}=1044$ kg/m³, $p_{\rm liq}=10^5$ Pa; NASG (Table 3).
- Each region nearly pure: volume fraction $1-10^{-8}$. $u=0$.
- 500 cells, CFL 0.4, $t=0.9$ ms.
- Fig. 1: $\rho,u,p,\alpha_1$ vs exact solution — good agreement; velocity and pressure invariance across the material interface preserved (no oscillations). Structure: left rarefaction, contact at $x\approx2$, right shock at $x\approx8$ ($u\approx2400$ m/s, $p\approx2.2\times10^{10}$ Pa in the star region, read from figure).

#### 8.1.2 Underwater explosion close to a rigid wall (2D)
- Variant of Xie et al. (2006, 2007) (they used a one-fluid cavitation model with Tait EOS + ideal gas).
- Gas bubble (JWL, Table 2) in liquid water (NASG, Table 3) near an upper flat wall.
- Domain $[-1.35,1.35]\times[-1.5,0]$ m², wall at $y=0$. Bubble centre $(x_b,y_b)=(0,-0.32)$ m, radius $r_b=0.05$ m.
- Inside: $p=4.6406\cdot10^{10}$ Pa, $\rho_{\rm gas}=2000$ kg/m³, $\alpha_{\rm gas}=1-10^{-6}$. Outside: $p=10^5$ Pa, $\rho_{\rm liq}=1044$ kg/m³, $\alpha_{\rm liq}=1-10^{-6}$ (same states as 8.1.1).
- Grid $721\times400$, CFL 0.4. Pressure fields at $t=0,0.03,0.06,0.08,0.21,0.3$ ms (Fig. 2; plotted on $x\in[-1,1]$). Fig. 3: pressure history at $(0,0)$ until $t=0.4$ ms.
- Results: circular shock, reflection at wall, rarefactions → strong pressure decrease near the wall, moderate recompression, then another expansion. Minimum pressure at the wall centre around $t=0.21$ ms; maximum of subsequent recompression at about $t=0.3$ ms. With this setup expansion regions near the wall do not reach very low pressures (no cavitation, unlike some cases of Xie et al. 2007). Peak wall pressure in Fig. 3 ≈ $2.2\times10^{10}$ Pa [read from plot, approximate].

### 8.2 Tests with thermo-chemical relaxation
Both instantaneous and finite-rate transfer; mainly finite-rate mass transfer under thermal equilibrium.

#### 8.2.1 Water two-phase cavitation tube (1D)
- Variant of Saurel et al. (2008), previously done in Pelanti & Shyue 2014 with SG EOS and instantaneous relaxation.
- Tube $[0,1]$ m, liquid water with uniform small vapour fraction $\alpha_{\rm vap}=10^{-2}$. $p=10^5$ Pa, $T=T_1=T_2=353$ K $<T_{\rm sat}(p)=372.3$ K. Phases in thermal equilibrium but **not** chemical equilibrium.
- Velocity discontinuity at $x=0.5$ m: $u=-2$ m/s left, $u=2$ m/s right.
- NASG for liquid water and vapour, Table 3.
- Five relaxation levels (techniques of Section 7.2):
  - *p-relax*: only instantaneous mechanical;
  - *pT(f)*: instantaneous mechanical + finite-rate thermal, $\vartheta=2000$ Pa/(s K);
  - *pT*: instantaneous mechanical + thermal;
  - *pTg(f)*: instantaneous mechanical + thermal + finite-rate chemical, $\nu=10^{-4}$ Pa kg²/(s J²);
  - *pTg*: full instantaneous thermodynamic relaxation.
  Phase transition only in the last two.
- Following Saurel et al. (2008), **chemical relaxation activated only when the equilibrium temperature exceeds the saturation temperature, $T>T_{\rm sat}(p)$** (happens in the middle of the tube). Because of this and the initial chemical-potential disequilibrium, the *pTg* solution does **not** correspond to a solution of the HEM (24).
- 5000 cells, CFL 0.5, $t=0.003$ s (Fig. 4: $p$, $u$, $\alpha_2$, $Y_2$; Fig. 5: $T$ and $g_k$).
- Results: two rarefactions moving apart lower pressure in the middle and raise vapour volume fraction. Without mass transfer $Y_{\rm vap}$ constant (purely mechanical cavitation). With mass transfer two evaporation waves appear, vapour mass fraction increases in the middle and pressure is driven to its saturation value (≈ $0.48\times10^5$ Pa plateau in *pTg* from Fig. 4 [read from plot]), whereas without mass transfer pressure keeps decreasing (min ≈ $0.05\times10^5$ Pa for *p-relax*, ≈ $0.13\times10^5$ for *pT*, ≈ $0.36\times10^5$ for *pTg(f)* [read from plot]).
- Fig. 5: liquid and vapour temperatures overlap for *pT*, *pTg(f)*, *pTg* → thermal equilibrium numerically imposed (T dip ≈ 352.94 K for *pTg*). Chemical potentials: region of equality in the middle for *pTg* (activation region $T>T_{\rm sat}$); with *pTg(f)* the $g$ difference is reduced relative to *pT*.
- Good agreement of *p-relax* with exact 5-eq. solution and of *pT* with exact 4-eq. solution → the method approximates relaxed limit models.
- Leading edges of rarefactions travel at speeds consistent with subcharacteristic condition (26), $c_{pT}\le c_p$. With phase transition, leading edges move at $c_{pT}$ (19), not $c_{pTg}$ (25) (much smaller), since chemical relaxation is not active where $T<T_{\rm sat}$.

#### 8.2.2 Dodecane liquid–vapour shock tube (1D)
- From Saurel et al. (2008); solved in P&S 2014 with SG EOS.
- Unit tube, discontinuity at $x=0.75$ m. Left: liquid dodecane; right: vapour dodecane; each region contains $\alpha=10^{-8}$ of the other phase.
- Left $p=10^8$ Pa, right $p=10^5$ Pa; $u=0$; $\rho_{\rm vap}=2$ kg/m³, $\rho_{\rm liq}=500$ kg/m³.
- NASG, Table 4. 2000 cells, CFL 0.5, $t=473$ μs.
- Two relaxation levels: *p-relax* and *pTg-relax*. In *pTg*, thermo-chemical relaxation activated only at interfaces, $\min(\alpha_{\rm vap},\alpha_{\rm liq})>\epsilon$, $\epsilon=10^{-4}$, and chemical relaxation activated only under $T>T_{\rm sat}$.
- Results: both cases have a leftward rarefaction, rightward contact discontinuity, and shock. With thermo-chemical effects an additional **evaporation wave** appears between rarefaction and contact; the evaporation front produces a liquid–vapour saturation region at higher speed. Left rarefaction runs in almost pure liquid ($\alpha_{\rm vap}\approx0$), where only mechanical relaxation is active, so phase temperatures differ; the vapour temperature there has no physical meaning and is unphysically high (order $10^6$ K), so Fig. 6 shows temperatures only up to 1200 K. *pTg*: liquid–vapour transition leads to **complete evaporation** (handled by the technique of Section 7.2.3, (83)); on the right end a region of nearly pure vapour ($\alpha_{\rm liq}\approx0$), again only mechanical relaxation. Small liquid–vapour mixture zone where thermo-chemical relaxation is active and phase temperatures coincide, ≈ $x=0.83$ m. Exact 5-eq. solution also shown; good agreement for mechanical relaxation only.

#### 8.2.3 Simões-Moreira & Shepherd dodecane evaporation tests
- Validation following Saurel et al. (2008): compare computed evaporation front velocities with experiments of Simões-Moreira & Shepherd (1999): liquid dodecane suddenly discharged into a low-pressure chamber (diaphragm broken), depressurising liquid and initiating evaporation. Initial reservoir pressure ≈ 1 mbar; eight test temperatures $T_{\rm test}$ from 180 °C to 300 °C.
- As in Saurel 2008, Zein 2010, Rodio–Abgrall 2015: instantaneous mass transfer; liquid and vapour dodecane modelled by a **stiffened gas** EOS with the dodecane parameters of Pelanti & Shyue (2014) (good match with experimental saturation curve in this range). $T_{\rm test}$ set everywhere initially. In the low-pressure chamber $p=1$ mbar as in the experiments, and in the other [continued]
