# diffuseIMEX-11 审稿意见（v11 PDF，2026-10-21）

五个角度审稿（审稿人视角、数学正确性、与代码一致性、内部一致性与表述、第 4 节算例设置），每条意见经独立代理尝试反驳核实，已删去被驳倒的条目。位置中的 `l.` 是 PDF 文本提取的行号，同时给出 `.tex` 行号。修改建议保留英文 LaTeX 原文。

## 总体评价

Major revision. The method is new and substantial in two ways. First, the pressure work and the mechanical, thermal and chemical relaxation are solved in one Newton pressure solve with the exact EOS closures, and a relaxation-limit theorem is proved. Second, the shock-energy defect is split between the phases using the shock relations of the model, and the weak-shock split is proved unique. In its present form the paper would not pass a JCP referee, for four reasons. (1) Four stated results are wrong as written: Proposition 3; the linearisation (34)/(36) on which Proposition 5 rests; G_k for NASG phases in Theorem 1, where the code has the same error; and the reason given for using the frozen sound speed in (41). (2) The main relaxation claim (Theorem 4, the finite-mu closure, 'arbitrary rates') is never checked numerically. (3) Two Table 1 problems are run entirely with backward-Euler start-up steps at acoustic Courant numbers up to about 38, which contradicts 'LSDIRK2, C = 0.2'. (4) The paper gives no cost data at all. Each of these can be fixed by restating the result, adding one 0D rate sweep that needs no new code, and adding a cost table that the existing run/pc logs almost fill. Known status, not counted here: Section 4 has no results yet, the front matter is placeholder, and citations print as '?'. Section 3 runs from p.7 to p.19; once results for 22 problems are added the paper will be long, so Section 3 should describe only what the code does.

## 最重要的几点

1. Proposition 3 is false as stated, and the abstract, introduction and conclusions repeat it. Part (i) holds only for SG/NASG phases or Mie–Grüneisen phases of uniform density. With thermal relaxation it also needs a uniform temperature, and with mass transfer it needs T <= T_sat(p). Part (ii) preserves rho_k and T_k only at uniform velocity; for any other velocity only e_k is preserved.
2. Eqs. (34)/(36) and Proposition 5 linearise about P*, the pressure of the explicit state, which already contains the explicit compression. The explicit compression is therefore counted twice, and the hypothesis P - P* = O(M^2) fails. Linearising about the transported state X^tr fixes both.
3. The main relaxation result has no numerical check. Theorem 4 (O(1/(tau*rate)) convergence, exact reduction at zero rate) is never tested, and the finite-mu closure (56) and its mu = 0 limit are never run. A rate sweep on the existing uniform cases 40 and 41, plus one finite-mu shock tube, would cover this.
4. The time stepping is not what Sec. 4 states, and the efficiency evidence is missing. Cases 81 and 85 consist of 25 and 27 backward-Euler start-up steps at acoustic Courant numbers up to about 23 and 38. There is no table of steps, Newton iterations per stage, acoustic Courant number or wall time. The Newton residual test also stalls the pressure convergence at p0 = 1e10.
5. Two formulas or justifications in the dissipation and shock parts are wrong. G_k = (gamma_k+1)/2 does not hold for NASG phases; the correct value is (gamma_k+1)/(2(1-rho_k b_k)), and the code (K_G) has the same error. The statement that the frozen sound speed makes the dissipation never smaller than that of the relaxed waves is false. The weight (41) is also derived for backward Euler, whereas LSDIRK2 has no O(dt) acoustic viscosity.

## 全部条目（按重要性排序）

### 1. [major] Sec. 3.3, Proposition 3 and proof (PDF p.16; paper11.txt l.1027-1037; 03_method.tex l.556-568); Abstract l.43-44; Introduction l.154-155 (01_introduction.tex l.123-124); Conclusions l.1499-1500 (05_conclusions.tex l.33)

Part (i) has three gaps.
- The proof uses (11), which is affine in (alpha_k, m_k) at fixed p only for SG/NASG phases. A Mie–Grüneisen phase whose rho_k varies leaves p after transport: for nitromethane, mixing 0.7 of rho = 1134 with 0.3 of rho = 1300 at 1e9 Pa gives an error of -2.6%.
- With T1 = T2 but T varying between cells, transport gives T1* != T2*. NASG water, cells at 300 K and 400 K: T1* = 313.5 K, T2* = 334.0 K, and the pT pressure is off by -6.1%.
- With T > T_sat(p), mass transfer is activated and the state changes.
Part (ii): Eq. (28) contains the strain term alpha^e*delta[u-bar^e], which the m_k and eps_k updates lack. So rho_k and T_k change under non-uniform u; only e_k is preserved.

**Fix.** Replace the statement by:
\begin{proposition}[Contacts and phase energies]\label{prop:wb}
(i) Let each phase be a stiffened or Noble--Abel stiffened gas, or a Mie--Gr\"uneisen phase of uniform density. A state with uniform $u$ and $p$ and arbitrary $\alpha_k$ is preserved exactly by the scheme with mechanical relaxation; with thermal or chemical relaxation active, this holds if in addition $T_1=T_2=T$ is uniform and, with mass transfer, $T\le T_{\rm sat}(p)$.
(ii) If each phase has uniform $p$ and $T_k$, transport preserves $e_k$ for any velocity, and $\rho_k$ and $T_k$ if in addition $u$ is uniform.
\end{proposition}
Changes to the proof:
- After (11), add 'for a Mie--Gr\"uneisen phase of uniform density, $m_k$ and $\varepsilon_k$ are proportional to $\alpha_k$'.
- Replace the parenthesis about the projection by 'With uniform $T$, transport keeps $T_k=T$ by (ii), so the $pT$ closure and the projection leave the state unchanged, and mass transfer is not activated for $T\le T_{\rm sat}(p)$.'
- Proof of (ii): '$\varepsilon_k=e_km_k$ on every face by the shared limiter, and $m_k$ and $\varepsilon_k$ are transported by the same Rusanov operator and increment; with uniform velocity so is $\alpha_k$.'
Abstract: 'We prove that the scheme preserves moving contacts between stiffened-gas and Noble--Abel stiffened-gas phases, at uniform temperature when thermal or chemical relaxation is active.'
Introduction l.123-124 and Conclusions l.33: replace 'the phase temperatures under transport' by 'the phase internal energies under transport', and qualify the contacts in the same way as the abstract.

_来源：referee-5, math-04, math-05, math-06, code-6, code-7, consistency-01_

### 2. [major] Sec. 3.2 'Derivation of the pressure wave equation': the sentence after (33) and Eqs. (34), (36) (PDF p.9-10; l.590-618; 03_method.tex l.154-171). Proposition 5 and proof (PDF p.17; l.1084-1105; 03_method.tex l.607-628)

X* is already compressed. Eqs. (25)-(26) advance m_k and eps_k with the full flux of u^n, while (28) only advects alpha. So P* already contains -K^c*dt*d_x u^n, and the incremental flux then adds -K^c*dt*d_x(u^{n+1}-u^n).
- Eq. (34) therefore drops the term +K^c*dt*d_x u^n, and (36) counts the explicit compression twice. For an ideal gas the scheme gives p^{n+1} - P* = -gamma*p*tau*d_x u^{n+1} + p*tau*d_x u^n.
- In Prop. 5, P* - P^n contains -tau*K^c*delta[u-bar^e]/dx, which does not scale with M^2, so the hypothesis P - P* = O(M^2) fails.
- The statement says 'For fixed dt', but the proof takes dt from (31).

**Fix.** l.157-158: 'To first order in $\Delta t$ the closure gives the cell pressure as the pressure $P^{\rm tr}$ of the transported state $X^{\rm tr}$, which to first order is the pressure at $t^n$ advected with the flow, minus the response of the mixture to the compression,'
(34): $\bar p^{\,n+1}=P^{\rm tr}-K^Ws$.
(36): right-hand side $P^{\rm tr}-\Delta t\,K^W\partial_xu^*$.
l.154: $p^{n+1}=p^{\rm tr}-\Delta t\,\rho c^2\partial_xu^{n+1}$.
Prop. 5 statement: 'For fixed $\Delta x$, $\Delta t$ given by \eqref{eq:cfl} and data well prepared at both time levels, $P=P_0+M^2P_2$, formally as $M\to0$: ...'. Drop the assumption.
Proof: 'Linearising the closure about $X^{\rm tr}$, ... $P_i-P^{\rm tr}_i=-K^W_is_i+\dots$, with $P^{\rm tr}$ the mixture pressure of $X^{\rm tr}$, which to first order is $P^n$ advected with the flow, so that $P_i-P^{\rm tr}_i=O(M^2)$ for well-prepared data.'

_来源：referee-7, math-01_

### 3. [major] Theorem 4 (PDF p.16; l.1040-1053; 03_method.tex l.573-587); finite-mu closure (56) (PDF p.14; l.906-925); Sec. 4 preamble l.1307; tests/jobs_all.txt

No test exercises the main relaxation result.
- No rate sweep checks the O(1/(tau*rate)) convergence of Theorem 4(i) or the exact reduction at zero rate of 4(ii).
- The finite-mu closure (56) and its frozen two-pressure limit mu = 0 are never run: the jobs use only relax=0..3, theta=2000, nu=1e-4 and no mu=.
Sec. 4 states 'Mechanical relaxation is instantaneous in all tests'. With 'arbitrary rates' in the title, this is the first thing a referee will ask for.

**Fix.** Add one compact test on the existing uniform non-equilibrium cases 40 and 41 (testRelax0D, testRelax0DChem, the SG water–vapour data of Pelanti–Shyue). These need no new code.
- Run tau*mu, tau*theta, tau*nu in {1e-2, ..., 1e4} and infinity.
- Plot the difference from the equilibrium-closure solution against tau*rate on log–log axes; the expected slope is -1.
- Report that rate 0 reproduces the scheme without that relaxation to round-off.
Add one finite-mu run of the water–air tube of Sec. 4.2.3, or of epoxy–spinel, at several mu, compared with mu = 0 and mu = infinity. Then replace l.1307 by 'Mechanical relaxation is instantaneous except in Sec. 4.x; ...'.

_来源：referee-2, consistency-15, referee-1(b)_

### 4. [major] Sec. 3.5 start-up rule (PDF p.17; l.1130-1131; 03_method.tex l.658-660); Sec. 4 preamble (PDF p.20; l.1302-1304); Sec. 4.2.2 Table 1 (PDF p.21): low-Mach water shock tube (case 81) and liquid–liquid shock tube (case 85)

In these two problems |u| is so small that the material step (31) is never reached (case 85) or reached only near the end (case 81). The run/pc/sch logs show:
- case 85: 'steps 27 (rk1 27)' and 'steps 18 (rk1 18)'
- case 81: 'steps 25 (rk1 23)' and 'steps 17 (rk1 17)'
Both problems are therefore computed with backward-Euler start-up stages. The acoustic Courant number grows by a factor 1.2 per step from 0.4 to about 38 (case 85) and 23 (case 81). This contradicts 'all computations use ... LSDIRK2, the material CFL condition (31) with C = 0.2', and the density error hardly decreases with N (1.995e-5 to 1.903e-5).
Also, the first step is the smaller of the material and acoustic steps, not always an acoustic step (vortex p0 = 1: rk1 0).

**Fix.** Choose one of two options.
Option 1: compute the two problems with the acoustic step. Add 'acdt=1 cfl=0.4' jobs for cases 81 and 85 at N = 200 and 1000; such runs already exist in data/Schropff2026_Riemann_problems/acoustic_cfl. In Sec. 4.2.2 write: 'For the low-Mach and the liquid--liquid shock tubes the step \eqref{eq:cfl} exceeds $t_{\rm end}$; these two problems are computed with the acoustic time step $\Delta t=0.4\,\Delta x/\max(|u|+c_f)$, without the start-up rule.'
Option 2: keep the runs and write: 'With \eqref{eq:cfl} these two problems consist of backward-Euler start-up steps only, whose acoustic Courant number grows from $0.4$ to about $23$ and $38$.'
In either case:
- Sec. 3.5: 'The first step is at most $0.4\big(\max(|u|+c_f)/\Delta x+\max(|v|+c_f)/\Delta y\big)^{-1}$ ...; reaching \eqref{eq:cfl} takes about $\log(0.5/M)/\log1.2$ steps.'
- Sec. 4: 'the first time step is limited to the acoustic Courant number $0.4$'.
- Report the number of start-up steps in the result tables.

_来源：referee-4, tests-2, code-9, consistency-17_

### 5. [major] Theorem 1 (PDF p.11-12; l.753-754; 03_method.tex l.318); Bethe–Weyl sentence (l.812-813; 03_method.tex l.373); src/spb6_physics.h K_G l.295

With K = gamma(p+pi)v/(v-b) along the isentrope (p+pi)(v-b)^gamma = const, dK/dp|_s = (gamma*v+b)/(v-b). Hence G_k = (gamma_k+1)/(2(1-rho_k*b_k)), not (gamma_k+1)/2. For the liquid water of Sec. 4.3 this is 3.23 instead of 1.09.
The code uses the same wrong value in the weak-shock share w_k and in the weak-shock dissipation (applied when dp < 1e-2 min K). The shooting branch uses the exact NASG Hugoniot, so for NASG phases the two branches disagree at the switch.

**Fix.** l.318: '...; for stiffened gases $\mathcal G_k=(\gamma_k+1)/2$ and for Noble--Abel stiffened gases $\mathcal G_k=(\gamma_k+1)/\big(2(1-\rho_kb_k)\big)$.'
l.373: '..., which holds for SG and NASG phases, $\mathcal G_k=(\gamma_k+1)/\big(2(1-\rho_kb_k)\big)>0$, ...'
Code: `if (!is_mg(k)) { G = 0.5 * (PH[k].g + 1.0) * v / (v - PH[k].b); return; }`.
Rerun the tests with NASG liquids, and check that the two branches agree at dp = 1e-2 min K.

_来源：referee-6, math-03, consistency-02_

### 6. [major] Sec. 3.2.1 'Acoustic weight', the sentence after (41) (PDF p.11; l.681-682; 03_method.tex l.241-242)

phi_f uses c_f >= c_eq, so phi_f <= (1 - tau*c_eq/dx)^+. For c_eq*tau < dx <= c_f*tau, phi_f = 0, and resolved relaxed waves keep only the implicit damping (1/2)c_eq^2*tau < (1/2)c_eq*dx.
Example, saturated water of Sec. 4.3 at alpha_v = 1e-2 and 1e5 Pa: c_f = 1569 m/s, c_pTg = 1.07 m/s, c_f*tau/dx = 46, c_pTg*tau/dx = 0.031.
So 'the dissipation is never smaller than that of the relaxed waves' is false. c_eq is also never defined.

**Fix.** Replace 'The frozen sound speed is used since ... never smaller than that of the relaxed waves.' by:
'The weight and the impedances use the frozen sound speed, the largest sound speed of the closures \cite{FlattenLund2011}. At a given $\varphi$ the frozen impedances give the waves of a relaxed closure, whose sound speed $c_{\rm eq}$ is $c_p$, $c_{pT}$ or $c_{pTg}$ of (9), the damping $\tfrac14\varphi\Delta x(c_{\rm eq}^2/c_f+c_f)\ge\tfrac12\varphi c_{\rm eq}\Delta x$, while waves with $c_{\rm eq}\tau<\Delta x\le c_f\tau$ receive only the damping of the implicit step.'

_来源：math-02_

### 7. [moderate] Abstract l.31-32; Sec. 3.2.4 l.1020-1021; Sec. 4 (all tests); Conclusions l.1501-1502

The paper gives no cost data: no steps, rejected or start-up steps, Newton iterations per stage, acoustic Courant number or wall time. Yet Sec. 3.2.4 and the conclusions refer to Sec. 4 for the Newton counts, and Schropff et al. report time-step ratios and CPU times (their Tables E.1-E.2).
For transonic 2D problems the gain is doubtful. Helium bubble: 3392 steps, average acoustic Courant number below 1, 2.11 Newton iterations per stage, 5224 s. UNDEX on the full grid takes about 14 h.

**Fix.** Add one table over all tests with these columns: N, steps, start-up steps, rejected steps, Newton iterations per stage, maximum acoustic Courant number $c_f\tau/\Delta x$, and wall time per cell and stage, with hardware and thread count. All of these except the acoustic Courant number are in the first line of each run/pc log; that one needs one extra print.
State that for transonic 2D problems such as the helium bubble the acoustic Courant number stays O(1), so the gain is at low Mach number.
Where the tests coincide, compare time steps with Schropff's Table E.1; no explicit runs are needed.
In Secs. 4.1 and 4.6 say that the Newton iterations per stage are reported.

_来源：referee-3, consistency-04_

### 8. [moderate] Sec. 3.2.4 stopping test (PDF p.15-16; l.1011-1024; 03_method.tex l.550); Sec. 4.6 l.1439-1440; claim of a Mach-independent Newton count at l.1020-1021

The residual is scaled by N = |P| + K + P_E + 1, about 3.4e10 at p0 = 1e10, so the tolerance 1e-13*N is about 3e-3 Pa, the size of the pressure error.
In Taylor–Green at p0 = 1e10, the Newton count is 0.75 and 0.71 iterations per stage and the P2 L2 error goes from 1.90e-3 (N128) to 2.27e-3 (N256). At p0 = 1e6 it converges, 1.94e-3 to 6.23e-4, at 1.00 iterations per stage.
The order study planned at p0 = 1e10 will therefore show a stalled pressure, and the paper presents the falling iteration count as Mach independence.

**Fix.** Either take at least one Newton step per stage, or scale the residual test by $\hat P_{i,j}$ or by $\rho u_0^2$, and rerun Sec. 4.6.
Otherwise state the floor in Sec. 3.2.4: 'the residual test resolves the pressure to $10^{-13}(p_0+K)$, which at $p_0=10^{10}$ is of the order of the discretisation error of $P_2$'. Then restrict the pressure-order study to $p_0\le10^6$, and replace 'independent of the Mach number' by 'does not grow as the Mach number decreases'.

_来源：referee-11_

### 9. [moderate] Theorem 4 and proof (PDF p.16; l.1040-1053; 03_method.tex l.573-587), against Sec. 3.2.2 l.856-858

The proof calls the residual of (60) 'a smooth function of P and zeta', but l.856-858 says it is only continuous and piecewise smooth. It has kinks at the break points delta = 0 and delta = Q of (53), at the donor switch in (40), at the bounds of Y, and at the switch at p_c. The implicit function theorem needs C^1 near the solution.
'Stiffly accurate and L-stable, so every stage lies on the equilibrium manifold' does not follow: each stage is a closure output, stiff accuracy gives U^{n+1} = U_2, and L-stability is not used. 'cell types' is never defined.

**Fix.** Hypotheses: '... and that at this solution no cell lies at a break point of \eqref{eq:Lambda}, at the switch of the donors in \eqref{eq:increment} or at a bound of $Y$.'
(i): replace 'For given cell types' by 'For closure types frozen as in Sec.~\ref{sec:relaxation}'.
Proof: 'The residual of \eqref{eq:global} is then smooth in $P$ and $\zeta$ near the solution, ...'. Replace the last sentence of (i) by 'Every stage value is the output of a closure and hence, for $\zeta=0$, an equilibrium state; since the tableau is stiffly accurate, $U^{n+1}=U_2$ is one too.'

_来源：referee-10, math-07, math-08, consistency-28_

### 10. [moderate] Introduction l.127-130 (01_introduction.tex l.92-94); Sec. 3.2.4 l.1005-1009 (03_method.tex l.535-539); Abstract l.29-30; Conclusions l.1478-1479

The Casulli–Zanolli theory needs a monotone diagonal nonlinearity plus a symmetric M-matrix coupling, and this system has neither.
- The coupling is nonlinear: the energy flux p_f*u_f, the kinetic terms u_f*M-hat and (1/2)u_f^2*M-hat in (rho e)^C, and the donor switch in (40).
- The Jacobian is not symmetric.
- The stated sign pattern is not proved, and it fails in paper tests. With the instrumented finite-difference Jacobian, case 5 has an interface-row west entry of -1.7e-4 at convergence, and case 61 has an off-diagonal entry of -3.06e7 against a diagonal of -4.2e6.

**Fix.** Introduction: '...; the implicit step is a nonlinear system for the cell pressures with a three-point Jacobian, five-point in 2D, whose nonlinearity from the equations of state and the relaxation is local to each cell, solved by Newton's method.'
Sec. 3.2.4: replace 'and the nonlinearity, ..., as in the mildly nonlinear systems of \cite{...}. The Jacobian is ... entries.' by 'and the nonlinearity of the equations of state and the relaxation is local to the cell. The Jacobian has the three-point pattern, five-point in 2D; its part from the strain and the pressure gradient has a negative diagonal and non-negative off-diagonal entries, as in the linearisation \eqref{eq:wave}.'
Abstract and Conclusions: replace 'mildly nonlinear system' by 'nonlinear system'.

_来源：referee-12, math-09, code-5_

### 11. [moderate] Sec. 3.2.1 'Acoustic weight', the derivation before (41) (PDF p.11; l.673-677; 03_method.tex l.231-234)

The term (1/2)c^2*tau is the O(dt) viscosity of one backward-Euler step. The LSDIRK2 tableau (62) has R(z) = (1+(1-2beta)z)/(1-beta*z)^2, and |R(iy)|^2 = 1 + O(y^4) for beta = 1-1/sqrt(2). LSDIRK2 therefore has no O(dt) acoustic viscosity, and the scheme as run carries only phi*(1/2)c*dx = (1-C_a)(1/2)c*dx, below the upwind value that (41) is meant to guarantee.

**Fix.** Write 'one backward-Euler step with \eqref{eq:faces} has ...' and add after 'gives \eqref{eq:weight}': 'The LSDIRK2 stages have no viscosity of order $\Delta t$, so that in the second-order scheme the Laplace terms alone carry the viscosity $\varphi\,\tfrac12c\Delta x$.'
Alternatively, derive $\varphi$ from the LSDIRK2 amplification factor.

_来源：math-12_

### 12. [moderate] Title (PDF p.1); Abstract l.21-28; Conclusions l.1476-1481; Sec. 2.2 l.316; spb6_physics.h closure(): 'if (mg_any() && in.type != CL_P) return false;'

The abstract and conclusions list SG, NASG and Mie–Grüneisen next to 'three relaxation processes, each at a finite or infinite rate'. But Mie–Grüneisen phases are used only with instantaneous mechanical relaxation, and the finite-rate closures (56), (58) and (59) are written for SG/NASG only. 'Unified' in the title appears nowhere in the text.
'At arbitrary rates' itself is acceptable: it is Pelanti's usage under the same ordering, and the abstract already qualifies it.

**Fix.** Abstract, first sentence, and Conclusions l.1476-1478: '... for stiffened-gas and Noble--Abel stiffened-gas phases with mechanical, thermal and chemical relaxation and for Mie--Gr\"uneisen phases with instantaneous mechanical relaxation'.
Title: either use 'unified' once in the text, or write '... with mechanical, thermal and chemical relaxation at arbitrary rates in one local closure'.

_来源：referee-1, consistency-15, consistency-23_

### 13. [moderate] Introduction l.116-121 (01_introduction.tex l.82-83); Sec. 4.2.3 l.1392; Sec. 4.8 l.1456; Sec. 4.9 l.1464; Sec. 4 preamble 'a local Rusanov speed on each face'

The paper criticises Schropff et al. for a per-test dissipation factor, yet switches to the global Rusanov speed per test itself (epoxy–spinel, Richtmyer–Meshkov, bubble compression). The code comments give tuning reasons, not a criterion: 18842 rejections, negative trace vapour mass, L1(u) 118 vs 124.
'Limited to ... instantaneous pressure relaxation' is inaccurate. Their scheme allows mu in [0, infinity]; only their tests use instantaneous relaxation, applied after the IMEX step, as all tests here do.

**Fix.** Replace the clause by 'it is one-dimensional, and its pressure relaxation is applied after the IMEX step and is not part of the implicit system'.
Keep 'and its explicit Rusanov speed is the material speed times a factor set for each test' only if one Rusanov setting is used for all tests here. Otherwise state in Sec. 3.1 a criterion for using the global speed, and drop that criticism.

_来源：referee-9_

### 14. [moderate] Sec. 4.8 l.1451-1456 and Sec. 4.10 l.1471-1472 (PDF p.22-23); paper/sections/04_results.tex l.202, l.222; tests/run_all.sh phase B l.15-16; data/README.md

The PDF grids, 300x100 and 361x200, match the code defaults and the cited papers. Everything else disagrees:
- run_all.sh runs 900x400 and 721x400, and the repository 04_results.tex states those.
- data/README.md contradicts its own notes.
- No SPB6_PC_IMEX run exists for RM, UNDEX or bubble compression. rm.log is empty, and data/ holds V21 SPB6_TV_IMEX runs, some with MOOD steps.
Each 2D problem is computed on one grid only, and UNDEX at half resolution. RM uses SPB's gas gamma = 1.8 on the Pelanti–Shyue grid, whose computation used gamma = 1.4.

**Fix.** Keep the PDF grids.
- In run_all.sh write '.|rm|case=71 snap=2e-3,3.1e-3,4.8e-3,6.4e-3,8.6e-3' and '.|undex|case=73 snap=... probe=0,-1e-6'.
- Align 04_results.tex and data/README.md with the PDF.
- Rerun phase B with SPB6_PC_IMEX.cpp before making the figures.
- Add one two-grid comparison, for example the UNDEX wall-pressure history on 361x200 and 721x400.
- In Sec. 4.8 write '300 $\times$ 100 cells, the grid of Pelanti and Shyue, whose gas has $\gamma=1.4$'.

_来源：tests-1, referee-16, code-8_

### 15. [moderate] Sec. 3.2.2 shock share (PDF p.12-13; l.803-812); Sec. 4.2.3 epoxy–spinel l.1390-1393; Conclusions l.1492 ('No calibrated heat exchange is needed')

The claim 'no calibrated heat exchange is needed' rests on the share heuristics: the pressure range taken from the stencil, the whole or part of the jump depending on the profile width, and a budget frozen per stage. Yet the only mixture-shock test is epoxy–spinel on 500 cells (spb52, case 13); the other shock tests are nearly pure phases with a trace partner.

**Fix.** Run case 13 on 500, 2000 and 8000 cells and add these jobs. Tabulate the errors of the post-shock $\alpha_k$, $\rho_k$ and $T_k$ against the exact five-equation solution. In Sec. 4.2.3 write: 'The post-shock phase states are compared with the exact solution on 500, 2000 and 8000 cells.'

_来源：referee-15_

### 16. [moderate] Conclusions l.1501 (05_conclusions.tex); Sec. 4.1 l.1330-1336

The archived vortex runs at p0 = 1e4 with cell-average data give pressure L2 errors of 9.86e-2, 4.23e-2 and 1.73e-2, i.e. orders 1.22 and 1.29. Second order (2.09) is reached only with well-prepared data, and only unlimited slopes are run.

**Fix.** Conclusions: 'second-order convergence of velocity and volume fraction at all Mach numbers on a smooth two-phase vortex, and of the pressure for well-prepared data'.
Sec. 4.1: say that both cell-average and well-prepared data are reported. Add a minmod row, or state that minmod reduces the order at smooth extrema.

_来源：referee-14_

### 17. [moderate] Sec. 3.2.1 'Second-order factors' (PDF p.11; l.694-697; 03_method.tex l.257-259); Introduction l.150-151 (01_introduction.tex l.118-119); Conclusions l.1496-1498 (05_conclusions.tex l.31)

Only psi^P is reduced by phi_s. psi^u = psi^mm[u] (Eq. 42; SPB6_PC_IMEX.cpp l.941-944) depends on the normal velocity alone. So at an extremum of rho*u or E with monotone u, the velocity term stays second order, although it modifies the momentum and energy fluxes. The sentence 'first order at a local extremum of any quantity whose flux they modify' contradicts (42).
The conclusions' 'limited to second order in smooth regions' reads as a cap on the order.

**Fix.** 03_method.tex l.257-259: 'The pressure term is first order where $P$, a volume-transported quantity or a velocity component has a local extremum, and the velocity term where the normal velocity has one; both are second order where the solution is smooth, which is the LED criterion of Jameson \cite{Jameson1995}.'
Introduction and Conclusions: '..., so that it vanishes at low Mach numbers; it is first order at local extrema and second order in smooth regions'.

_来源：code-1, consistency-07_

### 18. [moderate] Sec. 3.2.3 l.864-866 (03_method.tex l.425-426); finite-theta paragraph l.945-962 (03_method.tex l.484); Mie–Grüneisen p closure l.895-904; Theorem 4(ii) (PDF p.16)

The paper claims exact reduction at zero rate, but two cases break it.
- closure_theta rejects cells with N_k <= 0 for any finite tau*theta, including theta = 0 (spb6_physics.h l.667). closure_p has no such test, and N_k < 0 occurs for a trace gas at shocks, so theta = 0 is not exactly the p closure there.
- The MG p closure rejects a NASG partner with N_k <= 0 (l.526), which the paper does not say.
- For p >= p_c, the finite-nu closure at zeta = 0 picks the end value by the sign of g2 - g1, while the pTg closure picks the phase by comparing v with v_crit. So 'zeta = 0 gives exactly the corresponding equilibrium closure' fails above p_c.

**Fix.** Either change the code so the statements hold: skip the N_k test in closure_theta when tau*theta = 0, and use the pTg end-value rule above p_c. Or qualify the text:
- l.425-426: '... and, for $p<p_c$ and $N_k>0$, $\zeta=0$ gives exactly the corresponding equilibrium closure'.
- Finite theta: 'For finite $\theta$ the closure is used only where $N_k>0$ for both phases.'
- MG p closure: 'a stiffened or Noble--Abel stiffened partner phase requires $N_k>0$'.
- Theorem 4(ii): '$\theta=0$ the $p$ closure where $N_k>0$'.

_来源：code-4, referee-13_

### 19. [moderate] Sec. 3.2.3 p closure (PDF p.14; l.888-895; 03_method.tex l.444-446); finite-mu closure (l.919-921; 03_method.tex l.468-469)

p closure: Phi -> +infinity as p decreases to max_k(-pi_k) requires N_{k*} > 0 for the phase that attains this maximum. 'Since m_k b_k < B_k' fails under compression, for example for a trace NASG liquid with lambda_k|s| > alpha^tr_k - m_k b_k; what is actually needed is sum_k(m_k b_k - B_k)/gamma_k < 0.
mu closure: dp_1/dR < 0 and dp_2/dR > 0 require (gamma_k - 1)p_I + p_k + gamma_k*pi_k > 0, which fails when p_I is sufficiently negative. Yet 'F is increasing ... and the root is unique' is stated without conditions.

**Fix.** p closure: 'If $N_k>0$ for the phase with the largest $-\pi_k$ and $\sum_k(m_kb_k-B_k)/\gamma_k<0$, which holds for $s\ge0$, $\Phi$ changes sign and a root exists; it is unique when $N_k>0$ for both phases.'
MG: 'unique where $c_k^2>0$ on the bracket'.
mu closure: 'If $(\gamma_k-1)p_I+p_k+\gamma_k\pi_k>0$ for both phases, increasing $R$ ... and the root is unique; in general the root is bracketed between $0$ and $R_1$ and found by a safeguarded iteration.'

_来源：math-10, math-11_

### 20. [moderate] Abstract l.27-28; Introduction l.126 (01_introduction.tex l.92); Sec. 3.2.3 opening l.861-862; Conclusions l.1477-1478; Fig. 1 box 'one equation Phi(p) = 0 per cell'

For finite mu the closure solves F(R) = 0 in the volume increment, with a 2x2 linear solve (closure_mu). It is not a scalar equation in the cell pressure.

**Fix.** Sec. 3.2.3: 'which for each set of rates reduces to scalar root finding: $\Phi(p)=0$ in the common pressure, with an inner scalar equation for Mie--Gr\"uneisen phases and finite $\nu$, or $\mathcal F(R)=0$ in the volume increment for finite $\mu$'.
Abstract, Introduction and Conclusions: 'enter one local closure per cell, solved by scalar root finding'.
Fig. 1: 'one scalar equation per cell'.

_来源：consistency-08, consistency-29_

### 21. [moderate] Sec. 3.2.3 'Activation and projection' (PDF p.15; l.991-994; 03_method.tex l.520-522)

The second activation criterion belongs to the code's p–pTg level (RELAX_LEVEL == 3, spb6_physics.h l.809-812), which no Sec. 4 test uses. Sec. 4 defines only p, p–pT and p–pT–pTg, and p–pT–pTg uses the pT temperature. So '(the levels of Sec. 4)' points to a level that is not there.

**Fix.** Replace 'Mass transfer is activated where ... (the levels of Sec.~\ref{sec:results}).' by 'Mass transfer is activated where the temperature of the $pT$-equilibrium state exceeds $T_{\rm sat}(p)$.'

_来源：code-3, consistency-14_

### 22. [moderate] Sec. 2.3, last paragraph (PDF p.7; l.429-432; 02_governing.tex l.246), against (32)

Sec. 2.3 applies the increment X(u^{n+1}-u^n) also to rho*u and E. But (32), Sec. 3.6 and the code (face_flux, cell_residual) carry the momentum as u^{n+1}*M-hat and the kinetic energy as (1/2)(u^{n+1})^2*M-hat, both from the incremental mass flux.

**Fix.** 02_governing.tex l.246: '... $Xu^{n+1}=Xu^n+X(u^{n+1}-u^n)$ for $X\in\{m_1,m_2,\varepsilon_1,\varepsilon_2\}$, the momentum and the kinetic energy being carried with the new velocity by the incremental mass flux \eqref{eq:semidiscrete}: the first part is ...'

_来源：consistency-06, code-11_

### 23. [moderate] Eq. (47) (PDF p.11; l.737-743; 03_method.tex eq:heat)

Eq. (46) defines the dissipation per unit mass of phase k, which gives D_k = integral of (p - p-bar) dv_k. Eq. (47) puts Y_k on the left side only. The proof line 'Y_k D_k = v_k0 Y_k ...' and (48) use the corrected form.

**Fix.** Write (47) as $\mathcal D_k=\int_{\rm profile}(p-p_H)\,\mathrm dv_k$, with $p_H$ the mean shock pressure of rank 26.

_来源：math-13_

### 24. [moderate] Sec. 4.9 (PDF p.22; l.1458-1464)

'the left wall moves with velocity 100 m/s in +x (the grid is fixed, ...)' describes a moving wall that the code does not have. bc_apply mirrors 2u_w - u at a fixed x = 0, and the face carries the mass flux rho*u_w, so liquid enters the domain: dM_1/M = 0.129 at 1.2 ms. Pelanti–Shyue's piston would have advanced 0.12 m by then.

**Fix.** Replace by: 'The top, bottom and right boundaries are walls. The left boundary stays at $x=0$; its ghost cells mirror the adjacent cell with the normal velocity $2u_w-u$, $u_w=100$~m/s, so that liquid enters the fixed grid at 100~m/s in place of the moving piston of Pelanti and Shyue.'
In the first sentence, change 'compressed by a piston' to 'compressed by a liquid stream entering at 100~m/s'.

_来源：tests-4_

### 25. [moderate] Sec. 4.5 (PDF p.22; l.1424-1425); Sec. 4 preamble l.1311-1313

The 'exact' HEM reference (tests/exact_hem_bartak.py) neglects the initial vapour (alpha_v = 1e-3, Y0 = 6.3e-5), which the computation carries. It is exact at the probe only because the reflected wave arrives after t_end. The paper states neither point.

**Fix.** After 'which is computed exactly' add: 'The initial vapour is neglected in this reference; the rarefaction consists of the isentropic liquid wave down to the saturation pressure followed by the equilibrium two-phase fan. The wave reflected at the closed end flashes at once, moves at about 32~m/s and reaches the probe only after about 50~ms, so that the simple wave is exact at the probe up to $t_{\rm end}$.'

_来源：tests-5_

### 26. [moderate] Eqs. (45)-(50), Theorem 1 and proof, l.804 (PDF p.11-12; l.726-814); (46) 'P = p - Pi'; (34) and Prop. 5 'P*'; Secs. 4.1 and 4.6 'p0'

The pressure symbols are overloaded.
- p-bar is the mixture pressure alpha_1 p_1 + alpha_2 p_2 in Sec. 2, (32)-(33) and (60), but (p_0 + p_1)/2 in (45) and (47).
- p_1 is the phase-1 pressure and also the post-shock pressure.
- P is the Newton unknown in (38) and (60), and the total stress and Rayleigh line in (46).
- P* is never said to be p of (16).

**Fix.** - In (45)-(50), Theorem 1, its proof and l.804 ($p^-=p_{\rm lo}<p^+=p_{\rm hi}$), use $p^-$ and $p^+$ for the pre- and post-shock pressures and $p_H=\tfrac12(p^-+p^+)$ for their mean.
- Write the total stress as $\sigma=p-\Pi$ and the Rayleigh line as $\sigma(v)$ in (46).
- State once in Sec. 3.2 that $P_{i,j}$ is the Newton unknown, that $\bar p^{\rm cl}=P$ at convergence, and that $P^{\rm tr}$ is $p$ of (16) evaluated on $X^{\rm tr}$.

_来源：math-14, consistency-09_

### 27. [moderate] (18) K^c; (34) K^W; Sec. 3.2.1 l.653-654; Theorem 1 K_k, repeated at l.1060; Sec. 3.2.4 K_{i,j}, delta P_{i,j}, D_{i,j}; (63) K_i; (42) denominator '2D'; D_0 l.854; Prop. 5 'delta.', 'delta P_2', 'tau delta overline{P_2}/dx'

K has six meanings: K^c, K^W, K_k, a bare K in Sec. 3.2.1, K_{i,j}, and the stage slopes K_i. In Prop. 5, K_i^W carries a cell index i next to the stage index i.
delta stands for the energy defect, the face difference, the Newton step, and an undefined divergence 'delta.'.
In (42) the denominator '2D' reads as 'two-dimensional'.
In Prop. 5 the velocity correction lacks the density and the sign.

**Fix.** - Rename the stage slopes $H_i$ in (63) and Sec. 3.5.
- Write $K^W$ in Sec. 3.2.1 ('$\tau K^W\mathcal L(\eta^P,P)$, $K^W\eta^P=\tfrac12c\Delta x$') and $K^W_{i,j}$ in Sec. 3.2.4; define $K_k=\rho_kc_k^2$ once.
- Write the Newton step as $\Delta P_{i,j}$; replace $D$ in (42) by the face jump $[X]_f$; write $D_0\to\delta_0$ and $D_{i,j}\to J_{i,j}$.
- In Prop. 5 define $\mathrm{div}_h u_f=\delta_x[u_f]/\Delta x+\delta_y[v_f]/\Delta y$, write $\delta_xP_2$, and write the correction as $-\tau\,\delta_xp_{f,2}/(\rho\Delta x)$ with $p_{f,2}$ the weighted face average of $P_2$.

_来源：consistency-10, consistency-11, consistency-28_

### 28. [moderate] tau in (38) l.641 and (41) l.676, defined only at l.703 and in Sec. 3.5 l.1129; a_ii in Prop. 5; Y_min l.978 (defined l.987); S_k, h_k at the Gibbs relation (02_governing.tex l.53); c_eq l.681; (.)^+- in (40); Y^tr; e in the Y interval; p_sat in (65); SG, LSDIRK2, LED, and IMEX in the abstract

Sec. 3.2 uses dt in (32)-(36), then Sec. 3.2.1 switches to tau without explanation; tau is defined only later, in Sec. 3.2.2 and Sec. 3.5. Several other symbols and abbreviations are used before they are defined, or never defined.

**Fix.** - At the start of Sec. 3.2 write: 'The implicit step $\tau$ equals $\Delta t$ for a first-order step and $a_{ii}\Delta t$ in a stage of Sec.~\ref{sec:time}.' Use $\tau$ in (32)-(36).
- Define $S_k$ and $h_k$ at the Gibbs relation, $(a)^+=\max(a,0)$ and $(a)^-=\min(a,0)$ at (40), $Y^{\rm tr}$ as the mass fraction of $X^{\rm tr}$, $e=(\rho e)^C/\rho$, $p_{\rm sat}$ as the inverse of $T_{\rm sat}$, and $Y_{\min}$ at first use.
- Define $c_{\rm eq}$ as in rank 6.
- Spell out stiffened gas (SG), implicit--explicit (IMEX) in the abstract, and LSDIRK2 and LED at first use.

_来源：consistency-13_

### 29. [minor] Abstract (PDF p.1; l.21-45)

'converges to the corresponding equilibrium schemes as the relaxation rates tend to infinity' drops the fixed-mesh and non-degeneracy qualifiers of Theorem 4, which the introduction and conclusions keep. At scheme level, the equilibrium pT and pTg schemes also project U^n while the finite-rate ones do not.
IMEX is not defined. The abstract has about 336 words, and its second sentence alone about 91.

**Fix.** Write '..., and, at fixed mesh and under non-degeneracy assumptions, converges to the corresponding equilibrium schemes as the relaxation rates tend to infinity'. In the Conclusions, add 'for initial data on the equilibrium manifold'.
Write 'implicit--explicit (IMEX)'.
Split the second sentence after 'in a single implicit step.' and cut the Zha–Bilgen details, to reach about 250 words.

_来源：referee-19, math-17, consistency-27_

### 30. [minor] Sec. 3.2.3 finite chemical relaxation l.981-982; Sec. 4.5 l.1423-1424

Freezing the state-dependent rate (65) at U^n for both stages gives an O(dt) error in the rate. The Barták run is therefore first order in time where the transfer is not stiff. The same sentence also appears twice.

**Fix.** Keep the sentence once, in Sec. 4.5, and append ', which makes the step first order in time where the transfer is not stiff'. Alternatively, evaluate (65) on $U_{E,i}$ in each stage.

_来源：referee-20, math-18, consistency-20_

### 31. [minor] Introduction l.93-95 and l.132-134 (01_introduction.tex l.96-99); refs.bib

The all-regime Lagrange-projection schemes for two-phase models are not cited (Chalons, Girardin, Kokh, J. Comput. Phys. 335, 2017). The Sirianni & Re sentence does not mention that they treat concurrent finite-rate relaxation, which the ordering used here excludes.

**Fix.** Cite Chalons, Girardin and Kokh (2017) with the all-speed and preconditioning works.
After '... emulate it with large finite rates, respectively.' add: 'Sirianni and Re also treat concurrent finite-rate mechanical and thermo-chemical relaxation, which the ordering of Pelanti used here excludes.'

_来源：referee-8, referee-17_

### 32. [minor] Sec. 3.2.3 'Activation and projection' and Sec. 4 preamble; SPB6_PC_IMEX.cpp l.1300 versus l.1292-1293

The final state is projected (RelaxationActivation) before write_profile, but snapshots and probe values are written without projection. For runs with relax > 0, the final profile and the time histories are therefore not the same discrete state, and the paper does not say so.

**Fix.** Add 'The solution at the output times is the projected state.' Then either project before the snapshot and probe writes too, or drop the final projection.

_来源：code-2_

### 33. [minor] Sec. 3.6 step 2 (PDF p.19); SPB6_PC_IMEX.cpp face_x and face_y l.597-616, M_P l.683-692

The treatment of wall and transmissive faces in the pressure system is not described, although Secs. 4.5 and 4.7-4.10 depend on it.

**Fix.** Add to step 2: 'On a wall face $u_f$ is the wall velocity. On wall and transmissive faces the ghost pressure is the pressure of the adjacent cell, and the ghost velocity used in $D^*_f$ is that of the adjacent cell, with the normal component replaced by $2u_w-u$ at walls.'

_来源：code-17_

### 34. [minor] Sec. 3.2.4 l.1019-1023; SPB6_PC_IMEX.cpp GMRES_Solver l.741-757; line search l.828-832

After 4000 iterations, BiCGStab is accepted at a relative residual of 1e-10 (1e3 times the tolerance). This is an unstated fallback, against the author's rule. The Armijo factor 1 - 1e-4*lambda and the limit of 40 halvings are also unstated.

**Fix.** Remove the 1e3 acceptance, so that the linear solve either converges to 1e-13 or the step is rejected. Optionally state 'BiCGStab with at most 4000 iterations' and 'Armijo factor $1-10^{-4}\lambda$, at most 40 halvings'.

_来源：code-19_

### 35. [minor] Sec. 3.1 after (31) (PDF p.9; l.556-558; 03_method.tex l.114-116)

'with C = 0.2 unless stated otherwise' is leftover text: every Sec. 4 log reports CFL 0.2. 'all computations of Sec. 4 are advanced with it' does not fit the fixed-step runs of Sec. 4.6.

**Fix.** 'with $C=0.2$, and the first steps are limited by the start-up rule of Sec.~\ref{sec:time}. ... All computations of Sec.~\ref{sec:results} use time steps within \eqref{eq:cfl}, those of Sec.~4.6 a smaller fixed step.'

_来源：code-22, consistency-04_

### 36. [minor] tests/run_all.sh l.5-7; data/README.md l.3; tests/jobs_all.txt; src/bin_SPB6_PC_IMEX

The documented route to reproduce the results does not run the program the paper describes.
- run_all.sh defaults to SRC=SPB6_TV_IMEX.cpp, which updates alpha in the implicit step instead of by (28), and the README command passes no source.
- jobs_all.txt contains jobs that are not in the paper (cases 0-3, 6, 12, 87, 88, and relax=3) and lacks the acdt=1 runs.
- The binary is older than the last source commits.

**Fix.** Set SRC=${2:-SPB6_PC_IMEX.cpp*} and write 'bash tests/run_all.sh run/pc SPB6_PC_IMEX.cpp' in the README. Mark or move the jobs that are not in the paper. Add the jobs of ranks 3, 4 and 15. Rebuild and rerun before making the final figures.

_来源：code-20, tests-12_

### 37. [minor] Fig. 1 (PDF p.18; l.1138-1209)

Figure 1 is not cited anywhere in the text; the citing sentence was removed in this build.
- The box 'Algorithm: one time step' shows one stage and numbers its steps 0-3, while Sec. 3.6 numbers the same stage 1-6.
- 'Flux Vector Splitting' is a term the text never uses.
- The figure floats onto a page of its own.

**Fix.** Begin Sec. 3.6 with 'Fig.~\ref{fig:scheme} summarises the split and one stage of the scheme.'
Title the box 'Algorithm: one stage' and number its steps as in Sec. 3.6.
Use 'Zha--Bilgen splitting' in the box and the caption.
Place the figure with [t] near Sec. 3.6.

_来源：consistency-03, consistency-29_

### 38. [minor] Sec. 4.1 first sentence and Eq. (64) (PDF p.20; l.1316-1337)

The vortex is the author's single-velocity adaptation of the Battisti–Boscheri set-up: equal phase rho and p, and a passive alpha_2. The text presents it as their test, and it does not name the variables whose errors are measured.

**Fix.** 'This test is adapted from Battisti and Boscheri \cite{BattistiBoscheri2025}, Sec.~4.1: their domain, adiabatic exponents, translation velocity and background pressures, both phases carrying the density and pressure of the $\gamma=1.4$ vortex, as the single-pressure model requires, and $\alpha_2$ being a passive profile.'
Name the measured variables ($\rho$, $u$, $p$; and $\alpha_2\rho_2$, $\alpha_2\rho_2u$, $\alpha_2$ if they are compared with their Table 1).

_来源：tests-3_

### 39. [minor] Sec. 4.6 and Eq. (66) (PDF p.22; l.1428-1441)

Three points are not stated.
- In the two-phase case alpha_2 is transported, since u.grad(alpha_2) != 0.
- (66) is stationary for the compressible equations only up to O(M^2).
- The reported face-velocity divergence differs from Battisti–Boscheri's central difference of the cell velocities (their Table 4).

**Fix.** Add: 'In the two-phase case $\alpha_2$ is transported by the flow and is not compared; the compressible solution approaches \eqref{eq:tg} to $O(M^2)$.' State which divergence is reported, and add the central difference of the cell velocities if the results are compared with their Table 4.

_来源：tests-6_

### 40. [minor] Secs. 4.1, 4.2.1, 4.2.2, 4.4, 4.6, 4.9, 4.10 (PDF p.20-23)

Several small set-up items are missing or inexact.
- No relaxation level is stated for Secs. 4.1, 4.2.1, 4.2.2 and 4.6; the code uses p.
- Sec. 4.4 gives no domain.
- 'SG' is never defined.
- The 5-tuples in Sec. 4.9 have no labels.
- The UNDEX probe is the wall cell centred at (0, -dy/2), not the point (0, 0).

**Fix.** - Sec. 4 preamble: 'The relaxation level is $p$ unless stated otherwise.'
- Sec. 4.4: 'on $[0,1]$'.
- Sec. 2: 'the stiffened gas (SG), the NASG law with $b=0$'.
- Sec. 4.9: 'SG with $(\gamma,\pi,q,q',c_v)=\dots$'.
- Sec. 4.10: 'the pressure is recorded in the wall cell at $x=0$'.

_来源：tests-9_

### 41. [minor] Sec. 4.7 (PDF p.22; l.1443-1449)

The state of the trace helium behind the shock is set in the code (rho = 0.167 at the post-shock p and u) but not stated in the paper.

**Fix.** Add 'Behind the shock the trace helium has $\rho=0.167$ at the post-shock pressure and velocity.'

_来源：tests-7_

### 42. [minor] (2) l.198; (23) l.461; Theorem 1 '|_s'; Sec. 3.6 U^s; l.804-814; (32) u_n; (64) eps = 5 and T; (59) phi_k; l.920 R_1; Q_k vs script Q vs bold Q; Z_k vs Z_{i,j}; (65) p_crit and C_r; Sec. 4.3 theta and nu

Remaining symbol clashes:
- k is the phase index and the kinetic energy.
- v is the y-velocity and the specific volume.
- s is the strain, the path parameter, the entropy subscript and the update base.
- u_n collides with u^n.
- eps = 5 collides with eps_k, and T in (64) is not a temperature.
- phi_k in (59) has the glyph of the limiter phi_i.
- R_1 is the bracket end and a JWL parameter.
- p_crit in (65) is the critical pressure of real water, not p_c.
- C_r, theta and nu have non-obvious units.

**Fix.** - In (2) write $E=\varepsilon_1+\varepsilon_2+\tfrac12\rho|\mathbf u|^2$.
- In Sec. 3.2.2 say that $v$ denotes the specific volume.
- Write $|_S$ for the entropy subscript and $U^b$ for the update base; use $u_\perp$, $u_\parallel$ for the normal and tangential velocity.
- (64): $\rho_1=\rho_2=\big(1-\tfrac{25(\gamma-1)}{8\gamma\pi^2p_0}e^{1-r^2}\big)^{1/(\gamma-1)}$, $p_1=p_2=p_0\rho_1^\gamma$, $\gamma=\gamma_2$.
- Drop $\phi_k$ in (59).
- Rename the bracket end $R_{\max}$ and the shock budget $B^{\rm sh}_k$.
- After (65): '$p_{\rm crit}=22.064$~MPa, the critical pressure of water, not the end point $p_c$ of the NASG saturation curve; $C_r=1$~Pa~kg$^2$/(s~J$^2$)'.
- Sec. 4.3: $\theta=2000$~Pa/(s~K), $\nu=10^{-4}$~Pa~kg$^2$/(s~J$^2$).

_来源：consistency-12, tests-8_

### 43. [minor] Eq. (18), K^c (PDF p.6; l.361-381)

Row 7 needs kappa_k and Gamma_k at (rho_k, p), whereas rows 5-6 use them at the phase state (rho_k, p_k). The two differ whenever p_k != p.

**Fix.** After the definition of $K^c$ add 'with $\kappa_k$ and $\Gamma_k$ evaluated at $\rho_k$ and the pressure $p$'.

_来源：math-16_

### 44. [minor] Sec. 3.2.1, after (41) (PDF p.11; l.679; 03_method.tex l.239)

'for tau -> 0 (38) gives the acoustic Godunov face states' is wrong. Eq. (38) has arithmetic and density-weighted means, while the Godunov states use impedance-weighted means; only the diffusive terms coincide, which the previous paragraph already says.

**Fix.** Delete 'for $\tau\to0$ \eqref{eq:faces} gives the acoustic Godunov face states, and'. Begin the clause with 'For $c\tau\ge\Delta x$, ...'.

_来源：math-15_

### 45. [minor] Sec. 3.5 (PDF p.17; l.1132-1134; 03_method.tex l.660-663)

'non-positive density or volume' does not match the code. The volume test is vacuous, since alpha_1* + alpha_2* = 1. The density test covers U_{E,i} and X*. Closure and linear-solve failures also reject the step.

**Fix.** 'If in some stage the mixture density of $U_{E,i}$ or $X^*$ is non-positive in a cell, a local closure has no admissible solution, or the Newton iteration fails, the step is repeated with half the time step.'

_来源：code-10_

### 46. [minor] Eq. (39) and Eq. (42) (PDF p.10-11; l.655-697); Secs. 4.1 and 4.6

Several reconstruction details are not stated exactly.
- The code uses the velocities u*, v* of X* in (39).
- In 2D, psi^u is the minmod factor of the normal velocity and multiplies D*_f.
- phi-hat uses the factor after the cap (30), and the cap is skipped with unlimited slopes.
- 'previous Newton initial guess' suggests the guess of the previous stage, but the code uses the guess of the current stage.

**Fix.** - Write $u^*$, $v^*$ in (39).
- After (42): 'in two dimensions $\psi^u$ is evaluated on the normal velocity component and multiplies $D^*_f$'.
- '$\hat\phi$ the minimum of the capped shared factor (29)-(30) and ...'.
- Secs. 4.1 and 4.6: 'with unlimited central slopes, $\phi_i=1$ and (30) is not applied'.
- 'The slopes of $P$ are those of the Newton initial guess of the stage, the pressure of the previous stage, and those of $u$ are taken on $u^*$.'

_来源：code-12, code-13, code-21_

### 47. [minor] Sec. 3.2.2 and Theorem 2 weak-shock sentence (PDF p.12-13; l.789-815)

The text reads as if the restriction of Q_k to compressive cells were derived; the code imposes it. Q_k = 0 also where the stencil pressure does not vary or no admissible end state exists. w_k is written with alpha_k, while the code uses rho*Y_k*v_k0 with K_k and G_k evaluated at (p_0, v_k0).

**Fix.** '$Q_k$ is computed only in cells where the central difference of $u^*$ and $v^*$ has negative sum and the stencil pressure varies, a set that contains every admissible shock since ...; elsewhere, and where no admissible end state exists, $Q_k=0$.'
After '$\lambda_k=w_k$' add 'with $\alpha_k=\rho Y_kv_{k0}$ and $K_k$, $\mathcal G_k$ evaluated at $(p_0,v_{k0})$'.

_来源：code-14, code-16_

### 48. [minor] Sec. 3.2.3, MG inner solve (PDF p.14; l.903-904); outer equations of the pT, finite-theta, finite-nu and pTg closures

The MG inner solve takes Newton steps in alpha, with bisection in ln(alpha) only as a safeguard. The paper never says how Phi(p) = 0 is solved for the pT, finite-theta, finite-nu and pTg closures; the code uses Brent's method in ln(p - max_k(-pi_k)) (root_above).

**Fix.** Write 'a Newton iteration in $\alpha$, safeguarded by bisection in $\ln\alpha$', and add once: 'the equations $\Phi(p)=0$ of the $pT$, finite-$\theta$, finite-$\nu$ and $pTg$ closures are solved by Brent's method in $\ln(p-\max_k(-\pi_k))$'.

_来源：code-15_

### 49. [minor] Sec. 3.2.4 P_E (PDF p.15; l.1013-1014)

|m_k q_k| is not defined for MG phases; the code uses |m_k e_r| + |alpha_k p_r/Gamma_k| for them. The default K = max_k pi_k for an inadmissible X* is not stated.

**Fix.** Add 'with $|m_k\varepsilon_r|+|\alpha_kp_r/\Gamma_k|$ in place of $|m_kq_k|$ for Mie--Gr\"uneisen phases, and $K=\max_k\pi_k$ where $X^*$ is not admissible'.

_来源：code-18_

### 50. [minor] Introduction l.145 (01_introduction.tex l.113); Conclusions l.1493

The introduction says 'with the velocity' and the conclusions 'with the old velocity'. In the second stage the explicit fluxes use the velocity of U_{E,2}, not u^n.

**Fix.** Introduction: 'transports all quantities with the flow velocity'. Conclusions: 'transported explicitly with the velocity of the explicit state'.

_来源：consistency-05_

### 51. [minor] Conclusions l.1491-1492, against Proposition 2 (PDF p.13; l.823-844)

'for every defect' overstates Proposition 2: no admissible partition exists for delta <= -sum_j E_j.

**Fix.** '... the discretisation error with the unique linear share that keeps both phases admissible for every defect above $-\sum_j\mathcal E_j$'.

_来源：consistency-16_

### 52. [minor] Introduction road map l.160-162; Sec. 3 overview l.437-440; Sec. 3.2 l.569

- The road map omits the split form of Sec. 2.3.
- The Sec. 3 overview lacks 'and', omits Sec. 3.2.4, and places all proofs in Secs. 3.3-3.4, although Theorem 1 and Prop. 2 are in Sec. 3.2.2.
- Sec. 3.2 says Sec. 3.6 gives the fully discrete counterpart, although Secs. 3.2.1-3.2.4 are already fully discrete.

**Fix.** Road map: 'Section 2 presents the model, the equations of state and the split form'.
Overview: '..., and the implicit update of the pressure--relaxation sub-system with its dissipation, local closure, relaxation closures and pressure solver (Sec. 3.2)'.
Sec. 3.2: 'Sec.~\ref{sec:fully} collects the stage.'

_来源：consistency-19_

### 53. [minor] Sec. 3.6 last paragraph l.1284-1289; Sec. 4 preamble l.1304-1305; (43), (44), (56) omega_k; the zeta definitions in each closure

Repetitions, against the author's conciseness rule.
- The 'In summary ...' paragraph repeats Sec. 3.2.1.
- Sec. 4 repeats eps_I = 1e-4 and the 0.4 acoustic first step.
- omega_k = alpha*_k is a redundant symbol.
- zeta is redefined in every closure.

**Fix.** Delete the Sec. 3.6 summary paragraph and the Sec. 4 repeats. Write $\alpha^*_k$ for $\omega_k$, and drop the repeated $\zeta$ definitions. Keep the functionals $L^{(1)}$, $L^{(2)}$ in Sec. 3.1, since (30) needs them.

_来源：consistency-20_

### 54. [minor] Sec. 3.6 step 3 and last paragraph (PDF p.19; l.1246-1282)

u_t is defined as the tangential velocity, but the momentum fluxes use u_t on y-faces and v_t on x-faces for the same quantity. rho^{n+1} appears next to the superscript 'new'. 'H is the one-stage update given by steps 1-6' is wrong: in (61), H is a time derivative.

**Fix.** Write the donor velocities explicitly ($v^{L,R}$ on x-faces, $u^{L,R}$ on y-faces) and write $\rho^{\rm new}$.
'Steps 1--6 give the stage value $U_i$, and $H$ in (63) is $(U_i-\tilde U_{I,i})/\tau$; step 1 is its explicit part, evaluated on $U_{E,i}$, and steps 2--6 its implicit part.'

_来源：consistency-21, consistency-22_

### 55. [minor] Keywords (PDF p.1); Sec. 3.3 heading l.1026

The keyword 'All Mach number IMEX scheme' does not match the title's 'all-Mach'. The Sec. 3.3 heading 'Well-balanced property' names what Prop. 3 calls contact preservation.

**Fix.** Keyword 'All-Mach scheme'. Optionally rename Sec. 3.3 'Contact preservation'.

_来源：consistency-23_

### 56. [minor] (28) l.508; (41) l.677; l.1377-1378; l.945-946; l.1380; (8) l.259; l.79 vs l.719; headings l.1076, l.1084, l.1427; citations used as nouns l.105, l.614-615

Typographic slips:
- (28) and (41) end with a spaced colon, which renders as a ratio sign.
- '1000 and 100 cells' reverses the usual order.
- 'eps1 += tau Q, eps2 -= tau Q' is C syntax.
- Citation parentheses are nested or back-to-back.
- CFL in (8) but C in (31).
- 'mixture-energy-consistent' and 'mixture-energy consistent' both occur.
- Headings mix 'Low-Mach' and 'Low Mach number'.
- Citations are used as nouns.

**Fix.** - End (28) and (41) with a comma and continue '..., which ...'.
- 'with 100 and 1000 cells'.
- '$\varepsilon_1$ gains and $\varepsilon_2$ loses $\tau\mathcal Q$'.
- 'by the Cochran--Chan law (13) as in \cite[Sec.~4.2]{...}'.
- Use $C$ in (8); hyphenate one way; use one heading form; name the authors when a citation is a noun.

_来源：consistency-24_

### 57. [minor] Introduction l.100-102 and throughout Secs. 3-4, e.g. l.1347-1348, 1372, 1386-1387, 1418-1419, 1434-1436, 1455-1456, 1463-1464, 1471-1472

Many parentheses remain in running text, against the author's rule. Examples: the Dumbser–Casulli sentence '(the equation of state confines the nonlinearity to the diagonal)', '(five-point in 2D)' three times, '(frozen sound speed on X*)', '(M = |u|/c between ...)', '(single-phase case)', and a nested pair in Sec. 4.8.

**Fix.** - Dumbser–Casulli: 'this system is only mildly nonlinear, the nonlinearity being confined to the diagonal through the equation of state'.
- 'with a three-point Jacobian, five-point in 2D'.
- '; $M=|u|/c$ is $2.7\times10^{-4}$ in the air and $6\times10^{-4}$ in the water'.
- 'on $[0,1]$, and on $[-10,10]$ for the detonation products,'.
- '$\alpha_2=1/2$ in the single-phase case and $\alpha_2=\dots$ in the two-phase case'.
- '$p$--$pT$ without mass transfer, $p$--$pT$--$pTg$ with instantaneous mass transfer'.
- '361 $\times$ 200 cells, half the resolution of Pelanti \cite{Pelanti2022}'.

_来源：consistency-25, tests-10_

### 58. [minor] Conclusions l.1505; Sec. 4.2 heading; Sec. 4.2.1 heading

The conclusions list a 'bubble-collapse' problem, but Sec. 4.9 is the compression of a vapour bubble. The heading 'One-dimensional Riemann problems' for Sec. 4.2 also fits Secs. 4.3-4.4. Sec. 4.2.1 'Interface advection' contains a single-gas density step.

**Fix.** Conclusions: 'vapour-bubble compression'.
Sec. 4.2: 'One-dimensional problems with mechanical relaxation'.
Sec. 4.2.1: 'Advection of interfaces and contacts'.

_来源：consistency-26, tests-11_
