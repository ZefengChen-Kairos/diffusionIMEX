# flattener 的替代方案调查:36 个候选,没有一个更好

背景:L9 在 HLLEM 之上只加了一样东西 —— 把线性退化场的反扩散按
$\delta_k\leftarrow(1-s_p)\delta_k$ 衰减($s_p$ 是面邻域的相对压力跳,0.5 为饱和尺度)。
为了能回答 referee 的"为什么不用方法 X",对 36 个公开的替代构造做了系统评判,每个按
五个约束打分:AP 标度不被破坏、三阶不被破坏、与 $\eps$-split 的面分解相容、
MHD 可用、无可调常数。

**判决分布**

| 判决 | 个数 |
|---|---|
| equal_tradeoff | 1 |
| not_applicable | 9 |
| worse | 26 |

没有 `better`。唯一的 `equal_tradeoff` 是 ADC(见下),而它的致命问题是它限制的那个量
在本格式里含义不同。

## 逐条

| 判决 | 候选 | 为什么不用 |
|---|---|---|
| `equal_tradeoff` | Anti-Diffusion Control (ADC) for HLLEM — Simon & Mandal, IJNMF 2019 (arXiv:1806.03187) | ADC limits a quantity that does not mean here what it means in the paper. In classical HLLEM, delta_k = 1 - lam_-/sL - lam_+/sR measures an LD wave's speed against the ACOUSTIC envelope (sL ~ u-c_f, sR ~ u+c_f), and ADC's linear p… |
| `not_applicable` | Dumbser & Balsara's own HLLEM/HLLI formulation (JCP 304 (2016) 275-319) | It is not an alternative to the current cure -- it IS the current cure, and the scheme already cites it as such. /home/user/diffusionIMEX/blastwave_debug/paper/cas-refs-capitalised.bib:979 shows \cite{HLLEM} is this exact paper, p… |
| `not_applicable` | Entropy-stable matrix dissipation in entropy variables (Chandrashekar, CiCP 14(5) 2013, 1252-12… | The eps-split explicit sub-system provably admits NO strictly convex entropy in the stiff limit, by the author's own theorem. Theorem thm:hyp(iii) (paper.tex:863-875, proved at :900-972) states that at eps_M=0 the split primitive … |
| `not_applicable` | HLLC+ — antidissipation pressure fix (Chen, Lin, Li, Yan, SISC 42 (2020) B921–B950) | It acts on the wrong term. HLLC+'s antidissipation pressure fix modifies the pressure/acoustic part of the flux, but in this scheme the acoustic analogue — the magnetosonic family — carries NO anti-diffusion at all (Final_emb_omp_… |
| `not_applicable` | HLLC-type all-Mach / shock-robust solver (Xie group), IJNMF 89(10) 430-463, 2019 — "An accurate… | The HLLC star state does not exist for this scheme's explicit flux. HLLC closes the HLL fan by imposing, across each outer wave K with m_K = rho_K (S_K - u_K), the Rankine-Hugoniot conditions of the flux actually being solved toge… |
| `not_applicable` | Harten-Hyman entropy fix (HH1/HH2) on the Roe linearization | A two-sided pincer, both sides provable from the scheme's own flux. (1) On the fields HH is DEFINED for (genuinely nonlinear, i.e. magnetosonic) this scheme carries no anti-diffusion at all: delta_k is applied only to entropy + Al… |
| `not_applicable` | LHLLD — low-dissipation HLLD (theta shock-detecting factor on S_M + phi low-Mach total-pressure… | Both halves act on HLLD intermediate states that this scheme does not contain. The explicit solver is AP-HLLEM: a split-Einfeldt HLL base plus Roe-eigenvector anti-diffusion delta_C/delta_A+/delta_A- on the three linearly degenera… |
| `not_applicable` | MLAU — Multistate Low-dissipation AUSM for ideal MHD (Minoshima, Kitamura & Miyoshi, ApJS 248:1… | MLAU's low-speed machinery is dead, by its authors' own stated validity limit, in exactly the region where this scheme's carbuncle lives. The paper's own blast-wave diagnostics (paper.tex, Sec. "Ideal MHD blast wave problem", line… |
| `not_applicable` | Mass-flux / transverse-dissipation consistency condition for grid-aligned shocks (Xie et al., J… | Its only operative claim is factually inapplicable to this scheme. The leverage the task attributes to it is the negative result "entropy-wave dissipation alone cannot stabilise a strong shock", on the premise that "damping delta_… |
| `not_applicable` | SLAU — parameter-free simple low-dissipation AUSM-family mass flux (Shima & Kitamura, AIAA J 49… | SLAU cures a disease this scheme does not have, by a mechanism the scheme has already applied analytically. SLAU's shock stability IS Liou's criterion, and Liou's criterion is a statement about the mass flux. In this code the phys… |
| `worse` | ASHLLEM / ASHLLC (Xie et al., unified all-speed HLL-type construction, Computers & Fluids 2021) | The construction's premise does not exist in this scheme. The mass-flux-consistency condition is read off a linearisation around a planar shock that is a DISCRETELY STEADY state of ONE flux function — the physical Euler flux. Here… |
| `worse` | AUFSRR: rotated-hybrid AUFS + Roe (Tchuen, Tiam Kapen & Burtschell, IJNMHFF 26 (2016) 1310-1327… | The rotation destroys the exact staggered constrained transport, and it does so in proportion to the cure's own strength. In the current flux the normal-field row is identically zero by construction and this is load-bearing: f_con… |
| `worse` | Enhanced HLLEM / HLL-CPS with all-Mach anti-diffusion coefficients (Gogoi & Mandal, Phys. Fluid… | The half that is NEW here pushes dissipation the wrong way past a floor this code has already measured, and the half that would be the cure is only one of the two actions the cure needs. (1) Gogoi-Mandal recover low Mach by ADDING… |
| `worse` | Entropy-consistent flux with correct entropy production at shocks (Ismail & Roe, JCP 228(15) 20… | The method's defining condition cannot be written down for this scheme's flux. Ismail-Roe starts from a Tadmor entropy-CONSERVATIVE flux, i.e. a two-point flux satisfying (v_R - v_L).F* = psi_R - psi_L, where v and psi are the ent… |
| `worse` | Entropy-control flux term F_EC (Xie, Tian, Zhang, Yu, Ren/Yang, Comput. Math. Appl. 102 (2021) … | F_EC does not remove this scheme's amplifier, so it cannot replace the current cure — it can only be added on top of it. The instability here is not "insufficient entropy production inside a shock"; the cause is localised and writ… |
| `worse` | Entropy-production floor inside the numerical shock (Xie, Tian, Zhang, Yu, Yang, CAMWA 102 (202… | The floor is blind to the mode, because in this scheme the anti-diffusion acts ONLY on linearly degenerate fields, across which entropy production is zero at leading order. Verified in /home/user/diffusionIMEX/blastwave_debug/patc… |
| `worse` | Entropy-treatment flux aimed directly at the carbuncle (Ismail, Roe & Nishikawa, CFD 2006 / Spr… | The eps-split explicit sub-system has no entropy, so there is nothing for an entropy-consistent flux to be consistent with. This is provable in three lines from the paper's own Jacobian, /home/user/diffusionIMEX/blastwave_debug/pa… |
| `worse` | Fleischmann et al. (JCP 401, 2020, 109004) low-dissipation Roe — modified (limited) acoustic ei… | It has no authority over the term that causes the instability in THIS scheme, and the only setting in which it would acquire authority is a configuration this code has already measured as catastrophic. The instability here was loc… |
| `worse` | Genuinely multidimensional (vertex-based) Riemann solvers — Balsara family (multi-D HLLE/HLLC, … | It cannot remove the detector, which is the entire point of adopting it. Balsara's flattener chi (with its tuned kappa_1 = 0.4) is already in this scheme for a second, independent reason: it stabilises the IMPLICIT magnetic predic… |
| `worse` | H-correction (multidimensional / crossflow dissipation), Sanders, Morano & Druguet, JCP 145(2) … | eta is NOT a wave-speed *difference* that is small in smooth flow -- it is half the WIDTH OF THE RIEMANN FAN. Sanders et al. eq.(10), as implemented in Athena, is eta_face = 0.5*/(u_R + c_f,R) - (u_L - c_f,L)/, which in smooth uni… |
| `worse` | H-correction / multidimensional (crossflow) dissipation (Sanders, Morano & Druguet, JCP 145 (19… | An envelope/wave-speed floor is ALGEBRAICALLY INERT on the linearly degenerate fields in this scheme's HLLEM, which is exactly where the instability lives. In /home/user/diffusionIMEX/blastwave_debug/patch/Final_emb_omp_clean.cpp … |
| `worse` | H-correction as shipped by production MHD codes (Athena, MPI-AMRVAC, Cholla) | An envelope floor cannot reach the linearly degenerate fields at all, so in this flux the H-correction is a verified no-op on the field that carries the carbuncle. The H-correction's entire mechanism is max(/lambda^alpha/, eta): r… |
| `worse` | H-correction on unstructured / triangular grids (modified-H, and the mixed RoeVLPA entropy fix)… | In an HLLEM solver the H-correction's floor lands in the one slot where it is provably a no-op on the fields that carry the striping. The published action is /lambda_k/ -> max(/lambda_k/, eta^H), a floor on a Roe dissipation coeff… |
| `worse` | HLLC-LM (Fleischmann, Adami & Adams, JCP 423 (2020) 109762) — shock-stable HLLC by low-Mach red… | In an HLLEM flux, phi cannot change the dissipation on the linearly degenerate fields AT ALL — it is an exact algebraic no-op — and the author's own E17 run proves those fields carry 100% of this instability. Residual dissipation … |
| `worse` | HLLEMCC (continuous blend of HLLEM anti-diffusion driven by a Rankine-Hugoniot residual) | HLLEMCC has only ONE action (scale delta_k) and your cure has TWO. The repo's own ablation ladder (/home/user/diffusionIMEX/blastwave_debug/PATCHES.md, rows ref/L1-L8/none, blast 128^2, t=0.01) shows the delta_k action is the SMAL… |
| `worse` | Liou's mass-flux design criterion: no pressure-difference term in the numerical mass flux | Enforcing the prohibition removes the ONLY mass-row dissipation left at the carbuncle face instead of restoring the missing one, and the current cure already enforces the criterion in a strictly better way. In /home/user/diffusion… |
| `worse` | Refined H-correction with the current interface removed from its own stencil (Pandolfi–D'Ambros… | eta^H is a DIMENSIONAL speed drawn from the UNSPLIT physical spectrum, and this scheme's explicit time step is not. Two independent consequences, each fatal.  (1) CFL. caldt (Final_emb_omp_clean.cpp:3509-3535) sets dt = CFL/(a_max… |
| `worse` | Rotated Riemann solver for the 2D Euler equations (Levy-Powell-van Leer, JCP 106:201-214, 1993) | A rotation is a change of basis: it multiplies the wave-speed envelope by direction cosines, and cos(theta) * 0 = 0. The cure here must fix TWO failures, and the rotated solver can only reach one of them. Failure (a): on faces tan… |
| `worse` | Rotated Riemann solvers (rotate the 1D problem onto the shock normal) | For a GRID-ALIGNED shock the rotation is the identity, so it cannot touch the instability that was measured. At the faces that CROSS a grid-aligned shock, the velocity-difference vector IS the shock's normal velocity jump, so n1 =… |
| `worse` | Rotated Roe scheme (Ren 2003) | The scheme's numerical flux must be expressible as ONE HLL weight pair times the three physical flux parts, and a rotated solver destroys that. In HLLEM_Update (Final_emb_omp_clean.cpp:1682-1687, 1710-1722) the flux is decomposed … |
| `worse` | Rotated-RHLL / rotated-hybrid Roe-HLLE (Nishikawa & Kitamura 2008) | The whole detector-free premise is an Euler-only Rankine-Hugoniot identity, and it is FALSE for ideal MHD. "Across a grid-aligned shock dq is normal, so the flux becomes pure HLL" holds because the Euler jump conditions give [u_t]… |
| `worse` | Rotated-hybrid HLL/HLLC (Rusanov/HLLC) for relativistic hydrodynamics — Townsend, Könözsy & Jen… | It has no detector — stated as a virtue, but here it is the defect. The rotation is active at every face forever, and its direction is n1 = dv//dv/, where dv at a smooth face is O(h^3) after the WENO3-Z/MP reconstruction. So in sm… |
| `worse` | Rotated-hybrid solvers (Rotated-RHLL: HLL normal to the shock, complete-wave solver across it),… | The switch is a direction built from the reconstructed face jump Δv, and with this scheme's third-order characteristic WENO3-Z+MP reconstruction there is no regime in which that direction is meaningful. In smooth flow the reconstr… |
| `worse` | Rotationally biased upwind differencing (Davis 1984, JCP 56) | It never touches delta_k, which is where this instability actually lives. The measured mechanism in this code (src comment at Final_emb_omp_clean.cpp:1345 and REPORT.md section 6.1) is delta_C = clamp(1 - lam^-/s_L - lam^+/s_R) ->… |
| `worse` | Shock flattening / detector-switched reversion to HLL plus low-order (MinMod) reconstruction (P… | The half of the candidate that is not already your status quo — flag (b), dropping reconstruction to MinMod inside flagged zones — is measured IN YOUR OWN CODE to do the opposite of what the cure is for, and it is order-fatal at t… |
| `worse` | Tadmor's entropy-stable viscosity condition as a dissipation floor (Tadmor, Math. Comp. 49 (198… | The cap it produces is numerically inert at exactly the faces that carbuncle, and I can show this in closed form for THIS flux. Write the scheme's face flux in the form Tadmor's condition is stated in: F = w_L F_L + w_R F_R - /coe… |

## 两条引用上的行动项

1. **Simon & Mandal, "Anti-diffusion control for HLLEM", IJNMF 2019(arXiv:1806.03187)**
   —— 这是唯一的 `equal_tradeoff`,也是**对我们正在处理的那个机制的已发表分析**:
   他们用线性扰动分析把 HLLEM 的激波不稳定定位为"数值激波区内某些反扩散项的错误激活"
   (erroneous activation of certain antidiffusive terms within the numerical shock region)。
   **应当引用它作为机制的出处**,但不能照搬它的限制器:经典 HLLEM 里
   $\delta_k=1-\lambda^-/s_L-\lambda^+/s_R$ 是拿线性退化波速去比**声学**包络
   ($s_L\sim u-c_f$),而本格式比的是**分裂**包络,两者含义不同。

2. **Xie, Li, Li, Tian & Pan, JCP 350 (2017) 607** —— 之前记录它为"强激波附近退回 HLL"的
   现代出处。调查指出它其实是一篇**分析**(把网格对齐不稳定归因于横向数值耗散过小),
   而且它那条"仅靠熵波耗散不能稳定强激波"的否定结论**对本格式不适用**。
   引用时要按"分析/判据"引,不要当成"方法"引。

## 一条机制表述(来自调查,与源码注释一致)

$\delta_C=\mathrm{clamp}(1-\lambda^-_C/s_L-\lambda^+_C/s_R)\to1$ 发生在与激波**相切**的面上,
因为那里熵波速度 $\lambda_C=v_n\approx0$ —— 所以线性退化场的反扩散沿激波面完全没有被衰减。
这正是 flattener 要处理的东西,也是为什么衰减要加在 $\delta_k$ 上而不是加在包络上。
