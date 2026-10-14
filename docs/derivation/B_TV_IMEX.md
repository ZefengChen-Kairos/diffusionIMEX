# B. TV 分裂 IMEX：从要求到离散

目标：$\Delta t$ 只受物质速度限制，$\Delta t=\mathrm{CFL}\,\Delta x/\max|u|$。

- 声波（压力）部分隐式。
- 三种松弛也在同一个隐式闭合里完成。
- 速率 $\mu,\theta,\nu$ 是闭合的参数：
  - 速率 $\to\infty$ 时严格收敛到平衡闭合（渐近保持，AP）；
  - 速率 $\to0$ 时退化为冻结格式。

本章记号同 A 章。一个 IMEX-RK 级内的隐式步长记 $\tau=a_{ii}\Delta t$，$r=\tau/\Delta x$。

本页与 Notion 页 B 同步（2026-10-10，V12）。每一处修改的来由与被否定的做法见“C. 开发记录”（`docs/devlog.md`）。

---

## B0. 要求（每个离散选择都对应其中一条）

| | 要求 | 理由 |
|---|---|---|
| R1 | $m_k,\rho u,E$ 是通量形式 | 守恒；激波速度正确 |
| R2 | 均匀 $(p,u)$、任意 $\alpha_k,\rho_k$ 一步后仍均匀（WB） | 接触间断不被破坏 |
| R3 | 每一相的体积 $\alpha_k$、质量 $m_k$、内能 $\varepsilon_k$ 用**同一个**面通量算子输运 | 输运不改变 $\rho_k=m_k/\alpha_k$、$e_k=\varepsilon_k/m_k$，即不改变 $p_k,T_k$（温度一致性） |
| R4 | 相 $k$ 体积变化 $\Delta V_k$ 与压力功 $-p\,\Delta V_k$ 用同一个 $\Delta V_k$ | 离散 Gibbs 关系 $de+p\,dv=0$ |
| R5 | 隐式压力方程看到的刚度是松弛后的刚度 | 物质 CFL 下一步内 $\alpha$ 可变 $O(1)$；冻结刚度 + 事后松弛的分裂误差 $\propto\Delta t\,c_f$ |
| R6 | 混合内能与守恒总能一致 | SPB 的 “+1”；激波处相内能方程不可靠 |
| R7 | 同一个闭合，速率是参数，两端极限严格成立 | 创新点 5 |

---

## B1. TV 分裂

仿照 Toro–Vázquez-Cendón（`2DEuler_TV.cpp` 的 $f_{\rm conv}=(\rho u,\rho u^2,ku)$、$f_p=(0,p,(\rho e+p)u)$）。每个被输运量的通量 $Xu$ 写成

$$
Xu^{n+1}=\underbrace{X\,u^{e}}_{\text{显式}}+\underbrace{X\,(u^{n+1}-u^{e})}_{\text{隐式增量}},
$$

$u^e$ 是显式级态的速度。

| 方程 | 显式（对流，Rusanov） | 隐式（压力 + 闭合） |
|---|---|---|
| $\alpha_k$ | $\partial_x(\alpha_ku^e)$ | $\partial_x(\hat\alpha_k\Delta u)$，体积应变 $-\alpha_k\partial_xu$，松弛 |
| $m_k$ | $\partial_x(m_ku^e)$ | $\partial_x(\hat m_k\Delta u)$，传质 |
| $\varepsilon_k$ | $\partial_x(\varepsilon_ku^e)$ | $\partial_x(\hat\varepsilon_k\Delta u)$，压力功 $\alpha_kp_k\partial_xu$，松弛 |
| $\rho u$ | $\partial_x(\rho uu^e)$ | $\partial_x\bar p$，$\partial_x(u\,\hat M\Delta u)$ |
| $E$ | $\partial_x(Eu^e)$ | $\partial_x(\bar pu)$，$\partial_x(\hat E\Delta u)$ |

- 显式子系统的特征值全是 $u$。$\alpha,m,\varepsilon$ 是线性输运；$(\rho,\rho u)$ 是无压气体，特征值是二重的 $u$。所以 Rusanov 的耗散速度取物质速度就够了，与声速无关。
- 质量也带隐式增量：纯 TV 里质量只用 $u^e$ 输运，而内能用 $u^{n+1}$，那样 $e_k=\varepsilon_k/m_k$ 每步错 $O(r\Delta u)$，违反 R3。这就是创新点 4（增量通量形式与温度一致性）。

---

## B2. 体积输运形式与 path-conservative Rusanov 的等价

体积方程 $\partial_t\alpha+u\partial_x\alpha=0$ 的 path-conservative 离散用线性路径 $\Psi(s)=Q_L+s(Q_R-Q_L)$：

$$
\int_0^1u(\Psi)\,\frac{d\alpha(\Psi)}{ds}ds=\bar u_f\,(\alpha_R-\alpha_L),\qquad \bar u_f=\tfrac12(u_L+u_R)\quad(\text{面速度}).\tag{B2.1}
$$

恒等式 $\bar u\,\Delta\alpha=\Delta(\alpha u)-\bar\alpha\,\Delta u$（$\bar\alpha=\tfrac12(\alpha_L+\alpha_R)$）说明，PC-Rusanov 涨落

$$
D^\pm_f=\tfrac12\big[\Delta(\alpha u)-\bar\alpha_f\Delta u\big]\pm\tfrac12s\,\Delta\alpha
$$

可以拆成两部分：

1. **守恒的体积通量** $F_f[\alpha]=\tfrac12(\alpha_Lu_L+\alpha_Ru_R)-\tfrac12s\Delta\alpha$。它与 $m_k$、$\varepsilon_k$ 的 Rusanov 通量形式完全相同，满足 R3。
2. **体积应变** $-\bar\alpha_f\Delta u_f$，即离散的 $-\alpha\partial_xu$。

格式保留第 1 部分作为输运（显式 + 隐式增量，两部分用同一离散），把第 2 部分移入隐式闭合。在闭合里它与相能量方程中的压力功 $\alpha_kp_k\partial_xu$ 配对（R4），二者合起来成为离散 Gibbs 关系（B5）。

因此“显式 PC-Rusanov（线性路径，$\bar u$ 取面速度）、隐式部分对 $\alpha$ 用同一离散”与下面的体积输运形式是同一个格式。

---

## B3. 显式算子（MUSCL–Rusanov，耗散速度取物质速度）

在显式级态 $Q^e$ 上，对 $X\in\{\alpha_1,\alpha_2,m_1,m_2,\varepsilon_1,\varepsilon_2,\rho u,E\}$：

$$
F^e_{i+1/2}[X]=\tfrac12\big(X^{L}u^{L}+X^{R}u^{R}\big)-\tfrac12a_f\big(X^{R}-X^{L}\big),\qquad a_f=\max\big(|u^L_f|,|u^R_f|\big)\ \ \text{（默认）或}\ \ a_\Omega=\max_{f}\max\big(|u^L_f|,|u^R_f|\big) ,\tag{B3.1}
$$

**耗散速度：逐面局部（V14）或全场最大，按算例选（`TestSetup::rus`）。** 下面的保正条件是逐面的，两种取法都满足。全场最大是 BB (27) 与 `2DEuler_TV.cpp` 的 `Global_Diss`；它让整个区域都按最快的物质速度耗散：SPB 4.3.2（$\max u\approx3.2\times10^4$ m/s）的稀疏波里 $u$ 从 0 到 $2.2\times10^4$，全场最大时 $L_1(u)=84$，逐面局部时 55（再加 B4.1(b) 的 $\psi^u$ 为 44，加 B5.2 的负亏损分配后 45.7；显式 HLLC 17）；Gresho 涡的动能保持率 0.937 → 0.980。只有 SPB 5.2 用全场最大更好（$L_1(u)$ 118 对 124）。两者都是物质速度，低马赫行为相同（BB §3.6）。

$X^L=X_i+\tfrac12\sigma_{X,i}$，$X^R=X_{i+1}-\tfrac12\sigma_{X,i+1}$，$u^{L,R}=(\rho u)^{L,R}/\rho^{L,R}$。

**耗散速度取面态速度的最大值（单调性条件）。** 线性重构下格均值是两个面值的平均，$X_i=\tfrac12(X^R_{i-1/2}+X^L_{i+1/2})$，所以显式更新可以写成面值的组合：
$$
X^*_i=X^L_{i+1/2}\Big[\tfrac12-\tfrac r2\big(a+u^L_{i+1/2}\big)\Big]+X^R_{i-1/2}\Big[\tfrac12-\tfrac r2\big(a-u^R_{i-1/2}\big)\Big]+\tfrac r2\big(a-u^R_{i+1/2}\big)X^R_{i+1/2}+\tfrac r2\big(a+u^L_{i-1/2}\big)X^L_{i-1/2}.
$$
四个系数都非负时，$X^*$ 是非负面值的非负组合，微量相保正。充分条件是

- $a\ge|u^{L}_f|$ 且 $a\ge|u^{R}_f|$（每个面、两侧的面态），这一条也是必要的：$a<u^R_{i+1/2}$ 时第三个系数为负，取 $X^R_{i+1/2}$ 很大、其余面值很小就得到负值；
- $a\,r\le\tfrac12$。

MUSCL 的面速度 $(\rho u)^{L,R}/\rho^{L,R}$ 不是格速度的凸组合，可能超出 $\max_i|u_i|$，所以 $a_\Omega$ 必须对面态取最大，不能只对格均值取。只对格均值取时，下游面值 $X^R_{i+1/2}$ 以负系数进入本格：水–气界面处面速度只超出 0.34 m/s，界面格的微量气体质量就变负，后向 Euler（`rk=1`）因此卡住（C. 开发记录 V6）。时间步按格速度取（物质 CFL 0.2），$a\,r$ 远小于 $\tfrac12$。

- **体积输运量共用一个限制因子**。记 $\mathcal X_V=(\alpha_1,m_1,m_2,\varepsilon_1,\varepsilon_2)$，格 $i$ 的斜率为
  $$\sigma_{X,i}=\phi_i\,\tfrac12(X_{i+1}-X_{i-1}),\qquad \phi_i=\min_{X\in\mathcal X_V}\frac{\mathrm{minmod}(X_i-X_{i-1},X_{i+1}-X_i)}{\tfrac12(X_{i+1}-X_{i-1})}\in[0,1].$$
  - 对每个量，$|\sigma_X|$ 都不超过它自己的 minmod 斜率且同号，所以仍是 TVD。
  - 所有量的斜率是**同一个**系数乘中心差分，所以量之间的线性关系在面上原样成立。均匀 $p$ 时 $\varepsilon_k=c_k(p)\alpha_k+q_km_k$ 在面值上同样成立（(E1)），这是 WB 所需的。
  - 重构不经过 EOS，所以 $Q^e$ 只需 $\rho>0$；这一点对 IMEX-RK 的外推级态很重要（B10）。
  - 在模板内为常数（到舍入）的量不参与取最小（$0/0$）：$|X_i-X_{i-1}|,|X_{i+1}-X_i|\le10^{-12}\,S_X$，舍入尺度 $S_X$ 见下表；任一量在格 $i$ 取极值时 $\phi_i=0$。
- **速度单独重构（V11）。** 速度 $\mathbf u=(u,v)$ 是单速度模型里独立的原始变量，各分量单独 minmod：$u^{L,R}=u_i\pm\tfrac12\,\mathrm{minmod}(u_i-u_{i-1},u_{i+1}-u_i)$，$u_i=(\rho u)_i/\rho_i$。面上的动量与总能由面值组装：
  $$(\rho u)^{L,R}=\rho^{L,R}u^{L,R},\qquad E^{L,R}=\sum_k\varepsilon_k^{L,R}+\tfrac12\rho^{L,R}|\mathbf u^{L,R}|^2 .$$
  - WB：均匀 $\mathbf u$ 时速度斜率为 0（舍入阈值见下表），面速度就是 $\mathbf u$；加上上一条，(E1) 在面上成立。
  - 保正：minmod 面速度在相邻格速度之间；下面的条件 $a\ge|u^{L,R}_f|$ 照旧对面态取。
  - 为什么不共用：$\phi_i$ 取所有量的最小，任一量的极值使**所有**量降为一阶。二维涡里 $\rho u=-\rho v_\theta y/r$ 沿 $x=0$ 取极值，$v$ 恰在其梯度最大处被降为一阶。共用限制因子的理由（WB）只涉及经 EOS 相关联的 $\alpha,m,\varepsilon$；均匀速度对任何限制器都被精确重构。Gresho 涡（$M=0.01$，$80^2$，一圈）的动能保持率从 0.50 升到 0.94（C. 开发记录 V11）。一维算例基本不变（十二烷平台 $L_1(p)/p^*$ 0.044→0.056，水–气 $L_1(u)$ 2.00→1.87）。
  - 灵感来源：`2DEuler_TV.cpp` 的 `TVDReconstruction` 对原始变量逐个限制；BB (33)–(34) 对每个变量单独 minmod。
  - $\phi_i$ 存 $\min(\phi_i^{V},\ \text{速度的光滑比})$，用来判断“所有被输运量在格 $i$ 都光滑”（B4.1 的二阶声学耗散）。
- $\alpha_2$ 的面值取 $1-\alpha_1$，所以 $\sum_kF^e[\alpha_k]=\tfrac12(u^L+u^R)\equiv\bar u^e_f$，耗散两相抵消。

**限制因子的舍入尺度 $S_X$（V8、V11）。** $S_X$ 是格式产生该量舍入的尺度，不总是该量自身的大小（三个格求和）：

| 量 | $S_X$ | 理由 |
|---|---|---|
| $\alpha_1$ | 1 | $\sum\alpha=1$ 与闭合在 $O(1)$ 上算，微量相 $\alpha\sim10^{-8}$ 的舍入是 $\varepsilon_{\rm mach}$，不是 $\varepsilon_{\rm mach}\alpha$ |
| $u,v$ | $\lvert\mathbf u\rvert+c_f$ | 速度由压力差驱动，舍入 $\sim\varepsilon_{\rm mach}(\tau/h)p/\rho\sim\varepsilon_{\rm mach}c$，与当地 $\lvert u\rvert$ 无关；对向量取模，判据与坐标方向无关 |
| $\varepsilon_k$ | $\sum_k\lvert\varepsilon_k\rvert+\tfrac12\rho\lvert\mathbf u\rvert^2$ | 闭合在两相间搬运 $p\,d\alpha$，舍入相对于混合物内能 |
| $m_k$ | $\lvert m_k\rvert$ | 无传质时 $m_k$ 只被自身通量输运，舍入相对于自身；也是微量相质量保正所需 |

按各量自身大小判时，静止区（$u\approx0$）与二维的横向速度（一维中 $v\equiv0$，二维中是 $10^{-15}$ 量级的噪声）的舍入噪声会触发极值判据，使 $\phi$ 随舍入翻转；二维 y 向均匀的问题因此各行发散（C. 开发记录 V8）。保正要求微量 $\gg10^{-12}S_X$（算例中微量为 $10^{-6}$–$10^{-8}$）。

**中心部分为什么也用重构值（与 BB (35) 不同）。** BB 的中心部分用格均值、只有耗散用重构值。对线性输运 $X_t+uX_x=0$（$u>0$，$a_\Omega=u$），BB 的通量是
$$F=\tfrac12u(X_i+X_{i+1})-\tfrac12u(X^R-X^L)=uX_i+\tfrac14u(\sigma_i+\sigma_{i+1}),$$
含下游格的斜率 $\sigma_{i+1}$，不是 TVD；微量相（$\alpha\sim10^{-8}$）界面处会把 $\alpha$ 推成负值（数值上观察到）。中心部分也用重构值时 $F=u(X_i+\tfrac12\sigma_i)$，就是 MUSCL 迎风通量，TVD。两者都是二阶。

**与 BB 的对应：** BB (36)–(38) 要求，对 $\phi$ 加的人工耗散必须在内能方程里按 $\partial(\rho\phi e)/\partial\phi$ 缩放后同样加上，否则接触处出现压力振荡。本格式对 $\alpha_k$ 与 $\varepsilon_k$ 用同一个 $a_\Omega$ 和同一种重构；均匀 $p$ 时由 (E1) 有 $\Delta\varepsilon_k=\frac{p+\gamma\pi}{\gamma-1}\Delta\alpha_k+q\Delta m_k$，这正是 BB (38) 的缩放，而且对任意 $\alpha,m$ 分布都自动满足。

BB (39)–(40) 的压力跳跃耗散没有单独实现。它的系数是 $\sigma^{**}=\tfrac12\Delta x\,a_\Omega\,\partial(\rho\phi e)/\partial p$，仍是**物质速度**（BB §3.6 末：所有数值粘性都只与物质速度成正比，与 $\epsilon=M^2$ 无关）。用户论文 §3.4 (52)–(54) 与 `2DEuler_TV.cpp` 的 `MaterialSpeedp` 是同一思路：压力方程里加 $\chi|u|\Delta x$ 的隐式 Laplacian，$\chi$ 是激波指示子。本格式的对应物是 B4.1 的声学补足项：系数是声速（$c\Delta x$ 量级），只在时间步分辨声波时打开（$\varphi>0$），速度跳跃部分按散度作用（B4.1(c)）。V8 及以前这里写“与 BB 的原则一致”，并不准确：BB 的耗散是物质速度的；声学速度的耗散是否破坏 AP，要看它作用在什么上，见 B4.1(c)。

显式更新：$X^*_i=X^s_i-r\,\delta_i[F^e[X]]$，$Q^s$ 是级的更新基底。

---

## B4. 隐式部分：面速度、增量通量、输运态

$\rho^*_i=\sum_km^*_{k,i}$，$u^*_i=(\rho u)^*_i/\rho^*_i$。未知量是各格的混合压力 $P_i$。

### B4.1 面速度与面压力

**(a) 压力加速度与动量的一致性。** 格速度 $u_i$ 是主变量，面速度由两侧格速度加上压力加速度给出。面两侧各半格以同一加速度运动，半格动量平衡（质量 $\rho_i\Delta x/2$、$\rho_{i+1}\Delta x/2$）给出
$$
a_f=-\frac{P_{i+1}-P_i}{\bar\rho_f\Delta x},\quad \bar\rho_f=\tfrac12(\rho^*_i+\rho^*_{i+1}),\qquad
p_f=\frac{\rho^*_{i+1}P_i+\rho^*_iP_{i+1}}{\rho^*_i+\rho^*_{i+1}} .
$$
于是格的动量更新 $-r(p_{f,E}-p_{f,W})/\rho_i=\tfrac12\tau(a_E+a_W)$：格速度的压力增量是两侧面增量的平均，动量守恒，均匀 $P$ 时 $p_f=P$（WB）。（若 $p_f=\tfrac12(P_i+P_{i+1})$，轻格会以格密度承受面上的压力差，而面速度用的是面密度，二者不一致：1 GPa 水–气算例中气体格动能与面功不匹配，$(\rho e)^C<0$。）

**(b) 声学迎风耗散。** 对线性声学（$u_0$ 平流、$\rho,c$ 常数）做 von Neumann 分析（`tools/vonneumann_imex.py`；$C=u_0\tau/\Delta x$，$\nu=c\tau/\Delta x$）：

- 低波数：隐式步的等效耗散已经不低于声学 Godunov（例如 $C=0.2,M=0.5$：$d_{\rm eff}=0.16$ 对 $0.08$，单位 $\Delta x^2/\tau$）。所以“后向 Euler 只有 $c^2\tau/2$、$\nu<1$ 时不足”这个说法不对。
- 高波数才是缺口：奇偶速度模态 $u_i=(-1)^i$ 的面平均为零，不产生散度，压力看不见它；$p_f$ 是压力平均，也看不见它。所以中心隐式声学算子**对任何 $\nu$ 都不阻尼它**（$|\lambda(\pi)|=1$，$\nu$ 从 0.05 到 10 都如此）。只有显式的物质速度 Rusanov（BB 的耗散）阻尼它，每步因子 $|1-2C|$；从静止起步时 $a_\Omega=0$，完全不阻尼。
- 强间断处（起步的压力跳跃、阻抗跳跃很大的接触）这个模态被大量激发，表现为界面附近的色散过冲：十二烷算例液体速度过冲 13%、拉出张力、微量蒸汽空化（$p\to0$，与 $\Delta t$ 无关）；水–气算例 CFL $\le0.1$ 时界面格被抽空。

所以在隐式面量上加声学 Riemann 解（Godunov）的迎风项：速度跳跃项阻尼奇偶模态，压力扩散项使压力跳跃处的面速度单调（只加其一都不够：只加速度项时十二烷的空化仍在，只加压力项时水–气的气体过冲 7%）。权重取 $(1-\nu)^+$：
$$
\boxed{\begin{aligned}
u_f&=\bar u^*_f-\frac{\tau(P_{i+1}-P_i)}{\bar\rho_f\Delta x}-\varphi_f\,\psi^P_f\,\frac{P_{i+1}-P_i}{Z_i+Z_{i+1}},\\
p_f&=\frac{\rho^*_{i+1}P_i+\rho^*_iP_{i+1}}{\rho^*_i+\rho^*_{i+1}}-\varphi_f\,\psi^u_f\,\frac{Z_iZ_{i+1}}{Z_i+Z_{i+1}}\,D^*_f,
\end{aligned}}\qquad \varphi_f=\Big(1-\frac{\tau\min(c_i,c_{i+1})}{\Delta x}\Big)^+ ,\tag{B4.1}
$$
$Z=\rho^*c$，$\bar u^*_f=\tfrac12(u^*_i+u^*_{i+1})$。一维中 $D^*_f=u^*_{i+1}-u^*_i$；二维中 $D^*_f$ 是面两端两个顶点散度的平均乘 $\Delta x$（B4.1(c)，V11），一维与 y 向均匀时逐位退化为 $u^*_{i+1}-u^*_i$。

- $\tau\to0$：$\varphi=1$，$u_f-\bar u^*=-\Delta P/(Z_i+Z_{i+1})$，$p_f-\bar p=-\frac{Z_iZ_{i+1}}{Z_i+Z_{i+1}}\Delta u$，就是声学 Riemann 解（Godunov）。
- $c\tau\ge\Delta x$：$\varphi=0$，格式不变。固定 $\Delta t,\Delta x$ 令 $M\to0$ 时 $c\tau/\Delta x\to\infty$，补足项恒为零，极限格式与不加补足项时相同。
- $\varphi>0$ 当且仅当时间步分辨声波（$c\tau<\Delta x$，即 $M>a_{ii}\mathrm{CFL}$，或起步斜坡）。$a_{ii}\mathrm{CFL}<M<1$ 这一段与起步斜坡里它是声速量级的耗散，是否破坏低马赫精度见 (c)。
- 加上后奇偶模态在 $\theta=\pi$ 的放大因子为 $|(1-2\varphi\nu)(1-2C)|$（例如 $C=0.2,M=0.5$：0.31，Godunov 0.60）。
- 速度跳跃项是显式粘性，$\nu_{\rm add}\tau/\Delta x^2=\tfrac12\nu(1-\nu)\le\tfrac18$，稳定；压力项对 $P$ 隐式、线性。
- **$c$ 取冻结声速。** $\rho c_{\rm frozen}\ge\rho c$ 对所有松弛闭合（Wood、$pT$、$pTg$）成立（次特征条件），所以补足后的耗散对任何闭合都不低于其波的 Godunov 耗散。用 Wood 声速试过：数值混合区里 $Z$ 很小，$\Delta P/(Z_i+Z_{i+1})$ 很大，把界面推得过快（水–气 $\max u/u^*=1.11$，激波超前 20 格），舍弃。
- **二阶（V14 起两项分开）：** 压力项 $\psi^P_f=1-(1-\psi^{\rm mm}_f[P])\min(\phi_i,\phi_{i+1})$，速度项 $\psi^u_f=\psi^{\rm mm}_f[u]$。速度项是体粘性（B4.1(c)），只应由速度场决定：乘共用限制因子时，空化扇区里 $\alpha,m,\varepsilon$ 的极值（含微量相的舍入）使它在整个膨胀区变成一阶体粘性，把膨胀的动能变成内能（P-S ±500 m/s $pTg$ 中心压力 $7.8\times10^4\to2.6\times10^4$ Pa，SPB 4.5 扇区 $L_1(u)$ 0.066 → 0.045，C. 开发记录 V14）。下面关于接触的理由只涉及压力项。V13 及以前两项都用 $\psi_f=1-(1-\psi^{\rm mm}_f)\min(\phi_i,\phi_{i+1})$，$\psi^{\rm mm}_f=1-\frac{\sigma_i+\sigma_{i+1}}{2D_f}$（$D_f$ 为面跳跃，$\sigma$ 为 minmod 斜率，所以 $\psi^{\rm mm}\in[0,1]$，光滑区 $O(\Delta x)$）；$P$ 的斜率取上一次的压力场（Newton 初猜），$u$ 的取 $u^*$，残量对 $P$ 仍是线性的。$\phi$ 是 B3 的共用限制因子：只有所有被输运量都光滑时才降到二阶。接触间断处 $p,u$ 连续（$\psi^{\rm mm}\to0$）但阻抗跳跃，若只看 $p,u$，耗散在最需要的地方被关掉（水–气 CFL 0.05 时界面格气体质量变负，数值上观察到）。
- $c$、$Z$ 在显式预测态 $X^*$ 上计算，一级内不变（残量对 $P$ 仍是线性的）；$X^*$ 不可容许（某相 $\alpha^*$、$m^*\le0$ 或 $p+\pi\le0$）时取 $c=\infty$，即 $\varphi=0$，该面只剩半格动量平衡。固壁面 $u_f=0$，增量为零。
- **舍弃的另一条路：** 让两侧参与面动量平衡的质量取 $\rho\min(c\tau,\Delta x/2)$（声学影响域，不超过半格），$\tau$ 小时同样是 Godunov，$\tau$ 大时是带非弹性碰撞项的半格平衡。但大步长时它每一级都把格速度换成两侧面速度的平均（$\tfrac14\delta^2u$），耗散率 $\propto\Delta x^2/\tau$ 随 $a_{ii}$ 变化（LSDIRK2 是后向 Euler 的 3.4 倍），空化管的熵误差因此成倍增大。

**(c) 低马赫：两项的量级，以及速度项的散度形式（V11）。**

低马赫数 $M=U/C$，well-prepared 数据：压力脉动 $\delta p\sim\rho U^2$，$\nabla\cdot\mathbf u=O(M^2U/L)$（$p_0$ 不变的封闭或周期区域），光滑场的面跳跃 $\sim U\Delta x/L$。

- **压力扩散项**：$u_f$ 的修正 $\varphi\psi\,\Delta P/(Z_i+Z_{i+1})\sim\varphi\psi\,\rho U^2(\Delta x/L)/(2\rho C)=\varphi\psi\,\tfrac M2\tfrac{\Delta x}LU$，相对于 $U$ 是 $O(M)$，$M\to0$ 时消失。它不进动量方程，只抹平 $P$ 的跳跃（声波）。
- **速度跳跃项（逐维形式）**：$p_f$ 的修正 $\varphi\psi\,\tfrac Z2\Delta u\sim\varphi\psi\,\tfrac{\rho C}2U\tfrac{\Delta x}L$，与物理压力脉动 $\rho U^2$ 之比为 $\varphi\psi\,\frac{\Delta x/L}{2M}$，即 $O(1/M)$。这就是 Godunov/Roe 型格式的低马赫缺陷：动量方程里出现 $c\Delta x$ 量级的粘性，离散压力脉动是 $O(M)$ 而不是 $O(M^2)$。
- **一维里没有这个问题。** 一维的不可压极限是 $u$ 均匀，$\Delta u=O(M^2)U\Delta x/L$，速度跳跃项同样是 $O(M)$。问题只出在多维：逐维的法向跳跃 $u_{i+1,j}-u_{i,j}\approx\Delta x\,\partial_xu$ 对无散的涡不为零（$\partial_xu=-\partial_yv\neq0$）。
- **声学子系统只通过散度看见速度**（$\partial_t\mathbf u+\nabla p/\rho=0$，$\partial_tp+\rho c^2\nabla\cdot\mathbf u=0$）：涡（无散部分）在声学算子的核里，声学耗散不应作用在它上面。所以把速度跳跃换成散度。x 面 $(i+\tfrac12,j)$：
$$
D^*_{i+1/2,j}=\Delta x\cdot\tfrac12\big(\mathrm{div}_{i+1/2,j+1/2}+\mathrm{div}_{i+1/2,j-1/2}\big),\qquad
\mathrm{div}_{i+1/2,j+1/2}=\frac{u_{i+1,j}+u_{i+1,j+1}-u_{i,j}-u_{i,j+1}}{2\Delta x}+\frac{v_{i,j+1}+v_{i+1,j+1}-v_{i,j}-v_{i+1,j}}{2\Delta y},
$$
展开即
$$
D^*_{i+1/2,j}=\tfrac14\big[2(u_{i+1,j}-u_{i,j})+(u_{i+1,j+1}-u_{i,j+1})+(u_{i+1,j-1}-u_{i,j-1})\big]+\frac{\Delta x}{4\Delta y}\big[(v_{i,j+1}-v_{i,j-1})+(v_{i+1,j+1}-v_{i+1,j-1})\big],
$$
y 面对称。$\mathrm{div}_{i+1/2,j+1/2}$ 是四个格心值双线性插值在顶点处的散度；$u,v$ 取 $X^*$ 的格速度。

性质：

1. **一维不变。** 一维与 y 向均匀时 $D^*=u_{i+1}-u_i$；代码按 $\tfrac14\big(2d+(d+d)\big)$ 计算，逐位相等。所以所有一维结果与旧的补足项相同（V3 引入它的理由——奇偶模态、十二烷伪空化、水–气界面——都是一维的，照旧成立）。
2. **体粘性（grad-div）。** 记 $q_f=\varphi\psi\tilde Z_fD^*_f$，x 动量得到 $+r_x(q_{i+1/2,j}-q_{i-1/2,j})$。$\tilde Z$ 均匀、$\Delta x=\Delta y=h$、周期边界时分部求和：
$$\sum_{ij}\Big[u_{ij}\big(q_{i+1/2,j}-q_{i-1/2,j}\big)+v_{ij}\big(q_{i,j+1/2}-q_{i,j-1/2}\big)\Big]=-\varphi\psi\tilde Z\,h^2\sum_{\rm 顶点}\mathrm{div}_v^2\le0 ,$$
动能只按顶点散度的平方耗散，顶点意义下的无散场在它的核里。激波处 $\mathrm{div}$ 很负，它就是 von Neumann–Richtmyer 型的体粘性。
3. **低马赫量级（对所有 $M$ 一致）。** 光滑无散场的顶点散度是 $O(M^2U/L)$ 加上 $O(\Delta x^2U/L^3)$ 的截断误差，所以 $p_f$ 的修正与 $\rho U^2$ 之比为 $\varphi\psi\big[O(M\Delta x/L)+O\big((\Delta x/L)^3/M\big)\big]$。物质 CFL 定步长时 $\varphi=(1-a_{ii}\mathrm{CFL}/M)^+$，而
$$\sup_M\frac{\varphi}{M}=\sup_{M>m_0}\frac{M-m_0}{M^2}=\frac1{4m_0},\qquad m_0=a_{ii}\mathrm{CFL}$$
（$M=2m_0$ 处取到；LSDIRK2、CFL 0.2 时 $m_0=0.059$，$1/(4m_0)=4.3$）。所以这个比值对所有 $M$ 不超过 $\psi\big[O(M\Delta x/L)+4.3\,O((\Delta x/L)^3)\big]$，与 $M$ 无关地小。逐维形式（V10）是 $\varphi\psi\frac{\Delta x/L}{2M}\le\frac{\psi}{8m_0}\frac{\Delta x}L=2.1\,\psi\,\frac{\Delta x}L$：一阶、常数大，在 $M=2m_0\approx0.12$ 处最大，这就是 Gresho $M=0.1$ 的 0.36。
4. **奇偶模态照旧被阻尼。** $u_{ij}=(-1)^i$ 时 $D^*=\pm2$（法向紧致差分权 $\tfrac12$，两侧行各 $\tfrac14$）。

Gresho 涡（$80^2$，一圈，单相极限；C. 开发记录 V11）的动能保持率：

| | $M=0.1$ | $M=0.01$ | $M=0.001$ |
|---|---|---|---|
| 显式 HLLC（声学 CFL 0.4） | 0.83 | 0.47 | — |
| IMEX V10（逐维速度项、共用限制因子） | 0.36 | 0.50 | 0.50 |
| IMEX V11 | 0.93 | 0.94 | 0.94 |

$M=0.1$ 与 $M=0.01$ 之间剩下的 1% 来自压力扩散项（$O(M)$，见第一条）：把它也关掉时 $M=0.1$ 为 0.937。

**试过而舍弃的开关（V11）：**

| 速度项乘 | 压力项乘 | Gresho $M=0.1$ | 十二烷平台 $L_1(p)/p^*$（V10 0.044） | 其他一维 |
|---|---|---|---|---|
| 局部马赫数 $\theta=\min(1,\max\lvert\mathbf u\rvert/\min c)$ | 1 | 0.85 | 0.24 | SPB 4.5：1024 次拒绝 |
| 激波指示子 $\chi$（用户论文 (52)） | 1 | 0.93 | 0.21 | P-S ±500 m/s（只 $p$）：$\Delta t$ 下溢失败 |
| $\chi$ | $\chi$ | 0.94 | 0.39 | 同上 |
| $\max(\theta,\chi)$ | 1 | 0.85 | — | SPB 4.5：1024 次拒绝 |
| 两项都去掉 | | — | 0.64；最低 $p=3.8\times10^3$ | P-S ±500 m/s 失败；水–气 654 次拒绝、$\max u/u^*=1.055$ |

- $\theta$ 是 Dellacherie 的 all-Mach 修正与 Rieper 的 LMRoe 的做法；$\chi$ 是 Balsara 平坦器（用户论文 §3.4，`2DEuler_TV.cpp` 的 `ShockDetector`，$k_1=10^{-3}$，$3\times3$ 取最大）。
- 它们在二维 Gresho 上都有效，但在一维也把速度项关掉了。十二烷的失败在液体稀疏波尾之后的平台里：那里 $\nabla\cdot u\ge0$（$\chi=0$），$M\approx0.1$（$\theta\approx0.1$），$\varphi\approx0.5$。液体刚性，$\delta u=0.3$ m/s 就使平台压力（$1.87\times10^5$ Pa）降到零，这里需要声学耗散。
- 而一维里速度项本来就不破坏 AP（第三条）。所以开关关错了地方；散度形式只在需要的地方（多维的无散部分）去掉它。

灵感来源：Guillard–Viozat (1999) 与 Dellacherie (2010) 的低马赫分析（缺陷来自速度跳跃项，压力跳跃项无害）；Dellacherie–Omnes–Rieper (2010)（笛卡尔网格几何的作用）；Morton–Roe (2001) 的 vorticity-preserving 格式与 Barsukow (2019) 的 stationarity-preserving 格式（用顶点/格点平均的散度构造多维耗散，低马赫相容）；von Neumann–Richtmyer 体粘性；BB §3.6（耗散不随 $\epsilon$ 增长）；用户的提问与论文 §3.4（激波指示子，试过）。

**增量**（参考速度就是显式中心通量里用的那个面速度）：

$$
\Delta u_f=u_f-\bar u^e_f,\qquad \bar u^e_f=\tfrac12(u^L+u^R)^e .\tag{B4.2}
$$

### B4.2 增量通量：总通量对 $u_f$ 迎风

一步内面扫过的体积是 $u_f\tau$，其中的物质来自 $u_f$ 的上游格。显式部分已经按 $\bar u^e_f$ 从**它的**上游格扫过 $\bar u^e_f\tau$（Rusanov 在 $a_\Omega\ge|u|$ 时就是迎风加耗散），增量只需补足差额：
$$
I_f[X]=\hat X_{i}\big[(u_f)^+-(\bar u^e_f)^+\big]+\hat X_{i+1}\big[(u_f)^--(\bar u^e_f)^-\big],\qquad X\in\{\alpha_k,m_k,\varepsilon_k\},\tag{B4.2'}
$$
$(\cdot)^\pm=\max/\min(\cdot,0)$，$\hat X_j=X^*_j/\sum_l\alpha^*_{l,j}$ 是格 $j$ **单位体积**的内容。

- 同号（$u_f$ 与 $\bar u^e_f$ 同向）：$I_f=\hat X_{\rm up}\Delta u_f$。$\Delta u_f<0$ 时退回的是显式部分扫过来的**同一种**物质。
- 反向：退回显式扫过的全部，再从另一侧扫 $u_f\tau$。
- $\sum_k\hat\alpha_{k,j}=1$，两项系数之和为 $\Delta u_f$，所以 $\sum_kI_f[\alpha_k]=\Delta u_f$，(B4.4) 不变；同一系数用于 $\alpha_k,m_k,\varepsilon_k$，R3 不变；$I_f$ 对 $P$ 连续（分段线性）。
- **旧写法的错误：** 若按 $\Delta u_f$ 的符号取供体，$\Delta u_f<0$ 而 $u_f>0$ 时，增量从**下游**格取物质。显式部分把 L 的物质扫进 R，增量却从 R 退回 R 的物质，等于在面上交换两种物质：单速度模型里不该有的相间逆流（液–汽界面处把汽送进液体）。
- 显式通量中 $a_\Omega>|u|$ 的那部分是 Rusanov 的数值耗散，保持不变。

**为什么按单位体积归一。** 显式输运之后，格 $j$ 内的物质占体积 $\sum_l\alpha^*_{l,j}=1-r\delta_j[\bar u^e]\neq1$；增量通量搬走的是体积 $r\Delta u_f$ 的这种物质，所以按单位体积归一。于是 $\sum_k\hat\alpha_{k,j}=1$，(B4.4) 严格成立；三个量用同一个归一因子，比值 $\rho_k,e_k$ 不变（R3）；$\Delta u_f=O(\tau)$，一阶的供体内容在光滑区只带来 $O(\Delta x\,\tau)$ 的误差。

**动量与总能的增量。** 记 $M_{L}=\sum_k\hat m_{k,i}\big[(u_f)^+-(\bar u^e_f)^+\big]$、$M_R=\sum_k\hat m_{k,i+1}\big[(u_f)^--(\bar u^e_f)^-\big]$，$\hat M_f=M_L+M_R$。增量质量携带新时间层的法向面速度 $u_f$，切向速度取各自供体格的 $v^*$。

输运态为

$$
X^{\rm tr}_i=X^*_i-r\,\delta_i[I_f[X]],\qquad X\in\{\alpha_k,m_k,\varepsilon_k\}.\tag{B4.3}
$$

**体积应变**（显式中心通量与增量合起来恰好是 $u_f$；散度对压力完全隐式）：

$$
s_i\equiv1-\sum_k\alpha^{\rm tr}_{k,i}=r\,\delta_i[\bar u^e_f+\Delta u_f]=r\,\delta_i[u_f].\tag{B4.4}
$$

**动量与总能**（守恒，R1），$p_f$ 取 (B4.1)：

$$
\begin{aligned}
(\rho u)^{n+1}_i&=(\rho u)^*_i-r\,\delta_i[p_f]-r\,\delta_i[u_f\hat M_f],\\
E^{n+1}_i&=E^*_i-r\,\delta_i\Big[\textstyle\sum_kI_f[\varepsilon_k]+\tfrac12u_f^2\hat M_f+\tfrac12(v_L^{*2}M_L+v_R^{*2}M_R)+p_fu_f\Big].
\end{aligned}\tag{B4.5}
$$

- 增量质量携带新时间层面速度 $u_f$：在均匀加速度场中 $u$ 保持均匀。
- 守恒的混合内能：

$$
\rho^{n+1}_i=\sum_km^{\rm tr}_{k,i},\qquad (\rho e)^C_i=E^{n+1}_i-\frac{\big((\rho u)^{n+1}_i\big)^2}{2\rho^{n+1}_i}.\tag{B4.6}
$$

$X^{\rm tr}_i$、$s_i$、$(\rho e)^C_i$ 都只依赖 $(P_{i-1},P_i,P_{i+1})$。声学补足中被耗散的动能经 $E$ 进入 $(\rho e)^C$，由 B5.2 的 $\delta$ 分给各相（熵增，方向正确）。

---

## B5. 统一闭合

### B5.1 从模型到离散 Gibbs 关系

对 (A1.1) 的隐式部分在 $[0,\tau]$ 上做后向 Euler。体积应变 $s$ 在两相间按冻结比例 $\omega_k=\alpha^{\rm tr}_k/\sum_j\alpha^{\rm tr}_j$ 分配，这是 $\alpha_k\partial_xu$ 里的 $\alpha_k$；同一个 $\omega_k$ 也用在压力功 $\alpha_kp_k\partial_xu$ 里（R4）。记 $\chi_1=+1$，$\chi_2=-1$：

$$
\begin{aligned}
\alpha_k&=\alpha^{\rm tr}_k+\omega_ks+\chi_kR,\qquad R=\tau\big(\mu(p_1-p_2)+\mathcal M/\rho_I\big),\\
m_k&=m^{\rm tr}_k+\chi_k\tau\mathcal M,\\
\varepsilon^{\rm nc}_k&=\varepsilon^{\rm tr}_k-\omega_kp_k\,s-\chi_k\,\tau p_I\mu(p_1-p_2)+\chi_k\tau\big(\mathcal Q+e_I\mathcal M\big).
\end{aligned}\tag{B5.1}
$$

由 $\sum_k\omega_k=1$ 得 $\sum_k\alpha_k=\sum\alpha^{\rm tr}+s=1$（B4.4），体积约束自动满足。

### B5.2 混合能量一致性（“+1” 进入闭合，R6）

$\sum_k\varepsilon^{\rm nc}_k=\sum_k\varepsilon^{\rm tr}_k-s\sum_k\omega_kp_k$ 是非守恒的混合内能，与 (B4.6) 的守恒值 $(\rho e)^C$ 差

$$
\delta=(\rho e)^C-\sum_k\varepsilon^{\rm tr}_k+s\sum_k\omega_kp_k .\tag{B5.2}
$$

$\delta$ 与松弛项无关（$\chi_k$ 项两相相消）。光滑流里它是离散误差（正负交替的二阶小量），激波处是 $O(1)$ 的物理耗散（格式的数值粘性把动能变成内能）。分给两相：

$$
\varepsilon_k=\varepsilon^{\rm nc}_k+\Lambda_k(\delta),\qquad \sum_k\Lambda_k=\delta .\tag{B5.3}
$$

于是 $\sum_k\varepsilon_k=(\rho e)^C$ 恒成立，不需要级末投影，解出来的状态就是混合能量一致的（P-S 2014 的定义 3.1）。显式程序的力学松弛也用同一个闭合（$s=0$），不再做 SPB 3.4 的重初始化（A6）。$\Lambda_k$ 的形式由两条原则定（V16）：**物理耗散按模型自己的激波关系分，离散误差按不伤害任何一相的方式分。**

**(a) 分配决定格式的极限激波关系。** 6 方程/Kapila 模型的相能量方程是非守恒的。穿过激波时守恒律只定混合物的 Hugoniot $e-e_0+\bar p(v-v_0)=0$（$\bar p=\tfrac12(p_0+p_1)$）；耗散在两相之间怎么分，无粘模型本身不定，由激波内部的耗散机制决定——数值上就是由格式的粘性路径决定（路径守恒格式一般收敛不到所选路径的弱解，Abgrall–Karni 2010，原因相同）。模型的激波闭合是 SPB (II.5)（Saurel, Le Métayer, Massoni, Gavrilyuk 2007）：每相沿自己的 Hugoniot，

$$
e_k-e_{k0}+\bar p\,(v_k-v_{k0})=0\qquad(k=1,2),\tag{II.5}
$$

即 $(p,v_1,v_2)$ 空间里 DLM 直线路径的弱解；它用 epoxy–spinel 的实验激波速度验证过，`tests/exact5.py` 与 SPB 5.1/5.2 的精确解都用它。所以格式必须收敛到 (II.5)，而 (B5.3) 的 $\Lambda_k$ 是唯一的调节手段。只有一条正确的规则，就是让极限激波满足 (II.5) 的那条；下面推导它。

**(b) 剖面分析。** 激波系下的定常剖面：质量通量 $j=\rho u$，总正应力 $P=p-\sigma_{\rm num}$（$\sigma_{\rm num}$ 是格式数值粘性的应力）。动量与能量守恒给出 Rayleigh 线 $P(v)=p_0+j^2(v_0-v)$ 与混合物的 $\mathrm de=-P\,\mathrm dv$。相 $k$（单位相质量，$v_k=1/\rho_k$，$Y_k$ 质量分数）：

$$
\mathrm de_k=-p\,\mathrm dv_k+\mathrm dq_k,\qquad Y_k\,\mathrm dq_k=\lambda_k\,\mathrm dq,\qquad \mathrm dq=(p-P)\,\mathrm dv\ge0,
$$

$\lambda_k=\partial\Lambda_k/\partial\delta$ 是局部分配系数。相加得 $\mathrm de=-P\,\mathrm dv$：无论怎么分，混合物都落在自己的 Hugoniot 上（这就是 $\sum_k\varepsilon_k=(\rho e)^C$ 的作用），分配只决定终态的相分配 $(v_{11},v_{21})$。(II.5) 等价于

$$
Y_kq_k=\int_{\text{剖面}}(p-\bar p)\,\mathrm dv_k\qquad(k=1,2),\tag{B5.4}
$$

相 $k$ 得到的热 $=$ 它沿剖面的压缩功与 Hugoniot 功之差。
- **弱激波。** 相 $k$ 沿剖面近似等熵，$v_k(p)=v_{k0}\big[1-x/K_k+\mathcal G_kx^2/K_k^2+\dots\big]$，$x=p-p_0$，$K_k=\rho_kc_k^2$，$\mathcal G_k=\tfrac12(1+\partial K_k/\partial p|_s)$ 是基本导数（SG 为 $(\gamma_k+1)/2$）。代入 (B5.4)，每单位混合体积 $Y_kq_k\rho=\alpha_k\mathcal G_k\Delta p^3/(6K_k^2)$——正是相 $k$ 单相激波的熵产 $T_k\Delta s_k$，与路径无关。所以弱激波下与 (II.5) 相容的局部规则**唯一**：$\lambda_k\propto\alpha_k\mathcal G_k/K_k^2$。按质量 $Y_k$（V10–V15）、SPB 3.4 重初始化的 $\alpha_k/(\gamma_k-1)$、体积变化份额 $\alpha_k/K_k$ 都不满足。
- **强激波。** $q_k=\int T_k\,\mathrm ds_k$ 依赖剖面，局部规则不可能精确。把剖面 ODE 用常数份额 $\lambda_1=c$ 积分、解出使终态满足 (II.5) 的 $c$（`scratchpad/profile.py`），$c$ 被**终压处的 Hugoniot–等熵能隙**逼近得很好：

$$
Q_k=m_k\,(v^{\rm pre}_k-b_k)\,f_k(p_0\to p_1),\qquad
f_k=\frac{(p_1+\gamma_k\pi_k)(\tilde v_H-\tilde v_S)}{\gamma_k-1},\quad
\tilde v_H=\frac{(\gamma_k-1)x_1+(\gamma_k+1)x_0}{(\gamma_k+1)x_1+(\gamma_k-1)x_0},\quad
\tilde v_S=\Big(\frac{x_0}{x_1}\Big)^{1/\gamma_k},\tag{B5.5}
$$

$x_i=p_i+\pi_k$；$f_k$ 对 $v-b$ 齐次，$(v^{\rm pre}_k-b_k)f_k$ 是每单位质量的能隙，乘本格的质量 $m_k$ 得每单位本格体积（NASG）。MG 相用 (A2.5) 的 Hugoniot 与等熵线 $(p-p_r)v^{\Gamma+1}=$ 常数，各在 $(0,v_0]$ 上二分（`hug_gap_mg`）。$\lambda^{\rm II.5}_k=Q_k/\sum_jQ_j$，弱极限自动回到 $\alpha_k\mathcal G_k/K_k^2$。

剖面 ODE（SG，$p_0=10^5$ Pa）给出的第一相份额，$c$ 为精确值：

| 混合物 | $p_1$ (Pa) | $c$ | $Q_k$ 能隙 (B5.5) | 弱规则（$p_0$ 处） | 质量 $Y_k$ | SPB $\alpha_k/(\gamma_k-1)$ |
|---|---|---|---|---|---|---|
| epoxy–spinel | $10^8$ | 0.9983 | 0.9983 | 0.9983 | 0.325 | 0.389 |
| epoxy–spinel | $10^{10}$ | 0.9933 | 0.9923 | 0.9983 | 0.325 | 0.389 |
| epoxy–spinel | $2\times10^{11}$ | 0.820 | 0.846 | 0.9983 | 0.325 | 0.389 |
| 水 + 空气，$\alpha_{\rm a}=0.5$ | $10^8$ | 0.0019 | 0.0013 | 0.0000 | 0.999 | 0.105 |
| 水 + 空气，$\alpha_{\rm a}=0.5$ | $10^{10}$ | 0.145 | 0.113 | 0.0000 | 0.999 | 0.105 |
| 水 + 空气，$\alpha_{\rm a}=10^{-3}$ | $10^8$ | 0.582 | 0.559 | 0.0000 | 1.000 | 0.992 |
| 水 + 空气，$\alpha_{\rm a}=10^{-3}$ | $10^{10}$ | 0.994 | 0.992 | 0.0000 | 1.000 | 0.992 |

同一激波速度下各规则终态 $p_1$ 相对 (II.5) 的偏差：质量分配 epoxy–spinel $10^{10}$ 时 $+7.0\%$（$\rho_1/\rho_0$ 1.299 对 1.274）、$2\times10^{11}$ 时 $+9.8\%$，水–空气（0.5）$10^6$ 时已 $+9.9\%$、$10^8$ 时 $+11.7\%$；局部弱规则（$K_k$ 取剖面当地值）在 $p_1\le10^9$ 时偏差 $10^{-8}$–$10^{-3}$，但微量气体（$10^{-3}$，$10^{10}$）时 $+6\%$、气体温度比 (II.5) 高 50 倍（气体受热后体积不缩，份额 $\alpha_k\mathcal G_k/K_k^2$ 一直很大，正反馈）；能隙份额按构造与 $c$ 一致。

**(c) 局部激波强度。** $Q_k$ 需要 $p_0\to p_1$。数值激波宽 3–5 格，所以取 $\pm2$ 格模板（二维为十字）上 $X^*$ 混合压力 (A2.3) 的极值 $p_{\rm lo},p_{\rm hi}$，波前的相比容 $v^{\rm pre}_k=\alpha_k/m_k$ 取模板里压力最低的格、质量 $m_k$ 取本格自己的（激波穿过接触/混合区时被压缩的是本格的物质；组分均匀时与用波前格的 $\alpha_k$ 份额相同），只在压缩格（$u^*$ 的中心差分散度 $<0$）上算；否则 $Q_k=0$（`shock_budget`，`stencil_state`）。$Q_k$ 按级冻结（由 $X^*$ 算）。一次激波通过一格，各步亏损之和恰为物理耗散 $Q=\sum_kQ_k$（混合物落在 Hugoniot 上），所以每步 $\delta\le Q$；$\delta$ 超出 $Q$ 的部分一定是误差，光滑压缩区 $Q\propto\Delta p^3$ 几乎为零。

**(d) 离散误差的分配。** 按各相在可容许界之上的能量分（V14 的负亏损规则推广到两个符号）：

$$
\lambda^E_k=\frac{E_k}{\sum_jE_j},\qquad E_k=\varepsilon^{\rm tr}_k-m^{\rm tr}_kq_k-(\alpha^{\rm tr}_k-m^{\rm tr}_kb_k)\pi_k=\frac{(p_k+\pi_k)(\alpha_k-m_kb_k)}{\gamma_k-1}\Big|_{\rm tr}
$$

（MG 相取冷曲线之上的热能 $\varepsilon_k-m_k\varepsilon_r$）。各相的 $p_k+\pi_k$ 同乘 $1+\delta/\sum_jE_j$：两相的相对误差相同；$\delta<0$ 时扣完仍可容许（V14）；$\delta>0$ 时微量相不会被加热。后一条是必要的：水中 $10^{-6}$ 的空气，能量 $0.25$ J/m³，水 $\sim10^9$ J/m³，水能量 $10^{-9}$ 的相对误差就等于气体的全部能量，任何 $O(1)$ 的份额都会把它烤掉——把弱激波规则用于全部 $\delta\ge0$ 时，SPB 4.5 空化管的空气到 $3\times10^5$ K、4.3.1 界面速度过冲 5%（V16 §16.3，与 V10 的 SPB 重初始化同一机制）。质量分配保护微量相也是因为份额小，但它在激波处是错的（上表）。两相都在可容许界上（$\sum E=0$，退化）时按质量。

**(e) 合成。**

$$
\Lambda_k(\delta)=\begin{cases}\min(\delta,Q)\,\lambda^{\rm II.5}_k+\big(\delta-\min(\delta,Q)\big)\lambda^E_k,&\delta\ge0,\\ \delta\,\lambda^E_k,&\delta<0,\end{cases}\qquad Q=\sum_kQ_k .\tag{B5.6}
$$

$\Lambda_k$ 是 $\delta$ 的连续分段仿射函数（分段点 $\delta=0$、$\delta=Q$）。在 (B6.2) 里 $\lambda_k$ 只通过 $A_k+B_kp=\varepsilon^{\rm tr}_k-m_kq_k+p\alpha^{\rm tr}_k+\Lambda_k(\delta(p))$ 进入，$\delta(p)=D_0+ps$，$D_0=(\rho e)^C-\sum_k\varepsilon^{\rm tr}_k$；闭合在 $\delta$ 所在的分段上写 $\Lambda_k=\lambda_k\delta+\xi_k$（$\xi_k=Q_k-Q\lambda^E_k$ 只在第三段非零），即 $A_k=\varepsilon^{\rm tr}_k-m_kq_k+\lambda_kD_0+\xi_k$，$B_k=\alpha^{\rm tr}_k+\lambda_ks$，$p$ 闭合（B6.1）在每个试探 $p$ 处按 $\delta(p)$ 选段，所以闭合内 $\sum_k\varepsilon_k=(\rho e)^C$ 对任意 $p$ 严格成立；其余闭合按 $\delta(p^{\rm guess})$ 选段（$p^{\rm guess}$ 为本格当前的全局迭代值，收敛时相同）。所以 $G(P)$ 连续、分段光滑，B8 不再需要分支校正。各相的 $\alpha_k(p)$ 不必单调：激波处液格里分到 $\lambda\approx1$ 的微量气体在 $D_0\approx-ps$ 时 $N_k<0$，但 $\Phi=\sum\alpha_k-1$ 由液体主导仍递减；闭合只要求 $\Phi$ 有带区间的根、解处各相可容许（C. 开发记录 V16 §16.8）。（V15 的不连续来自按 $D_0$ 而不是 $\delta$ 的符号分支：$D_0=0$ 处 $\delta=ps\ne0$，$\lambda_k\delta$ 跳变。按 $\delta$ 分支时跳点处 $\lambda_k\delta=0$。）`defect_split` 返回 $(\lambda_k,\xi_k)$；有限 $\mu$ 的两压力闭合 (B6.3) 用同一分段，$\xi_k$ 加在右端。

**(f) 结果（V16，C. 开发记录 §16.5）。** SPB 5.2（$2\times10^{11}$ Pa）激波后密度 3966 → 3780.0 kg/m³（精确 3766.7，误差 5.3% → 0.35%），$L_1(u)$ 118 → 25.4；SPB 5.1 2533.7 → 2528.6（精确 2527.5）。质量分配的 5% 误差就是分配规则的误差。其他一维算例的指标只在第 3 位以后变化。两相等熵涡 $64^2$：Newton 2.0 次/级、$\max|G|$ 到 $6\times10^{-16}$（V15：2.99 次、10503 次分支翻转、$10^{-13}$）。

**(g) 试过而舍弃（V16）。**
- 局部弱规则 $\alpha_k\mathcal G_k/K_k^2$（$K_k$ 取当地值）用于全部 $\delta\ge0$：5.2 修好，但 4.3.1 过冲 1.048、4.5 空气 $3\times10^5$ K（(d)）。
- 弱规则冻结在波前（非局部）：剖面 ODE 里微量气体 $10^{10}$ Pa 时 $p_1$ 偏 $+38\%$——强激波下气体份额必须随压力降低。
- 体积变化份额 $\lambda_k=Y_k\,\mathrm dv_k/\mathrm dv\propto\alpha_k/K_k$（把数值粘性应力作用在各相自己的压缩上，$\mathrm de_k=-P\,\mathrm dv_k$）：只在两相同步压缩时等于 (II.5)，弱极限就是错的（熵产 $\propto\mathcal G/K^2$，不是 $1/K$）；ODE 里水–空气 $10^{10}$ 时 $-4.4\%$。
- 按 $T_k\Delta s_k$（$T_k$ 取 $p_0$、$\bar p$ 或 $p_1$ 处的等熵温度）：气–液强激波下份额错一个量级（$10^{10}$ 时 0.95–0.99 对 $c=0.145$）。
- 质量分配（V10–V15）：由格式的公共数值运动粘性推出（$\nu=a\Delta x/2$，粘性加热 $\propto m_k$），描述的是格式不是模型。
- $D_0$ 符号分支 + active set（V15）：解决的是错误分支变量造成的不连续。
- 灵感来源：Saurel et al. 2007 的 (II.5)；弱激波熵产 $\Delta s=\mathcal Gv\Delta p^3/(6TK^2)$（Landau–Lifshitz §86）；Abgrall–Karni 2010 关于路径守恒格式极限解由数值粘性决定的分析。

### B5.3 闭合的输出

闭合给出 $\alpha_k,m_k,\varepsilon_k,p_k,(T_k)$，以及混合压力 $\bar p^{\rm cl}=\sum_k\alpha_kp_k$。全局方程（B8）是

$$
G_i(P)\equiv\bar p^{\rm cl}_i\big(s_i,X^{\rm tr}_i,(\rho e)^C_i;\ \mu,\theta,\nu\big)-P_i=0 .\tag{B5.4}
$$

含义：格子的压力等于它的物质对本步体积应变的热力学响应。这是 $p^{n+1}=p^*-\tau\rho c^2\partial_xu^{n+1}$ 的非线性、带松弛的推广。

---

## B6. 闭合的各个实例（按速率）

假设与 Pelanti (2022) 相同的时间尺度顺序：$\theta>0$ 或 $\nu>0$ 时 $\mu=\infty$；$\nu>0$ 时 $\theta=\infty$。

### B6.1 $\mu=\infty$，$\theta=\nu=0$：$p$ 闭合

$p_1=p_2=p$，$p_I=p$。把 $\omega_ks+\chi_kR=\alpha_k-\alpha^{\rm tr}_k$ 代入 (B5.1)：

$$
\boxed{\varepsilon_k=\varepsilon^{\rm tr}_k-p\,(\alpha_k-\alpha^{\rm tr}_k)+\lambda_k\delta},\qquad \delta=D_0+p\,s,\quad D_0=(\rho e)^C-\sum_j\varepsilon^{\rm tr}_j .\tag{B6.1}
$$

这就是离散 Gibbs 关系，与 $\omega_k$、$p_I$ 都无关。用 (A2.2)（性质 E2）：

$$
\alpha_k(p)=\frac{(\gamma_k-1)\big[A_k+B_kp\big]+(p+\gamma_k\pi_k)m_kb_k}{\gamma_k(p+\pi_k)},\qquad A_k=\varepsilon^{\rm tr}_k-m_kq_k+\lambda_kD_0,\quad B_k=\alpha^{\rm tr}_k+\lambda_ks .\tag{B6.2}
$$

局部方程是 $\Phi(p)=\sum_k\alpha_k(p)-1=0$。

- **单调性：** $\alpha_k'(p)=-\dfrac{(\gamma_k-1)N_k}{\gamma_k(p+\pi_k)^2}$，$N_k=A_k-(B_k-m_kb_k)\pi_k$。若 $\varepsilon^{\rm tr}_k$ 是某压力 $p^{\rm tr}_k>-\pi_k$ 下的状态，则 $\varepsilon^{\rm tr}_k-m_kq_k-(\alpha^{\rm tr}_k-m_kb_k)\pi_k=(\alpha^{\rm tr}_k-m_kb_k)\frac{p^{\rm tr}_k+\pi_k}{\gamma_k-1}>0$；可容许条件 $N_k>0$ 就是这个量加上 $\lambda_k(D_0-s\pi_k)$ 后仍为正。
- **存在唯一：** $p\downarrow\max_k(-\pi_k)$ 时 $\Phi\to+\infty$；$p\to\infty$ 时 $\Phi\to\sum_k\frac{(\gamma_k-1)B_k+m_kb_k}{\gamma_k}-1<\sum_kB_k-1=0$（因为 $m_kb_k<B_k$）。所以根存在且唯一，用带区间的 Newton 求解。
- 输出 $\bar p^{\rm cl}=p$。

**MG 相（V12，A2.1）。** (B6.1) 不变，只是 $\varepsilon_k$ 换成 (A2.5)。对 MG 相，$\alpha_k(p)$ 是

$$
f_k(\alpha)=m_k\varepsilon_r(m_k/\alpha)+\alpha\frac{p-p_r(m_k/\alpha)}{\Gamma_k}+p\,\alpha-A_k-B_kp=0,\qquad A_k=\varepsilon^{\rm tr}_k+\lambda_kD_0
$$

的根（$B_k$ 同 (B6.2)）。用 $\varepsilon_r'=p_r/\rho^2$、$d\rho/d\alpha=-\rho/\alpha$：

$$
f_k'(\alpha)=\frac{(\Gamma_k+1)(p-p_r)+\rho_kp_r'}{\Gamma_k}=\frac{\rho_kc_k^2}{\Gamma_k}>0 ,
$$

所以只要 $c_k^2>0$，$f_k$ 严格递增，根唯一，带区间的 Newton（在 $\ln\alpha$ 上）求解。外层 $\Phi(p)=\sum\alpha_k(p)-1$：

$$
\alpha_k'(p)=-\frac{\partial f_k/\partial p}{f_k'}=-\frac{\Gamma_k\big[\alpha_k(1+1/\Gamma_k)-B_k\big]}{\rho_kc_k^2},
$$

$s=0$、$\delta=0$ 时 $B_k=\alpha^{\rm tr}_k$，在 $\alpha_k=\alpha^{\rm tr}_k$ 附近为 $-\alpha_k/(\rho_kc_k^2)<0$（等熵压缩率）。$\Phi$ 一般不再有 SG 那样的全局单调证明，所以外层用带区间的 Newton，区间从初猜向两侧几何扩张；某个 $p$ 下 $\alpha_k(p)$ 无解（$c_k^2\le0$）当作 $\Phi=+\infty$（$p$ 太低）。B8 的残量尺度里 $K$ 用 Wood 体积模量（同前），$P_E$ 的分子把 $|m_kq_k|$ 换成 $|m_k\varepsilon_r|+|\alpha_kp_r/\Gamma_k|$（$\varepsilon_k$ 舍入的来源）。

### B6.2 有限 $\mu$，$\theta=\nu=0$：两压力闭合

记 $a_k=\alpha^{\rm tr}_k+\omega_ks=\alpha^{\rm tr}_k/(1-s)$，则 $\alpha_1=a_1+R$，$\alpha_2=a_2-R$；$p_I=w_1p_1+w_2p_2$，$w_1=Z_2/(Z_1+Z_2)$（$Z_k=\rho_kc_k$ 冻结在显式预测态 $X^*$ 的单位体积内容上，Newton 迭代中不变）。给定 $R$ 时，(B5.1)–(B5.3) 对 $(p_1,p_2)$ 是 2×2 线性方程组（$k=1,2$）：

$$
\Big(\frac{\alpha_k-m_kb_k}{\gamma_k-1}+\omega_ks\Big)p_k+\chi_kR\,(w_1p_1+w_2p_2)-\lambda_ks\,(\omega_1p_1+\omega_2p_2)
=\varepsilon^{\rm tr}_k-m_kq_k-\frac{\gamma_k\pi_k(\alpha_k-m_kb_k)}{\gamma_k-1}+\lambda_kD_0 .\tag{B6.3}
$$

再解标量方程

$$
\mathcal F(R)=R-\tau\mu\big(p_1(R)-p_2(R)\big)=0 .
$$

$R$ 增大时 $\alpha_1$ 变大、$p_1$ 下降、$p_2$ 上升，所以 $\mathcal F$ 在可容许的连通区间上递增，根唯一。输出 $\bar p^{\rm cl}=\alpha_1p_1+\alpha_2p_2$。

- **根的区间。** 不能取 $\big(-(a_1-m_1b_1),\,a_2-m_2b_2\big)$：被压缩相的功按终态 $p_I$ 计（后向 Euler），2×2 方程组在该相体积到零之前就奇异（SG、$p_I\approx p_2$ 时在 $R=a_2/\gamma_2$），$\mathcal F$ 在区间端点不定号。松弛使 $|p_1-p_2|$ 减小，所以根在 $0$ 与 $R_1=\tau\mu(p_1-p_2)|_{R=0}$ 之间（$\mathcal F(0)$ 与 $\mathcal F(R_1)$ 异号）。$R_1$ 不可容许时，奇点处 $\mathcal F$ 趋于被压缩相一侧的无穷大（符号正确），向 $0$ 二分即可得到可容许的区间。

- $\mu=0$：$R=0$，每相按自己的 $p_k$ 做 Gibbs 功，$\varepsilon_k-\varepsilon^{\rm tr}_k=-p_k(\alpha_k-\alpha^{\rm tr}_k)+\lambda_k\delta$，体积按 $\omega_k$ 各向同性压缩：**冻结格式**（与 Schropff 等 2026 的冻结 IMEX 同类）。

### B6.3 $\mu=\theta=\infty$，$\nu=0$：$pT$ 闭合

$m_k=m^{\rm tr}_k$。混合能量 $\sum_km_ke_k(p,T)=(\rho e)^C$ 给出（E3）

$$
T(p)=\frac{(\rho e)^C-\sum_km_kq_k}{\sum_km_kc_{vk}\frac{p+\gamma_k\pi_k}{p+\pi_k}},\qquad
\Phi(p)=\sum_km_kv_k\big(p,T(p)\big)-1=0 .\tag{B6.4}
$$

- 与 $\alpha^{\rm tr},s,\lambda$ 无关：$\alpha$ 是从属变量，这正是 4 方程 $pT$ 模型。
- 根唯一（P-S (C.3) 的可容许根）。
- 输出 $\alpha_k=m_kv_k$，$\varepsilon_k=m_ke_k$，$\bar p^{\rm cl}=p$。

### B6.4 $\mu=\infty$，有限 $\theta$，$\nu=0$

在 (B6.1) 上加传热 $\varepsilon_1\mathrel{+}=\Delta Q$，$\varepsilon_2\mathrel{-}=\Delta Q$，$\Delta Q=\tau\theta(T_2-T_1)$。

由 (B6.2)，$\alpha_k$ 对 $\Delta Q$ 线性：$\alpha_k=\alpha_k^0(p)+\chi_k\frac{(\gamma_k-1)\Delta Q}{\gamma_k(p+\pi_k)}$。于是

$$
T_k=\frac{(\alpha_k/m_k-b_k)(p+\pi_k)}{(\gamma_k-1)c_{vk}}=T^0_k(p)+\chi_k\frac{\Delta Q}{C_{pk}},\qquad C_{pk}=\gamma_kc_{vk}m_k ,
$$

即等压下吸热 = 焓增。所以

$$
\boxed{\Delta Q(p)=\frac{\tau\theta\,\big(T^0_2(p)-T^0_1(p)\big)}{1+\tau\theta\,(C_{p1}^{-1}+C_{p2}^{-1})}}\tag{B6.5}
$$

是闭式，局部方程仍是 $\Phi(p)=\sum_k\alpha_k-1=0$。$\theta=0$ 时回到 B6.1。

### B6.5 $\mu=\theta=\infty$，有限 $\nu$

$\rho=\sum m^{\rm tr}_k$ 不变。给定 $p$ 与 $Y=Y_2$：

$$
T(p,Y)=\frac{(\rho e)^C/\rho-(1-Y)q_1-Yq_2}{(1-Y)c_{v1}\phi_1+Yc_{v2}\phi_2},\quad \phi_k=\frac{p+\gamma_k\pi_k}{p+\pi_k};\qquad
\mathcal R(Y)=\rho Y-m^{\rm tr}_2+\tau\nu\,(g_2-g_1)\big(p,T(p,Y)\big)=0 .\tag{B6.6}
$$

- **单调：** 能量固定时 $Y\uparrow$ 使 $T\downarrow$（$e_2(p,T)>e_1(p,T)$）；$\partial_T(g_2-g_1)=-(s_2-s_1)<0$。所以 $\mathcal R'(Y)=\rho+\tau\nu(s_2-s_1)|\partial_YT|>0$，$Y$ 在可容许区间上唯一。
- **可容许区间。** 能量固定时 $T(p,Y)>0\iff e-q_1-Y(q_2-q_1)>0$，$q_2>q_1$ 时 $Y<(e-q_1)/(q_2-q_1)$（随比内能 $e$ 而变，例如 0D 算例 41 为 0.507）。所以区间是 $[Y_{\min},\min(1-Y_{\min},(e-q_1)/(q_2-q_1))]$，而不是 $[Y_{\min},1-Y_{\min}]$。$T\to0^+$ 时 $g_k\to q_k+b_kp$，$\mathcal R$ 在该端为正（$q_2>q_1$），与单调性一致。
- 若 $\mathcal R(Y_{\min})>0$ 或 $\mathcal R(1-Y_{\min})<0$，$Y$ 取界值：完全凝结/蒸发。这是速率方程的解落在可容许域边界，不是额外规则。
- 然后解 $\Phi(p)=\rho\sum_kY_kv_k(p,T)-1=0$。
- 能量方程只用混合量，$\rho_I,e_I$ 不出现。

### B6.6 $\mu=\theta=\nu=\infty$：$pTg$ 闭合

$T=T_{\rm sat}(p)$，$Y$ 由能量定（A7.2），$\Phi(p)=\rho\sum_kY_kv_k-1$。$Y$ 出界时钳住，$T$ 由能量定，与 B6.5 的界值分支相同，所以 $\Phi$ 连续。$p\ge p_c$（SG 饱和线的端点，A2）时没有两相平衡，按 $1/\rho$ 与 $v_c=v_1(p_c,T_{\rm sat}(p_c))$ 的大小取单相（$Y=Y_{\min}$ 或 $1-Y_{\min}$），$T$ 由能量定。

### B6.7 激活（按论文）与投影

$\theta,\nu$ 只在 $\epsilon_I\le\alpha_1\le1-\epsilon_I$ 的格子非零；$\nu$ 还要 $T>T_{\rm sat}(p)$（P-S (10)、Pelanti 8.2.1）。

**激活在每步开始时按 $U^n$ 判断一次，并把 $U^n$ 投影到该闭合的平衡流形上；本步各级冻结闭合类型。** 理由：

- 瞬时松弛（速率 $=\infty$）在激活时是状态的跳跃。若在级内激活，第一级把状态投影到流形上，这个变化是 $O(1)$、与 $\Delta t$ 无关，于是 $K_1=(U_1-U^n)/\tau=O(1/\Delta t)$。
- IMEX-RK 的级组合会外推这个跳跃。LSDIRK2 第二级的基底是 $Q^s_2=U^n+\frac{1-\gamma}{\gamma}(U_1-U^n)$，系数 2.41：微量相被过量移走，含量变负。数值上观察到：十二烷 $pT$ 中 $10^6$ K 的微量蒸汽热平衡后能量从 39 降到 0.02 J/m³，外推得 $-55$；$p$-$pTg$ 中微量液体完全蒸发，外推得 $\alpha_1=-1.5\times10^{-4}$。两者都使闭合无可容许根，$\Delta t$ 减半也无效。
- 刚性精确 IMEX-RK 的 AP 性质（B7）本来就要求初值在流形上（well-prepared）。投影后 $U_1-U^n=O(\Delta t)$，级组合不再外推跳跃。
- 投影 = $s=0$ 的同一个闭合：守恒 $m_k$（$pT$）或 $\rho$（$pTg$）、$\rho u$、$E$。按 $U^n$ 判断激活就是 P-S 的算子顺序（双曲步之后判断、松弛），只是松弛放在下一步开头。输出前同样投影。
- 有限速率的类型不投影：速率有限，状态以有限速度松弛，由各级闭合处理。
- 投影（$s=0$ 的闭合）无解的格子，本步退回力学闭合（不做热/化学松弛）。

## B7. 两个极限（创新点 5）

**定理 1（速率 $\to\infty$，AP）。** 设 B6.1 / B6.3 / B6.6 的局部方程在其解处 Jacobian 非奇异（B6 的单调性保证）。则：

- **$\mu$：** 令 $\kappa=1/(\tau\mu)$，(B6.3) 等价于 $p_1-p_2=\kappa R$。$\kappa=0$ 时代回 (B5.1)、(B5.2) 就是 (B6.1)，因为 $\omega_ks\pm R=\alpha_k-\alpha^{\rm tr}_k$，$p_I=p$，$\delta=D_0+ps$。由隐函数定理，局部解对 $\kappa$ 光滑，$|\text{解}(\kappa)-\text{解}(0)|=O(\kappa)$，常数与 $\Delta x$ 无关（逐格）。再由全局 Jacobian 非奇异（B8），全局离散解收敛到瞬时格式的离散解，误差 $O(1/(\tau\mu))$。
- **$\theta$：** (B6.5) 给出 $T_1-T_2=O(1/(\tau\theta))$，极限方程组就是 (B6.4)；由 $pT$ 根的唯一性，极限就是 $pT$ 闭合。
- **$\nu$：** (B6.6) 给出 $|g_2-g_1|=|\rho Y-m^{\rm tr}_2|/(\tau\nu)=O(1/(\tau\nu))$（未钳位时），极限就是 $pTg$ 闭合。

**极限格式是平衡模型的相容离散。** 以 $\mu$ 为例：把 (B6.1) 线性化，每相等熵且等压，$d\alpha_k=s\,\frac{\alpha_k/(\rho_kc_k^2)}{\sum_j\alpha_j/(\rho_jc_j^2)}$。代入 $s=\Delta t\,\partial_xu$：

$$
\partial_t\alpha_1+u\partial_x\alpha_1=\alpha_1\Big(\frac{\rho c_p^2}{\rho_1c_1^2}-1\Big)\partial_xu=\frac{\alpha_1\alpha_2(\rho_2c_2^2-\rho_1c_1^2)}{\alpha_1\rho_2c_2^2+\alpha_2\rho_1c_1^2}\partial_xu ,
$$

这正是 Kapila 模型的体积方程 (Pelanti (9a))。

**时间上：** 隐式表是刚性精确、L 稳定的 SDIRK（B10），每一级都在平衡流形上，$U^{n+1}=U_s$。所以 $\mu\to\infty$ 与 $\Delta t\to0$ 可交换，格式是 AP 的。

**定理 2（速率 $\to0$，冻结）。** $\mu=0$ 时 $R=0$（B6.2）；$\theta=0$ 时 $\Delta Q=0$；$\nu=0$ 时 $Y=Y^{\rm tr}$。闭合退化为：

- 无松弛的 6 方程 TV-IMEX：每相 $\varepsilon_k-\varepsilon^{\rm tr}_k=-p_k(\alpha_k-\alpha^{\rm tr}_k)+\lambda_k\delta$，体积 $\alpha_k=\alpha^{\rm tr}_k/(1-s)$；
- 或依次退化为 $p$ 闭合、$pT$ 闭合。

这些退化是**逐位**的：代码里速率为 0 与“关掉松弛”走同一组公式。

**定理 3（低马赫极限，V11）。** 混合物的低马赫极限（BB §3.6 的意义：一个混合物马赫数 $M$）。无量纲化：压力以 $\rho_rc_r^2$ 为尺度，$\Delta t,\Delta x$ 固定，well-prepared 初值 $P=P_0+M^2P_2$，$P_0$ 为常数。

- **压力方程。** 隐式级的方程 $G_i(P)=\bar p^{\rm cl}_i-P_i=0$ 线性化后是 $P_i-P^*_i=-\tau K_i(\delta\cdot u_f)_i+\dots$，$K=\rho c_W^2=O(1)$（无量纲），而面速度里 $P$ 的系数是 $\tau/(M^2\bar\rho\Delta x)$。代入 $P=P_0+M^2P_2$：$O(M^{-2})$ 与 $O(1)$ 两阶给出 $(\delta\cdot u_f)_i=O(M^2)$ 与
  $$\delta\cdot\Big(\bar u^*_f-\frac{\tau}{\bar\rho_f}\frac{\delta P_2}{\Delta x}\Big)=O(M^2) ,$$
  即面速度的离散投影（紧致 Laplacian 的 Poisson 方程）。
- **动量。** $(\rho\mathbf u)^{n+1}=(\rho\mathbf u)^*-\tau\,\delta p_f/\Delta x$，$p_f$ 是 $P$ 的密度加权平均（$P_0$ 抵消），剩下 $\tau\,\delta(\overline{P_2})/\Delta x$：格速度用宽模板梯度修正（近似投影）。
- **补足项。** $c\tau/\Delta x=a_{ii}\mathrm{CFL}\,c/a_\Omega=O(1/M)$，$M<a_{ii}\mathrm{CFL}$ 后 $\varphi=0$ 恒成立，补足项不出现在极限格式里。
- **显式部分。** 只有物质速度 $a_\Omega$ 的 Rusanov 耗散（B3），与 $M$ 无关。
- **松弛。** 闭合在每一级给出平衡流形上的状态（定理 1），与 $M$ 无关。

所以极限格式是混合物不可压 Euler 方程的相容离散（显式的对流 + 面速度投影），格式在低马赫极限下是 AP 的。

**一致性（不只是极限）。** AP 只说 $M\to0$ 的极限格式相容，还要求中间马赫数与起步斜坡里没有 $O(1/M)$ 的耗散（Dellacherie 2010，BB §3.6 末）。逐项：

| 项 | 作用在 | 量级 |
|---|---|---|
| 显式 Rusanov | 所有量 | $a_\Omega\Delta x$（物质速度） |
| 隐式步本身 | $P$ | $\tau c^2$：L 稳定地阻尼声波 |
| 补足项的压力扩散 | $u_f$ | 相对 $O(M)$（B4.1(c)） |
| 补足项的速度项，逐维形式（V10 及以前） | 动量 | 相对 $O(1/M)$，$\varphi>0$ 时 |
| 补足项的速度项，散度形式（V11） | 动量 | 相对 $O(M)+O((\Delta x/L)^3/M)$ |

V10 的逐维速度项在 $a_{ii}\mathrm{CFL}<M<1$ 和起步斜坡（$\varphi\approx0.9$，约 $\log(1/M)/\log1.2$ 步）里是声速量级的粘性。Gresho 涡 $M=0.1$ 时动能保持率 0.36。$M=0.01$ 时为 0.50：补足项只在斜坡里作用，只占其中不到 0.1%（速度项乘 $\theta$ 时为 0.503）；其余损失来自共用限制因子（B3）。V11 改成散度形式、速度单独重构后，三个马赫数都是 0.93–0.94（B4.1(c) 的表）。

**直接检验（固定 $\Delta t$）。** AP 要求的是 $\Delta t$ 与 $M$ 无关时一致稳定、并收敛到极限。Gresho 涡用同一个 $\Delta t=0.0025$（503 步，无起步斜坡），$M=10^{-1}\ldots10^{-6}$（算例 50–55，`tests/run_gresho_ap.sh`）：
- 每个 $M$ 都是 0 次拒绝，每级约 1 次 Newton；
- 动能保持率 $M\le0.01$ 时到四位相同（0.9372），$M=0.1$ 为 0.930；
- $|\rho-1|\le0.11M^2$；$p_2$ 与不可压精确解之差 $1.4\times10^{-2}$（离散误差），与 $M$ 无关；
- 速度的 Cauchy 差（$M$ 与 $M/10$）为 $1.8\times10^{-3}$、$1.9\times10^{-6}$、$1.6\times10^{-7}$、$2.0\times10^{-9}$（$M=0.1\ldots10^{-4}$）。

做到这一点，求解器有三处必须与 $M$ 无关（B8）：未知量取压力扰动 $P-P_{\rm ref}$（否则舍入噪声 $\nu^2\varepsilon_{\rm mach}p_0$，以 $\rho U^2$ 为单位是 $O(\varepsilon_{\rm mach}/M^4)$）；差分步长与步长判据按 $P$ 的自然尺度 $s=N/D$；线性解法用两级预条件。离散格式本身不变。$M\le10^{-5}$ 时剩下的是双精度的底：$\rho U^2/p_0=\gamma M^2$，相对绝对压力的 Newton 容差以 $\rho U^2$ 为单位是 $2.4\times10^{-13}/M^2$，$M=10^{-6}$ 时 $p_2$ 有约 $10^{-2}$ 的噪声。详见 C. 开发记录 V11 §11.8，图 `docs/figures/gresho_ap.png`。

---

## B8. 全局非线性方程组与求解

$G_i(P)=\bar p^{\rm cl}_i-P_i$ 只依赖 $(P_{i-1},P_i,P_{i+1})$（二维 5 点），Jacobian 三对角（二维 5 对角）。

**单调性：**

- $\partial s_i/\partial P_i=r^2(1/\bar\rho_{i+1/2}+1/\bar\rho_{i-1/2})+r\sum_f\varphi_f\psi^P_f/(Z_i+Z_{i\pm1})>0$（B4.1 的补足项只增加对角占优），且 $\partial\bar p^{\rm cl}/\partial s<0$：膨胀降压（E5）。
- $\partial E_i/\partial P_i\approx-r^2\,2\bar p/\bar\rho<0$。

所以 $\partial G_i/\partial P_i<-1$，相邻项 $O(r^2\rho c^2/\bar\rho)$，矩阵是 $I-\tau^2\partial_x(c^2\partial_x)$ 型：对角占优的 M 矩阵。

**求解：**

- Newton。Jacobian 用逐格差分：扰动三个（五个）压力，各算一次局部闭合。
- **$P$ 的自然尺度（V11）。** $s_i=N_i/D_i$。$N_i=|P_i|+K_i+P_{E,i}$ 是残量尺度（见下面的收敛判据）；$D_i=1+(K_i+P_{E,i})\sum_f(\tau/h_f)k_f$ 是 $\partial G_i/\partial P_i$ 的尺度（$k_f$ 是面速度里压力差的系数 $\tau/(\bar\rho_fh_f)+\varphi\psi/(Z_L+Z_R)$，故 $D\approx1+2d\,\nu^2$，$d$ 为维数）。$s$ 是残量所能分辨的压力变化：非刚性时 $s\approx N$，低马赫时 $s\approx N/(2d\,\nu^2)$。
- Jacobian 的差分步长 $h_i=10^{-7}s_i$（V11）；Newton 至多 60 次。
  - 原来 $h=10^{-7}(|P|+1)$。扰动 $h$ 在本格引起的体积变化约 $h\sum_f(\tau/h_f)k_f\approx10^{-7}\nu^2/\gamma$（$|P|\approx K/\gamma$ 时），低马赫时 $\nu=a_{ii}\mathrm{CFL}/M$，$M=10^{-5}$ 时已是 $O(1)$：差分 Jacobian 失效，Newton 失败，Gresho $M=10^{-5}$ 每步被拒绝两次（4019 次拒绝）。
  - 按 $N/D$ 取步长：体积变化不超过 $10^{-7}N/K$，舍入引起的导数相对误差约 $\varepsilon_{\rm mach}N/(hD)=\varepsilon_{\rm mach}/10^{-7}$，两者都与 $M$ 无关。这是差分 Jacobian 按变量的典型尺度取步长的通常做法（Dennis–Schnabel 的 typical size），这里典型尺度是“残量尺度 / Jacobian 尺度”。
- 一维非周期用 Thomas；二维与周期边界用 BiCGStab（`2DEuler_TV` 风格，但用相对容差 $10^{-13}$），预条件是两级的（V11）：
  $$M^{-1}r=D^{-1}r+\mathbf 1\,\frac{\mathbf 1^{\rm T}r}{\mathbf 1^{\rm T}J\mathbf 1},$$
  即 Jacobi（用本格扰动的那个元素）加常数向量的粗空间校正（Nicolaides 1987 的两级预条件）。
  - 理由：面速度与动量只依赖压力差，所有边界条件下（周期、固壁、透射虚格 = 本格）$P$ 的常数平移只使 $G$ 变化 $O(1)$（只经能量通量与 $-P$ 项），而对角元约 $1+2d\,\nu^2$。Jacobi 预条件后常数模态的特征值约 $1/(2d\,\nu^2)$，其余模态的谱与 $M$ 无关。
  - $M=10^{-6}$ 时这个特征值约 $10^{-10}$，BiCGStab 发散（相对残量 $10^{129}$），37 次拒绝。粗空间校正把它移到 $[1,2]$，其余模态不变：预条件后的谱与 $M$ 无关。
  - 粗空间只加一个标量，各格加的是同一个数：y 向均匀的问题各行仍逐位相同。
- 边界虚格（透射、固壁）的压力等于本格压力。它在 Jacobian 里是单独的元素（列号为本格），不并进对角元：对角元只扰动本格，虚格的导数单独差分；每行按“本格、W、E、S、N”的固定次序存 5 个元素，同列的也不合并；Jacobi 预条件用本格的那个元素，Thomas 把同列元素相加。这样数学上 Jacobian 不变，而边界行与内部行的浮点运算相同：y 向均匀的问题各行逐位相同（V8；并进对角元时边界行的舍入不同，再经共用限制因子的开关放大成截断误差大小的行间差）。
- **能量亏损的分配（V16）。** B5.2 的 $\Lambda_k(\delta)$ 是 $\delta$ 的连续分段仿射函数，分段按本格当前迭代值处的 $\delta=D_0+Ps$ 取（收敛时与解处一致），所以 $G(P)$ 连续、分段光滑（$\delta=0$ 与 $\delta=Q$ 处导数跳变），差分 Jacobian 直接用。V15 的分支冻结与 active-set 校正（分支按 $D_0$ 的符号取，$\lambda_k$ 又乘 $s$，$G$ 在 $D_0=0$ 处跳 $(\lambda^E_k-Y_k)ps$）因此去掉：$\lambda_k$ 只通过 $\lambda_k\delta$ 进入 (B6.2)，按 $\delta$ 分支时跳点处为 0。两相等熵涡 $64^2$：Newton 2.0 次/级、$\max|G|$ 到 $6\times10^{-16}$（V15：2.99 次、10503 次翻转、$10^{-13}$）。
  - 两相等熵涡：MOOD 回退 145 → 0；空化 RM：拒绝 41 → 0、MOOD 162 → 25；一维算例的指标在第 5 位以后变化。代价是每级 Newton 次数约 +1。
- 线搜索（步长减半，至多 40 次）：接受条件是所有格的局部闭合有解（可容许），且残量范数满足 Armijo 下降 $\|G\|_{\rm new}<(1-10^{-4}\lambda)\|G\|$。
- 收敛判据：$\max_i|G_i|/(|P_i|+K_i+P_{E,i})<10^{-13}$，$K=\rho c_W^2$（Wood 体积模量，$G/K$ 就是体积约束残量 $\sum\alpha-1$），$P_E=\sum_k(|\varepsilon_k|+|m_kq_k|)/\sum_k\partial\varepsilon_k/\partial p$（$X^*$ 处；代码里分母再加 1 Pa 防零）。理由：闭合给出的压力只精确到 $\varepsilon_{\rm mach}(p+\pi)$，守恒总能的舍入 $\varepsilon_{\rm mach}|E|$ 又经 $\lambda_k\delta$ 进入各相压力，所以 $G$ 的可达精度是 $\varepsilon_{\rm mach}\max(K,P_E)$。按 $|P|$ 归一时 $p\ll\pi$ 的格子永远达不到容差：环氧–尖晶石（$\pi=1.41\times10^{11}$）算例因此每步都被拒绝（3659 次拒绝、13906 步），改后 0 次拒绝、122 步，每级 3 次 Newton。
- **扰动未知量（V11）。** Newton 的未知量是 $P-P_{\rm ref}$，$P_{\rm ref}$ 是初始混合压力的最小值（常数）。
  - $p_f=w_LP_L+w_RP_R-p_v$ 的权和为 1，所以面压力、面速度（$P_R-P_L$）、动量（$p_{f,e}-p_{f,w}$）都只含扰动量；只有能量通量 $p_fu_f$ 与闭合需要绝对压力，在那里加回 $P_{\rm ref}$：$(P_{\rm ref}+p_f')u_f$，闭合输入与残量 $G=(\bar p^{\rm cl}-P_{\rm ref})-P'$。精确算术下与原格式完全相同。
  - 浮点下：存下来的 $P$ 带 $\varepsilon_{\rm mach}|P|$ 的舍入，它经面速度、体积通量与闭合放大 $\nu^2$ 倍（见下一条）进入残量与解。低马赫时 $|P|\approx p_0=1/(\gamma M^2)$、$\nu\sim1/M$，这个噪声是 $\nu^2\varepsilon_{\rm mach}p_0\sim\varepsilon_{\rm mach}/M^4$（相对 $\rho U^2$）。Gresho $M=10^{-4}$ 时约 $3\times10^{-3}$：压力的 Cauchy 差（相邻两个马赫数之差）从 $M=10^{-3}$ 的 $3.9\times10^{-4}$ 反而升到 $2.8\times10^{-3}$，$\max|\rho-1|/M^2$ 从 0.10 升到 0.24。
  - 用扰动量后放大的只是 $\varepsilon_{\rm mach}|P-P_{\rm ref}|\sim\varepsilon_{\rm mach}\rho U^2$，噪声变为 $\nu^2\varepsilon_{\rm mach}\rho U^2\sim\varepsilon_{\rm mach}/M^2$。剩下的舍入底是总能形式本身的 $\varepsilon_{\rm mach}p_0$（$p=(\gamma-1)(E-\tfrac12\rho u^2)$），它不再被放大。
  - $P_{\rm ref}$ 取最小值而不是平均值：$|P-P_{\rm ref}|$ 不超过当地压力，多相问题（压力跨几个量级）里残量的舍入底不会因参考值变大。
  - 这是低马赫可压缩求解器的标准做法：Klein 1995 的多尺度展开 $p=p_0+M^2p_2$；Sesterhenn–Müller–Thomann 1999 以压力扰动为未知量。
- **步长判据（V11）。** 完整 Newton 步满足
  $$\max_i\frac{|\delta P_i|}{|P_i-P_{\rm ref}|+s_i}<10^{-13}$$
  时接受 $P+\delta P$：Newton 步相对未知量本身的大小与残量的分辨率都已在容差内，$P$ 的误差已在容差内。这是有限精度下 Newton 法的标准判据（步长与残量两个判据）。
  - **来由。** $P$ 的舍入 $\varepsilon_{\rm mach}|P|$ 经面速度（$\tau/(\rho\Delta x)$）、体积通量（$\tau/\Delta x$）与闭合（$K=\rho c^2$）进入残量，放大 $K\tau^2/(\rho\Delta x^2)=\nu^2$ 倍（$\nu=c\tau/\Delta x$）：残量有一个约 $\nu^2\varepsilon_{\rm mach}|P|$ 的舍入底。低马赫时 $\nu=a_{ii}\mathrm{CFL}/M$ 很大，Gresho $M=0.001$（$\nu\approx60$，未知量还是绝对压力时）这个底约 $|G|\approx6\times10^{-7}$ Pa（$P=7.1\times10^5$），归一后 $(1\text{–}2)\times10^{-13}$，高于容差。Newton 一步就把残量从 $10^{-5}$ 降到这个底，线搜索再也降不下去，只有残量判据时判为失败：拒绝、减半、斜坡重启，1495 次拒绝，7205 个一阶步。加步长判据后 0 次拒绝、911 步。
  - **尺度。** 最初用与残量相同的尺度 $|P_i|+K_i+P_{E,i}$。低马赫时它比 Newton 步能改变的量 $s$ 大约 $\nu^2$ 倍：$M=10^{-6}$ 时第一步（$|\delta P|\approx3\times10^{-4}$ Pa，$N\approx2.4\times10^{12}$）就被接受，而一步之后残量的二次项还有 $\nu^4|\delta P|^2/K\approx10$ Pa（归一后 $6\times10^{-12}$，高于容差）：判据接受了一个没有收敛的迭代。改用 $|P-P_{\rm ref}|+s$：舍入底的情形里 Newton 步就是未知量的舍入 $\varepsilon_{\rm mach}|P-P_{\rm ref}|$ 或残量舍入折回的 $\varepsilon_{\rm mach}s$，判据照样成立；真正的修正（$M=10^{-6}$ 第一步的 $3\times10^{-4}$ Pa 对 $s\approx170$ Pa）不再被当成收敛，终残量不超过容差。非刚性时 $s\approx N$，与原判据相同。
  - 用扰动未知量之后，低马赫时舍入底变为 $\nu^2\varepsilon_{\rm mach}|P-P_{\rm ref}|$，不再随 $p_0$ 增大；判据仍保留，因为压力跨度大的问题（$|P-P_{\rm ref}|$ 大的格子）里同样的底仍可能高于容差。
  - 日志里的 `maxres` 是各级接受时的残量，可以略高于容差（它记的是舍入底）。

**步长控制（不是补丁，而是可容许性条件）：** 以下任一情况发生时，拒绝本步，$\Delta t\leftarrow\Delta t/2$ 重算（至多 30 次，否则报告失败）：
- $Q^e$ 出现 $\rho\le0$（显式通量只需要这一条，B3）；
- 闭合输入不可容许（$\alpha^{\rm tr}_k\le m^{\rm tr}_kb_k$、$m^{\rm tr}_k\le0$ 或 $N_k\le0$）；
- Newton 不收敛。

这和显式格式的正性 CFL 条件是同一类东西。所有拒绝都计数并报告。$\Delta t<10^{-12}t_{\rm end}$ 时报告失败并停止（V6：此前 $\Delta t$ 下溢为 0 时会无限循环）。

---

## B9. 性质证明

### B9.1 均匀 $(p,u)$ 保持（R2，创新点 1）

设 $u\equiv u_0$、$p_k\equiv p$，$\alpha_k,\rho_k$ 任意（各相温度也可以任意）。

1. **显式：** 由 (E1)，格均值满足 $\varepsilon_k=c_k(p)\alpha_k+q_km_k$、$\rho u=u_0\rho$；共用限制因子的重构保持这些线性关系，所以面值也满足。(B3.1) 对 $X$ 线性，所以 $\varepsilon^*_k=c_k\alpha^*_k+q_km^*_k$，$(\rho u)^*=u_0\rho^*$，$E^*=\sum\varepsilon^*_k+\tfrac12\rho^*u_0^2$。
2. **隐式：** 取 $P\equiv p$。$\bar u^*_f=u_0=\bar u^e_f$，(B4.1) 的两个压力项与速度跳跃项都为零，所以 $u_f=u_0$、$p_f=p$，$\Delta u_f=0$，(B4.2') 的增量为零，$s=0$，$X^{\rm tr}=X^*$。$(\rho u)^{n+1}=u_0\rho^*$，$E^{n+1}=E^*-r\delta[pu_0]=E^*$，所以 $(\rho e)^C=\sum\varepsilon^*_k$，$D_0=0$，$\delta=0$。
3. **闭合：**
   - $p$ 闭合：$p$ 处 $\alpha_k=\alpha^*_k$ 满足 (B6.1) 与 EOS，$\sum\alpha^*_k=1$，由根唯一，$\bar p^{\rm cl}=p$，$G=0$。
   - $\mu$ 有限：$R=0$，$p_1=p_2=p$ 满足 (B6.3)。
   - 热闭合：只在两相 $T$ 相等时才是平衡态，此时 $(p,T)$ 满足 (B6.4)，同理成立。
4. 所以 $P\equiv p$ 是全局解，$u,p$ 一步后仍均匀。证明只用了 (E1) 与 R3。

**二阶空间、二阶时间：** 共用限制因子的重构保持面值之间的线性关系（B3）。IMEX-RK 的 $Q^e,Q^s$ 是均匀 $(u,p)$ 态的线性组合，由 (E1) 仍是均匀 $(u,p)$ 态。所以每一级都满足上面的论证。

### B9.2 温度一致性（R3，创新点 4）

设各相 $(p,T_k)$ 都均匀（$\rho_k,e_k$ 均匀），$u$ 均匀。则 $m_k=\rho_k\alpha_k$、$\varepsilon_k=\rho_ke_k\alpha_k$ 都与 $\alpha_k$ 成正比。同一个线性算子输运后比例不变，所以 $\rho_k,e_k,T_k$ 不变。

若质量用 $u^e$、能量用 $u^{n+1}$（纯 TV），比例被破坏，$e_k$ 的误差为 $O(r|\Delta u|)$。对力学松弛，SG 的压力与 $\rho_k$ 无关，看不出来；对热/化学闭合，$T$ 直接决定饱和压力，就会出错。

### B9.3 守恒

$m_k$（含传质时 $\rho$）、$\rho u$、$E$ 都是通量形式，内部守恒到舍入误差。

### B9.4 混合能量一致性

$\sum_k\varepsilon_k=(\rho e)^C=E-\tfrac12\rho u^2$，且 $\bar p=\sum\alpha_kp_k$ 由相 EOS 给出。所以 P-S 定义 3.1 的 (i)(ii) 都满足（$\mu=\infty$ 时 $\bar p=p$ 就是混合 EOS (A2.3) 给出的压力）。

---

## B10. 时间积分（`2DEuler_TV` 的半隐式 IMEX-RK，Boscarino–Filbet–Russo）

`2DEuler_TV` 的两张表（`RKs=2`、`RKs=3`）都可用：

$$
\text{LSDIRK2(2,2,2)}:\ \tilde A=\begin{pmatrix}0&0\\ \hat c&0\end{pmatrix},\ A=\begin{pmatrix}\gamma&0\\ 1-\gamma&\gamma\end{pmatrix},\ \gamma=1-\tfrac1{\sqrt2},\ \hat c=\tfrac1{2\gamma};
$$

$$
\text{SSP2(3,3,2)}:\ \tilde A=\begin{pmatrix}0&0&0\\ \tfrac12&0&0\\ \tfrac12&\tfrac12&0\end{pmatrix},\ A=\begin{pmatrix}\tfrac14&0&0\\ 0&\tfrac14&0\\ \tfrac13&\tfrac13&\tfrac13\end{pmatrix}.
$$

两张表都满足：$b=A$ 的末行（刚性精确），$\tilde b=b$，$\sum b=1$，$b\cdot\tilde c=b\cdot c=\tfrac12$（二阶），$|R(-\infty)|=0$（L 稳定，已数值核对）。

区别在显式级态的外推：
- LSDIRK2 的 $\tilde c_2=1/(2\gamma)=1.707$，$Q^e-U^n=\frac1{2\gamma^2}(U_1-U^n)$，是第一级增量的 5.83 倍，在界面附近 $Q^e$ 容易出可容许域；
- SSP2(3,3,2) 的 $\tilde c\in[0,1]$，外推因子为 2。

默认 `rk=2`（LSDIRK2）；水–气激波管上 `rk=3` 与之结果相近（B13）。

**不可容许时的处理（V14）。** 第二级的显式预测 $X^*=Q^s-r\,\delta F(Q^e)$ 中 $Q^s\neq Q^e$（LSDIRK2 的 $Q^e$ 是外推），B3 的保正论证（$X^*$ 是非负面值的凸组合）不成立，界面附近微量相的 $\alpha^*$ 可以为负，任何压力下闭合都不可容许（十二烷，C. 开发记录 V14 §14.4）。所以一步不可容许时，先在**同一步长**上用一阶（单级后向 Euler，$Q^s=Q^e=U^n$，B3 保正）重算这一步；仍不可容许才减半步长（报告为拒绝）。这是 a posteriori 降阶（MOOD，Clain–Diot–Loubère 2011），只影响这一步，日志记为 `mood`。

**时间步：物质 CFL 按算例（V14，`TestSetup::cfl`，默认 0.2）。** IMEX 的步长按物质速度取，声波不被时间分辨。声波本身是要解的结构时（低马赫以外），这在声学 Courant 数接近 1 时给出 LSDIRK2 的色散误差：十二烷激波管里液体中 $10^8$ Pa 的稀疏波（$c\approx1400$ m/s）在物质 CFL 0.2 下声学 Courant 数约 0.9，平台压力误差 0.17；物质 CFL 0.05 时 0.054。所以十二烷（SG）与 SPB 5.1（激波与稀疏波都是声学结构，$u\ll c$）取 0.05，其余 0.2。低马赫算例照旧用大步长（AP 不受影响）。

**收敛阶试验（V14）。** 光滑解的阶用 `limiter=0`：体积输运量与速度都用不限制的中心斜率，补足项的 $\psi$ 用同一个中心重构（$\psi D=X^R-X^L$，截到 $[0,1]$）。minmod 在光滑极值处降为一阶，测出来的是限制器的阶。

**有限速率松弛的刚性中间区。** 对松弛 ODE $y'=\lambda y$（$z=\lambda\Delta t<0$），LSDIRK2 的稳定函数是 $R(z)=\frac{1+(1-2\gamma)z}{(1-\gamma z)^2}$：$z>-2.41$ 时 $R>0$（单调），$z<-2.41$ 时 $R<0$（越过平衡；最小值 $-0.207$，在 $z\approx-8.2$），$z\to-\infty$ 时 $R\to0$（L 稳定）。SSP2(3,3,2) 同类（最小值 $-0.212$，在 $z\approx-8.3$）。二阶方法不可能对所有 $z$ 都单调（Bolley–Crouzeix），后向 Euler（$R=1/(1-z)$）是单调的。

- 一般情况下越过平衡的量会在后续各步衰减。
- P-S/Pelanti 的化学松弛只在 $T>T_{\rm sat}$ 时激活（只蒸发）：越过饱和线后激活条件不再成立，越过的量被冻结。0D 算例 41：$\nu=1$（$z\approx-14$）时 $g_1-g_2$ 从 $7.9\times10^4$ 变到 $-1.3\times10^4$ 后不再变化。
- $\nu\to\infty$ 时偏差按 $O(1/\nu)$ 消失（AP 成立）。后向 Euler 在任何 $\nu$ 下都单调到达饱和。
- 所以：速率刚性但有限（$z\in(-10^3,-2.4)$）时，用瞬时闭合（AP 保证其为极限）或 `rk=1`。Pelanti 8.2.1 的 $\nu=10^{-4}$ 在物质 CFL 步长下 $|z|\sim10^{-3}$，不受影响。

第 $i$ 级：

1. $Q^e=U^n+\Delta t\sum_{j<i}\tilde a_{ij}K_j$，$Q^s=U^n+\Delta t\sum_{j<i}a_{ij}K_j$。
2. $U_i=\Phi(Q^e,Q^s;\tau=a_{ii}\Delta t)$：B3 的显式更新，加 B4–B8 的隐式闭合。
3. $K_i=(U_i-Q^s)/\tau$。

隐式表刚性精确，所以 $U^{n+1}$ 就是末级 $U_s$。$K_i$ 包含松弛源项的贡献，松弛是隐式表的一部分（BB §3.5 / Pelanti 的源项在这里由同一个闭合给出）。

**起步：** 从静止或间断初值起算时，首步用声学步长（`cflac`，默认 0.4），之后每步至多放大 1.2 倍，直到物质 CFL（`2DEuler_TV` 的 `full_step` 与 BB 的做法）；被 1.2 倍限制的步用一阶（单级后向 Euler）。`dt=` 给定固定步长（时间收敛阶试验用，不做斜坡）。斜坡里 $\varphi\approx0.9$，补足项的速度部分按散度作用（B4.1(c)），所以斜坡不在低马赫流里引入声速量级的粘性（定理 3）。

---

## B11. 二维

- 显式：$x,y$ 两个方向的面各算一次 (B3.1)。
- 面速度 (B4.1) 只修正法向分量；增量质量通量携带新的法向面速度 $u_f$，以及按供体取的切向速度 $\hat v$。
- 补足项的速度部分用顶点散度 $D^*_f$（B4.1(c)），用到面两侧行（列）的格速度，虚格按边界条件取值。
- 显式重构：体积输运量共用限制因子，$u,v$ 各自 minmod（B3）。
- $s_{ij}=r_x\delta_x[u_f]+r_y\delta_y[v_f]$。
- 动量：$x$ 方向 $-r_x\delta_x[p_f]-r_x\delta_x[u_f\hat M]-r_y\delta_y[\hat u\hat M]$，$y$ 方向对称。
- 总能的增量动能取 $\tfrac12(u_f^2+\hat v^2)\hat M$。
- 全局 Jacobian 5 对角。
- **二维检查（V8）：** 一维算例放在 $N\times4$ 网格上（y 向透射、均匀），各行应与一维结果相同、$v\equiv0$。12 个算例（SPB 4.1、4.3.1、4.3.2、4.5、5.1，P-S ±2、±500 m/s，十二烷 SG/NASG 的 $p$ 与 $pTg$，Pelanti $pTg(f)$，光滑解）各行逐位相同，$v\equiv0$。固定步长时第 0 行与一维结果之差：SPB 5.1 为 $4\times10^{-9}$，4.3.1 为 $10^{-6}$（一维用 Thomas、二维用 BiCGStab，舍入不同；欠分辨的 P-S 中心尖峰处经限制器开关放大到 $10^{-2}$，仍远小于离散误差）。显式程序同样各行逐位相同。真正二维的算例见下面（V13）。
- **边界条件（V13）。** 每一边单独给（`spb6_physics.h` 的 `bc_apply`，两个程序共用）：透射（虚格复制边界格）、固壁（镜像，法向速度取 $2u_w-u$，$u_w=0$ 时就是取负；总能随动能改）、周期。固壁面上隐式面速度 $u_f=u_w$（静止壁为 0，体积增量为 0；活塞 $u_w\neq0$ 时体积通量从虚格按 (B4.2') 迎风进入）。网格不动：活塞面固定在边界上，被推动的液体从该面进入（P-S 5.4 用移动边界，这里是固定网格上的近似，记录在 C 页 V13）。
- **初值：** 圆形区域的格按 $4\times4$ 个子点上守恒量的平均初始化（`TestSetup::sub`）。
- **并行：** 两个程序的逐格循环与 BiCGStab 用 OpenMP；BiCGStab 的内积与粗空间求和按固定大小的块求和，结果与线程数无关、逐位可重复；各格的算术与串行相同，一维输出逐字节不变。

---

## B12. 与 handoff v23 的差别（换路径的地方）

| | handoff v23 | 本推导 |
|---|---|---|
| 混合能量一致性 | 闭合用非守恒相能量，级末再用总能重初始化 | $\delta$ 进入闭合（B5.2），没有级末投影；动量方程看到的就是守恒能量下的压力 |
| 松弛 | 瞬时（$p$、$pT$、$pTg$）在闭合里；有限速率没有 | 一个闭合，$\mu,\theta,\nu$ 是参数（B6），两端极限严格（B7） |
| 全局残量 | $\sum\alpha-1$ 或能量残量，按闭合类型不同 | $G=\bar p^{\rm cl}-P$，所有闭合同一形式 |
| 传质判据 | 上一级的标记（滞后） | 每步开始按 $U^n$ 判断并投影到平衡流形，本步各级冻结（B6.7） |
| 显式耗散 | 中心部分也用重构值 | 同样中心部分也用重构值（BB 的格均值中心部分不是 TVD，B3）；耗散速度取物质速度，且不小于面态速度（B3 的保正条件）；体积输运量共用限制因子，速度单独重构（V11） |
| 声学耗散 | 无（只有后向 Euler） | 声学 CFL < 1 时补足到 Godunov（B4.1），否则不变；速度部分按散度作用，低马赫一致（B4.1(c)，定理 3） |
| 增量供体 | 按 $\Delta u$ 的符号 | 总通量对 $u_f$ 迎风（B4.2）；按 $\Delta u$ 取会产生相间逆流 |
| 不可容许时 | 逐格限制 $Q^e$、退回一阶、压力下限 | 同一步长整步退回一阶（MOOD），仍不行再减半 $\Delta t$；计数报告（B10，V14） |
| 代码结构 | 一维专用 | `2DEuler_TV` 结构（`Grid`、`Qe/Qs/Hs`、二维）；显式与 IMEX 是两个程序（`src/SPB6_explicit.cpp`、`src/SPB6_TV_IMEX.cpp`），共用 `src/spb6_physics.h`（V9 起直接在 `2DEuler_TV.cpp` 上改）；所有算例同一套格式与参数（V7） |

---

## B13. 数值验证记录（与精确解 / 论文结果比较）

**设置。**

- $N=1000$，二阶。
- IMEX：物质 CFL 0.2，LSDIRK2（`rk=2`）。
- 显式参考（MODE 0，SPB 原格式）：声学 CFL 0.4。
- 精确解：
  - `tests/exact_sg.py`：两相 SG Euler；
  - `tests/exact5.py`：5 方程，SPB (II.5) 激波关系，有真空时给出真空；
  - `tests/exactfan.py`：$pT\to pTg$ 平衡扇（P-S 的 SG 水）；
  - `tests/exact_relaxed.py`：4 方程（$pT$ 松弛）的扇与带激活的 $pT\to pTg$ 扇，SG 与 NASG（V8）；`tests/exact5.py` 也加了 NASG 余容。
- 复现：
  - `bash tests/run_suite.sh OUT all`（编译并运行 `src/spb6_ex` 与 `src/spb6_imex`）；
  - `python3 tests/metrics.py CASE files`；
  - `python3 tests/relax0d.py`；
  - 论文网格上的对照图：`bash tests/run_paper_figs.sh run/figs`，`python3 tests/plot_paper_figs.py run/figs docs/figures`。

各数字来自当前代码（V11，见 C. 开发记录）。与论文各图一一对应的结果（论文网格、变量与论文相同、显式与 IMEX 并列）见 D. 测试结果（`docs/results.md`）。

**参考解的选取。**

- SPB 4.5 空化管的 5 方程精确解有真空：逃逸速度 5.79 m/s < 100 m/s，真空区 $x\in[0.326,0.674]$。$N=4000$ 的显式解只是数值解（中心 $p=332$ Pa，精确为 0）。
- SPB 4.3.1 的微量气体（$\alpha=10^{-6}$）在稀疏波中膨胀约 233 倍，所以用 5 方程精确解（$p^*=4.802\times10^5$，$u^*=492.14$），不用纯相解。
- P-S ±2 m/s 的中心尖峰：精确解（`tests/exactfan.py`）就是 P-S 的 $p$–$pT$–$pTg$ 模型（带 $T>T_{\rm sat}$ 激活）的解：先是 $pT$ 稀疏波降到 $p_{\rm sat}(T)$，再接 $pTg$ 蒸发扇，扇的特征速度只有几 m/s，所以 $x=0.5$ 处是宽约 1 cm 的尖峰。IMEX 的中心值随网格收敛到精确值；显式 HLLC 收敛得慢得多，P-S 图 1 的 $N=5000$ 值也明显偏低。这是低马赫耗散的差别：HLLC 对这个慢结构的耗散 $\sim c\Delta x$（$c\approx1500$ m/s），IMEX 只有物质速度耗散 $\sim|u|\Delta x$（$|u|\le2$ m/s）。
- Pelanti 8.2.1：5 方程精确解（含 NASG 余容）中心 4795 Pa，4 方程（$pT$）13683 Pa，带激活的 $pT\to pTg$ 扇 47724 Pa（$\alpha_2=0.698$，$Y_2=7.02\times10^{-4}$，$T=352.632$ K）。论文图 5 的 $pTg$ 中心温度约 352.94 K，对应约 48.3 kPa：论文的 5000 格解同样抹低了这个蒸发尖峰。
- Pelanti 8.2.2（NASG 十二烷）：$pTg$ 的液体平台就是液体温度（593 K）对应的饱和压力。饱和曲线对表 4 只印到四位的参数极其敏感：液体 $\gamma$ 取 1.2055–1.2065（印出来都是 1.206）时 $p_{\rm sat}(593\ {\rm K})$ 从 $2.3\times10^5$ 变到 $4.8\times10^5$ Pa。表 4 原样给出 $3.33\times10^5$ Pa，我们的平台 $3.37\times10^5$；论文 $4.7\times10^5$。$\gamma_l=1.2065$ 时 IMEX 给出 $4.78\times10^5$ Pa、$\max u=303$ m/s，与论文（约 295 m/s）一致。代码仍用表 4 原样的参数（D 页图 6）。

| 算例 | 指标 | 显式 | IMEX |
|---|---|---|---|
| SPB 4.1 界面平移 | $\max\lvert p/p_0-1\rvert$ | $2.9\times10^{-11}$ | $3.8\times10^{-9}$（Newton 容差） |
| SPB 4.3.1 水–气 | $\max u/u^*$；$L_1(u)$ | 1.0019；1.82 | 1.0015；1.87 |
| SPB 4.3.1，其他时间积分 | $\max u/u^*$ | — | `rk=3` 1.0015；`rk=1` 1.0063；CFL 0.05：1.0017 |
| SPB 4.3.2 水–气 1 TPa | $L_1(u)$ | 16.8 | 84.1（未解决） |
| SPB 4.5 空化 | 扇区 $L_1(u)$；中心 $p$（精确 0） | 0.21；0.18 Pa | 0.066；0.32 Pa |
| P-S 十二烷（只 $p$） | 平台 $L_1(p)/p^*$；$L_1(u)$ | 0.059；0.48 | 0.056；1.03 |
| SPB 5.1 环氧–尖晶石 | $L_1(u)$；$\rho^*$（精确 2527.5） | 3.78；2531 | 5.50；2547；0 次拒绝 |
| SPB 5.2 | $\rho^*$（精确 3767） | 3942 | 3978；两者都不收敛到精确解，与 SPB 5.2 的结论一致 |
| P-S ±2 m/s，$p$ / $pT$ / $pTg$ | 中心 $p$（精确 4010 / 11847 / 50306） | 3414 / 11331 / 50924 | 4067 / 11817 / 50531 |
| P-S ±2 m/s，$pTg$ | 中心尖峰 $\alpha_v$ / $Y_v$（精确扇 0.735 / $7.65\times10^{-4}$；P-S 图 1，$N=5000$：约 0.33 / $1.3\times10^{-4}$） | $N=1000$：0.38 / $1.7\times10^{-4}$（V8 时 $N=2000$、5000：0.51、0.68） | $N=1000$：0.66 / $5.5\times10^{-4}$（V8 时 $N=2000$、5000：0.72、0.74） |
| P-S ±2 m/s，$N=5000$（V8） | 中心 $p$：$p$ / $pT$（精确 4010 / 11847） | 4021 / 11854 | 4012 / 11847 |
| 十二烷 $pTg$ / $p$-$pTg$ | 对 $N=4000$ 显式：$L_1(u)$，$\max u$（参考 317.3） | 1.65，317.5 | V11：$\max u$ 314.7（V8：4.21，315.2） |
| P-S ±500 m/s，$pTg$ | 中心 $p$（带激活的 HEM 精确扇 664；P-S 图 2 约 $1.1\times10^4$） | $N=1000$：$3.0\times10^4$；$N=5000$（V8）：$9.0\times10^3$ | $N=1000$：$7.8\times10^4$（V10 为 $9.3\times10^4$）；$N=5000$（V8）：$2.7\times10^4$（起步加热，未解决） |
| P-S ±500 m/s，$p$ | 中心 $p$（精确为真空） | $N=1000$：0.09 Pa；V10 起没有起步尖峰（V8 时 $N=5000$ 在 $x\approx0.11$ 处有 $1.8\times10^5$ Pa 的尖峰，C 页 V8 §8.7、V10） | $N=1000$：0.71 Pa（V8 时 $N=5000$ 为 5390 Pa） |

以上为当前代码（显式 V10 起未改，IMEX V11），$N=1000$；标注 V8 的是当时的结果。

**Gresho 涡（$80^2$，转一圈，单相极限；B4.1(c)、定理 3）。** 动能保持率：

> **注（V18 后）：** Gresho 涡算例（50–55）已从代码与结果中删除（你的决定：SPB、P-S、Pelanti、BB 都没有这个算例）。本节引用的 Gresho 数字是推导当时的检验记录；低马赫渐近保持现在由 Battisti–Boscheri §4.4 的 Taylor–Green 涡（算例 57、58，C 页 V17）检验。

| | $M=0.1$ | $M=0.01$ | $M=0.001$ |
|---|---|---|---|
| 显式 HLLC（声学 CFL 0.4） | 0.83 | 0.47 | — |
| IMEX（物质 CFL 0.2，带起步斜坡） | 0.926 | 0.937 | 0.937（0 次拒绝，911 步） |

$M=0.001$ 时 $\rho-1\in[-1.1\times10^{-7},3.6\times10^{-8}]$，即 $O(M^2)$。

同一个固定 $\Delta t=0.0025$（503 步）下 $M=10^{-1}\ldots10^{-6}$（`tests/run_gresho_ap.sh`、`tests/plot_gresho_ap.py`，图 `docs/figures/gresho_ap.png`）：

| $M$ | $10^{-1}$ | $10^{-2}$ | $10^{-3}$ | $10^{-4}$ | $10^{-5}$ | $10^{-6}$ |
|---|---|---|---|---|---|---|
| 拒绝次数 | 0 | 0 | 0 | 0 | 0 | 0 |
| KE/KE$_0$ | 0.9303 | 0.9372 | 0.9372 | 0.9372 | 0.9372 | 0.9372 |
| $\max\lvert\rho-1\rvert/M^2$ | 0.113 | 0.102 | 0.101 | 0.103 | 0.035 | 0.106 |
| $\lVert p_2-p_2^{\rm inc}\rVert_2$ | $1.45\times10^{-2}$ | $1.39\times10^{-2}$ | $1.40\times10^{-2}$ | $1.41\times10^{-2}$ | $1.41\times10^{-2}$ | $1.71\times10^{-2}$ |
| $\lVert u_M-u_{M/10}\rVert_2$ | $1.8\times10^{-3}$ | $1.9\times10^{-6}$ | $1.6\times10^{-7}$ | $2.0\times10^{-9}$ | $3.8\times10^{-7}$ | |

$M\le10^{-5}$ 的数字受双精度的底限制（C 页 V11 §11.8）。

**Mie–Grüneisen 激波管（V12，$N=500$，与一般 EOS 的精确 Euler 解比）。** SPB 4.2 硝基甲烷（Cochran–Chan，20 GPa | 0.2 MPa，67 μs）：精确 $p^*=9.827$ GPa、$u^*=1770.6$ m/s；$L_1(p)/p^*$ 显式 0.0030、IMEX 0.0062（V14：0.0049），$L_1(u)/u^*$ 0.0031 / 0.0063（V14：0.0049）。Pelanti 8.1.1（JWL 爆轰气体 | NASG 水，0.9 ms）：精确 $p^*=2.191\times10^{10}$ Pa、$u^*=2397$ m/s；$L_1(p)/p^*$ 0.0035 / 0.0064（V14：0.0045），$L_1(u)/u^*$ 0.0025 / 0.0054（V14：0.0037）。V12 时记录的 $L_1(p)$ 按左态压力而不是 $p^*$ 归一，偏小约一半（V14 更正）。两者接触处都没有振荡（图 `spb6.png`、`spb7.png`、`pel1.png`）。

**Pelanti 8.2.1（NASG 水，±2 m/s，$t=3$ ms，$N=5000$，V8 时的结果，即 V10 之前）。** 中心压力，单位 Pa：

| | 论文（读图） | 显式 | IMEX（178 步） |
|---|---|---|---|
| $p$ | $\approx5\times10^3$ | 4803 | 4795 |
| $pT(f)$，$\vartheta=2000$ | — | 8886 | 8849 |
| $pT$ | $\approx1.3\times10^4$ | 13688 | 13683 |
| $pTg(f)$，$\nu=10^{-4}$ | $\approx3.6\times10^4$ | 35874 | 35855 |
| $pTg$ | $\approx4.8\times10^4$ | 47894 | 47727 |

精确解：$p$ 4795（5 方程），$pT$ 13683（4 方程），$pTg$ 47724（带激活的 $pT\to pTg$ 扇）。$pTg$ 中心尖峰 $\alpha_2$ / $Y_2$：显式 0.638 / $5.4\times10^{-4}$，IMEX 0.706 / $7.27\times10^{-4}$，精确 0.698 / $7.02\times10^{-4}$。

当前代码 $N=1000$：IMEX（58 步）4879 / 8820 / 13677 / 35678 / 47954，显式 4264 / 8687 / 13321 / 35911 / 48277。有限速率的显式参考在 V6 中改正了松弛分裂（A3），修正后 $N=1000$ 为 $pT(f)$ 8639、$pTg(f)$ 35911。

**0D 松弛（算例 40、41，`tests/relax0d.py`）。**

| | 结果 |
|---|---|
| 冻结极限 $\mu=0$ | 状态逐位不变；$\theta=0$ 时两相温度保持不等 |
| AP 极限 | $\mu\ge10^{-2}$、$\theta\ge10^{10}$ 时到舍入误差；$\nu=10,10^2,10^3,10^4,10^6$ 时 $\lvert p-p_\infty\rvert=1.9\times10^3,200,20,2,0.02$ Pa，即 $O(1/\nu)$ |
| 有限 $\mu=10^{-5}$ | 对 Radau 参考：`rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.00） |
| 有限 $\theta=5\times10^6$（$\mu=\infty$） | `rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.03、2.02、2.00、1.95） |
| 有限 $\nu=10^{-3}$（$\mu=\theta=\infty$） | `rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.00） |
| $\nu$ 刚性中间区 | 见 B10（LSDIRK2 越过饱和后被只蒸发的激活条件冻结；`rk=1` 单调） |

**光滑解的观测阶（V14，`limiter=0`，内部 $0.15<x<0.85$，对同一格式 $N=3200$，$N=400\to800$）。**

| | $L_1$：$p$ / $u$ / $\alpha_2$ | $L_2$：$p$ / $u$ / $\alpha_2$ |
|---|---|---|
| 气（算例 20），IMEX | 2.07 / 2.07 / 2.07 | 2.08 / 2.08 / 2.07 |
| 气（算例 20），显式 | 2.06 / 2.06 / 2.07 | 2.06 / 2.07 / 2.07 |
| 水（算例 21），IMEX | 2.04 / 2.05 / 2.08 | 2.04 / 2.05 / 2.08 |
| 水（算例 21），显式 | 1.96 / 1.87 / 2.00 | 1.87 / 1.81 / 1.98 |

V13 及以前的“IMEX 1.2–1.4”是限制器开着测的（minmod 在光滑极值处一阶，共用限制因子又把补足项在任一量的极值处变成一阶）。

**未解决（V14 之后）：**

- SPB 4.3.2 的尾部误差：$L_1(u)$ 84 → 46（V14），显式 HLLC 17；剩下的是 Rusanov 型物质速度耗散与 HLLC 的差别；
- $M\le10^{-5}$ 的双精度底（不修）；活塞为固定网格近似；
- NASG 十二烷平台窗口误差 0.036 → 0.070（全场 $L_1(p)$ 2.44 → 2.00）。

V14 解决：光滑区阶（试验方法）、±500 m/s 起步加热、SPB 4.5 扇区误差、十二烷微量汽相拒绝（C. 开发记录 V14）。

已解决：SPB 4.3.1 的界面速度过冲（V10，1.0014）；SPB 4.5 中心压力（V10）。
