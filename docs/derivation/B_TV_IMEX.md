# B. TV 分裂 IMEX：从要求到离散

目标：$\Delta t$ 只受物质速度限制，$\Delta t=\mathrm{CFL}\,\Delta x/\max|u|$。

- 声波（压力）部分隐式。
- 三种松弛也在同一个隐式闭合里完成。
- 速率 $\mu,\theta,\nu$ 是闭合的参数：
  - 速率 $\to\infty$ 时严格收敛到平衡闭合（渐近保持，AP）；
  - 速率 $\to0$ 时退化为冻结格式。

本章记号同 A 章。一个 IMEX-RK 级内的隐式步长记 $\tau=a_{ii}\Delta t$，$r=\tau/\Delta x$。

本页与 Notion 页 B 同步（2026-10-10，V6）。每一处修改的来由与被否定的做法见“C. 开发记录”（`docs/devlog.md`）。

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
F^e_{i+1/2}[X]=\tfrac12\big(X^{L}u^{L}+X^{R}u^{R}\big)-\tfrac12a_\Omega\big(X^{R}-X^{L}\big),\qquad a_\Omega=\max_{f}\max\big(|u^L_f|,|u^R_f|\big) ,\tag{B3.1}
$$

$X^L=X_i+\tfrac12\sigma_{X,i}$，$X^R=X_{i+1}-\tfrac12\sigma_{X,i+1}$，$u^{L,R}=(\rho u)^{L,R}/\rho^{L,R}$。

**耗散速度取面态速度的最大值（单调性条件）。** 线性重构下格均值是两个面值的平均，$X_i=\tfrac12(X^R_{i-1/2}+X^L_{i+1/2})$，所以显式更新可以写成面值的组合：
$$
X^*_i=X^L_{i+1/2}\Big[\tfrac12-\tfrac r2\big(a+u^L_{i+1/2}\big)\Big]+X^R_{i-1/2}\Big[\tfrac12-\tfrac r2\big(a-u^R_{i-1/2}\big)\Big]+\tfrac r2\big(a-u^R_{i+1/2}\big)X^R_{i+1/2}+\tfrac r2\big(a+u^L_{i-1/2}\big)X^L_{i-1/2}.
$$
四个系数都非负时，$X^*$ 是非负面值的非负组合，微量相保正。充分条件是

- $a\ge|u^{L}_f|$ 且 $a\ge|u^{R}_f|$（每个面、两侧的面态），这一条也是必要的：$a<u^R_{i+1/2}$ 时第三个系数为负，取 $X^R_{i+1/2}$ 很大、其余面值很小就得到负值；
- $a\,r\le\tfrac12$。

MUSCL 的面速度 $(\rho u)^{L,R}/\rho^{L,R}$ 不是格速度的凸组合，可能超出 $\max_i|u_i|$，所以 $a_\Omega$ 必须对面态取最大，不能只对格均值取。只对格均值取时，下游面值 $X^R_{i+1/2}$ 以负系数进入本格：水–气界面处面速度只超出 0.34 m/s，界面格的微量气体质量就变负，后向 Euler（`rk=1`）因此卡住（C. 开发记录 V6）。时间步按格速度取（物质 CFL 0.2），$a\,r$ 远小于 $\tfrac12$。

- **重构在被输运的守恒量上做，各量共用一个限制因子**。记 $\mathcal X=(\alpha_1,m_1,m_2,\varepsilon_1,\varepsilon_2,\rho u,E)$，格 $i$ 的斜率为
  $$\sigma_{X,i}=\phi_i\,\tfrac12(X_{i+1}-X_{i-1}),\qquad \phi_i=\min_{X\in\mathcal X}\frac{\mathrm{minmod}(X_i-X_{i-1},X_{i+1}-X_i)}{\tfrac12(X_{i+1}-X_{i-1})}\in[0,1].$$
  - 对每个量，$|\sigma_X|$ 都不超过它自己的 minmod 斜率且同号，所以仍是 TVD。
  - 所有量的斜率是**同一个**系数乘中心差分，所以量之间的线性关系在面上原样成立。均匀 $p,u$ 时 $\varepsilon_k=c_k(p)\alpha_k+q_km_k$、$\rho u=u\rho$、$E=\sum\varepsilon_k+\tfrac12u^2\rho$ 在面值上同样成立（(E1)），这是 WB 所需的。
  - 重构不经过 EOS，所以 $Q^e$ 只需 $\rho>0$；这一点对 IMEX-RK 的外推级态很重要（B10）。
  - 在模板内为常数（到舍入）的量不参与取最小（$0/0$）：$|X_i-X_{i-1}|,|X_{i+1}-X_i|\le10^{-12}\,S_X$，舍入尺度 $S_X$ 见下表；任一量在格 $i$ 取极值时 $\phi_i=0$。
  - $\phi_i$ 也用来判断“所有被输运量在格 $i$ 都光滑”（B4.1 的二阶声学耗散）。
- $\alpha_2$ 的面值取 $1-\alpha_1$，所以 $\sum_kF^e[\alpha_k]=\tfrac12(u^L+u^R)\equiv\bar u^e_f$，耗散两相抵消。

**限制因子的舍入尺度 $S_X$（V8）。** $S_X$ 是格式产生该量舍入的尺度，不总是该量自身的大小（三个格求和）：

| 量 | $S_X$ | 理由 |
|---|---|---|
| $\alpha_1$ | 1 | $\sum\alpha=1$ 与闭合在 $O(1)$ 上算，微量相 $\alpha\sim10^{-8}$ 的舍入是 $\varepsilon_{\rm mach}$，不是 $\varepsilon_{\rm mach}\alpha$ |
| $\rho u,\rho v$ | $\lvert\rho\mathbf u\rvert+\rho c_f$ | 动量由压力差驱动，舍入 $\sim\varepsilon_{\rm mach}(\tau/h)p\sim\varepsilon_{\rm mach}\rho c$，与当地 $\lvert\rho u\rvert$ 无关；对向量取模，判据与坐标方向无关 |
| $\varepsilon_k$、$E$ | $\sum_k\lvert\varepsilon_k\rvert+\tfrac12\rho\lvert\mathbf u\rvert^2$ | 闭合在两相间搬运 $p\,d\alpha$，舍入相对于混合物内能 |
| $m_k$ | $\lvert m_k\rvert$ | 无传质时 $m_k$ 只被自身通量输运，舍入相对于自身；也是微量相质量保正所需 |

按各量自身大小判时，静止区（$\rho u\approx0$）与二维的横向动量（一维中 $\rho v\equiv0$，二维中是 $10^{-15}$ 量级的噪声）的舍入噪声会触发极值判据，使 $\phi$ 随舍入翻转；二维 y 向均匀的问题因此各行发散（C. 开发记录 V8）。保正要求微量 $\gg10^{-12}S_X$（算例中微量为 $10^{-6}$–$10^{-8}$）。

**中心部分为什么也用重构值（与 BB (35) 不同）。** BB 的中心部分用格均值、只有耗散用重构值。对线性输运 $X_t+uX_x=0$（$u>0$，$a_\Omega=u$），BB 的通量是
$$F=\tfrac12u(X_i+X_{i+1})-\tfrac12u(X^R-X^L)=uX_i+\tfrac14u(\sigma_i+\sigma_{i+1}),$$
含下游格的斜率 $\sigma_{i+1}$，不是 TVD；微量相（$\alpha\sim10^{-8}$）界面处会把 $\alpha$ 推成负值（数值上观察到）。中心部分也用重构值时 $F=u(X_i+\tfrac12\sigma_i)$，就是 MUSCL 迎风通量，TVD。两者都是二阶。

**与 BB 的对应：** BB (36)–(38) 要求，对 $\phi$ 加的人工耗散必须在内能方程里按 $\partial(\rho\phi e)/\partial\phi$ 缩放后同样加上，否则接触处出现压力振荡。本格式对 $\alpha_k$ 与 $\varepsilon_k$ 用同一个 $a_\Omega$ 和同一种重构；均匀 $p$ 时由 (E1) 有 $\Delta\varepsilon_k=\frac{p+\gamma\pi}{\gamma-1}\Delta\alpha_k+q\Delta m_k$，这正是 BB (38) 的缩放，而且对任意 $\alpha,m$ 分布都自动满足。

BB (39) 的压力跳跃耗散（高马赫时打开）没有单独实现：它的作用由 B4.1 的声学迎风项承担，后者同样只在声波被时间步分辨（$c\tau<\Delta x$，即高马赫或起步）时打开。

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
p_f&=\frac{\rho^*_{i+1}P_i+\rho^*_iP_{i+1}}{\rho^*_i+\rho^*_{i+1}}-\varphi_f\,\psi^u_f\,\frac{Z_iZ_{i+1}}{Z_i+Z_{i+1}}\big(u^*_{i+1}-u^*_i\big),
\end{aligned}}\qquad \varphi_f=\Big(1-\frac{\tau\min(c_i,c_{i+1})}{\Delta x}\Big)^+ ,\tag{B4.1}
$$
$Z=\rho^*c$，$\bar u^*_f=\tfrac12(u^*_i+u^*_{i+1})$。

- $\tau\to0$：$\varphi=1$，$u_f-\bar u^*=-\Delta P/(Z_i+Z_{i+1})$，$p_f-\bar p=-\frac{Z_iZ_{i+1}}{Z_i+Z_{i+1}}\Delta u$，就是声学 Riemann 解（Godunov）。
- $c\tau\ge\Delta x$：$\varphi=0$，格式不变。低马赫时 $c\tau/\Delta x=\mathrm{CFL}/M\gg1$，所以 AP 性质（B7）与低马赫精度不受影响。
- $\varphi>0$ 当且仅当时间步分辨声波（$c\tau<\Delta x$，即 $M>\mathrm{CFL}$，或起步阶段）。这与 BB 的原则一致：声学型耗散 (39) 只在高马赫时打开，低马赫时只有物质速度耗散。耗散率不超过 $c\Delta x/2$。
- 加上后奇偶模态在 $\theta=\pi$ 的放大因子为 $|(1-2\varphi\nu)(1-2C)|$（例如 $C=0.2,M=0.5$：0.31，Godunov 0.60）。
- 速度跳跃项是显式粘性，$\nu_{\rm add}\tau/\Delta x^2=\tfrac12\nu(1-\nu)\le\tfrac18$，稳定；压力项对 $P$ 隐式、线性。
- **$c$ 取冻结声速。** $\rho c_{\rm frozen}\ge\rho c$ 对所有松弛闭合（Wood、$pT$、$pTg$）成立（次特征条件），所以补足后的耗散对任何闭合都不低于其波的 Godunov 耗散。用 Wood 声速试过：数值混合区里 $Z$ 很小，$\Delta P/(Z_i+Z_{i+1})$ 很大，把界面推得过快（水–气 $\max u/u^*=1.11$，激波超前 20 格），舍弃。
- **二阶：** $\psi_f=1-(1-\psi^{\rm mm}_f)\min(\phi_i,\phi_{i+1})$，$\psi^{\rm mm}_f=1-\frac{\sigma_i+\sigma_{i+1}}{2D_f}$（$D_f$ 为面跳跃，$\sigma$ 为 minmod 斜率，所以 $\psi^{\rm mm}\in[0,1]$，光滑区 $O(\Delta x)$）；$P$ 的斜率取上一次的压力场（Newton 初猜），$u$ 的取 $u^*$，残量对 $P$ 仍是线性的。$\phi$ 是 B3 的共用限制因子：只有所有被输运量都光滑时才降到二阶。接触间断处 $p,u$ 连续（$\psi^{\rm mm}\to0$）但阻抗跳跃，若只看 $p,u$，耗散在最需要的地方被关掉（水–气 CFL 0.05 时界面格气体质量变负，数值上观察到）。
- $c$、$Z$ 在显式预测态 $X^*$ 上计算，一级内不变（残量对 $P$ 仍是线性的）；$X^*$ 不可容许（某相 $\alpha^*$、$m^*\le0$ 或 $p+\pi\le0$）时取 $c=\infty$，即 $\varphi=0$，该面只剩半格动量平衡。固壁面 $u_f=0$，增量为零。
- **舍弃的另一条路：** 让两侧参与面动量平衡的质量取 $\rho\min(c\tau,\Delta x/2)$（声学影响域，不超过半格），$\tau$ 小时同样是 Godunov，$\tau$ 大时是带非弹性碰撞项的半格平衡。但大步长时它每一级都把格速度换成两侧面速度的平均（$\tfrac14\delta^2u$），耗散率 $\propto\Delta x^2/\tau$ 随 $a_{ii}$ 变化（LSDIRK2 是后向 Euler 的 3.4 倍），空化管的熵误差因此成倍增大。

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

$\delta$ 与松弛项无关（$\chi_k$ 项两相相消）。光滑流里它是离散误差，激波处是 $O(1)$。把它按等容加热的比例分给两相：

$$
\varepsilon_k=\varepsilon^{\rm nc}_k+\lambda_k\delta,\qquad
\lambda_k=\frac{\partial\varepsilon_k/\partial p|_{\alpha,m}}{\sum_j\partial\varepsilon_j/\partial p|_{\alpha,m}}=\frac{(\alpha^{\rm tr}_k-m^{\rm tr}_kb_k)/(\gamma_k-1)}{\sum_j(\alpha^{\rm tr}_j-m^{\rm tr}_jb_j)/(\gamma_j-1)} .\tag{B5.3}
$$

于是 $\sum_k\varepsilon_k=(\rho e)^C$ 恒成立。

- SPB 的重初始化是在 $\alpha$ 固定时整体平移 $p$，也就是按同样的 $\lambda_k$ 分配 $\delta$；区别只是 SPB 分两步做，这里在闭合里一起解。
- 所以不需要级末投影，解出来的状态就是混合能量一致的（P-S 2014 的定义 3.1）。
- $\lambda_k$ 在一次闭合求解中冻结在输运态。

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

---

## B8. 全局非线性方程组与求解

$G_i(P)=\bar p^{\rm cl}_i-P_i$ 只依赖 $(P_{i-1},P_i,P_{i+1})$（二维 5 点），Jacobian 三对角（二维 5 对角）。

**单调性：**

- $\partial s_i/\partial P_i=r^2(1/\bar\rho_{i+1/2}+1/\bar\rho_{i-1/2})+r\sum_f\varphi_f\psi^P_f/(Z_i+Z_{i\pm1})>0$（B4.1 的补足项只增加对角占优），且 $\partial\bar p^{\rm cl}/\partial s<0$：膨胀降压（E5）。
- $\partial E_i/\partial P_i\approx-r^2\,2\bar p/\bar\rho<0$。

所以 $\partial G_i/\partial P_i<-1$，相邻项 $O(r^2\rho c^2/\bar\rho)$，矩阵是 $I-\tau^2\partial_x(c^2\partial_x)$ 型：对角占优的 M 矩阵。

**求解：**

- Newton。Jacobian 用逐格差分：扰动三个（五个）压力，各算一次局部闭合。
- Jacobian 的差分步长 $h=10^{-7}(|P|+1)$；Newton 至多 60 次。
- 一维非周期用 Thomas；二维与周期边界用 BiCGStab + Jacobi 预条件（`2DEuler_TV` 风格，但用相对容差 $10^{-13}$）。
- 边界虚格（透射、固壁）的压力等于本格压力。它在 Jacobian 里是单独的元素（列号为本格），不并进对角元：对角元只扰动本格，虚格的导数单独差分；每行按“本格、W、E、S、N”的固定次序存 5 个元素，同列的也不合并；Jacobi 预条件用本格的那个元素，Thomas 把同列元素相加。这样数学上 Jacobian 不变，而边界行与内部行的浮点运算相同：y 向均匀的问题各行逐位相同（V8；并进对角元时边界行的舍入不同，再经共用限制因子的开关放大成截断误差大小的行间差）。
- 线搜索（步长减半，至多 40 次）：接受条件是所有格的局部闭合有解（可容许），且残量范数满足 Armijo 下降 $\|G\|_{\rm new}<(1-10^{-4}\lambda)\|G\|$。
- 收敛判据：$\max_i|G_i|/(|P_i|+K_i+P_{E,i})<10^{-13}$，$K=\rho c_W^2$（Wood 体积模量，$G/K$ 就是体积约束残量 $\sum\alpha-1$），$P_E=\sum_k(|\varepsilon_k|+|m_kq_k|)/\sum_k\partial\varepsilon_k/\partial p$（$X^*$ 处；代码里分母再加 1 Pa 防零）。理由：闭合给出的压力只精确到 $\varepsilon_{\rm mach}(p+\pi)$，守恒总能的舍入 $\varepsilon_{\rm mach}|E|$ 又经 $\lambda_k\delta$ 进入各相压力，所以 $G$ 的可达精度是 $\varepsilon_{\rm mach}\max(K,P_E)$。按 $|P|$ 归一时 $p\ll\pi$ 的格子永远达不到容差：环氧–尖晶石（$\pi=1.41\times10^{11}$）算例因此每步都被拒绝（3659 次拒绝、13906 步），改后 0 次拒绝、122 步，每级 3 次 Newton。

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

**起步：** 从静止或间断初值起算时，首步用声学步长（`cflac`，默认 0.4），之后每步至多放大 1.2 倍，直到物质 CFL（`2DEuler_TV` 的 `full_step` 与 BB 的做法）；被 1.2 倍限制的步用一阶（单级后向 Euler）。`dt=` 给定固定步长（时间收敛阶试验用，不做斜坡）。

---

## B11. 二维

- 显式：$x,y$ 两个方向的面各算一次 (B3.1)。
- 面速度 (B4.1) 只修正法向分量；增量质量通量携带新的法向面速度 $u_f$，以及按供体取的切向速度 $\hat v$。
- $s_{ij}=r_x\delta_x[u_f]+r_y\delta_y[v_f]$。
- 动量：$x$ 方向 $-r_x\delta_x[p_f]-r_x\delta_x[u_f\hat M]-r_y\delta_y[\hat u\hat M]$，$y$ 方向对称。
- 总能的增量动能取 $\tfrac12(u_f^2+\hat v^2)\hat M$。
- 全局 Jacobian 5 对角。
- **二维检查（V8）：** 一维算例放在 $N\times4$ 网格上（y 向透射、均匀），各行应与一维结果相同、$v\equiv0$。12 个算例（SPB 4.1、4.3.1、4.3.2、4.5、5.1，P-S ±2、±500 m/s，十二烷 SG/NASG 的 $p$ 与 $pTg$，Pelanti $pTg(f)$，光滑解）各行逐位相同，$v\equiv0$。固定步长时第 0 行与一维结果之差：SPB 5.1 为 $4\times10^{-9}$，4.3.1 为 $10^{-6}$（一维用 Thomas、二维用 BiCGStab，舍入不同；欠分辨的 P-S 中心尖峰处经限制器开关放大到 $10^{-2}$，仍远小于离散误差）。显式程序同样各行逐位相同。真正二维的算例（如 SPB 4.7、P-S 的气泡压缩）尚未运行。

---

## B12. 与 handoff v23 的差别（换路径的地方）

| | handoff v23 | 本推导 |
|---|---|---|
| 混合能量一致性 | 闭合用非守恒相能量，级末再用总能重初始化 | $\delta$ 进入闭合（B5.2），没有级末投影；动量方程看到的就是守恒能量下的压力 |
| 松弛 | 瞬时（$p$、$pT$、$pTg$）在闭合里；有限速率没有 | 一个闭合，$\mu,\theta,\nu$ 是参数（B6），两端极限严格（B7） |
| 全局残量 | $\sum\alpha-1$ 或能量残量，按闭合类型不同 | $G=\bar p^{\rm cl}-P$，所有闭合同一形式 |
| 传质判据 | 上一级的标记（滞后） | 每步开始按 $U^n$ 判断并投影到平衡流形，本步各级冻结（B6.7） |
| 显式耗散 | 中心部分也用重构值 | 同样中心部分也用重构值（BB 的格均值中心部分不是 TVD，B3）；耗散速度取物质速度，且不小于面态速度（B3 的保正条件） |
| 声学耗散 | 无（只有后向 Euler） | 声学 CFL < 1 时补足到 Godunov（B4.1），否则不变 |
| 增量供体 | 按 $\Delta u$ 的符号 | 总通量对 $u_f$ 迎风（B4.2）；按 $\Delta u$ 取会产生相间逆流 |
| 不可容许时 | 逐格限制 $Q^e$、退回一阶、压力下限 | 拒绝本步并减半 $\Delta t$，计数报告 |
| 代码结构 | 一维专用 | `2DEuler_TV` 结构（`Grid`、`Qe/Qs/Hs`、二维）；显式与 IMEX 是两个程序（`src/SPB6_explicit.cpp`、`src/SPB6_TV_IMEX.cpp`），共用 `src/spb6_common.h`；所有算例同一套格式与参数（V7） |

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

各数字来自当前代码（V8，见 C. 开发记录）。与论文各图一一对应的结果（论文网格、变量与论文相同、显式与 IMEX 并列）见 D. 测试结果（`docs/results.md`）。

**参考解的选取。**

- SPB 4.5 空化管的 5 方程精确解有真空：逃逸速度 5.79 m/s < 100 m/s，真空区 $x\in[0.326,0.674]$。$N=4000$ 的显式解只是数值解（中心 $p=332$ Pa，精确为 0）。
- SPB 4.3.1 的微量气体（$\alpha=10^{-6}$）在稀疏波中膨胀约 233 倍，所以用 5 方程精确解（$p^*=4.802\times10^5$，$u^*=492.14$），不用纯相解。
- P-S ±2 m/s 的中心尖峰：精确解（`tests/exactfan.py`）就是 P-S 的 $p$–$pT$–$pTg$ 模型（带 $T>T_{\rm sat}$ 激活）的解：先是 $pT$ 稀疏波降到 $p_{\rm sat}(T)$，再接 $pTg$ 蒸发扇，扇的特征速度只有几 m/s，所以 $x=0.5$ 处是宽约 1 cm 的尖峰。IMEX 的中心值随网格收敛到精确值；显式 HLLC 收敛得慢得多，P-S 图 1 的 $N=5000$ 值也明显偏低。这是低马赫耗散的差别：HLLC 对这个慢结构的耗散 $\sim c\Delta x$（$c\approx1500$ m/s），IMEX 只有物质速度耗散 $\sim|u|\Delta x$（$|u|\le2$ m/s）。
- Pelanti 8.2.1：5 方程精确解（含 NASG 余容）中心 4795 Pa，4 方程（$pT$）13683 Pa，带激活的 $pT\to pTg$ 扇 47724 Pa（$\alpha_2=0.698$，$Y_2=7.02\times10^{-4}$，$T=352.632$ K）。论文图 5 的 $pTg$ 中心温度约 352.94 K，对应约 48.3 kPa：论文的 5000 格解同样抹低了这个蒸发尖峰。
- Pelanti 8.2.2（NASG 十二烷）：$pTg$ 的液体平台就是液体温度（593 K）对应的饱和压力。饱和曲线对表 4 只印到四位的参数极其敏感：液体 $\gamma$ 取 1.2055–1.2065（印出来都是 1.206）时 $p_{\rm sat}(593\ {\rm K})$ 从 $2.3\times10^5$ 变到 $4.8\times10^5$ Pa。表 4 原样给出 $3.33\times10^5$ Pa，我们的平台 $3.37\times10^5$；论文 $4.7\times10^5$。$\gamma_l=1.2065$ 时 IMEX 给出 $4.78\times10^5$ Pa、$\max u=303$ m/s，与论文（约 295 m/s）一致。代码仍用表 4 原样的参数（D 页图 6）。

| 算例 | 指标 | 显式 | IMEX |
|---|---|---|---|
| SPB 4.1 界面平移 | $\max\lvert p/p_0-1\rvert$ | $10^{-9}$ | $2.2\times10^{-8}$（Newton 容差） |
| SPB 4.3.1 水–气 | $\max u/u^*$；激波位置误差 | 1.001；+4 格 | 1.051；+10 格（$N=4000$：1.005） |
| SPB 4.3.1，其他时间积分 | $\max u/u^*$ | — | `rk=3` 1.051；`rk=1` 1.075；CFL 0.05：1.040 |
| SPB 4.3.2 水–气 1 TPa | $L_1(u)$ | 17 | 128（未解决） |
| SPB 4.5 空化 | 扇区 $L_1(u)$；中心 $p$（精确 0） | 0.090；2211 Pa | 0.049；2254 Pa |
| P-S 十二烷（只 $p$） | 平台 $L_1(p)/p^*$；$L_1(u)$ | 0.058；0.48 | 0.043；0.86 |
| SPB 5.1 环氧–尖晶石 | $L_1(u)$；$\rho^*$（精确 2527.5） | 3.70；2530 | 5.40；2547；0 次拒绝 |
| SPB 5.2 | $\rho^*$（精确 3767） | 3953 | 3990；两者都不收敛到精确解，与 SPB 5.2 的结论一致 |
| P-S ±2 m/s，$p$ / $pT$ / $pTg$ | 中心 $p$（精确 4010 / 11847 / 50306） | 3583 / 11331 / 50924 | 4075 / 11815 / 50608 |
| P-S ±2 m/s，$pTg$ | 中心尖峰 $\alpha_v$ / $Y_v$（精确扇 0.735 / $7.64\times10^{-4}$；P-S 图 1，$N=5000$：约 0.33 / $1.3\times10^{-4}$） | $N=1000$：0.38 / $1.7\times10^{-4}$；$N=2000$：0.51 / $2.9\times10^{-4}$ | $N=1000$：0.65 / $5.1\times10^{-4}$；$N=2000$：0.72 / $7.2\times10^{-4}$；$N=5000$：0.74 / $7.8\times10^{-4}$ |
| 十二烷 $pTg$ / $p$-$pTg$ | 对 $N=4000$ 显式：$L_1(u)$，$\max u$（参考 317.3） | 1.65，317.5 | 4.21，315.2 |
| P-S ±500 m/s，$pTg$ | 中心 $p$（HEM 精确 664） | $3.0\times10^4$ | $9.3\times10^4$（未解决） |

**Pelanti 8.2.1（NASG 水，±2 m/s，$t=3$ ms，$N=5000$）。** 中心压力，单位 Pa：

| | 论文（读图） | 显式 | IMEX（178 步） |
|---|---|---|---|
| $p$ | $\approx5\times10^3$ | 4803 | 4795 |
| $pT(f)$，$\vartheta=2000$ | — | 8886 | 8849 |
| $pT$ | $\approx1.3\times10^4$ | 13688 | 13683 |
| $pTg(f)$，$\nu=10^{-4}$ | $\approx3.6\times10^4$ | 35874 | 35855 |
| $pTg$ | $\approx4.8\times10^4$ | 47894 | 47727 |

精确解：$p$ 4795（5 方程），$pT$ 13683（4 方程），$pTg$ 47724（带激活的 $pT\to pTg$ 扇）。$pTg$ 中心尖峰 $\alpha_2$ / $Y_2$：显式 0.638 / $5.4\times10^{-4}$，IMEX 0.706 / $7.27\times10^{-4}$，精确 0.698 / $7.02\times10^{-4}$。

$N=1000$（IMEX 58 步）为 4880 / 8818 / 13676 / 35689 / 47929。有限速率的显式参考在 V6 中改正了松弛分裂（A3），修正后 $N=1000$ 为 $pT(f)$ 8639、$pTg(f)$ 35911。

**0D 松弛（算例 40、41，`tests/relax0d.py`）。**

| | 结果 |
|---|---|
| 冻结极限 $\mu=0$ | 状态逐位不变；$\theta=0$ 时两相温度保持不等 |
| AP 极限 | $\mu\ge10^{-2}$、$\theta\ge10^{10}$ 时到舍入误差；$\nu=10,10^2,10^3,10^4,10^6$ 时 $\lvert p-p_\infty\rvert=1.9\times10^3,200,20,2,0.02$ Pa，即 $O(1/\nu)$ |
| 有限 $\mu=10^{-5}$ | 对 Radau 参考：`rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.00） |
| 有限 $\theta=5\times10^6$（$\mu=\infty$） | `rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.03、2.02、2.00、1.95） |
| 有限 $\nu=10^{-3}$（$\mu=\theta=\infty$） | `rk=1` 一阶，`rk=2`、`rk=3` 二阶（2.00） |
| $\nu$ 刚性中间区 | 见 B10（LSDIRK2 越过饱和后被只蒸发的激活条件冻结；`rk=1` 单调） |

**光滑解（算例 20，对同一格式 $N=3200$）的观测阶。**

| | $p$、$u$ | $\alpha_2$ |
|---|---|---|
| 显式 | 2.0 | 2.0 |
| IMEX | 1.4 | 1.2（未解决） |

**未解决：**

- IMEX 光滑区 $\alpha$ 的阶只有 1.2–1.4，$p$、$u$ 约 1.4（声学补足项的 $\psi$ 在其激活处降阶）；
- SPB 4.3.1 的界面速度过冲约 5%：偏差随时间增大，随网格一阶收敛；
- SPB 4.3.2 的尾部误差；
- ±500 m/s $pTg$ 远离 HEM 精确扇。
