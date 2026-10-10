# A. 6+1 方程模型、状态方程与全显式格式（SPB 2009 + 三种松弛）

记号：相 $k=1$（液），$k=2$（汽/气）。$\alpha_k$ 体积分数，$\alpha_1+\alpha_2=1$；$m_k=\alpha_k\rho_k$ 相质量；
$\varepsilon_k=\alpha_k\rho_ke_k$ 相内能密度（**含**生成焓 $q_k$）；$\rho=\sum_km_k$，$Y_k=m_k/\rho$；
$u$ 公共速度；$E=\sum_k\varepsilon_k+\tfrac12\rho u^2$；$\bar p=\sum_k\alpha_kp_k$ 混合压力；
一维记 $r=\Delta t/\Delta x$，$\delta_i[\phi_f]=\phi_{i+1/2}-\phi_{i-1/2}$。

---

## A1. 模型

Saurel–Petitpas–Berry (2009) 的 6 方程模型 (II.6)，加冗余的混合总能 (II.7)，再加 Pelanti–Shyue (2014) / Pelanti (2022) 的传热传质源项：

$$
\begin{aligned}
&\partial_t\alpha_1+u\,\partial_x\alpha_1=\mathcal P+\mathcal M/\rho_I,\\
&\partial_tm_1+\partial_x(m_1u)=\mathcal M,\qquad \partial_tm_2+\partial_x(m_2u)=-\mathcal M,\\
&\partial_t(\rho u)+\partial_x(\rho u^2+\bar p)=0,\\
&\partial_t\varepsilon_1+\partial_x(\varepsilon_1u)+\alpha_1p_1\,\partial_xu=-p_I\mathcal P+\mathcal Q+e_I\mathcal M,\\
&\partial_t\varepsilon_2+\partial_x(\varepsilon_2u)+\alpha_2p_2\,\partial_xu=+p_I\mathcal P-\mathcal Q-e_I\mathcal M,\\
&\partial_tE+\partial_x\big((E+\bar p)u\big)=0 .
\end{aligned}\tag{A1.1}
$$

$$
\mathcal P=\mu(p_1-p_2),\qquad \mathcal Q=\theta(T_2-T_1),\qquad \mathcal M=\nu(g_2-g_1),\qquad
p_I=\frac{Z_2p_1+Z_1p_2}{Z_1+Z_2},\ Z_k=\rho_kc_k .\tag{A1.2}
$$

- $\mu$ [1/(Pa·s)]、$\theta$ [Pa/(s·K)]、$\nu$ [Pa·kg²/(s·J²)] 是三个速率。$\mathcal M>0$ 表示相 2 → 相 1（$g_2>g_1$ 时凝结）。
- 最后一个方程是前面方程的推论（冗余），但离散后它是**守恒**的，用来修正非守恒的相内能（SPB 的 “+1”）。
- $\rho_I,e_I$ 在本文用到的所有闭合（压力、温度瞬时平衡下的传质）中都会消掉，不需要给定（P-S 2014 §2.3；下面 B5 节再验证）。
- 熵：$\alpha_k\rho_kT_k\,ds_k/dt=(p_k-p_I)(\pm\mathcal P)+\dots$；取阻抗加权的 $p_I$ 时，力学松弛产生的熵 $\ge0$。

**层级极限（Pelanti 2022 §3）：** $\mu\to\infty$ 给 5 方程 Kapila 模型（$p$ 平衡）；再加 $\theta\to\infty$ 给 4 方程 $pT$ 模型；再加 $\nu\to\infty$ 给 HEM（$pTg$）。声速满足次特征条件 $c_{pTg}\le c_{pT}\le c_{p}\le c_f$。

---

## A2. 状态方程（NASG；$b=0$ 时为 stiffened gas）

每相参数 $(\gamma,\pi,q,q',c_v,b)$：

$$
\begin{aligned}
&\rho e=\frac{(p+\gamma\pi)(1-\rho b)}{\gamma-1}+\rho q,\qquad p(\rho,e)=\frac{(\gamma-1)\rho(e-q)}{1-\rho b}-\gamma\pi,\\
&v(p,T)=\frac{(\gamma-1)c_vT}{p+\pi}+b,\qquad e(p,T)=c_vT\,\frac{p+\gamma\pi}{p+\pi}+q,\qquad T(p,v)=\frac{(v-b)(p+\pi)}{(\gamma-1)c_v},\\
&h=\gamma c_vT+bp+q,\qquad s=c_v\ln\frac{T^\gamma}{(p+\pi)^{\gamma-1}}+q',\qquad
g=h-Ts=(\gamma c_v-q')T-c_vT\ln\frac{T^\gamma}{(p+\pi)^{\gamma-1}}+bp+q,\\
&c^2=\frac{\gamma(p+\pi)}{\rho(1-\rho b)} .
\end{aligned}\tag{A2.1}
$$

**按相体积写（关键）：** 一个格子里相 $k$ 占体积 $\alpha_k$、质量 $m_k$，在压力 $p$ 下的内能是

$$
\varepsilon_k(\alpha_k,m_k,p)=\frac{(p+\gamma_k\pi_k)(\alpha_k-m_kb_k)}{\gamma_k-1}+m_kq_k .\tag{A2.2}
$$

状态方程在本算法里承担五个角色：

| 编号 | 性质 | 用在哪里 |
|---|---|---|
| (E1) | 固定 $p$ 时 $\varepsilon_k$ 对 $(\alpha_k,m_k)$ **仿射** | 任何对 $(\alpha_k,m_k,\varepsilon_k)$ 同一个线性输运算子都保持均匀 $p$：接触间断不被破坏（WB，B9 节） |
| (E2) | (A2.2) 对 $\alpha_k$ 可显式反解 | 离散 Gibbs 闭合是代数的（A5、B5 节） |
| (E3) | $v(p,T)$、$e(p,T)$ 显式，$e$ 对 $T$ 线性 | $pT$、$pTg$、有限速率 $\theta,\nu$ 闭合中 $T$ 显式（B6 节） |
| (E4) | $g_1(p,T)=g_2(p,T)$ 给出饱和线 $T_{\rm sat}(p)$ | $pTg$ 闭合、传质判据 |
| (E5) | 凸性：$c_k^2>0$、$\partial v/\partial p\vert _T<0$、$c_p>0$ | 所有局部闭合的残量对 $p$ 单调、根唯一（B6 节） |

**混合压力（SPB (II.4) 的 NASG 推广）：** 给定 $\alpha_k,m_k$ 和混合内能 $\rho e=\sum\varepsilon_k$，按 (A2.2) 对 $p$ 线性：

$$
p=\frac{\rho e-\sum_km_kq_k-\sum_k\gamma_k\pi_k\frac{\alpha_k-m_kb_k}{\gamma_k-1}}{\sum_k\frac{\alpha_k-m_kb_k}{\gamma_k-1}} .\tag{A2.3}
$$

**声速：** 冻结 $c_f^2=\sum_kY_kc_k^2$；Wood $\dfrac{1}{\rho c_p^2}=\sum_k\dfrac{\alpha_k}{\rho_kc_k^2}$；$c_{pT}$、$c_{pTg}$ 见 Pelanti (2022) 式 (19)(25)。

**饱和线：** $g_1(p,T)=g_2(p,T)$ 在 $T$ 上用带区间的 Newton 解出 $T_{\rm sat}(p)$。SG 下这条曲线只在 $v_2(p,T_{\rm sat})>v_1(p,T_{\rm sat})$ 的一段上是真的相平衡；端点 $p_c$ 处 $v_1=v_2$（十二烷 SG：$p_c\approx1.39\times10^7$ Pa）。


### A2.1 Mie–Grüneisen 相（V12：Cochran–Chan、JWL）

SPB §4.2 的硝基甲烷（Cochran–Chan）与 Pelanti §4.1–4.2 的 JWL 都是 MG 形式（Pelanti (27)）：

$$
p(\rho,e)=\Gamma\rho\big(e-\varepsilon_r(\rho)\big)+p_r(\rho),\qquad \Gamma=\text{常数} .\tag{A2.4}
$$

- Cochran–Chan（SPB §4.2）：$x=\rho/\rho_{\rm ref}$，$p_r=A_1x^{E_1}-A_2x^{E_2}$，$\varepsilon_r=\dfrac{A_1x^{E_1-1}}{\rho_{\rm ref}(E_1-1)}-\dfrac{A_2x^{E_2-1}}{\rho_{\rm ref}(E_2-1)}$。
- JWL（Pelanti (29)）：$p_r=a\,e^{-r_1\rho_0/\rho}+b\,e^{-r_2\rho_0/\rho}$，$\varepsilon_r=\dfrac{a}{r_1\rho_0}e^{-r_1\rho_0/\rho}+\dfrac{b}{r_2\rho_0}e^{-r_2\rho_0/\rho}-\varepsilon_0$。
- 两者的参考曲线都是等熵线：$\varepsilon_r'(\rho)=p_r(\rho)/\rho^2$（直接求导可验）。下面所有公式只用这一条。NASG 也是 MG 形式（$\Gamma=(\gamma-1)/(1-\rho b)$，$\varepsilon_r=q$，$p_r=-\gamma\pi$），但 $\Gamma$ 依赖 $\rho$；代码里 NASG 保留 (A2.2) 的闭式，逐位不变。

**按相体积写：** $\rho_k=m_k/\alpha_k$，

$$
\varepsilon_k(\alpha_k,m_k,p)=m_k\varepsilon_r(\rho_k)+\alpha_k\frac{p-p_r(\rho_k)}{\Gamma_k} .\tag{A2.5}
$$

**声速：** $c^2=(\partial p/\partial\rho)_e+\frac{p}{\rho^2}(\partial p/\partial e)_\rho$，用 $\varepsilon_r'=p_r/\rho^2$：

$$
c^2=\frac{(\Gamma+1)(p-p_r(\rho))}{\rho}+p_r'(\rho) .\tag{A2.6}
$$

可容许条件是 $c^2>0$（代替 SG 的 $p+\pi>0$）。

**混合压力（Pelanti (28)）：** (A2.5) 在固定 $\alpha_k,m_k$ 时对 $p$ 线性，

$$
p=\frac{\rho e-\sum_k\big(m_k\varepsilon_r(\rho_k)-\alpha_kp_r(\rho_k)/\Gamma_k\big)}{\sum_k\alpha_k/\Gamma_k},\qquad \frac{\partial\varepsilon_k}{\partial p}\Big|_{\alpha,m}=\frac{\alpha_k}{\Gamma_k} .\tag{A2.7}
$$

**五个角色（A2 的表）对 MG 的变化：**
- (E1) 不再成立：固定 $p$ 时 $\varepsilon_k$ 对 $(\alpha_k,m_k)$ 不仿射（$\rho_k$ 进 $\varepsilon_r,p_r$）。但若每相的密度 $\rho_k$ 在界面两侧相同（初值里微量相取与主相同一个 $\rho_k$），输运只改变 $\alpha_k$、$m_k$ 的比例而不改变 $\rho_k$，$\varepsilon_k$ 对 $\alpha_k$ 线性，WB 仍成立。SPB §3.2.2 末与 §4.2 正是这样做的（同一种流体当两相）。
- (E2) 不再显式：(A2.5) 对 $\alpha_k$ 的反解要一维求根（A6、B6.1）。
- (E3)(E4) 不用：MG 相只做力学松弛（SPB 4.2、Pelanti 8.1 都只有力学松弛）。代码里 MG 相与热、化学松弛的组合会报错。
- (E5) 见下：力学闭合的单调性由 $c^2>0$ 保证。

**相 Hugoniot（显式 HLLC 的相能量，SPB (II.8)）：** $e^*-e+\tfrac12(p^*+p)(v^*-v)=0$，代入 (A2.4) 对 $p^*$ 仍线性：

$$
p^*\Big(\frac{v^*}{\Gamma}+\frac{v^*-v}{2}\Big)=e-\varepsilon_r(\rho^*)+\frac{p_r(\rho^*)v^*}{\Gamma}-\frac{p\,(v^*-v)}{2} .\tag{A2.8}
$$

---

## A3. 全显式格式：总流程（SPB §3.5 + P-S §4 / Pelanti §5）

每个时间步 $\Delta t$（声学 CFL：$\Delta t=\mathrm{CFL}\,\Delta x/\max(|u|+c_f)$）：

1. **双曲步**：在每个面上解 HLLC Riemann 问题（A4），用 Godunov 型格式更新全部 7 个量（A5）。输出记上标 $0$。
2. **力学松弛**（$\mu\to\infty$，全场）：解 (A6.3) 得 $\alpha_k^*,p^*$。
3. **混合能量一致（“+1”，V10）**：第 2 步的闭合直接用守恒的混合内能 $\rho e=E-\tfrac12(\rho u)^2/\rho$，它与相内能之和的差 $\delta$ 按质量分数分给两相（A6）。SPB §3.4 原文是在 $\alpha$ 固定下用 (A2.3) 重算 $p$（重初始化）；V10 起不用，理由见 A6。
4. **热松弛**（界面格 $\epsilon_I\le\alpha_1\le1-\epsilon_I$）：$pT$ 平衡（A7.1），或有限速率 $\theta$（A7.3）。
5. **化学松弛**（界面格，且液体亚稳 $T>T_{\rm sat}(p)$）：$pTg$ 平衡（A7.2），或有限速率 $\nu$（A7.4）。

二阶：原始变量 $(\alpha_1,\rho_1,\rho_2,u,p_1,p_2)$ 的 MUSCL（minmod）+ Heun。瞬时松弛（$\mu,\theta,\nu=\infty$）是投影，幂等，每一级之后都做；有限速率松弛是对时间的积分，每个时间步只在 Heun 平均之后做一次（$\tau=\Delta t$）。若每一级之后都做，对松弛 ODE $y'=-\lambda y$，Heun 的一步是 $R\big(\tfrac12(y+R(R(y)))\big)=(1-2\lambda\Delta t)\,y+O(\Delta t^2)$（$R$ 为步长 $\Delta t$ 的后向 Euler），速率被放大约 2 倍（0D 算例实测 1.93 倍，见 C. 开发记录 V6）。

代码中 4、5 两步合成一次判断：`relax=2` 时先算该格的 $pT$ 平衡态，若其 $T>T_{\rm sat}(p)$ 就直接做 $pTg$（或有限 $\nu$）闭合，否则做 $pT$（或有限 $\theta$）。这与先 $pT$ 再 $pTg$ 等价：$pTg$ 与有限 $\nu$ 的闭合只用到守恒的 $\rho$、$\rho e$ 与 $m_2$，而 $pT$ 步不改变它们。

---

## A4. HLLC（SPB §3.1.2）

面两侧状态 $L,R$；混合量 $\rho,u,\bar p$，冻结声速 $c_f$。

$$
S_L=\min(u_L-c_L,u_R-c_R),\quad S_R=\max(u_L+c_L,u_R+c_R),\quad
S_M=\frac{\bar p_R-\bar p_L+\rho_Lu_L(S_L-u_L)-\rho_Ru_R(S_R-u_R)}{\rho_L(S_L-u_L)-\rho_R(S_R-u_R)} .\tag{A4.1}
$$

$K\in\{L,R\}$ 的星态：

$$
m_{k,K}^*=m_{k,K}\frac{S_K-u_K}{S_K-S_M},\quad (\rho u)^*_K=\rho^*_KS_M,\quad
\bar p^*=\bar p_K+\rho_K(u_K-S_K)(u_K-S_M),\quad
E^*_K=\frac{E_K(u_K-S_K)+\bar p_Ku_K-\bar p^*S_M}{S_M-S_K}.\tag{A4.2}
$$

$\alpha_K^*=\alpha_K$（声波不改变体积分数）；相密度 $\rho^*_{k,K}=\rho_{k,K}(S_K-u_K)/(S_K-S_M)$；**相内能沿各相 Hugoniot (II.8)** $e_k^*-e_k+\tfrac12(p_k^*+p_k)(v_k^*-v_k)=0$，SG 闭式：

$$
p_k^*+\pi_k=(p_k+\pi_k)\frac{(\gamma_k+1)\rho_k^*-(\gamma_k-1)\rho_k}{(\gamma_k+1)\rho_k-(\gamma_k-1)\rho_k^*}\quad(\text{NASG：把 }\rho\text{ 换成 }1/(v-b)).\tag{A4.3}
$$

面通量：$F=F_K+S_K(U_K^*-U_K)$（亚声速）或 $F_K$（超声速），给出守恒通量 $(m_ku)^*,(\rho u^2+\bar p)^*,((E+\bar p)u)^*$，以及 $u^*=S_M$、$(u\alpha_1)^*$、$(\varepsilon_ku)^*$。

---

## A5. Godunov 型更新（SPB §3.2.1）

$$
\begin{aligned}
U_i^0&=U_i^n-r\,\delta_i[F^*],\qquad U=(m_1,m_2,\rho u,E),\\
\alpha_{1,i}^0&=\alpha_{1,i}^n-r\Big(\delta_i[(u\alpha_1)^*]-\alpha_{1,i}^n\,\delta_i[u^*]\Big),\\
\varepsilon_{k,i}^0&=\varepsilon_{k,i}^n-r\Big(\delta_i[(\varepsilon_ku)^*]+(\alpha_kp_k)_i^n\,\delta_i[u^*]\Big).
\end{aligned}\tag{A5.1}
$$

**均匀 $(p,u)$ 保持（SPB §3.2.2）：** $u$、$p$ 均匀时 $u^*=u$，$\delta[u^*]=0$；$\alpha^*,\varepsilon_k^*$ 取同一侧（按 $S_M$ 迎风），由 (E1) $\varepsilon_k=c_k(p)\alpha_k+q_km_k$，所以 $\varepsilon_k^0=c_k(p)\alpha_k^0+q_km_k^0$，压力不变。**这要求 $\alpha$ 与 $\varepsilon_k$ 用同一组面值和同一组 $u^*$。**IMEX 必须继承这一点（B9）。

---

## A6. 力学松弛与重初始化（SPB §3.3–3.4）

松弛 ODE（$\mu\to\infty$）里 $m_k,\rho u,E$ 不变，由相能量方程与体积方程得 $\partial_te_k+p_I\partial_tv_k=0$；积分：

$$
e_k-e_k^0+\hat p_I\,(v_k-v_k^0)=0 .\tag{A6.1}
$$

混合能量守恒要求两相用同一个 $\hat p_I$。乘 $m_k$，用 $m_kv_k=\alpha_k$：

$$
\boxed{\varepsilon_k-\varepsilon_k^0=-\hat p_I\,(\alpha_k-\alpha_k^0)}\qquad(\text{离散 Gibbs 关系 }de+p\,dv=0)\tag{A6.2}
$$

取 SPB 允许的 $\hat p_I=p$（松弛后压力），代入 (A2.2)：

$$
\alpha_k(p)=\frac{(\gamma_k-1)\big[\varepsilon_k^0-m_kq_k+p\,\alpha_k^0\big]+(p+\gamma_k\pi_k)m_kb_k}{\gamma_k(p+\pi_k)},\qquad
\sum_k\alpha_k(p)=1 .\tag{A6.3}
$$

（取 $\hat p_I=p_I^0$ 时为 SPB (III.4)：$v_k-b_k=(v_k^0-b_k)\dfrac{p_k^0+\gamma_k\pi_k+(\gamma_k-1)\hat p_I}{p+\gamma_k\pi_k+(\gamma_k-1)\hat p_I}$。）$\sum_k\alpha_k(p)$ 对 $p$ 严格递减（(E5)，B6.1 证明），根唯一，带区间的 Newton 解。

**混合能量一致（V10 起代替 SPB §3.4 的重初始化）。** 双曲步后守恒的 $\rho e=E-\tfrac12(\rho u)^2/\rho$ 与 $\sum_k\varepsilon^0_k$ 差 $\delta$（主要是数值耗散掉的动能）。把 $\delta$ 按质量分数 $Y_k$ 加进 (A6.2)：

$$
\varepsilon_k-\varepsilon_k^0=-p\,(\alpha_k-\alpha_k^0)+Y_k\delta ,\qquad \delta=\rho e-\sum_j\varepsilon^0_j ,\tag{A6.4}
$$

即 (A6.3) 中 $\varepsilon^0_k$ 换成 $\varepsilon^0_k+Y_k\delta$。由 $\sum\alpha_k=\sum\alpha_k^0=1$，$\sum_k\varepsilon_k=\rho e$。这就是 B5 的统一闭合在 $s=0$ 时的特例（B5.2、B6.1）。

- SPB §3.4 的重初始化（$\alpha$ 固定，用 (A2.3) 重算 $p$，再按 (A2.2) 重置 $\varepsilon_k$）等价于按 $\partial\varepsilon_k/\partial p|_\alpha\propto\alpha_k/(\gamma_k-1)$ 分配 $\delta$：耗散的动能主要进了可压缩性最大、质量极小的气或汽，微量相被加热到 $10^5$–$10^6$ K（C. 开发记录 V10）。
- 在用相总能写的模型里（P-S 2014 §2.2）第 $k$ 相的动能是 $\tfrac12\alpha_k\rho_k|\mathbf u|^2$，所以按质量分配。这也满足 P-S 2014 §3.1 的条件 (ii)。
- 效果（$N=1000$）：SPB 4.3.1 显式 $L_1(u)$ 2.91 → 1.82；P-S ±500 m/s、只 $p$ 松弛的起步尖峰消失（P-S 图 2 没有）；SPB 4.5 中心压力约 2200 Pa → 小于 1 Pa（SPB §4.5 说中心“压力约为 0”，5 方程精确解是真空）。
- 灵感来源：P-S 2014 §2.2、§3.1。

---

## A7. 热、化学松弛（P-S 2014 §4、Pelanti 2022 §7）

所有松弛步中 $\rho,\rho u,E$ 不变，所以混合内能 $\rho e$ 不变。

### A7.1 $pT$ 平衡（$\mu,\theta\to\infty$）

未知 $p,T$，$m_k$ 不变：

$$
\sum_km_k\,e_k(p,T)=\rho e,\qquad \sum_km_k\,v_k(p,T)=1 .\tag{A7.1}
$$

由 (E3) 第一式给出 $T(p)=\dfrac{\rho e-\sum m_kq_k}{\sum_km_kc_{vk}\frac{p+\gamma_k\pi_k}{p+\pi_k}}$；代入第二式得关于 $p$ 的标量方程（SG 下是 P-S (C.3) 的二次式）。然后 $\alpha_k=m_kv_k(p,T)$。

### A7.2 $pTg$ 平衡（$\mu,\theta,\nu\to\infty$）

未知 $p,T,Y_2$，$\rho$ 不变：$T=T_{\rm sat}(p)$，

$$
Y_2=\frac{e-e_1(p,T_{\rm sat})}{e_2(p,T_{\rm sat})-e_1(p,T_{\rm sat})},\qquad
F(p)=\rho\big[(1-Y_2)v_1+Y_2v_2\big](p,T_{\rm sat})-1=0 .\tag{A7.2}
$$

$Y_2\notin(0,1)$ 时（完全蒸发/凝结）$Y_2$ 钳在界上，$T$ 改由 (A7.1) 的能量式定。在界上 $T=T_{\rm sat}$，所以 $F$ 连续（handoff §11.3 的延拓）。

### A7.3 有限速率热松弛 $\theta$（$\mu=\infty$）

Pelanti (2022) §7.2.2 在 5 方程模型的 ODE 上冻结系数，得指数解 $\Delta T^{**}=\Delta T^*e^{-K_T\Delta t}$，再用混合能量关系 (70) 重算 $p$。

本文的统一闭合用后向 Euler（B6.4 节，取 $s=0$）：在等压下相 $k$ 吸热 $\Delta Q$，焓增 $\Delta Q$，有 $T_k=T_k^0\pm\Delta Q/C_{pk}$，$C_{pk}=\gamma_kc_{vk}m_k$，于是

$$
\Delta Q=\frac{\tau\theta\,(T_2^0-T_1^0)}{1+\tau\theta\,(C_{p1}^{-1}+C_{p2}^{-1})} .
$$

### A7.4 有限速率化学松弛 $\nu$（$\mu,\theta=\infty$）

Pelanti §7.2.3 用 4 方程模型的 ODE（71），冻结系数，得 $\Delta g^\otimes=\Delta g^{**}e^{-K_g\Delta t}$。

统一闭合（B6.5，$s=0$）：给定 $p$，$T(p,Y)$ 由能量显式定，$Y$ 由单调方程 $\rho Y_2-m_2^0+\tau\nu(g_2-g_1)(p,T)=0$ 定，再由体积约束定 $p$。

**激活判据（按原文，不自加）：** 热、化学松弛只在 $\epsilon_I\le\alpha_1\le1-\epsilon_I$ 的格子；化学松弛还要 $T>T_{\rm sat}(p)$（P-S (10)；Pelanti 8.2.1）。P-S 的 “p–pTg” 档跳过 3a：$\epsilon_I\le\alpha_1\le1-\epsilon_I$ 且任一相 $T_k>T_{\rm sat}(p)$ 时直接做 $pTg$。

---

## A8. 显式参考格式与 IMEX 的关系

显式格式的松弛步就是 B5 节统一闭合在“**没有体积应变**”（$s=0$）、输运态取双曲步结果时的特例：A6 是 $\mu=\infty$ 的力学闭合，A7.1/A7.2 是 $\theta,\nu=\infty$，A7.3/A7.4 是有限速率。

SPB 原文把“松弛 → 用总能重初始化”分两步做，相当于按 $\lambda_k\propto\alpha_k/(\gamma_k-1)$ 分配 $\delta$。V10 起显式程序 `src/SPB6_explicit.cpp` 与 IMEX 都用同一个闭合（B5.3）、按质量分数分配 $\delta$（A6）；显式程序的力学松弛就是 $s=0$ 的 $p$ 闭合（B6.1）。
