# B. TV 分裂 IMEX：从要求到离散

目标：$\Delta t$ 只受物质速度限制，$\Delta t=\mathrm{CFL}\,\Delta x/\max|u|$。

- 声波（压力）部分隐式。
- 三种松弛也在同一个隐式闭合里完成。
- 速率 $\mu,\theta,\nu$ 是闭合的参数：
  - 速率 $\to\infty$ 时严格收敛到平衡闭合（渐近保持，AP）；
  - 速率 $\to0$ 时退化为冻结格式。

本章记号同 A 章。一个 IMEX-RK 级内的隐式步长记 $\tau=a_{ii}\Delta t$，$r=\tau/\Delta x$。

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

## B3. 显式算子（BB 2025 的耗散方式）

在显式级态 $Q^e$ 上，对 $X\in\{\alpha_1,\alpha_2,m_1,m_2,\varepsilon_1,\varepsilon_2,\rho u,E\}$：

$$
F^e_{i+1/2}[X]=\tfrac12\big(X_iu_i+X_{i+1}u_{i+1}\big)^e-\tfrac12a_\Omega\big(X^{R}_{i+1/2}-X^{L}_{i+1/2}\big)^e,\qquad a_\Omega=\max_\Omega|u^e| .\tag{B3.1}
$$

- **中心部分**用格均值；**耗散**用重构的面值（BB (33)–(35)：minmod，中心部分不动）。
- **重构在被输运的守恒量上做，各量共用一个限制因子**。记 $\mathcal X=(\alpha_1,m_1,m_2,\varepsilon_1,\varepsilon_2,\rho u,E)$，格 $i$ 的斜率为
  $$\sigma_{X,i}=\phi_i\,\tfrac12(X_{i+1}-X_{i-1}),\qquad \phi_i=\min_{X\in\mathcal X}\frac{\mathrm{minmod}(X_i-X_{i-1},X_{i+1}-X_i)}{\tfrac12(X_{i+1}-X_{i-1})}\in[0,1].$$
  - 对每个量，$|\sigma_X|$ 都不超过它自己的 minmod 斜率且同号，所以仍是 TVD。
  - 所有量的斜率是**同一个**系数乘中心差分，所以量之间的线性关系在面上原样成立。均匀 $p,u$ 时 $\varepsilon_k=c_k(p)\alpha_k+q_km_k$、$\rho u=u\rho$、$E=\sum\varepsilon_k+\tfrac12u^2\rho$ 在面值上同样成立（(E1)），这是 WB 所需的。
  - 重构不经过 EOS，所以 $Q^e$ 只需 $\rho>0$；这一点对 IMEX-RK 的外推级态很重要（B10）。
- $\alpha_2$ 的面值取 $1-\alpha_1$，所以 $\sum_kF^e[\alpha_k]=\tfrac12(u_i+u_{i+1})^e\equiv\bar u^e_f$，耗散两相抵消。

**与 BB 的对应：** BB (36)–(38) 要求，对 $\phi$ 加的人工耗散必须在内能方程里按 $\partial(\rho\phi e)/\partial\phi$ 缩放后同样加上，否则接触处出现压力振荡。本格式对 $\alpha_k$ 与 $\varepsilon_k$ 用同一个 $a_\Omega$ 和同一种重构；均匀 $p$ 时由 (E1) 有 $\Delta\varepsilon_k=\frac{p+\gamma\pi}{\gamma-1}\Delta\alpha_k+q\Delta m_k$，这正是 BB (38) 的缩放，而且对任意 $\alpha,m$ 分布都自动满足。

BB (39) 的压力跳跃耗散（高马赫）作为可选项，按 BB 原式显式加在相能量上，默认关。

显式更新：$X^*_i=X^s_i-r\,\delta_i[F^e[X]]$，$Q^s$ 是级的更新基底。

---

## B4. 隐式部分：面速度、增量通量、输运态

$\rho^*_i=\sum_km^*_{k,i}$，$u^*_i=(\rho u)^*_i/\rho^*_i$。未知量是各格的混合压力 $P_i$。

**面速度**（面上的动量方程，紧致模板）：

$$
u_f=\bar u^*_f-r\,\frac{P_{i+1}-P_i}{\bar\rho_f},\qquad \bar u^*_f=\tfrac12(u^*_i+u^*_{i+1}),\quad \bar\rho_f=\tfrac12(\rho^*_i+\rho^*_{i+1}).\tag{B4.1}
$$

**增量**（参考速度就是显式中心通量里用的那个面速度）：

$$
\Delta u_f=u_f-\bar u^e_f,\qquad \bar u^e_f=\tfrac12(u^e_i+u^e_{i+1}).\tag{B4.2}
$$

**增量通量：** $\hat X_f\Delta u_f$。$\hat X_f$ 是供体格（按 $\Delta u_f$ 的符号取上游格 $d$）**单位体积**的内容：

$$
\hat X_f=\frac{X^*_d}{\sum_j\alpha^*_{j,d}},\qquad X\in\{\alpha_k,m_k,\varepsilon_k\}.
$$

- 显式输运之后，供体格内的物质占体积 $\sum_j\alpha^*_{j,d}=1-r\delta_d[\bar u^e]\neq1$。增量通量搬走的是体积 $r\Delta u_f$ 的这种物质，所以要按单位体积归一。
- 于是 $\sum_k\hat\alpha_{k,f}=1$，(B4.4) 严格成立。
- 三个量用同一个归一因子，比值 $\rho_k,e_k$ 不变（R3）。
- 增量 $\Delta u_f=O(\tau)$，所以供体格的一阶面值在光滑区只带来 $O(\Delta x\,\tau)$ 的误差。

输运态为

$$
X^{\rm tr}_i=X^*_i-r\,\delta_i[\hat X_f\Delta u_f],\qquad X\in\{\alpha_k,m_k,\varepsilon_k\}.\tag{B4.3}
$$

**体积应变**（显式中心通量与增量合起来恰好是 $u_f$；散度对压力完全隐式）：

$$
s_i\equiv1-\sum_k\alpha^{\rm tr}_{k,i}=r\,\delta_i[\bar u^e_f+\Delta u_f]=r\,\delta_i[u_f].\tag{B4.4}
$$

**动量与总能**（守恒，R1）。面压力取 $p_f=\tfrac12(P_i+P_{i+1})$，$\hat M_f=\sum_k\hat m_{k,f}\Delta u_f$：

$$
\begin{aligned}
(\rho u)^{n+1}_i&=(\rho u)^*_i-r\,\delta_i[p_f]-r\,\delta_i[u_f\hat M_f],\\
E^{n+1}_i&=E^*_i-r\,\delta_i\Big[\textstyle\sum_k\hat\varepsilon_{k,f}\Delta u_f+\tfrac12u_f^2\hat M_f+p_fu_f\Big].
\end{aligned}\tag{B4.5}
$$

- 增量质量携带新时间层面速度 $u_f$：在均匀加速度场中 $u$ 保持均匀。
- 动量压力项 $-r\delta[p_f]=-\tfrac r2(P_{i+1}-P_{i-1})$。
- 守恒的混合内能：

$$
\rho^{n+1}_i=\sum_km^{\rm tr}_{k,i},\qquad (\rho e)^C_i=E^{n+1}_i-\frac{\big((\rho u)^{n+1}_i\big)^2}{2\rho^{n+1}_i}.\tag{B4.6}
$$

$X^{\rm tr}_i$、$s_i$、$(\rho e)^C_i$ 都只依赖 $(P_{i-1},P_i,P_{i+1})$。

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

记 $a_k=\alpha^{\rm tr}_k+\omega_ks=\alpha^{\rm tr}_k/(1-s)$，则 $\alpha_1=a_1+R$，$\alpha_2=a_2-R$；$p_I=w_1p_1+w_2p_2$，$w_1=Z_2/(Z_1+Z_2)$（$Z$ 冻结在输运态）。给定 $R$ 时，(B5.1)–(B5.3) 对 $(p_1,p_2)$ 是 2×2 线性方程组（$k=1,2$）：

$$
\Big(\frac{\alpha_k-m_kb_k}{\gamma_k-1}+\omega_ks\Big)p_k+\chi_kR\,(w_1p_1+w_2p_2)-\lambda_ks\,(\omega_1p_1+\omega_2p_2)
=\varepsilon^{\rm tr}_k-m_kq_k-\frac{\gamma_k\pi_k(\alpha_k-m_kb_k)}{\gamma_k-1}+\lambda_kD_0 .\tag{B6.3}
$$

再解标量方程

$$
\mathcal F(R)=R-\tau\mu\big(p_1(R)-p_2(R)\big)=0,\qquad R\in\big(-(a_1-m_1b_1),\ a_2-m_2b_2\big).
$$

$R$ 增大时 $\alpha_1$ 变大、$p_1$ 下降、$p_2$ 上升，所以 $\mathcal F$ 递增，根唯一。输出 $\bar p^{\rm cl}=\alpha_1p_1+\alpha_2p_2$。

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

- **单调：** 能量固定时 $Y\uparrow$ 使 $T\downarrow$（$e_2(p,T)>e_1(p,T)$）；$\partial_T(g_2-g_1)=-(s_2-s_1)<0$。所以 $\mathcal R'(Y)=\rho+\tau\nu(s_2-s_1)|\partial_YT|>0$，$Y$ 在 $[Y_{\min},1-Y_{\min}]$ 上唯一。
- 若 $\mathcal R(Y_{\min})>0$ 或 $\mathcal R(1-Y_{\min})<0$，$Y$ 取界值：完全凝结/蒸发。这是速率方程的解落在可容许域边界，不是额外规则。
- 然后解 $\Phi(p)=\rho\sum_kY_kv_k(p,T)-1=0$。
- 能量方程只用混合量，$\rho_I,e_I$ 不出现。

### B6.6 $\mu=\theta=\nu=\infty$：$pTg$ 闭合

$T=T_{\rm sat}(p)$，$Y$ 由能量定（A7.2），$\Phi(p)=\rho\sum_kY_kv_k-1$。$Y$ 出界时钳住，$T$ 由能量定，与 B6.5 的界值分支相同，所以 $\Phi$ 连续。

### B6.7 激活（按论文）

$\theta,\nu$ 只在 $\epsilon_I\le\alpha_1\le1-\epsilon_I$ 的格子非零；$\nu$ 还要 $T>T_{\rm sat}(p)$（P-S (10)、Pelanti 8.2.1）。

判据在每一级开始时用**显式预测态**（$\Delta u_f=0$ 时的输运态 $X^*$ 的 $pT$ 平衡）计算，并在该级的 Newton 中冻结，使残量对 $P$ 连续。这与 P-S 把判据作用在双曲步之后的状态上是同一个意思。

---

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

- $\partial s_i/\partial P_i=r^2(1/\bar\rho_{i+1/2}+1/\bar\rho_{i-1/2})>0$，且 $\partial\bar p^{\rm cl}/\partial s<0$：膨胀降压（E5）。
- $\partial E_i/\partial P_i\approx-r^2\,2\bar p/\bar\rho<0$。

所以 $\partial G_i/\partial P_i<-1$，相邻项 $O(r^2\rho c^2/\bar\rho)$，矩阵是 $I-\tau^2\partial_x(c^2\partial_x)$ 型：对角占优的 M 矩阵。

**求解：**

- Newton。Jacobian 用逐格差分：扰动三个（五个）压力，各算一次局部闭合。
- 一维用 Thomas；二维用 BiCGStab + Jacobi 预条件（`2DEuler_TV` 风格，但用相对容差）。
- 线搜索保证可容许性：$X^{\rm tr}$ 为正、局部闭合有解。
- 收敛判据：$\max_i|G_i|/(|P_i|+p_{\rm ref})<10^{-10}$。

**步长控制（不是补丁，而是可容许性条件）：** 以下任一情况发生时，拒绝本步，$\Delta t\leftarrow\Delta t/2$ 重算：
- $Q^e$ 出现 $\rho\le0$（显式通量只需要这一条，B3）；
- 闭合输入不可容许（$\alpha^{\rm tr}_k\le m^{\rm tr}_kb_k$、$m^{\rm tr}_k\le0$ 或 $N_k\le0$）；
- Newton 不收敛。

这和显式格式的正性 CFL 条件是同一类东西。所有拒绝都计数并报告。

---

## B9. 性质证明

### B9.1 均匀 $(p,u)$ 保持（R2，创新点 1）

设 $u\equiv u_0$、$p_k\equiv p$，$\alpha_k,\rho_k$ 任意（各相温度也可以任意）。

1. **显式：** 由 (E1)，格均值满足 $\varepsilon_k=c_k(p)\alpha_k+q_km_k$、$\rho u=u_0\rho$；共用限制因子的重构保持这些线性关系，所以面值也满足。(B3.1) 对 $X$ 线性，所以 $\varepsilon^*_k=c_k\alpha^*_k+q_km^*_k$，$(\rho u)^*=u_0\rho^*$，$E^*=\sum\varepsilon^*_k+\tfrac12\rho^*u_0^2$。
2. **隐式：** 取 $P\equiv p$。$\bar u^*_f=u_0=\bar u^e_f$，所以 $\Delta u_f=0$，$u_f=u_0$，$s=0$，$X^{\rm tr}=X^*$。$(\rho u)^{n+1}=u_0\rho^*$，$E^{n+1}=E^*-r\delta[pu_0]=E^*$，所以 $(\rho e)^C=\sum\varepsilon^*_k$，$D_0=0$，$\delta=0$。
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

默认表由数值试验决定，两者都报告。

第 $i$ 级：

1. $Q^e=U^n+\Delta t\sum_{j<i}\tilde a_{ij}K_j$，$Q^s=U^n+\Delta t\sum_{j<i}a_{ij}K_j$。
2. $U_i=\Phi(Q^e,Q^s;\tau=a_{ii}\Delta t)$：B3 的显式更新，加 B4–B8 的隐式闭合。
3. $K_i=(U_i-Q^s)/\tau$。

隐式表刚性精确，所以 $U^{n+1}$ 就是末级 $U_s$。$K_i$ 包含松弛源项的贡献，松弛是隐式表的一部分（BB §3.5 / Pelanti 的源项在这里由同一个闭合给出）。

**起步：** 从静止或间断初值起算时，前若干步用声学步长且用一阶（单级后向 Euler），之后每步至多放大 1.2 倍，直到物质 CFL（`2DEuler_TV` 的 `full_step` 与 BB 的做法）。

---

## B11. 二维

- 显式：$x,y$ 两个方向的面各算一次 (B3.1)。
- 面速度 (B4.1) 只修正法向分量；增量质量通量携带新的法向面速度 $u_f$，以及按供体取的切向速度 $\hat v$。
- $s_{ij}=r_x\delta_x[u_f]+r_y\delta_y[v_f]$。
- 动量：$x$ 方向 $-r_x\delta_x[p_f]-r_x\delta_x[u_f\hat M]-r_y\delta_y[\hat u\hat M]$，$y$ 方向对称。
- 总能的增量动能取 $\tfrac12(u_f^2+\hat v^2)\hat M$。
- 全局 Jacobian 5 对角。

---

## B12. 与 handoff v23 的差别（换路径的地方）

| | handoff v23 | 本推导 |
|---|---|---|
| 混合能量一致性 | 闭合用非守恒相能量，级末再用总能重初始化 | $\delta$ 进入闭合（B5.2），没有级末投影；动量方程看到的就是守恒能量下的压力 |
| 松弛 | 瞬时（$p$、$pT$、$pTg$）在闭合里；有限速率没有 | 一个闭合，$\mu,\theta,\nu$ 是参数（B6），两端极限严格（B7） |
| 全局残量 | $\sum\alpha-1$ 或能量残量，按闭合类型不同 | $G=\bar p^{\rm cl}-P$，所有闭合同一形式 |
| 传质判据 | 上一级的标记（滞后） | 本级显式预测态，级内冻结 |
| 显式耗散 | 中心部分也用重构值 | BB 的做法：中心部分用格均值，耗散用重构值 |
| 不可容许时 | 逐格限制 $Q^e$、退回一阶、压力下限 | 拒绝本步并减半 $\Delta t$，计数报告 |
| 代码结构 | 一维专用 | `2DEuler_TV` 结构（`Grid`、`Qe/Qs/Hs`、二维） |
