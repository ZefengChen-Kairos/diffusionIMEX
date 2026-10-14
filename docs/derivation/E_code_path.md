# E. 唯一代码路径（最终版，V16）：从输入到输出，按本页即可复现

本章把最终固定下来的**唯一**代码路径从头到尾写一遍，每一步给出**全部公式与常数**，并注明代码位置与推导节号（A、B 章是理由，本章是做法）。程序里没有对照开关：所有算例走同一条路径，只有算例自带的物理参数、区域、边界、网格、终止时间、松弛档位不同（E9）。

- IMEX：`src/SPB6_TV_IMEX.cpp`（推导 B 章）。显式参考：`src/SPB6_explicit.cpp`（推导 A 章，E6）。两者共用 `src/spb6_physics.h`（状态方程、闭合、激波预算、算例、边界条件）。
- 两个程序都直接由 `reference/code/2DEuler_TV.cpp` 改来（函数名、`Grid`/`Array`、`Qe/Qs/Hs` 的 IMEX-RK 循环、VTK 输出保留）。
- 每处选择的来由、试过而放弃的做法见 C. 开发记录（`docs/devlog.md`）。

记号：格 $i$（二维 $(i,j)$），面 $f=i+\tfrac12$；$h$ 为该方向的网格宽（$\Delta x$ 或 $\Delta y$）；$r=\tau/h$；$\delta_i[F_f]=F_{i+1/2}-F_{i-1/2}$；$(\cdot)^\pm=\max/\min(\cdot,0)$；$\chi_1=+1,\chi_2=-1$。

---

## E0. 固定的设置

| 项 | IMEX | 显式参考 |
|---|---|---|
| 时间积分 | LSDIRK2(2,2,2)，$\gamma=1-1/\sqrt2$（`rk=2`）；起步斜坡内的步与 MOOD 回退的步用单级后向 Euler | Heun（二阶 SSP-RK） |
| 时间步 | 物质 CFL $C$：$\Delta t=C\big/\big(\max_i\lvert u_i\rvert/\Delta x+\max_i\lvert v_i\rvert/\Delta y\big)$，$C=0.2$；十二烷 SG（算例 3）与 SPB 5.1（算例 12）$C=0.05$ | 声学 CFL 0.4：$\Delta t=0.4\big/\big(\max(\lvert u\rvert+c_f)/\Delta x+\dots\big)$ |
| 起步 | 首步 $\Delta t=0.4\,\Delta t_{\rm ac}/1.2$ 再放大 1.2（即声学 CFL 0.4），之后每步至多放大 1.2 倍；被限制的步用后向 Euler | — |
| 空间 | 二阶 MUSCL：体积输运量共用 minmod 限制因子，速度单独 minmod（E3.1） | 原始变量 minmod（SPB 附录 C） |
| 对流通量 | Rusanov，耗散速度为物质速度：逐面 $\max(\lvert u^L_f\rvert,\lvert u^R_f\rvert)$，或全场最大（算例 13、71、72）（E3.2） | HLLC（A4） |
| 压力/声学 | 隐式（E4），面上加按 $(1-c\tau/h)^+$ 加权的声学 Godunov 耗散（E4.2） | 显式（HLLC） |
| 松弛 | 统一闭合，在隐式压力方程里（E5） | 每级之后的瞬时投影 + 每步一次有限速率（A3） |
| 能量亏损分配 | 物理耗散按 SPB (II.5) 的相 Hugoniot 能隙份额，误差按可容许能量份额（E5.0，V16） | 同一规则（E6） |
| Newton | 残量容差 $10^{-13}$（相对 $\lvert P\rvert+K+P_E+1$），至多 60 次；步长判据同容差 | — |
| 线性解法 | 一维非周期 Thomas；其余 BiCGStab + 两级预条件，相对容差 $10^{-13}$，至多 4000 次 | — |
| 界面判据 $\epsilon_I$ | $10^{-4}$（热、化学松弛只在 $\epsilon_I\le\alpha_1\le1-\epsilon_I$ 的格） | 同左 |
| 质量分数下限 $Y_{\min}$ | $10^{-10}$（$pTg$、有限 $\nu$ 的钳位） | 同左 |

命令行 `key=value`：`case= N= ny= cfl= rus= relax= mu= theta= nu= epsI= rk= cflac= ramp= tol= dt= limiter= order= tend= out= snap= probe= prec= p0= wellprep= vtk= verbose=`。收敛阶试验另加 `limiter=0`（不限制的中心斜率）。

---

## E1. 变量与状态方程（A1、A2）

**守恒变量**（每格 8 个，次序 `IA1 IM1 IM2 IMU IMV IE1 IE2 IEN`）：$\alpha_1,\ m_1=\alpha_1\rho_1,\ m_2=\alpha_2\rho_2,\ \rho u,\ \rho v,\ \varepsilon_1=\alpha_1\rho_1e_1,\ \varepsilon_2,\ E$。$\alpha_2=1-\alpha_1$，$\rho=m_1+m_2$，$E=\varepsilon_1+\varepsilon_2+\tfrac12\rho\lvert\mathbf u\rvert^2$ 守恒；相内能 $\varepsilon_k$ 非守恒（模型的 $\alpha_kp_k\nabla\cdot\mathbf u$ 项）。

**NASG / SG**（`Phase{g, pi, q, qp, cv, b}`，`kind = 0`；$b=0$ 为 SG）。以单位体积内容 $(\alpha,m,\varepsilon)$ 写：

$$
\varepsilon_k=\frac{(p+\gamma_k\pi_k)(\alpha_k-m_kb_k)}{\gamma_k-1}+m_kq_k,\qquad
p_k=\frac{(\varepsilon_k-m_kq_k)(\gamma_k-1)}{\alpha_k-m_kb_k}-\gamma_k\pi_k,\qquad
T_k=\frac{(\alpha_k/m_k-b_k)(p_k+\pi_k)}{(\gamma_k-1)c_{vk}},
$$

$$
v_k(p,T)=\frac{(\gamma_k-1)c_{vk}T}{p+\pi_k}+b_k,\quad e_k(p,T)=c_{vk}T\frac{p+\gamma_k\pi_k}{p+\pi_k}+q_k,\quad
s_k=c_{vk}\big[\gamma_k\ln T-(\gamma_k-1)\ln(p+\pi_k)\big]+q'_k,\quad g_k=\gamma_kc_{vk}T+b_kp+q_k-Ts_k,
$$

$$
c_k^2=\frac{\gamma_k(p+\pi_k)}{\rho_k(1-\rho_kb_k)},\qquad \partial\varepsilon_k/\partial p|_{\alpha,m}=\frac{\alpha_k-m_kb_k}{\gamma_k-1}.
$$

可容许：$p_k+\pi_k>0$，$\alpha_k>m_kb_k$。

**Mie–Grüneisen**（`kind = 1` Cochran–Chan，`kind = 2` JWL；参数 `G, rr, c1, c2, e1, e2, e0`）：$p=\Gamma\rho(e-\varepsilon_r(\rho))+p_r(\rho)$，$\Gamma$ 常数，$\varepsilon_r'(\rho)=p_r/\rho^2$。记 $x=\rho/\rho_r$：

- CC：$p_r=c_1x^{e_1}-c_2x^{e_2}$，$p_r'=(e_1c_1x^{e_1}-e_2c_2x^{e_2})/\rho$，$\varepsilon_r=\big(c_1x^{e_1}/(e_1-1)-c_2x^{e_2}/(e_2-1)\big)/\rho$；
- JWL：$t_i=c_ie^{-e_i/x}$，$p_r=t_1+t_2$，$p_r'=(e_1t_1+e_2t_2)\rho_r/\rho^2$，$\varepsilon_r=t_1/(e_1\rho_r)+t_2/(e_2\rho_r)-e_0$。

$\varepsilon_k=m_k\varepsilon_r(\rho_k)+\alpha_k(p-p_r)/\Gamma_k$，$p_k=\Gamma_k(\varepsilon_k-m_k\varepsilon_r)/\alpha_k+p_r$，$c_k^2=(\Gamma_k+1)(p-p_r)/\rho_k+p_r'$，$T_k=\alpha_k(p-p_r)/(\Gamma_km_kc_{vk})$，$\partial\varepsilon_k/\partial p|_{\alpha,m}=\alpha_k/\Gamma_k$。可容许：$c_k^2>0$。MG 相只做力学松弛（$p$ 闭合）。

**混合量。** 混合压力（给定 $\alpha_k,m_k$ 与混合内能 $\rho e$，A2.3/A2.7；对 $p$ 线性）：

$$
p_{\rm mix}=\frac{\rho e-\sum_kN_k}{\sum_k\partial\varepsilon_k/\partial p|_{\alpha,m}},\qquad
N_k=\begin{cases}m_kq_k+\gamma_k\pi_k\dfrac{\alpha_k-m_kb_k}{\gamma_k-1}&\text{NASG}\\[2pt] m_k\varepsilon_r-\alpha_kp_r/\Gamma_k&\text{MG}\end{cases}
$$

冻结声速 $c_f^2=\sum_kY_kc_k^2$（$Y_k=m_k/\rho$）；Wood 体积模量 $K=\big[\sum_k\alpha_k/(\rho_kc_k^2)\big]^{-1}$；$\bar p=\sum_k\alpha_kp_k$（`pbar_of`）。

**饱和线 $T_{\rm sat}(p)$**（`build_sat`，只在 `relax>0` 且无 MG 相时建）：解 $g_1(p,T)=g_2(p,T)$，取 $s_2>s_1$ 的分支（$f=g_1-g_2$ 对 $T$ 递增）。在 $p\in[0.1,5\times10^9]$ Pa 的对数网格（1600 点）上逐点求根：第一点在 $[20,5000]$ K 上扫描找变号再二分，其后用上一点为初猜做 Newton（$T'=T-f/(s_2-s_1)$，步长限在 $[\tfrac12T,2T]$）；当 $v_v(p,T)\le v_l(p,T)$ 或无根时停止，端点 $p_c$ 用几何二分细化，$v_c=v_l(p_c,T_{\rm sat}(p_c))$。`Tsat(p)`：表上线性插值作初猜，再 Newton 抛光；$p\ge p_c$ 或 $p\le0$ 返回 NaN。

---

## E1.5 怎么分裂（TV 分裂，B1–B2）

模型（单速度、两压力的 6+1 方程，SPB 2009）一维写作

$$
\partial_t\alpha_1+u\partial_x\alpha_1=\mu(p_1-p_2),\quad \partial_tm_k+\partial_x(m_ku)=\chi_k\dot m,\quad
\partial_t(\rho u)+\partial_x(\rho u^2+\bar p)=0,
$$
$$
\partial_t\varepsilon_k+\partial_x(\varepsilon_ku)+\alpha_kp_k\partial_xu=-\chi_k\,p_I\mu(p_1-p_2)+\chi_k(\mathcal Q+e_I\dot m),\quad
\partial_tE+\partial_x\big((E+\bar p)u\big)=0 .
$$

仿照 Toro–Vázquez-Cendón（`2DEuler_TV.cpp`），每个被输运量的通量 $Xu^{n+1}$ 拆成 $X\,u^e$（显式级态的速度，**显式**）加 $X\,(u^{n+1}-u^e)$（**隐式增量**）：

| 方程 | 显式（对流，E3：MUSCL–Rusanov，耗散速度 = 物质速度） | 隐式（E4–E5：压力、声学、闭合，一个全局压力方程） |
|---|---|---|
| $\alpha_k$ | $\partial_x(\alpha_ku^e)$ | 增量 $\partial_x(\hat\alpha_k\Delta u)$；体积应变 $-\alpha_k\partial_xu$ 与 $\mu$ 松弛进闭合 |
| $m_k$ | $\partial_x(m_ku^e)$ | 增量 $\partial_x(\hat m_k\Delta u)$；传质 $\nu$ 进闭合 |
| $\varepsilon_k$ | $\partial_x(\varepsilon_ku^e)$ | 增量 $\partial_x(\hat\varepsilon_k\Delta u)$；压力功 $\alpha_kp_k\partial_xu$、松弛、能量亏损分配进闭合 |
| $\rho u$ | $\partial_x(\rho u\,u^e)$ | $\partial_x p_f$、$\partial_x(u_f\hat M\Delta u)$ |
| $E$ | $\partial_x(Eu^e)$ | $\partial_x(p_fu_f)$、$\partial_x(\hat E\Delta u)$ |

- 显式子系统的特征值全是 $u$（$\alpha,m,\varepsilon$ 线性输运；$(\rho,\rho u)$ 是无压气体），所以 Rusanov 的耗散速度取物质速度，时间步只受物质 CFL 限制；声波全在隐式部分。
- 体积方程的 path-conservative Rusanov（线性路径）$=$ 守恒的体积通量 $F[\alpha]$ $+$ 体积应变 $-\bar\alpha\Delta u$（B2）：前者与 $m,\varepsilon$ 同形式输运，后者移入闭合与压力功配对，成为离散 Gibbs 关系（E5）。
- 质量与内能都带隐式增量、用同一组系数，所以 $e_k=\varepsilon_k/m_k$ 与 $\rho_k$ 在输运中一致（温度一致，B9.2）。
- 隐式部分的未知量只有每格一个混合压力 $P$；面速度、面压力对 $P$ 线性（E4.2），格的输运态、动量、总能对 $P$ 分段线性（E4.3），闭合（E5）把它们映射成压力，残量 $G=\bar p^{\rm cl}-P$。所有松弛（三种、有限或无限速率）都在这个闭合里，与压力同时隐式求解。

---

## E2. 一个时间步（`main` 的时间循环）

**E2.1 初始化。** 网格 $\Delta x=(x_1-x_0)/n_x$，$\Delta y$ 同；初值按 `test_state`（E9）的原始变量 $(\alpha_1,\rho_1,\rho_2,u,v,p_1,p_2)$ 换算成守恒量（圆形区域的格取 $4\times4$ 子点平均，`TestSetup::sub`）；填虚格（E7）。压力参考值 $P_{\rm ref}=\min_i\bar p_i$（初值的混合压力最小值；B8），Newton 的未知量是 $P'_i=P_i-P_{\rm ref}$，初值 $P'_i=\bar p_i-P_{\rm ref}$。

**E2.2 步长。** $\Delta t=C\cdot\big[\max_i\lvert u_i\rvert/\Delta x+\max_i\lvert v_i\rvert/\Delta y\big]^{-1}$（格速度，不含声速）。起步斜坡：记上一步 $\Delta t_{\rm prev}$（首步前取 $0.4\,\Delta t_{\rm ac}/1.2$，$\Delta t_{\rm ac}$ 用 $\lvert u\rvert+c_f$ 算）；若 $\Delta t>1.2\,\Delta t_{\rm prev}$ 则 $\Delta t=1.2\,\Delta t_{\rm prev}$ 且本步标记为斜坡步（用后向 Euler）。再截到终止时间与下一个 `snap=` 时刻。`dt=` 给定时用固定步长（无斜坡）。$\Delta t<10^{-12}t_{\rm end}$ 时报告失败。

**E2.3 一步的尝试循环**（至多 30 次）：

1. 选表：斜坡步或 MOOD 回退步用后向 Euler（$a_{11}=1$），否则 LSDIRK2
   $$\tilde A=\begin{pmatrix}0&0\\ \tfrac1{2\gamma}&0\end{pmatrix},\quad A=\begin{pmatrix}\gamma&0\\ 1-\gamma&\gamma\end{pmatrix},\quad\gamma=1-\tfrac1{\sqrt2}$$
   （`rk=3`：SSP2(3,3,2)，$\tilde A=\big(\begin{smallmatrix}0&0&0\\ \frac12&0&0\\ \frac12&\frac12&0\end{smallmatrix}\big)$，$A=\big(\begin{smallmatrix}\frac14&0&0\\0&\frac14&0\\ \frac13&\frac13&\frac13\end{smallmatrix}\big)$）。两张表刚性精确（$b=A$ 末行），$U^{n+1}=U_s$。
2. **激活与投影**（`relax>0`，E5.7）：按 $U^n$ 逐格判断闭合类型（$p$ / 有限 $\mu$ / $pT$ / 有限 $\theta$ / 有限 $\nu$ / $pTg$），瞬时类型（$pT$、$pTg$）的格把 $U^n$ 投影到该闭合的平衡流形上（$s=0$、$\tau=0$ 的同一个闭合，守恒 $m_k$ 或 $\rho$、$\rho\mathbf u$、$E$；投影失败的格退回力学类型）。本步各级冻结这些类型（`STEP_TYPE`）。
3. **各级** $s=1,\dots$：
   1. $Q^e=U^n+\Delta t\sum_{r<s}\tilde a_{sr}H_r$，$Q^s=U^n+\Delta t\sum_{r<s}a_{sr}H_r$（含虚格），填虚格；$\tau=a_{ss}\Delta t$。
   2. 在 $Q^e$ 上重构（E3.1）。
   3. 显式部分（E3.2）：$X^*=Q^s-\tau\,\nabla\!\cdot F^e(Q^e)$，并给出 $\alpha_2^*$ 与面速度 $\bar u^e_f$。
   4. 隐式部分（E4）：由 $X^*$ 解全局压力方程，得级终态 $U_s$ 与压力场 $P'$。
   5. $H_s=(U_s-Q^s)/\tau$。
4. 任一级失败（$Q^e$ 的 $\rho\le0$、$X^*$ 某格 $\sum\alpha^*\le0$ 或 $\rho^*\le0$、闭合不可容许、Newton 不收敛）：
   - 第一次失败且不是斜坡步：**同一步长**改用后向 Euler 重算（MOOD 回退，计 `mood`）；
   - 仍失败：拒绝，$\Delta t\leftarrow\Delta t/2$（计 `rejects`），从第 1 项重来。
5. 成功：$U^{n+1}=U_s$，$P'$ 留作下一步的 Newton 初猜；$t\leftarrow t+\Delta t$。

**E2.4 输出。** `snap=` 时刻写剖面；`probe=` 每步记录最近格的压力；结束时做一次激活投影（P-S：输出松弛后的态）再写剖面与统计行（E10）。

---

## E3. 显式部分（B3，`TVDReconstruction`、`Update`）

### E3.1 重构（在 $Q^e$ 上，每个方向各一次）

对格 $i$，三点 $a=i-1$、$c=i$、$b=i+1$。

- **体积输运量** $\mathcal X_V=(\alpha_1,m_1,m_2,\varepsilon_1,\varepsilon_2)$ 共用一个限制因子：
  $$\phi_i=\min_{X\in\mathcal X_V}\frac{\mathrm{minmod}(X_c-X_a,\ X_b-X_c)}{\tfrac12(X_b-X_a)}\in[0,1],$$
  其中：某量两侧差分都满足 $\lvert\cdot\rvert\le10^{-12}S_X+10^{-300}$（模板内为常数到舍入）时该量不参与；某量两侧差分异号（极值）时 $\phi_i=0$。舍入尺度 $S_X$：$\alpha_1$ 为 3；$m_k$ 为 $\lvert X_a\rvert+\lvert X_c\rvert+\lvert X_b\rvert$；$\varepsilon_k$ 为 $\sum_{a,c,b}\big(\lvert\varepsilon_1\rvert+\lvert\varepsilon_2\rvert+\tfrac12\rho\lvert\mathbf u\rvert^2\big)$。面值 $X^{L}_i=X_i-\tfrac12\sigma_i$（格的左面）、$X^{R}_i=X_i+\tfrac12\sigma_i$（格的右面），$\sigma_i=\phi_i\cdot\tfrac12(X_b-X_a)$。
- **速度** $u,v$（$u_i=(\rho u)_i/\rho_i$）各分量单独 minmod：$\sigma^u_i=\mathrm{minmod}(u_c-u_a,u_b-u_c)$（两侧差分都 $\le10^{-12}S_u$ 时为 0，$S_u=\sum_{a,c,b}(\lvert\mathbf u\rvert+c_f)$），$u^{L,R}_i=u_i\mp\tfrac12\sigma^u_i$。面上组装 $(\rho u)^{L,R}=\rho^{L,R}u^{L,R}$，$E^{L,R}=\varepsilon_1^{L,R}+\varepsilon_2^{L,R}+\tfrac12\rho^{L,R}\lvert\mathbf u^{L,R}\rvert^2$（$\rho^{L,R}=m_1^{L,R}+m_2^{L,R}$）。
- 存 $\phi_i\leftarrow\min(\phi_i,\ \text{速度的光滑比}\ \sigma^u_i/\tfrac12(u_b-u_a))$（极值处 0），给 E4.2 的 $\psi^P$ 用。
- `limiter=0`：所有斜率取中心差分 $\tfrac12(X_b-X_a)$（只用于光滑解的收敛阶）。`order=1`：无重构。

### E3.2 Rusanov 通量与更新

面 $f=i+\tfrac12$ 的左右态取 $Q_L=Q^{R}_i$、$Q_R=Q^{L}_{i+1}$（上面的面值）。对流通量 $F(Q)=u\,(\alpha_1,m_1,m_2,\rho u,\rho v,\varepsilon_1,\varepsilon_2,E)$（**不含压力**，TV 分裂 B1）：

$$
F^e_f=\tfrac12\big(F(Q_L)+F(Q_R)\big)-\tfrac12a_f(Q_R-Q_L),\qquad
a_f=\begin{cases}\max(\lvert u_L\rvert,\lvert u_R\rvert)&\text{逐面（默认）}\\ a_\Omega=\max\big(\max_i\lvert u_i\rvert,\ \max_f\max(\lvert u_L\rvert,\lvert u_R\rvert)\big)&\text{全场（算例 13、71、72）}\end{cases}
$$

（$x$、$y$ 方向各自取 $a$；全场最大先对格速度、再对所有面态速度取最大。）更新：

$$
X^*_i=Q^s_i-\frac{\tau}{\Delta x}\big(F^e_{i+1/2}-F^e_{i-1/2}\big)-\frac{\tau}{\Delta y}\big(\cdot\big)_y ,\qquad
\bar u^e_f=\tfrac12(u_L+u_R),\qquad
\alpha_2^*=(1-\alpha_1)^s_i-\frac{\tau}{\Delta x}\delta_i\big[\bar u^e_f-F^e_f[\alpha_1]\big]-\dots
$$

（$\alpha_2$ 的面值取 $1-\alpha_1$ 的面值，所以 $F^e[\alpha_2]=\bar u^e_f-F^e[\alpha_1]$，两相的 Rusanov 耗散抵消，$\sum_kF^e[\alpha_k]=\bar u^e_f$）。$\rho^e\le0$ 的格使本级失败。之后填 $X^*$ 与 $\alpha_2^*$ 的虚格。

---

## E4. 隐式部分（B4–B8，`PressureSubsystem`）

### E4.1 级数据（由 $X^*$，`StageData`，一级内冻结）

- 单位体积内容 $\hat c_i=(\hat\alpha_1,\hat\alpha_2,\hat m_1,\hat m_2,\hat\varepsilon_1,\hat\varepsilon_2)=(\alpha_1^*,\alpha_2^*,m_1^*,m_2^*,\varepsilon_1^*,\varepsilon_2^*)/(\alpha_1^*+\alpha_2^*)$；供体速度 $u^*_i=(\rho u)^*_i/\rho^*_i$，$v^*_i$ 同（$\rho^*_i=m_1^*+m_2^*$）。$\alpha_1^*+\alpha_2^*\le0$ 或 $\rho^*\le0$ 的格使本级失败。
- 冻结声速 $c_i$：两相都可容许时 $c_i^2=(m_1^*c_1^2+m_2^*c_2^2)/\rho^*$（各相用 $p_k(X^*)$），否则 $c_i=\infty$（该格相邻面不加声学项）。
- 残量尺度 $K_i+P_{E,i}$：$K$ 为 Wood 体积模量（不可容许时取 $\max_k\pi_k$），$P_E=\sum_k\big(\lvert\varepsilon^*_k\rvert+\text{eps\_round}_k\big)\big/\sum_k\partial\varepsilon_k/\partial p$，$\text{eps\_round}_k=\lvert m_kq_k\rvert$（NASG）或 $\lvert m_k\varepsilon_r\rvert+\lvert\alpha_kp_r/\Gamma_k\rvert$（MG）。
- 闭合类型：`relax>0` 时取本步冻结的 `STEP_TYPE`，否则 $p$（$\mu=\infty$）或有限 $\mu$。有限 $\mu$ 的格另存阻抗 $Z_k=\rho_kc_k$（由 $\hat c$ 算）。
- 激波预算 $Q_{k,i}$（E4.6）。

### E4.2 面量（B4.1，`set_face`）

对每个 $x$ 面 $(i+\tfrac12,j)$（$y$ 面对称，把 $u\to v$、$\Delta x\to\Delta y$）：记 $L=i$、$R=i+1$，$\rho_{L,R}=\rho^*$，$\bar\rho_f=\tfrac12(\rho_L+\rho_R)$，

$$
\boxed{\ u_f=\bar u^*_f-k_f\,(P_R-P_L),\qquad p_f=w_LP_L+(1-w_L)P_R-p^v_f\ }
$$

$$
\bar u^*_f=\tfrac12(u^*_L+u^*_R),\quad w_L=\frac{\rho_R}{\rho_L+\rho_R},\quad
k_f=\frac{\tau}{\bar\rho_f\,h}+\varphi_f\psi^P_f\frac1{Z_L+Z_R},\quad
p^v_f=\varphi_f\psi^u_f\frac{Z_LZ_R}{Z_L+Z_R}D^*_f,\quad Z=\rho^*c .
$$

- 声学权 $\varphi_f=\big(1-\tau\min(c_L,c_R)/h\big)^+$；$c_L$ 或 $c_R$ 为 $\infty$ 时 $\varphi_f=0$（$k_f$ 只剩半格动量平衡项，$p^v_f=0$）。
- 二阶因子（面的四点 $m=L-1,L,R,p=R+1$ 都在网格内时才算，否则为 1）：
  $$\psi_{\rm face}(x_m,x_L,x_R,x_p)=1-\frac{\mathrm{minmod}(x_L-x_m,D)+\mathrm{minmod}(D,x_p-x_R)}{2D},\quad D=x_R-x_L\ (D=0\Rightarrow1),$$
  压力项 $\psi^P_f=1-\big(1-\psi_{\rm face}[P^{\rm old}]\big)\min(\phi_L,\phi_R)$（$P^{\rm old}$ 为本级 Newton 的初猜压力场，$\phi$ 为 E3.1 存的限制因子：只有所有被输运量都光滑时才降到二阶）；速度项 $\psi^u_f=\psi_{\rm face}[u^*]$。`limiter=0` 时 $\psi_{\rm face}=\mathrm{clip}\big(1-\tfrac14[(x_R-x_m)+(x_p-x_L)]/D,\,0,\,1\big)$；一阶空间时 $\psi=1$。
- 速度跳跃 $D^*_f$（顶点散度形式，B4.1(c)）：一维，或 $j$ 行在边界上时 $D^*=u^*_{R}-u^*_{L}$；二维内部
  $$D^*_{i+1/2,j}=\tfrac14\Big[2(u^*_{i+1,j}-u^*_{i,j})+(u^*_{i+1,j+1}-u^*_{i,j+1})+(u^*_{i+1,j-1}-u^*_{i,j-1})\Big]+\tfrac14\frac{\Delta x}{\Delta y}\Big[(v^*_{i,j+1}-v^*_{i,j-1})+(v^*_{i+1,j+1}-v^*_{i+1,j-1})\Big].$$
- **固壁面**：$u_f=u_w$（壁速度，静止壁 0），$p_f$ 照上式（虚格 $P=$ 本格）。
- 增量 $\Delta u_f=u_f-\bar u^e_f$（$\bar u^e_f$ 来自 E3.2）。

$k_f,w_L,p^v_f,\bar u^*_f$ 一级内冻结；$u_f,p_f$ 对 $P$ 线性。

### E4.3 增量通量与格更新（B4.2–B4.6，`face_x/y`、`inc`、`face_flux`、`cell_residual`）

面上的增量权重

$$
w^L_f=(u_f)^+-(\bar u^e_f)^+,\qquad w^R_f=(u_f)^--(\bar u^e_f)^-\qquad(\text{固壁且 }u_w=0\text{ 时 }w^L=w^R=0),
$$

增量通量（同一组系数输运六个体积量，供体为单位体积内容）

$$
I_f[X]=\hat X_Lw^L_f+\hat X_Rw^R_f,\quad X\in\{\alpha_1,\alpha_2,m_1,m_2,\varepsilon_1,\varepsilon_2\};\qquad
M^{L}_f=(\hat m_1+\hat m_2)_Lw^L_f,\ M^{R}_f=(\hat m_1+\hat m_2)_Rw^R_f,\ \hat M_f=M^L_f+M^R_f .
$$

给定本格压力 $P_i$ 与四个邻格压力（虚格按 E7 取本格或周期），格 $i$ 的输运态与守恒量：

$$
\begin{aligned}
X^{\rm tr}_i&=X^*_i-r_x\,\delta_i[I_f[X]]-r_y\,\delta_i[I_f[X]]_y ,\qquad s_i=1-\alpha^{\rm tr}_{1,i}-\alpha^{\rm tr}_{2,i},\\
(\rho u)_i&=(\rho u)^*_i-r_x\,\delta_i[p_f]-r_x\,\delta_i[u_f\hat M_f]-r_y\,\delta_i[\hat u_f\hat M_f]_y ,\\
(\rho v)_i&=(\rho v)^*_i-r_y\,\delta_i[p_f]_y-r_y\,\delta_i[v_f\hat M_f]_y-r_x\,\delta_i[\hat v_f\hat M_f],\\
E_i&=E^*_i-r_x\,\delta_i\Big[I_f[\varepsilon_1]+I_f[\varepsilon_2]+\tfrac12u_f^2\hat M_f+\tfrac12\big(v^{*2}_LM^L_f+v^{*2}_RM^R_f\big)+(P_{\rm ref}+p_f)u_f\Big]-r_y[\cdot]_y ,
\end{aligned}
$$

切向动量通量 $\hat v_f\hat M_f\equiv v^*_LM^L_f+v^*_RM^R_f$（切向速度按供体取）；$y$ 面的切向是 $u^*$。守恒的混合内能 $(\rho e)^C_i=E_i-\big[(\rho u)_i^2+(\rho v)_i^2\big]/(2\rho^{\rm tr}_i)$，$\rho^{\rm tr}_i=m^{\rm tr}_{1,i}+m^{\rm tr}_{2,i}$。

### E4.4 闭合输入与残量（B5.3）

闭合（E5）的输入：$\alpha^{\rm tr}_k,m^{\rm tr}_k,\varepsilon^{\rm tr}_k,s,(\rho e)^C,\tau$、类型、$Z_k$、$p^{\rm guess}=P_{\rm ref}+P'_i$、$Q_k$。输出 $\alpha_k,m_k,\varepsilon_k,p_k,T_k$ 与 $\bar p^{\rm cl}=\sum_k\alpha_kp_k$。格残量

$$
G_i(P')=\big(\bar p^{\rm cl}_i-P_{\rm ref}\big)-P'_i .
$$

闭合无解（不可容许）时 $G_i$ 未定义：Jacobian 构造失败则本级失败；线搜索中则该步长不被接受。

### E4.5 全局 Newton（B8，`Newton_Solver`）

- 残量范数 $\lVert G\rVert=\max_i\lvert G_i\rvert/N_i$，$N_i=\lvert P_{\rm ref}+P'_i\rvert+K_i+P_{E,i}+1$；收敛判据 $\lVert G\rVert<10^{-13}$。
- $P$ 的自然尺度 $s_i=N_i/D_i$，$D_i=1+(K_i+P_{E,i})\big[r_x(k_{i+1/2}+k_{i-1/2})+r_y(k_{j+1/2}+k_{j-1/2})\big]$。
- Jacobian：逐格差分，步长 $h_i=10^{-7}s_i$；行 $i$ 存 5 个元素（本格、W、E、S、N，一维 3 个），邻格扰动用邻格自己的步长；同列（边界虚格 = 本格）也分开存不合并。
- 线性解：一维非周期 Thomas（同列元素相加）；否则 BiCGStab（函数名沿用 `GMRES_Solver`），预条件 $M^{-1}r=D^{-1}r+\mathbf 1\,(\mathbf 1^{\rm T}r)/(\mathbf 1^{\rm T}J\mathbf 1)$（Jacobi 用本格扰动的那个元素 + 常数向量粗空间），相对容差 $10^{-13}$，至多 4000 次；内积按固定块求和。
- 步长判据：完整 Newton 步满足 $\max_i\lvert\delta P_i\rvert/(\lvert P'_i\rvert+s_i)<10^{-13}$ 时直接接受 $P'+\delta P$（重算残量后返回）。
- 否则线搜索：$\lambda=1,\tfrac12,\dots$（至多 40 次），接受条件：所有格闭合有解且 $\lVert G_{\rm new}\rVert<(1-10^{-4}\lambda)\lVert G\rVert$（或已低于容差）。
- 至多 60 次 Newton；失败 → 本级失败（E2.3 第 4 项）。日志 `maxres` 是各级接受时的残量最大值。

### E4.6 激波预算（B5.2 (c)，`shock_budget`、`stencil_state`、`phys_dissipation`）

对每个内部格，在 $\pm2$ 格的一维模板（二维为十字：$i\pm1,i\pm2$ 与 $j\pm1,j\pm2$）上用 $X^*$ 的状态：

- 每个模板格的压力 $p^*=p_{\rm mix}(\alpha_1^*,m_1^*,m_2^*,\varepsilon_1^*+\varepsilon_2^*)$（$\alpha_2$ 取 $1-\alpha_1$），两相在该 $p^*$ 下都可容许才计入，否则跳过；$p_{\rm lo}=\min p^*$（记该格为波前格），$p_{\rm hi}=\max p^*$。
- 压缩判据：$\mathrm{div}=(u^*_{i+1}-u^*_{i-1})+(v^*_{j+1}-v^*_{j-1})<0$。
- 压缩且 $p_{\rm hi}>p_{\rm lo}$ 时，本格的物质（**本格自己的**质量 $m_k$）按 (II.5) 应得的耗散（每单位本格体积）为 $Q_k=m_k\,d_k$，波前比容 $v_k^{\rm pre}=\alpha_k/m_k$ 取波前格（相在波前格中不存在时取本格）；否则 $Q_k=0$。每单位质量的能隙 $d_k$：

$$
d_k=\begin{cases}(v^{\rm pre}_k-b_k)\,f_k,\quad f_k=\dfrac{(p_{\rm hi}+\gamma_k\pi_k)(\tilde v_H-\tilde v_S)}{\gamma_k-1},\ \ \tilde v_H=\dfrac{(\gamma_k-1)x_1+(\gamma_k+1)x_0}{(\gamma_k+1)x_1+(\gamma_k-1)x_0},\ \ \tilde v_S=\Big(\dfrac{x_0}{x_1}\Big)^{1/\gamma_k},\ x_0=p_{\rm lo}+\pi_k,\ x_1=p_{\rm hi}+\pi_k&\text{NASG}\\[6pt]
e_k(p_{\rm hi},v_H)-e_k(p_{\rm hi},v_S),\quad e_k(p,v)=\varepsilon_r(1/v)+v\,(p-p_r)/\Gamma_k；v_H\text{ 解 }e_k(p_{\rm hi},v)-e_k(p_{\rm lo},v_0)+\tfrac12(p_{\rm lo}+p_{\rm hi})(v-v_0)=0，&\text{MG}\\ \quad v_S\text{ 解 }p_r(1/v)+(p_{\rm lo}-p_r(1/v_0))(v_0/v)^{\Gamma_k+1}=p_{\rm hi}\text{（等熵线）；}v_0=v^{\rm pre}_k\text{，各在 }(0,v_0]\text{ 上按 }0.9\text{ 倍缩小找变号再二分 100 次}\end{cases}
$$

用本格自己的质量而不是波前格的 $\alpha_k$：激波穿过接触/混合区时，被压缩的是本格的物质；组分均匀时两种取法给出同一份额。（$d_k\ge0$；MG 求根失败则 $Q_k=0$。）$Q_k$ 是模型激波闭合 (II.5) 下相 $k$ 在 $p_{\rm lo}\to p_{\rm hi}$ 的激波中应得的耗散，弱激波时 $\to\alpha_k\mathcal G_k\Delta p^3/(6K_k^2)$。

### E4.7 终态

用解出的 $P'$ 对每格再算一次 E4.3–E4.4：写回 $\alpha_1=\alpha_1^{\rm cl}$、$m_k$、$(\rho u)$、$(\rho v)$、$\varepsilon_k=\varepsilon_k^{\rm cl}$、$E$（$E$ 是 E4.3 的守恒值，与 $\sum\varepsilon_k+\tfrac12\rho\lvert\mathbf u\rvert^2$ 一致到舍入）；填虚格。记录 $T_{\max}$（$\alpha_k>10^{-9}$ 的格）。

---

## E5. 统一闭合（B5、B6，`closure` in `spb6_physics.h`）

输入见 E4.4（显式程序与投影：$s=0$）。先算共用量，再按类型求解。所有一维根都带区间（`root_above`：在 $\ln(p-p_{\rm lo})$ 上从初猜向两侧几何扩张找变号，再 Brent；`brent`：标准 Brent）。

### E5.0 能量亏损的分配（B5.2，V16，`defect_split`、`p_coef`）

$$
D_0=(\rho e)^C-\varepsilon^{\rm tr}_1-\varepsilon^{\rm tr}_2,\qquad \delta(p)=D_0+p\,s\quad(\text{真实亏损}),\qquad \delta_g=\delta(p^{\rm guess}).
$$

可容许界之上的能量 $E_k=\max\big(0,\ \varepsilon^{\rm tr}_k-m_kq_k-(\alpha^{\rm tr}_k-m_kb_k)\pi_k\big)$（MG：$\max(0,\varepsilon^{\rm tr}_k-m_k\varepsilon_r(m_k/\alpha^{\rm tr}_k))$），$\lambda^E_k=E_k/\sum_jE_j$（$\sum E=0$ 时按质量 $m_k/\sum m$）；$Q=Q_1+Q_2$，$\lambda^{\rm II.5}_k=Q_k/Q$。分配函数

$$
\Lambda_k(\delta)=\begin{cases}\delta\,\lambda^E_k,&\delta<0\ \text{或}\ Q=0,\\ \delta\,\lambda^{\rm II.5}_k,&0\le\delta<Q,\\ \delta\,\lambda^E_k+\big(Q_k-Q\lambda^E_k\big),&\delta\ge Q,\end{cases}
$$

连续、分段仿射；按 $\delta_g$ 所在的段取 $\Lambda_k=\lambda_k\delta+\xi_k$（$\xi_k$ 只在第三段非零）。闭合系数

$$
A_k=\varepsilon^{\rm tr}_k-m_kq_k+\lambda_kD_0+\xi_k,\qquad B_k=\alpha^{\rm tr}_k+\lambda_ks,\qquad N_k=A_k-(B_k-m_kb_k)\pi_k\ (>0\ \text{才可容许}),
$$

于是 $A_k+B_kp=\varepsilon^{\rm tr}_k-m_kq_k+p\alpha^{\rm tr}_k+\Lambda_k(\delta(p))$：相 $k$ 得到自己的压缩功 $-p(\alpha_k-\alpha^{\rm tr}_k)$ 加分到的亏损 $\Lambda_k$。`nodelta`（只用于激活判据）时 $\lambda_k=\xi_k=0$。

### E5.1 $p$ 闭合（$\mu=\infty$，$\theta=\nu=0$；B6.1，`closure_p`）

NASG：

$$
\alpha_k(p)=\frac{(\gamma_k-1)(A_k+B_kp)+(p+\gamma_k\pi_k)m_kb_k}{\gamma_k(p+\pi_k)},\qquad \alpha_k'(p)=-\frac{(\gamma_k-1)N_k}{\gamma_k(p+\pi_k)^2}<0,
$$

解 $\Phi(p)=\alpha_1(p)+\alpha_2(p)-1=0$，$p>p_{\rm lo}=\max_k(-\pi_k)$：先从 $p^{\rm guess}$ 找区间（$\Phi>0$ 的下端、$\Phi<0$ 的上端，向上 ×4 或向下 ×1/4 扩张），再带区间的 Newton（步出区间时取几何/算术中点），终止于相对 $2\times10^{-16}$。输出 $\alpha_k=\alpha_k(p)$，$m_k=m^{\rm tr}_k$，$\varepsilon_k$ 由 E1 的 EOS，$p_k=p$，$T_k$ 由 E1；要求 $\alpha_k>m_kb_k$。

MG 相（`closure_p_mg`，`alpha_mg`）：$\alpha_k(p)$ 是 $f_k(\alpha)=m_k\varepsilon_r(m_k/\alpha)+\alpha(p-p_r(m_k/\alpha))/\Gamma_k+p\alpha-A_k-B_kp=0$ 的根（$f_k'=\rho_kc_k^2/\Gamma_k>0$；在 $\alpha$ 上从初猜几何扩张找区间，再带区间的 Newton，$A_k$ 中无 $m_kq_k$ 项）；外层 $\Phi(p)$ 用带区间的 Newton，$\alpha_k(p)$ 无解当作 $\Phi=+\infty$。

### E5.2 有限 $\mu$（$\theta=\nu=0$；B6.2，`closure_mu`）

$\omega_k=\alpha^{\rm tr}_k/\sum_j\alpha^{\rm tr}_j$，$a_k=\omega_k$（即 $\alpha^{\rm tr}_k+\omega_ks$），$\alpha_1=a_1+R$，$\alpha_2=a_2-R$；界面压力 $p_I=w_1p_1+w_2p_2$，$w_1=Z_2/(Z_1+Z_2)$（$Z_k$ 冻结）。给定 $R$，解 2×2 线性方程组（$k=1,2$）

$$
\Big(\frac{\alpha_k-m_kb_k}{\gamma_k-1}+\omega_ks\Big)p_k+\chi_kR\,p_I-\lambda_ks\,(\omega_1p_1+\omega_2p_2)
=\varepsilon^{\rm tr}_k-m_kq_k-\frac{\gamma_k\pi_k(\alpha_k-m_kb_k)}{\gamma_k-1}+\lambda_kD_0+\xi_k ,
$$

再解 $\mathcal F(R)=R-\tau\mu\,(p_1-p_2)=0$：根在 $0$ 与 $R_1=-\mathcal F(0)$ 之间；若 $R_1$ 不可容许（$\alpha_k\le m_kb_k$ 或 $p_k+\pi_k\le0$）向 0 二分收缩区间，再 Brent。$\tau\mu=0$ 时 $R=0$（冻结）。输出 $\alpha_k$、$p_k$、$\varepsilon_k=\varepsilon_k(\alpha_k,m_k,p_k)$、$\bar p^{\rm cl}=\alpha_1p_1+\alpha_2p_2$。

### E5.3 $pT$ 闭合（$\mu=\theta=\infty$，$\nu=0$；B6.3，`closure_pT`）

$m_k=m^{\rm tr}_k$。$T(p)=\dfrac{(\rho e)^C-\sum_km_kq_k}{\sum_km_kc_{vk}(p+\gamma_k\pi_k)/(p+\pi_k)}$，解 $\Phi(p)=\sum_km_kv_k(p,T(p))-1=0$（`root_above`）。输出 $\alpha_k=m_kv_k$，$\varepsilon_k=m_ke_k(p,T)$，$p_k=p$，$T_k=T$。与 $\alpha^{\rm tr},s,\Lambda$ 无关。

### E5.4 有限 $\theta$（$\mu=\infty$，$\nu=0$；B6.4，`closure_theta`）

在 E5.1 的 $\alpha^0_k(p)$、$T^0_k(p)$ 上加传热 $\Delta Q$：

$$
\Delta Q(p)=\frac{\tau\theta\,\big(T^0_2(p)-T^0_1(p)\big)}{1+\tau\theta\,(1/C_{p1}+1/C_{p2})},\quad C_{pk}=\gamma_kc_{vk}m_k;\qquad
\alpha_k=\alpha^0_k(p)+\chi_k\frac{(\gamma_k-1)\Delta Q}{\gamma_k(p+\pi_k)} ,
$$

解 $\Phi(p)=\alpha_1+\alpha_2-1=0$（`root_above`）；输出同 E5.1（$\varepsilon_k$ 由 $\alpha_k$ 与 $p$）。$\theta=0$ 时逐位回到 E5.1。

### E5.5 有限 $\nu$（$\mu=\theta=\infty$；B6.5，`closure_nu`）

$\rho=\sum m^{\rm tr}_k$，$e^C=(\rho e)^C/\rho$。给定 $p$：

$$
T(p,Y)=\frac{e^C-(1-Y)q_1-Yq_2}{(1-Y)c_{v1}\phi_1+Yc_{v2}\phi_2},\ \phi_k=\frac{p+\gamma_k\pi_k}{p+\pi_k};\qquad
\mathcal R(Y)=\rho Y-m^{\rm tr}_2+\tau\nu\,\big(g_2-g_1\big)\big(p,T(p,Y)\big)=0
$$

在 $Y\in[Y_{\min},1-Y_{\min}]$ 与 $T>0$ 的交集上（$q_2>q_1$ 时上界 $(e^C-q_1)/(q_2-q_1)$，$q_2<q_1$ 时下界）：$\mathcal R$ 在下端 $\ge0$ 取下端，上端 $\le0$ 取上端，否则 Brent。外层 $\Phi(p)=\rho\sum_kY_kv_k(p,T)-1=0$。输出 $m_1=\rho(1-Y)$，$m_2=\rho Y$，$\alpha_k=m_kv_k$，$\varepsilon_k=m_ke_k$。

### E5.6 $pTg$ 闭合（$\mu=\theta=\nu=\infty$；B6.6，`closure_pTg`）

给定 $p$：$T_s=T_{\rm sat}(p)$；若有限，$Y=\dfrac{e^C-e_1(p,T_s)}{e_2(p,T_s)-e_1(p,T_s)}$；$Y\in(Y_{\min},1-Y_{\min})$ 时 $T=T_s$，否则钳到界值并由能量定 $T=T(p,Y)$（E5.5 的公式）。$p\ge p_c$（$T_{\rm sat}$ 无定义）：$1/\rho<v_c$ 取 $Y=Y_{\min}$，否则 $1-Y_{\min}$，$T$ 由能量定。外层 $\Phi(p)=\rho\sum_kY_kv_k(p,T)-1=0$。输出同 E5.5。

### E5.7 类型选择与松弛档位（B6.7，`closure_type_of`，`RelaxationActivation`）

速率：$\mu$（`mu=`，默认 $\infty$）、$\theta$（`theta=`）、$\nu$（`nu=`）；`relax=` 档位：0 只有力学；1 $p$–$pT$；2 $p$–$pT$–$pTg$；3 $p$–$pTg$。档位 $\ge1$ 时 $\theta$ 默认 $\infty$，档位 2、3 时 $\nu$ 默认 $\infty$；命令行给有限值即为有限速率。**有限速率与无限速率走同一个闭合**：无限速率是 E5.1/E5.3/E5.6（平衡流形上的投影），有限速率是 E5.2/E5.4/E5.5（后向 Euler 的松弛 ODE，$\tau=a_{ss}\Delta t$ 为级的步长；速率 $\to\infty$ 时连续地趋于平衡闭合，$=0$ 时逐位退化为冻结，B7）。

每格类型（输入该格的 $\alpha_1,m_k,\varepsilon_k$ 与比内能）：
1. 力学类型 $\mathcal M$ = $p$（$\mu=\infty$）或有限 $\mu$；档位 0 → $\mathcal M$。
2. 非界面格（$\alpha_1<\epsilon_I$ 或 $>1-\epsilon_I$）→ $\mathcal M$。
3. 档位 3：以 $s=0$、不含 $\delta$ 的 $p$ 闭合算力学平衡态，任一相 $T_k>T_{\rm sat}(p)$ → 化学（$pTg$ 或有限 $\nu$），否则 $\mathcal M$。
4. 档位 1、2：热类型 $\mathcal T$ = $pT$（$\theta=\infty$）或有限 $\theta$；档位 2 时先算 $pT$ 平衡态，其 $T>T_{\rm sat}(p)$ → 化学，否则 $\mathcal T$。

类型按 $U^n$ 在每步开始判断一次，$pT$、$pTg$ 类型的格把 $U^n$ 投影到流形上（$s=0$ 的同一闭合；无解则该格退回 $\mathcal M$ 且不投影），本步各级冻结；输出前再做一次。MG 相的算例只允许 $p$ 类型。

---

## E6. 显式参考程序的路径（A3–A8，`SPB6_explicit.cpp`）

1. 步长：声学 CFL 0.4。
2. Heun 每一级：原始变量 $(\alpha_1,\rho_1,\rho_2,u,v,p_1,p_2)$ 的 minmod 重构 → HLLC（A4；相能量用各相 Hugoniot，MG 相用 A2.8）→ SPB 的 Godunov 型更新（A5）→ 瞬时松弛投影：先填虚格，按 E4.6 在 $\pm2$ 格模板上算每格的 $Q_k$（压力用各格的 $p_{\rm mix}$、散度用格速度的中心差分），再对每格做 E5 的闭合（$s=0$，$\tau=$ 本步 $\Delta t$，E5.0 的分配同 IMEX）。
3. 松弛起点：输运态不可容许的相（$p_k+\pi_k\le0$ 或 MG 的 $c^2\le0$）先按 SPB §3.4 取混合压力下的能量（差值进 $\delta$），再做闭合（A6，V13）。
4. 热/化学类型按 E5.7 判断（在闭合之后，`closure_type_of`），有限速率松弛每步只在 Heun 平均之后做一次（$\tau=\Delta t$，A3，V6）。
5. 拒绝与失败处理同 IMEX。

---

## E7. 边界、初值、输出、并行（B11，`bc_apply`）

- 两层虚格。每一边单独（左、右、下、上）：**透射**——虚格复制最近的内部格；**固壁**——镜像（第 $q$ 层虚格取第 $q$ 层内部格），法向动量取 $2\rho u_w-(\rho u)$（$u_w=0$ 时取负），总能按动能的变化修正；**周期**——绕回。一维（$n_y=1$ 且 $y$ 非周期）时 $y$ 向只复制。隐式部分：固壁面 $u_f=u_w$（E4.2）；虚格压力 $P=$ 本格（周期时绕回），Jacobian 里该列单独存。
- 初值：圆形区域的格按 $4\times4$ 子点上守恒量的平均初始化（`sub=4` 的算例）。
- 输出列（`write_profile`）：`x y p u v rho alpha2 Y2 rho1 rho2 c_f c_W T1 T2 p1 p2 g1-g2`（`prec=` 有效位数，默认 16）；`snap=` 写 `<out>_s<k>.dat`；`probe=` 写 `<out>.probe`。
- OpenMP：逐格循环、Jacobian、BiCGStab；内积按固定块求和，结果与线程数无关；一维输出与串行逐字节相同。

---

## E8. 这条路径上的已知限制（C 页“仍未解决”）

SPB 4.3.2 的起步误差（激波超前约 15 格，恒定、一阶收敛；V16 后 $L_1(u)$ 480）；$M\le10^{-5}$ 的双精度底（不修）；活塞为固定网格近似；NASG 十二烷平台窗口误差（0.070）；Pelanti 8.2.2 $pTg$ 的少量拒绝。

---

## E9. 各算例（`test_setup` / `test_state`；数据与图见 `data/`、`docs/figures/`）

SG/NASG 参数记为 $(\gamma,\pi,q,q',c_v,b)$，省略的 $q=q'=b=0$、$c_v=1000$；MG 记为 (类型; $\Gamma,\rho_r,c_1,c_2,e_1,e_2,e_0$)。区域默认 $[0,1]$，边界默认透射；初值“$L\mid R$”为界面两侧的 $(\alpha_1,\rho_1,\rho_2,u,p)$，微量相 $\alpha=10^{-6}$（SPB）或 $10^{-8}$（十二烷、MG）。

| `case` | 论文 / 图 | 相 1 | 相 2 | 初值 | $t_{\rm end}$ | 网格 | CFL | Rusanov | 松弛 |
|---|---|---|---|---|---|---|---|---|---|
| 4 | SPB 4.1 / 图 5 | 水 (4.4, 6e8) | 空气 (1.4, 0) | $x<0.5$：$(1-10^{-8},1000,10,100,10^5)$，右侧 $\alpha_1=10^{-8}$ | 2.79e-3 | 200 | 0.2 | 局部 | $p$ |
| 60 | SPB 4.2 / 图 6–7 | 硝基甲烷 CC（1.19, 1134, 0.819181e9, 1.50835e9, 4.52969, 1.42144） | 同一流体 | $x<0.5$：$p=2\times10^{10}$，$\rho=1134$；右 $p=2\times10^5$，$\rho=1200$ | 67e-6 | 500 | 0.2 | 局部 | $p$ |
| 5 | SPB 4.3.1 / 图 8–9 | 水 (4.4, 6e8) | 空气 (1.4, 0) | $x<0.75$：水 $p=10^9$；右：空气 $\rho=1$，$10^5$ | 240e-6 | 1000、100 | 0.2 | 局部 | $p$ |
| 6 | SPB 4.3.2 / 图 10 | 同上 | 同上 | $x<0.6$（论文图 10 反推，C 页 §15.8）：水 $p=10^{12}$；右：空气 $\rho=10$，$10^5$ | 8.3e-6 | 1000 | 0.2 | 局部 | $p$ |
| 0 | SPB 4.5 / 图 12 | 水 | 空气 | $\alpha_{\rm air}=10^{-2}$，$\rho=(1000,1)$，$u=\mp100$，$p=10^5$ | 1.85e-3 | 1000 | 0.2 | 局部 | $p$ |
| 12 | SPB 5.1 / 图 19 | 环氧 (2.43, 5.3e9) | 尖晶石 (1.62, 141e9) | $\alpha_1=0.5954$，$\rho=(1185,3622)$；$x<0.6$：$p=10^{10}$，右 $10^5$ | 80e-6 | 500 | **0.05** | 局部 | $p$ |
| 13 | SPB 5.2 / 图 20 | 同上 | 同上 | 同上，$p_L=2\times10^{11}$ | 34e-6 | 500 | 0.2 | **全场** | $p$ |
| 1, 2 | P-S 2014 5.1 / 图 1、2 | 水 (2.35, 1e9, −1167e3, 0, 1816) | 汽 (1.43, 0, 2030e3, −23.4e3, 1040) | $\alpha_v=10^{-2}$，$T=354.728$ K，$p=10^5$，$u=\mp2$ / $\mp500$ | 3.2e-3 / 0.58e-3 | 5000 | 0.2 | 局部 | 0–3 |
| 3 | P-S 2014 5.2 / 图 3 | 液十二烷 (2.35, 4e8, −775.269e3, 0, 1077.7) | 汽 (1.025, 0, −237.547e3, −24.4e3, 1956.45) | $x<0.75$：液 $\rho=500$，$p=10^8$；右：汽 $\rho=2$，$10^5$ | 473e-6 | 2000 | **0.05** | 局部 | 0、2、3 |
| 30 | Pelanti 8.2.1 / 图 4–5 | 水 NASG (1.187, 7028e5, −1177788, 0, 3610, 6.61e-4) | 汽 (1.467, 0, 2077616, 14317, 955) | $\alpha_v=10^{-2}$，$T=353$ K，$p=10^5$，$u=\mp2$ | 3e-3 | 5000 | 0.2 | 局部 | $p$、$pT(f)$ $\theta=2000$、$pT$、$pTg(f)$ $\nu=10^{-4}$、$pTg$ |
| 31 | Pelanti 8.2.2 / 图 6 | 十二烷 NASG (1.206, 1681e5, −996054, 0, 2532, 7.51e-4) | 汽 (1.021, 0, −384592, −4301, 2274) | 同算例 3 | 473e-6 | 2000 | 0.2 | 局部 | $p$、$pTg$ |
| 61 | Pelanti 8.1.1 / 图 1 | 爆轰气 JWL（0.25, 1840, 854.5e9, 20.5e9, 4.6, 1.35, 0） | 水 NASG（同 30） | $[-10,10]$，$x<0$：$\rho=2000$，$p=4.6406\times10^{10}$；右：$\rho=1044$，$10^5$ | 0.9e-3 | 500 | 0.2 | 局部 | $p$ |
| 70 | SPB 4.6 / 图 13–14 | 空气 (1.4, 0) | 氦 (1.67, 0) | $[0,0.3]\times[0,0.08]$，上下固壁；氦泡心 $(0.2,0.04)$ 半径 0.02，$\rho=0.167$；空气 $1.29$ kg/m³、101325 Pa；$x>0.221$ 为 Mach 1.5 波后态（向左） | 600e-6 | $600\times160$ | 0.2 | 局部 | $p$ |
| 71 | SPB 4.7 / 图 15–17 | 水 | 气 (1.8, 0) | $[0,3]\times[0,1]$，左、上、下固壁；水在 $x<1.2+\sqrt{0.36-(y-0.5)^2}$，$\rho=(1000,100)$，$u=-200$，$p=101325$ | 8.6e-3 | $300\times100$ | 0.2 | **全场** | $p$ |
| 72 | P-S 2014 5.4 / 图 5 | 十二烷（同 3） | 汽（同 3） | 单位正方形，四周固壁，左壁活塞 $u_w=100$；汽泡心 $(0.5,0.5)$ 半径 0.2；$p=10^5$，$T=600$ K | 1.2e-3 | $200\times200$ | 0.2 | **全场** | 2 |
| 73 | Pelanti 8.1.2 / 图 2–3 | JWL（同 61） | 水 NASG | $[-1.35,1.35]\times[-1.5,0]$，上边固壁；气泡心 $(0,-0.32)$ 半径 0.05，内 $\rho=2000$，$4.6406\times10^{10}$ Pa；水 $1044$，$10^5$ | 0.4e-3 | $361\times200$（论文 $721\times400$） | 0.2 | 局部 | $p$ |
| 50–55 | Gresho 涡，$M=10^{-1}\ldots10^{-6}$ | 理想气体 (1.4, 0) | 同 | $[-0.5,0.5]^2$ 周期；$\rho=1$，$p_0=1/(\gamma M^2)$，Gresho 剖面 | $0.4\pi$ | $80^2$ | `dt=0.0025` | 局部 | $p$ |
| 56 | Battisti–Boscheri §4.1 两相等熵涡 | (1.35, 0) | (1.4, 0) | $[-5,5]^2$ 周期；Hu–Shu 涡 $\epsilon=5$，$T=p_0+\delta T$，$\rho=(T/p_0)^{1/0.4}$，$p=\rho T$，$\alpha_2=0.5+0.25e^{(1-r^2)/2}$，平移 $(2,2)$；`p0=` 背景压力 | 5 | $40^2$（试验 32–512） | 0.2 | 局部 | $p$ |
| 20, 21 | 光滑解（时间收敛阶） | (1.4, 0) / 水 | (1.4, 0) | $u=50e^{-(x-0.5)^2/0.01}$，$\alpha_1=0.5+0.2\sin2\pi x$，$p=10^5(1+0.2e^{\cdots})$ | 5e-4 / 1e-4 | — | 0.2 | 局部 | $p$ |
| 40, 41 | 0D 松弛（周期） | 水（同 1） | 汽（同 1） | 40：$\alpha_1=0.8$，$(T,p)=(350,1.2\times10^5)\mid(380,0.8\times10^5)$；41：$T=380$，$p=0.8\times10^5$ | 1e-4 | — | 0.2 | 局部 | 有限速率 |

选择 CFL 0.05 的理由：声波本身是要解的结构，物质 CFL 0.2 时其声学 Courant 数接近 1，LSDIRK2 的色散使平台偏低（C 页 V14 §14.3）。选择全场 Rusanov 速度的理由：空化 RM、汽泡压缩用局部速度时不稳健（C 页 V14 §14.5、§14.9）；SPB 5.2 在 V14 时全场更好（$L_1(u)$ 118 对 124）。

---

## E10. 日志行与诊断

`case N mode 1 order 2 rk 2 N nx CFL C | steps S (rk1 R1, mood M) rejects R | stages St newton it I (I/St) maxres r Tmax T1 T2 | closure fails F | dM1/M dM2/M dE/E | wall s`：`rk1` 为用后向 Euler 的步数（斜坡 + MOOD），`mood` 为 MOOD 回退的步数，`rejects` 为减半次数，`newton it` 为 Newton 迭代总数，`maxres` 为各级接受时残量的最大值，`Tmax` 为两相温度在整个计算中的最大值（$\alpha_k>10^{-9}$ 的格），`closure fails` 为显式程序闭合失败数，最后三项是质量与总能的守恒误差（有开边界时含流出）。算例 56 另打印 $L_2$ 误差行。
