# E. 唯一代码路径（最终版，V13）：从输入到输出

本章把最终固定下来的**唯一**代码路径从头到尾写一遍：每一步用哪个公式、在代码哪个函数、为什么这样选（推导节号）。程序里没有对照开关：所有算例走同一条路径，只有算例自带的物理参数、区域、边界、网格、终止时间不同。

- IMEX：`src/SPB6_TV_IMEX.cpp`（推导 B 章）。显式参考：`src/SPB6_explicit.cpp`（推导 A 章）。两者共用 `src/spb6_physics.h`（状态方程、闭合、算例、边界条件）。
- 两个程序都直接由 `reference/code/2DEuler_TV.cpp` 改来（函数名、`Grid`/`Array`、`Qe/Qs/Hs` 的 IMEX-RK 循环、VTK 输出保留）。
- 每处选择的来由、试过而放弃的做法见 C. 开发记录（`docs/devlog.md`）。

---

## E0. 固定的设置（所有算例相同）

| 项 | IMEX | 显式参考 |
|---|---|---|
| 时间积分 | LSDIRK2(2,2,2)（$\gamma=1-1/\sqrt2$，`rk=2`）；起步斜坡内的步用后向 Euler | Heun（二阶 SSP-RK） |
| 时间步 | 物质 CFL：$\Delta t=0.2/(\max\lvert u\rvert/\Delta x+\max\lvert v\rvert/\Delta y)$ | 声学 CFL：$\Delta t=0.4/(\max(\lvert u\rvert+c_f)/\Delta x+\dots)$ |
| 起步 | 首步取声学 CFL 0.4 的步长，之后每步至多放大 1.2 倍，被限制的步用一阶 | — |
| 空间 | 二阶 MUSCL（minmod），B3 | 原始变量 minmod（SPB 附录 C） |
| 通量 | Rusanov，耗散速度为全局物质速度 $a_\Omega$（BB (27)） | HLLC（SPB 3.1.2） |
| 压力/声学 | 隐式（B4–B8） | 显式（HLLC） |
| 松弛 | 统一闭合，在隐式压力方程里（B5、B6） | 每级之后的瞬时投影 + 每步一次有限速率（A3） |
| Newton 容差 | $10^{-13}$（相对 $\lvert P\rvert+K+P_E$），至多 60 次 | — |
| 线性解法 | 一维非周期 Thomas；其余 BiCGStab + 两级预条件，相对容差 $10^{-13}$，至多 4000 次 | — |
| 界面判据 $\epsilon_I$ | $10^{-4}$（热、化学松弛只在 $\epsilon_I\le\alpha_1\le1-\epsilon_I$ 的格） | 同左 |

CFL 是两个程序各自的唯一值（IMEX 0.2、显式 0.4），所有算例相同；命令行 `cfl=` 只用于收敛阶等专门试验。

---

## E1. 变量与状态方程（A1、A2、A2.1）

守恒变量（每格 8 个）：$\alpha_1,\ m_1=\alpha_1\rho_1,\ m_2,\ \rho u,\ \rho v,\ \varepsilon_1=\alpha_1\rho_1e_1,\ \varepsilon_2,\ E$。$\alpha_2=1-\alpha_1$；$E$ 守恒（SPB 的 “+1”），相内能非守恒。

每相的状态方程（`Phase::kind`）：

- **NASG / SG**（`kind = 0`）：$\varepsilon_k=(p+\gamma_k\pi_k)(\alpha_k-m_kb_k)/(\gamma_k-1)+m_kq_k$（A2.2），$c^2=\gamma(p+\pi)/(\rho(1-\rho b))$，$T$、$s$、$g$ 显式（A2.1）。可容许：$p+\pi>0$。
- **Mie–Grüneisen**（`kind = 1` Cochran–Chan，`kind = 2` JWL）：$p=\Gamma\rho(e-\varepsilon_r(\rho))+p_r(\rho)$，$\Gamma$ 常数，$\varepsilon_r'=p_r/\rho^2$（A2.4）；$\varepsilon_k=m_k\varepsilon_r(\rho_k)+\alpha_k(p-p_r)/\Gamma_k$（A2.5），$c^2=(\Gamma+1)(p-p_r)/\rho+p_r'$（A2.6）。可容许：$c^2>0$。只做力学松弛。
- 混合压力（A2.3、A2.7）对 $p$ 线性；冻结声速 $c_f^2=\sum Y_kc_k^2$；Wood 体积模量 $K=1/\sum_k\alpha_k/(\rho_kc_k^2)$。

---

## E2. 一个时间步（IMEX，`main` 的时间循环）

1. **步长。** `caldt`（物质速度）× CFL 0.2；起步斜坡（首步为声学步长 / 1.2，之后每步至多 ×1.2）；截到终止时间与输出时刻（`snap=`）。
2. **激活与投影（B6.7，`RelaxationActivation`）。** 按 $U^n$ 判断每格的闭合类型（$p$ / 有限 $\mu$ / $pT$ / 有限 $\theta$ / 有限 $\nu$ / $pTg$；判据：界面格 $\epsilon_I\le\alpha_1\le1-\epsilon_I$，化学松弛另要 $T>T_{\rm sat}(p)$），瞬时类型把 $U^n$ 投影到其平衡流形（$s=0$ 的同一个闭合）。本步各级冻结这些类型。
3. **各级**（LSDIRK2 两级；Butcher 表 `Ae`、`Ai` 为 2DEuler 的形式）：
   1. $Q^e=U^n+\Delta t\sum_{r<s}a^e_{sr}H_r$，$Q^s=U^n+\Delta t\sum_{r<s}a^i_{sr}H_r$，填虚格；$\tau=a^i_{ss}\Delta t$。
   2. 重构（E3.1）在 $Q^e$ 上。
   3. 显式部分（E3.2，`Update`）：$X^*=Q^s-\tau\,\nabla\cdot F^{\rm Rus}(Q^e)$。
   4. 隐式部分（E4，`PressureSubsystem`）：由 $X^*$ 解全局压力方程，得本级终态 $U_s$。
   5. $H_s=(U_s-Q^s)/\tau$。
4. **接受或拒绝。** 任一级出现可容许性失败（$Q^e$ 的 $\rho\le0$、闭合输入不可容许、Newton 不收敛）时拒绝本步，$\Delta t\leftarrow\Delta t/2$ 重算，至多 30 次；$\Delta t<10^{-12}t_{\rm end}$ 时报告失败（B8 的步长控制）。接受时 $U^{n+1}=U_{\rm last}$（刚性精确），压力场 $P$ 留作下一步的 Newton 初猜。
5. 输出：`snap=` 时刻的剖面、`probe=` 的压力历史；结束时写剖面与统计行（步数、一阶步、拒绝、每级 Newton 次数、终残量最大值、守恒误差、用时）。

---

## E3. 显式部分（B3）

### E3.1 重构（`TVDReconstruction`）

- **体积输运量** $(\alpha_1,m_1,m_2,\varepsilon_1,\varepsilon_2)$ **共用一个限制因子**：$\phi=\min_q(\text{minmod 斜率}/\text{中心斜率})$，极值处 $\phi=0$。“常数到舍入”的量不参与（阈值 $10^{-12}\times$ 舍入尺度：$\alpha$ 用 3，$m_k$ 用自身，能量用混合内能；V8）。各量面值 $=x_c\pm\tfrac12\phi\cdot\tfrac12(x_b-x_a)$。共用的理由是 WB：(E1) 要求 $\alpha,m,\varepsilon$ 用同一组系数（B9.1）。
- **速度**单独 minmod（舍入尺度 $\lvert u\rvert+c_f$），面上组装 $\rho u=\rho^fu^f$，$E=\sum\varepsilon^f+\tfrac12\rho^f\lvert\mathbf u^f\rvert^2$（V11）。
- 面上的限制因子 $\phi_f$ 留给隐式部分的二阶化（$\psi$）。

### E3.2 Rusanov 通量（`RusanovC_x/y`、`Update`）

- 中心部分与耗散都用重构面值：$F=\tfrac12(F(Q_L)+F(Q_R))-\tfrac12a_\Omega(Q_R-Q_L)$，**不含压力**（TV 分裂，B1）：对流通量 $(\alpha_1u,\ m_ku,\ \rho u^2,\ \rho uv,\ \varepsilon_ku,\ Eu)$。
- 耗散速度 $a_\Omega$：全局的物质速度最大值，且不小于所有面态速度（V6，Rusanov 单调性条件）；x、y 方向各一个。
- $\alpha_2$ 单独输运（`a2_star`），面体积通量和 $\bar u^e_f=\tfrac12(u_L+u_R)$ 记下（`ube`），给隐式增量用。

---

## E4. 隐式部分（B4–B8，`PressureSubsystem`）

### E4.1 级数据（`StageData`）

由 $X^*$：每格单位体积内容 $(\alpha_k,m_k,\varepsilon_k)/\sum\alpha$、供体速度 $u^*,v^*$、冻结声速 $c_f$、阻抗 $Z=\rho c$、残量尺度 $K+P_E$（$P_E=\sum(\lvert\varepsilon_k\rvert+\text{eps\_round}_k)/\sum\partial\varepsilon_k/\partial p$）、本步的闭合类型。

### E4.2 面量（B4.1）

每个面（x 面 $i+\tfrac12$，y 面对称）：

$$
u_f=\bar u^*_f-k_f(P_R-P_L),\qquad k_f=\frac{\tau}{\bar\rho_fh}+\frac{\varphi\psi^P}{Z_L+Z_R},\qquad
p_f=w_LP_L+(1-w_L)P_R-\varphi\psi^u\frac{Z_LZ_R}{Z_L+Z_R}D^*_f ,
$$

- $\bar\rho_f=\tfrac12(\rho_L+\rho_R)$，$w_L=\rho_R/(\rho_L+\rho_R)$（密度加权面压力，半格动量平衡）。
- 声学补足的权 $\varphi=(1-\tau\min(c_L,c_R)/h)^+$：时间步分辨声波时才补到 Godunov 的耗散，低马赫 / 大步长时为 0（B4.1(b)）。
- 二阶化 $\psi=1-(1-\psi_{\rm face})\phi_f$：$\psi_{\rm face}$ 由上一级压力（或 $u^*$）的 minmod 给出，只在共用限制因子允许二阶处（$\phi_f>0$）降低耗散。
- 速度跳跃 $D^*_f$ 用面两端两个顶点散度的平均乘 $h$（B4.1(c)，V11）：一维与 y 向均匀时就是 $u_{i+1}-u_i$；二维中无散场在它的核里，低马赫一致。
- 固壁面：$u_f=u_w$（静止壁 0；活塞时体积通量从虚格进入）。

### E4.3 增量与格残量（B4.2'、B4.3–B4.6）

- 面增量：总体积通量对 $u_f$ 迎风，显式部分已扫过 $\bar u^e_f\tau$，增量补到 $u_f\tau$：$I_f[X]=\hat X_Lw_L+\hat X_Rw_R$，$w_L=(u_f)^+-(\bar u^e_f)^+$，$w_R=(u_f)^--(\bar u^e_f)^-$。同一组系数输运 $\alpha_k,m_k,\varepsilon_k$（温度一致，B9.2）。
- 格内：$X\leftarrow X^*-r\,\delta[I_f]$（6 个量），体积应变 $s=1-\sum\alpha$；动量 $-r\,\delta[p_f]-r\,\delta[u_f\hat M]$（切向用供体速度）；总能 $-r\,\delta[\sum I_f[\varepsilon]+\tfrac12u_f^2\hat M+\tfrac12\hat v^2\hat M+(P_{\rm ref}+p_f)u_f]$。
- 闭合（E5）给出 $\bar p^{\rm cl}$，残量 $G=(\bar p^{\rm cl}-P_{\rm ref})-P'$。

### E4.4 全局 Newton（B8）

- 未知量是压力扰动 $P'=P-P_{\rm ref}$，$P_{\rm ref}$ 为初始混合压力的最小值（V11.1）。
- Jacobian：逐格差分，每行“本格、W、E、S、N”5 个元素固定次序、同列不合并（V8）；步长 $h=10^{-7}s$，$s=N/D$ 为 $P$ 的自然尺度（$N=\lvert P\rvert+K+P_E+1$，$D=1+(K+P_E)\sum_f(\tau/h_f)k_f$）。
- 线性解：一维非周期 Thomas；其余 BiCGStab，预条件 $M^{-1}r=D^{-1}r+\mathbf 1(\mathbf 1^{\rm T}r)/(\mathbf 1^{\rm T}J\mathbf 1)$（Jacobi + 常数向量粗空间，V11.1）；内积按固定块求和（与线程数无关）。
- 线搜索：步长减半至多 40 次，接受条件为所有格闭合可容许且 Armijo 下降。
- 收敛：残量 $\max\lvert G\rvert/N<10^{-13}$；或完整 Newton 步 $\max\lvert\delta P\rvert/(\lvert P'\rvert+s)<10^{-13}$（步长判据，舍入底）。
- 终态：用解出的 $P$ 再算一次每格（闭合的 $\alpha_k,m_k,\varepsilon_k$、动量、总能），写回 $U_s$。

---

## E5. 统一闭合（B5、B6，`closure` in `spb6_physics.h`）

所有闭合共用离散 Gibbs 关系 (B6.1) $\varepsilon_k=\varepsilon^{\rm tr}_k-p(\alpha_k-\alpha^{\rm tr}_k)+\lambda_k\delta$，能量亏损 $\delta=D_0+ps$ 按质量分数 $\lambda_k=Y_k$ 分给两相（V10，B5.2）。按格的类型：

| 类型 | 速率 | 方程 | 节 |
|---|---|---|---|
| $p$ | $\mu=\infty$ | $\sum_k\alpha_k(p)=1$；NASG 的 $\alpha_k(p)$ 显式，MG 的为单调一维根 | B6.1 |
| 有限 $\mu$ | $\theta=\nu=0$ | $R-\tau\mu(p_1-p_2)=0$，$p_k$ 由 2×2 线性方程 | B6.2 |
| $pT$ | $\mu=\theta=\infty$ | $\sum m_ke_k(p,T)=\rho e$，$\sum m_kv_k=1$ | B6.3 |
| 有限 $\theta$ | $\mu=\infty$ | $p$ 闭合 + 后向 Euler 传热 $\Delta Q$ | B6.4 |
| 有限 $\nu$ | $\mu=\theta=\infty$ | $\rho Y-m^{\rm tr}_2+\tau\nu(g_2-g_1)=0$ | B6.5 |
| $pTg$ | 全瞬时 | $T=T_{\rm sat}(p)$，$Y$ 由能量，$p$ 由体积 | B6.6 |

速率 $\to\infty$ 时各闭合连续地趋于平衡闭合（AP，定理 1）；速率 $=0$ 时逐位退化为冻结格式（定理 2）。所有根都带区间求解。

---

## E6. 显式参考程序的路径（A3–A8）

1. 步长：声学 CFL 0.4。
2. Heun 每一级：原始变量 $(\alpha_1,\rho_1,\rho_2,u,v,p_1,p_2)$ 的 minmod 重构 → HLLC（相能量用各相 Hugoniot，MG 相用 (A2.8)）→ SPB 的 Godunov 更新（A5）→ 瞬时松弛投影（每格的统一闭合，$s=0$）。
3. 松弛起点：输运态不可容许的相（$p_k+\pi_k\le0$ 或 MG 的 $c^2\le0$）先按 SPB §3.4 取混合压力下的能量（差值进 $\delta$），再做闭合（A6，V13）。
4. 有限速率松弛每步只在 Heun 平均之后做一次（$\tau=\Delta t$，A3，V6）。
5. 拒绝与失败处理同 IMEX。

---

## E7. 边界、初值、输出、并行（B11）

- 每一边单独：透射、固壁（法向速度 $2u_w-u$，可带活塞速度）、周期（`bc_apply`，两个程序共用）。
- 初值：圆形区域 $4\times4$ 子点平均。
- 每个算例的区域、状态方程、各边边界、终止时间、论文网格、松弛档位都在 `test_setup` / `test_state`；切换算例改 `Test test = ...;` 或 `case=N`（清单见程序顶部）。
- OpenMP：逐格循环与 BiCGStab；结果与线程数无关，一维输出与串行逐字节相同。
- 输出列：`x y p u v rho alpha2 Y2 rho1 rho2 c_f c_W T1 T2 p1 p2 g1-g2`（`prec=` 有效位数，默认 16）。

---

## E8. 这条路径上的已知限制（C 页“仍未解决”）

光滑区收敛阶（IMEX 约 1.2–1.4）；SPB 4.3.2 的尾部误差；±500 m/s $pTg$ 的起步加热；十二烷微量汽相的拒绝；SPB 4.5 扇区误差；$M\le10^{-5}$ 的双精度底；活塞为固定网格近似。
