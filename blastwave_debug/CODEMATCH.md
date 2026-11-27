# 论文 ↔ 代码 逐条核对表

核对对象：`blastwave_debug/patch/Final_emb_omp_patched.cpp`（= 原始
`Final_emb_omp.cpp` + 修复 A–H）。**我的修复 A/C/D 不在论文里，不参与比对**；
凡用到它们的行会标注。论文行号指
`scratchpad/paper/paper.tex`（即 `cas-sc-template.tex`）。

状态记号：✓ 一致 · ✗ 不一致 · ⚠ 一致但措辞需要补强 · ○ 未核

---

## 1. 分裂参数 $\varepsilon_M,\varepsilon_B$

| 项 | 论文 | 代码 | |
|---|---|---|---|
| $M_0$ | $0.15$ (1153) | `M0 = 0.15` (112) | ✓ |
| $M_1$ | $0.4$ (1153) | `M1 = 0.4` (112) | ✓ |
| $\alpha$ | $14$ (1153) | `alpha = 14` (112) | ✓ |
| 形式 | $\varepsilon=1-s(M)$，含 $(1-M_0^\alpha)$ | 同 (112) | ✓ |
| 面值 | 取相邻单元最大 | `face value = max of neighbours` (112) | ✓ |
| 物理地板 | 通量不设地板 | `eps_M_min = 0` (176) | ✓ |
| 重构地板 | 「基底远离零」 | `g_eps_basis_min = 1e-2` (fix 3) | ⚠ 论文未给数值 |

**数值验证**：field loop ($A_0=1$, 128×64) 的输出场给出
$\varepsilon_M\equiv0$、$\varepsilon_B\equiv1$。解析核对：$M=|u|/c_s=2.236/408=0.0055$
→ $M^{14}\approx10^{-31.6}$ 下溢为 0；$M_A=|u|/c_B=2.236/0.282=7.93>M_1$ → $\varepsilon_B=1$。
**这个算例同时实现 $\varepsilon_M=0$ 与 $\varepsilon_B=1$**，正是导言所说单参数无法表示
的混合区（$\beta=2.5\times10^6$ 高、$M$ 低、超 Alfvén）。

## 2. CFL 条件

| 项 | 论文 | 代码 | |
|---|---|---|---|
| 只用分裂谱 | `eq:cfl` (1112) | `caldt` 用 $\tilde\lambda^c$ | ✓ |
| 两约束可独立解除 | 1118–1123 | $\varepsilon_M,\varepsilon_B$ 独立进入 | ✓ |

## 3. HLLEM 反扩散

| 项 | 论文 | 代码 | |
|---|---|---|---|
| $\delta_k=1-\tilde\lambda_k^-/s_L-\tilde\lambda_k^+/s_R$ | `eq:deltak` | 1468–1470 | ✓ 逐字 |
| clamp 到 $[0,1]$ | 「restricted to $[0,1]$」 | `HLLEM_clamp01` (1229) | ✓ |
| $\lambda^\pm$ 定义 | 正/负部 | `0.5(lam∓fabs(lam))` (1227–1228) | ✓ |
| 波族 $k\in\{\rm ent,A^\pm\}$ | `eq:hllem` (1315) | `delta_C, delta_Ap, delta_Am` | ✓ |

## 4. 激波激活稳定化

| 项 | 论文 | 代码 | |
|---|---|---|---|
| $\kappa_1$ | $0.4$ (`eq:chi` 后) | `shock_kappa1 = 0.4` (fix H) | ✓ 修复后 |
| | | 原始码为 `0.05` | ✗ **原始码与论文差 8 倍** |
| $3\times3$ 取最大 | `eq:chi` | `max_{|k|,|l|<=1}` | ✓ |
| $c^{\min}$ 为 $3\times3$ 最小 | ✓ | ✓ | ✓ |

## 5. 磁子系统

| 项 | 论文 | 代码 | |
|---|---|---|---|
| 动量代入感应 | `eq:mag_mom` | `PreMagSubsystem_Pressure_HF` | ✓ |
| 磁压切线线性化 | `eq:mag_mom` | 同 | ✓ |
| 一次 GMRES | 「solved once with GMRES」 | 头注 104「one GMRES solve」 | ✓ |
| $\tilde B$ 只用于 Lorentz 力 | `eq:mag_inter` | ✓ | ✓ |
| 稳定化 Laplacian | `eq:mag_stab` | `Magnetic_Sweep_FV` 内 | ✓ |

## 6. 约束输运 (CT)

| 项 | 论文 | 代码 | |
|---|---|---|---|
| $\partial_tA+u\cdot\nabla A=0$ | `eq:Aadvect` | 4168–4171 | ✓ |
| HJ-WENO5-Z 上风 | `eq:Aupdate` | ✓ | ✓ |
| 四阶 curl | `eq:curl4` | 4304–4305 逐字 | ✓ |
| $B_\perp$ **替换**显式更新 | 「replaces the in-plane components」 | 赋值，非累加 | ✓ |
| 不加 $\varepsilon_B$ 权 | 「no weighting by $\epsB$」 | 面内感应被 `if(!g_use_ct)` 守卫 (6023, 6031) | ✓ |
| 初场由 `eq:curl4` 定义 | ✓ | 6828–6838 | ✓ |
| **为什么这是对的** | 原本缺 → 已补（一致性论证） | — | ⚠→✓ |

见 `STORY.md` 第 10 节：势函数方程特征速度只有 $|u|$，不含 Alfvén 刚性，
故无需分裂；刚性全在 Lorentz 力里，而那一项**是**被 $(1-\varepsilon_B)$ 送进隐式的。

## 7. 压力子系统

| 项 | 论文 | 代码 | |
|---|---|---|---|
| `eq:p_mom`,`eq:Bz`,`eq:Umag` | ✓ | ✓ | ✓ |
| `eq:klin` 动能/磁能线性化 | ✓ | `kemb_lin` | ✓ |
| $\Phi$ 对 $m$ **仿射** | 1786–1790 | 头注 6136「Phi is AFFINE in m」 | ✓ |
| 单一线性压力方程 | `eq:pwave` | 「is LINEAR in p: one GMRES solve」 | ✓ |
| **每 stage 一次压力解** | 「without a second pressure solve」 | `ipic < (emb_on ? 1 : g_npicard)` (6393)，`emb_on` 恒真 | ✓ |
| | | `g_npicard = 3` 在默认配置下是**死代码** | ✓ 不影响 |
| 左端每次求值算无散场 | 1795–1798 | `emb_eval` 调 `Unstaggered_CT` | ✓ |
| 算子非对称 | ✓ | ✓ | ✓ |

## 8. 压力保正（嵌套 Newton）

| 项 | 论文 | 代码 | |
|---|---|---|---|
| $p_\varepsilon$ | $10^{-7}$ | `p_floor = 1e-7` (221) | ✓ |
| $\theta$ | $0.9$ | `theta_floor = 0.9` (223) | ✓ |
| 二阶 M 矩阵 $T_2$ | 紧致五点 + 对流上风 | ✓ | ✓ |
| CT/$B_z$/磁能**不在** $T_2$ 内，只进 defect | 1864–1867 | 头注 105–106 逐字 | ✓ |
| defect correction `eq:defect` | ✓ | ✓ | ✓ |
| $M$ 大斜率取值 | 未给（只说 steep） | `Mbig_fac=1e9` + Cor.1 下界 | ✓ 无冲突 |

> §sec:pressure 里「左端每次求值都算无散场」与「CT 由 defect 承担而非由被求逆的
> 算子承担」**不矛盾**：被求逆的是预条件子 $T_2$，整个左端（含 CT）是被求解的算子。
> `eq:defect` 把这件事写清楚了，前向引用有效。

## 9. 时间离散 SA-DIRK(3,4,3)

论文给公式，代码给硬编码数值，两边 round-trip 一致，阶条件全部满足。

| 量 | 论文 | 实算 / 代码 | |
|---|---|---|---|
| $\delta$（中间根） | $0.4358665215$ | $0.4358665215084597$ | ✓ |
| $c_3=(1+\delta)/2$ | ✓ | $0.7179332607542298$ | ✓ |
| $b_2=-(6\delta^2-16\delta+1)/4$ | ✓ | $1.2084966491760119$ | ✓ |
| $b_3=(6\delta^2-20\delta+5)/4$ | ✓ | $-0.6443631706844715$ | ✓ |
| $\tilde a_{32}=c_3-\tilde a_{31}$ | $\tilde a_{31}=1.4377449738$ | 代码 `-0.71981171300978564`，反推 $\tilde a_{31}=1.4377449737640$ | ✓ |
| $\tilde a_{43}=\tfrac12-\tilde a_{41}$ | $\tilde a_{41}=0.9169932984$ | 代码 `-0.41699329835202015`，反推 $\tilde a_{41}=0.9169932983520$ | ✓ |

阶条件（共享 $b$）：

| | 实算 | 应为 |
|---|---|---|
| $\sum b_i$ | $1.0000000000000000$ | $1$ ✓ |
| $\sum b_ic_i$ | $0.5000000000000003$ | $1/2$ ✓ |
| $\sum b_ic_i^2$ | $0.3333333333333335$ | $1/3$ ✓ |
| $\sum b_i\tilde a_{ij}c_j$ | $0.1666666666617603$ | $1/6$ ✓（残差 $5\times10^{-12}$ 来自论文把 $\tilde a_{31}$ 截到 10 位） |
| $\sum b_ia_{ij}c^I_j$ | $0.1666666666666667$ | $1/6$ ✓ |
| 两表行和 $=c$ | ✓ | ✓ |
| stiffly accurate（$A_I$ 末行 $=b$） | True | ✓ |

两个负元是论文公式的**直接后果**，不是笔误：$c_3-\tilde a_{31}<0$ 且
$\tfrac12-\tilde a_{41}<0$。显式表带负元不影响相容性与三阶，但它意味着
显式部分不是 SSP——任何基于「显式系数非负」的保正论证在这里都不能用。
（这正是 fix G 的 IDP 构造必须把低阶通量从 **stage 的 cell averages** 重算、
而不能依赖凸组合系数的原因。）

## 10. E1 / E2：两个决定性实验（2026-11-27）

Field loop，$A_0=1$，$128\times64$，CFL 0.2，senv=2，跑满一个周期 $t_f=1$。
精确解 = 初值（周期域、$u=(2,1)$、域 $2\times1$，正好一周期）。
脚本 `blastwave_debug/e2cmp.py`，驱动 `blastwave_debug/loop2.sh`。

实测 $\varepsilon_M\equiv0$、$\varepsilon_B\equiv1$（VTK 输出场），
与论文 §Field loop advection 所述一致（$M\approx5.5\times10^{-3}$、$M_A\approx7.9$）。

| | base ($\varepsilon_B=1$，$\delta_k$ 开) | noad ($\delta_k$ 全关) | epsB0 ($\varepsilon_B\equiv0$) |
|---|---|---|---|
| 步数 | 1155 | 1155 | **976** |
| $\Delta t/\Delta t_{\rm expl}$ 中位 | 229.69 | 229.69 | **272.75** |
| 墙钟 | 10.4 min | 10.4 min | **8.9 min** |
| $E_{\rm mag}(1)/E_{\rm mag}(0)$ | 0.922029 | 0.922029 | 0.922105 |
| $\max\lvert B_\perp\rvert$ 保留 | 0.960036 | 0.960036 | 0.960384 |
| $L_1$ 误差 $B_x$ | 8.8176e-03 | 8.8176e-03 | 8.8215e-03 |
| $L_1$ 误差 $B_y$ | 9.4006e-03 | 9.4006e-03 | 9.3868e-03 |

逐点差（舍入量级为 $10^{-14}$）：

| | $\max\lvert\Delta\rvert/\max\lvert{\rm base}\rvert$ $B_x$ | $B_y$ | $u$ | $v$ |
|---|---|---|---|---|
| noad vs base | 1.21e-06 | 1.43e-06 | 1.31e-07 | 2.48e-07 |
| epsB0 vs base | 2.67e-03 | 2.53e-03 | 4.36e-04 | 9.96e-04 |

### E2 判定（HLLEM 反扩散到不到得了 2D 面内场）

**到得了，但 $10^{-6}$，比离散误差低四个数量级。** 不是舍入（$10^{-6}\gg10^{-14}$），
所以间接通路 $\delta_k\to$ 动量 $\to u\to A\to B_\perp$ 真实存在；但质量指标逐位相同。
→ **改措辞，不改设计。** 论文不能声称 Alfvén 反扩散锐化 2D 面内场的旋转结构；
该声称 1D（$B_y,B_z$ 本身是有限体积行，RP7 即孤立旋转间断）、2D 的 $B_z$、动量/能量行。
已写入 §sec:CT。

限制：这个 field loop 是纯输运，**本来没有旋转间断**，所以它测的是"到不到得了"，
不是"真有旋转间断时能锐化多少"。后者由结构论证封死（面内通量被丢弃，唯一通路是 $u$）。

### E1 判定（单参数 vs 双参数）

**单参数（安全选择）在这里步长大 19%、墙钟快 14%、精度不可区分**
（四个误差指标里两个单参数更好，差异都在第四位）。
所以论文原来那句「the admissible time step is the same ... what the second
parameter buys there is *resolution*」**不成立**。

原因是定义的后果，不是偶然：

$$\varepsilon_M\to0\Rightarrow\lvert u\rvert\ll a,\qquad
\varepsilon_B\to1\Rightarrow c_B\lesssim\lvert u\rvert
\;\Longrightarrow\; c_B\ll a\;\Longrightarrow\;\beta\gg1.$$

**只要分裂把磁通量设显式而声学设隐式，磁场必然是高 $\beta$ 弱扰动**
（此例 $\beta=2.5\times10^6$，$\lvert B\rvert^2/2=0.04$ 对 $p=10^5$）。
再加上 2D 面内场不走有限体积通量，$\varepsilon_B$ 对 $B_\perp$ 的影响只能经由 $u$，
而 $u$ 在高 $\beta$ 下由声学主导。

反向区 $\varepsilon_M\to1,\varepsilon_B\to0$ 必然 $\beta\ll1$，磁通量主导——
**这才是 $\varepsilon_B$ 值钱的地方**，而且 blast wave 里这个环带占 41% 的网格
（$\varepsilon_M\approx1$、$\varepsilon_B\approx0$，$r\in[0.139,0.463]$，
快激波在 $r=0.288$）。

→ 导言已改成**不对称**的故事：$\varepsilon_B$ 的收益在强磁化区（步长，数量级），
在高 $\beta$ 区反而要付约 19% 的步长，上界为因子 2。

## 11. blast wave 的 $\varepsilon$ 分布（256², $\kappa_1=0.4$, $t=0.01$）

| | $\varepsilon_M$ | $\varepsilon_B$ |
|---|---|---|
| 恰为 0 | 45.95% | 77.66% |
| 恰为 1 | 38.45% | 2.56% |
| 均值 | 0.4603 | 0.0450 |
| 两者相等 | 47.60% | |
| $\lvert\varepsilon_M-\varepsilon_B\rvert>0.5$ | **41.26%**（该处 $\varepsilon_M$ 均值 0.988、$\varepsilon_B$ 均值 0.013） | |

已写入 §Ideal MHD blast wave（替掉 `% CHECK: describe eps_M and eps_B`），
并在 .tex 里留了注释说明这些数来自 256²，投稿前需用 1000² 生产运行刷新。

## 12. 特征结构（§Eigenstructure）

| 断言 | 论文 | 验证 | |
|---|---|---|---|
| 把 $F^p+F^B$ 全给隐式 → 显式部分只弱双曲 | 导言 304–307 | `jordan_check.txt`：$F^c$ 的特征值 $u$ 代数重数 4、几何重数 3 → **defective** | ✓ |
| | | $F^B$ 单独：特征值 0 代数 4、几何 3 → defective | ✓ |
| 分裂系统谱为实 | `eq:ms_eig` 四次方程 | `hyperbolicity_check.py`：**0 / 200000** 随机态有复根（相对虚部 > 1e-8），含 $\varepsilon_M\ne\varepsilon_B$ 与 $\varepsilon_M=0$ | ✓ |
| 基底在 $\varepsilon_M\to0$ 退化 | 「the closed-form basis degenerates as $\epsM\to0$」 | 扫描在 $\varepsilon_M=\varepsilon_B=0$ 处仍为实根（实但 defective），与之一致 | ✓ |

## 13. 实测 $\Delta t/\Delta t_{\rm explicit}$（中位数）

直接回答 .tex 里结论段那条 `% CHECK: quote the measured ratio to the explicit step`。

| 算例 | 网格 | 中位 | 最大 | 说明 |
|---|---|---|---|---|
| field loop $A_0=1$ | 128×64 | **229.69** | 229.74 | 真低 Mach，$M=5.5\times10^{-3}$ |
| field loop $\varepsilon_B\equiv0$ | 128×64 | 272.75 | 273.02 | 单参数代理，比双参数**大 19%** |
| blast wave | 256² | **1.084** | 4.409 | 激波主导 |
| rotor | 256² | **1.000** | — | 盘内超声速超 Alfvén → $\varepsilon_M=\varepsilon_B=1$，全域最大值取到完整 MHD 谱 |
| RP5 | 800×1 | **1.000** | 1.0005 | 1D 激波管 |

**格局很清楚：所有激波主导的算例收益基本为零，收益只在真正的低 Mach 算例。**
这是「resolution not time step」改写的直接数值依据。

## 14. 待核

- ○ §Eigenstructure of the explicit sub-system（已在 `REVIEW_IMEX_MHD-3.md` 逐式推过）
- ○ §Time discretisation（`eq:88` 表已核；其余未核）
- ○ §Finite volume operators / WENO3-Z 重构
- ○ 各算例初值与论文表格
