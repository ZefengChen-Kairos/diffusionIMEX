# 论文叙事审查：卖点和实际贡献的错配

来源：作者 2026-11-27 提出的五条批评。本文件逐条对着 `.tex` 源码核实，
记下我同意/修正的地方、建议的重写方向，以及两个能把悬而未决的问题一次测掉的实验。

**总判断：五条基本都成立。其中第 2 条和第 4 条是我上一轮逐式复核时应该抓到却没抓到的
——我核了每个公式对不对，没核故事讲不讲得通。**

---

## 1. 双参数 ε：卖错了点

### 核实

卖时间步的地方，三处：

| 位置 | 原文 |
|---|---|
| Highlights L211 | "time step freed from the acoustic and Alfvén constraints where they are stiff" |
| §3.1 L1082–1084 | "either constraint may be removed on its own: the time step is freed from the acoustic restriction wherever $M\ll1$ and from the Alfvénic one wherever $M_A\ll1$" |
| 结论 L3407 | "The time step is governed by the material velocity wherever the corresponding Mach number is small" |

ε 的二值性，实测（$\alpha=14$，$M_0=0.15$，$M_1=0.4$）：

| $M$ | 0.07 | 0.15 | 0.30 | 0.40 | 0.60 | 1.00 |
|---|---|---|---|---|---|---|
| $\varepsilon_M$ | 1.1e-16 | 2.9e-12 | **0.430** | 0.9992 | 0.99999 | 1 |

过渡带 $[0.15,0.4]$ 里确实有中间值（$M=0.3$ 时 0.43），但带很窄，两侧各自贴死 0 和 1。
"continuous adaptation" 偏大，"a smooth but narrow switch" 才是实情。

### 我的修正：时间步的论断要看和谁比

批评里说"Δt 只会相等或更小"——这取决于单参数由什么驱动，必须分情况：

| 单参数由什么驱动 | 结果 |
|---|---|
| $\max(M,M_A)$（保守选择） | 任一通道刚性就全隐式 → 显式谱永远 $\mathcal O(\|u\|)$ → **双参数永不占优** |
| 快磁声马赫数 $M_f=\|u\|/c_f$ | $c_f\ge\max(a,c_B)$ 所以 $M_f\le\min(M,M_A)$ → 同上，**永不占优** |
| **$M$ 单独驱动**（气体动力学格式的朴素推广，即 [32–35] 的实际做法） | 超声速但强磁化时（$M\gtrsim1$，$M_A\ll1$）单参数把 $\mathbf F^B$ 留在显式 → $\Delta t\sim\Delta x/c_B$ 极小；双参数 $\varepsilon_B=0$ → $\Delta t\sim\Delta x/(\|u\|+a)$。**双参数大幅占优** |

所以正确的说法是**两面的**，而且比现在写的更有力，因为它点名了真正的竞争对手：

> 相对于由两个马赫数中更严格者驱动的单参数格式，双参数不买到时间步，买到的是
> **非刚性通道重新获得迎风激波捕捉，即分辨率**。相对于仅由声学马赫数驱动的单参数
> 格式（即现有工作对 MHD 的朴素推广），双参数在磁刚性占主导的区域买到的是**数量级的
> 时间步**。

### 论文现有的证据不支持这个论断

Fig. 10 量的是 $\Delta t/\Delta t_{\rm MHD}$，分母是**全显式**。那是最容易的对比。
Fig. 1/2 里的 "Standard IMEX" 是 $\varepsilon_M=\varepsilon_B\equiv0$ 的固定选择，
也不是单参数自适应格式。

**缺的正是 dual-vs-single-adaptive 的那条曲线。**

### 实验 E1（不需要改代码）

Field loop：$M\approx5.5\times10^{-3}$（声学刚性），$M_A\approx7.9\times10^3\ldots7.9$
（磁非刚性）。单参数由 $M$ 驱动会给出 $\varepsilon=0$，即**磁通量也全隐式**。

```bash
./mhd cfl=0.2 tag=dual                 # 双参数：eps_M≈0, eps_B=1
./mhd cfl=0.2 epsBmax=0 tag=single_M   # 模拟"单参数由 M 驱动"：强制 eps_B=0
```

`epsBmax=0` 把 $\varepsilon_B$ 钉死在 0，正好就是单参数格式在这个算例下的行为。
比较环的形状保持、$\|B\|$ 的扩散量、以及 $\Delta t$。这一张图能把第 1 条的两面论断
同时立住。

---

## 2. 特征结构：数学最硬，因果链断了

### 核实：确实是直接自相矛盾

```
L320–322 (Intro):
  "The decisive structural consequence is the restoration of strong hyperbolicity:
   whenever eps_M > 0 the explicit Jacobian has eight real eigenvalues, a complete
   set of right eigenvectors... This makes HLLEM applicable"

L1119–1122 (§3.1):
  "The Riemann solver does not require eps_M > 0, since it uses only the regularised
   Alfven projection and the classical entropy strength (Section 3.2). The
   characteristic basis does: ... evaluated with eps_M bounded from below by 10^-2"
```

一个说定理让 HLLEM 可用，另一个说 HLLEM 不需要它。审稿人读一遍就能看出来。

### HLLEM 里到底用了多少 split 结构

| HLLEM 的组成 | 用的是什么 |
|---|---|
| 波速包络 (52) | **split** 四次式的极端根 ✓ |
| 两道 Alfvén 波 (56)(57) | **split** 2×2 块 (41) + Corollary 2.4 的正则化特征对 ✓ |
| 熵波 (58) | **classical** $\Delta\rho-\Delta p/\tilde a^2$ ✗ |
| $\delta_k$ (59) | split 特征值 ✓ |

熵振幅被 $\delta_k$ 的结构逼回了经典形式——§3.2 自己写了原因："the anti-diffusive
term of (54) weights the fields with different $\delta_k$, so that the cancellation
is lost"。

所以"split 特征结构让 HLLEM 适用"要打折：**真正用上 split 结构的是两道 Alfvén 波
和波速包络**，熵波没有。

### 真正的亮点被埋了

Lemma 2.3(iii) + Corollary 2.4 合起来说的是：

> Alfvén 波张成二维不变子空间 $S_A=\mathrm{span}\{\mathbf e_v,\mathbf e_B\}$，
> 其对偶泛函 $\ell_v(\delta\mathbf Q)=\beta_z\delta v-\beta_y\delta w$ 和
> $\ell_B(\delta\mathbf Q)=\beta_z\delta B_y-\beta_y\delta B_z$ **在不含于 $S_A$ 的
> 每一个特征向量上取零**。因此 Alfvén 投影有解析闭式，且在 $B_\perp\to0$ 和
> $\varepsilon\to0$ 两个极限下都正则——**对任意 $\varepsilon$ 成立，包括 $\varepsilon=0$**。

这才是结构部分的主角：**完整特征基在刚性极限下退化，但 Alfvén 扇区不退化**，
所以在全特征基失效的地方仍然能做一个分辨波的 Riemann 求解器。

现在它是一个 Lemma 的第三款加一个 Corollary。应该升级成命名的结果，并且让 Intro 的
因果链改成：

> ~~强双曲 → HLLEM 可用~~
> **定理 2.1（完整特征向量，$\varepsilon_M>0$）→ 特征 WENO 重构可用（且需 $\varepsilon_M\ge10^{-2}$）**
> **Lemma 2.3 + Corollary 2.4（$S_A$ 的无奇异解析投影，任意 $\varepsilon$）→ HLLEM 可用**

两条腿各管各的，矛盾消失，而且第二条比第一条更不平凡。

---

## 3. CT 写进压力算子："无迭代"是借来的

### 核实：论文自己已经写了，只是摘要没认

§3.4.4 原文：

> "The constrained transport, the induction of $B_z$ and the magnetic energy flux are
> not part of $\underline{\mathbf T}_2$; they enter only through the defect."

而 §3.4.3 收尾说 "without any iteration between the two sub-systems"。字面成立
（B 和 p 之间确实没有外层 Picard），但 defect correction 本身就是迭代，而且 CT、$B_z$、
磁能通量全部滞后在里面。

实测（1000²，CFL 0.2，1318 步 × 4 级 = 5272 级）：

```
pressure solves 5272        →  每级 1 次压力系统 ✓
defect-correction sweeps 28066  →  每级 5.32 次 sweep
```

**所谓"无迭代"，是把外层迭代挪进了压力求解内部。** 这不是错，但摘要和 Highlights
的措辞要认这件事。（"two linear solves per stage" 已在本轮删除。）

低 β 下 defect correction 的收敛性是这一点的命门，结论里已经诚实写了
"supported here by numerical evidence only" ✓。

三维那句话要软化：2D 用标量 $A_z$，3D 要换向量势并处理规范，不是"same construction"。

---

## 4. 最尖锐的一条：两个创新点互相抵消

### 核实：论文自己把证据写在 L1671

```
L1669–1671:
  "Since (73) is the full induction equation, the field obtained in this way
   replaces the in-plane components produced by the explicit update, and
   no weighting by eps_B is involved."
```

### 这意味着什么

在 2D，面内 $\mathbf B_\perp$ 的终值 $=\mathrm{curl}_4A$，而 $A$ 按**物质速度**迎风输运，
整个过程不带 $\varepsilon_B$。所以：

| $\varepsilon_B$ 还管什么 | |
|---|---|
| 动量方程里的磁应力 (67)(70)(94) | ✓ |
| $B_z$ 的感应 (78)(96) | ✓ |
| 磁能通量 (79)(97) | ✓ |
| 磁预测步 $\tilde{\mathbf B}$ (69) | ✓ |
| 时间步 | ✓ |
| **面内感应方程** | ✗ **被 CT 整个覆盖** |

而 HLLEM 的 Alfvén 反扩散作用在显式通量的横向速度和横向磁场上——**横向磁场随后被
$\mathrm{curl}_4A$ 覆盖掉**。所以在 2D，Alfvén 反扩散对面内场只剩一条间接路径：
通过 $v_\perp$ 影响 $A$ 的输运速度。

1D 测试里没有 CT（$B_x$ 是常数，横向场按通量形式更新），那里展示的 Alfvén 锐利度
**不能直接推到 2D 的面内场**。

### 实验 E2（决定性，不需要改代码）

Field loop 是纯面内场输运，$\varepsilon_B=1$，Alfvén 分辨率应该最该起作用的地方：

```bash
./mhd cfl=0.2 tag=hllem_on              # 默认：带 delta_k 反扩散
./mhd cfl=0.2 hllemoff=1 tag=plainHLL   # 关掉全部反扩散，退回 plain HLL
```

比较两者的 $B_x,B_y$：

- **若逐点相同到舍入** → Alfvén 反扩散在 2D 根本到不了面内场，论文里关于 Alfvén
  分辨率的卖点必须限定在 1D 和 $B_z$，2D 面内场另说。
- **若有可见差异** → 机制经由速度真实存在，可以量化并写成一句有依据的话。

这个实验的结果决定第 4 条是"改写措辞"还是"改设计"。

### 如果结论是"到不了"，设计上有一条出路

现在 $A$ 的输运是非守恒形式 $\partial_tA+\mathbf u\cdot\nabla A=0$，里面没有地方放
$\varepsilon_B$。但 Rossmanith 的方法也支持用 **EMF 形式**驱动 $A$——如果 EMF 取自
HLLEM 的数值通量（它带 $\varepsilon_B$ 和 Alfvén 反扩散），那么两个创新点就不再抵消，
Alfvén 分辨率会真正进入面内场。

**这是目前最值得先想清楚的设计问题**，而且正好和手上的零守卫改造在同一片代码里。

---

## 5. 三阶是规格不是创新 —— 而且先行工作没引

### 核实

```
grep -c 'ChenWuXiong2023' paper.tex  →  0
```

`.bib` 里躺着：

> Chen, Wu, Xiong, *High order asymptotic preserving finite difference WENO schemes
> with constrained transport for MHD equations in all sonic Mach numbers*,
> JCP 488 (2023) 112240.

**高阶 + AP + 无散 + MHD，四个属性齐全，一次都没引。** 这既坐实了"三阶是规格"，
也是一个会被直接点名的遗漏。

Proposition 3.1（MP3 在光滑极值处保持设计阶）是真结果，但它是 Suresh–Huynh 限制器的
性质，不是这个格式的性质，属于附带产物。

---

## 建议的重写骨架

把贡献从三条改成三条，但换掉内容：

| | 现在卖的 | 应该卖的 |
|---|---|---|
| 1 | 双参数 → 时间步自由 | 双参数 → **非刚性通道恢复迎风分辨率**；对声学驱动的单参数，在磁刚性区**另有数量级时间步** |
| 2 | 强双曲 → HLLEM 可用 | **$S_A$ 在任意 $\varepsilon$ 下的无奇异解析投影** → 在全特征基退化处仍能分辨波；定理 2.1 单独服务于特征重构 |
| 3 | CT 嵌入压力算子 → 每级两次线性解、无迭代 | CT 嵌入压力算子 → **解决"CT 覆盖 B 之后能量守恒与 E–p 一致二选一"的老问题**；代价是 defect correction 承担了迭代，低 β 收敛性目前只有数值证据 |
| 4 | 三阶（作为卖点） | 三阶（作为**规格**），并引 Chen–Wu–Xiong 2023 |

另外必须补的一句话（第 4 条）：在 2D，面内感应由 CT 承担，$\varepsilon_B$ 管的是
Lorentz 力、磁能通量、$B_z$ 和时间步；Riemann 求解器的 Alfvén 分辨率经由输运速度
进入面内场——**这句话的强弱取决于 E2 的结果**。

---

## 待办

| | 内容 | 依赖 |
|---|---|---|
| E1 | field loop：`epsBmax=0` vs 默认，立"分辨率 vs 时间步"的两面论断 | 无，可直接跑 |
| E2 | field loop：`hllemoff=1` vs 默认，测 Alfvén 反扩散在 2D 到不到得了面内场 | 无，可直接跑 |
| W1 | 重写 Intro 的因果链（第 2 条），把 $S_A$ 提成命名结果 | 无 |
| W2 | 摘要/Highlights/§3.1/结论：分辨率优先于时间步（第 1 条） | E1 更好 |
| W3 | 认 defect correction 是迭代；软化 3D 那句（第 3 条） | 无 |
| W4 | 引 Chen–Wu–Xiong 2023，三阶降格为规格（第 5 条） | 无 |
| W5 | 2D 面内场由 CT 承担那句话（第 4 条） | **必须等 E2** |

---

## 10. 第 4 条批评的代码级判定（2026-11-27）

用户的批评 4 原话：「两个贡献在 2D 里部分抵消，因为 CT 覆盖了面内场」。
**方向对，机制说错了，而且改对之后对论文更有利。** 以下每条都在
`blastwave_debug/patch/Final_emb_omp_patched.cpp` 里核过，给行号。

### 10.1 已核实的代码事实

1. **`Unstaggered_CT` (4163–4308)** 把矢势 z 分量 $A$ 按**纯输运方程**推进
   $$\partial_t A + u\,\partial_x A + v\,\partial_y A = 0,$$
   用 HJ-WENO5-Z 对 $\nabla A$ 做上风重构，**显式**处理。源码注释 4195–4199
   自己写明了理由：*"A is advected at the MATERIAL speed, so an explicit
   treatment costs nothing: the Alfven stiffness sits in the Lorentz force,
   not here."*

2. 紧接着 **赋值** $B_x,B_y=\mathrm{curl}_4 A$（4300–4306，四阶中心差）：
   ```cpp
   u_con(i,j,5) =  ( -A(i,j+2,2) + 8.0*A(i,j+1,2) - 8.0*A(i,j-1,2) + A(i,j-2,2) ) / (12.0*dy);
   u_con(i,j,6) = -( -A(i+2,j,2) + 8.0*A(i+1,j,2) - 8.0*A(i-1,j,2) + A(i-2,j,2) ) / (12.0*dx);
   ```
   是赋值不是累加，于是 $\mathrm{div}_4(\mathrm{curl}_4 A)\equiv 0$ 逐点成立。

3. **`Magnetic_Sweep_FV` (6007+)** 打开 CT 时**只**更新能量行 (4) 和 $B_z$ 行 (7)。
   面内感应通量和第 5、6 行的更新被 `if(!g_use_ct)` 守卫（6023、6031），
   即**只有 1D 才算、才施加**：
   ```cpp
   if(!g_use_ct){ FY(i,j,1) = v*b5 - u*b6; FX(i,j,2) = u*b6 - v*b5; }
   if(!g_use_ct && g_Bx_st){ u_con(i,j,5) = ...; u_con(i,j,6) = ...; }
   ```

4. `ct_now()` 在 6158（star 态）、6395（每个 Picard/embedded 迭代）、6539（最终）
   都调用，`emb_eval` 内部也调用——**每一次状态成形，$B_\perp$ 都被 $A$ 重新赋值一遍**。

结论：2D 里 $B_x,B_y$ 的有限体积磁通量**整条被丢弃**，不是"部分抵消"。

### 10.2 这对 $\varepsilon_B$ 的 claim 意味着什么

$\varepsilon_B$ 加权 $F^B$ 的显式份额。$F^B$ 出现在：

| 行 | 2D 里是否被 $\varepsilon_B$ 分裂 |
|---|---|
| 动量（Lorentz 力 $-\nabla\!\cdot\!(BB)+\nabla(\lvert B\rvert^2/2)$） | 是 |
| 能量 | 是（`eps_B_grid_x/y` 是 `Magnetic_Sweep_FV` 的实参） |
| $B_z$ | 是 |
| 面内感应 $B_x,B_y$ | **2D 里根本不存在这一项** |

所以正确的说法不是"少分裂了一块"，而是：

> **$\varepsilon_B$ 作用在 Alfvén 刚性真正出现的每一行上。** 势函数形式下的面内感应
> 方程特征速度只有 $\lvert u\rvert$，**本身不含 Alfvén 刚性**，因此无需分裂，
> 全程显式推进即可。

这把第 4 条从"设计缺陷"翻成了**一致性论证**，而且读者能自己验证——
$\partial_t A+u\cdot\nabla A=0$ 的特征速度就是 $\lvert u\rvert$，一眼可见。
顺带也解释了为什么势函数 CT 在这个 IMEX 框架里是**必需**的而不是可选的：
交错 CT 要把面 $E_z$ 拆成三个 $\varepsilon$ 份额做通量对偶记账（源码注释
4180–4184 明确提到 "the flux-dual bookkeeping needed by the staggered variant
does not arise here"），势函数形式把这个记账问题整个消掉。

### 10.3 这对 HLLEM 反扩散的 claim 意味着什么

$\delta_A$ 反扩散只活在**显式** HLLEM 通量里。于是

* **1D**：直接作用在 $B_y,B_z$ 上。
* **2D**：作用在动量、能量、$B_z$ 上；到面内场只能**经由 $u$**
  （$\delta_A$ 改动量 → $u$ → HJ-WENO5 输运 $A$ → $B_\perp=\mathrm{curl}_4A$）。

因此论文**不能**直说"Alfvén 反扩散锐化 2D 面内场的旋转间断"。必须写成：
锐化 Riemann 解算器所推进的那些行里的 Alfvén 波；面内场经由速度场继承这一效果。
E2（`hllemoff=1` vs 默认）量化的就是这条间接通路的强度——
若 $B_x,B_y$ 逐点相同到舍入，措辞必须限定在 1D 与 $B_z$。
