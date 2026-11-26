# IMEX_MHD-3.pdf 逐式复核报告

对 `IMEX_MHD-3.pdf`（39 页）做了一遍完整复核：第 2、3 节的每一个公式都独立重推了一遍，
第 4 节的每一个表格数字都重算了一遍，所有 claim 都对照代码 `Final_emb_omp.cpp` 核对过。

**总结：数学是对的。** (1)–(115) 里我没有找到一个推导错误。问题集中在排版、交叉引用、
图注、以及几处"论文写的"和"代码跑的"不一致。

---

## A. 必须修（投稿前的 blocker）

### A1. 至少 13 处交叉引用没解析，全是 `??`

| 位置 | 原文 | 应指向 |
|---|---|---|
| p.15 §3.3 Admissibility | "unrelated to the pressure floor $p_\varepsilon$ of Section **??**" | §3.4.4 |
| p.16 式(72)下 | "the largest magnetic eigenvalue (**??**) in direction $d$" | (15) |
| p.21 §4.1 | "listed in Table **??** for $\rho_0=10^{-k}$" | 新表 |
| p.21 §4.1 | "collected in Table **??**" | Table 1 |
| p.21 §4.1 | "Fig. **??** shows the time step permitted by the present scheme" | 缺图 |
| p.21 §4.2 | "Having established the convergence order in Section **??**" | §4.1 |
| p.26 §4.3 | "The initial data are (**??**)–(**??**)" | (101)–(105) |
| p.26 §4.3 | "an exact steady state of (**??**)" | (1) |
| p.26 §4.3 | "Fig. **??** shows the normalised magnetic field" / "Fig. **??** a cut" | Fig. 3 / Fig. 4 |
| p.27 Fig.4 图注 | "the CT method described in Section **??**" | §3.4.2 |
| p.27 §4.4 | "are displayed in Fig. **??**" / "Fig. **??** shows the corresponding distributions" | Fig. 5 |
| p.29 §4.6 | "The splitting parameters at the final time are reported in Fig. **??**" | Fig. 7 |
| p.29 §4.7 | "is reported in Fig. **??** and remains at the level of round-off" | 缺图 |

其中 p.21 的 "Table ??"（$M_c$ 和 $c_a$ 的列表）、"Fig. ??"（时间步比值）和 p.29 的
"Fig. ??"（field loop 的 $\mathrm{div}\,B$）看起来是**图表本身就没放进来**，不只是引用坏了。

### A2. 图 5 和图 7 的图注是错的（都是图 6 的复制）

- **Figure 5**（p.28）画的是 **rotor**（§4.4）：$\rho$ 最大 10.5、域 $[-0.5,0.5]^2$、
  六联图 $\rho,p,|u|^2,m,\varepsilon_M,\varepsilon_B$。图注却写 "for the ideal Orszag–Tang
  vortex at times $t=0.5/2\pi,\,2.0/2\pi,\,3/2\pi,\,5/2\pi$"。
- **Figure 7**（p.31）画的是 **blast wave**（§4.6）：$\log_{10}\rho$、$\log_{10}p$、沿 $x$
  拉长的结构。图注是和图 5 **一模一样**的那句 Orszag–Tang。
- 图 6（p.30）的图注才是真正属于 Orszag–Tang 的。

两张主结果图标错了算例名，这是审稿人一眼就看到的。另外 §4.4 正文说"Fig. ?? 显示
$\rho,p,m$"和"Fig. ?? 显示两个分裂参数"，是两张图；实际只有图 5 一张六联图，正文描述
也要跟着改。

### A3. 式 (46) 印出来是**不连续**的 —— 真正的公式错误

论文印的是

$$s(\mu)=\exp\!\Big(1-\tfrac{1}{1-(\frac{\mu-M_0}{M_1-M_0})^2}\Big)\big[\,\boldsymbol{(1-M_0)^\alpha}-(1-M_1)^\alpha\big]+(1-M_1)^\alpha,\qquad M_0<\mu<M_1$$

在 $\mu=M_0$ 处 bump $=1$，中间支给出 $(1-M_0)^\alpha$；第一支给出 $1-M_0^\alpha$。代入
$M_0=0.15,\ \alpha=14$：

```
第一支  1 - M0^α   = 0.999999999997
中间支  (1-M0)^α   = 0.102770
跳变                = 0.897
```

**$\varepsilon_M$ 会在 $M=0.15$ 处从 $\approx 0$ 跳到 $\approx 0.897$。**

代码是对的，用的是 $(1-M_0^\alpha)$：

```cpp
const double A = (1.0 - std::pow(SW_EPS0, SW_ALPHA))      //  1 - M0^α   ← 论文印成 (1-M0)^α
               -  std::pow(1.0 - SW_EPS1, SW_ALPHA);      //  (1-M1)^α   ← 与论文一致
const double B =  std::pow(1.0 - SW_EPS1, SW_ALPHA);
return bump * A + B;
```

所以是**排版把 $1-M_0^\alpha$ 错排成 $(1-M_0)^\alpha$**。其余三项（$1-\mu^\alpha$、
$(1-M_1)^\alpha$、$(1-\mu)^\alpha$）论文和代码一致。改掉这一项即可。

### A4. Blast wave 初值：论文写锐跳变，实际跑的是 erf 平滑

式 (114) 给的是

$$p_0=\begin{cases}1000,&r\le R\\ 0.1,&r>R\end{cases}$$

而代码当前默认 `g_blast_smooth = 1`，用的是一个单元对角宽的 erf 光滑：

```cpp
double eps = 1 * sqrt(dx*dx + dy*dy);
q[4] = 0.5*(prsR+prsL) + 0.5*(prsR-prsL) * erf((rad - rD)/eps);
```

论文必须说清楚跑的是哪一个。建议在 §4.6 加一句说明用了一个单元宽的 erf 光滑，并给出
`eps` 的定义 —— 这是可复现性问题，不是风格问题。

---

## B. 审稿人大概率会问（应该改）

### B1. Table 1 的 EOC 全部自洽，但"attained third order"的说法和数据对不上

我把 90 个 EOC 全部从表里的 $L^2$ 误差重算了一遍，**全部吻合**（最大偏差 < 0.025，
都在三位有效数字的舍入范围内）。表本身没有算错。

但是：**90 个里有 31 个超过 3.5**，而且低密度两块根本没有向 3 收敛的迹象：

| $\rho_0$ | $N$=256 那一行的 EOC（$\rho,u,p,B_x,A_z$） |
|---|---|
| $10^{-4}$ | 3.94, 2.83, **4.46**, 3.93, 3.82 |
| $10^{-5}$ | **4.59**, 3.39, **4.60**, **4.60**, **4.50** |

最细的网格还在 4.5 左右。正文写 "The formal third order of accuracy is attained for all
variables and across all regimes considered"，结论又重复一遍 "confirm third-order
accuracy" —— 读起来像是测到了 $\approx 3$，但实际测到的是 $3.3\!-\!4.9$。

这其实**对论文有利**（隐式空间算子本来就是四阶，低密度下几乎全隐式，所以表现出 ~4 阶），
但应当照实写，例如"至少三阶；低密度算例呈现约四阶，与隐式算子的四阶空间精度一致"。
现在的写法会被审稿人抓住说"你的 EOC 不收敛到设计阶"。

另外 $\rho_0=10^{-5},\ N=32$ 一行：$u$ 的 $L^2$ 误差是 **5.37e-01**，而速度扰动的幅值
只有 $\tilde v/2\pi=\sqrt{2\pi}/2\pi=0.399$ —— 这个网格什么都没解出来，不应该用来锚定
EOC。建议该算例从 $N=64$ 起算，或者加一句说明。

### B2. $\varepsilon_M$ 在 $M=0.4$ 就已经 $\approx 1$，正文的 bullet 没说清

用修正后的 (46)：

| $M$ | 0.07 | 0.15 | 0.30 | **0.40** | 0.60 | 1.00 |
|---|---|---|---|---|---|---|
| $\varepsilon_M$ | 1.1e-16 | 2.9e-12 | 0.430 | **0.9992** | 0.99999 | 1 |

整个过渡发生在 $[M_0,M_1]=[0.15,0.4]$ 之内，$M\ge 0.4$ 就已经是完全显式。正文的 bullet
写 "Acoustically supersonic ($M\ge 1$): $\varepsilon_M=1$"，容易让人以为显式区从 $M=1$
才开始。建议直接写出过渡窗口。

这一点有实际后果：在激波主导的算例上，解活跃的区域 $\varepsilon$ 基本都是 1，IMEX 相对
全显式几乎没有时间步收益。我在 $1000^2$ blast wave 上实测 $\Delta t/\Delta t_{\rm MHD}$
中位数只有 **1.175**（最小 0.417，最大 4.393）—— 和这个结论一致。论文在 blast wave 一节
没有给 $\Delta t$ 比值，如果给了反而能把"$\varepsilon$ 自适应是局部的、该显式就显式"这件事
讲圆。

### B3. (46) 在 $\mu=M_1$ 只有 $C^0$

bump 支在 $M_1^-$ 的导数为 0，$(1-\mu)^\alpha$ 支在 $M_1$ 的导数是 $-\alpha(1-M_1)^{\alpha-1}
=-0.0183$。所以 $s$（以及 $\varepsilon$）在 $M=0.4$ 有 **0.0183 的斜率跳变**，数值通量对
状态只 Lipschitz、不 $C^1$。

在 $\mu=1$ 处反而没问题（$(1-\mu)^{14}$ 在 $\mu=1$ 是 $C^{13}$），在 $\mu=M_0$ 跳变只有
$3\times10^{-10}$，可忽略。所以唯一的 kink 在 $M_1$。

这和三阶的 claim 有张力，而且 Table 1 的 $\rho_0=10^{-1}$ 算例正好穿过它：
$a=\sqrt{\gamma p_0/\rho_0}=4.08$，$|u|$ 在 $1.0$–$1.8$ 之间 → $M\in[0.25,0.44]$，跨过
$M_1=0.4$。该算例的 EOC 确实是全表最干净的 2.65–3.14，所以经验上不碍事，但正文值得加
一句说明（或者直接用一个在 $M_1$ 也 $C^1$ 的衔接）。

### B4. RP1（Brio–Wu）的 $\gamma$

Table 2 没有 $\gamma$ 列，§4 说 "unless otherwise stated, $\gamma=5/3$"。代码里 `test1`
（= Brio–Wu，$B_x=0.75$、$B_{y,L}=1$、$B_{y,R}=-1$）**没有**单独设 $\gamma$，用的是文件
默认 $5/3$。

**Brio–Wu 的经典设置是 $\gamma=2$。** 熟悉这个算例的审稿人会拿标准图去比，中间态对不上。
要么给 RP1 用 $\gamma=2$，要么在 Table 2 加一列 $\gamma$ 并明说 RP1 用的是 Brio–Wu 的
初值但 $\gamma=5/3$。

### B5. Rotor 的转速只有标准的一半

式 (110)：$(\rho_0,u_0,v_0)=(10,\,-10y,\,10x)$，在 $r=0.1$ 处 $|u|=1$。

Balsara–Spicer [53]（以及 Tóth 2000）的第一个 rotor 是 $|u|=2$ 在 $r_0=0.1$，即
$(-20y,\,20x)$。代码和论文一致（都是 10），所以不是笔误而是有意为之 —— 但引用 [53] 就
应该说明这是修改过的版本，否则看起来像抄错了。

### B6. (51) 叫 "Cargo–Gallice Roe mean" 但少了 $X$ 项

Cargo–Gallice 的 MHD Roe 平均里，声速要带一个磁场跳变的修正：

$$X=\frac{\|\Delta\mathbf B\|^2}{2(\sqrt{\rho_L}+\sqrt{\rho_R})^2},\qquad
\tilde a^2=(\gamma-1)\Big(\tilde H-\tfrac12\|\tilde{\mathbf u}\|^2-\tfrac{\|\tilde{\mathbf B}\|^2}{\tilde\rho}\Big)-(\gamma-2)X$$

(51) 里没有 $X$。对 HLLEM 来说**无害**（守恒由 HLL 主体保证，Roe 均值只进反扩散项和包络），
但"Cargo–Gallice Roe mean"这个称呼就不准确了。建议改成"Cargo–Gallice 型的 $\sqrt\rho$ 加权
平均（略去使其成为精确 Roe 矩阵的 $X$ 项，本格式不需要）"。

### B7. (92) 第二项没做协方差修正

$$\mathcal D^x\Big[(1-\varepsilon^x_M)\Big((1-\varepsilon^x_M)\overline{\mathcal I^x[h^n]\,\mathcal G^x[p]}-\mathcal I^x[h^n]\,\mathcal I^x[p]\,\partial_x\varepsilon^x_M\Big)\Big]$$

第一个乘积用了 (91) 的面平均协方差修正（上横线），第二个乘积 $\mathcal I^x[h^n]\mathcal I^x[p]$
没有。严格说该项在面上只有二阶。因为它带 $\Delta t^2$ 因子，整体应该仍是三阶 —— 但论文
宣称隐式算子是四阶，这里值得补一句说明为什么不影响。

---

## C. 小问题

- **p.7 定理 2.1 证明**："Step 1: symmetrisation. Step 1: symmetrisation." 重复了一遍。
- **p.21 §4.1**："the periodic-free square $\Omega=[-5,5]^2$" —— 应为 "periodic square"。
- **式 (108)** 用 $\tilde\mu=\sqrt{4\pi}$ 参数化压强，而正文说 $\tilde v=\tilde B=1$。
  两者等价（$\frac{1}{8\pi}(\frac{\tilde\mu}{2\pi})^2=\frac{1}{8\pi^2}=\frac12(\frac{\tilde B}{2\pi})^2$，
  即 $\tilde\mu=\sqrt{4\pi}\,\tilde B$），但同一个量用两个符号，建议统一成 $\tilde B$。
- **§4.8**："exceeds four orders of magnitude for $M_x=10^{-4}$" —— 图 10 显示该比值从
  $\approx1.5\times10^4$ 降到 $\approx5\times10^3$，只有开头超过四个量级。图 10 自己的
  图注 "reaches up to four orders of magnitude" 才是准确的说法。
- **参考文献**：
  - [27] `alfv\'en` —— LaTeX 转义泄漏到了正文，应为 "Alfvén"
  - [62] 标题重复："...on general meshes: K. wu, c.-w. shu, Numerische Mathematik..."
  - [49] Ryu–Jones、[50] Falle 都引的是 arXiv preprint，但都已正式发表
    （Ryu & Jones, ApJ 442 (1995) 228）；Elsevier 一般要求引正式版
  - 大量标题被 bib 转成了小写："mhd", "imex", "hll", "euler", "riemann", "mach",
    "alfven", "navier- stokes"（还带一个断字残留的空格）
- **p.1 作者块是空的**：`a,∗,1 , b,2`，单位是 `a . . . .`、`b . . . .`，通讯作者邮箱、
  ORCID 全空；Highlights 页标题下有一个孤零零的逗号。

---

## D. 复核通过的部分（逐条验证）

### 第 2 节

| 式 | 验证方式 | 结果 |
|---|---|---|
| (1) | 四个方程逐个对标准理想 MHD | ✅ |
| (2)(4)(5) | 定义 | ✅ |
| (7) | 1D 通量 8 个分量逐个核对 | ✅ |
| (8)(9)(10) | 八波特征值、$c_a,c_{f,s},c_B$、两个 Mach 数 | ✅ |
| (11)(12) | 逐行相加 $=\mathbf F$；第 5 行用 $E+p+\frac12\|B\|^2=\rho h+k+2m$ 验证 | ✅ |
| (13) | 自己算 $\partial\mathbf F^c/\partial\mathbf U$：$\det=\lambda^4\det(A-\lambda I)$，$A$ 四重根全为 $u$；$A-uI$ 秩 1 → 几何重数 3 → **确实缺一个特征向量** | ✅ |
| (14) | 只有第 2、5 行非零，$2\times2$ 约化块 ${\rm tr}=u$、$\det=-a^2$ → $\lambda^2-u\lambda-a^2=0$ | ✅ |
| (15) | 从原始变量完整推了磁子系统（顺带验证 $p_t\equiv 0$），$5\times5$ 给出 $\lambda^2-u\lambda-c_B^2=0$ 与 $\lambda^2-u\lambda-c_a^2=0$，加 4 个零根 | ✅ |
| (16)(17) | $\tilde{\mathbf F}^c+\tilde{\mathbf F}^p+\tilde{\mathbf F}^B=\mathbf F$；每个分量核对 | ✅ |
| **(18)** | 从 (17) 完整推导原始变量形式，八行逐行核对。**第 5 行**做了完整的能量→压强变换，$\rho_x,v_x,w_x,{B_y}_x,{B_z}_x$ 五组系数全部精确抵消，得 $p_t+\varepsilon_M\rho a^2u_x+\varepsilon_M u\,p_x=0$ | ✅ |
| (19)(20) | $\beta_z\!\cdot\!(v) -\beta_y\!\cdot\!(w)$ 等；$\partial_xu$ 耦合消失因 $\beta_zB_y-\beta_yB_z\equiv0$ | ✅ |
| (21)(22) | 行列式展开；判别式 $(1+\varepsilon_B)^2-4\varepsilon_B=(1-\varepsilon_B)^2$ | ✅ |
| (23)(24)(25)(26) | 两个 $2\times2$ 块、$D_s,D_A$ 消元、磁声四次式 | ✅ |
| (27)(28) | $u=0$ 双二次（用 $c_a^2+b_\perp^2=c_B^2$）、求根公式 | ✅ |
| (29) | $\varepsilon\to0$ 退回 $\{0^{(4)},u^{(4)}\}$；$\varepsilon\to1$ 退回全 MHD 谱 | ✅ |
| **定理 2.1** | $S={\rm diag}(\rho,\rho,\rho,(\rho a^2)^{-1},1,1)$ 使 $SC$ 对称 —— 15 个非对角配对逐个核对；Step 2 的提升公式 $(\rho r_u/(\mu-u),\mathbf r)$；Step 3 的 (32)–(35) 逐行；$\varepsilon_M^2\rho a^2(1-\varepsilon_B)+(1-\varepsilon_M)\varepsilon_B^2B_\perp^2=0$ 的推导（包括 $\varepsilon_B=1$ 的退化情形也成立） | ✅ |
| **(36)** | 自己解了左特征向量方程 $\boldsymbol\ell_r(C-uI)=-\rho\mathbf e_u^{\mathsf T}$，得到 $\ell_1=-\kappa_M$、$\ell_4=-1/(\varepsilon_Ma^2)$、$\ell_{2,3}=-\kappa_MB_{y,z}/B_x$，与论文的 $\kappa_M=\rho(1-\varepsilon_M)u/(\varepsilon_M^2a^2)$ **完全一致** | ✅ |
| (37) | 代回 (18) 八行全部为 0（用 $\alpha(\varepsilon_Bu-\lambda)=\varepsilon_B^2c_a^2$） | ✅ |
| (38) | 从 $\delta u=1$ 逐分量解出，八个分量全部吻合；最后一行自洽性退回 (26) | ✅ |
| (39) | $\tilde{\mathbf F}^p=(1-\varepsilon_M)\mathbf F^p$ ⟹ 特征值乘 $(1-\varepsilon_M)$ | ✅ |
| (40)(41) | $A_5$、$A_2$ 逐行对照我自己推的原始变量方程 | ✅ |
| (42) | $\hat{\mathbf r}_A^\pm=\sigma_v\mathbf e_v+\alpha_A^\pm\mathbf e_B$ 恰好等于 (37) 乘 $-\alpha_A^\pm/B_\perp$；$\boldsymbol\ell_\pm\!\cdot\!\mathbf r_\pm=\varepsilon_B^2c_a^2+(\alpha_A^\pm)^2$ | ✅ |

### 第 3 节

| 式 | 结果 |
|---|---|
| (43)–(45) | ✅ |
| (46) | ❌ 见 A3（印刷错误）；代码版本 ✅ |
| (47)(48)(49)(50) | ✅ |
| (51) | ⚠️ 见 B6（缺 $X$ 项）；$\tilde a^2$ 的形式本身由 $H=h+\frac12\|u\|^2+\|B\|^2/\rho$ 验证 ✅ |
| (52) | ✅ 并且我**数值验证了一个论文没说的性质**：30 万个随机状态下，磁声包络 $[s_L,s_R]$ 始终包含分裂 Alfvén 速度 (22)（最大越界 $0.000$）。这个性质是 HLLEM 的 $\delta_k\in[0,1]$ 成立的前提，论文应该作为一个 remark 写出来 |
| (53) | Fujiwara 界 $2\max(|c_3|,|c_2|^{1/2},|c_1|^{1/3},|c_0/2|^{1/4})$ ✅，与代码一致 |
| (54)(55) | 标准 HLLEM / HLL ✅ |
| (56)(57) | 左右特征向量投影与 Corollary 2.4 的归一化一致 ✅ |
| (58) | $\mathbf R_{\rm ent}=(1,u,v,w,\frac12\|u\|^2,0,0,0)$ 由 $\delta p=\delta\mathbf u=\delta\mathbf B=0$ 验证 ✅ |
| (59) | $\frac12(\lambda\pm|\lambda|)=\max/\min(\lambda,0)$，标准 HLLEM 反扩散系数 ✅ |
| (60)(61) | WENO3-Z，$d_0=1/3,d_1=2/3$，$\tau=|\beta_0-\beta_1|$，指数 $p=2$ ✅（$p=2$ 是临界点保三阶所必需的，论文用对了） |
| (62)(63)(64) | Suresh–Huynh MP 的 $q^{\rm UL},q^{\rm MD},q^{\rm LC}$ 系数 $2$ 和 $4/3$ ✅ |
| **命题 3.1** | Case A 的 $q^W-q_i=\frac12hq'+O(h^2)$ ✅；Case B 的 $\kappa(k\vartheta+\frac12k^2)$、$q(x_{i+1/2})-q_i=\kappa(\frac12\vartheta+\frac18-\frac1{24})$、$d_M=\kappa$、$q^W-q_i=\kappa(\frac12\vartheta+\frac1{12})$、$q^{\rm UL},q^{\rm MD},q^{\rm LC}$、(65) 的 $q^W=q^{\rm MD}+\frac\kappa3=q^{\rm LC}-\kappa$ 全部重推 ✅；下界的 $\max$ 在 $\vartheta=\frac12$ 取到 $\frac13$ ✅，上界在 $\vartheta=-\frac12$ 取到 $\frac16$ ✅ |
| (66) | Zhang–Shu 收缩：代入 $\theta$ 恰得 $\tilde x=x_{\min}$ ✅ |
| (67)–(72) | 线性化形式；(70) 的 $\mathbf F^B_{\rho u}=\frac12\|B\|^2\mathbf e_x-B_x\mathbf B$ 与 (12) 一致 ✅；(72) Balsara flattener 的斜坡区间 $[\kappa_1c^{\min},2\kappa_1c^{\min}]$ ✅ |
| **(73)(74)(75)** | $\partial_tA+\mathbf u\!\cdot\!\nabla A=0$ ✅；(74) 的四阶单元平均乘积（$\Delta x^2/12$ 协方差）✅；(75) 的四阶中心一阶导 ✅，**并且 $D_4^x[\partial_y^{D_4}A]+D_4^y[-\partial_x^{D_4}A]\equiv0$ 因为两个差分算子对易 —— 这就是 $\mathrm{div}\,B$ 到机器精度的原因** ✅ |
| (76)–(84) | 压强方程的完整消元重推了一遍，(82) 的每一项都对 ✅；(81) 的动能/磁能切线化满足 $H(y,y)=$ 原非线性项 ✅；(84) 的守恒性论证成立 ✅ |
| (85)–(88) | SA-DIRK(3,4,3)：$\delta=0.4358665215$ 确为 $x^3-3x^2+\frac32x-\frac16$ 的**中间根**（三根 $\approx0.1588,\,0.4359,\,2.4053$）✅；$b_2+b_3+\delta=1.000000$ ✅；每行行和 $=c_i$ ✅；刚性精确（$A$ 末行 $=b$）✅ |
| (89)(90) | $\bar q=q+\frac1{24}(\delta_x^2+\delta_y^2)q$ ✅；面平均 $\frac1{12}(-\bar q_{i-1}+7\bar q_i+7\bar q_{i+1}-\bar q_{i+2})$ ✅；面法向导数 $\frac{\bar q_{i-1}-15\bar q_i+15\bar q_{i+1}-\bar q_{i+2}}{12\Delta x}$ —— 用原函数 $Q$ 的四阶二阶导模板验证 ✅ |
| (91) | 面平均乘积的横向协方差 ✅ |
| (92) | ⚠️ 见 B7 |
| (93)–(100) | 与 (77)–(84) 的半离散形式逐条对应 ✅ |

### 第 4 节

| 项 | 结果 |
|---|---|
| (101)–(105) | $\delta A_z=\frac{\tilde B}{2\pi}\varphi$ ⟹ $\mathbf B=\frac{\tilde B}{2\pi}r\varphi(-\sin\theta,\cos\theta)$，与 (104) **精确吻合** ✅ |
| **(106)** | 方位角磁场的径向平衡 $\rho v_\theta^2/r=\partial_rp+B_\theta\partial_rB_\theta+B_\theta^2/r$，重排即得 ✅ |
| **(107)** | 从 (106) 解析积分（$\int_r^\infty s\varphi^2=\frac12\varphi^2$，$\int_r^\infty s(2-s^2)\varphi^2=\frac12\varphi^2(1-r^2)$），**精确得到论文的表达式** ✅ |
| 代码 vs 论文的 vortex 初值 | $q[1],q[2],q[4],q[5],q[6]$ 逐项吻合，$\tilde v=\sqrt{2\pi}$、$\tilde B=1$、$p_0=1$、$\mathbf u_0=(1,1,0)$、$[-5,5]^2$、$t_f=1$ 全部一致 ✅ |
| Table 1 | 90 个 EOC 全部从 $L^2$ 重算吻合 ✅（但见 B1） |
| Table 2 RP5/6/7 | 与代码 `test5/6/7` 初值逐项核对 ✅（但 RP1 的 $\gamma$ 见 B4） |
| (108) | 与 (107) 等价（$\tilde\mu=\sqrt{4\pi}\tilde B$）✅ |
| (109)(110) | 与代码一致 ✅（但转速见 B5） |
| (111)(112) | 标准 Orszag–Tang ✅ |
| (113)(114) | $\beta_{\rm amb}=2p/\|B\|^2=0.2/795.77=2.51\times10^{-4}$，"of order $10^{-4}$" ✅ |
| §4.6 的定性描述 | "the lower bound on the pressure is attained in the cells behind the fast shock **in the direction across the field**" —— 我独立测到低压格点集中在 $x\in[-0.096,0.092]$（即 $y$ 轴上，垂直于 $\mathbf B\parallel\hat x$），**完全印证** ✅ |
| (115) | $\|\mathbf B\|=A_0/\sqrt{4\pi}$ 在环内为常数 ✅ |
| §4.7 的 Mach 数 | $M=\sqrt5/\sqrt{\gamma p/\rho}=2.236/408.25=5.48\times10^{-3}$ ✅；$M_A$ 从 $7.93\times10^3$ 到 $7.93$ ✅；$t_f=1$ 确为 $x$、$y$ 各一个周期 ✅ |
| (116)(117) | $\eta$ 在 $y=\pm7/32,\pm9/32$ 四处连续**且** $C^1$（正弦斜坡两端导数为 0）✅；$a=\sqrt{1.4\times1/1.4}=1$ ✅；$M_A=1/(0.1/\sqrt{1.4})=11.83$ 且与 $M_x$ 无关 ✅；$t_f=0.8/M_x=\frac16\cdot4.8/M_x$ ✅ |
| 图 10 | $\Delta t/\Delta t_{\rm MHD}$ 与 $1/M_x$ 成正比（15, 150, 1500, 15000）✅ |

---

## E. 一个额外的数值发现（不是错误，但值得知道）

磁声四次式 (26) 的根**数学上一定是实的**（定理 2.1(i)，我用 $A_5$ 的谱在 20 万个随机
状态上验证了 0 次失败）。但直接对四次式的系数做通用求根，在 **2.4%** 的状态下会返回
带虚部的结果 —— 这些是根接近重合、系数病态的情形。

论文 §3.2 选的是"从 Fujiwara 界出发的单调 Newton"，这**恰好避开了这个问题**（对实根
多项式，从最大根之上出发的 Newton 单调下降）。代码注释里也记了之前用 Cauchy 界时
"blast-like states: 0.6% of the interfaces, envelope up to 20% too wide"。

建议在 §3.2 加半句说明**为什么**用 Newton 而不是闭式/通用求根 —— 这是个真实的鲁棒性
理由，写出来是加分项。

---

## F. 建议的修改优先级

1. **A1** 补全 13 处 `??`，并确认 p.21 的 Table、p.21 和 p.29 的两张 Fig 是不是真的漏了
2. **A2** 改图 5、图 7 的图注（现在都标成了 Orszag–Tang）
3. **A3** 式 (46) 中间支改成 $(1-M_0^\alpha)$
4. **A4** §4.6 写明 blast wave 用的是 erf 光滑还是锐跳变
5. **B1** 改写 §4.1 和结论里关于 EOC 的说法
6. **B4** Table 2 加 $\gamma$ 列（RP1 尤其）
7. 其余 B、C 类按时间允许处理
