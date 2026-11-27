# L9 干净版 ↔ 论文:仍不一致的七项

核对对象：`blastwave_debug/patch/Final_emb_omp_clean.cpp`（L9 级，七步逐位验证）
与 `blastwave_debug/paper/paper.tex`。本文件只列**仍然不一致**的；已核对一致的
在 `CODEMATCH.md`。

| # | 代码里有 | 论文 | 处理 |
|---|---|---|---|
| **1** | flattener：$\delta_k\leftarrow(1-\max(\chi,s_p))\,\delta_k$ | **完全没有** | **必须写**——这是 L9 在 HLLEM 之上唯一加的东西 |
| **2** | $s_p$：面邻域的相对压力跳，$\min(1,\max_{\rm nbr}\lvert\Delta p\rvert/\min p\,/\,0.5)$ | 没有 | 随 1 一起写 |
| **3** | `senv`：HLL 包络的物质速度地板（**本次新增的 fix E**，默认 2 = 局部 $\max\lvert u_n\rvert$ over $\{L,{\rm Roe},R\}$） | `eq:einfeldt` 是裸的 Einfeldt 包络，无地板 | **是新增的**：要保留就写进论文，否则删 |
| 4 | FV 修正的二阶差分用 3 候选 minmod 限制 | §`sec:fv_ops` 写了 covariance（`eq:covariance`），没写这个限制 | 一句话补上 |
| 5 | 完备性**数值**自检 + 分量式回退，容差 `g_char_fallback_tol = 1e-6` | 911 行的 "criterion for completeness"（`eq:criterion`）是双曲性证明里的**解析**判据，不是这个 | 要写或说明 |
| 6 | 反卷积的密度点值正性缩放，`RHO_PT_FLOOR = 0.5` | 没有 | 要写（实测在 blast 与 OT 上**从不触发**，可一句话带过） |
| 7 | blast 初值的 erf 平滑（一个单元对角宽） | `eq:blast_ic2` 写锐跳变 | **作者已决定保持此分歧**，论文不写 |

不是分歧：`g_eps_basis_min = 1e-2` 的重构基底地板，论文 §3.1 与 §3.3 **已写**。

## 第 1 项的一处张力

`eq:einfeldt` 之后有一句：

> "The anti-diffusion coefficients $\delta_k$ of~\eqref{eq:deltak} therefore lie in
> $[0,1]$ for all three linearly degenerate fields **without further clipping**."

原意是"不需要额外 clamp 到 $[0,1]$"，但 flattener 恰恰在进一步乘 $(1-\sigma)$。
审稿人会把这句读成"没有别的东西碰 $\delta_k$"。写 flattener 时这句必须改。

## 原始码与论文不符、而干净版已修好的两项

说明论文原本描述是对的，代码没做到：

| | 原始码 | 论文 | 干净版 |
|---|---|---|---|
| $\kappa_1$ | `0.05` | `eq.(72)` 是 `0.4` | **0.4** ✓ |
| 激波激活压力 Laplacian | 算了 `lap_p` 但**从未加进 Ax**（死代码） | `eq.(45)/(85)` 描述它生效 | **真的加了** ✓ |

## 第 1、2 项的引用（已查证）

| 块 | 出处 | 论文里有吗 |
|---|---|---|
| $\delta_k$ 机器 | Dumbser & Balsara, *A new efficient formulation of the HLLEM Riemann solver*, JCP **304** (2016) 275–319, DOI 10.1016/j.jcp.2015.10.014 | ✓ 已引（`\cite{HLLEM}`） |
| $\chi$ | Balsara (2012) | ✓ 已引（`eq:chi`） |
| "强激波附近把低耗散通量退回 HLL，开关用压力" | **Quirk** 最早（*A contribution to the great Riemann solver debate*，细节待核）；现代版 **Xie, Li, Li, Tian & Pan**, JCP **350** (2017) 607–637, DOI 10.1016/j.jcp.2017.08.063 | ✗ **需新增** |
| $s_p$ 的逐字形式 | 未找到逐字出处；"压力开关"正是 Quirk 的做法，可作为该构造的一个实例 | — |

**没有逐字出处的只有 $s_p$ 的具体形式**，而它可以写成被引构造的实例，不必声称新颖。
