# L9 干净版 ↔ 论文:仍不一致的七项

核对对象：`blastwave_debug/patch/Final_emb_omp_clean.cpp`（L9 级，七步逐位验证）
与 `blastwave_debug/paper/paper.tex`。本文件只列**仍然不一致**的；已核对一致的
在 `CODEMATCH.md`。

| # | 代码里有 | 论文 | 处理 |
|---|---|---|---|
| **1** | flattener：$\delta_k\leftarrow(1-\max(\chi,s_p))\,\delta_k$ | **完全没有** | **必须写**——这是 L9 在 HLLEM 之上唯一加的东西 |
| **2** | $s_p$：面邻域的相对压力跳，$\min(1,\max_{\rm nbr}\lvert\Delta p\rvert/\min p\,/\,0.5)$ | 没有 | 随 1 一起写 |
| **3** | `senv`：HLL 包络的物质速度地板（**本次新增的 fix E**，默认 2 = 局部 $\max\lvert u_n\rvert$ over $\{L,{\rm Roe},R\}$） | `eq:einfeldt` 是裸的 Einfeldt 包络，无地板 | **是新增的**：要保留就写进论文，否则删 |
| **4** | 单元均值↔点值的转换有**两族**：未限制的 $d_2/24$（`fv_cell`，3 处）与限制版（`fv_cell_lim_to` 等，3 处）。限制版 = 三候选 minmod + 门 $\kappa hS$（**$\kappa=0.5$**）+ 按 $\varepsilon_M$ 混合 $d_2=\varepsilon_M d_2^{\rm lim}+(1-\varepsilon_M)d_2^{\rm lin}$ | §`sec:fv_ops` 只写了 covariance（`eq:covariance`），这三样都没有 | **比"一句话"大**：三个机制 + 一个常数 |
| 5 | 完备性**数值**自检：固定探针向量的 $\lVert RL\,\text{probe}-\text{probe}\rVert>10^{-6}$ 时**把特征基底整个换成单位矩阵**（重构退化为分量式） | 911 行的 "criterion for completeness"（`eq:criterion`）是双曲性证明里的**解析**判据，不是这个 | 要写或说明 |
| ~~6~~ | 反卷积的密度点值正性缩放，`RHO_PT_FLOOR = 0.5` | 没有 | **不写**（作者决定）：实测在 blast 128²/256²（三种 senv）与 OT 128² 上**从不触发**，是断言而非方法，论文不写自己的 `assert`；留在代码里，计数器会报告，万一在生产运行里触发了它才变成真机制 |
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


---

## 三项机制的确切定义（核过源码）

### $s_p$：0.5 是饱和尺度，不是阈值

```
d  = max over 四个面邻居 of |p_nbr - p_0| / min(p_nbr, p_0)
s_p = min(1, d / 0.5)
```

相对压力跳达到 **50%** 时 $s_p=1$（完全退回 HLL），从 0 线性爬升，**没有死区**。
对比被删掉的 $\Theta(J_\rho)=\min(1,\max(0,(J-0.5)/0.5))$：那个在 50% 以下恒为 0。

后果：光滑流里 $d=O(h)$，所以 $s_p=O(h)$ 而非恒零。**无害**——它乘的那个面跳在三阶
重构后是 $O(h^3)$，扰动 $O(h^4)$；低 Mach 涡（$\rho_0=10^{-2}$）的 EOC 在它开着时
仍是 3.52/2.89/3.29/3.52/3.64，实测确认。

### 限制版二阶差分：三候选 + 门 + $\varepsilon_M$ 混合

```
a = q_i   - 2 q_{i-1} + q_{i-2}      // 中心在 i-1
b = q_{i+1} - 2 q_i   + q_{i-1}      // 中心在 i
c = q_{i+2} - 2 q_{i+1} + q_i        // 中心在 i+1
if (|a|,|b|,|c| 全 < thr) return b;  // 未限制
else                      return mm3(a,b,c);
thr = kappa * h * S,  kappa = 0.5,  S = 该变量的全域 max|q|
d2  = eps_M * d2lim + (1 - eps_M) * d2lin
```

$\varepsilon_M$ 混合有测量理由（源码注释）：minmod 在欠分辨的低 Mach 剪切层上
在相邻面之间开开关关会注入网格尺度动量——KH 64×32、$t=81$ 实测密度漂移
限制版 $2.8\times10^{-3}$ 对线性版 $2.8\times10^{-4}$。

**$\kappa=0.5$ 是一个我先前没数到的可调常数。**

### 完备性自检：运行时的 $RL=I$ 检验

```
probe = (1, -0.7, 0.4, -0.3, 0.9, 0.6, -0.5, 0.2) / ||.||     固定的任意单位向量
S.project(probe, cc);  S.reconstruct(cc, bb);
err = ||bb - probe||
if (!isfinite(err) || err > g_char_fallback_tol)   ->  R := I  （分量式重构）
```

源码自带的实测：$B_n=O(1)$ 时 err $=6\times10^{-14}$（$\varepsilon_M=1$）到
$3\times10^{-12}$（$\varepsilon_M=0.01$）；只有 $B_n\sim10^{-6}$（比横向场低六个
数量级）时升到 $2\times10^{-6}$。所以回退只在那个角落触发。

**注释与代码不符**：注释描述一个 `g_weno3_no_fallback` 开关，说那个"过保守"的
$10^{-6}$ 判据已退休、只守 NaN——**该变量全文零声明**，判据无条件生效。
