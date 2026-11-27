# L9 干净版 ↔ 论文:仍不一致的项

核对对象：`blastwave_debug/patch/Final_emb_omp_clean.cpp`（L9 级，七步逐位验证）
与 `blastwave_debug/paper/paper.tex`。本文件只列**仍然不一致**的；已核对一致的
在 `CODEMATCH.md`。

| # | 代码里有 | 论文 | 处理 |
|---|---|---|---|
| **1** | flattener：$\delta_k\leftarrow(1-\max(\chi,s_p))\,\delta_k$ | **完全没有** | **必须写**——这是 L9 在 HLLEM 之上唯一加的东西 |
| **2** | $s_p$：面邻域的相对压力跳，$\min(1,\max_{\rm nbr}\lvert\Delta p\rvert/\min p\,/\,0.5)$ | 没有 | 随 1 一起写 |
| ~~3~~ | `senv`：HLL 包络的物质速度地板 | ~~`eq:einfeldt` 是裸的 Einfeldt 包络，无地板~~ | **已解决**：论文改了（commit `6b37333`），见下节。代码固定为全局版、删开关 |
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


---

## 第 3 项:已通过改论文解决(2026-11-27, commits `6b37333` → `5d81ee8`)

作者决定:不去迁就论文,而是把论文改成代码做的事。界**保留**(它是必需的),但取在
**界面自己的三个态上**,不是全域。

### 为什么不能"不写"

这不是 Boscheri 那种"大家都知道显式部分用 Rusanov,所以不必写"的情况。论文**正面断言了相反的行为**,
而且当卖点写(原第 1463–1465 行):

> Where $u$ has one sign across the interface, $s_Ls_R=0$ and the flux reduces to the upwind value.

界开着时那个分支永远不触发。按论文原文重新实现的人得到的是另一个格式。

### 为什么界是必需的(不是可选的补丁)

$\varepsilon_M,\varepsilon_B\to0$ 时磁声四次式 $\to(\lambda-u)^4$,分裂 Einfeldt 区间塌成单点 $u$;
$u$ 在界面上同号时 HLL 部分**零耗散**,格式**收敛到错的弱解**:

| RP6 | 无界 | 精确 |
|---|---|---|
| 后激波 $\rho$ | 3.884 | 3.982 |
| 激波位置 | −0.383 | −0.371 |

**从 250 格到 2000 格不变** —— 不是分辨率问题,是错的弱解。所以"不设界"根本不是选项。

### 为什么是面局部而不是全域(这一条我先改错过一次)

两者都能闭合退化,但界是通过

$$\delta_{\rm ent}=1-\frac{|u|}{u^{\max}}$$

进入反扩散的。**取全域**:$|u|\ll\max_\Omega|u|$ 处 $\delta\to1$,把三个线性退化场的耗散
**抽走**——与意图相反。**取面局部**:均匀界面上 $u^{\max}_{i+1/2}=|u|$ 故 $\delta=0$,
只在真有速度跳的界面上为正,反扩散自己关掉。

L9 配置、blast 128²、相对幅值 > 阈值的局部极值格数占比:

| 界 | ρ>1% | p>1% | **ρ>5%** | p>5% | ≤8格条纹% |
|---|---|---|---|---|---|
| 无 | 2.20 | 6.68 | 0.10 | 5.59 | 1.2 |
| 全域 | 6.30 | 9.08 | **1.03** | 7.73 | 0.6 |
| **面局部** | **1.32** | **5.39** | **0.10** | **4.43** | 1.0 |

面局部版的 1D 残差是**局域化的接触后振铃**:RP5 的 21 个振荡格里 19 个落在含接触面的
那一个 0.05 宽 bin 里(中位幅值 1.68%),与全域版的 L1 差 0.33%。

一个曾误用的代理指标:`idp=1` 下"U^L 不可容许压强"计数 1062(无界)/ 909(局部)/ 164(全域),
曾被读成"局部几乎无效"。那个计数器量的是**重构界面态的容许性**,不是弱解正确性;RP6 的
收敛结果才是决定性的。该组数在 `hllemcap=0.25 idp=1` 下测,不是 L9 配置。

### 出处

面局部的 $\max_{\{L,\tilde{\ },R\}}|u|$ 在刚性极限($|\tilde\lambda^c_{\max}|\to|u|$)**就是**
`eq:rusanov` 自己的 $s^{\max}$,所以不需要自己的引用,`Rusanov1962` 已在 bib。不引 Lax,
避免加条目。

### 论文改了哪几处

| 处 | 现在 |
|---|---|
| `eq:einfeldt` | 加 $\mp u^{\max}_{i+1/2}$ 界;新 `eq:umaxface` 定义它取在 $\{L,\tilde{\ },R\}$ 上 |
| 同处 | RP6 反例作为理由;给 1D RP 小节加了 `\label{sec:riemann1d}` 以便引用 |
| AP 段落 | 推导刚性极限:HLL 部分**正好等于** `eq:rusanov`(新 `eq:aplimit`);写明 $\delta_{\rm ent}$ 自关断,以及全域版为何失效 |
| upwind 分支(§1394) | 注明"除非三个态都静止,界使其不可达" |
| $\delta_k\in[0,1]$(§1374) | 界只**加宽**包络,故不影响该包含关系 |

改完比原文**更强**:原文声称零耗散(referee 可攻"那不稳定",而且 RP6 证明确实错),
新文本说 $\eps$-split 把弱双曲性强迫的 Rusanov 基线作为刚性极限恢复出来,并在三个线性
退化场上超过它。

### 代码

`make_clean.py` 新增 step 4e:固定面局部界,删 `g_senv` 开关、`senv` 运行时键、banner 字段,
以及 `g_senv_u`(全局 max,现在无人读)。step 4b 学会了 `apply_patch.py` 为 χ 贡献测量加的
`hllemflat=3` 隔离钩子,把它折叠回 L9 形式。

干净版现在的运行时键:`b0 cfl gamma novtk pin rho0 tag tend` —— 方法层面只剩 `cfl`。

### 仍待定

你明确说过 "RP5 要用全局最大速度做耗散"。固定面局部后 RP5 会得到上述接触面附近那一簇
(约 21/800 点,中位 1.68%)而不是 0%。L9 配置下的 RP5 三路在测。如果那簇点在 1D 图上
不可接受,就得保留 per-test 选择,论文需说明 1D 图用的是全域变体。
