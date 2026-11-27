# L9 干净版 ↔ 论文:仍不一致的项

核对对象：`blastwave_debug/patch/Final_emb_omp_clean.cpp`（L9 级，七步逐位验证）
与 `blastwave_debug/paper/paper.tex`。本文件只列**仍然不一致**的；已核对一致的
在 `CODEMATCH.md`。

| # | 代码里有 | 论文 | 处理 |
|---|---|---|---|
| ~~1~~ | flattener：$\delta_k\leftarrow(1-\sigma)\,\delta_k$ | ~~完全没有~~ | **已写**：新 `eq:flattener`，见下节 |
| ~~2~~ | $s_p$：面邻域的相对压力跳，$\min(1,\max_{\rm nbr}\lvert\Delta p\rvert/\min p\,/\,0.5)$ | ~~没有~~ | **已写**：新 `eq:sp` |
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
| 激波激活压力 Laplacian | 算了 `lap_p` 但**从未加进 Ax**（死代码） | paper **1** 的 eq.(45)/(85) 描述它生效 | **真的加了** ✓ ——但见下面第 12 项：paper **2** 的 `eq:pwave` 里没有这一项 |

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

$\varepsilon_M,\varepsilon_B\to0$ 时 `SplitMagnetosonicCoeffs` 给 $c_3=-2u,\;c_2=u^2,\;c_1=c_0=0$,即
磁声四次式 $\to\lambda^2(\lambda-u)^2$,根 $\to\{0,0,u,u\}$,包络 $\to[\min(0,u),\max(0,u)]$。
$u$ 在界面上同号时**一个端点恰好为 0**,$s_Ls_R=0$,HLL 部分**完全零耗散** —— 不只是在孤立
滞止点,而是在光滑低马赫流的几乎每一个界面上。格式**收敛到错的弱解**:

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

### 出处(已在原文核实)

Paper 1(`MHDIMEX` = Chen, Demattè, Boscheri & Millmore, arXiv:2608.15837)eq.(24)–(25):

$$\mathbf F^c_{i+\frac12,j}=\tfrac12\Big(\mathbf F^c(\mathbf U^R)+\mathbf F^c(\mathbf U^L)-s^{\max}_{i+\frac12,j}\,\mathsf P\,(\mathbf U^R-\mathbf U^L)\Big),\qquad
\mathsf P=\mathrm{diag}(1,1,1,1,1,0,0,0)$$
$$s^{\max}_{i+\frac12,j}=\max\big(|u^L_{i+\frac12,j}|,|u^R_{i+\frac12,j}|\big)$$

**速度是局部的**,取两个界面态的 max;paper 1 正文也写 "the dissipative component of the
explicit Rusanov flux, which **scales with the local flow velocity** $u$"。全域 max 在 paper 1
里只出现在 CFL 时间步 eq.(26),**从不出现在任何通量里**。

Roe 态是凸组合,故 $|\tilde u|\le\max(|u_L|,|u_R|)$,于是 $\max_{\{L,\tilde{\ },R\}}|u|$
**恰好等于** paper 1 的 eq.(25)。所以这个界引入的速度是基线格式本来就在用的那一个:
零新机制、零新引用。

顺带发现 paper 1 的 eq.(24) 带一个投影 $\mathsf P=\mathrm{diag}(1,1,1,1,1,0,0,0)$ ——
Rusanov 耗散**不作用在 $\mathbf B$ 上**。paper 2 的 `eq:rusanov` 没写这个 $\mathsf P$,
对照基线时值得核一下是不是有意省略。

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


---

## 2026-11-27 的通量路径审计:新确认的分歧

对 `paper.tex` 与 `Final_emb_omp_clean.cpp` 做了一次显式通量路径与"全局量渗进局部式"的
系统核对,每条都经过一次对抗性复核(目标是**推翻**它)。下面只列**没被推翻**的,
并标出哪些是我亲自在源码里复验过的。

| # | 代码里有 | 论文 | 严重度 |
|---|---|---|---|
| **8** | 第 1 步的 dt 用**完整 MHD 谱**(`if(timestep < full_step+1){ a_max_x = ch_x; }`,`full_step=1`,:3505),之后每步增长被**压到前一步的 1.5 倍**(`if(dt_prev>0 && dt>1.5*dt_prev) dt=1.5*dt_prev;`,:6632) | `eq:cfl` 前一句写 "**always** subject to a CFL-type condition based \emph{only} on the characteristic speeds of the explicit sub-system" | **blocks-reproduction**(我已复验两处代码;每个 run log 里 `dt/dt_MHD` = 1, 1.50004, 2.25025 … 都看得见) |
| **9** | 显式通量的"面中点→面均值"转换是**按 ε 份额分别做**的:通量先用 HLL 基权拆成 `Pc/Pp/Pb` + 余项 `FxR = f-(Pc+eM*Pp+eB*Pb)`,各自卷积,再用该面的 ε 重组(:1677–1678 / :1696–1697) | §`sec:fv_ops`(:2717–2720)说的是对**数值通量本身**做一次转换 | blocks-reproduction |
| **10** | 横向的面均值↔面中点转换是 **minmod 门 + ε_M 混合**的二阶差分(`fv_transverse_lim`, :1518–1541),两个方向都是;外侧两圈**不转换**。门限 `thr = 0.5*h*max_Ω|q_k|` —— 单元级的判断依赖**全域**量 | `eq:convdeconv`(:2690–2700)是无条件的线性式 | misleading(扩展已有第 4 项:新信息是"两个方向"、"外两圈不转换"、"门限含全域 max") |
| **11** | 喂给 Einfeldt 包络的声速被**地板化**:`a2_t = fmax(a2_t, 1e-12)`(:1172)并把该值传给 Roe 态的 `SplitFastSpeeds`(:1195);左右态用 `fmax(a2L,0.0)` | `eq:einfeldt` 与 $\tilde a^2$ 的定义没有地板;论文的退化条款只覆盖**反扩散项** | misleading |
| **12** | 压力算子里加了一整个**激波激活 Laplacian**,系数是**全域** $\max_\Omega\lvert u\rvert$(`g_plap` 默认 1) | 论文写了**磁**的那一个(§1798–1802,`eq:mag_stab`,连全域 max 一起写了),但 `eq:pwave`(:1952–1958)**没有任何耗散项** | **blocks-reproduction**(我已复验论文两处) |
| **13** | CT 里 $\nabla A$ 的斜率限制器由**全域** $\max\lvert\nabla A\rvert$ 开门(`mm2g`, :4023),并按 $\varepsilon_M$ 混合 | 论文只说 "limited differences" | misleading |
| 14 | $\beta$ 旋转在 $\lvert\tilde B_\perp\rvert\le10^{-12}$ 以下回退到固定向量(:1237–1241);同一个 `eq:beta` 在两处用了**两个不一致**的门限 | 论文未写 | cosmetic |

### 被推翻的(**不是**分歧)

| 曾怀疑 | 复核结论 |
|---|---|
| Alfvén 右特征向量论文给的是原始变量、代码按守恒变量用 | 论文在 :994 写了 $\mathbf M=\partial\mathbf U/\partial\mathbf Q$。只是 $\mathbf M$ 从未逐项写出,属行文问题 |
| 代码把 $\delta_k$ clamp 到 $[0,1]$,而论文说 "without further clipping" | 同一个量在相距 94 行的两处被描述为 "without further clipping"(:1383–1385)与 "restricted to $[0,1]$"(:1479)。是措辞不一致,不是分歧 |
| nested Newton 的陡斜率 $M$ 用了**全域**最小值,论文称 $V$ 是单元局部的 | $M$ 的闭式论文从未写出(代码 :5306–5311)。是"论文没写"的可复现性缺口,不是局部性被破坏 |
| `eq:einfeldt` 写全域界而代码默认局部 | 读的是中间状态。现在两边都是面局部,一致 |

### 我自己写进论文后被这次审计抓到的两个错误(已改,commit `fbb6a60` 与本次)

1. **$(\lambda-u)^4$** —— 错。$\varepsilon=0$ 时 `SplitMagnetosonicCoeffs` 给
   $c_3=-2u,c_2=u^2,c_1=c_0=0$,即 $\lambda^2(\lambda-u)^2$,根 $\{0,0,u,u\}$。
   改正后论证**更强**:包络 $\to[\min(0,u),\max(0,u)]$,$u$ 同号时一个端点恰为 0,
   零耗散发生在光滑低马赫流的**几乎每一个界面**上,不只是孤立滞止点。
2. **"about ten scalar iterations suffice"** —— 在刚性极限是错的。我自己复现了那个
   Newton 循环:

   | 情形 | 迭代 | 得到的精度 |
   |---|---|---|
   | 分离良好($\varepsilon_M=1,a=u=1$) | 13 | 精确 |
   | 刚性极限($\varepsilon=0$) | 33 | 相对 **1.8e-9** |

   双根使 $P$ 在 $\lvert\lambda-\lambda^\pm_f\rvert<\sqrt{\epsilon_{\rm mach}}$ 时成为纯噪声,
   循环靠**符号守卫**退出而非容差。论文已改为实测值,并指出这在该极限下无关紧要 ——
   因为 `eq:einfeldt` 的下界正好在那里把两个速度整个替换掉。


---

## 第 1、2 项:已写进论文(2026-11-27)

新增 `\paragraph{Shock stabilisation of the anti-diffusion.}`,在 AP 标度段之后,含两个新式子
`eq:flattener` 与 `eq:sp`。内容:

1. **机制** —— 与激波**横向**的界面(沿激波面而非跨过它)上法向速度小,熵波速度
   $\tilde\lambda_{\rm ent}=\tilde u$ 离 $[s_L,s_R]$ 两端都远,`eq:deltak` 给出
   $\delta_{\rm ent}\to1$,沿激波面的反扩散完全没被衰减 → 网格对齐激波不稳定(carbuncle)。
   两个 Alfvén 系数同理。
2. **做法** —— $\delta_k\leftarrow(1-\sigma)\delta_k$,$\sigma$ 取两侧单元 $s^p$ 的 max。
3. **三条性质**:
   - $s^p$ 对接触面和旋转间断**完全不响应**(两者压力连续),所以衰减只打激波,
     恰好放过反扩散存在的理由;
   - $\varsigma=1/2$ 是**饱和尺度不是阈值**:线性上升、无死区、无求解器间的跳变;
     饱和处就是裸 HLL,而刚性极限下那正是 `eq:aplimit` 的 Rusanov 基线,所以永不低于
     弱双曲性本来就会强迫的那个通量;
   - **不破坏阶**:光滑流里 $\sigma=O(h)$,乘的是重构后 $O(h^3)$ 的反扩散跳变,
     扰动 $O(h^4)$。低马赫涡($\rho_0=10^{-2}$)实测 3.52/2.89/3.29/3.52/3.64。
4. **$\chi$ 不在里面** —— 实测它可删(见下节),所以论文里 $\sigma=s^p$,只有一个传感器、
   一个常数。

同时把 `eq:einfeldt` 之后那句 "lie in $[0,1]$ ... without further clipping" 改成
"without being clipped there; the shock sensor of eq:flattener is the only further
modification they receive" —— 原句会被读成"没有别的东西碰 $\delta_k$"。

### 一条待核的引用

`\cite{SimonMandal2019}` = Simon & Mandal, *Strategies to cure numerical shock instability
in the HLLEM Riemann solver*。它是**对这个机制的已发表线性扰动分析**,把不稳定归因于
"数值激波区内横向界面上反扩散项的虚假激活",正是我们处理的东西。

**但它的出版信息我没能核实**:`arxiv.org`、`doi.org`、`api.crossref.org` 都被本会话的
网络策略挡住,两个第三方目录给的信息还互相矛盾(一个给 IJNMF 89(12) 533–569 (2019),
另一个把 DOI `10.1016/j.jcp.2018.11.022` 挂在同一个 preprint 上,后者更可能是同作者的
HLLC 姊妹篇)。bib 条目已加入,**并在 .bib 里用显著注释标出投稿前必须核对**
journal / volume / issue / pages / DOI。
