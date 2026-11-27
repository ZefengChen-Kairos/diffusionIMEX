# 六个补丁能不能全去掉

清单不是靠记忆，是从**源码重建**的：代码自己在收尾横幅里报告每个保护的触发次数
（`Final_emb_omp_original.cpp:6835`），那个报告就是清单，加上两个未计数的常数。

| # | 保护 | 位置 | 可调常数 |
|---|---|---|---|
| 1 | density point-value guard | `fv_cell`/`fvN`/`rho_pt_guard`（1572, 1605, 3934） | `0.5` |
| 2 | entropy anti-diffusion limiter | HLLEM_x（1368） | `0.5`, `1` |
| 3 | envelope widening | HLLEM_x（1225） | `0.5`, `1` |
| 4 | `g_eps_basis_min = 1e-2` | `build_split_char_x`（2010） | `1e-2` |
| 5 | `g_char_fallback_tol = 1e-6` | 完备性自检（2106） | `1e-6` |
| 6 | 压力地板（嵌套 Newton） | `sec:nestednewton` | **按作者决定保留** |

---

## 决定性测量：IDP 路线还有没有命

完整 128² blast（618 步、2472 次压力解、约 $10^7$ 次单元更新，senv=2, $\kappa_1=0.4$）：

```
U^L inadmissible: rho in 0 cells, p in 909 cells
density point-value guard applied 0 times
entropy anti-diffusion limited at 194832 faces
envelope widened at 0 faces
```

**读法：**

1. **$U^L$ 密度 0 个不可容许** → 低阶更新从不产生负密度，支撑 IDP 的引理经验成立。
   **密度方向的凸限制 $U=U^L+\theta(U^H-U^L)$ 可行。**
2. **$U^L$ 压力 909 个不可容许** → 往一个不可容许的状态混合不可能恢复可容许性。
   **IDP 对压力无效**，压力必须留给保护 6。这不是退让，是测出来的边界。
3. **保护 1 与保护 3 在这个算例上零触发。** 保护 3 的零触发值得注意：
   `g_env_wide` 只在加宽**真的改变** $s_L,s_R$ 时才自增，所以 senv=2（fix E2 的
   局部物质速度地板）可能已经把包络撑够宽，**把保护 3 吞掉了**。正在用 senv=0/1 对照验证。
4. **保护 2 触发 194832 面**（约占全部面-步的 1%）→ **承重**。

---

## 逐条判决

### 保护 2 + 3：真正 ad-hoc 的一对，而且 IDP 原理上救不了

两者用同一个 $J=\lvert\rho_R-\rho_L\rvert/\min(\rho_L,\rho_R)$ 和同两个常数，
源码注释写明是为**同一个失效**存在的：低 Mach 区滞止点（OT 的 X 点，$\varepsilon_M=0$，
$u\to0$）分裂包络是 $O(\lvert u\rvert)\to0$，通量无耗散，而隐式算子是中心型的——
**质量方程在那里没有任何耗散来源**，网格尺度密度结构一直长到点值失效。

**为什么 IDP 救不了**：IDP 往 $U^L$ 混合，而在那里 plain HLL 的包络同样是
$O(\lvert u\rvert)=0$，$s_L=s_R\to0$ 连 HLL 通量本身都退化。$U^L$ 不比 $U^H$ 更耗散。
2+3 做的是**制造两者都没有的耗散**，这超出限制器的能力范围。

**有原则的替代（论文里已有这个机制）**：论文 `eq:mag_stab` 是 $\chi$ 门控的 Laplacian。
它盖不到这个情形的原因很具体——$\chi>0$ 要求
$$h\,\nabla\!\cdot\!\mathbf u<-\kappa_1 c^{\min},\qquad c^{\min}=\sqrt{a^2+c_B^2},$$
而低 Mach 区 $a\gg\lvert u\rvert$，所以 $h\lvert\nabla\!\cdot\!\mathbf u\rvert\sim\lvert u\rvert$
永远打不过 $\kappa_1a$，**$\chi$ 在那里从不触发**。把 $\chi$ 改为门控在
**分裂谱半径 $\lvert\tilde\lambda^c\rvert_{\max}$** 上（显式算子真正解析的那个速度），
则 $\kappa_1\lvert\tilde\lambda^c\rvert\sim0.4\lvert u\rvert$ 与
$h\lvert\nabla\!\cdot\!\mathbf u\rvert$ 同量级，**恰好在那里触发**。

一行改动；用论文自己的机制和自己的 $\kappa_1=0.4$ 替掉两个 $J$-keyed 补丁和它们的
两个常数；而且它比现状**更自洽**：指示器本该门控在它所稳定的那个算子的速度上，
而不是门控在一个被送进隐式求解器的速度上。

### 保护 1：可以换成有名字、无可调常数的限制器

```cpp
if(!(rho_pt >= 0.5*rho_bar)) rho_pt = rho_bar;
```

位置在**反卷积算子内部**（四阶 FV 用 $\bar q\pm d^2\bar q/24$ 在均值与点值间转换），
不是更新里。IDP 碰不到它：IDP 保证**均值**可容许，而这个管的是从可容许均值
插出来的**点值**——即使均值完美，$\bar\rho+d^2\bar\rho/24$ 在 $\rho$ 剧变处仍可为负。
这是插值假象。

**替代**：连续缩放而非 0/1 开关，
$$\rho_{\rm pt}=\bar\rho+\theta\,\frac{d^2\bar\rho}{24},\qquad
\theta=\min\Big(1,\ \frac{\bar\rho}{\lvert d^2\bar\rho/24\rvert}\Big)\ \text{当增量为负},$$
即 Zhang–Shu 缩放限制器。**那个 0.5 直接消失**，只缩放必要的量（设计阶保得更好），
而且可引文献。本次 blast 里它零触发，所以不是承重项。

### 保护 4 + 5：不是可容许性补丁

- **4 `g_eps_basis_min=1e-2`**：**任意选择的重构基底**的条件数地板。源码有实测
  （`Final_emb_omp_original.cpp:168-173`）：$\varepsilon_M=10^{-2}$ 时完备性误差 $10^{-13}$、
  $\lvert L_5\rvert\sim10^{2..3}$；$10^{-6}$ 时误差 $10^{-7}$、$\lvert L_5\rvert\sim10^{7..8}$；
  $\le10^{-10}$ 时误差 $O(1)$、$\lvert L_5\rvert\sim10^{15..19}$（基底数值上已死）。
  $\kappa(R)\sim\varepsilon^{-4}$。**论文里已经写了这条**（§3.1 与 §3.3）。
  原理上可去：$\varepsilon_M\to0$ 时分裂通量就是纯输运，任何基底都可用于重构，
  换成一致良条件的基底（例如非分裂 MHD 基底）即可完全去掉地板。
- **5 `g_char_fallback_tol=1e-6`**：完备性自检，在**真实退化**（$B_n$ 比 $B_\perp$ 低约六个
  数量级）时回退到分量式重构。实测（2084–2088）：$B_n=O(1)$ 时误差 6e-14…3e-12，
  只在那个角落到 2e-6。标准做法。
  - **顺带一个代码问题**：2100 行注释描述了一个 `g_weno3_no_fallback` 开关，说那个
    "过保守"的 1e-6 判据已经退休、只守 NaN/Inf——但**这个变量在源码里根本不存在**，
    2106 行的 `err > g_char_fallback_tol` 仍在无条件生效。注释与代码不一致。

---

## 结论

**六个里：**

- **3 个能真正去掉**：保护 1（无参缩放限制器）、保护 2+3（重新门控 $\chi$，一个机制替两个补丁）。
- **2 个不是补丁**：保护 4、5 守的是真实数值退化与基底条件数，有实测支撑，4 已写进论文。
- **1 个你决定保留**：保护 6。

**"零可调可容许常数"可以达到。"零保护"不行，也不该**——4 和 5 守的是数值退化，不是建模漏洞。

**而 fix G（IDP 凸限制）不是那个"一个方法论解决全部"的答案**：它对密度可行（$U^L$ 零次不可容许），
对压力无效（909 次），对保护 2+3 所解决的问题原理上无能（$U^L$ 与 $U^H$ 同样无耗散）。
它能替掉的是"密度正性"这一条，而那一条在本次 blast 里**本来就没触发过**。
