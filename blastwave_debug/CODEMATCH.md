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

## 9. 待核

- ○ §Eigenstructure of the explicit sub-system（已在 `REVIEW_IMEX_MHD-3.md` 逐式推过）
- ○ §Time discretisation（`eq:88` 表已核；其余未核）
- ○ §Finite volume operators / WENO3-Z 重构
- ○ 各算例初值与论文表格
