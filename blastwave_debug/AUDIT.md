# 叙事审计（4 reader × 对抗性验证，49 agent / 4.8M token）

44 条候选 → 5 条通过对抗性验证。其中 **3 条在审计期间已被我后来的改写消掉**
（审计原话：*"fixed better than the proposed fixes"*）。另外它抓到了**我自己新引入的
两个错误**，都往好看的方向偏——下面标 ⚠ 的两条。

我独立复核了审计的两个承重数字，**都准确**：

| 审计的数 | 我的复核 | |
|---|---|---|
| 表里 90 个 EOC 有 31 个 < 3.00 | 90 个，31 个 < 3.00（34%），min 2.610、median 3.085、max 4.880、无 < 2.5 | ✓ |
| $\max[1+\varepsilon_B/M_A]=3.627$ @ $M_A=0.3755$ | 用论文 `eq:alfven_eig` 精确式：$\max=3.6213$ @ $M_A=0.3768$（$\varepsilon_B=0.9897$）；并证明 $1+\varepsilon_B/M_A$ 是严格上界（差 $4\varepsilon(1-\varepsilon)/M\ge0$） | ✓ |

另独立验了 KH 的 $M_A$：代码 4103 行 `Az_val = 0.1*KH_MX*y` → $B_x=0.1M_x$，
$\rho=\gamma=1.4$ → $c_B=0.0845M_x$，$\lvert u\rvert=M_x$ → $M_A=11.83$，因子 1.0845 ✓

## 已改（14 处）

| # | 位置 | 问题 |
|---|---|---|
| 1 | §3.2 `sec:hllem` 开头 | **被判为"论文里最糟的一处"**。「Whenever $\epsM>0$ ... This makes the HLLEM scheme ... applicable」与 §3.1 仅隔 80 行的「**The Riemann solver does not require $\epsM>0$**」正面矛盾；也与 `eq:hllem`（只对 ent 与 $A^\pm$ 求和）、`eq:walf`（只含 $\varepsilon_B$）、`eq:went`（不含 $\varepsilon$）矛盾；也与运行点矛盾（$\varepsilon_M$ 在 $M\lesssim0.07$ 就低于机器精度，而那里照用 HLLEM）。改成「HLLEM 只需线性退化场及其投影，不需完整基底」，并把 $\varepsilon_M\ge10^{-2}$ 的界明确限定在重构上 |
| 2 | 同段末 | "supplied by the split basis" → "by the block decomposition"（前者在完整基底与块分解之间含混） |
| 3 | 导言 leg-two | 「regular both as $B_\perp\to0$ **and as $\epsB B_x\to0$**」超出 `cor:alfven`（只证了 $B_\perp\to0$）；正文说那里**范数化本身消失**并在 $10^{-8}$ 以下换用解耦分支。改成"闭形式在所有参数值下可得 + 解耦极限是参数化的退化而非系统的退化" |
| 4 ⚠ | 导言 cost 段 | **我的错**：「a factor bounded by two, since a non-stiff channel has a wave speed below the material speed」。错在把"非刚性"当成 $M_A\ge1$，而 `eq:switch` 在 $M_A\approx0.4$ 就已返回 $\varepsilon_B\approx0.999$。真实上界 **3.63**，在 $M_A\approx0.38$ 取到（那里 $c_B$ 仍是 $2.7\lvert u\rvert$）。同时把无出处的 1.19 / "four digits" 换成有出处的解析因子 1.13（field loop）与 1.08（KH） |
| 5 ⚠ | 结论 dt 句 | **我的错**：$2.3\times10^{2}$ 与 1.08 在第 4 节里查不到（审稿人可以直接 grep）。改成引用论文**已有**的 KH 结果（`fig:kh_dtratio`，四个数量级） |
| 6 | 摘要 | 两条腿写反：强双曲性被断言却**不给任何后果**，而 characteristic reconstruction 被归给块分解（导言与结论都归给定理的完整性）。并补上 $\varepsilon_M\ge10^{-2}$ 地板 |
| 7 | 摘要 + 结论 + §4.1 | "at least third-order" 对着 31/90 低于 3.00 的表。改成"design third order，低密度处从下方逼近、最刚性处超出（粗网格尚未进入渐近区）" |
| 8 | Highlights | 从未被改写触及，现在是**全文唯一无条件断言完整性的地方**；"no iteration" 未限定；时间步仍被当作卖点而 resolution 完全缺席。四条全部重写 |
| 9 | 导言动机 | 同一混淆：「precludes ... Roe or HLLEM」。改成"依赖显式 Jacobian 特征向量的解算器"，并指向 §3.2 |
| 10 | §2.3 末 | 同上：「neither a complete approximate Riemann solver nor a characteristic reconstruction can be built」 |
| 11 | Standard IMEX baseline | 同上：改成"整个磁通量在隐式侧时 Alfvén 场在显式谱里不留痕迹，所以没有东西供中间波修正作用" |
| 12 | §sec:CT | MF4 建议的 "reported after eq:CTmap" 实际写反了——测量在 1796 行、`eq:CTmap` 在 1811 行。改为 `Section~\ref{sec:CT}` |
| 13 | KH 节 | 补 `\label{sec:KH}`（MF4 的交叉引用需要它） |
| 14 | Field loop 节 | 补 `\label{sec:fieldloop}`（见 CODEMATCH §6） |

## 未改，待你定（审计的 should_fix / still_unfixed 里仍活的）

1. **"accurate and efficient across all Mach and Alfvén number regimes"**（结论）——
   全文没有任何计时、迭代次数、defect correction 容差或成本对比，而每步成本是
   磁预测器一次 GMRES 加压力的嵌套 Newton + defect correction。要么给数据，要么去掉 "efficient"。
2. **能量守恒是双重条件的**——正文里取决于 defect correction 收敛**且**压力地板未激活，
   这两个条件都没进摘要/导言/结论。
3. **嵌套 Newton 在摘要/导言/结论里完全没出现**，尽管它在每个 defect 步里都在跑。
4. **`thm:hyp`(i)（所有参数值下八个实特征值）没有被说明它许可什么**——
   它支撑 `eq:cfl` 和 Newton 包络，但现在读起来像是也依赖 $\varepsilon_M>0$。
   审计认为这是比条件性的完整性**更强更稳**的陈述，被低估了。
5. **`prop:mp3`（MP 限制器下光滑极值处保持设计阶）是一个被证明的结果，是贡献而非规格**——
   "third order is a specification" 的降格把它一起降了。
6. **双向权衡仍未被测量**：第 4 节没有本方案的单自适应参数运行；唯一的基线是
   固定 $\varepsilon\equiv0$ + Rusanov，而那条基线**现在的代码跑不出来**（见 CODEMATCH §15）。
7. **低 $\beta$ 中等 Mach 侧只有 blast wave 考验**，而它的 $\varepsilon$ 统计来自 256²（已标注待刷新）。
8. 摘要说旋转间断被解析得 "sharply"，结论对 resolution 什么都不说；
   鉴于 2D 面内场那条 $1.2\times10^{-6}$ 的测量，"sharply" 需要限定到 1D。
