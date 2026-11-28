# 干净版源码的清理清单

目标：单一 `.cpp`，决定写进代码，被替换的补丁**删掉**而非关掉，死代码清空。
以下与 Orszag–Tang 实验结果**无关**的部分现在就可以定。

## A. 死代码（声明/注释存在，实际零引用）

| 项 | 证据 | 处理 |
|---|---|---|
| `RusanovK_x/y`, `RusanovK_Update`, `RusanovK_Pressure_x/y`, `RusanovK_Update_Pressure` | 六个前向声明，**无定义、无调用**（未使用的声明不需定义，所以照样链接） | 删 |
| `double epsilon_min = 0.1; double epsilon_max = 1.0;`（3358–3359） | **全文零引用**。连带上方那段注释（"eps==0 remains a valid configuration ONLY with the Rusanov flux path"）描述的也是没在用的变量 | 删 |
| `double Bx_amb = 70.0/sqrt(4π);`（206） | 零引用，$B_0=70$ 时代的残留；现在 $B_0=100$，留着会误导 | 删 |
| `const int g_npicard = 3;` 及 6393 的 Picard 循环 | `emb_on = (magfield==UnStaggerCT && g_ct_A)`，而 `g_ct_A` 在 7047 行**无条件**赋值、`magfield` 默认 `UnStaggerCT` → `emb_on` 恒真 → `for(ipic=1; ipic<1; ...)` 零次 | 删循环与常量 |
| `g_pst_warned` | 只有声明+一处 | 查，大概可删 |
| 22 处 `// [removed: ...]` 注释残桩 | `RusanovC_x/y/Update`, `Magnetic_Subsystem_potential_AX`, `Pressure_NestedNewton` 等 | 删 |

**注意**：删掉 `epsilon_min = 0.1` 会消除一处真实的混淆。实际生效的两个地板是
`eps_M_min = 0`（通量，**无地板**）与 `g_eps_basis_min = 1e-2`（重构基底），
这与论文「while the flux itself is not [bounded away from zero]」一致；
但代码里那个未使用的 `0.1` 会让读者以为通量有 0.1 的地板。

## B. 注释与代码不符

| 位置 | 问题 | 处理 |
|---|---|---|
| 2100 行附近 | 注释描述一个 `g_weno3_no_fallback` 开关，说那个"过保守"的 `g_char_fallback_tol = 1e-6` 判据已退休、只守 NaN/Inf；**该变量全文零声明**，2106 行的 `err > g_char_fallback_tol` 仍无条件生效 | 改注释以匹配代码（或真加那个开关，二选一） |
| 6729 行 | `// match CASTOR: rho/v/p Dirichlet(ambient) + B transmissive` 描述的是 `x_Dirichlet` 分支，而 blast 实际选的是 `x_Transmissive`（零阶外推） | 改注释 |

## C. 运行时开关：只留两个 + 算例数据

作者规格：**`cfl` 与 `senv` 保留为可选项**（`senv` = 物质速度地板取 global 还是 local，
确实按算例选）；**其余一切单一路径，写死，不留分支**。

**保留（两个真正的选项）**

| 键 | 为什么留 |
|---|---|
| `cfl` | 论文各算例用不同 CFL（0.2 / 0.8） |
| `senv` | global $\max\lvert u\rvert$ 还是 local $\max\lvert u_n\rvert$，按算例选 |

**保留（算例数据，不是算法分支；论文的参数扫描需要）**

`tag` `novtk` · `rho0` `vL`（Table 1 需要 $\rho_0=10^0\dots10^{-5}$）·
`b0` `pin` `tend`（blast）· `gamma`（各算例默认不同）

field loop 的 $A_0$ 目前是编译期常量 `LOOP_A0`，论文要四个振幅 → 也应提为算例数据。

**删除（写死为选定值，不留开关）**

| 键 | 写死为 | 依据 |
|---|---|---|
| `kappa1` | `0.4` | 论文 eq.(72) |
| `swm0` `swm1` `swalpha` | `0.15` `0.4` `14` | 论文 eq.(switch)，已核对一致 |
| `epsMmin` `epsMmax` `epsBmin` `epsBmax` | `0` `1` `0` `1` | 论文的 $\varepsilon\in[0,1]$；`epsBmax=0` 只是我的 E1 实验 |
| `plap` | `1`（开） | 论文 eq.(45)/(85) 的激波激活 Laplacian |
| `hllemflat` | 定值 | 取 E16/E17/E18 定下的值 |
| `hllemcap` | 定值 | 同上 |
| `betafloor` | 定值 | 取 $\kappa\beta_{\rm amb}$ 研究定下的值 |
| `epsmax3` | `0` | 实验用 |
| `pswitch` | `0` | 实验用 |
| `bsmooth` | 单一路径 | **需你定**：代码默认 erf 平滑，论文 eq.(107) 写锐跳变 |
| `chisplit` `rhoscale` `rhothr` `jguard` `idp` `hllemoff` | 不留 | 我的诊断开关；决定一旦写进代码就不需要了 |

## D. 待 Orszag–Tang 结果决定

- 若 C（`jguard=0 chisplit=1`）追平 A：**删掉** HLLEM_x 里两个 $J$-keyed 块
  （熵波反扩散限制器 1496 行、包络加宽 1346 行）及其两个常数 0.5/1，
  把 $\chi$ 的门控改为 $\sqrt{\varepsilon_Ma^2+\varepsilon_Bc_A^2}$ 写死。
- 若 C 没追平：两个块**保留**，并在论文里如实描述（它们目前完全没写）。
- `rhoscale`：阈值缩放版写死（无论 OT 结果如何，它都是严格改进）。

---

## 执行记录

### STEP 1 死代码删除：**逐位验证通过** ✓

`make_clean.py` 第一步，128² blast，1212 步，终态 VTK **逐字节相同**。
所以那六个 `RusanovK_*` 声明、`epsilon_min/max` 及其注释块、`Bx_amb`、
22 处 `[removed: ...]` 残桩确实全部不可达。删掉 2406 字符 / 47 行。

### 两个开关不需要任何决定

| 开关 | 默认 | 结论 |
|---|---|---|
| `g_betafloor` | **0（关）** | 本次所有参考运行都没开过 → **不进干净版**，连常数一起省掉 |
| `g_pswitch` | **0（关）** | 同上 |

所以我之前问的"betafloor 用绝对值还是相对形式"是**伪问题**——生产配置里它不存在。
这比两个选项中的任何一个都好，符合"尽量少"。

### 干净版的最终内容

**进（有论文依据）**

| 项 | 值 | 依据 |
|---|---|---|
| 压力 Laplacian | 开 | 论文 eq.(45)/(85)（上篇） |
| $\kappa_1$ | 0.4 | 论文 eq.(72) |
| $M_0,M_1,\alpha$ | 0.15, 0.4, 14 | 论文 eq:switch，已核对 |
| $\varepsilon$ 界 | $[0,1]$ | 论文 |
| 密度点值保护 | 缩放形式，`RHO_PT_FLOOR` | Zhang–Shu 缩放限制器 |
| fix L 统一 flattener | 待阶梯定 | Dumbser–Balsara + Einfeldt + Balsara，三者均已引 |
| 压力地板（嵌套 Newton） | 开 | 作者上篇文章 |
| 基底地板 $10^{-2}$ | 开 | 论文 §3.1/§3.3 |
| 完备性回退 $10^{-6}$ | 开 | 标准做法 |

**保留为运行时选项（作者规格）**：`cfl`、`senv`

**算例数据（论文表格需要）**：`tag` `novtk` `rho0` `vL` `b0` `pin` `tend` `gamma`，
外加应从编译期常量提出来的 field loop 振幅 `LOOP_A0`

**不进**：`idp`（诊断）、`hllemoff`（诊断）、`chisplit`（**已证伪**：OT 全 NaN）、
`jguard`（开关去掉，代码被 fix L 吸收）、`rhoscale`/`rhothr`（写死）、
`epsmax3`（实验）、`pswitch`（默认关）、`betafloor`（默认关）


---

## 干净版：三步全部逐位验证通过

`make_clean.py` 输入**补丁版**（机器齐全），输出 `Final_emb_omp_clean.cpp`。
每步的检验方式相同：两边各跑 128² blast 618 步，终态 VTK 逐字节比。

| 步 | 内容 | 删除 | 验证 |
|---|---|---|---|
| 1 | 死代码 | 2406 字符 | **BIT-IDENTICAL ✓** |
| 2 | fix I（已证伪）+ 探针开关，fix J 写死 | 1191 字符 | 含在下面 |
| 3 | flattener 写死为 L5 级 | 735 字符 | **BIT-IDENTICAL ✓** 对 `flat1=4 hllemcap=0 rhoscale=1` |

共 4332 字符 / 69 行。具名常数：`RHO_PT_FLOOR = 0.5`、`J0_FLAT = 0.5, J1_FLAT = 1.0`。

### 更正一处计数

L5 跑的时候 `hllemflat=2` 是默认值，而 `g_face_chi` 取自 `hchi`，后者在
`hllemflat>=2` 时是 $\max(\chi,s_p)$ —— $s_p$ 是 fix C 的单元级压力跳传感器，
自带一个 0.5。所以 L5 的 $\sigma_\delta=\max(\chi,\ s_p,\ \Theta(J_\rho))$，
**压力传感器仍在**，经由 `g_face_chi` 进来。

我先前说"L5 删掉压力跳"只对 $\sigma$ 层面的 $\Theta(J_p)$ 成立，对 $s_p$ 不成立。
准确计数是 **3 个常数**（$s_p$ 的 0.5，$\Theta$ 的 0.5 和 1），不是 2 个。
L7（`hllemflat=1`，只留 $\chi$）与 L8（`hllemflat=0`，连 $\chi$ 也去掉）正在测，
结果决定最终是 3 个还是 2 个。

### 第四步：尚未做

剩余开关全局量与它们的使用点：

| | 行 | 处理 |
|---|---|---|
| `g_hllemflat` | 1226, 1777, 1779, 1804, 1833 | `hchi` 机器变无条件；1779 的 $s_p$ 块按 L7/L8 结果决定 |
| `g_hllemcap` | 1227 + 键 + banner | 用它的代码已删，去掉全局量 |
| `g_hllemoff` | 1226 + banner | 诊断 |
| `g_idp` | 1223, 1815, 1844, 1884 + banner | fix G 整块，连 `g_hllem_lo`、`Ulo`、`FxL/FyL`、计数器 |
| `g_epsmax3` | 3552, 3657 + 键 + banner | 实验，默认关 |
| `g_pswitch` | 3553, 3641 + 键 + banner | 实验（fix B），默认关 |
| `g_betafloor` | 328, 6538, 6545 + 键 + banner | 默认关，不进 |
| `kappa1`/`swm0,m1,alpha`/`plap` 的键 | 6594 | 值烘死，键删 |

---

## step 4e:包络下界固定为面局部(2026-11-27)

`make_clean.py` 新增的第十步。它把三路 `g_senv` 开关换成无条件的面局部界,并删掉
`senv` 运行时键、banner 字段、`g_senv` 声明,以及 `g_senv_u`(全局 max,只被已删的
全局分支读过)。生成的干净版运行时键降到:

```
b0 cfl gamma novtk pin rho0 tag tend
```

方法层面只剩 `cfl`,其余全是各 test 自己的数据。

### 为什么这一步要单独验证

这是**唯一一个不与上一级逐位相同**的步骤 —— 旧默认是 `senv=2`,新代码无条件走同一段
算术,语义上在 `senv=2` 下是 no-op(生成源的 diff 确认:只有分支外壳、死写入、注释和键
的变化),但 `-O3 -march=native` 下**不**逐位相同:

| | 128² blast |
|---|---|
| 4 步后 | ρ 相对 5.2e-12,p 3.9e-11,B_X 1.9e-14,**`eps_M` 逐位相同** |
| 618 步后 | ρ 相对 4.5e-4 |

4 步就有 1e-12、618 步长到 1e-4,是混沌放大,不是语义差。诊断确认:

```
g++ -O3 -march=native -ffp-contract=off ...
```

下干净版与 patched 跑 `senv=2` **逐位一致**。所以差异来自去掉分支后编译器把
`fmax(fabs(vxL), fmax(fabs(vx_t), fabs(vxR)))` 这条链收缩成了原先没用的 FMA。

**与 step 4c 的 `const` 注记同一类问题。** 不要靠"加回一个分支"去消除它 ——
验证方式改成 `-ffp-contract=off`,不要改代码。

**完整验证已通过**:`-ffp-contract=off` 下 128² blast 跑完全部 618 步,干净版与
patched `senv=2` 的最终 VTK **逐位一致**(2026-11-27 19:01)。

---

## step 4f:删掉转换与 CT 的限制器(2026-11-28)

删掉 `d2lim`/`mm3`/`mm2g`、两个 $\kappa=0.5$、两个全域规约($\verb|grid_absmax|$ 两处与
$S_{gx},S_{gy}$)、整套 $\varepsilon_M$ 混合管线(`g_epsM_cell`/`g_epsM_fx`/`g_epsM_fy`,
3 个全局量 + 4 个赋值点 + 4 个读取点)、以及 `fv_cell_lim_to` 的两趟重试。转换因此
**逐字变回论文 `eq:convdeconv`**。保留 `fv_cell_lim_to` 的正性 guard。共 4644 字符。

同一次提交里把 **fix J 回退**:密度 guard 的动作改回原始的"丢弃增量"
($\rho_{\rm pt}\leftarrow\bar\rho$)。

### 决定依据(限制器开/关)

| 测试 | 开 | 关 |
|---|---|---|
| 低马赫涡 $\rho_0=10^{-2}$,$N=32/64/128$ | — | **逐位相同** |
| 涡 $\rho_0=1$,$N=128$ | — | **逐位相同** |
| 涡 $\rho_0=1$,64→128 阶 | 3.13/2.94/2.94/3.01/3.00 | 3.03/2.95/2.94/3.02/3.02 |
| 涡 $\rho_0=1$,ρ 的 32→64 阶 | 2.48 | **2.67** |
| blast 128² ρ>5% / ρ>1% | 0.101 / 1.316 | 0.101 / 1.373 |
| RP5 ρ>1% 格数 / 中位幅值 | 9 / 1.92% | 8 / 2.18% |

两条结构性事实解释了为什么可以删:

1. **门限的标度比被限制的量少一个 $h$**($\kappa hS$ 对二阶差分、$0.5\max\lvert G\rvert$
   对一阶差分),所以限制器随分辨率自动失效 —— $N=128$ 起逐位相同就是直接证据。
2. **$\varepsilon_M$ 混合在低马赫区已经把它关掉了**,所以论文 Table 1 的低马赫行从来
   没经过它。

### 一个被逐位验证抓到的 bug

把 `for(int pass = 0; pass < 2; pass++){` 换成普通块时,**循环末尾的 `break;` 留了下来**,
改绑到外层的 `i` 循环,导致每行遇到第一个内部单元就跳出。结果:**能编译、能跑、不报错**,
但第 1 步就 NaN、密度 guard 触发 1.44 亿次。删除语句的锚点要用正性 guard 的回退行,
因为 `}`/`break;`/`}` 这个序列在源里有 4 处。

### 验证的三个层级

| 层 | 排除 | 4e | 4f |
|---|---|---|---|
| 1 直接比生产标志 `-O3 -march=native` | 粗错 | ✗ | ✗ |
| 2 两边加 `-ffp-contract=off` | + FMA 收缩 | **够** | ✗ |
| 3 再让参照算**同一个表达式** | + 表达式重结合 | — | **需要** |

4e 只去掉分支、没动表达式形状,第 2 层即可。4f 删掉了 $w\,b+(1-w)\,b$ 的混合,
第 2 层下仍差 ρ 相对 1.7e-4(618 步放大后),**那测的是重结合不是算法**。
参照改成在 `convlim=0` 时把混合权强制为 1(`1.0*b + 0.0*b` 在 IEEE 下精确等于 `b`)后,
**618 步全程逐位一致** ✓。

**生产运行不用 `-ffp-contract=off`**,只有逐位验证用;所有物理测量都来自
`-O3 -march=native -fopenmp` 的二进制。
