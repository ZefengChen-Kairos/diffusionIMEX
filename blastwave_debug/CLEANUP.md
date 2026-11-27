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
