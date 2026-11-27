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

## C. 运行时开关：逐个决定

当前 26 个。干净版里应分三类。

**保留（算例参数，复现论文表格需要）**
`cfl` `tag` `novtk` `b0` `pin` `tend` `bsmooth` `gamma` `rho0` `vL`

**保留（方法参数，论文里有的）**
`kappa1`（论文 eq.72 = 0.4） `epsMmin/epsMmax/epsBmin/epsBmax`（论文的 $\varepsilon$ 界）
`swm0/swm1/swalpha`（论文 $M_0,M_1,\alpha$） `plap`（论文 eq.45/85 的激波激活 Laplacian）

**删除（我的诊断/实验开关，决定一旦写进代码就不需要了）**
`chisplit` `rhoscale` `rhothr` `jguard` `idp` `hllemoff` `epsmax3` `pswitch`
`senv` `betafloor` `hllemflat` `hllemcap` —— 后四个取定值后写死

## D. 待 Orszag–Tang 结果决定

- 若 C（`jguard=0 chisplit=1`）追平 A：**删掉** HLLEM_x 里两个 $J$-keyed 块
  （熵波反扩散限制器 1496 行、包络加宽 1346 行）及其两个常数 0.5/1，
  把 $\chi$ 的门控改为 $\sqrt{\varepsilon_Ma^2+\varepsilon_Bc_A^2}$ 写死。
- 若 C 没追平：两个块**保留**，并在论文里如实描述（它们目前完全没写）。
- `rhoscale`：阈值缩放版写死（无论 OT 结果如何，它都是严格改进）。
