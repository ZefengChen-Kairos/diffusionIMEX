# 跑 1000×1000 blast wave

## 文件

`Final_emb_omp_clean.cpp` —— 单文件,无需 Makefile。默认配置**就是**论文 §4.6 的
blast wave,1000×1000。

## 编译

```bash
g++ -O3 -march=native -fopenmp -std=c++17 -I/usr/include/eigen3 Final_emb_omp_clean.cpp -o mhd
```

唯一的外部依赖是 **Eigen 3**(纯头文件)。AutoDL 上若没有:

```bash
apt-get install -y libeigen3-dev
```

或者下载解包后把 `-I` 指到解包目录。`-march=native` 可按机器换成具体架构;
换掉会让末位结果不同(FMA 收缩),但不影响精度和稳定性。

## 运行

```bash
OMP_NUM_THREADS=<核数> OMP_WAIT_POLICY=passive ./mhd
```

`OMP_WAIT_POLICY=passive` 在核数少于线程数时很重要,否则自旋等待会让速度下降一个数量级。

## 可改的东西

**编译期**(文件顶部 `RUN SETTINGS` 那一段,唯一需要编辑的地方):

```cpp
Test test = testblastwave;      // 换 test
const int nxCells = 1000;       // 网格
const int nyCells = 1000;
double CFL = 0.8;               // 也可用运行时 cfl= 覆盖
```

**运行时**(`key=value`,空格分隔,任意顺序):

| 键 | 默认 | 含义 |
|---|---|---|
| `cfl=` | 0.8 | CFL 数 |
| `tend=` | 0.01 | blast 的终止时间(论文 §4.6 用 0.01) |
| `b0=` | 100 | $\mathbf B_0=(b0,0,0)/\sqrt{4\pi}$ |
| `pin=` | 1000 | 中心压强 |
| `gamma=` | 各 test 自带 | 绝热指数 |
| `rho0=` | 1 | 仅 testVortex 的背景密度 |
| `novtk=1` | 关 | 不写 VTK / state_*.dat(只要日志时用) |
| `tag=` | 空 | 输出文件名的后缀,便于区分多次运行 |

方法本身**没有可调选项** —— flattener、包络下界、$\varepsilon$-split 的参数
($M_0=0.15$、$M_1=0.4$、$\alpha=14$、$\kappa_1=0.4$)全部烘死在源码里,与论文一致。

## 输出

- `plot/<testname>_<N>x<N>_CFL<cfl>_<tag>/plt_*.vtk` —— 各输出时刻的 VTK(ASCII,RECTILINEAR_GRID)
- 同目录 `state_*.dat` —— 全状态二进制
- 同目录 `steplog_*.csv` —— 每步一行:step, t, dt, dt_MHD, maxB, maxB/maxB0, divB_L2, divB_Linf
- stdout 每步一行,末行是汇总计数器(压力求解次数、defect 扫数、未收敛数、flattener 触发面数)

blast 的输出时刻由 `tend` 决定,终止时写一次。

## 跑之前值得确认的两件事

1. **CFL**:论文第 4 节开头写 "Unless otherwise stated ... $\mathrm{CFL}=0.8$",
   blast 小节没有另行声明,所以 0.8 与论文一致。若要和旧结果对比,用 `cfl=` 指定。
2. **初值的平滑**:源码默认把压强跳变做一个单元对角宽的 erf 平滑(与原始代码一致),
   而论文 `eq:blast_ic2` 写的是锐跳变。这是已知且作者决定保留的分歧。
