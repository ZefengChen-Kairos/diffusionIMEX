# latte — 拉花杯内液体 FVM 求解器（V3 第一阶段）

Robotic Latte Art 仿真的杯内液体部分。规格来源：Notion《V3 · 拉花物理模型》§1–§17
（`v3_proposal_2026-10-28`）。本目录与仓库其余部分（可压缩多相 IMEX 论文代码）完全独立。

* 设计、离散方法、所有自加假设：`docs/stage1_design.md`
* 第一阶段验证结果：`reports/stage1/README.md`（人工总结）与 `reports/stage1/RESULTS.md`（自动生成表格与图）
* 配置（物理、数值、阈值、各测试参数、随机种子）：`configs/stage1/suite.yaml`

```
pip install numpy scipy pyamg pyyaml matplotlib pytest
cd latte
python -m pytest -q tests/            # 快速回归，数秒
python -m latte.stage1.suite          # 第一阶段全部测试（CPU 约 40 分钟）
python -m latte.stage1.suite --only T11,T12
```

| 模块 | 内容 |
|---|---|
| `latte/mesh.py` | 笛卡尔列 + σ 层网格、拓扑、几何（深度场 D 决定的层厚、坡度、面积） |
| `latte/operators.py` | 度量一致、按构造对称半正定的压力/黏性算子 |
| `latte/transport.py` | 体积、组分、质量、动量的同一份面通量输运；联合选择性 FCT |
| `latte/inlet.py` | 交接截面 Γ_J、求积样本、输入契约（零/未知/无效/出杯/事件） |
| `latte/solver.py` | 时间步（V3 §7.4）、自由面、隐式黏性、共同压力投影、回滚与缩步 |
| `latte/remap.py` | 守恒有界垂向重映射（T10） |
| `latte/runio.py` | 运行目录、manifest、metrics（V3 §15 字段） |
| `latte/stage1/` | 第一阶段测试与运行器 |
