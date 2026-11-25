#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RP5/RP6/RP7 的 Density / Pressure / B_y 对比图 (exact vs standard IMEX vs dual eps-IMEX).

    python3 plot_rp.py                 # 画全部 3 个算例 x 3 个变量 = 9 张
    python3 plot_rp.py RP6             # 只画一个算例

所有图共用同一套字号和图幅, 符号间隔固定, 便于并排放进论文.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ 配置 ============================
OUTDIR = os.path.dirname(os.path.abspath(__file__)) + "/fig"

# exact: 10 列 x rho u v w p Bx By Bz E
# dual : 12 列 x rho vx vy vz p Bx By Bz E eps_M eps_B
# std  : 3 列  x y value (每个变量一个文件)
# CFL02 = True 用 CFL 0.2 的运行; False 用 CFL 0.8 (论文 §4 的默认值)
CFL02 = True
_S = "_c02_senv1" if CFL02 else "_senv1"
U = "/root/.claude/uploads/fa41aa46-3215-5a47-a5d7-66870c4d191b"
D = os.path.dirname(os.path.abspath(__file__)) + "/data"
CASES = {
    "RP5": dict(exact=U+"/01b3c9db-RP5.dat", dual=D+"/rp5"+_S+".dat",
                std={"rho": U+"/3612a0da-RK5_rho.dat", "p": U+"/c7d25eb0-RK5_p.dat", "By": U+"/ca49d3c3-RK5_By.dat"}),
    "RP6": dict(exact=U+"/89c34670-RP6.dat", dual=D+"/rp6"+_S+".dat",
                std={"rho": U+"/94186312-RK6_rho.dat", "p": U+"/17835d25-RK6_p.dat", "By": U+"/15ad3d0d-RK6_By.dat"}),
    "RP7": dict(exact=U+"/0df4d866-RP7.dat", dual=D+"/rp7"+_S+".dat",
                std={"rho": U+"/b8931b67-RK7_rho.dat", "p": U+"/ed44c5d9-RK7_p.dat", "By": U+"/ff3ed7a3-RK7_By.dat"}),
}
VARS = [("rho", 1, "Density"), ("p", 5, "Pressure"), ("By", 7, r"$\mathbf{B_y}$")]

STRIDE   = 10          # 符号每隔几格画一个 (2000 格 -> 200 个符号)
FIGSIZE  = (7.2, 6.6)      # 统一拉高, 给右上角的 legend 留位置
DPI      = 200
FS_LABEL = 22          # 坐标轴标签
FS_TICK  = 16          # 刻度
FS_LEG   = 16          # 图例
C_EXACT, C_STD, C_DUAL = "black", "tab:orange", "tab:blue"
MS, MEW = 5.0, 1.1     # 符号大小 / 线宽
XLIM = (-0.5, 0.5)
LEGLOC = "upper right"
# 固定 y 范围: (算例, 变量) -> (ymin, ymax); 不在表里的用数据范围.
# 若 legend 会压到数据, 下面会自动把上限抬高 (下限不动).
YLIM = {("RP7","rho"): (0.0700, 0.0900),
        ("RP7","p")  : (0.80, 1.20),
        ("RP7","By") : (-0.30, 0.30)}
# =============================================================

plt.rcParams.update({
    "font.size": FS_TICK, "axes.labelsize": FS_LABEL,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_LEG, "mathtext.default": "bf",
    "axes.linewidth": 1.2, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
})


def load(case):
    c = CASES[case]
    for k in ("exact", "dual"):
        if not os.path.exists(c[k]):
            sys.exit(f"{case}: 缺少 {k} 文件\n  {c[k]}")
    E = np.loadtxt(c["exact"])
    Dn = np.loadtxt(c["dual"])
    if E.shape[0] != Dn.shape[0]:
        sys.exit(f"{case}: exact {E.shape[0]} 行 vs dual {Dn.shape[0]} 行, 网格不一致")
    if np.abs(E[:, 0] - Dn[:, 0]).max() > 1e-9:
        sys.exit(f"{case}: exact 与 dual 的 x 不一致")
    S = {}
    for v, f in c["std"].items():
        if os.path.exists(f):
            a = np.loadtxt(f)
            if np.abs(a[:, 0] - E[:, 0]).max() > 1e-9:
                sys.exit(f"{case}/{v}: standard IMEX 的 x 与 exact 不一致")
            S[v] = a[:, 2]
    return E, Dn, S


def _make_room(fig, ax, leg, pad=0.015, max_grow=3.0):
    """legend 固定在右上角; 若它压到数据, 就抬高 y 上限 (下限不动) 直到不压为止.
    图幅对所有图是同一个, 所以字号的相对大小保持一致."""
    y0, y1 = ax.get_ylim()
    span0 = y1 - y0
    for _ in range(60):
        fig.canvas.draw()
        bb = leg.get_window_extent().transformed(ax.transData.inverted())
        hit = False
        for ln in ax.lines:
            xd, yd = ln.get_xdata(), ln.get_ydata()
            if len(xd) == 0:
                continue
            m = (xd >= bb.x0) & (xd <= bb.x1)
            if m.any() and np.nanmax(yd[m]) > bb.y0 - pad*(ax.get_ylim()[1]-ax.get_ylim()[0]):
                hit = True; break
        if not hit:
            return
        cur0, cur1 = ax.get_ylim()
        if (cur1 - cur0) > max_grow*span0:
            return
        ax.set_ylim(cur0, cur1 + 0.02*span0)


def one(case, E, Dn, S, var, col, ylab):
    x = E[:, 0]
    sl = slice(STRIDE // 2, None, STRIDE)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(x, E[:, col], "-", color=C_EXACT, lw=1.6, label="Exact", zorder=3)
    if var in S:
        ax.plot(x[sl], S[var][sl], "o", mfc="none", mec=C_STD, ms=MS, mew=MEW,
                label="Standard IMEX", zorder=4)
    ax.plot(x[sl], Dn[sl, col], "o", mfc="none", mec=C_DUAL, ms=MS, mew=MEW,
            label=r"Dual $\varepsilon$-IMEX", zorder=5)
    ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax.set_ylabel(ylab, fontweight="bold")
    ax.set_xlim(*XLIM)
    ax.grid(True, color="0.85", lw=0.8)
    ax.set_axisbelow(True)
    if (case, var) in YLIM:
        ax.set_ylim(*YLIM[(case, var)])
    leg = ax.legend(loc=LEGLOC, framealpha=1.0,
                    prop={"weight": "bold", "size": FS_LEG},
                    labelspacing=0.3, handlelength=1.6, handletextpad=0.6,
                    borderpad=0.45, borderaxespad=0.4)
    fig.tight_layout()
    _make_room(fig, ax, leg)
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    out = f"{OUTDIR}/{case}_{var}.png"
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    lo, hi = float(np.min(Dn[:, col])), float(np.max(Dn[:, col]))
    print(f"  {case} {var:3s} -> {os.path.basename(out)}   dual 范围 [{lo:.5g}, {hi:.5g}]"
          + ("" if var in S else "   (没有 standard IMEX 数据, 只画两条)"))


def main():
    cases = sys.argv[1:] or list(CASES)
    for c in cases:
        if c not in CASES:
            sys.exit(f"未知算例 {c}; 可用: {', '.join(CASES)}")
        E, Dn, S = load(c)
        print(f"{c}: {E.shape[0]} 格, standard IMEX 变量 {sorted(S) or '无'}")
        for v, col, ylab in VARS:
            one(c, E, Dn, S, v, col, ylab)


if __name__ == "__main__":
    main()
