#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_blast.py —— blast wave 六格图, 改好 FILE 一行, 直接运行.

    python3 plot_blast.py                      # 用下面的 FILE
    python3 plot_blast.py 某个.vtk             # 或命令行给路径

依赖只有 numpy + matplotlib (不需要 pyvista, 内置 ASCII legacy VTK 解析器).

关于压力地板 (BETA_FLOOR):
  本算例 B0 = 70/sqrt(4*pi), 磁压 m = |B|^2/2 = 194.96, 环境热压 0.1,
  即 beta = 2p/|B|^2 = 5.129e-4, 热能只占总能的 0.128%. 于是
  p = (gamma-1)(E - KE - |B|^2/2) 是灾难性相消: 快激波脚处 CT 场与守恒通量
  之间不到一个网格的相位失配 (|B| 偏 0.056%) 就能把环境压力从 0.1 打到 0.006,
  对数图的下界被这一圈薄环拽到 -2.2.
  BETA_FLOOR 在取 log10 之前做 p <- max(p, BETA_FLOOR * |B|^2/2), 和求解器的
  betafloor= 同一个判据, 于是已经算完的数据不必重跑也能出图.
  脚本会打印被抬起的格子数 —— 这是作图处理不是物理, 图注里要写明这个数字.
  设 BETA_FLOOR = None 就关掉, 看原始数据.
"""

import os
import sys

import numpy as np
import matplotlib
if not os.environ.get("DISPLAY") and sys.platform != "darwin":
    matplotlib.use("Agg")          # 无显示环境: 只存文件
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import ImageGrid

# ========================== 配 置 区 ==========================

FILE = "/Users/chenzefeng/Downloads/plt_blast_wave_1000x1000_CFL0.4_fix_2531.vtk"
OUT  = None                 # None = 存到 FILE 同目录的 panel_blast.png

# 真实计算域 (文件里的 X/Y_COORDINATES 存的是偏移过的值, 一律按这里重建)
DOMAIN = ((-0.5, 0.5), (-0.5, 0.5))

# 取 log10 之前的压力地板, 见文件开头说明
BETA_FLOOR = 5.0e-4         # p >= BETA_FLOOR * |B|^2/2 ; None = 关掉
P_FLOOR    = None           # 另加一个绝对下界, 两者取大; None = 不加

# 六格布局
PANEL_FIELDS = [["log_rho", "log_p"],
                ["velMag2", "pmag"],
                ["eps_M",   "eps_B"]]

CBAR_LABELS = {
    "rho": r"$\rho$", "p": r"$p$", "pmag": r"$\mathbf{m}$",
    "log_rho": r"$\log_{10}\,\rho$", "log_p": r"$\log_{10}\,p$",
    "velMag": r"$|\mathbf{u}|$", "velMag2": r"$|\mathbf{u}|^2$",
    "Bmag": r"$\|\mathbf{B}\|$", "eps_M": r"$\varepsilon_M$", "eps_B": r"$\varepsilon_B$",
}

# 固定色标. 不给的 field 用数据自身的范围.
# 注意: log_p 的下界不要取正好 -1.0 —— 远场约 30% 的格子停在 0.09999 这种
# "比环境值低一丝"的位置, -1.0 会把它们整片判为超范围 (实测 128^2 刷掉 31%).
# 加了地板之后 log_p 的自然范围就是 [-1.01, 2.30], 本来就不需要固定.
FIELD_RANGE  = {"eps_M": (0.0, 1.0), "eps_B": (0.0, 1.0)}
FIELD_LEVELS = {"eps_M": 0, "eps_B": 0}      # 0 = 不画黑色等值线

CLIP_UNDER, CLIP_OVER, CLIP_BAD = "magenta", "lime", "white"   # 超范围/无定义的格子

CMAP, LEVELS = "jet", 30
LINE_COLOR, LINE_WIDTH = "k", 0.35
SMOOTH = True               # True = gouraud 平滑; 查单格缺陷时设 False
TICK_STEP = 0.2
PANEL_W, PANEL_H = 6.0, 5.2
PAD_W, PAD_H = 1.25, 0.85
DPI, SHOW, SAVE = 150, True, True

GAMMA_FALLBACK = 1.4        # 反推不出来时用它 (blast/rotor 1.4, OT 5/3)
SW_M0, SW_M1, SW_ALPHA = 0.15, 0.40, 14.0      # 和 C++ switchS 一致
EPS_MIN, EPS_MAX = 0.0, 1.0

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "axes.titlesize": 22,
                     "xtick.labelsize": 14, "ytick.labelsize": 14})
plt.rcParams["mathtext.default"] = "bf"

# ======================= VTK 读取 (内置) =======================

_HEADER = {"DATASET", "DIMENSIONS", "ORIGIN", "SPACING", "X_COORDINATES",
           "Y_COORDINATES", "Z_COORDINATES", "POINT_DATA", "CELL_DATA",
           "FIELD", "SCALARS", "VECTORS", "LOOKUP_TABLE", "METADATA", "POINTS"}


def _is_header(s):
    tok = s.split()[0]
    if tok.upper() in _HEADER:
        return True
    try:
        float(tok); return False
    except ValueError:
        return True


def read_vtk(path):
    """ASCII legacy RECTILINEAR_GRID + CELL_DATA/FIELD. 容忍文件被截断."""
    if not os.path.exists(path):
        sys.exit(f"文件不存在:\n  {path}")
    txt = open(path, errors="replace").read().split("\n")
    i, n = 0, len(txt)
    dims, fields, t, truncated = None, {}, float("nan"), False

    def nextline():
        nonlocal i
        while i < n:
            s = txt[i].strip(); i += 1
            if s:
                return s
        return None

    def floats(k):
        """读 k 个数; 返回 (数组, 是否读满)."""
        nonlocal i
        out, got = np.empty(k), 0
        while got < k and i < n:
            s = txt[i].strip()
            if not s:
                i += 1; continue
            if _is_header(s):
                return out[:got], False
            a = np.fromstring(s, dtype=float, sep=" ")
            i += 1
            if a.size == 0:
                return out[:got], False
            m = min(a.size, k - got)
            out[got:got + m] = a[:m]; got += m
        return (out, True) if got == k else (out[:got], False)

    while True:
        s = nextline()
        if s is None:
            break
        tok = s.split(); key = tok[0].upper()
        if key == "DIMENSIONS":
            dims = tuple(int(v) for v in tok[1:4])
        elif key == "FIELD":
            for _ in range(int(tok[2]) if len(tok) > 2 else 0):
                h = nextline()
                if h is None:
                    break
                hp = h.split()
                if len(hp) < 4:
                    i -= 1; break
                name, ncomp, ntup = hp[0], int(hp[1]), int(hp[2])
                vals, ok = floats(ncomp * ntup)
                if not ok:
                    truncated = True
                    print(f"[警告] 数组 '{name}' 不完整 ({vals.size}/{ncomp*ntup}), 丢弃 —— 文件被截断")
                    continue
                if name in ("TIME", "CYCLE"):     # FieldData 头部, 不是网格数组
                    if name == "TIME":
                        t = float(vals[0])
                elif ncomp == 1:
                    fields[name] = vals
                else:
                    a = vals.reshape(ntup, ncomp)
                    for c, sfx in zip(range(ncomp), "_x _y _z".split()):
                        fields[name + sfx] = a[:, c]
        elif key in ("X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"):
            # 坐标一律按 DOMAIN 重建 (writer 存的是偏移过的 cell center, 而且
            # 声明的个数常和实际写出的差 1), 所以这里只是跳过, 不判断完整性.
            floats(int(tok[1]))
        elif key == "SCALARS":
            name = tok[1]
            nx2 = nextline()
            if nx2 and not nx2.upper().startswith("LOOKUP_TABLE"):
                i -= 1
            ncell = (dims[0]-1)*(dims[1]-1) if dims else 0
            vals, ok = floats(ncell)
            if ok:
                fields[name] = vals
            else:
                truncated = True

    if dims is None:
        sys.exit(f"没读到 DIMENSIONS, 这不像一个 legacy VTK 文件:\n  {path}")
    nx, ny = dims[0] - 1, dims[1] - 1
    ncell = nx * ny
    bad = [k for k, v in fields.items() if v.size != ncell]
    for k in bad:
        print(f"[警告] 数组 '{k}' 长度 {fields[k].size} != {ncell}, 丢弃"); fields.pop(k)
    if not fields:
        sys.exit(f"没读到任何和网格对得上的数组 (文件 {os.path.getsize(path)/1e6:.1f} MB, "
                 f"可能被截断):\n  {path}")
    print(f"读入 {os.path.basename(path)}: {nx}x{ny}, t = {t:.6g}, "
          f"{len(fields)} 个数组{'  [文件被截断, 只恢复了完整的部分]' if truncated else ''}")
    return fields, nx, ny, t


# ======================= 派生场 + 地板 =======================

def _B(fields):
    def pick(c):
        for k in c:
            if k in fields:
                return k
        return None
    return (pick(["B_X", "B_x", "Bx"]), pick(["B_Y", "B_y", "By"]), pick(["B_Z", "B_z", "Bz"]))


def detect_gamma(fields):
    """从守恒量反推: p == (gamma-1)(U - KE - |B|^2/2)."""
    if not all(k in fields for k in ("p", "U", "rho", "mom_X", "mom_Y", "mom_Z")):
        return GAMMA_FALLBACK
    bx, by, bz = _B(fields)
    if bx is None or by is None:
        return GAMMA_FALLBACK
    with np.errstate(divide="ignore", invalid="ignore"):
        ke = 0.5*(fields["mom_X"]**2 + fields["mom_Y"]**2 + fields["mom_Z"]**2) \
             / np.maximum(fields["rho"], 1e-300)
    mb = 0.5*(fields[bx]**2 + fields[by]**2 + (fields[bz]**2 if bz else 0.0))
    e = fields["U"] - ke - mb
    scale = max(1.0, float(np.nanmax(np.abs(fields["p"]))))
    best, err = GAMMA_FALLBACK, np.inf
    for g in (1.4, 5.0/3.0, 2.0, 1.2, 1.1):
        d = float(np.nanmax(np.abs(fields["p"] - (g - 1.0)*e)))
        if d < err:
            best, err = g, d
    if err < 1e-8*scale:
        print(f"[gamma] 从守恒量反推 gamma = {best:.6g} (残差 {err:.3g})")
        return best
    print(f"[gamma] 反推失败 (最小残差 {err:.3g}), 用 {GAMMA_FALLBACK:.6g}")
    return GAMMA_FALLBACK


def _switch_s(M):
    """C++ switchS 的向量化版本, 分支顺序一致."""
    M = np.asarray(M, float); a, m0, m1 = SW_ALPHA, SW_M0, SW_M1
    with np.errstate(all="ignore"):
        u = (M - m0)/(m1 - m0)
        bump = np.exp(1.0 - 1.0/(1.0 - u*u))
        A = (1.0 - m0**a) - (1.0 - m1)**a
        B = (1.0 - m1)**a
        return np.select([M <= 0.0, M >= 1.0, M <= m0, M >= m1],
                         [1.0, 0.0, 1.0 - M**a, (1.0 - M)**a], default=bump*A + B)


def apply_floor(p, fields):
    """取 log 之前给压力下界. 返回新数组, 不动 fields['p']."""
    floor = None
    if BETA_FLOOR:
        bx, by, bz = _B(fields)
        if bx is None or by is None:
            print("[地板] 找不到 B 分量, BETA_FLOOR 跳过")
        else:
            mb = 0.5*(fields[bx]**2 + fields[by]**2 + (fields[bz]**2 if bz else 0.0))
            floor = BETA_FLOOR*np.asarray(mb, float)
    if P_FLOOR:
        pf = np.full_like(np.asarray(p, float), float(P_FLOOR))
        floor = pf if floor is None else np.maximum(floor, pf)
    if floor is None:
        print(f"[地板] 关闭, 用原始数据 (p_min = {np.nanmin(p):.6g})")
        return np.asarray(p, float)
    out = np.maximum(p, floor)
    lift = p < floor
    if lift.any():
        print(f"[地板] {int(lift.sum())} 个格子 ({100.0*lift.mean():.3f}%) 被抬到下界: "
              f"p_min {np.nanmin(p):.6g} -> {np.nanmin(out):.6g}"
              + (f"   (beta >= {BETA_FLOOR:g}, 即 p >= {np.nanmin(floor):.4g}..{np.nanmax(floor):.4g})"
                 if BETA_FLOOR else ""))
        print(f"       ^ 作图地板, 不是物理. 图注里写这个数: {100.0*lift.mean():.2f}% 的网格")
    else:
        print(f"[地板] 没有格子低于下界 (p_min = {np.nanmin(p):.6g}), 图是原始数据")
    return out


def derive(fields):
    gamma = detect_gamma(fields)
    bx, by, bz = _B(fields)
    p_eff = None
    if bx and by:
        Bz = fields[bz] if bz else 0.0
        B2 = fields[bx]**2 + fields[by]**2 + Bz**2
        fields["pmag"] = 0.5*B2
        fields["Bmag"] = np.sqrt(B2)
    if "p" in fields:
        p_eff = apply_floor(np.asarray(fields["p"], float), fields)
    if all(k in fields for k in ("velX", "velY")):
        vz = fields.get("velZ", 0.0)
        v2 = fields["velX"]**2 + fields["velY"]**2 + vz**2
        fields["velMag2"] = v2
        fields["velMag"] = np.sqrt(v2)
        # eps: 文件里有就用文件里的 (求解器写的就是这两个)
        if "rho" in fields and "p" in fields:
            with np.errstate(divide="ignore", invalid="ignore"):
                v = np.sqrt(v2)
                if "eps_M" not in fields:
                    M = v/np.sqrt(gamma*p_eff/fields["rho"])   # 用打过地板的 p
                    fields["eps_M"] = np.clip(1.0 - _switch_s(M), EPS_MIN, EPS_MAX)
                    print("[eps] 文件里没有 eps_M, 用公式重算")
                if "eps_B" not in fields and bx and by:
                    vA = np.maximum(np.sqrt(fields["Bmag"]**2/fields["rho"]), 1e-12)
                    fields["eps_B"] = np.clip(1.0 - _switch_s(v/vA), EPS_MIN, EPS_MAX)
                    print("[eps] 文件里没有 eps_B, 用公式重算")
    for src, dst in (("rho", "log_rho"), ("p", "log_p")):
        if src not in fields:
            continue
        v = p_eff if (src == "p" and p_eff is not None) else np.asarray(fields[src], float)
        bad = ~(v > 0.0)
        if bad.any():
            print(f"[警告] {src}: {int(bad.sum())} 个格子 ({100.0*bad.mean():.3f}%) <= 0, "
                  f"min = {np.nanmin(v):.6g} -> log10 无定义, 画成 {CLIP_BAD}")
        with np.errstate(divide="ignore", invalid="ignore"):
            fields[dst] = np.where(bad, np.nan, np.log10(np.where(bad, 1.0, v)))
    return fields


# ========================= 画 图 =========================

def _cmap_for(name):
    try:
        cm = plt.get_cmap(CMAP).copy()
    except AttributeError:
        import copy; cm = copy.copy(plt.get_cmap(CMAP))
    cm.set_bad(CLIP_BAD)
    if name in FIELD_RANGE:
        cm.set_under(CLIP_UNDER); cm.set_over(CLIP_OVER)
    return cm


def _range_levels(name, v):
    if name in FIELD_RANGE:
        vmin, vmax = FIELD_RANGE[name]
    else:
        fin = np.isfinite(v)
        vmin, vmax = (float(v[fin].min()), float(v[fin].max())) if fin.any() else (0.0, 1.0)
    if vmax <= vmin:
        pad = 0.5*abs(vmin) if vmin else 0.5
        vmin, vmax = vmin - pad, vmax + pad
    k = FIELD_LEVELS.get(name, LEVELS)
    return vmin, vmax, (np.linspace(vmin, vmax, k) if k and k > 1 else None)


def _extend(name, v, vmin, vmax):
    v = np.asarray(v, float); fin = np.isfinite(v)
    if (~fin).any():
        print(f"      {int((~fin).sum())} 个格子 NaN/Inf -> {CLIP_BAD}")
    if name not in FIELD_RANGE:
        return "neither"
    lo, hi = int((v[fin] < vmin).sum()), int((v[fin] > vmax).sum())
    if lo or hi:
        print(f"      固定色标裁掉: {lo} 个 < {vmin:g} ({CLIP_UNDER}), {hi} 个 > {vmax:g} ({CLIP_OVER})")
    return "both" if (lo and hi) else ("min" if lo else ("max" if hi else "neither"))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else FILE
    fields, nx, ny, t = read_vtk(path)
    fields = derive(fields)

    (x0, x1), (y0, y1) = DOMAIN
    xc = x0 + (np.arange(nx) + 0.5)*(x1 - x0)/nx
    yc = y0 + (np.arange(ny) + 0.5)*(y1 - y0)/ny
    X, Y = np.meshgrid(xc, yc)
    dx, dy = (x1 - x0)/nx, (y1 - y0)/ny
    extent = (x0, x1, y0, y1)

    rows, cols = len(PANEL_FIELDS), max(len(r) for r in PANEL_FIELDS)
    fig = plt.figure(figsize=(PANEL_W*cols, PANEL_H*rows))
    grid = ImageGrid(fig, 111, nrows_ncols=(rows, cols), axes_pad=(PAD_W, PAD_H),
                     cbar_mode="each", cbar_location="right", cbar_pad=0.10,
                     cbar_size="4.5%", label_mode="all", aspect=True)

    for r, row in enumerate(PANEL_FIELDS):
        for c, name in enumerate(row):
            k = r*cols + c
            ax, cax = grid[k], grid.cbar_axes[k]
            if name not in fields:
                print(f"[跳过] 没有 field '{name}'. 可用: {', '.join(sorted(fields))}")
                ax.set_axis_off(); cax.set_axis_off(); continue
            v = np.asarray(fields[name], float)
            vmin, vmax, levels = _range_levels(name, v)
            fin = np.isfinite(v)
            print(f"  {name}: data [{v[fin].min():.4g}, {v[fin].max():.4g}], "
                  f"色标 [{vmin:.4g}, {vmax:.4g}]")
            V = v.reshape(ny, nx)
            im = ax.pcolormesh(X, Y, V, shading="gouraud" if SMOOTH else "nearest",
                               cmap=_cmap_for(name), vmin=vmin, vmax=vmax)
            if levels is not None:
                cs = ax.contour(X, Y, np.nan_to_num(V, nan=vmin), levels=levels,
                                colors=LINE_COLOR, linewidths=LINE_WIDTH)
            ax.set_xlabel("X", fontweight="bold", fontsize=18)
            ax.set_ylabel("Y", fontweight="bold", fontsize=18)
            ax.tick_params(axis="both", labelsize=12)
            ax.set_aspect("equal")
            ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
            if TICK_STEP:
                ax.xaxis.set_major_locator(ticker.MultipleLocator(TICK_STEP))
                ax.yaxis.set_major_locator(ticker.MultipleLocator(TICK_STEP))
            cb = fig.colorbar(im, cax=cax, extend=_extend(name, v, vmin, vmax))
            cb.locator = ticker.MaxNLocator(nbins=8)
            fmt = ticker.ScalarFormatter(useMathText=True); fmt.set_powerlimits((-2, 3))
            cb.formatter = fmt; cb.update_ticks()
            cb.ax.tick_params(labelsize=15)
            cax.set_title(CBAR_LABELS.get(name, name), fontsize=18, fontweight="bold", pad=10)
    for k in range(rows*cols, len(grid)):
        grid[k].set_axis_off(); grid.cbar_axes[k].set_axis_off()

    if SAVE:
        out = OUT or os.path.join(os.path.dirname(os.path.abspath(path)), "panel_blast.png")
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        print(f"已保存 {out}")
    if SHOW:
        plt.show()


if __name__ == "__main__":
    main()
