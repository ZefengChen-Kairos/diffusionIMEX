#!/usr/bin/env python3
"""Patch plot_vtk.py: floor the pressure before log10, and stop the figure from
hiding what it cannot plot.

    python3 apply_plot_patch.py /path/to/plot_vtk.py          # in place, keeps a .bak
    python3 apply_plot_patch.py /path/to/plot_vtk.py --dry    # only report

Four changes (background: blastwave_debug/REPORT.md 6.4):

 1. PRESSURE FLOOR before the log, so an existing run can be plotted without being
    re-run.  At beta = 5.1e-4 the thermal energy is 0.128% of the total, so
    p = (g-1)(E - KE - |B|^2/2) is a catastrophic cancellation and the ambient
    pressure just ahead of the fast shock undershoots.  BETA_FLOOR applies
    p <- max(p, BETA_FLOOR * |B|^2/2), the same criterion as the solver's
    betafloor= option, and PRINTS how many cells it lifted.  That count belongs in
    the figure caption: it is a plotting floor, not physics.  Only log_p is floored;
    a raw "p" panel still shows the unfloored data.

 2. log10(|p|) -> masked log10(p).  The script took the ABSOLUTE VALUE before the
    log, so p = -0.005 was drawn as if it were +0.005.  Paper 1's code puts a QUARTER
    of this test's far field below zero and abs() renders every one of those cells as
    a plausible dark blue.  Non-positive cells now become NaN, are drawn white, and
    are counted.

 3. No fixed colour range is imposed, because the floor already does the job: with the
    floor on, log_p lands at [-1.01, 2.30] on its own, so the bar bottom IS the ambient
    pressure.  Commented-out entries are left in FIELD_RANGE for cross-figure
    comparability, with a warning not to set the log_p bottom at exactly -1.0 -- about
    30% of the far field sits at 0.09999, a hair below ambient, and a bound of -1.0
    throws all of it out of range (31% of the cells at 128^2).

 4. A fixed range clips silently - pcolormesh clamps out-of-range cells to the end
    colour, where they look exactly like ambient, which would hide the very thing
    items 1-3 exist to expose.  Below-range cells are now magenta, above-range lime,
    NaN white, the colour bar grows the matching arrow, and the counts are printed.
    A clean figure is one with no magenta, not one whose colour bar was cropped.

 5. gamma was hard-coded to 5/3 while this test runs at 1.4, biasing Mach and any
    recomputed eps.  It is now read back from the conservative variables, with the
    hard-coded value as a fallback.  (eps_M / eps_B are stored in the VTK, so those
    panels were already correct.)
"""
import os, re, shutil, sys

LIT = 'literal'
RGX = 'regex'

EDITS = [
# ------------------------------------------------------------------ 1. gamma
(LIT, """GAMMA = 5.0 / 3.0    # ε_M 和 Mach 都用它: blast / rotor 是 1.4, OT 是 5/3""",
      """GAMMA = 5.0 / 3.0    # 仅作后备: 现在优先从守恒量反推 (detect_gamma); blast/rotor 1.4, OT 5/3""",
 1),

# -------------------------------------------- 2. floor + fixed colour ranges
(LIT, """FIELD_RANGE  = {"eps_M": (0.0, 1.0), "eps_B": (0.0, 1.0)}
FIELD_LEVELS = {"eps_M": 0, "eps_B": 0}""",
      """# ---------------- 取 log 之前的压力地板 ----------------
# 低 beta 下 p = (gamma-1)(E - KE - |B|^2/2) 是灾难性相消: 本算例 beta = 5.129e-4,
# 热能只占总能的 0.128%, 快激波脚处 CT 场与守恒通量之间亚格的相位失配就能把环境
# 压力从 0.1 打到 0.006 (REPORT.md 6.4). 下面在取 log10 之前给 p 一个下界, 这样
# 已经算完的数据可以直接出图, 不必重跑.
#   BETA_FLOOR : p <- max(p, BETA_FLOOR * |B|^2/2)  —— 和求解器 betafloor= 同一个判据.
#                5e-4 略低于本算例的环境 beta (5.129e-4), 即远场被钉在 ~0.0975.
#   P_FLOOR    : 另加一个绝对下界 (两者取较大的那个). 不需要就留 None.
#   两个都 None : 不做地板.
# 只作用于 log_p; 直接画 "p" 那一格仍是未处理的原始数据.
# 被抬起的格子数会打印出来 —— 这是作图处理不是物理, 图注里要写明.
BETA_FLOOR = 5.0e-4
P_FLOOR    = None

# 色标范围. 加了上面的地板之后, log_p 的自然范围就落在 [-1.01, 2.30] —— 底端
# 已经是环境压力, 不需要再固定. 所以默认不给 log_p / log_rho 固定范围.
#
# 如果要跨图/跨分辨率严格可比, 可以打开下面两行; 但注意 log_p 的下界不要正好取
# -1.0: 远场有约 30% 的格子就停在 0.09999 这种"比环境值低一丝"的位置,
# log_p = -1.0000x 会整片落到范围外 (实测 128^2 上会刷掉 31% 的格子).
# 要固定就留点余量, 例如 (-1.05, 2.35); log_rho 的 (-0.6, 0.6) 也会裁掉约 8%.
#   FIELD_RANGE["log_p"]   = (-1.05, 2.35)
#   FIELD_RANGE["log_rho"] = (-0.70, 0.60)
# 被裁掉的格子由 _cmap_for() 画成洋红/亮绿, 不会偷偷藏起来 —— 干净的图是
# "没有洋红"的图, 不是"把 colorbar 下界裁掉"的图.
FIELD_RANGE  = {"eps_M": (0.0, 1.0), "eps_B": (0.0, 1.0)}
FIELD_LEVELS = {"eps_M": 0, "eps_B": 0}

CLIP_UNDER = "magenta"    # 低于 vmin 的格子
CLIP_OVER  = "lime"       # 高于 vmax 的格子
CLIP_BAD   = "white"      # NaN/Inf, 例如 p <= 0 取不了 log10""",
 1),

# ------------------------------- 3. honest log fields (floor + mask, no abs)
(LIT, """    if "rho" in fields:
        with np.errstate(divide="ignore", invalid="ignore"):
            fields["log_rho"] = np.log10(np.abs(fields["rho"]))
    if "p" in fields:
        with np.errstate(divide="ignore", invalid="ignore"):
            fields["log_p"] = np.log10(np.abs(fields["p"]))
    return fields""",
      """    # ---- log 场: 先打地板, 再 log; 绝不用 abs() ----
    # 原来是 np.log10(np.abs(p)): 负压会被画成同样大小的正压, 等于把失效藏起来.
    for src, dst in (("rho", "log_rho"), ("p", "log_p")):
        if src not in fields:
            continue
        v = np.asarray(fields[src], dtype=float)
        if src == "p":
            v = _apply_p_floor(v, fields)
        bad = ~(v > 0.0)
        if bad.any():
            print(f"[警告] {src}: {int(bad.sum())} 个格子 ({100.0*bad.mean():.3f}%) <= 0, "
                  f"min = {np.nanmin(v):.6g} -> log10 无定义, 画成 {CLIP_BAD}")
        with np.errstate(divide="ignore", invalid="ignore"):
            fields[dst] = np.where(bad, np.nan, np.log10(np.where(bad, 1.0, v)))
    return fields""",
 1),

# --------------------- 4. helpers + gamma from the file, before the big entry
(LIT, """def add_derived_fields(fields, x, y, shape, gamma=GAMMA):""",
      """def _apply_p_floor(p, fields):
    \"\"\"取 log 之前给压力一个下界. 返回新数组, 不动 fields["p"].\"\"\"
    floor = None
    if BETA_FLOOR:
        bx, by, bz = _find_B_components(fields)
        if bx is None or by is None:
            print("[地板] 找不到 B 分量, BETA_FLOOR 跳过")
        else:
            mb = 0.5*(fields[bx]**2 + fields[by]**2 + (fields[bz]**2 if bz else 0.0))
            floor = BETA_FLOOR * np.asarray(mb, dtype=float)
    if P_FLOOR:
        pf = np.full_like(np.asarray(p, dtype=float), float(P_FLOOR))
        floor = pf if floor is None else np.maximum(floor, pf)
    if floor is None:
        return p
    out = np.maximum(p, floor)
    lifted = p < floor
    if lifted.any():
        msg = (f"[地板] {int(lifted.sum())} 个格子 ({100.0*lifted.mean():.3f}%) 被抬到下界: "
               f"p_min {np.nanmin(p):.6g} -> {np.nanmin(out):.6g}")
        if BETA_FLOOR:
            msg += f"   (beta >= {BETA_FLOOR:g}, 即 p >= {np.nanmin(floor):.4g}..{np.nanmax(floor):.4g})"
        print(msg)
        print("       ^ 这是作图地板, 不是物理: 图注里要写明, 数量就是这个.")
    else:
        print(f"[地板] 没有格子低于下界 (p_min = {np.nanmin(p):.6g}), 图是原始数据.")
    return out


def detect_gamma(fields, fallback=GAMMA):
    \"\"\"从守恒量反推 gamma: 应当有 p == (gamma-1)(U - KE - |B|^2/2).
    GAMMA 填错 (本算例 1.4 却填 5/3) 会让 Mach 和重算的 eps 整体偏;
    eps_M/eps_B 若已存在 VTK 里则不受影响.\"\"\"
    if not all(k in fields for k in ("p", "U", "rho", "mom_X", "mom_Y", "mom_Z")):
        return fallback
    bx, by, bz = _find_B_components(fields)
    if bx is None or by is None:
        return fallback
    with np.errstate(divide="ignore", invalid="ignore"):
        ke = 0.5*(fields["mom_X"]**2 + fields["mom_Y"]**2 + fields["mom_Z"]**2) \\
             / np.maximum(fields["rho"], 1e-300)
    mb = 0.5*(fields[bx]**2 + fields[by]**2 + (fields[bz]**2 if bz else 0.0))
    e = fields["U"] - ke - mb
    scale = max(1.0, float(np.nanmax(np.abs(fields["p"]))))
    best, berr = fallback, np.inf
    for g in (1.4, 5.0/3.0, 2.0, 1.2, 1.1):
        d = float(np.nanmax(np.abs(fields["p"] - (g - 1.0)*e)))
        if d < berr:
            best, berr = g, d
    if berr < 1e-8 * scale:
        if abs(best - fallback) > 1e-12:
            print(f"[gamma] 从守恒量反推 gamma = {best:.6g} (残差 {berr:.3g}), "
                  f"覆盖脚本里的 GAMMA = {fallback:.6g}")
        return best
    print(f"[gamma] 反推失败 (最小残差 {berr:.3g}), 沿用 GAMMA = {fallback:.6g}")
    return fallback


def _cmap_for(name):
    \"\"\"固定色标时把超范围/无定义的格子画成刺眼颜色, 而不是压到端点色 ——
    否则 FIELD_RANGE 就等于"把缺陷裁掉", 和直接改 colorbar 下界没区别.\"\"\"
    try:
        cm = plt.get_cmap(CMAP).copy()
    except AttributeError:                      # 老 matplotlib
        import copy as _copy
        cm = _copy.copy(plt.get_cmap(CMAP))
    if CLIP_BAD:
        cm.set_bad(CLIP_BAD)
    if name in FIELD_RANGE:
        if CLIP_UNDER:
            cm.set_under(CLIP_UNDER)
        if CLIP_OVER:
            cm.set_over(CLIP_OVER)
    return cm


def _extend_for(name, vals, vmin, vmax):
    \"\"\"colorbar 箭头 + 打印被裁掉多少, 一眼确认图里没有东西被藏起来.\"\"\"
    v = np.asarray(vals, dtype=float)
    fin = np.isfinite(v)
    nbad = int((~fin).sum())
    if nbad:
        print(f"      {nbad} 个格子 NaN/Inf -> {CLIP_BAD}")
    if name not in FIELD_RANGE:
        return "neither"
    nu = int((v[fin] < vmin).sum())
    no = int((v[fin] > vmax).sum())
    if nu or no:
        print(f"      固定色标裁掉: {nu} 个 < {vmin:g} ({CLIP_UNDER}), "
              f"{no} 个 > {vmax:g} ({CLIP_OVER})")
    return "both" if (nu and no) else ("min" if nu else ("max" if no else "neither"))


def add_derived_fields(fields, x, y, shape, gamma=None):
    gamma = detect_gamma(fields) if gamma is None else gamma""",
 1),

# ------------------------------- 5. wire the honest cmap / colorbar into plots
(LIT, """cmap=CMAP,""", """cmap=_cmap_for(name),""", None),
# 缩进无关: 只认 panel/grid 里那两处 (>= 8 空格); _add_cbar 里的 4 空格版本不在 name/vals 作用域内
(RGX, r'^([ \t]{8,})cb = fig\.colorbar\(im, cax=cax\)[ \t]*$',
      lambda m: (f"{m.group(1)}cb = fig.colorbar(im, cax=cax,\n"
                 f"{m.group(1)}                  extend=_extend_for(name, vals, vmin, vmax))"),
 2),
]


def patch(text):
    report = []
    for kind, old, new, want in EDITS:
        if kind == LIT:
            n = text.count(old)
            head = old.strip().splitlines()[0][:66]
        else:
            n = len(re.findall(old, text, re.M))
            head = old[:66]
        if n == 0:
            raise SystemExit("找不到要替换的片段 (你的 plot_vtk.py 可能已改过, 或缩进不同):\n"
                             f"---\n{old[:240]}\n---")
        if want is not None and n != want:
            raise SystemExit(f"片段出现 {n} 次, 期望 {want} 次 —— 先确认是哪几处:\n"
                             f"---\n{old[:240]}\n---")
        text = text.replace(old, new) if kind == LIT else re.sub(old, new, text, flags=re.M)
        report.append((n, head))
    return text, report


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        raise SystemExit(__doc__)
    path = args[0]
    src = open(path, encoding="utf-8").read()
    out, report = patch(src)
    for n, head in report:
        print(f"  x{n}  {head}")
    if "--dry" in sys.argv:
        print(f"\n--dry: {path} 未改动 ({len(src)} -> {len(out)} 字节)")
    else:
        shutil.copy2(path, path + ".bak")
        open(path, "w", encoding="utf-8").write(out)
        print(f"\n已改写 {path}   (备份: {os.path.basename(path)}.bak)")
