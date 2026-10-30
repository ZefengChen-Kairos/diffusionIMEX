"""Figure of the pour inputs: flow Q(t), spout height h(t) and the spout path, one row per pattern.

    python -m actions.plot [out.png]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from actions.scripts import ALL, describe          # noqa: E402

INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"


def plot(out):
    names = list(ALL)
    fig, axs = plt.subplots(len(names), 3, figsize=(12, 2.1 * len(names)),
                            gridspec_kw=dict(width_ratios=[2.2, 2.2, 1.0]))
    for r, n in enumerate(names):
        a = ALL[n]()
        t = a.tg[a.tg <= a.T_end]
        k = len(t)
        xy = np.array([a.pos_mm(tt) for tt in t[::10]])
        d = describe(a)
        for c, (ax, y, lab) in enumerate(((axs[r, 0], a.val["Q"][:k], "Q [ml/s]"),
                                          (axs[r, 1], a.val["h"][:k] * 100, "h [cm]"))):
            ax.plot(t, y, color=INK, lw=1.6)
            for lbl, ts in a.strokes:
                ax.axvline(ts, color=GRID, lw=1, zorder=0)
                if c == 0:
                    ax.text(ts + 0.05, 13.2, lbl, fontsize=7, color=MUTED,
                            va="top")
            ax.set_ylabel(lab, fontsize=8)
            ax.tick_params(labelsize=7, colors=MUTED)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.set_xlim(0, a.T_end)
        axs[r, 0].set_ylim(0, 14)
        axs[r, 1].set_ylim(0, 8)
        axs[r, 0].set_title(f"{n}: {d['T_pour']:.1f} s continuous pour, {d['V_ml']:.0f} ml, "
                            f"max |dQ/dt| {d['dQdt_max']:.0f} ml/s$^2$", fontsize=8, loc="left", color=INK)
        axs[r, 1].set_title(f"max |dh/dt| {d['dhdt_max']:.2f} m/s", fontsize=8, loc="left", color=INK)
        ax = axs[r, 2]
        ax.add_patch(plt.Circle((0, 0), 42, fill=False, color=GRID, lw=1))
        on = a.val["Q"][:k][::10] > 0
        ax.plot(xy[:, 0], xy[:, 1], color=GRID, lw=0.8)
        ax.scatter(xy[on, 0], xy[on, 1], s=1.5 + 0.25 * a.val["Q"][:k][::10][on] * (a.val["h"][:k][::10][on] < 0.02),
                   color=INK, lw=0)
        ax.set_xlim(-45, 45)
        ax.set_ylim(-45, 45)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("spout path (pitcher at -y)", fontsize=8, color=INK)
    for ax in axs[-1, :2]:
        ax.set_xlabel("t [s]", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "docs", "figures", "pours.png")
    plot(out)
    print(out)
