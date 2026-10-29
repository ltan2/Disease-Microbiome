"""Box-and-whisker plot of each group x group block of a sample x sample matrix (all samples, no subsampling).

Same pairs as 05_summary.py: within-group blocks use off-diagonal pairs (each pair once), between-group
blocks use all pairs. Box = quartiles, line = median, diamond = mean, whiskers = 2.5th-97.5th percentile
(--whis lo,hi to change); points beyond the whiskers are not drawn (blocks hold millions of pairs).
Statistics are computed exactly with numpy and drawn with Axes.bxp. --log draws a log y axis (pairs <= 0
are left out of that plot and their share noted under the block).
Usage: 38_block_boxplot.py MATRIX.npy GROUPS.tsv OUT.png TITLE YLABEL [--whis LO,HI] [--log]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import common

ap = argparse.ArgumentParser()
for a in ["matrix", "groups", "out", "title", "ylabel"]:
    ap.add_argument(a)
ap.add_argument("--whis", default="2.5,97.5")
ap.add_argument("--log", action="store_true")
args = ap.parse_args()
lo, hi = map(float, args.whis.split(","))

groups = pd.read_csv(args.groups, sep="\t").group.to_numpy()
M = np.load(args.matrix, mmap_mode="r")
SHORT = {"disease": "D", "control": "C", "other_control": "O"}
WITHIN, BETWEEN = "#2a78d6", "#eb6834"


def pairs(n):
    """Compact pair count (e.g. 4.85M) so the tick labels do not run into each other."""
    return f"{n / 1e6:.3g}M" if n >= 1e6 else f"{n / 1e3:.3g}K" if n >= 1e3 else str(n)


# within-group blocks first, then between-group blocks
blocks = [(g, g) for g in common.GROUPS] + [(a, b) for i, a in enumerate(common.GROUPS)
                                             for b in common.GROUPS[i + 1:]]
stats, colors, notes = [], [], []
for a, b in blocks:
    block = M[np.ix_(groups == a, groups == b)]
    v = block[np.triu_indices(len(block), 1)] if a == b else block.ravel()
    v = v.astype(np.float64)
    n_all = v.size
    if args.log:
        v = v[v > 0]
    q = np.percentile(v, [lo, 25, 50, 75, hi])
    stats.append({"whislo": q[0], "q1": q[1], "med": q[2], "q3": q[3], "whishi": q[4], "mean": v.mean(),
                  "fliers": [], "label": f"{SHORT[a]} x {SHORT[b]}\n{pairs(n_all)} pairs"})
    colors.append(WITHIN if a == b else BETWEEN)
    notes.append(f"{1 - v.size / n_all:.1%} <= 0" if args.log else "")
    print(f"{a} x {b}: n={n_all} median={q[2]:.4g} mean={v.mean():.4g} IQR=[{q[1]:.4g}, {q[3]:.4g}]")

fig, ax = plt.subplots(figsize=(9, 5))
bp = ax.bxp(stats, showmeans=True, showfliers=False, patch_artist=True, widths=0.55,
            medianprops={"color": "black", "lw": 2},
            meanprops={"marker": "D", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": 6},
            whiskerprops={"lw": 1.2}, capprops={"lw": 1.2})
for box, c in zip(bp["boxes"], colors):
    box.set(facecolor=c, alpha=0.75, edgecolor="black", lw=1)
if args.log:
    ax.set_yscale("log")
    for i, note in enumerate(notes, 1):
        ax.annotate(note, (i, 0), xycoords=("data", "axes fraction"), xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8, color="#52514e")
ax.axvline(3.5, color="#c3c2b7", lw=1, ls="--")
ax.set_ylabel(args.ylabel)
ax.set_title(args.title)
ax.grid(axis="y", color="#e5e4df", lw=0.8)
ax.set_axisbelow(True)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.legend(handles=[Patch(facecolor=WITHIN, alpha=0.75, label="within group"),
                   Patch(facecolor=BETWEEN, alpha=0.75, label="between groups"),
                   plt.Line2D([], [], color="black", lw=2, label="median"),
                   plt.Line2D([], [], marker="D", ls="", markerfacecolor="white", markeredgecolor="black",
                              label="mean")],
          loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
fig.text(0.01, -0.02, f"D = disease, C = control, O = other control. Box = IQR, whiskers = {lo:g}th-{hi:g}th "
         "percentile; all samples, outliers not drawn.", fontsize=8, color="#52514e")
fig.savefig(args.out, dpi=200, bbox_inches="tight")
