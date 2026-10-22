"""Per-sample count and count per Mbase against sequencing depth, one point per sample, colored by group.

Left panel: the count saturates with depth (a rarefaction-like curve), so depth drives raw counts.
Right panel: dividing by depth does not remove it; count per Mbase falls as depth rises.
Lines are each group's median in log-spaced depth bins (bins with >= 5 samples); the legend gives
each group's Spearman rho with depth, also written to OUT_spearman.tsv (panel, group, n, rho, p).
Samples with a count of 0 are left out (log axes).

Usage: plot_depth_curve.py IN.tsv COLUMN OUT.png TITLE [--query EXPR] [--ylabel TEXT] [--xlabel TEXT]
                           [--panels raw per_mbase]
  IN.tsv    from 36_per_mbp.py: acc, group, mbases, COLUMN, COLUMN_per_mbp
  --query   pandas query to select rows first (e.g. "coverage == 0.0625" for taxa)
  --panels  which panels to draw (default both; per_mbase alone gives the one-panel plot)
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import common

# same group colors as plot_group_hist.py and 18_metadata_eda.py
GROUP_COLOR = dict(zip(common.GROUPS, ["#2a78d6", "#eb6834", "#1baf7a"]))
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--query")
ap.add_argument("--ylabel")
ap.add_argument("--xlabel", default="sequencing depth (Mbases, before human removal)")
ap.add_argument("--panels", nargs="+", choices=["raw", "per_mbase"], default=["raw", "per_mbase"])
args = ap.parse_args()

df = pd.read_csv(args.tsv, sep="\t")
if args.query:
    df = df.query(args.query)
df = df[(df[args.column] > 0) & (df.mbases > 0)]
ylabel = args.ylabel or args.column
bins = np.logspace(np.log10(df.mbases.min()), np.log10(df.mbases.max()), 16)

rows = []
panels = {"raw": (args.column, ylabel, "raw count"), "per_mbase": (f"{args.column}_per_mbp", f"{ylabel} / Mbase", "per Mbase")}
fig, axes = plt.subplots(1, len(args.panels), figsize=(5.5 * len(args.panels) + (1 if len(args.panels) == 1 else 0), 4.4),
                         squeeze=False)
for ax, (col, lab, name) in zip(axes[0], [panels[p] for p in args.panels]):
    # other_control (largest) first so the smaller groups sit on top
    for grp in reversed(common.GROUPS):
        d = df[df.group == grp]
        ax.scatter(d.mbases, d[col], s=6, alpha=0.25 if grp == "other_control" else 0.5,
                   color=GROUP_COLOR[grp], linewidths=0, rasterized=True)
    for grp in common.GROUPS:
        d = df[df.group == grp]
        rho, p = spearmanr(d.mbases, d[col]) if len(d) > 2 else (np.nan, np.nan)
        rows.append({"panel": "per_mbase" if col != args.column else "raw", "group": grp, "n": len(d), "rho": rho, "p": p})
        b = d.groupby(pd.cut(d.mbases, bins), observed=True)
        med = pd.DataFrame({"x": b.mbases.median(), "y": b[col].median(), "n": b.size()})
        med = med[med.n >= 5]
        ax.plot(med.x, med.y, color=GROUP_COLOR[grp], lw=2,
                label=f"{grp} (n={len(d)}, Spearman rho={rho:.2f})")
    # what dividing by depth assumes: count proportional to depth, at the median count per Mbase
    k, xs = df[f"{args.column}_per_mbp"].median(), np.array([df.mbases.min(), df.mbases.max()])
    ax.plot(xs, k * xs if col == args.column else [k, k], color=MUTED, lw=1.5, ls="--",
            label=f"proportional to depth (what per-Mbase assumes; {k:.3g} / Mbase)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(args.xlabel, color=INK)
    ax.set_ylabel(lab, color=INK)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=7, frameon=False, labelcolor=INK)
    if len(args.panels) > 1:
        ax.set_title(name, color=INK, fontsize=10)
fig.suptitle(args.title, color=INK)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
pd.DataFrame(rows).to_csv(args.out.removesuffix(".png") + "_spearman.tsv", sep="\t", index=False, float_format="%.3g")
