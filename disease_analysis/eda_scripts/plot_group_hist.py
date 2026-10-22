"""Grouped histogram per group of one or more per-observation value columns, one panel per column.

Each group's bars are the proportion of that group's rows (e.g. sample x taxon observations) in
each bin, so groups of different sizes are comparable. With --per-sample each row is weighted by
1 / (rows of its sample), so every sample counts equally and a bar is the mean over the group's
samples of the share of that sample's rows in the bin. Bins are shared by all groups of a panel.
The legend gives samples (distinct --sample-col values) and observations (rows) per group.

Usage: plot_group_hist.py IN.tsv[.gz] OUT.png TITLE --columns COL [COL ...] [--labels TEXT ...]
                          [--bins N] [--linear] [--per-sample] [--sample-col COL] [--ylabel TEXT]
  --columns  value columns, one panel each (rows need a 'group' column with common.GROUPS values)
  --labels   x-axis label per column (default: the column name)
  --bins     number of bins (default 20)
  --linear   evenly spaced bins; default is log10-spaced bins, for values spanning orders of
             magnitude (values <= 0 are then dropped and counted in the legend)
  --per-sample  weight rows so each sample counts equally (default: each row counts equally)
  --sample-col  sample id column (default acc)
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common

# same group colors as 18_metadata_eda.py
GROUP_COLOR = dict(zip(common.GROUPS, ["#2a78d6", "#eb6834", "#1baf7a"]))
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

ap = argparse.ArgumentParser()
for a in ["tsv", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--columns", nargs="+", required=True)
ap.add_argument("--labels", nargs="+")
ap.add_argument("--bins", type=int, default=20)
ap.add_argument("--linear", action="store_true")
ap.add_argument("--per-sample", action="store_true")
ap.add_argument("--sample-col", default="acc")
ap.add_argument("--ylabel")
args = ap.parse_args()
labels = args.labels or args.columns
if len(labels) != len(args.columns):
    ap.error("give one --labels entry per --columns entry")

ylabel = args.ylabel or ("mean proportion per sample" if args.per_sample else "proportion of observations")

df = pd.read_csv(args.tsv, sep="\t", usecols=["group", args.sample_col, *args.columns])
groups = [g for g in common.GROUPS if g in set(df.group)]

fig, axes = plt.subplots(1, len(args.columns), figsize=(6.5 * len(args.columns), 4.2), squeeze=False)
for ax, col, label in zip(axes.flat, args.columns, labels):
    v = df[col] if args.linear else df[col][df[col] > 0]
    lo, hi = v.min(), v.max()
    edges = np.linspace(lo, hi, args.bins + 1) if args.linear else np.geomspace(lo, hi, args.bins + 1)
    # bars side by side within each bin, on an index axis so log bins get equal width
    width = 0.8 / len(groups)
    for i, g in enumerate(groups):
        x = v[df.group[v.index] == g]
        n_all = (df.group == g).sum()
        sample = df[args.sample_col][x.index]
        w = 1 / sample.map(sample.value_counts()) if args.per_sample else np.ones(len(x))
        prop = np.histogram(x, edges, weights=w)[0] / w.sum()
        dropped = "" if args.linear or n_all == len(x) else f", {n_all - len(x):,} <= 0 dropped"
        ax.bar(np.arange(args.bins) + (i - (len(groups) - 1) / 2) * width, prop, width * 0.9,
               color=GROUP_COLOR[g], label=f"{g} ({sample.nunique():,} samples, {len(x):,} obs{dropped})")
    ticks = np.arange(0, args.bins + 1, max(1, args.bins // 5))
    ax.set_xticks(ticks - 0.5)
    ax.set_xticklabels([f"{edges[t]:.2g}" for t in ticks], color=MUTED, fontsize=8)
    ax.set_xlim(-0.6, args.bins - 0.4)
    ax.set_xlabel(f"{label} ({'linear' if args.linear else 'log'}-spaced bins)", color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.tick_params(axis="y", colors=MUTED, labelsize=8)
    ax.tick_params(length=0)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.legend(fontsize=7, frameon=False)
fig.suptitle(args.title, color=INK)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
