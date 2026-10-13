"""One density curve (Gaussian KDE) per group for a per-sample count column.

Usage: plot_group_density.py IN.tsv COLUMN OUT.png TITLE [--log] [--facet COL] [--xlabel TEXT]
  --log    plot log10(value + 1), for counts spanning orders of magnitude
  --facet  one panel per value of COL (e.g. coverage)
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import common

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--log", action="store_true")
ap.add_argument("--facet")
ap.add_argument("--xlabel")
args = ap.parse_args()

df = pd.read_csv(args.tsv, sep="\t", dtype={args.facet: str} if args.facet else None)
df["x"] = np.log10(df[args.column] + 1) if args.log else df[args.column]
panels = list(df.groupby(args.facet, sort=False)) if args.facet else [(None, df)]

ncol = min(len(panels), 4)
nrow = -(-len(panels) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.8 * nrow), squeeze=False)
for ax, (key, d) in zip(axes.flat, panels):
    grid = np.linspace(d.x.min(), d.x.max(), 500)
    for grp in common.GROUPS:
        v = d.x[d.group == grp]
        if v.nunique() > 1:  # KDE needs some spread
            ax.plot(grid, gaussian_kde(v)(grid), label=f"{grp} (n={len(v)}, median={d[args.column][d.group == grp].median():,.0f})")
    xlabel = args.xlabel or args.column
    ax.set_xlabel(f"log10({xlabel} + 1)" if args.log else xlabel)
    ax.set_ylabel("density")
    if key is not None:
        ax.set_title(f"{args.facet} = {key}")
    ax.legend(fontsize=7)
for ax in axes.flat[len(panels):]:
    ax.remove()
fig.suptitle(args.title)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
