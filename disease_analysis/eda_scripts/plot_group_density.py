"""One density curve (Gaussian KDE) per group for a per-sample count column.

Usage: plot_group_density.py IN.tsv COLUMN OUT.png TITLE [--log] [--pseudo P] [--facet COL] [--xlabel TEXT]
  --log    plot log10(value + P), for counts spanning orders of magnitude
  --pseudo pseudocount P added before the log (default 1); with 0, values <= 0 are left out
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
ap.add_argument("--pseudo", type=float, default=1)
ap.add_argument("--facet")
ap.add_argument("--xlabel")
args = ap.parse_args()

df = pd.read_csv(args.tsv, sep="\t", dtype={args.facet: str} if args.facet else None)
if args.log:
    df = df[df[args.column] + args.pseudo > 0]
    df["x"] = np.log10(df[args.column] + args.pseudo)
else:
    df["x"] = df[args.column]
panels = list(df.groupby(args.facet, sort=False)) if args.facet else [(None, df)]

ncol = min(len(panels), 4)
nrow = -(-len(panels) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.8 * nrow), squeeze=False)
for ax, (key, d) in zip(axes.flat, panels):
    grid = np.linspace(d.x.min(), d.x.max(), 500)
    for grp in common.GROUPS:
        v = d.x[d.group == grp]
        if v.nunique() > 1:  # KDE needs some spread
            # median on the plotted axis, so with --log it reads off the x axis directly
            med = v.median()
            med = f"{med:.2f}" if args.log else f"{med:,.0f}" if med >= 100 else f"{med:.3g}"
            ax.plot(grid, gaussian_kde(v)(grid), label=f"{grp} (n={len(v)}, median={med})")
    xlabel = args.xlabel or args.column
    ax.set_xlabel((f"log10({xlabel} + {args.pseudo:g})" if args.pseudo else f"log10({xlabel})") if args.log else xlabel)
    ax.set_ylabel("density")
    if key is not None:
        ax.set_title(f"{args.facet} = {key}")
    ax.legend(fontsize=7)
for ax in axes.flat[len(panels):]:
    ax.remove()
fig.suptitle(args.title)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
