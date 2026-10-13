"""Heatmap of a sample x sample matrix with borders around each group x group block.

Rows/columns follow GROUPS.tsv, which is ordered by group. --n-sub keeps a random (seed 0)
subsample of up to N samples per group. The color scale is capped at the 99th percentile of
off-diagonal values so outliers do not wash out the blocks.
Usage: 04_heatmap.py MATRIX.npy GROUPS.tsv OUT.png TITLE COLORBAR_LABEL [--n-sub N]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common

ap = argparse.ArgumentParser()
for a in ["matrix", "groups", "out", "title", "label"]:
    ap.add_argument(a)
ap.add_argument("--n-sub", type=int)
args = ap.parse_args()

groups = pd.read_csv(args.groups, sep="\t").group.to_numpy()
rng = np.random.default_rng(0)
idx = []
for g in common.GROUPS:
    ids = np.flatnonzero(groups == g)
    if args.n_sub and len(ids) > args.n_sub:
        ids = np.sort(rng.choice(ids, args.n_sub, replace=False))
    idx.append(ids)
sizes = [len(i) for i in idx]
idx = np.concatenate(idx)
M = np.load(args.matrix, mmap_mode="r")[np.ix_(idx, idx)]
vmax = np.percentile(M[~np.eye(len(M), dtype=bool)], 99)

fig, ax = plt.subplots(figsize=(10, 9))
im = ax.imshow(M, cmap="viridis", vmin=0, vmax=vmax, interpolation="antialiased")
edges = np.concatenate([[0], np.cumsum(sizes)])
for e in edges[1:-1]:  # group borders split the matrix into the 3 x 3 blocks
    ax.axhline(e - 0.5, color="white", lw=1.5)
    ax.axvline(e - 0.5, color="white", lw=1.5)
mid = (edges[:-1] + edges[1:]) / 2 - 0.5  # cell i is centered on i
labels = [f"{g}\n(n={s})" for g, s in zip(common.GROUPS, sizes)]
ax.set_xticks(mid, labels)
ax.set_yticks(mid, labels)
ax.set_title(args.title)
fig.colorbar(im, ax=ax, label=args.label, shrink=0.8)
fig.savefig(args.out, dpi=200, bbox_inches="tight")
