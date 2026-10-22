"""Forest plot: one panel per contrast, one row per label, effect with its 95% CI.

Reads a table with label, contrast, effect, ci_lo, ci_hi, q (e.g. 25_taxa_assoc.py output). Rows
keep the table's label order, unless --top picks the N labels with the smallest q in
--sort-contrast, then ordered by their effect there. Filled points: q < --alpha; hollow: not.
When some labels start with whitespace (members listed under a parent, as in markers_*.tsv),
the unindented parent labels are bold with a rule above each.
Usage: plot_forest.py IN.tsv[.gz] OUT.png TITLE [--query EXPR] [--top N] [--sort-contrast C]
                      [--alpha A] [--xlabel TEXT]
  --query  pandas query to select rows first, e.g. "rank == 'species'"
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

ap = argparse.ArgumentParser()
for a in ["tsv", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--query")
ap.add_argument("--top", type=int)
ap.add_argument("--sort-contrast", default="disease vs control")
ap.add_argument("--alpha", type=float, default=0.05)
ap.add_argument("--xlabel", default="log2 fold change of geometric mean abundance (95% CI)")
args = ap.parse_args()

df = pd.read_csv(args.tsv, sep="\t")
if args.query:
    df = df.query(args.query)
if args.top:
    s = df[df.contrast == args.sort_contrast].dropna(subset=["q"])
    s = s.assign(abs_eff=s.effect.abs()).sort_values(["q", "abs_eff"], ascending=[True, False]).head(args.top)
    labels = list(s.sort_values("effect", ascending=False).label)
else:
    labels = list(dict.fromkeys(df.label))
contrasts = list(dict.fromkeys(df.contrast))

if not labels:
    fig, ax = plt.subplots(figsize=(6, 1.5))
    ax.text(0.5, 0.5, "no rows to plot", ha="center", va="center", color=MUTED)
    ax.axis("off")
    fig.suptitle(args.title, color=INK)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    raise SystemExit(0)

pos = {lab: i for i, lab in enumerate(labels)}
nested = any(lab[:1].isspace() for lab in labels)
parents = [i for i, lab in enumerate(labels) if nested and not lab[:1].isspace()]
fig, axes = plt.subplots(1, len(contrasts), figsize=(3.6 * len(contrasts) + 3, 0.26 * len(labels) + 1.6),
                         sharey=True, squeeze=False)
for ax, c in zip(axes.flat, contrasts):
    d = df[(df.contrast == c) & df.label.isin(pos)]
    y = d.label.map(pos)
    sig = d.q < args.alpha
    ax.hlines(y, d.ci_lo, d.ci_hi, color=INK, lw=1.2)
    ax.scatter(d.effect[sig], y[sig], s=28, color=INK, zorder=3, label=f"q < {args.alpha}")
    ax.scatter(d.effect[~sig], y[~sig], s=28, facecolor="white", edgecolor=INK, lw=1.2, zorder=3,
               label=f"q ≥ {args.alpha}")
    ax.axvline(0, color=MUTED, lw=1, ls="--")
    for i in parents[1:]:
        ax.axhline(i - 0.5, color=GRID, lw=1)
    ax.set_title(c, color=INK, fontsize=10)
    ax.set_xlabel(args.xlabel, color=INK, fontsize=8)
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
ax0 = axes.flat[0]
ax0.set_yticks(range(len(labels)))
ax0.set_yticklabels([lab.strip() for lab in labels], color=INK)
for i in parents:
    ax0.get_yticklabels()[i].set_fontweight("bold")
ax0.set_ylim(len(labels) - 0.5, -0.5)
fig.suptitle(args.title, color=INK)
fig.tight_layout()
fig.legend(*axes.flat[-1].get_legend_handles_labels(), fontsize=8, frameon=False, ncol=2,
           loc="lower right", bbox_to_anchor=(1, 1))
fig.savefig(args.out, dpi=200, bbox_inches="tight")
