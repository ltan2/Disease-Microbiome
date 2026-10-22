"""Baseline (adversary off) vs DANN (adversary on) from 34_dann.py, for one disease, abundance column, rank and CV.

Reads results/<disease>/ml/dann/{auroc,predictions}_TAG_<baseline|dann>, TAG = <COLUMN>_<RANK>_<CV>; a mode
that has not been run is left out.

Writes results/<disease>/ml/dann/compare_TAG.png, one panel: disease vs control AUROC and, for stratified CV,
the study discriminator's balanced accuracy on the test folds (how much study signal the latent features
still carry). Both are on a 0-1 scale; dots are repeats (stratified) or held-out studies (loso), the bar is
the mean, and each metric has its own dashed chance line (AUROC 0.5, study 1 / #studies). Leave-one-study-out
has no study discriminator score (the held-out study is not among the training domains), so it shows AUROC only.
Usage: 35_plot_dann.py DISEASE COLUMN RANK [CV=stratified]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common

MODES = {"baseline": "Baseline (adversary off)", "dann": "DANN (adversary on)"}
MODE_COLOR = {"baseline": "#6b6a64", "dann": "#2a78d6"}
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"
W = 0.36  # offset between the modes within a metric

disease, column, rank = sys.argv[1], sys.argv[2], sys.argv[3]
cv_mode = sys.argv[4] if len(sys.argv) > 4 else "stratified"
out = common.out_dir(disease, "ml/dann")
tag = f"{column}_{rank}_{cv_mode}"

auroc, n_studies = {}, None
for m in MODES:
    if os.path.exists(f"{out}/auroc_{tag}_{m}.tsv"):
        auroc[m] = pd.read_csv(f"{out}/auroc_{tag}_{m}.tsv", sep="\t")
        n_studies = pd.read_csv(f"{out}/predictions_{tag}_{m}.tsv.gz", sep="\t", usecols=["study"]).study.nunique()
if not auroc:
    raise SystemExit(f"no 34_dann.py output for {tag} in {out}")
modes = list(auroc)
metrics = [("auroc", "Disease vs control\nAUROC", 0.5)]
if cv_mode == "stratified":
    metrics.append(("study_bacc", "Study discriminator\nbalanced accuracy", 1 / n_studies))

fig, ax = plt.subplots(figsize=(3.2 * len(metrics) + 1.6, 4.8))
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color(MUTED)
ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=9)
ax.grid(axis="y", color=GRID, lw=0.8)
ax.set_axisbelow(True)
rng = np.random.default_rng(0)
for j, (col, _, chance) in enumerate(metrics):
    ax.plot([j - 0.45, j + 0.45], [chance] * 2, color=MUTED, ls="--", lw=1)
    ax.text(j + 0.45, chance, f"chance {chance:.2f}", va="bottom", ha="right", fontsize=8, color=MUTED)
    for k, m in enumerate(modes):
        x = j + (k - (len(modes) - 1) / 2) * W
        v = auroc[m][col]
        ax.scatter(x + rng.uniform(-0.05, 0.05, len(v)), v, s=40, color=MODE_COLOR[m], edgecolor="white", lw=1,
                   zorder=3, label=MODES[m] if j == 0 else None)
        ax.plot([x - 0.12, x + 0.12], [v.mean()] * 2, color=INK, lw=2, zorder=4)
        ax.text(x + 0.14, v.mean(), f"{v.mean():.2f}", va="center", fontsize=8, color=INK)
ax.set_xticks(range(len(metrics)), [label for _, label, _ in metrics])
ax.set(xlim=(-0.6, len(metrics) - 0.4), ylim=(0, 1.02), ylabel="score")
ax.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK)
unit = "repeats of 5-fold stratified CV" if cv_mode == "stratified" else "held-out studies"
ax.set_title(f"{disease} ({rank}, {column})\ndots: {unit}; bar: mean", color=INK, fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(f"{out}/compare_{tag}.png", dpi=200)
print(f"wrote {out}/compare_{tag}.png")
