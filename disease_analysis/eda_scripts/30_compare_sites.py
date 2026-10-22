"""Before/after the body-site filter: all samples vs each body site, for one disease.

Reads, for 'all' (results/<disease>/) and each SITE (results/<disease>/site_<SITE>/, from 29_site_filter.py):
ml/auroc_TAG.tsv and ml/ranked_TAG.tsv (27_classify.py, 28_rank_taxa.py) and taxa/assoc_<COLUMN>.tsv.gz
(25_taxa_assoc.py), TAG = <COLUMN>_<RANK>.

Writes results/<disease>/site_compare/:
  compare_TAG.tsv          one row per setting: disease and control samples, LASSO and RF AUROC (mean, SD
                           over repeats), features, clades passing each test and all three, literature
                           markers passing >= 2 tests
  compare_species_TAG.tsv  clades passing all three tests in any setting, plus the literature markers:
                           agree, univariate log2FC and q, lasso_freq and lasso_dir, rf_shadow_hits per
                           setting (empty: not a feature there, i.e. below MIN_PREV)
  compare_TAG.png          AUROC per setting and model; literature markers' disease vs control log2FC
                           (95% CI) per setting, filled when q < 0.05
Usage: 30_compare_sites.py DISEASE COLUMN RANK SITE [SITE ...]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import common

MODELS = {"lasso": "LASSO", "rf": "Random forest"}
MODEL_COLOR = {"lasso": "#2a78d6", "rf": "#eb6834"}
SETTING_COLOR = ["#6b6a64", "#1baf7a", "#6250d6", "#eda100"]  # all (before) in neutral, then sites
SETTING_MARKER = ["o", "s", "D", "^"]
ALPHA = 0.05
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

disease, column, rank, sites = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
if not sites:
    raise SystemExit("give at least one SITE")
tag = f"{column}_{rank}"
settings = ["all"] + sites
out = common.out_dir(disease, "site_compare", site="")


def read(setting):
    site = "" if setting == "all" else setting
    ml, taxa = common.out_dir(disease, "ml", site=site), common.out_dir(disease, "taxa", site=site)
    assoc = pd.read_csv(f"{taxa}/assoc_{column}.tsv.gz", sep="\t")
    assoc = assoc[(assoc["rank"] == rank) & (assoc.contrast == "disease vs control")].set_index("label")
    return (pd.read_csv(f"{ml}/auroc_{tag}.tsv", sep="\t"), pd.read_csv(f"{ml}/ranked_{tag}.tsv", sep="\t",
            keep_default_na=False, na_values=[""]), assoc)


data = {s: read(s) for s in settings}

# ---- per-setting summary ----
rows = []
for s, (auroc, ranked, _) in data.items():
    first = auroc[auroc.model == "lasso"].iloc[0]  # every repeat predicts every sample once
    r = {"setting": s, "n_disease": first.n_disease, "n_control": first.n_control, "features": len(ranked)}
    for m in MODELS:
        a = auroc[auroc.model == m].auroc
        r[f"{m}_auroc"], r[f"{m}_auroc_sd"] = a.mean(), a.std()
    for test in ["lasso", "rf", "univariate"]:
        r[f"pass_{test}"] = ranked[f"pass_{test}"].sum()
    r["agree_3"] = (ranked.agree == 3).sum()
    lit = ranked[ranked.literature.notna()]
    r["literature_features"], r["literature_agree_ge2"] = len(lit), (lit.agree >= 2).sum()
    rows.append(r)
summary = pd.DataFrame(rows)[["setting", "n_disease", "n_control", "lasso_auroc", "lasso_auroc_sd", "rf_auroc",
                              "rf_auroc_sd", "features", "pass_lasso", "pass_rf", "pass_univariate", "agree_3",
                              "literature_features", "literature_agree_ge2"]]
summary.to_csv(f"{out}/compare_{tag}.tsv", sep="\t", index=False, float_format="%.4g")
print(f"{disease}: all samples vs {', '.join(sites)} ({column}, {rank})")
print("  " + summary.to_string(index=False, float_format="%.3f").replace("\n", "\n  "))

# ---- per-species, wide over settings ----
keep = set()
for _, ranked, _ in data.values():
    keep |= set(ranked.feature[(ranked.agree == 3) | ranked.literature.notna()])
cols = ["agree", "univariate_log2fc", "univariate_q", "lasso_freq", "lasso_dir", "rf_shadow_hits"]
wide = []
for s, (_, ranked, _) in data.items():
    wide.append(ranked.set_index("feature").reindex(sorted(keep))[cols].add_suffix(f"_{s}"))
lit = pd.concat([r.set_index("feature").literature for _, r, _ in data.values()]).dropna()
species = pd.concat(wide, axis=1).rename_axis("feature").reset_index()
species.insert(1, "literature", species.feature.map(lit[~lit.index.duplicated()]))
order = [f"agree_{s}" for s in settings]
species = species.sort_values(order, ascending=False, na_position="last")
species.to_csv(f"{out}/compare_species_{tag}.tsv", sep="\t", index=False, float_format="%.4g")
n3 = {s: set(r.feature[r.agree == 3]) for s, (_, r, _) in data.items()}
for s in sites:
    print(f"  passing all three: all {len(n3['all'])}, {s} {len(n3[s])}, in both {len(n3['all'] & n3[s])}")


# ---- figure ----
def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=8)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


markers = sorted(lit.index.unique(), key=lambda f: (lit[f] if isinstance(lit[f], str) else lit[f].iloc[0], f))
fig, axes = plt.subplots(1, 2, figsize=(13, max(4.5, 0.42 * len(markers) + 1.8)),
                         gridspec_kw={"width_ratios": [1, 1.6]})

ax = axes[0]
style(ax)
ax.grid(axis="x", visible=False)
ax.axhline(0.5, color=MUTED, ls="--", lw=1)
x = np.arange(len(settings))
for k, (m, label) in enumerate(MODELS.items()):
    mean, sd = summary[f"{m}_auroc"], summary[f"{m}_auroc_sd"].fillna(0)
    ax.errorbar(x + (k - 0.5) * 0.22, mean, yerr=sd, fmt="o", ms=8, color=MODEL_COLOR[m], mec="white", mew=1.5,
                ecolor=MODEL_COLOR[m], elinewidth=2, capsize=0, label=label)
    for xi, v in zip(x + (k - 0.5) * 0.22, mean):
        ax.text(xi, v + 0.025, f"{v:.2f}", ha="center", va="bottom", fontsize=7, color=INK)
ax.set_xticks(x, [f"{s}\n{d} D / {c} C" for s, d, c in zip(summary.setting, summary.n_disease, summary.n_control)])
ax.set(ylim=(0.4, 1.02), ylabel="AUROC (mean ± SD over CV repeats)", xlim=(-0.5, len(settings) - 0.5))
ax.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK)
ax.set_title("Disease vs control classification", color=INK, fontsize=10, loc="left")

ax = axes[1]
style(ax)
ax.grid(axis="y", visible=False)
ax.axvline(0, color=MUTED, lw=1)
pos = {f: i for i, f in enumerate(markers[::-1])}
dy = np.linspace(0.25, -0.25, len(settings)) if len(settings) > 1 else [0]
for j, s in enumerate(settings):
    a = data[s][2].reindex(markers).dropna(subset=["effect"])
    yy = a.index.map(pos).to_numpy() + dy[j]
    sig = (a.q < ALPHA).to_numpy()
    ax.hlines(yy, a.ci_lo, a.ci_hi, color=SETTING_COLOR[j], lw=2)
    ax.scatter(a.effect[sig], yy[sig], marker=SETTING_MARKER[j], s=40, color=SETTING_COLOR[j], edgecolor="white",
               lw=1, zorder=3)
    ax.scatter(a.effect[~sig], yy[~sig], marker=SETTING_MARKER[j], s=40, facecolor="white", edgecolor=SETTING_COLOR[j],
               lw=1.5, zorder=3)
ax.set_yticks(list(pos.values()), [f"{f}  ({lit[f] if isinstance(lit[f], str) else lit[f].iloc[0]})" for f in pos])
ax.set_xlabel("disease vs control log2 fold change (95% CI); filled: q < 0.05")
handles = [Line2D([], [], marker=SETTING_MARKER[j], ls="-", lw=2, ms=7, color=SETTING_COLOR[j], mec="white",
                  label=f"{s} ({summary.n_disease[j]} D / {summary.n_control[j]} C)") for j, s in enumerate(settings)]
ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, 1.0), ncol=len(settings),
          labelcolor=INK, borderaxespad=0)
ax.set_ylim(-0.6, len(pos) + 0.4)  # room for the legend above the top marker
ax.set_title("Literature markers (expected direction)", color=INK, fontsize=10, loc="left")
fig.suptitle(f"{disease}: all samples vs body site ({rank}, {column})", color=INK, fontsize=11)
fig.tight_layout()
fig.savefig(f"{out}/compare_{tag}.png", dpi=200)
print(f"  -> {out}/compare_{tag}.png")
