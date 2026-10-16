"""One ranked table of the clades that separate disease from control, combining three lines of evidence.

Per clade at RANK (species, genus or family), from files already written:
  lasso       27_classify.py features_<COLUMN>_<RANK>.tsv: selection frequency over the cv folds and
              mean coefficient (> 0: higher in disease). Passes when lasso_freq >= LASSO_FREQ.
  rf          same file: random forest mean impurity importance (rf_rank 1 = most important) and the
              fraction of folds where it beat the best shuffled shadow feature. Passes when
              rf_shadow_hits >= RF_HITS, i.e. it is more important than noise in most folds.
  univariate  25_taxa_assoc.py assoc_<COLUMN>.tsv.gz, disease vs control: log2 fold change of geometric
              means (> 0: higher in disease) and BH q. Passes when q < ALPHA.
agree = number of the three that pass (0-3). direction_agrees: LASSO coefficient and univariate effect
have the same sign (empty when the LASSO coefficient is 0). Literature markers (literature_markers.tsv
rows for DISEASE; a genus or species name also covers GTDB's suffixed splits, as in 25_taxa_assoc.py)
are flagged with their expected direction, but do not count toward agree.
Sorted by agree, then lasso_freq, then rf_importance. Only clades used as classifier features (present
in >= MIN_PREV of the samples) appear.

Writes results/<disease>/ml/ranked_<COLUMN>_<RANK>.tsv
Usage: 28_rank_taxa.py DISEASE COLUMN RANK
"""
import sys

import numpy as np
import pandas as pd

import common

LASSO_FREQ = 0.5
RF_HITS = 0.5
ALPHA = 0.05

disease, column, rank = sys.argv[1], sys.argv[2], sys.argv[3]
if rank not in ("species", "genus", "family"):
    raise SystemExit("RANK must be species, genus or family (the ranks 25_taxa_assoc.py tests)")
ml, taxa = common.out_dir(disease, "ml"), common.out_dir(disease, "taxa")

feat = pd.read_csv(f"{ml}/features_{column}_{rank}.tsv", sep="\t")
feat["rf_rank"] = feat.rf_importance.rank(ascending=False, method="min").astype(int)
uni = pd.read_csv(f"{taxa}/assoc_{column}.tsv.gz", sep="\t")
uni = uni[(uni["rank"] == rank) & (uni.contrast == "disease vs control")]
uni = uni.set_index("label")[["effect", "q"]].rename(columns={"effect": "univariate_log2fc", "q": "univariate_q"})
t = feat.join(uni, on="feature")

lit = pd.read_csv("literature_markers.tsv", sep="\t")
lit = lit[lit.disease == disease]
t["literature"], t["marker"] = "", ""
for clade, expected in zip(lit.clade, lit.expected):
    hit = t.feature == clade
    if rank in ("genus", "species"):
        hit |= t.feature.str.startswith(clade + "_")
    t.loc[hit, "literature"] = expected
    t.loc[hit, "marker"] = clade

t["lasso_dir"] = np.select([t.lasso_coef > 0, t.lasso_coef < 0], ["up", "down"], "")
t["univariate_dir"] = np.select([t.univariate_log2fc > 0, t.univariate_log2fc < 0], ["up", "down"], "")
t["direction_agrees"] = np.where(t.lasso_dir == "", "", np.where(t.lasso_dir == t.univariate_dir, "yes", "no"))
t["pass_lasso"] = t.lasso_freq >= LASSO_FREQ
t["pass_rf"] = t.rf_shadow_hits >= RF_HITS
t["pass_univariate"] = t.univariate_q < ALPHA
t["agree"] = t[["pass_lasso", "pass_rf", "pass_univariate"]].sum(axis=1)
t = t.sort_values(["agree", "lasso_freq", "rf_importance"], ascending=False)

cols = ["feature", "agree", "pass_lasso", "pass_rf", "pass_univariate", "lasso_freq", "lasso_coef", "lasso_dir",
        "rf_importance", "rf_rank", "rf_shadow_hits", "univariate_log2fc", "univariate_q", "univariate_dir", "direction_agrees",
        "literature", "prev_disease", "prev_control", "mean_abundance_disease", "mean_abundance_control"]
path = f"{ml}/ranked_{column}_{rank}.tsv"
t[cols].to_csv(path, sep="\t", index=False)

print(f"{disease}: {len(t)} {rank} features ({column}); lasso_freq >= {LASSO_FREQ}, rf_shadow_hits >= {RF_HITS}, q < {ALPHA}")
print("  clades by agree: " + ", ".join(f"{k}/3: {v}" for k, v in t.agree.value_counts().sort_index(ascending=False).items()))
print(f"  of agree 3, LASSO and univariate direction disagree: {((t.agree == 3) & (t.direction_agrees == 'no')).sum()}")
m = t[t.literature != ""]
print(f"  literature markers with a feature: {m.marker.nunique()} of {len(lit)} ({len(m)} features, GTDB splits included)")
if len(m):
    print("    " + m[["feature", "literature", "agree", "lasso_dir", "univariate_dir", "rf_rank", "rf_shadow_hits"]].to_string(index=False)
          .replace("\n", "\n    "))
print(f"  top 20 -> {path}")
print("    " + t.head(20)[["feature", "agree", "lasso_freq", "lasso_dir", "rf_rank", "rf_shadow_hits", "univariate_log2fc", "univariate_q",
                          "literature"]].to_string(index=False, float_format="%.3g").replace("\n", "\n    "))
