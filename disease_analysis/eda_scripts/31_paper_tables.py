"""Publication tables for the disease vs control classifiers, built from the raw per-fold outputs of
27_classify.py. Nothing is copied from its printed report; every number is recomputed here.

  Table 1  performance per model: samples, features in, AUROC mean ± SD over the cv repeats (each repeat:
           out-of-fold predictions of all samples) and a 95% CI from a stratified bootstrap of samples
           (BOOT resamples) on each sample's out-of-fold P(disease) averaged over the repeats; balanced
           accuracy, sensitivity and specificity at P(disease) >= 0.5 (mean ± SD over repeats).
           The per-repeat AUROC is checked against auroc_TAG.tsv.
  Table 2  LASSO: features selected (coef != 0) in >= MIN_FREQ of the cv folds. Odds ratio per 10-fold higher
           abundance = exp(coef_per_10x), median and IQR over the folds that selected the feature.
  Table 3  random forest: features beating the best shadow feature in >= MIN_FREQ of the folds, with mean
           impurity importance ± SD over folds and its rank among all features. Importance has no sign;
           direction is that of the univariate disease vs control log2 fold change (28_rank_taxa.py).
  Table S  samples per BioProject and group that entered the classifier. Every fold mixes studies, so a
           study with only cases or only controls lets the model separate the groups by study.
Tables 2 and 3 give prevalence (fraction of samples with the clade) in disease and control.

Writes results/<disease>/ml/paper/{table1_performance,table2_lasso,table3_rf,tableS_samples}_TAG.tsv
(full precision) and tables_TAG.md (rounded, for reading).
Usage: 31_paper_tables.py DISEASE COLUMN RANK [MIN_FREQ=0.5]
"""
import os
import sys

import duckdb
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

import common

BOOT = 2000
SEED = 0
MODELS = {"lasso": "LASSO", "rf": "Random forest"}

disease, column, rank = sys.argv[1], sys.argv[2], sys.argv[3]
min_freq = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
ml, taxa = common.out_dir(disease, "ml"), common.out_dir(disease, "taxa")
tag = f"{column}_{rank}"
out = f"{ml}/paper"
os.makedirs(out, exist_ok=True)

preds = pd.read_csv(f"{ml}/predictions_{tag}.tsv.gz", sep="\t", dtype={"acc": str})
folds = pd.read_csv(f"{ml}/lasso_folds_{tag}.tsv.gz", sep="\t")
feat = pd.read_csv(f"{ml}/features_{tag}.tsv", sep="\t").set_index("feature")
ranked = pd.read_csv(f"{ml}/ranked_{tag}.tsv", sep="\t").set_index("feature")
reported = pd.read_csv(f"{ml}/auroc_{tag}.tsv", sep="\t").set_index(["model", "repeat"]).auroc
name = lambda f: f.split("__", 1)[-1]


def pm(m, s, d=3):
    return f"{m:.{d}f} ± {s:.{d}f}"


# ---- Table 1: performance ----
rng = np.random.default_rng(SEED)
perf = []
for model, label in MODELS.items():
    p = preds[preds.model == model]
    reps = []
    for rep, d in p.groupby("repeat"):
        yy, call = (d.group == "disease").to_numpy(), d.p.to_numpy() >= 0.5
        auc = roc_auc_score(yy, d.p)
        if not np.isclose(auc, reported[(model, rep)], rtol=0, atol=1e-9):
            raise SystemExit(f"{model} repeat {rep}: AUROC {auc} recomputed, {reported[(model, rep)]} in auroc_{tag}.tsv")
        reps.append({"auroc": auc, "balanced_accuracy": balanced_accuracy_score(yy, call),
                     "sensitivity": call[yy].mean(), "specificity": (~call[~yy]).mean()})
    reps = pd.DataFrame(reps)
    mean_p = p.groupby(["acc", "group"]).p.mean().reset_index()
    yy = (mean_p.group == "disease").to_numpy()
    cases, ctrls = np.flatnonzero(yy), np.flatnonzero(~yy)
    boot = [roc_auc_score(np.r_[np.ones(len(cases)), np.zeros(len(ctrls))],
                          mean_p.p.to_numpy()[np.r_[rng.choice(cases, len(cases)), rng.choice(ctrls, len(ctrls))]])
            for _ in range(BOOT)]
    row = {"model": label, "n_disease": int(yy.sum()), "n_control": int((~yy).sum()), "features_in": len(feat),
           "cv": f"{folds.fold.nunique()}-fold stratified x {reps.shape[0]} repeats",
           "auroc_mean": reps.auroc.mean(), "auroc_sd": reps.auroc.std(),
           "auroc_ci_lo": np.quantile(boot, 0.025), "auroc_ci_hi": np.quantile(boot, 0.975)}
    for k in ["balanced_accuracy", "sensitivity", "specificity"]:
        row[f"{k}_mean"], row[f"{k}_sd"] = reps[k].mean(), reps[k].std()
    row["features_stable"] = int((feat.lasso_freq >= min_freq).sum() if model == "lasso" else (feat.rf_shadow_hits >= min_freq).sum())
    perf.append(row)
perf = pd.DataFrame(perf)
perf.to_csv(f"{out}/table1_performance_{tag}.tsv", sep="\t", index=False)

# ---- Table 2: LASSO ----
n_folds = folds.groupby(["repeat", "fold"]).ngroups
sel = folds[folds.coef != 0].groupby("feature").coef_per_10x
lasso = pd.DataFrame({"folds_selected": sel.size(), "or_per_10x_median": np.exp(sel.median()),
                      "or_per_10x_q1": np.exp(sel.quantile(0.25)), "or_per_10x_q3": np.exp(sel.quantile(0.75))})
lasso["selection_freq"] = lasso.folds_selected / n_folds
if not np.allclose(lasso.selection_freq, feat.lasso_freq.reindex(lasso.index)):
    raise SystemExit(f"selection frequency from lasso_folds_{tag}.tsv.gz disagrees with features_{tag}.tsv")
lasso = lasso[lasso.selection_freq >= min_freq].join(feat[["prev_disease", "prev_control"]])
lasso["direction"] = np.where(lasso.or_per_10x_median > 1, "higher in disease", "lower in disease")
lasso = lasso.assign(_a=np.abs(np.log(lasso.or_per_10x_median))).sort_values(["selection_freq", "_a"], ascending=False)
lasso = lasso.drop(columns="_a").rename_axis("feature").reset_index()
lasso = lasso[["feature", "direction", "selection_freq", "folds_selected", "or_per_10x_median", "or_per_10x_q1",
               "or_per_10x_q3", "prev_disease", "prev_control"]]
lasso.to_csv(f"{out}/table2_lasso_{tag}.tsv", sep="\t", index=False)

# ---- Table 3: random forest ----
feat["rf_rank"] = feat.rf_importance.rank(ascending=False, method="min").astype(int)
rf = feat[feat.rf_shadow_hits >= min_freq].join(ranked[["univariate_log2fc", "univariate_q"]])
rf["direction"] = np.select([rf.univariate_log2fc > 0, rf.univariate_log2fc < 0], ["higher in disease", "lower in disease"], "")
rf = rf.sort_values("rf_importance", ascending=False).rename_axis("feature").reset_index()
rf = rf[["feature", "direction", "rf_rank", "rf_importance", "rf_importance_sd", "rf_shadow_hits", "univariate_log2fc",
         "univariate_q", "prev_disease", "prev_control"]]
rf.to_csv(f"{out}/table3_rf_{tag}.tsv", sep="\t", index=False)

# ---- Table S: samples per study ----
used = preds[["acc", "group"]].drop_duplicates()
bp = duckdb.connect(common.METADATA_DB, read_only=True).sql(
    "SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata").df().set_index("acc").bioproject
used["bioproject"] = bp.reindex(used.acc).fillna("none").to_numpy()
study = pd.crosstab(used.bioproject, used.group).reindex(columns=["disease", "control"], fill_value=0)
study["design"] = np.select([(study.disease > 0) & (study.control > 0), study.disease > 0], ["case+control", "case only"], "control only")
study = study.sort_values(["design", "disease", "control"], ascending=[True, False, False]).reset_index()
study.to_csv(f"{out}/tableS_samples_{tag}.tsv", sep="\t", index=False)


# ---- markdown, rounded ----
def md(df):
    df = df.astype(str)
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * df.shape[1]]
    return "\n".join(lines + ["| " + " | ".join(r) + " |" for r in df.to_numpy()])


pct = lambda x: f"{x:.0%}"
t1 = pd.DataFrame({
    "Model": perf.model, "Disease / control": perf.n_disease.astype(str) + " / " + perf.n_control.astype(str),
    "AUROC (mean ± SD)": [pm(m, s) for m, s in zip(perf.auroc_mean, perf.auroc_sd)],
    "AUROC 95% CI": [f"{lo:.3f}–{hi:.3f}" for lo, hi in zip(perf.auroc_ci_lo, perf.auroc_ci_hi)],
    "Balanced accuracy": [pm(m, s) for m, s in zip(perf.balanced_accuracy_mean, perf.balanced_accuracy_sd)],
    "Sensitivity": [pm(m, s) for m, s in zip(perf.sensitivity_mean, perf.sensitivity_sd)],
    "Specificity": [pm(m, s) for m, s in zip(perf.specificity_mean, perf.specificity_sd)],
    f"Stable features (≥ {min_freq:.0%} of folds)": perf.features_stable.astype(str) + f" of {len(feat)}"})
t2 = pd.DataFrame({
    "Species": lasso.feature.map(name), "Direction": lasso.direction,
    "Selected in folds": [f"{f} / {n_folds} ({pct(s)})" for f, s in zip(lasso.folds_selected, lasso.selection_freq)],
    "OR per 10× abundance, median (IQR)": [f"{m:.2f} ({a:.2f}–{b:.2f})" for m, a, b in
                                           zip(lasso.or_per_10x_median, lasso.or_per_10x_q1, lasso.or_per_10x_q3)],
    "Prevalence disease": lasso.prev_disease.map(pct), "Prevalence control": lasso.prev_control.map(pct)})
t3 = pd.DataFrame({
    "Species": rf.feature.map(name), "Direction": rf.direction, "Importance rank": rf.rf_rank,
    "Importance (mean ± SD)": [f"{m:.4f} ± {s:.4f}" for m, s in zip(rf.rf_importance, rf.rf_importance_sd)],
    "Beats shadow (folds)": rf.rf_shadow_hits.map(pct),
    "Univariate log2FC (q)": [f"{e:.2f} ({q:.2g})" for e, q in zip(rf.univariate_log2fc, rf.univariate_q)],
    "Prevalence disease": rf.prev_disease.map(pct), "Prevalence control": rf.prev_control.map(pct)})
n_studies = (study.design == "case+control").sum()
one_class = study.loc[study.design != "case+control", ["disease", "control"]].to_numpy().sum()
doc = f"""# {disease}: disease vs control classification ({column}, rank {rank})

Features: log10(abundance + pseudocount) of {len(feat)} {rank} clades present in ≥ 10% of samples. Cross-validation pools
all studies; AUROC 95% CI by stratified bootstrap ({BOOT} resamples) of out-of-fold P(disease) averaged over repeats;
accuracy measures at P(disease) ≥ 0.5.

**Confounding check:** {n_studies} of {len(study)} studies have both cases and controls; {one_class} of {len(used)} samples come
from single-group studies (Table S), so the models can separate the groups by study.

## Table 1. Classification performance

{md(t1)}

## Table 2. LASSO: species selected in ≥ {min_freq:.0%} of {n_folds} cv folds

OR > 1: higher odds of disease with higher abundance (per 10-fold increase; median and IQR over the folds selecting it).

{md(t2)}

## Table 3. Random forest: species more important than the best shadow (noise) feature in ≥ {min_freq:.0%} of folds

{md(t3)}

## Table S. Samples per BioProject

{md(study)}
"""
open(f"{out}/tables_{tag}.md", "w").write(doc)
print(f"{disease} ({tag}): {n_studies} of {len(study)} studies case+control, {one_class} of {len(used)} samples from single-group studies")
print(md(t1))
print(f"LASSO {len(lasso)} species, RF {len(rf)} species -> {out}/")
