"""Disease vs control classification from YACHT taxa abundance: LASSO (L1 logistic regression) and
random forest (after Wirbel et al. 2019, Nat Med / SIAMCAT).

Samples: disease and control (common.split_groups) with YACHT taxa, pooled across studies: BioProject is
not used anywhere, the only label is disease vs control. other_control is left out (it is not a control
for this disease). Samples with no taxa (relative abundance undefined) are left out.

Features: COLUMN of results/<disease>/taxa/abundance.tsv.gz summed per GTDB r232 clade at RANK (sum
of its genomes' abundances, as in 25_taxa_assoc.py; RANK = genome keeps the YACHT genomes), 0 where
not detected. A feature is kept when present in >= MIN_PREV of the samples, counted over disease and
control together so the filter does not use the labels. Transform: log10(a + pc), pc = half the
smallest nonzero abundance (same pc as 25_taxa_assoc.py).

Models (balanced class weights, since disease:control is often far from 1:1):
  lasso  standardize, then L1 logistic regression; C picked from CS by inner stratified FOLDS-fold CV
         on AUROC within each training set
  rf     random forest, N_TREES trees, sqrt(#features) per split, grown fully
Shadow-feature test for the random forest (after Boruta, Kursa & Rudnicki 2010): in each cv fold a
second forest is fit on the training samples with a shadow copy of every feature added, each shadow
column shuffled independently across samples so it is noise with the same distribution as its feature.
A feature beats the shadows in that fold when its importance exceeds the largest shadow importance,
i.e. more than any pure-noise feature got by chance.
Evaluation: REPEATS x stratified FOLDS-fold CV over all samples, stratified on disease/control only.
Balanced accuracy calls P(disease) >= 0.5 disease.

Writes results/<disease>/ml/, with TAG = <COLUMN>_<RANK>:
  predictions_TAG.tsv.gz  out-of-fold P(disease): acc, group, repeat, model, p
  auroc_TAG.tsv           AUROC and balanced accuracy per model and repeat
  features_TAG.tsv        per feature over the cv folds: LASSO mean coefficient (per SD of log10
                          abundance; > 0 = higher in disease) and selection frequency (fraction of
                          folds with a nonzero coefficient), RF mean impurity importance, fraction of
                          folds where the feature beat the best shadow feature (rf_shadow_hits),
                          prevalence and mean abundance per group
  lasso_folds_TAG.tsv.gz  raw LASSO fit of every cv fold, one row per (fold, feature): repeat, fold,
                          C (chosen by the inner CV), inner_auroc (its inner-CV AUROC), intercept,
                          feature, coef (log-odds of disease per SD of log10 abundance; the scale the
                          model is fit on, 0 = dropped) and coef_per_10x (log-odds per 10-fold higher
                          abundance, = coef / that fold's SD; exp() of it is the odds ratio per 10x)
  roc_TAG.png, features_TAG.png, report_TAG.txt
Usage: 27_classify.py DISEASE COLUMN RANK MIN_PREV PROCS [REPEATS=5]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import common

FOLDS = 5
CS = np.logspace(-3, 1, 13)  # LASSO inverse penalty strengths searched by the inner CV
N_TREES = 500
SEED = 0
TOP = 25  # features per panel in features_TAG.png
MODELS = {"lasso": "LASSO", "rf": "Random forest"}
MODEL_COLOR = {"lasso": "#2a78d6", "rf": "#eb6834"}
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

disease, column, rank, min_prev, procs = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5])
repeats = int(sys.argv[6]) if len(sys.argv) > 6 else 5
if rank != "genome" and rank not in common.RANK_PREFIX.values():
    raise SystemExit(f"RANK must be genome or one of {list(common.RANK_PREFIX.values())}")
taxa = common.out_dir(disease, "taxa")
out = common.out_dir(disease, "ml")
tag = f"{column}_{rank}"
report = []
log = lambda s: (report.append(s), print(s, flush=True))

# ---- samples ----
ab = pd.read_csv(f"{taxa}/abundance.tsv.gz", sep="\t", usecols=["acc", "group", "taxon", column], dtype={"acc": str})
pc = ab[column][ab[column] > 0].min() / 2  # over all groups, as in 25_taxa_assoc.py
groups = pd.read_csv(f"{taxa}/groups.tsv", sep="\t", dtype={"acc": str})
groups = groups[groups.group.isin(["disease", "control"])]

log(f"{disease}: disease vs control classification ({column}, rank {rank}, pseudocount {pc:.3g})")
log(f"  labeled samples with a YACHT result: {groups.group.value_counts().to_dict()}")
no_taxa = ~groups.acc.isin(set(ab.acc))
log(f"  excluded, no taxa detected: {groups[no_taxa].group.value_counts().to_dict()}")
samples = groups[~no_taxa].reset_index(drop=True)
y = (samples.group == "disease").to_numpy().astype(int)
log(f"  samples used: {samples.group.value_counts().to_dict()}")
if min(y.sum(), len(y) - y.sum()) < FOLDS:
    raise SystemExit(f"need at least {FOLDS} disease and {FOLDS} control samples")

# ---- features: samples x clades, log10(a + pc) ----
ab = ab[ab.acc.isin(set(samples.acc))]
if rank == "genome":
    ab["feature"] = ab.taxon
else:
    lineage = common.read_gtdb().set_index("taxon")[rank]
    ab["feature"] = lineage.reindex(ab.taxon).to_numpy()
    if ab.feature.isna().any():
        raise SystemExit(f"{ab[ab.feature.isna()].taxon.nunique()} taxa missing from {common.GTDB_TAXONOMY}")
A = ab.groupby(["acc", "feature"])[column].sum().unstack(fill_value=0).reindex(samples.acc, fill_value=0)
prev = (A > 0).mean()
A = A.loc[:, prev >= min_prev]
log(f"  features: {len(prev)} detected at rank {rank}, {A.shape[1]} present in >= {min_prev:.0%} of samples")
if A.shape[1] == 0:
    raise SystemExit("no feature passes MIN_PREV")
X = np.log10(A.to_numpy() + pc)


# ---- models ----
def lasso():
    lr = LogisticRegression(l1_ratio=1, solver="liblinear", class_weight="balanced", max_iter=5000)
    return GridSearchCV(make_pipeline(StandardScaler(), lr), {"logisticregression__C": CS}, scoring="roc_auc",
                        cv=StratifiedKFold(FOLDS, shuffle=True, random_state=SEED), n_jobs=procs)


def rf():
    return RandomForestClassifier(N_TREES, max_features="sqrt", class_weight="balanced_subsample", n_jobs=procs,
                                  random_state=SEED)


preds, coefs, imps, best_c, shadow_hits, lasso_folds = [], [], [], [], [], []
rng = np.random.default_rng(SEED)


def beats_shadow(train):
    """Per feature: importance above the best shadow (independently shuffled copy) feature's, one fold."""
    shadow = rng.permuted(X[train], axis=0)  # each column shuffled on its own
    imp = rf().fit(np.hstack([X[train], shadow]), y[train]).feature_importances_
    return imp[:X.shape[1]] > imp[X.shape[1]:].max()


def fit_predict(rep, fold, train, test):
    """Fit both models on train, predict test; keep the folds' coefficients and importances."""
    for name, make in [("lasso", lasso), ("rf", rf)]:
        m = make().fit(X[train], y[train])
        p = m.predict_proba(X[test])[:, 1]
        preds.append(samples.iloc[test][["acc", "group"]].assign(repeat=rep, model=name, p=p))
        if name == "lasso":
            scaler, lr = m.best_estimator_
            coefs.append(lr.coef_[0])
            best_c.append(m.best_params_["logisticregression__C"])
            lasso_folds.append(pd.DataFrame({
                "repeat": rep, "fold": fold, "C": best_c[-1], "inner_auroc": m.best_score_, "intercept": lr.intercept_[0],
                "feature": A.columns, "coef": lr.coef_[0], "coef_per_10x": lr.coef_[0] / scaler.scale_}))
        else:
            imps.append(m.feature_importances_)
            shadow_hits.append(beats_shadow(train))


cv = RepeatedStratifiedKFold(n_splits=FOLDS, n_repeats=repeats, random_state=SEED)
for i, (train, test) in enumerate(cv.split(X, y)):
    print(f"  cv fold {i + 1}/{FOLDS * repeats}", flush=True)
    fit_predict(i // FOLDS, i % FOLDS, train, test)
preds = pd.concat(preds, ignore_index=True)
pd.concat(lasso_folds, ignore_index=True).to_csv(f"{out}/lasso_folds_{tag}.tsv.gz", sep="\t", index=False,
                                                float_format="%.6g")
preds[["acc", "group", "repeat", "model", "p"]].to_csv(
    f"{out}/predictions_{tag}.tsv.gz", sep="\t", index=False)


# ---- scores ----
def score(d):
    yy = (d.group == "disease").to_numpy().astype(int)
    return {"n_disease": yy.sum(), "n_control": len(yy) - yy.sum(), "auroc": roc_auc_score(yy, d.p),
            "balanced_accuracy": balanced_accuracy_score(yy, d.p >= 0.5)}


auroc = pd.DataFrame([{"model": model, "repeat": rep, **score(d)} for (model, rep), d in preds.groupby(["model", "repeat"])])
auroc.to_csv(f"{out}/auroc_{tag}.tsv", sep="\t", index=False)

log(f"  LASSO C chosen in the cv folds: median {np.median(best_c):.3g} (range {min(best_c):.3g}-{max(best_c):.3g})")
for model, label in MODELS.items():
    a = auroc[auroc.model == model]
    log(f"  {label:13s} AUROC {a.auroc.mean():.3f} ± {a.auroc.std():.3f}, balanced accuracy "
        f"{a.balanced_accuracy.mean():.3f} ± {a.balanced_accuracy.std():.3f} (mean ± SD over {len(a)} repeats)")

# ---- features ----
coefs, imps, shadow_hits = np.array(coefs), np.array(imps), np.array(shadow_hits)
feat = pd.DataFrame({"feature": A.columns, "lasso_coef": coefs.mean(0), "lasso_freq": (coefs != 0).mean(0),
                     "rf_importance": imps.mean(0), "rf_importance_sd": imps.std(0),
                     "rf_shadow_hits": shadow_hits.mean(0)})
for g in ["disease", "control"]:
    a = A.to_numpy()[(samples.group == g).to_numpy()]
    feat[f"prev_{g}"] = (a > 0).mean(0)
    feat[f"mean_abundance_{g}"] = a.mean(0)
feat = feat.assign(abs_coef=feat.lasso_coef.abs()).sort_values(["lasso_freq", "abs_coef"], ascending=False)
feat = feat.drop(columns="abs_coef")
feat.to_csv(f"{out}/features_{tag}.tsv", sep="\t", index=False)
log(f"  LASSO features selected in >= 50% of folds: {(feat.lasso_freq >= 0.5).sum()} of {len(feat)}")
log(f"  RF features beating the best shadow feature in >= 50% of folds: {(feat.rf_shadow_hits >= 0.5).sum()} of {len(feat)}")


# ---- plots ----
def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=8)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


title = f"{disease}: disease vs control from {rank} abundance ({column})"
fig, ax = plt.subplots(figsize=(5.6, 5.2))
style(ax)
ax.plot([0, 1], [0, 1], color=MUTED, ls="--", lw=1)
for model, label in MODELS.items():
    a = auroc[auroc.model == model]
    for rep, r in preds[preds.model == model].groupby("repeat"):
        fpr, tpr, _ = roc_curve(r.group == "disease", r.p)
        ax.plot(fpr, tpr, color=MODEL_COLOR[model], lw=1.2, alpha=0.5,
                label=f"{label}  AUROC {a.auroc.mean():.2f} ± {a.auroc.std():.2f}" if rep == 0 else None)
ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal", xlabel="False positive rate", ylabel="True positive rate")
ax.set_title(f"{FOLDS}-fold CV x {repeats} repeats", color=INK, fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK)
fig.suptitle(title, color=INK, fontsize=11)
fig.tight_layout()
fig.savefig(f"{out}/roc_{tag}.png", dpi=200)
plt.close(fig)

# top features: LASSO mean coefficient (selection frequency printed) and RF importance
fig, axes = plt.subplots(1, 2, figsize=(14, 0.28 * TOP + 1.6))
sel = feat[feat.lasso_coef != 0].head(TOP)[::-1]
ax = axes[0]
style(ax)
ax.grid(axis="y", visible=False)
ax.barh(range(len(sel)), sel.lasso_coef, color=INK, height=0.7)
ax.axvline(0, color=MUTED, lw=1)
for i, (c, f) in enumerate(zip(sel.lasso_coef, sel.lasso_freq)):
    ax.text(c, i, f" {f:.0%} ", ha="left" if c >= 0 else "right", va="center", fontsize=7, color=MUTED)
ax.set_yticks(range(len(sel)), sel.feature)
ax.set_xlabel("mean LASSO coefficient (> 0: higher in disease); label = fraction of folds selected")
ax.set_title(f"LASSO, top {len(sel)} by selection frequency", color=INK, fontsize=10, loc="left")
top = feat.nlargest(TOP, "rf_importance")[::-1]
ax = axes[1]
style(ax)
ax.grid(axis="y", visible=False)
ax.barh(range(len(top)), top.rf_importance, xerr=top.rf_importance_sd, color=INK, height=0.7,
        error_kw={"ecolor": MUTED, "lw": 1})
ax.set_yticks(range(len(top)), top.feature)
ax.set_xlabel("mean impurity importance (± SD over folds)")
ax.set_title(f"Random forest, top {len(top)} by importance", color=INK, fontsize=10, loc="left")
fig.suptitle(title, color=INK, fontsize=11)
fig.tight_layout()
fig.savefig(f"{out}/features_{tag}.png", dpi=200)
plt.close(fig)

open(f"{out}/report_{tag}.txt", "w").write("\n".join(report) + "\n")
