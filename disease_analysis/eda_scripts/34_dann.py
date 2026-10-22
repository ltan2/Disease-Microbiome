"""Disease vs control classification from YACHT taxa abundance with a domain-adversarial neural network
(DANN; Ganin et al. 2016, JMLR), where the domain is the study (BioProject). The adversary is a switch:

  ADVERSARY=off  baseline: encoder + disease classifier only. The study discriminator head is still
                 trained, but on detached features (gradient reversal weight 0), so it only measures how
                 much study information the latent features carry and never changes them
  ADVERSARY=on   DANN: the study discriminator's gradient is reversed into the encoder (weight ramped
                 from 0 to LAMBDA over training, alpha = LAMBDA * (2 / (1 + exp(-10 p)) - 1), p = training
                 progress), pushing the encoder toward features that predict disease but not study

Samples, features and pseudocount as in 27_classify.py: disease and control (common.split_groups) with
YACHT taxa, COLUMN of results/<disease>/taxa/abundance.tsv.gz summed per GTDB r232 clade at RANK, kept
when present in >= MIN_PREV of the samples, log10(a + pc), then standardized with the training fold's
mean and SD. Domains: BioProject (SRA metadata); studies with fewer than MIN_STUDY_N samples, and samples
without a BioProject, are pooled into one 'other' domain.

Evaluation (CV):
  stratified  REPEATS x stratified FOLDS-fold CV on disease/control (same splits as 27_classify.py, so
              AUROC is comparable to LASSO/RF). Also reports the study discriminator's balanced accuracy
              on the test fold (chance = 1 / #domains): lower with the adversary on = less study signal
  loso        leave-one-study-out: each study with >= MIN_STUDY_N samples of both disease and control is
              held out in turn and the model is trained on all other samples (the held-out study is not
              seen at all, not even unlabeled). This is the test of whether removing study signal helps
              the model transfer to a new study. REPEATS is ignored

Model (SpeciesDANN): encoder Linear(n, HIDDEN)-BN-ReLU-Dropout-Linear(HIDDEN, LATENT)-BN-ReLU; disease
head Linear(LATENT, 16)-ReLU-Dropout-Linear(16, 1); study head the same with one output per domain.
Losses: BCE with pos_weight = #control / #disease (training fold) and class-weighted cross-entropy over
domains (each study counts equally); Adam, EPOCHS epochs, batch BATCH.

Writes results/<disease>/ml/dann/, with TAG = <COLUMN>_<RANK>_<CV>_<baseline|dann>:
  predictions_TAG.tsv.gz  out-of-fold P(disease): acc, group, study, repeat, fold, p
  auroc_TAG.tsv           stratified: AUROC and balanced accuracy per repeat, plus mean test-fold study
                          discriminator balanced accuracy; loso: the same per held-out study
  history_TAG.tsv.gz      per fold and epoch: alpha, disease loss, study loss (training means)
  features_TAG.tsv        which species the trained model relies on, from its held-out samples in each fold:
                          auroc_drop     permutation importance: test AUROC minus test AUROC with that
                                         feature shuffled across the test samples (mean of N_PERM shuffles),
                                         averaged over folds (_sd: SD over folds); > 0 = the model needs it
                          frac_folds_pos fraction of folds with auroc_drop > 0
                          direction      mean gradient of the disease logit w.r.t. the feature (per SD of
                                         log10 abundance) over test samples and folds; > 0 = higher abundance
                                         pushes toward disease
                          study_drop     stratified only: drop in study discriminator balanced accuracy with
                                         the feature shuffled (> 0 = the feature carries study signal; with the
                                         adversary on it should shrink)
                          prevalence and mean abundance per group; sorted by auroc_drop. Shuffling one of
                          several correlated species understates each of them (the others stand in for it)
  report_TAG.txt
Usage: 34_dann.py DISEASE COLUMN RANK MIN_PREV ADVERSARY(on|off) [CV=stratified] [REPEATS=5] [LAMBDA=1]
"""
import sys

import duckdb
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from torch.autograd import Function

import common

FOLDS = 5
SEED = 0  # same CV splits as 27_classify.py
MIN_STUDY_N = 10
HIDDEN, LATENT = 64, 16
EPOCHS, BATCH, LR, WEIGHT_DECAY = 100, 64, 1e-3, 1e-4
N_PERM = 5  # shuffles per feature for the permutation importance

disease, column, rank, min_prev, adversary = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
cv_mode = sys.argv[6] if len(sys.argv) > 6 else "stratified"
repeats = int(sys.argv[7]) if len(sys.argv) > 7 else 5
lam = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0
if adversary not in ("on", "off"):
    raise SystemExit("ADVERSARY must be on or off")
if cv_mode not in ("stratified", "loso"):
    raise SystemExit("CV must be stratified or loso")
if rank != "genome" and rank not in common.RANK_PREFIX.values():
    raise SystemExit(f"RANK must be genome or one of {list(common.RANK_PREFIX.values())}")
use_adversary = adversary == "on"
taxa = common.out_dir(disease, "taxa")
out = common.out_dir(disease, "ml/dann")
tag = f"{column}_{rank}_{cv_mode}_{'dann' if use_adversary else 'baseline'}"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
report = []
log = lambda s: (report.append(s), print(s, flush=True))

# ---- samples ----
ab = pd.read_csv(f"{taxa}/abundance.tsv.gz", sep="\t", usecols=["acc", "group", "taxon", column], dtype={"acc": str})
pc = ab[column][ab[column] > 0].min() / 2  # over all groups, as in 25_taxa_assoc.py
groups = pd.read_csv(f"{taxa}/groups.tsv", sep="\t", dtype={"acc": str})
groups = groups[groups.group.isin(["disease", "control"])]

log(f"{disease}: disease vs control DANN, adversary {adversary} (lambda {lam if use_adversary else 0}), "
    f"{cv_mode} CV ({column}, rank {rank}, pseudocount {pc:.3g}, device {device})")
no_taxa = ~groups.acc.isin(set(ab.acc))
log(f"  excluded, no taxa detected: {groups[no_taxa].group.value_counts().to_dict()}")
samples = groups[~no_taxa].reset_index(drop=True)
y = (samples.group == "disease").to_numpy().astype(int)
log(f"  samples used: {samples.group.value_counts().to_dict()}")
if min(y.sum(), len(y) - y.sum()) < FOLDS:
    raise SystemExit(f"need at least {FOLDS} disease and {FOLDS} control samples")

# ---- domains: BioProject, small studies and missing BioProject pooled into 'other' ----
bioproject = duckdb.connect(common.METADATA_DB, read_only=True).sql(
    "SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata").df().set_index("acc").bioproject
samples["study"] = bioproject.reindex(samples.acc).fillna("other").to_numpy()
size = samples.study.value_counts()
samples.loc[samples.study.isin(size.index[size < MIN_STUDY_N]), "study"] = "other"
log(f"  studies (domains): {samples.study.nunique()} ({(samples.study == 'other').sum()} samples in 'other', from "
    f"studies with < {MIN_STUDY_N} samples or no BioProject)")
log("  samples per study:\n" + samples.groupby(["study", "group"]).size().unstack(fill_value=0).to_string())

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
X = np.log10(A.to_numpy() + pc).astype(np.float32)


# ---- model ----
class GradientReversalFn(Function):
    @staticmethod
    def forward(ctx, x, alpha):
        ctx.alpha = alpha
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.alpha, None  # alpha = 0: the encoder gets no study gradient


class SpeciesDANN(nn.Module):
    def __init__(self, n_features, n_domains):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_features, HIDDEN), nn.BatchNorm1d(HIDDEN), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(HIDDEN, LATENT), nn.BatchNorm1d(LATENT), nn.ReLU())
        self.disease_head = nn.Sequential(nn.Linear(LATENT, 16), nn.ReLU(), nn.Dropout(0.2), nn.Linear(16, 1))
        self.domain_head = nn.Sequential(nn.Linear(LATENT, 16), nn.ReLU(), nn.Dropout(0.2), nn.Linear(16, n_domains))

    def forward(self, x, alpha):
        z = self.encoder(x)
        return self.disease_head(z).squeeze(1), self.domain_head(GradientReversalFn.apply(z, alpha))


def importance(model, xte, yte, dte, domains):
    """Per feature on the test samples: AUROC drop and study balanced-accuracy drop when the feature is
    shuffled (mean of N_PERM shuffles), and the mean gradient of the disease logit w.r.t. the feature."""
    def scores(x):
        with torch.no_grad():
            d_logit, s_logit = model(x, 0.0)
        return (roc_auc_score(yte, d_logit.cpu().numpy()),
                balanced_accuracy_score(dte, domains[s_logit.argmax(1).cpu().numpy()]) if dte is not None else np.nan)

    auc0, study0 = scores(xte)
    gen = torch.Generator(device=device).manual_seed(SEED)
    drop = np.zeros((xte.shape[1], 2))
    for j in range(xte.shape[1]):
        for _ in range(N_PERM):
            x = xte.clone()
            x[:, j] = x[torch.randperm(len(x), device=device, generator=gen), j]
            auc, study = scores(x)
            drop[j] += [auc0 - auc, study0 - study]
    x = xte.clone().requires_grad_(True)
    model(x, 0.0)[0].sum().backward()
    return drop / N_PERM, x.grad.mean(0).cpu().numpy()


def fit_predict(train, test, fold_id):
    """Train on train, return (P(disease) on test, predicted domain on test, history rows, feature importance)."""
    torch.manual_seed(SEED)
    mu, sd = X[train].mean(0), X[train].std(0) + 1e-6
    xtr = torch.tensor((X[train] - mu) / sd, device=device)
    xte = torch.tensor((X[test] - mu) / sd, device=device)
    ytr = torch.tensor(y[train], dtype=torch.float32, device=device)
    domains, dtr = np.unique(samples.study.to_numpy()[train], return_inverse=True)
    dtr = torch.tensor(dtr, device=device)

    model = SpeciesDANN(X.shape[1], len(domains)).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    disease_loss = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((1 - y[train]).sum() / y[train].sum(), dtype=torch.float32, device=device))
    # balanced: each study weighs the same in the study loss, however many samples it has
    dom_w = len(dtr) / (len(domains) * torch.bincount(dtr, minlength=len(domains)).float())
    domain_loss = nn.CrossEntropyLoss(weight=dom_w)

    n = len(train)
    n_batches = max(n // BATCH, 1)  # drop the last partial batch: BatchNorm fails on a batch of 1
    history = []
    for epoch in range(EPOCHS):
        model.train()
        perm = torch.randperm(n, device=device)
        sums = np.zeros(2)
        for b in range(n_batches):
            idx = perm[b * BATCH:(b + 1) * BATCH] if n >= BATCH else perm
            p = (b + epoch * n_batches) / (EPOCHS * n_batches)
            alpha = lam * (2 / (1 + np.exp(-10 * p)) - 1) if use_adversary else 0.0
            d_logit, s_logit = model(xtr[idx], alpha)
            l_d, l_s = disease_loss(d_logit, ytr[idx]), domain_loss(s_logit, dtr[idx])
            opt.zero_grad()
            (l_d + l_s).backward()  # with alpha = 0 l_s only trains the study head
            opt.step()
            sums += [l_d.item(), l_s.item()]
        history.append({"fold": fold_id, "epoch": epoch + 1, "alpha": alpha, "disease_loss": sums[0] / n_batches,
                        "study_loss": sums[1] / n_batches})
        if epoch == 0 or (epoch + 1) % 25 == 0:
            print(f"    epoch {epoch + 1}/{EPOCHS}  alpha {alpha:.3f}  disease loss {sums[0] / n_batches:.4f}  "
                  f"study loss {sums[1] / n_batches:.4f}", flush=True)

    model.eval()
    with torch.no_grad():
        d_logit, s_logit = model(xte, 0.0)
    pred_domain = domains[s_logit.argmax(1).cpu().numpy()]
    # study drop only where the test studies are among the training domains (stratified CV)
    dte = samples.study.to_numpy()[test] if cv_mode == "stratified" else None
    drop, grad = importance(model, xte, y[test], dte, domains)
    imp = pd.DataFrame({"fold": fold_id, "feature": A.columns, "auroc_drop": drop[:, 0], "study_drop": drop[:, 1],
                        "direction": grad})
    return torch.sigmoid(d_logit).cpu().numpy(), pred_domain, history, imp


# ---- CV ----
preds, histories, study_acc, imps = [], [], [], []
if cv_mode == "stratified":
    splits = [(i // FOLDS, i % FOLDS, tr, te)
              for i, (tr, te) in enumerate(RepeatedStratifiedKFold(n_splits=FOLDS, n_repeats=repeats,
                                                                   random_state=SEED).split(X, y))]
else:
    by_study = samples[samples.study != "other"].groupby(["study", "group"]).size().unstack(fill_value=0)
    held = by_study.index[(by_study.reindex(columns=["disease", "control"], fill_value=0) >= MIN_STUDY_N).all(axis=1)]
    log(f"  held-out studies (>= {MIN_STUDY_N} disease and control): {len(held)}")
    if len(held) == 0:
        raise SystemExit("no study has enough disease and control samples to hold out")
    idx = np.arange(len(samples))
    splits = [(0, s, idx[samples.study != s], idx[samples.study == s]) for s in held]

for k, (rep, fold, train, test) in enumerate(splits):
    print(f"  fold {k + 1}/{len(splits)} (repeat {rep}, fold {fold})", flush=True)
    p, pred_domain, history, imp = fit_predict(train, test, k)
    imps.append(imp)
    preds.append(samples.iloc[test][["acc", "group", "study"]].assign(repeat=rep, fold=fold, p=p))
    histories += history
    if cv_mode == "stratified":
        study_acc.append({"repeat": rep, "study_bacc": balanced_accuracy_score(samples.study.iloc[test], pred_domain)})
preds = pd.concat(preds, ignore_index=True)
preds.to_csv(f"{out}/predictions_{tag}.tsv.gz", sep="\t", index=False)
pd.DataFrame(histories).to_csv(f"{out}/history_{tag}.tsv.gz", sep="\t", index=False, float_format="%.6g")


# ---- scores ----
def score(d):
    yy = (d.group == "disease").to_numpy().astype(int)
    return {"n_disease": yy.sum(), "n_control": len(yy) - yy.sum(), "auroc": roc_auc_score(yy, d.p),
            "balanced_accuracy": balanced_accuracy_score(yy, d.p >= 0.5)}


label = "DANN" if use_adversary else "Baseline"
if cv_mode == "stratified":
    auroc = pd.DataFrame([{"repeat": rep, **score(d)} for rep, d in preds.groupby("repeat")])
    auroc = auroc.merge(pd.DataFrame(study_acc).groupby("repeat").study_bacc.mean().reset_index(), on="repeat")
    log(f"  {label}: AUROC {auroc.auroc.mean():.3f} ± {auroc.auroc.std():.3f}, balanced accuracy "
        f"{auroc.balanced_accuracy.mean():.3f} ± {auroc.balanced_accuracy.std():.3f} (mean ± SD over {len(auroc)} repeats)")
    log(f"  study discriminator balanced accuracy on test folds: {auroc.study_bacc.mean():.3f} "
        f"(chance {1 / samples.study.nunique():.3f}; lower = less study signal in the latent features)")
else:
    auroc = pd.DataFrame([{"study": s, **score(d)} for s, d in preds.groupby("study")])
    log(f"  {label}: per held-out study AUROC mean {auroc.auroc.mean():.3f} ± {auroc.auroc.std():.3f} "
        f"(median {auroc.auroc.median():.3f}, {len(auroc)} studies); balanced accuracy "
        f"{auroc.balanced_accuracy.mean():.3f} ± {auroc.balanced_accuracy.std():.3f}")
    log("  per study:\n" + auroc.to_string(index=False, float_format="%.3f"))
auroc.to_csv(f"{out}/auroc_{tag}.tsv", sep="\t", index=False)

# ---- features the model relies on ----
imps = pd.concat(imps, ignore_index=True)
feat = imps.groupby("feature").agg(auroc_drop=("auroc_drop", "mean"), auroc_drop_sd=("auroc_drop", "std"),
                                   frac_folds_pos=("auroc_drop", lambda v: (v > 0).mean()),
                                   direction=("direction", "mean"), study_drop=("study_drop", "mean"))
if cv_mode != "stratified":
    feat = feat.drop(columns="study_drop")
for g in ["disease", "control"]:
    a = A.to_numpy()[(samples.group == g).to_numpy()]
    feat[f"prev_{g}"] = pd.Series((a > 0).mean(0), index=A.columns)
    feat[f"mean_abundance_{g}"] = pd.Series(a.mean(0), index=A.columns)
feat = feat.sort_values("auroc_drop", ascending=False).reset_index()
feat.to_csv(f"{out}/features_{tag}.tsv", sep="\t", index=False, float_format="%.6g")
top = feat.head(15).assign(direction=lambda d: np.where(d.direction > 0, "disease", "control"))
log(f"  top 15 features by test AUROC drop when shuffled (direction: higher abundance pushes toward):\n"
    + top[["feature", "auroc_drop", "frac_folds_pos", "direction"]].to_string(index=False, float_format="%.4f"))
open(f"{out}/report_{tag}.txt", "w").write("\n".join(report) + "\n")
