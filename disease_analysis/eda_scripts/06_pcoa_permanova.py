"""PERMANOVA and PERMDISP between groups, plus a PCoA plot, for any sample x sample matrix.

PERMANOVA tests whether group centroids differ; PERMDISP tests whether group spreads differ (a
significant PERMANOVA with a significant PERMDISP may come from spread alone). Both are run on
all three groups and on each pair.
PERMANOVA's pseudo-F is the same as skbio's, but the permutations are computed in batches as matrix
products (onehot' D^2 onehot), which is much faster on large matrices.
PERMDISP is the centroid version (as R vegan::betadisper type="centroid"): each sample's distance
to its group centroid, taken straight from the distance matrix, compared with a permuted ANOVA F.
Writes OUT.tsv (test, groups, sample size, F, p-value) and OUT.png.
Usage: 06_pcoa_permanova.py MATRIX.npy GROUPS.tsv OUT.tsv OUT.png TITLE PERMUTATIONS [--similarity]
  --similarity  MATRIX holds similarities (e.g. Jaccard); distance = 1 - similarity
"""
import argparse
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa

import common

ap = argparse.ArgumentParser()
for a in ["matrix", "groups", "out_tsv", "out_png", "title"]:
    ap.add_argument(a)
ap.add_argument("permutations", type=int)
ap.add_argument("--similarity", action="store_true")
args = ap.parse_args()

g = pd.read_csv(args.groups, sep="\t", index_col=0).group  # first column: sample id
D = np.load(args.matrix).astype(np.float64)
if args.similarity:
    D = 1 - D
np.fill_diagonal(D, 0)
dm = DistanceMatrix(D, g.index, validate=False)  # symmetric by construction; validating 15k x 15k is slow


def anova_f(z, codes):
    """One-way ANOVA F of values z across groups coded 0..k-1."""
    k, n = codes.max() + 1, len(z)
    cnt, s = np.bincount(codes, minlength=k), np.bincount(codes, z, minlength=k)
    between = (s ** 2 / cnt).sum() - z.sum() ** 2 / n
    within = (z ** 2).sum() - (s ** 2 / cnt).sum()
    return (between / (k - 1)) / (within / (n - k))


def permanova(pos, labels, perms, batch=100):
    """Pseudo-F = (SS_between / (k-1)) / (SS_within / (n-k)), with SS_within = sum over groups of
    (sum of squared within-group distances) / group size; p-value from label permutations."""
    codes = pd.factorize(labels)[0]
    k, n = codes.max() + 1, len(codes)
    D2 = D[np.ix_(pos, pos)] ** 2
    ss_total = D2.sum() / (2 * n)
    size = np.bincount(codes)

    def f_stats(perm_codes):  # perm_codes: P x n
        onehot = np.zeros((n, len(perm_codes) * k))
        onehot[np.arange(n)[None, :], perm_codes + k * np.arange(len(perm_codes))[:, None]] = 1
        q = (onehot * (D2 @ onehot)).sum(axis=0).reshape(-1, k)  # 1_g' D^2 1_g per permutation and group
        ss_within = (q / (2 * size)).sum(axis=1)
        return ((ss_total - ss_within) / (k - 1)) / (ss_within / (n - k))

    F = f_stats(codes[None, :])[0]
    rng = np.random.default_rng(0)
    null = np.concatenate([f_stats(np.array([rng.permutation(codes) for _ in range(min(batch, perms - i))]))
                           for i in range(0, perms, batch)])
    return {"sample size": n, "test statistic": F, "p-value": ((null >= F).sum() + 1) / (perms + 1)}


def permdisp(pos, labels, perms):
    """Distance to group centroid: z_i^2 = mean_j d_ij^2 - mean_jk d_jk^2 / 2 over the group's samples
    (abs() for non-Euclidean distances, as betadisper); then F and its permutation p-value."""
    codes = pd.factorize(labels)[0]
    z = np.empty(len(pos))
    for c in np.unique(codes):
        m = pos[codes == c]
        d2 = D[np.ix_(m, m)] ** 2
        z[codes == c] = np.sqrt(np.abs(d2.mean(axis=1) - d2.mean() / 2))
    F = anova_f(z, codes)
    rng = np.random.default_rng(0)
    null = np.array([anova_f(z, rng.permutation(codes)) for _ in range(perms)])
    return {"sample size": len(z), "test statistic": F, "p-value": ((null >= F).sum() + 1) / (perms + 1)}


res = []
for t in [tuple(common.GROUPS)] + list(itertools.combinations(common.GROUPS, 2)):
    keep = g.isin(t).to_numpy()
    pos, labels = np.flatnonzero(keep), g[keep].to_numpy()
    for name, r in [("PERMANOVA", permanova(pos, labels, args.permutations)),
                    ("PERMDISP", permdisp(pos, labels, args.permutations))]:
        res.append({"test": name, "groups": " vs ".join(t), "sample size": r["sample size"],
                    "F": r["test statistic"], "p-value": r["p-value"]})
res = pd.DataFrame(res)
res.to_csv(args.out_tsv, sep="\t", index=False)
print(res.to_string(index=False))

ord_ = pcoa(dm, method="fsvd", dimensions=2, seed=0)  # fast randomized PCoA: only two axes are plotted
pc, var = ord_.samples, ord_.proportion_explained
fig, ax = plt.subplots(figsize=(7, 6))
for grp in common.GROUPS:
    m = (g == grp).to_numpy()
    ax.scatter(pc.iloc[m, 0], pc.iloc[m, 1], s=4, alpha=0.5, label=f"{grp} (n={m.sum()})")
ax.set_xlabel(f"PC1 ({var.iloc[0]:.1%})")
ax.set_ylabel(f"PC2 ({var.iloc[1]:.1%})")
F, p = res.iloc[0][["F", "p-value"]]
Fd, pd_ = res.iloc[1][["F", "p-value"]]
ax.set_title(f"{args.title}\nPERMANOVA F={F:.2f}, p={p:.3g}; PERMDISP F={Fd:.2f}, p={pd_:.3g} (3 groups)")
ax.legend(markerscale=3)
fig.savefig(args.out_png, dpi=200, bbox_inches="tight")
