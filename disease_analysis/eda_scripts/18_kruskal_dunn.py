"""Kruskal-Wallis test across the three groups, then pairwise Dunn tests, for a per-sample count.

Dunn z uses the mean ranks of the pooled data with the tie correction; two-sided p-values are
Holm-adjusted over the three pairs. With --facet, everything is run separately per facet value.
Writes OUT.tsv: one Kruskal-Wallis row (groups = all three) and one Dunn row per pair, with
n, the test statistic (H or z), p-value and Holm-adjusted p-value.
Usage: 18_kruskal_dunn.py IN.tsv COLUMN OUT.tsv [--facet COL]
"""
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.stats import kruskal, norm, rankdata

import common

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "out"]:
    ap.add_argument(a)
ap.add_argument("--facet")
args = ap.parse_args()


def holm(p):
    """Holm step-down adjusted p-values."""
    order = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * np.asarray(p)[order]).clip(max=1)
    return adj[np.argsort(order)]


def tests(d):
    x = {g: d.loc[d.group == g, args.column].to_numpy() for g in common.GROUPS}
    H, p = kruskal(*x.values())
    rows = [{"test": "Kruskal-Wallis", "groups": " vs ".join(common.GROUPS), "n": len(d), "statistic": H,
             "p-value": p, "p-adj (Holm)": p}]

    r = rankdata(d[args.column])
    N = len(r)
    _, t = np.unique(r, return_counts=True)
    var = N * (N + 1) / 12 - (t ** 3 - t).sum() / (12 * (N - 1))  # rank variance with tie correction
    mean_rank = pd.Series(r).groupby(d.group.to_numpy()).mean()
    dunn = []
    for a, b in itertools.combinations(common.GROUPS, 2):
        z = (mean_rank[a] - mean_rank[b]) / np.sqrt(var * (1 / len(x[a]) + 1 / len(x[b])))
        dunn.append({"test": "Dunn", "groups": f"{a} vs {b}", "n": len(x[a]) + len(x[b]), "statistic": z,
                     "p-value": 2 * norm.sf(abs(z))})
    for row, padj in zip(dunn, holm([row["p-value"] for row in dunn])):
        row["p-adj (Holm)"] = padj
    return pd.DataFrame(rows + dunn)


df = pd.read_csv(args.tsv, sep="\t", dtype={args.facet: str} if args.facet else None)
if args.facet:
    res = pd.concat(tests(d).assign(**{args.facet: k}) for k, d in df.groupby(args.facet, sort=False))
    res = res[[args.facet] + list(res.columns[:-1])]
else:
    res = tests(df)
res.to_csv(args.out, sep="\t", index=False)
print(res.to_string(index=False))
