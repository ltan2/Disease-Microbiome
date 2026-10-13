"""Disease vs control per-sample counts within BioProjects (other_control never shares a project, so it
is left out).

Counts are the ones plotted: unique hashes (sketch), distinct KOs (func) and distinct taxa (taxa, per
coverage). Only BioProjects with at least MIN_PER_GROUP disease and MIN_PER_GROUP control samples are used.
  per project  Mann-Whitney U (two-sided), p Holm-adjusted over the projects of each data row
  combined     van Elteren test: Wilcoxon rank sums computed within each project, weighted 1/(n+1)
Cliff's delta = P(disease > control) - P(disease < control): >0 means disease counts are higher.
Writes results/<disease>/within_project.tsv.
Usage: 21_within_project.py DISEASE [MIN_PER_GROUP=10]
"""
import os
import sys

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, rankdata

import common

disease = sys.argv[1]
min_n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
out = os.path.dirname(common.out_dir(disease, "func"))

con = duckdb.connect(common.METADATA_DB, read_only=True)
bioproject = con.sql("SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata").df().set_index("acc").bioproject

# data -> (file, count column, facet column)
inputs = {"Sketch (unique hashes)": (f"{out}/sketch/hash_counts.tsv", "n_unique_hashes", None),
          "KO (distinct KOs)": (f"{out}/func/distinct_per_sample.tsv", "n_distinct", None),
          "Taxa (distinct taxa)": (f"{out}/taxa/distinct_per_sample.tsv", "n_distinct", "coverage")}


def holm(p):
    """Holm step-down adjusted p-values."""
    order = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * np.asarray(p)[order]).clip(max=1)
    return adj[np.argsort(order)]


def van_elteren(d, col):
    """z and two-sided p of the van Elteren stratified rank-sum test (strata = bioproject)."""
    T = E = V = 0.0
    for _, s in d.groupby("bioproject"):
        r = rankdata(s[col])
        N, m = len(r), (s.group == "disease").sum()
        _, t = np.unique(r, return_counts=True)
        T += r[(s.group == "disease").to_numpy()].sum() / (N + 1)
        E += m / 2
        V += m * (N - m) * ((N + 1) - (t ** 3 - t).sum() / (N * (N - 1))) / (12 * (N + 1) ** 2)
    z = (T - E) / np.sqrt(V)
    return z, 2 * norm.sf(abs(z))


def tests(d, col):
    rows = []
    for proj, s in d.groupby("bioproject"):
        x, y = s.loc[s.group == "disease", col], s.loc[s.group == "control", col]
        U, p = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({"bioproject": proj, "test": "Mann-Whitney U", "disease": len(x), "control": len(y),
                     "median disease": x.median(), "median control": y.median(), "statistic": U,
                     "cliffs delta": 2 * U / (len(x) * len(y)) - 1, "p-value": p})
    for row, padj in zip(rows, holm([r["p-value"] for r in rows])):
        row["p-adj (Holm)"] = padj
    z, p = van_elteren(d, col)
    x, y = d.loc[d.group == "disease", col], d.loc[d.group == "control", col]
    rows.append({"bioproject": "all (combined)", "test": "van Elteren", "disease": len(x), "control": len(y),
                 "median disease": x.median(), "median control": y.median(), "statistic": z,
                 "cliffs delta": np.nan, "p-value": p, "p-adj (Holm)": p})
    return pd.DataFrame(rows)


res = []
for data, (path, col, facet) in inputs.items():
    if not os.path.exists(path):
        print(f"skipping {data}: no {path}")
        continue
    df = pd.read_csv(path, sep="\t", dtype={"acc": str, "coverage": str})
    df = df[df.group.isin(["disease", "control"])].assign(bioproject=lambda d: bioproject.reindex(d.acc).to_numpy())
    counts = pd.crosstab(df.drop_duplicates("acc").bioproject, df.drop_duplicates("acc").group)
    keep = counts.index[(counts.get("disease", 0) >= min_n) & (counts.get("control", 0) >= min_n)]
    df = df[df.bioproject.isin(keep)]
    if df.empty:
        print(f"skipping {data}: no BioProject with >= {min_n} disease and control samples")
        continue
    for level, d in df.groupby(facet, sort=False) if facet else [("", df)]:
        res.append(tests(d, col).assign(data=data, coverage=level))

res = pd.concat(res)
res = res[["data", "coverage"] + [c for c in res.columns if c not in ("data", "coverage")]]
res.to_csv(f"{out}/within_project.tsv", sep="\t", index=False)
print(res.to_string(index=False))
