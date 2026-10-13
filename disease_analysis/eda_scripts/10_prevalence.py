"""Prevalence of each feature (taxon or KO) in each group, at each presence level.

prevalence(f, g) = 100 * (# samples in g with f present) / (# samples in g).
Reads results/<disease>/<kind>/presence_<group>.tsv.gz (see common.read_presence).
Writes results/<disease>/<kind>/prevalence.tsv: one row per (level, feature) present in any group.
Usage: 10_prevalence.py DISEASE taxa|func
"""
import sys

import pandas as pd

import common

disease, kind = sys.argv[1], sys.argv[2]
out = common.out_dir(disease, kind)
totals = pd.read_csv(f"{out}/groups.tsv", sep="\t").group.value_counts()

pos = {}
for g in common.GROUPS:
    df, feature, levels = common.read_presence(out, g)
    pos[g] = df.groupby(feature)[levels].sum()
pos = pd.concat(pos, axis=1).fillna(0).astype(int)  # columns: (group, level)

level = common.LEVEL[kind] or "level"
rows = []
for lv in levels:
    df = pd.DataFrame({level: lv}, index=pos.index)
    for g in common.GROUPS:
        n = pos[(g, lv)]
        df[f"{g}: positive/total"] = n.astype(str) + f"/{totals.get(g, 0)}"
        df[f"{g} %"] = (100 * n / totals.get(g, 0)).round(2)
    rows.append(df[pos.xs(lv, axis=1, level=1).sum(axis=1) > 0].sort_values("disease %", ascending=False))
res = pd.concat(rows).reset_index()
if common.LEVEL[kind] is None:
    res = res.drop(columns=level)
res.to_csv(f"{out}/prevalence.tsv", sep="\t", index=False)
