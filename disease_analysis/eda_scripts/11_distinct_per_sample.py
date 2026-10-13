"""Number of distinct features (taxa or KOs) per sample at each presence level (0 if none).

Writes results/<disease>/<kind>/distinct_per_sample.tsv (acc, group, [level], n_distinct);
plotted by plot_group_density.py.
Usage: 11_distinct_per_sample.py DISEASE taxa|func
"""
import sys

import pandas as pd

import common

disease, kind = sys.argv[1], sys.argv[2]
out = common.out_dir(disease, kind)
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str})

counts = []
for g in common.GROUPS:
    df, _, levels = common.read_presence(out, g)
    counts.append(df.groupby("acc")[levels].sum())
counts = pd.concat(counts).reindex(groups.acc, fill_value=0)  # samples with nothing get 0

level = common.LEVEL[kind] or "level"
df = groups.merge(counts.reset_index().melt("acc", var_name=level, value_name="n_distinct"), on="acc")
if common.LEVEL[kind] is None:
    df = df.drop(columns=level)
df.to_csv(f"{out}/distinct_per_sample.tsv", sep="\t", index=False)
