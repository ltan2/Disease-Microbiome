"""A per-sample count divided by the sample's sequencing depth (SRA mbases, before human removal).

Samples without a positive mbases value are left out.
Writes OUT.tsv: the input columns, mbases and COLUMN_per_mbp; plotted by plot_group_density.py.
Usage: 36_per_mbp.py COUNTS.tsv DEPTH.tsv COLUMN OUT.tsv
"""
import sys

import pandas as pd

counts, depth, column, out = sys.argv[1:]
df = pd.read_csv(counts, sep="\t", dtype={"acc": str})
dep = pd.read_csv(depth, sep="\t", dtype={"acc": str})[["acc", "mbases"]]
df = df.merge(dep[dep.mbases > 0], on="acc")
df[f"{column}_per_mbp"] = df[column] / df.mbases
df.to_csv(out, sep="\t", index=False)
print(f"{df.acc.nunique()} samples with mbases > 0")
