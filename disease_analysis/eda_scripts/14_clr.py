"""Func abundance step 3: CLR transform of the KO abundance matrix.

KOs present in fewer than MIN_PREV (fraction) of the samples are dropped. Zeros are then filled by
multiplicative replacement (skbio multi_replace), which keeps each sample a composition, and
clr(x) = log(x) - mean(log(x)) per sample. Samples with none of the kept KOs are dropped.
Writes results/<disease>/func/clr.tsv.gz (KO x sample) and clr_groups.tsv (acc, group) in the
same sample order, which is also the order of aitchison.npy.
Usage: 14_clr.py DISEASE MIN_PREV
"""
import sys

import pandas as pd
from skbio.stats.composition import clr, multi_replace

import common

disease, min_prev = sys.argv[1], float(sys.argv[2])
out = common.out_dir(disease, "func")
ab = pd.read_csv(f"{out}/abundance.tsv.gz", sep="\t", index_col=0)
ab = ab[(ab > 0).mean(axis=1) >= min_prev]
ab = ab.loc[:, ab.sum() > 0]

values = clr(multi_replace(ab.T.to_numpy())).T  # skbio works on samples x features
pd.DataFrame(values, index=ab.index, columns=ab.columns).to_csv(f"{out}/clr.tsv.gz", sep="\t", float_format="%.6g")

groups = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str}).set_index("acc").loc[ab.columns]
groups.reset_index().to_csv(f"{out}/clr_groups.tsv", sep="\t", index=False)
print(f"{ab.shape[0]} KOs with prevalence >= {min_prev}; {ab.shape[1]} samples kept")
print(groups.group.value_counts().reindex(common.GROUPS).to_string())
