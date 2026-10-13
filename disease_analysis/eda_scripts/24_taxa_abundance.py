"""Per-sample taxa relative abundance proxies from YACHT, for one disease's three groups.

Containment-based: a_i^(C) = C_i / sum_j C_j, with C_i the containment of taxon i (fraction of
its exclusive k-mers detected; see 23_yacht_containment.py).
Detection-score-based: with the coverage thresholds ordered highest to lowest (k = 1..K, i.e.
common.COVERAGES) and weights w_k, k_i* is the highest threshold taxon i passes; detection there
implies detection at every lower threshold, so
    s_i = sum_{k=k_i*}^{K} w_k / sum_{k=1}^{K} w_k,   a_i^(S) = s_i / sum_j s_j.
Both are normalized over the taxa detected in a sample (present at any common.COVERAGES level).

Writes results/<disease>/taxa/abundance.tsv.gz: acc, group, taxon, containment,
containment_abundance, detection_score, detection_abundance. Samples with no taxa have no rows.
Usage: 24_taxa_abundance.py DISEASE [WEIGHTS]   (WEIGHTS: comma-separated w_1..w_K, default all 1)
"""
import sys

import numpy as np
import pandas as pd

import common

disease = sys.argv[1]
w = np.array(sys.argv[2].split(","), dtype=float) if len(sys.argv) > 2 else np.ones(len(common.COVERAGES))
if len(w) != len(common.COVERAGES) or (w < 0).any() or w.sum() == 0:
    raise SystemExit(f"need {len(common.COVERAGES)} non-negative weights, one per coverage {common.COVERAGES}")

out = common.out_dir(disease, "taxa")
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str})

# detection score: tail[k] = sum_{k'>=k} w_k' / sum w; k* = first (highest) coverage passed
pres = pd.concat(common.read_presence(out, g)[0] for g in common.GROUPS)
passed = pres[common.COVERAGES].to_numpy().astype(bool)
tail = w[::-1].cumsum()[::-1] / w.sum()
pres["detection_score"] = tail[passed.argmax(axis=1)]
pres = pres[passed.any(axis=1)].drop(columns=common.COVERAGES)

cont = pd.read_csv(common.YACHT_CONTAINMENT, sep="\t", dtype={"acc": str},
                   usecols=["acc", "taxon", "containment"])
cont = cont[cont.acc.isin(groups.acc)]  # merging on all samples' rows is very slow
df =groups.merge(pres, on="acc").merge(cont, on=["acc", "taxon"], how="left")
if df.containment.isna().any():
    raise SystemExit(f"{df.containment.isna().sum()} detected (sample, taxon) pairs lack a containment; "
                     f"re-run 23_yacht_containment.py")

for score, abund in [("containment", "containment_abundance"), ("detection_score", "detection_abundance")]:
    df[abund] = df[score] / df.groupby("acc")[score].transform("sum")
cols = ["acc", "group", "taxon", "containment", "containment_abundance", "detection_score", "detection_abundance"]
df[cols].to_csv(f"{out}/abundance.tsv.gz", sep="\t", index=False)
print(df.groupby("group").acc.nunique().reindex(common.GROUPS).to_string())
