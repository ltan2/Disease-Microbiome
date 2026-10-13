"""Func step 2: KO presence/absence per group from the funprofiler KO profiles (unitigs).

A KO is present in a sample if it is listed in <acc>.unitigs_ko_profiles.csv (header only = no KOs).
Writes results/<disease>/func/presence_<group>.tsv.gz (acc, ko) and matrix_<group>.tsv.gz,
the KO x sample 0/1 matrix.
Usage: 12_ko_presence.py DISEASE PROCS
"""
import sys
from multiprocessing import Pool

import pandas as pd

import common


if __name__ == "__main__":
    disease, procs = sys.argv[1], int(sys.argv[2])
    out = common.out_dir(disease, "func")
    groups = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str})
    with Pool(procs) as p:
        pres = pd.concat(p.map(common.read_ko, groups.acc, chunksize=64)).merge(groups, on="acc")

    for g in common.GROUPS:
        sub = pres[pres.group == g]
        sub[["acc", "ko"]].to_csv(f"{out}/presence_{g}.tsv.gz", sep="\t", index=False)
        accs = groups.acc[groups.group == g]
        matrix = pd.crosstab(sub.ko, sub.acc).reindex(columns=accs, fill_value=0).clip(upper=1)
        matrix.to_csv(f"{out}/matrix_{g}.tsv.gz", sep="\t")
    print(groups.group.value_counts().reindex(common.GROUPS).to_string())
