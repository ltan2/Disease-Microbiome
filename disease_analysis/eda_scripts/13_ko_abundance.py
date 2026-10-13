"""Func abundance step 2: KO x sample relative-abundance matrix (funprofiler, unitigs).

Columns are the samples of func/groups.tsv in group order; a KO missing from a sample is 0.
Writes results/<disease>/func/abundance.tsv.gz.
Usage: 13_ko_abundance.py DISEASE PROCS
"""
import sys
from multiprocessing import Pool

import pandas as pd

import common

if __name__ == "__main__":
    disease, procs = sys.argv[1], int(sys.argv[2])
    out = common.out_dir(disease, "func")
    accs = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str}).acc
    with Pool(procs) as p:
        long = pd.concat(p.map(common.read_ko, accs, chunksize=64))
    ab = long.pivot(index="ko", columns="acc", values="abundance").reindex(columns=accs).fillna(0)
    ab.to_csv(f"{out}/abundance.tsv.gz", sep="\t", float_format="%.6g")
    print(f"{ab.shape[0]} KOs x {ab.shape[1]} samples")
