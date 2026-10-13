"""Taxa abundance step 0: YACHT containment of every detected taxon in every sample.

containment = num_matches / num_exclusive_kmers_to_genome, the fraction of a taxon's exclusive
k-mers found in the sample. num_matches does not depend on min_coverage and the sheets are nested,
so the lowest common.COVERAGES sheet holds every taxon kept in results/yacht_presence.tsv.gz.
Writes results/yacht_containment.tsv.gz (acc, taxon, num_matches, num_exclusive_kmers_to_genome,
containment).
Usage: 23_yacht_containment.py PROCS
"""
import glob
import os
import re
import sys
from multiprocessing import Pool

import pandas as pd

import common

COLS = ["organism_name", "in_sample_est", "num_matches", "num_exclusive_kmers_to_genome"]


def read(path):
    if os.path.getsize(path) == 0:
        return None
    df = pd.read_excel(path, sheet_name=f"min_coverage{common.COVERAGES[-1]}", usecols=COLS)
    df = df.loc[df.in_sample_est.astype(bool)].drop(columns="in_sample_est")  # empty sheet: object dtype
    df.insert(0, "taxon", df.pop("organism_name").str.removesuffix("_genomic"))
    df.insert(0, "acc", re.search(r"_ANI_(.+)_unitigs\.xlsx$", path).group(1))
    df["containment"] = df.num_matches / df.num_exclusive_kmers_to_genome
    return df


if __name__ == "__main__":
    with Pool(int(sys.argv[1])) as p:
        res = p.map(read, glob.glob(f"{common.YACHT}/result_*_ANI_*_unitigs.xlsx"), chunksize=16)
    df = pd.concat([d for d in res if d is not None])
    df.to_csv(common.YACHT_CONTAINMENT, sep="\t", index=False)
    print(f"{df.acc.nunique()} samples, {len(df)} (sample, taxon) rows -> {common.YACHT_CONTAINMENT}")
