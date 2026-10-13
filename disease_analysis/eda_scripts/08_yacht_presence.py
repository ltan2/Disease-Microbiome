"""Taxa step 0: YACHT presence/absence of every taxon in every sample at each coverage level.

Each result_*_ANI_<acc>_unitigs.xlsx has one sheet per min_coverage listing the taxa YACHT calls
present. Writes results/yacht_presence.tsv.gz (acc, taxon, one 0/1 column per common.COVERAGES),
keeping taxa present at any of those levels, and results/yacht_samples.txt with every sample that
has a result file (empty file = YACHT found nothing, so the sample has zero taxa).
Usage: 08_yacht_presence.py PROCS
"""
import glob
import os
import re
import sys
from multiprocessing import Pool

import pandas as pd

import common


def read(path):
    acc = re.search(r"_ANI_(.+)_unitigs\.xlsx$", path).group(1)
    if os.path.getsize(path) == 0:
        return acc, None
    sheets = pd.read_excel(path, sheet_name=[f"min_coverage{c}" for c in common.COVERAGES],
                           usecols=["organism_name", "in_sample_est"])
    present = {c: set(s.organism_name[s.in_sample_est]) for c, s in zip(common.COVERAGES, sheets.values())}
    taxa = sorted(set.union(*present.values()))
    df = pd.DataFrame({c: [int(t in present[c]) for t in taxa] for c in common.COVERAGES})
    df.insert(0, "taxon", [t.removesuffix("_genomic") for t in taxa])
    df.insert(0, "acc", acc)
    return acc, df


if __name__ == "__main__":
    with Pool(int(sys.argv[1])) as p:
        res = p.map(read, glob.glob(f"{common.YACHT}/result_*_ANI_*_unitigs.xlsx"), chunksize=16)
    pd.concat([df for _, df in res if df is not None]).to_csv(common.YACHT_PRESENCE, sep="\t", index=False)
    with open(common.YACHT_SAMPLES, "w") as f:
        f.write("\n".join(sorted(acc for acc, _ in res)) + "\n")
    print(f"{len(res)} YACHT samples -> {common.YACHT_PRESENCE}")
