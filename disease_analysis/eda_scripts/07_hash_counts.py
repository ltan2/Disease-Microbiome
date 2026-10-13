"""Step 7: number of unique k-mer hashes per sample, after human hashes were removed.

Writes results/<disease>/sketch/hash_counts.tsv (acc, group, n_unique_hashes); plotted by plot_group_density.py.
Usage: 07_hash_counts.py DISEASE PROCS
"""
import sys

import pandas as pd

import common

if __name__ == "__main__":
    disease, procs = sys.argv[1], int(sys.argv[2])
    out = common.out_dir(disease, "sketch")
    df = pd.read_csv(f"{out}/groups.tsv", sep="\t")
    df["n_unique_hashes"] = [len(h) for h in common.load_nohuman(df.acc, procs)]
    df.to_csv(f"{out}/hash_counts.tsv", sep="\t", index=False)
