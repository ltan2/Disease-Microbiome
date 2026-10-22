"""Sequencing depth after human removal, estimated per sample from its sketches.

The non-human fraction is the share of the sample's unitig hashes left after GRCh38 hashes were removed
(01_remove_human.py), so non-human mbases = SRA mbases x (human-free hashes / all hashes). This is a share
of distinct k-mers, not of reads: once a genome is covered, more reads of it add no hashes, so a sample
with deeply covered human DNA has its human read share, and so its human-free depth, misestimated.
Writes results/<disease>/<kind>/depth_nonhuman.tsv (acc, group, mbases = human-free estimate,
mbases_raw, frac_nonhuman); samples without mbases > 0 are left out.
Usage: 37_nonhuman_depth.py DISEASE sketch|taxa|func PROCS
"""
import os
import sys
from multiprocessing import Pool

import pandas as pd
import sourmash

import common


def n_hashes(path):
    if not os.path.exists(path):
        return 0
    return len(next(iter(sourmash.load_file_as_signatures(path))).minhash)


def counts(acc):
    return n_hashes(f"{common.SKETCHES}/{acc}.unitigs.k31.sig.zip"), n_hashes(common.nohuman_path(acc))


if __name__ == "__main__":
    disease, kind, procs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    out = common.out_dir(disease, kind)
    df = pd.read_csv(f"{out}/depth.tsv", sep="\t", dtype={"acc": str})
    df = df[df.mbases > 0].rename(columns={"mbases": "mbases_raw"})
    with Pool(procs) as p:
        n = p.map(counts, df.acc, chunksize=64)
    df["n_all"], df["n_nonhuman"] = zip(*n)
    df = df[df.n_all > 0]  # no original sketch: fraction unknown
    df["frac_nonhuman"] = df.n_nonhuman / df.n_all
    df["mbases"] = df.mbases_raw * df.frac_nonhuman
    df[["acc", "group", "mbases", "mbases_raw", "frac_nonhuman"]].to_csv(
        f"{out}/depth_nonhuman.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"{len(df)} samples; median non-human fraction by group:")
    print(df.groupby("group").frac_nonhuman.median().round(3).to_string())
