"""Step 1: remove GRCh38 hashes from every unitig sketch in all-res.

Writes results/nohuman_unitigs/ACC.sig.zip (abundances kept). Existing outputs are skipped, and
samples left with no hashes are not written.
Usage: 01_remove_human.py PROCS
"""
import glob
import os
import sys
from multiprocessing import Pool

import numpy as np
import sourmash
from sourmash.save_load import SaveSignaturesToLocation

import common

HUMAN = common.load_hashes(common.HUMAN)  # loaded before the fork, shared by all workers


def run(path):
    acc = os.path.basename(path).split(".")[0]
    out = common.nohuman_path(acc)
    if os.path.exists(out):
        return
    mh = next(iter(sourmash.load_file_as_signatures(path))).minhash
    abund = mh.hashes
    h = np.fromiter(abund.keys(), dtype=np.uint64)
    is_human = HUMAN[np.searchsorted(HUMAN, h).clip(max=len(HUMAN) - 1)] == h
    if is_human.all():
        return
    # rebuild from the kept hashes; MinHash.remove_many is quadratic on large sketches
    mh = mh.copy_and_clear()
    mh.set_abundances({int(x): abund[int(x)] for x in h[~is_human]})
    tmp = f"{common.NOHUMAN}/.{acc}.tmp.sig.zip"  # rename at the end so partial files never look done
    with SaveSignaturesToLocation(tmp) as save:
        save.add(sourmash.SourmashSignature(mh, name=acc))
    os.replace(tmp, out)


if __name__ == "__main__":
    os.makedirs(common.NOHUMAN, exist_ok=True)
    paths = glob.glob(f"{common.SKETCHES}/*.unitigs.k31.sig.zip")
    with Pool(int(sys.argv[1])) as p:
        p.map(run, paths, chunksize=64)
    print(f"{len(os.listdir(common.NOHUMAN))} human-free sketches in {common.NOHUMAN}")
