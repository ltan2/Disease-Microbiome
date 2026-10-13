"""Step 3: all-vs-all Jaccard similarity of the human-free sketches in groups.tsv order.

Jaccard = |A & B| / |A | B|, with intersections from a sparse sample x hash matrix product.
Writes results/<disease>/sketch/jaccard.npy (float32, n x n).
Usage: 03_jaccard.py DISEASE PROCS
"""
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import common

BLOCK = 1000  # rows per worker task


def block(i):
    inter = (X[i:i + BLOCK] @ X.T).toarray()
    return inter / (n[i:i + BLOCK, None] + n[None, :] - inter)


if __name__ == "__main__":
    disease, procs = sys.argv[1], int(sys.argv[2])
    out = common.out_dir(disease, "sketch")
    accs = pd.read_csv(f"{out}/groups.tsv", sep="\t").acc
    hashes = common.load_nohuman(accs, procs)

    # binary sample x hash matrix; columns are the distinct hashes across all samples
    n = np.array([len(h) for h in hashes], dtype=np.float32)
    _, cols = np.unique(np.concatenate(hashes), return_inverse=True)
    rows = np.repeat(np.arange(len(hashes)), n.astype(int))
    X = csr_matrix((np.ones(len(cols), np.float32), (rows, cols)))

    with Pool(procs) as p:  # workers inherit X and n through fork
        J = np.vstack(p.map(block, range(0, len(accs), BLOCK))).astype(np.float32)
    np.save(f"{out}/jaccard.npy", J)
    print(f"{J.shape} Jaccard matrix -> {out}/jaccard.npy")
