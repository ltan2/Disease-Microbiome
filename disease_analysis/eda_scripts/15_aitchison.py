"""Func abundance step 4: Aitchison distance between every pair of samples.

Aitchison distance = Euclidean distance between CLR profiles, i.e. the norm of their CLR
difference. Computed as sqrt(|a|^2 + |b|^2 - 2 a.b) so it runs as one matrix product.
Writes results/<disease>/func/aitchison.npy (float32, n x n, order of clr_groups.tsv).
Usage: 15_aitchison.py DISEASE
"""
import sys

import numpy as np
import pandas as pd

import common

out = common.out_dir(sys.argv[1], "func")
X = pd.read_csv(f"{out}/clr.tsv.gz", sep="\t", index_col=0).T.to_numpy()  # samples x KOs
sq = (X ** 2).sum(axis=1)
D = np.sqrt(np.clip(sq[:, None] + sq[None, :] - 2 * X @ X.T, 0, None))  # clip float round-off
np.fill_diagonal(D, 0)
np.save(f"{out}/aitchison.npy", D.astype(np.float32))
print(f"{D.shape} Aitchison distance matrix -> {out}/aitchison.npy")
