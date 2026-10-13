"""Median / mean of each group x group block of a sample x sample matrix.

Within-group blocks use only off-diagonal pairs (each pair once); between-group blocks use all pairs.
Usage: 05_summary.py MATRIX.npy GROUPS.tsv OUT.tsv
"""
import itertools
import sys

import numpy as np
import pandas as pd

import common

matrix, groups, out = sys.argv[1:4]
groups = pd.read_csv(groups, sep="\t").group.to_numpy()
M = np.load(matrix, mmap_mode="r")

rows = []
for a, b in itertools.combinations_with_replacement(common.GROUPS, 2):
    block = M[np.ix_(groups == a, groups == b)]
    v = block[np.triu_indices(len(block), 1)] if a == b else block.ravel()
    rows.append((f"{a} x {b}", v.size, np.median(v), v.mean()))

df = pd.DataFrame(rows, columns=["block", "n_pairs", "median", "mean"])
df.to_csv(out, sep="\t", index=False, float_format="%.4g")
print(df.to_string(index=False))
