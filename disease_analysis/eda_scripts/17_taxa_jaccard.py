"""Taxa step: Jaccard similarity of YACHT taxon presence between samples, one matrix per coverage.

At each coverage, samples with no taxa present are left out (Jaccard is undefined for them).
Writes results/<disease>/taxa/jaccard_<coverage>.npy (float32, n x n) and groups_<coverage>.tsv
(acc, group: the matrix's row order).
Usage: 17_taxa_jaccard.py DISEASE
"""
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import common

out = common.out_dir(sys.argv[1], "taxa")
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t")
pres = pd.concat(common.read_presence(out, g)[0] for g in common.GROUPS)

for cov in common.COVERAGES:
    p = pres[pres[cov] == 1]
    gr = groups[groups.acc.isin(p.acc)].reset_index(drop=True)
    rows = p.acc.map(pd.Series(gr.index, gr.acc)).to_numpy()
    cols = p.taxon.astype("category").cat.codes.to_numpy()
    X = csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(gr), cols.max() + 1))
    inter = (X @ X.T).toarray()
    n = np.asarray(X.sum(axis=1)).ravel()
    np.save(f"{out}/jaccard_{cov}.npy", (inter / (n[:, None] + n[None, :] - inter)).astype(np.float32))
    gr.to_csv(f"{out}/groups_{cov}.tsv", sep="\t", index=False)
    print(f"coverage {cov}: {len(gr)} samples, {X.shape[1]} taxa")
