"""Disease vs control PERMANOVA before and after accounting for BioProject (other_control never shares a
project, so it is left out).

Only samples that finished every step (human-free sketch, KO profile and YACHT result) are used, and
only BioProjects with at least MIN_PER_GROUP disease and MIN_PER_GROUP control of them; for taxa,
samples with no taxa at a coverage are also dropped (Jaccard is undefined), so n can differ per row.
On the same samples, two models:
  before  distance ~ group               labels permuted freely
  after   distance ~ bioproject + group  sequential SS (as vegan adonis2), labels permuted only
                                         within each project (adonis2 strata = bioproject)
Run on the sketch Jaccard, KO Aitchison and taxa Jaccard (per coverage) matrices of the disease.
Writes results/<disease>/stratified_permanova.tsv and results/<disease>/stratified_projects.tsv.
Usage: 20_stratified_permanova.py DISEASE PERMUTATIONS [MIN_PER_GROUP=10]
"""
import os
import sys

import duckdb
import numpy as np
import pandas as pd

import common

disease, perms = sys.argv[1], int(sys.argv[2])
min_n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
out = os.path.dirname(common.out_dir(disease, "func"))

# data -> (matrix, groups, matrix holds similarities)
inputs = {"Sketch (Jaccard)": (f"{out}/sketch/jaccard.npy", f"{out}/sketch/groups.tsv", True),
          "KO (Aitchison)": (f"{out}/func/aitchison.npy", f"{out}/func/clr_groups.tsv", False)}
inputs |= {f"Taxa (Jaccard, cov {c})": (f"{out}/taxa/jaccard_{c}.npy", f"{out}/taxa/groups_{c}.tsv", True)
           for c in common.COVERAGES}

# samples with a sketch, a KO profile (in the Aitchison matrix) and a YACHT result
complete = set.intersection(*(set(pd.read_csv(f, sep="\t").iloc[:, 0]) for f in
                              [f"{out}/sketch/groups.tsv", f"{out}/func/clr_groups.tsv", f"{out}/taxa/groups.tsv"]))

con = duckdb.connect(common.METADATA_DB, read_only=True)
bioproject = con.sql("SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata").df().set_index("acc").bioproject


def stratified_permanova(D, y, strata, perms, batch=100):
    """Pseudo-F of the 0/1 group y after the strata, from the Gower-centred matrix G = -J D^2 J / 2.

    With the strata fitted first, group's SS is r'Gr / r'r, where r is y centred within each stratum;
    permuting y within strata keeps r'r and the strata SS fixed, so only r'Gr changes.
    A single stratum gives the ordinary one-way PERMANOVA.
    """
    n = len(y)
    D2 = D ** 2
    G = -(D2 - D2.mean(axis=0) - D2.mean(axis=1)[:, None] + D2.mean()) / 2
    codes = pd.factorize(strata)[0]
    onehot = np.eye(codes.max() + 1)[codes]
    ss_total = np.trace(G)
    ss_strata = ((onehot.T @ G @ onehot).diagonal() / onehot.sum(axis=0)).sum()
    df_res = n - onehot.shape[1] - 1

    def f_stats(Y):  # Y: n x P group labels
        R = Y - onehot @ ((onehot.T @ Y) / onehot.sum(axis=0)[:, None])
        ss_group = (R * (G @ R)).sum(axis=0) / (R ** 2).sum(axis=0)
        return ss_group / ((ss_total - ss_strata - ss_group) / df_res), ss_group

    F, ss_group = (v[0] for v in f_stats(y[:, None]))
    rng = np.random.default_rng(0)
    idx = [np.flatnonzero(codes == c) for c in range(onehot.shape[1])]

    def shuffle():
        z = y.copy()
        for i in idx:
            z[i] = y[rng.permutation(i)]
        return z

    null = np.concatenate([f_stats(np.column_stack([shuffle() for _ in range(min(batch, perms - i))]))[0]
                           for i in range(0, perms, batch)])
    return F, ss_group / ss_total, ss_strata / ss_total, ((null >= F).sum() + 1) / (perms + 1)


rows, projects = [], []
for data, (matrix, groups, similarity) in inputs.items():
    if not os.path.exists(matrix):
        continue
    g = pd.read_csv(groups, sep="\t", index_col=0).group  # first column: sample id
    g = g.to_frame().assign(bioproject=bioproject.reindex(g.index).to_numpy(), pos=np.arange(len(g)))
    g = g[g.group.isin(["disease", "control"]) & g.index.isin(complete)]
    counts = pd.crosstab(g.bioproject, g.group)
    keep = counts[(counts.disease >= min_n) & (counts.control >= min_n)]
    g = g[g.bioproject.isin(keep.index)]
    projects.append(keep.reset_index().assign(data=data))

    D = np.load(matrix, mmap_mode="r")[np.ix_(g.pos, g.pos)].astype(np.float64)
    if similarity:
        D = 1 - D
    np.fill_diagonal(D, 0)
    y = (g.group == "disease").to_numpy(float)
    F0, r2_0, _, p0 = stratified_permanova(D, y, np.zeros(len(y)), perms)
    F1, r2_1, r2_proj, p1 = stratified_permanova(D, y, g.bioproject.to_numpy(), perms)
    rows.append({"data": data, "n": len(g), "disease": int(y.sum()), "control": int(len(y) - y.sum()),
                 "projects": len(keep),
                 "before F": round(F0, 2), "before R2": round(r2_0, 4), "before p": p0,
                 "after F": round(F1, 2), "after R2 (group)": round(r2_1, 4),
                 "after R2 (bioproject)": round(r2_proj, 4) + 0, "after p": p1})

res = pd.DataFrame(rows)
res.to_csv(f"{out}/stratified_permanova.tsv", sep="\t", index=False)
print(f"{len(complete)} samples finished every step\n" + res.to_string(index=False))
pd.concat(projects)[["data", "bioproject", "disease", "control"]].to_csv(
    f"{out}/stratified_projects.tsv", sep="\t", index=False)
