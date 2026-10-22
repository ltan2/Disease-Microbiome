"""Per-taxon association of YACHT relative abundance with the groups, for forest plots.

Taxa are aggregated to GTDB r232 species, genus and family (sum of the genomes' abundances per
sample). For each clade and each pair of groups (CONTRASTS), y = log2(abundance + pc) per sample,
with 0 for samples where the clade was not detected and pc = half the smallest nonzero abundance;
effect = mean_A(y) - mean_B(y), a log2 fold change of geometric means (> 0: higher in A), with a
Welch t 95% CI and p, and BH q within each (rank, contrast). Samples with no taxa at all are left
out. Clades are tested when present in at least MIN_PREV of the samples of some group.

Literature markers (literature_markers.tsv rows for DISEASE: clade with its GTDB rank prefix,
expected direction in disease up/down/variable) are tested the same way; a genus or species name also
covers GTDB's suffixed splits of it (g__Roseburia_A, s__Faecalibacterium prausnitzii_I). Each
marker is followed by its MEMBERS most prevalent species (unless its only species is itself). A
marker detected in no sample keeps an empty 'not detected' row. verdict (disease vs each control
group): consistent / opposite (q < 0.05) or n.s.; for 'variable': enriched / depleted / n.s.

Writes results/<disease>/taxa/assoc_<COLUMN>.tsv.gz (all ranks) and markers_<COLUMN>.tsv, both
with label, contrast, effect, ci_lo, ci_hi, p, q, prevalence and mean abundance per group (the
format read by plot_forest.py).
Usage: 25_taxa_assoc.py DISEASE COLUMN MIN_PREV [MEMBERS]   (COLUMN: an abundance.tsv.gz column)
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

import common

CONTRASTS = [("disease", "control"), ("disease", "other_control"), ("control", "other_control")]
RANKS = ["species", "genus", "family"]
ALPHA = 0.05

disease, column, min_prev = sys.argv[1], sys.argv[2], float(sys.argv[3])
members = int(sys.argv[4]) if len(sys.argv) > 4 else 3
out = common.out_dir(disease, "taxa")

ab = pd.read_csv(f"{out}/abundance.tsv.gz", sep="\t", usecols=["acc", "group", "taxon", column])
ab = ab.rename(columns={column: "a"}).merge(common.read_gtdb(), on="taxon", how="left")
if ab.species.isna().any():
    raise SystemExit(f"{ab[ab.species.isna()].taxon.nunique()} taxa missing from {common.GTDB_TAXONOMY}")
n = ab.groupby("group").acc.nunique().reindex(common.GROUPS)  # samples with any taxon
y0 = np.log2(pc := ab.a[ab.a > 0].min() / 2)


def assoc(obs):
    """obs: acc, group, label, a (genomes of one clade per label). One row per (label, contrast)."""
    obs = obs.groupby(["label", "group", "acc"]).a.sum().reset_index()  # clade abundance per sample
    obs["y"] = np.log2(obs.a + pc)
    obs["y2"] = obs.y ** 2
    s = obs.groupby(["label", "group"]).agg(n_det=("y", "size"), s1=("y", "sum"), s2=("y2", "sum"), sa=("a", "sum"))
    s = s.unstack("group", fill_value=0)  # a group with no detection gets zero columns, not a missing level
    s = s.reindex(columns=pd.MultiIndex.from_product([s.columns.levels[0], common.GROUPS]), fill_value=0)
    g = {}
    for grp in common.GROUPS:
        N, z = n[grp], n[grp] - s[("n_det", grp)]
        mean = (s[("s1", grp)] + z * y0) / N
        var = ((s[("s2", grp)] + z * y0 ** 2 - N * mean ** 2) / (N - 1)).clip(lower=0)
        g[grp] = (N, mean, var)
    base = pd.DataFrame({f"prev_{grp}": s[("n_det", grp)] / n[grp] for grp in common.GROUPS})
    for grp in common.GROUPS:
        base[f"mean_abundance_{grp}"] = s[("sa", grp)] / n[grp]
    rows = []
    for A, B in CONTRASTS:
        (na, ma, va), (nb, mb, vb) = g[A], g[B]
        se = np.sqrt(va / na + vb / nb).replace(0, np.nan)
        df = se ** 4 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
        r = base.copy()
        r.insert(0, "contrast", f"{A} vs {B}")
        r.insert(1, "effect", ma - mb)
        half = stats.t.ppf(0.975, df) * se
        r.insert(2, "ci_lo", r.effect - half)
        r.insert(3, "ci_hi", r.effect + half)
        r.insert(4, "p", 2 * stats.t.sf(np.abs(r.effect / se), df))
        rows.append(r)
    return pd.concat(rows).reset_index()


def bh(df, by):
    df["q"] = np.nan
    for _, idx in df.groupby(by).groups.items():
        p = df.loc[idx, "p"].dropna()
        df.loc[p.index, "q"] = stats.false_discovery_control(p)
    return df


# ---- every clade at each rank ----
res = []
for rank in RANKS:
    r = assoc(ab[["acc", "group", "a"]].assign(label=ab[rank]))
    r.insert(0, "rank", rank)
    res.append(r[r.filter(like="prev_").max(axis=1) >= min_prev])
res = bh(pd.concat(res, ignore_index=True), ["rank", "contrast"])
res = res[["rank", "label", "contrast", "effect", "ci_lo", "ci_hi", "p", "q"]
          + [c for c in res if c.startswith(("prev_", "mean_abundance_"))]]
res.to_csv(f"{out}/assoc_{column}.tsv.gz", sep="\t", index=False)

# ---- literature markers and their most prevalent species ----
lit = pd.read_csv("literature_markers.tsv", sep="\t")
lit = lit[lit.disease == disease]
mk = []
for clade, expected in zip(lit.clade, lit.expected):
    rank = common.RANK_PREFIX[clade[:3]]
    in_clade = (ab[rank] == clade) | ab[rank].str.startswith(clade + "_") if rank in ("genus", "species") \
        else ab[rank] == clade
    arrow = {"up": "expected ↑", "down": "expected ↓"}.get(expected, expected)
    if not in_clade.any():
        print(f"{clade}: not detected in any sample")
        mk.append(pd.DataFrame({"label": f"{clade} ({arrow}; not detected)", "contrast": [f"{a} vs {b}" for a, b in CONTRASTS],
                                "rank": rank, "marker": clade, "expected": expected, "kind": "clade"}))
        continue
    sub = ab[in_clade]
    r = assoc(sub[["acc", "group", "a"]].assign(label=f"{clade} ({arrow})"))
    r["rank"], r["marker"], r["expected"], r["kind"] = rank, clade, expected, "clade"
    if set(sub.species) == {clade}:  # a species without GTDB splits: no member rows
        mk.append(r)
        continue
    sp = res[(res["rank"] == "species") & res.label.isin(set(sub.species)) & (res.contrast == "disease vs control")]
    top = sp.assign(m=sp.filter(like="prev_").max(axis=1)).nlargest(members, "m").label
    m = res[(res["rank"] == "species") & res.label.isin(top)].assign(marker=clade, expected=expected, kind="member")
    m["label"] = "    " + m.label  # indented under the marker in the plot
    m = m.set_index("label").loc[[f"    {t}" for t in top]].reset_index()
    mk += [r, m]
cols = ["marker", "expected", "kind", "rank", "label", "contrast", "effect", "ci_lo", "ci_hi", "p", "q"]
if mk:
    mk = pd.concat(mk, ignore_index=True)
    # marker clades get their own BH within each contrast; member species keep their q from assoc
    clades = mk.kind == "clade"
    mk.loc[clades, "q"] = bh(mk[clades].copy(), "contrast").q
    sign = np.where(mk.expected == "up", 1, -1) * np.sign(mk.effect)
    mk["verdict"] = np.where(mk.q.isna(), "", np.where(mk.q >= ALPHA, "n.s.", np.where(sign > 0, "consistent", "opposite")))
    var = (mk.expected == "variable") & (mk.q < ALPHA)
    mk.loc[var, "verdict"] = np.where(mk.effect[var] > 0, "enriched", "depleted")
    mk.loc[mk.contrast == "control vs other_control", "verdict"] = ""
    mk = mk[cols + ["verdict"] + [c for c in mk if c.startswith(("prev_", "mean_abundance_"))]]
else:
    mk = pd.DataFrame(columns=cols + ["verdict"])
mk.to_csv(f"{out}/markers_{column}.tsv", sep="\t", index=False)

print(f"pseudocount {pc:.3g}; samples with taxa: {n.to_dict()}")
print(res.groupby(["rank", "contrast"]).q.apply(lambda q: f"{(q < ALPHA).sum()}/{len(q)} q<{ALPHA}").unstack().to_string())
if len(mk):
    show = mk[mk.kind == "clade"]
    print(show[["label", "contrast", "effect", "q", "verdict"]].round(4).to_string(index=False))
