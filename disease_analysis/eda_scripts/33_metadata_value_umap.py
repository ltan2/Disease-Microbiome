"""Group the CRC / Crohn's taxa-analysis samples by their metadata values for chosen key clusters of 32_metadata_key_umap.py.

A key group is the union of the key clusters that hold the given key labels. The keys of
results/metadata_keys/keys.tsv are re-embedded and re-clustered (average linkage on 1 - cosine, as in 32) at the
group's key cut: 0.5 reproduces 32's clusters; a coarser cut merges neighbours, e.g. diet at 0.7 = the food / diet
type / food-frequency / caffeine super-cluster plus the vitamin-carotenoid and the meat super-clusters.
environment = the 0.5 clusters of "env biome", "env material", "env broad scale", "env feature". For each group a
sample keeps its non-missing values of the group's keys and is dropped if it has none. Sample x sample similarity:

  text   (text values only: environment) each "key: value" is embedded with BioLORD-2023, a sample is the mean of
         its embeddings, similarity = cosine. The key is part of the text so "no" for "exclude dairy" and "no" for
         "kosher" differ, while different studies' names for one field (env biome / environment biome) match.
  gower  (numbers present: ulcer scores, diet nutrient amounts and frequency answers) per key both samples have,
         |a - b| / key range for two numbers, 1 - BioLORD cosine of the values for two texts (Daily vs Regularly),
         1 for a number vs a text (e.g. Healthy_control); similarity = 1 - mean over shared keys, 0 when none shared.

--mode auto takes gower when at least 10% of a group's values are numbers. Samples are clustered (average linkage, cut at
--cut in 1 - similarity) and projected with UMAP on 1 - similarity.

Writes results/metadata_values/<group>/:
  samples.tsv     one row per kept sample: diseases, status, its "key: value" text, sample cluster, UMAP x/y
  clusters.tsv    per sample cluster: size, samples per disease x status, most common values
  similarity.npy  sample x sample similarity (rows in samples.tsv order)
  umap.png        UMAP colored by disease x status and by sample cluster
Usage: 33_metadata_value_umap.py [DISEASE ...] [--group NAME=LABEL[;LABEL...][@KEYCUT] ...] [--mode auto|text|gower] [--cut 0.2]
"""
import argparse
import os
from collections import Counter

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import umap
from sentence_transformers import SentenceTransformer
from sklearn.cluster import AgglomerativeClustering

import common

# categorical slots 1-8 (blue, orange, aqua, yellow, magenta, green, violet, red); other clusters in grey
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#b4b2a9"
# group: (key labels whose clusters form it, key-tree cut)
GROUPS = {
    "environment": (["env biome", "env material", "env broad scale", "env feature"], 0.5),
    "ulcers": (["presence and size of ulcers ses cd rectum"], 0.5),
    "diet": (["diet", "beta carotene provitamin a carotenoid", "beef"], 0.7),
}
MISSING = {"", "na", "n/a", "nan", "none", "missing", "not applicable", "not collected", "not provided",
           "unknown", "not available", "uncalculated", "-", "missing: not provided", "null"}
BOOL = {"true": "yes", "y": "yes", "yes": "yes", "false": "no", "n": "no", "no": "no"}

p = argparse.ArgumentParser()
p.add_argument("diseases", nargs="*", default=["Colorectal cancer", "Crohn's disease"])
p.add_argument("--groups-from", default="disease,control", help="taxa-analysis groups to take samples from")
p.add_argument("--group", action="append", default=[],
               help="NAME=LABEL[;LABEL...][@KEYCUT]: key group made of the key clusters (key tree cut at KEYCUT, "
                    "default 0.5) holding these key labels (default: environment, ulcers, diet)")
p.add_argument("--mode", default="auto", choices=["auto", "text", "gower"])
p.add_argument("--cut", type=float, default=0.2, help="1 - similarity at which the sample dendrogram is cut")
p.add_argument("--neighbors", type=int, default=15, help="UMAP n_neighbors")
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
groups = {}
for g in args.group:
    name, spec = g.split("=", 1)
    labels, _, cut = spec.partition("@")
    groups[name] = (labels.split(";"), float(cut or 0.5))
groups = groups or GROUPS

keys = pd.read_csv(f"{common.RESULTS}/metadata_keys/keys.tsv", sep="\t")
label_of = {r: l for l, rs in zip(keys.label, keys.raw_keys) for r in rs.split(",")}

# ---- samples and their disease x status ----
member = []
for d in args.diseases:
    ab = pd.read_csv(f"{common.out_dir(d, 'taxa', site='')}/abundance.tsv.gz", sep="\t", usecols=["acc", "group"])
    ab = ab[ab.group.isin(args.groups_from.split(","))].drop_duplicates("acc")
    member += [(a, d, g) for a, g in zip(ab.acc, ab.group)]
member = pd.DataFrame(member, columns=["acc", "disease", "status"])
info = member.groupby("acc").agg(diseases=("disease", " + ".join), status=("status", " + ".join))
# one color class per sample: its first disease x status
info["set"] = member.drop_duplicates("acc").set_index("acc").apply(lambda r: f"{r.disease}: {r.status}", axis=1)
sets = [f"{d}: {s}" for d in args.diseases for s in args.groups_from.split(",")]

con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("s", pd.DataFrame({"acc": info.index}))
model = SentenceTransformer("FremyCompany/BioLORD-2023", device="cuda" if torch.cuda.is_available() else "cpu")
key_emb = model.encode(keys.label.tolist(), batch_size=256, normalize_embeddings=True)
key_dist = np.clip(1 - key_emb @ key_emb.T, 0, None)
key_dist = (key_dist + key_dist.T) / 2
np.fill_diagonal(key_dist, 0)


def key_clusters(cut):
    return AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="average",
                                   distance_threshold=cut).fit_predict(key_dist)


def number(v):
    try:
        return float(v)
    except ValueError:
        return np.nan


def gower(wide):
    """1 - mean over shared keys of |a - b| / range (numbers), 1 - BioLORD cosine (texts) or 1 (number vs text);
    0 where nothing is shared."""
    n = len(wide)
    total, shared = np.zeros((n, n)), np.zeros((n, n))
    for k in wide:
        v = wide[k].to_numpy(dtype=object)
        have = np.array([x is not None and x == x for x in v])
        x = np.array([number(a) if h else np.nan for a, h in zip(v, have)])
        num = ~np.isnan(x)
        rng = np.nanmax(x) - np.nanmin(x) if num.any() else 0
        d = np.abs(x[:, None] - x[None, :]) / (rng or 1)
        both_num = num[:, None] & num[None, :]
        text = [str(a) for a in v]
        uniq = sorted(set(t for t, h, nm in zip(text, have, num) if h and not nm))
        if uniq:
            e = model.encode(uniq, batch_size=256, normalize_embeddings=True)
            row = {t: i for i, t in enumerate(uniq)}
            idx = np.array([row.get(t, 0) for t in text])
            sem = np.clip(1 - e[idx] @ e[idx].T, 0, 1)
        else:
            sem = np.ones((n, n))
        both_text = (have & ~num)[:, None] & (have & ~num)[None, :]
        d = np.where(both_num, d, np.where(both_text, sem, 1.0))
        both = have[:, None] & have[None, :]
        total += np.where(both, d, 0)
        shared += both
    return np.where(shared > 0, 1 - total / np.maximum(shared, 1), 0.0)


def plot(df, out, name, mode):
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    for s, c in zip(sets, PALETTE):
        m = df.set == s
        axes[0].scatter(df.umap_x[m], df.umap_y[m], s=14, c=c, edgecolors="white", linewidths=0.3,
                        label=f"{s} ({m.sum()})")
    axes[0].legend(loc="best", fontsize=9, frameon=False)
    axes[0].set_title("disease x status")
    big = df.cluster < len(PALETTE)
    axes[1].scatter(df.umap_x[~big], df.umap_y[~big], s=14, c=GREY, edgecolors="white", linewidths=0.3,
                    label="other clusters")
    for c in range(min(len(PALETTE), df.cluster.max() + 1)):
        m = df.cluster == c
        axes[1].scatter(df.umap_x[m], df.umap_y[m], s=14, c=PALETTE[c], edgecolors="white", linewidths=0.3,
                        label=f"{c}: {Counter(df.text[m]).most_common(1)[0][0][:70]} ({m.sum()})")
        axes[1].text(df.umap_x[m].median(), df.umap_y[m].median(), str(c), fontsize=12, fontweight="bold",
                     color="#2b2b2b")
    axes[1].legend(loc="best", fontsize=8, frameon=False, title="sample cluster: most common values")
    axes[1].set_title("sample clusters")
    for ax in axes:
        ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle(f"{name}: {len(df)} samples by metadata values ({mode} similarity)")
    fig.tight_layout()
    fig.savefig(f"{out}/umap.png", dpi=200)
    plt.close(fig)


for name, (labels, key_cut) in groups.items():
    out = f"{common.RESULTS}/metadata_values/{name}"
    os.makedirs(out, exist_ok=True)
    missing = [l for l in labels if l not in set(keys.label)]
    if missing:
        raise SystemExit(f"{name}: key labels not in keys.tsv: {missing}")
    kc = key_clusters(key_cut)
    cl = set(kc[keys.label.isin(labels)])
    raw = [r for rs in keys.raw_keys[np.isin(kc, list(cl))] for r in rs.split(",")]

    # ---- non-missing values of the group's keys ----
    con.register("k", pd.DataFrame({"k": raw}))
    val = con.sql("SELECT acc, k, any_value(v) AS v FROM attributes JOIN s USING (acc) JOIN k USING (k) "
                  "GROUP BY acc, k").df()
    val.v = val.v.str.strip()
    val = val[~val.v.str.lower().isin(MISSING)]
    val["v"] = val.v.map(lambda v: BOOL.get(v.lower(), v))
    val["label"] = val.k.map(label_of)
    val = val.drop_duplicates(["acc", "label"])
    mode = args.mode
    if mode == "auto":
        mode = "gower" if val.v.map(number).notna().mean() >= 0.1 else "text"
    val["pair"] = val.label + ": " + val.v
    df = info.loc[sorted(set(val.acc))].copy()
    df["text"] = val.sort_values("label").groupby("acc").pair.agg(" | ".join)
    print(f"{name}: {len(cl)} key clusters at key cut {key_cut}, {len(raw)} raw keys, {len(df)} of {len(info)} samples with a value, "
          f"{mode} similarity")

    # ---- sample similarity ----
    if mode == "text":
        pairs = sorted(set(val.pair))
        emb = model.encode(pairs, batch_size=256, normalize_embeddings=True)
        row = pd.Series(range(len(pairs)), index=pairs)
        vec = np.stack([emb[row[g.pair]].mean(0) for _, g in val.set_index("acc").loc[df.index].groupby(level=0, sort=False)])
        vec /= np.linalg.norm(vec, axis=1, keepdims=True)
        sim = vec @ vec.T
    else:
        sim = gower(val.pivot(index="acc", columns="label", values="v").loc[df.index])
    np.save(f"{out}/similarity.npy", sim)
    dist = np.clip(1 - sim, 0, None)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)

    c = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="average",
                                distance_threshold=args.cut).fit_predict(dist)
    rank = pd.Series(c).value_counts().index
    df["cluster"] = pd.Series(c).map({x: i for i, x in enumerate(rank)}).values
    xy = umap.UMAP(n_neighbors=min(args.neighbors, len(df) - 1), metric="precomputed",
                   random_state=args.seed).fit_transform(dist)
    df["umap_x"], df["umap_y"] = xy[:, 0], xy[:, 1]
    df.sort_values("cluster").to_csv(f"{out}/samples.tsv", sep="\t")

    summary = pd.crosstab(df.cluster, df.set).reindex(columns=[s for s in sets if s in set(df.set)], fill_value=0)
    summary.insert(0, "size", summary.sum(axis=1))
    summary["top_values"] = df.groupby("cluster").text.agg(
        lambda t: " || ".join(f"{v} ({n})" for v, n in Counter(t).most_common(3)))
    summary.to_csv(f"{out}/clusters.tsv", sep="\t")
    print(f"  {len(summary)} sample clusters at cut {args.cut}; largest:")
    print(summary.head(8).to_string(max_colwidth=110))
    plot(df, out, name, mode)
print(f"wrote {common.RESULTS}/metadata_values/")
