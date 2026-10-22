"""Group the SRA/BioSample metadata keys of the taxa-analysis samples by meaning (BioLORD-2023) and map them with UMAP.

Samples: the disease and control samples of each disease's taxa analysis (results/<disease>/taxa/abundance.tsv.gz,
unfiltered by site). Keys: every key of each sample's sra_metadata.jattr JSON; when a value is a dictionary (also
inside a list) its keys count as metadata keys too, recursively. Keys that differ only by the SRA suffixes
(_sam, _run, _calc, _s_dpl<N>, ...) are merged and embedded as plain words ("env_broad_scale_sam" -> "env broad scale").

The cleaned keys are embedded with BioLORD-2023; cosine similarity between them gives the distance (1 - similarity)
for UMAP and for average-linkage clustering cut at --cut. Keys whose value differs in every sample having them
(sample ids, timestamps; only judged for keys in >= --unique-min samples) are left out, unless every value is a
number (a measurement: age, BMI, nutrient intake) and the key's name is not an id / label / name.

Clusters are ranked by closeness to the disease label: the anchor is the sample-weighted mean embedding of the keys
the samples' disease labels were matched on (selection.matched_keys, disease_health keys only), and a cluster's
score is the cosine similarity of its mean embedding to that anchor. Every point is labeled; the --top clusters
closest to the disease label are colored.

Writes results/metadata_keys/:
  keys.tsv            one row per cleaned key: label, raw keys, samples having it (all and per disease), distinct
                      values, cluster, UMAP x/y
  unique_keys.tsv     keys left out because every sample has a different non-numeric value or it is a numeric id
  clusters.tsv        per cluster: size, samples, disease-label similarity and rank, label keys it holds, member
                      labels (largest first)
  disease_clusters.tsv  the --top clusters closest to the disease label, one row per member key
  similarity.npy      key x key cosine similarity (rows in keys.tsv order)
  umap.png umap.pdf   labeled UMAP (zoom into the PDF for the labels)
Usage: 32_metadata_key_umap.py [DISEASE ...] [--groups disease,control] [--min-samples 1] [--unique-min 2] [--cut 0.5] [--top 10]
"""
import argparse
import json
import os
import re
from collections import defaultdict

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import umap
from adjustText import adjust_text
from sentence_transformers import SentenceTransformer
from sklearn.cluster import AgglomerativeClustering

import common

# categorical slots 1-8 (blue, orange, aqua, yellow, magenta, green, violet, red); other clusters in grey
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#b4b2a9"
# highlighted clusters by rank: 5 hues, then the same hues with the next marker (ranks 1-5 o, 6-10 ^, 11-15 s)
N_HUE = 5
MARKERS = ["o", "^", "s"]


def style(r):
    return PALETTE[(r - 1) % N_HUE], MARKERS[(r - 1) // N_HUE]

SUFFIX = re.compile(r"(_(sam|run|calc|exp))?(_s+_dpl\d+)?$")
# words that make a numeric key an identifier rather than a measurement
ID_WORD = re.compile(r"(^| )(id|ids|sampleid|seqid|label|name|accession|barcode)( |$)")

p = argparse.ArgumentParser()
p.add_argument("diseases", nargs="*", default=["Colorectal cancer", "Crohn's disease"])
p.add_argument("--groups", default="disease,control", help="taxa-analysis groups to take samples from")
p.add_argument("--min-samples", type=int, default=1, help="keep keys present in at least this many samples")
p.add_argument("--unique-min", type=int, default=2,
               help="drop keys in at least this many samples whose value differs in every one of them")
p.add_argument("--cut", type=float, default=0.5, help="cosine distance at which the key dendrogram is cut")
p.add_argument("--top", type=int, default=10, help="clusters closest to the disease label to highlight")
p.add_argument("--neighbors", type=int, default=15, help="UMAP n_neighbors")
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
if args.top > N_HUE * len(MARKERS):
    p.error(f"--top is at most {N_HUE * len(MARKERS)} (hue x marker combinations)")
out = f"{common.RESULTS}/metadata_keys"
os.makedirs(out, exist_ok=True)

# ---- samples of each disease's taxa analysis ----
groups = args.groups.split(",")
samples = {}
for d in args.diseases:
    ab = pd.read_csv(f"{common.out_dir(d, 'taxa', site='')}/abundance.tsv.gz", sep="\t", usecols=["acc", "group"])
    samples[d] = set(ab.acc[ab.group.isin(groups)])
    print(f"{d}: {len(samples[d])} samples ({args.groups})")
accs = sorted(set().union(*samples.values()))


# ---- every metadata key per sample with its value, dictionary keys included ----
def walk(v, kv):
    if isinstance(v, dict):
        for k, x in v.items():
            kv[k].add(json.dumps(x, sort_keys=True))
            walk(x, kv)
    elif isinstance(v, list):
        for x in v:
            walk(x, kv)


con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("s", pd.DataFrame({"acc": accs}))
rows = con.sql("SELECT DISTINCT ON (acc) acc, jattr FROM sra_metadata JOIN s USING (acc)").fetchall()
acc_keys = {}
for acc, jattr in rows:
    acc_keys[acc] = defaultdict(set)
    walk(json.loads(jattr) if jattr else {}, acc_keys[acc])
print(f"{len(acc_keys)} of {len(accs)} samples have metadata")


def clean(k):
    return re.sub(r"_+", " ", SUFFIX.sub("", k)).strip().lower() or k


raw = defaultdict(set)                          # cleaned label -> raw keys
values = defaultdict(lambda: defaultdict(set))  # cleaned label -> sample -> its values
for acc, kv in acc_keys.items():
    for k, v in kv.items():
        raw[clean(k)].add(k)
        values[clean(k)][acc] |= v
have = {l: set(v) for l, v in values.items()}
keys = pd.DataFrame({"label": sorted(raw)})
keys["raw_keys"] = keys.label.map(lambda l: ",".join(sorted(raw[l])))
keys["n_samples"] = keys.label.map(lambda l: len(have[l]))
keys["n_values"] = keys.label.map(lambda l: len({frozenset(v) for v in values[l].values()}))


def is_number(v):
    v = json.loads(v)
    try:
        return all(float(x) == float(x) for x in (v if isinstance(v, list) else [v]))
    except (TypeError, ValueError):
        return False


keys["numeric"] = keys.label.map(lambda l: all(is_number(x) for v in values[l].values() for x in v))
for d in args.diseases:
    keys[f"n_{d.replace(' ', '_').replace(chr(39), '')}"] = keys.label.map(lambda l: len(have[l] & samples[d]))
# identifier-like keys: a different value in every sample that has them (sample ids, timestamps); numeric ones
# are measurements (age, BMI, nutrient intake) and kept unless the name says id / label / name
ident = ((keys.n_values == keys.n_samples) & (keys.n_samples >= args.unique_min)
         & (~keys.numeric | keys.label.str.contains(ID_WORD)))
keys[ident].sort_values("n_samples", ascending=False).to_csv(f"{out}/unique_keys.tsv", sep="\t", index=False)
keys = keys[~ident & (keys.n_samples >= args.min_samples)].reset_index(drop=True)
print(f"{len(raw)} distinct keys after merging suffixes; dropped {ident.sum()} identifier-like ones with a unique "
      f"value per sample (in >= {args.unique_min} samples); {len(keys)} left in >= {args.min_samples} samples")

# ---- BioLORD embedding, similarity, clusters, UMAP ----
model = SentenceTransformer("FremyCompany/BioLORD-2023", device="cuda" if torch.cuda.is_available() else "cpu")
emb = model.encode(keys.label.tolist(), batch_size=256, show_progress_bar=True)
sim = model.similarity(emb, emb).numpy()
np.save(f"{out}/similarity.npy", sim)
dist = np.clip(1 - sim, 0, None)
np.fill_diagonal(dist, 0)
dist = (dist + dist.T) / 2

cl = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="average",
                             distance_threshold=args.cut).fit_predict(dist)
# renumber clusters largest first (ties by summed samples)
size = pd.DataFrame({"c": cl, "n": keys.n_samples}).groupby("c").agg(size=("n", "size"), samples=("n", "sum"))
order = size.sort_values(["size", "samples"], ascending=False).index
keys["cluster"] = pd.Series(cl).map({c: i for i, c in enumerate(order)})

xy = umap.UMAP(n_neighbors=min(args.neighbors, len(keys) - 1), metric="precomputed",
               random_state=args.seed).fit_transform(dist)
keys["umap_x"], keys["umap_y"] = xy[:, 0], xy[:, 1]
keys.sort_values(["cluster", "n_samples"], ascending=[True, False]).to_csv(f"{out}/keys.tsv", sep="\t", index=False)

# ---- clusters closest to the disease label ----
# anchor: keys the disease labels were matched on, weighted by the samples each labeled
disease_keys = {k for (k,) in con.sql("SELECT DISTINCT k FROM disease_bmi_fields WHERE category = 'disease_health'")
                .fetchall()}
matched = con.sql("SELECT acc, matched_keys FROM selection JOIN s USING (acc)").fetchall()
weight = pd.Series([clean(k) for _, ks in matched if ks for k in ks.split(";") if k in disease_keys]).value_counts()
weight = weight[weight.index.isin(keys.label)]
print("disease-label keys (samples labeled):", weight.to_dict())
unit = emb / np.linalg.norm(emb, axis=1, keepdims=True)
row = pd.Series(range(len(keys)), index=keys.label)
anchor = (unit[row[weight.index]] * weight.values[:, None]).sum(0)
anchor /= np.linalg.norm(anchor)
centroid = pd.DataFrame(unit).groupby(keys.cluster).mean()
score = centroid.values @ anchor / np.linalg.norm(centroid.values, axis=1)

clusters = (keys.sort_values("n_samples", ascending=False).groupby("cluster")
            .agg(size=("label", "size"), samples=("n_samples", "sum"),
                 label_keys=("label", lambda l: " | ".join(x for x in l if x in weight.index)),
                 labels=("label", " | ".join)))
clusters.insert(2, "disease_similarity", score)
clusters.insert(3, "disease_rank", clusters.disease_similarity.rank(ascending=False, method="first").astype(int))
clusters.to_csv(f"{out}/clusters.tsv", sep="\t")
top = clusters.nsmallest(args.top, "disease_rank")
(keys[keys.cluster.isin(top.index)].assign(disease_rank=lambda k: k.cluster.map(top.disease_rank))
 .sort_values(["disease_rank", "n_samples"], ascending=[True, False])
 .to_csv(f"{out}/disease_clusters.tsv", sep="\t", index=False))
print(f"{len(clusters)} clusters at cosine distance {args.cut}; closest {args.top} to the disease label:")
print(top.drop(columns="label_keys").to_string(max_colwidth=120))

# ---- labeled UMAP ----
side = max(10, np.sqrt(len(keys)) * 1.1)
fig, ax = plt.subplots(figsize=(side, side))
rank = keys.cluster.map(clusters.disease_rank)
area = 8 + 40 * np.log10(keys.n_samples) / max(np.log10(keys.n_samples.max()), 1)
rest = rank > args.top
ax.scatter(keys.umap_x[rest], keys.umap_y[rest], s=area[rest], c=GREY, edgecolors="white", linewidths=0.5, zorder=2)
for r in range(1, args.top + 1):
    m = rank == r
    color, marker = style(r)
    ax.scatter(keys.umap_x[m], keys.umap_y[m], s=area[m] * 2, c=color, marker=marker,
               edgecolors="white", linewidths=0.5, zorder=2)
texts = [ax.text(x, y, l, fontsize=8 if r <= args.top else 6, color="#2b2b2b",
                 fontweight="bold" if r <= args.top else "normal", zorder=3)
         for x, y, l, r in zip(keys.umap_x, keys.umap_y, keys.label, rank)]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4), expand=(1.05, 1.2))
for r, (c, t) in enumerate(top.iterrows(), 1):
    members = keys[keys.cluster == c].nlargest(3, "n_samples").label
    color, marker = style(r)
    ax.scatter([], [], s=30, c=color, marker=marker,
               label=f"{r}. {', '.join(members)} ({t['size']} keys, similarity {t.disease_similarity:.2f})")
ax.scatter([], [], s=30, c=GREY, label="other clusters")
ax.legend(loc="upper left", fontsize=9, frameon=False,
          title=f"{args.top} clusters closest to the disease label (top keys by samples)")
ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2",
       title=f"BioLORD-2023 metadata keys of {len(acc_keys)} {' + '.join(args.diseases)} samples "
             f"({len(keys)} keys; point size = samples with the key)")
for s in ax.spines.values():
    s.set_visible(False)
fig.tight_layout()
fig.savefig(f"{out}/umap.png", dpi=200)
fig.savefig(f"{out}/umap.pdf")
print(f"wrote {out}/")
