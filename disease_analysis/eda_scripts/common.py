"""Paths and helpers shared by the sketch-similarity and taxa scripts."""
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import sourmash

SKETCHES = "/scratch/lbt5343/all-res/sketches/unitigs"
YACHT = "/scratch/lbt5343/all-res/yacht_results/unitigs"
KO = "/scratch/lbt5343/all-res/ko_profiles/unitigs"
HUMAN = "/home/grads/lbt5343/Disease_Inv/GRCh38.k31.zip"
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NOHUMAN = f"{RESULTS}/nohuman_unitigs"
METADATA_DB = "/scratch/lbt5343/disease-sra/metadata/disease_sra_metadata.duckdb"
LABELS = f"{RESULTS}/disease_labels.tsv"
YACHT_PRESENCE = f"{RESULTS}/yacht_presence.tsv.gz"
YACHT_SAMPLES = f"{RESULTS}/yacht_samples.txt"
YACHT_CONTAINMENT = f"{RESULTS}/yacht_containment.tsv.gz"
GROUPS = ["disease", "control", "other_control"]
# YACHT min_coverage levels used; 0 is left out since any single matching k-mer passes it
COVERAGES = ["1.0", "0.5", "0.25", "0.125", "0.0625", "0.03125", "0.015625"]
# name of the presence-level column per data kind (None: plain presence, no levels)
LEVEL = {"taxa": "coverage", "func": None}


# body site from the SRA organism (sra_metadata.organism); human metagenome / metagenome do not say
SITE_ORGANISMS = {
    "gut": ["human gut metagenome", "gut metagenome", "feces metagenome", "human feces metagenome"],
    "oral": ["human oral metagenome", "human saliva metagenome", "oral metagenome"],
}
UNKNOWN_ORGANISMS = ["human metagenome", "metagenome"]


def body_site(organism):
    """gut / oral / unknown (no organism, or one that names no site) / other (any other site or host)."""
    for site, organisms in SITE_ORGANISMS.items():
        if organism in organisms:
            return site
    return "unknown" if pd.isna(organism) or organism in UNKNOWN_ORGANISMS else "other"


def out_dir(disease, kind, site=None):
    """results/<disease>/<kind>; disease spaces -> '_' and apostrophes dropped (same rule as the Makefile).

    With a body site (site=, else the SITE environment variable) results/<disease>/site_<site>/<kind>,
    where 29_site_filter.py writes that site's samples; site="" forces the unfiltered directory."""
    site = os.environ.get("SITE", "") if site is None else site
    d = f"{RESULTS}/{disease.replace(' ', '_').replace(chr(39), '')}" + (f"/site_{site}" if site else "") + f"/{kind}"
    os.makedirs(d, exist_ok=True)
    return d


def split_groups(disease):
    """{group: set of accs} for one disease.

      disease        case for DISEASE
      control        control for DISEASE (and not a case of it)
      other_control  control for some other disease, never a case of any disease, no DISEASE label
    """
    lab = pd.read_csv(LABELS, sep="\t")
    if disease not in set(lab.disease):
        raise SystemExit(f"unknown disease {disease!r}; labels are listed in {LABELS}")
    this = lab[lab.disease == disease]
    case = set(this.acc[this.status == "case"])
    return {
        "disease": case,
        "control": set(this.acc[this.status == "control"]) - case,
        "other_control": set(lab.acc[lab.status == "control"]) - set(lab.acc[lab.status == "case"]) - set(this.acc),
    }


def nohuman_path(acc):
    return f"{NOHUMAN}/{acc}.sig.zip"


def ko_path(acc):
    return f"{KO}/{acc}.unitigs_ko_profiles.csv"


def read_ko(acc):
    """KO profile of one sample as (acc, ko, abundance) rows; header-only file = no KOs."""
    df = pd.read_csv(ko_path(acc)).rename(columns={"ko_id": "ko"})
    df.insert(0, "acc", acc)
    return df


def read_presence(out, group):
    """presence_<group>.tsv.gz as (df, feature column, level columns).

    Columns are acc, feature, then one 0/1 column per level; without level columns every row is
    a plain presence, given a single 'present' level.
    """
    df = pd.read_csv(f"{out}/presence_{group}.tsv.gz", sep="\t", dtype={"acc": str})
    if df.shape[1] == 2:
        df["present"] = 1
    return df, df.columns[1], list(df.columns[2:])


def load_hashes(path):
    """Sorted uint64 hashes of the single signature in `path`."""
    sig = next(iter(sourmash.load_file_as_signatures(path)))
    return np.sort(np.fromiter(sig.minhash.hashes, dtype=np.uint64))


def load_nohuman(accs, procs):
    with Pool(procs) as p:
        return p.map(load_hashes, [nohuman_path(a) for a in accs], chunksize=64)


GTDB_TAXONOMY = [f"{RESULTS}/gtdb_r232/{d}_taxonomy_r232.tsv.gz" for d in ["bac120", "ar53"]]
RANK_PREFIX = {"d__": "domain", "p__": "phylum", "c__": "class", "o__": "order", "f__": "family", "g__": "genus",
               "s__": "species"}


def read_gtdb():
    """GTDB lineage per YACHT taxon (genome accession without the RS_/GB_ prefix): taxon, then one
    column per rank holding the prefixed name (e.g. genus = 'g__Fusobacterium')."""
    tax = pd.concat(pd.read_csv(p, sep="\t", header=None, names=["taxon", "lineage"]) for p in GTDB_TAXONOMY)
    tax["taxon"] = tax.taxon.str[3:]
    ranks = tax.pop("lineage").str.split(";", expand=True)
    ranks.columns = list(RANK_PREFIX.values())
    return pd.concat([tax, ranks], axis=1)
