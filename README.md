# Disease-Microbiome

Tools for studying disease-associated microbiome differences using SRA metadata and Logan metagenome assemblies. The project selects samples with disease, health, or BMI metadata; processes assemblies into sequence sketches, KEGG Orthology (KO) profiles, and YACHT taxon detections; and compares disease and control groups through exploratory data analysis (EDA).

The repository contains code and configuration. Input databases, reference sketches, accession lists, and generated analysis results are not included.

## Project layout

| File or folder | Purpose |
| --- | --- |
| `sql/` | Select relevant SRA samples and assign disease-specific case/control labels. |
| `nextflow/` | Download and process Logan unitigs and contigs, recording progress for each sample. |
| `disease_analysis/eda_scripts/` | Build comparison groups, calculate microbiome metrics, run statistical tests, and produce figures and tables. |
| `README.md` | Project overview, file guide, and usage notes. |
| `.gitignore` | Excludes common Python caches, environments, build artifacts, and local configuration from Git. |
| `LICENSE` | MIT license, copyright 2026 ltan2. |

## Metadata SQL files

Paths in this table are relative to `sql/`.

| File | What it does |
| --- | --- |
| `build_disease_db.sql` | Builds a DuckDB database from an external SRA metadata database, a WGS-metagenome accession manifest, and a BMI accession list. Stores selected metadata, attributes, selection flags, and summary views; calls `select_disease_accessions.sql`. |
| `select_disease_accessions.sql` | Filters BioSample/run attributes for disease, health, and BMI information, excluding unrelated fields and missing-value placeholders. Creates intermediate accession-selection tables. |
| `disease_groups.sql` | Converts disease flags and free-text metadata into case/control labels. Uses BioProject context to resolve generic labels and records supporting evidence in `disease_labels`. |

## Nextflow pipeline files

Paths in this table are relative to `nextflow/`.

| File | What it does |
| --- | --- |
| `main.nf` | Pipeline entry point. Reads a todo TSV, routes requested downloads and processing stages, and passes existing unitig sketches directly to YACHT when available. |
| `nextflow.config` | Active configuration: reference/output paths, tool parameters, retries, resource scheduling, download backpressure, traces, and reports. Currently schedules up to 150 declared CPUs and an 880 GB Nextflow memory pool. |
| `nextflow.config.bak_100cpu` | Previous configuration with a 100-CPU scheduling pool and fewer concurrent YACHT tasks; retained as a backup. |
| `subworkflows/analyze.nf` | Connects decompression, sketching, KO profiling, YACHT, and cleanup. Runs only the stages requested for each sequence type. |
| `modules/fetch.nf` | Downloads requested assemblies from the public Logan S3 bucket. Handles missing/tiny files and limits the number of assemblies awaiting analysis. |
| `modules/decompress.nf` | Decompresses `.fa.zst` assemblies into shared FASTA files and removes the consumed compressed files from the work directory. |
| `modules/sketch.nf` | Creates abundance-aware sourmash DNA sketches and writes completion records. |
| `modules/funprofiler.nf` | Runs the `funcprofiler` command against KO reference sketches, producing KO profiles, prefetch CSVs, and completion records. |
| `modules/yacht.nf` | Runs YACHT taxon detection on sketches with memory requests based on sketch size. Publishes Excel results and records empty-sketch/no-hit outcomes. |
| `modules/cleanup.nf` | Deletes decompressed FASTA files after sketching and KO profiling finish, then releases their download-backpressure tokens. |
| `modules/ledger.nf` | Writes progress records for events such as downloads skipped because the assembly is empty. |
| `bin/write_ledger_row.sh` | Writes one CSV progress record containing accession, sequence type, stage, status, shard ID, and UTC timestamp. |

The pipeline sketches and functionally profiles both unitigs and contigs; YACHT runs on unitigs only. Published outputs are organized under `outdir` into `sketches/`, `ko_profiles/`, `yacht_results/`, and `ledger/`.

## Analysis files

Paths in this table are relative to `disease_analysis/eda_scripts/`. The numbered scripts form several analysis branches; they are not intended to be run as one numerical sequence. Both files numbered `18` are separate analyses.

| File | What it does / main output |
| --- | --- |
| `Makefile` | Coordinates dependencies and exposes targets for sketch, taxon, KO, metadata, and statistical analyses. |
| `common.py` | Shared data paths, disease-group definitions, coverage thresholds, and helpers for reading sketches and profiles. |
| `00_labels.sh` | Runs the disease-label SQL against a read-only metadata database and exports `results/disease_labels.tsv`. |
| `01_remove_human.py` | Removes hashes shared with the GRCh38 reference from unitig sketches, preserving abundance information. Writes nonempty sketches to `results/nohuman_unitigs/`. |
| `02_split_groups.py` | Creates `groups.tsv` for a disease and data type, retaining only samples with the required data. |
| `03_jaccard.py` | Calculates pairwise Jaccard similarities between human-filtered sketches; writes `sketch/jaccard.npy`. |
| `04_heatmap.py` | Plots a similarity/distance matrix with group boundaries and optional subsampling. |
| `05_summary.py` | Summarizes within-group and between-group matrix values with counts, means, and medians. |
| `06_pcoa_permanova.py` | Produces principal coordinates analysis (PCoA) plots and runs PERMANOVA and PERMDISP for all groups and pairwise comparisons. |
| `07_hash_counts.py` | Counts unique hashes per sample after human-hash removal; writes `sketch/hash_counts.tsv`. |
| `08_yacht_presence.py` | Parses YACHT workbooks into a shared taxon presence table at each coverage threshold and a list of samples with results. |
| `09_taxa_presence.py` | Splits taxon presence records into disease, control, and other-control tables. |
| `10_prevalence.py` | Calculates the percentage of samples containing each taxon or KO within each group. |
| `11_distinct_per_sample.py` | Counts distinct taxa or KOs per sample, including zero counts for samples without detections. |
| `12_ko_presence.py` | Converts unitig KO profiles into per-group presence records and binary KO-by-sample matrices. |
| `13_ko_abundance.py` | Combines KO profiles into a KO-by-sample abundance matrix, filling absent KOs with zero. |
| `14_clr.py` | Filters rare KOs, replaces zeros, and applies the centered log-ratio (CLR) transform; saves transformed values and matching sample groups. |
| `15_aitchison.py` | Calculates pairwise Aitchison distances from CLR-transformed KO profiles; writes `func/aitchison.npy`. |
| `16_seq_depth.py` | Extracts sequencing depth in megabases from SRA metadata for samples in an analysis group. |
| `17_taxa_jaccard.py` | Calculates taxon-presence Jaccard similarities separately at each YACHT coverage threshold, excluding samples with no detected taxa at that threshold. |
| `18_kruskal_dunn.py` | Tests per-sample counts using Kruskal–Wallis and pairwise Dunn tests with Holm correction; optionally repeats tests by coverage level. |
| `18_metadata_eda.py` | Summarizes group sizes, data availability, BioProjects, geography, sequencing characteristics, age, sex, and BMI; produces tables, plots, and confounding checks. |
| `19_permanova_table.py` | Combines KO and taxon PERMANOVA/PERMDISP results into a TSV and prints a Markdown summary with explained-variation estimates. |
| `20_stratified_permanova.py` | Compares disease versus control before and after accounting for BioProject, using samples with all three data types and qualifying projects. |
| `21_within_project.py` | Compares per-sample counts within BioProjects using Mann–Whitney tests, Cliff's delta, and a combined van Elteren test. |
| `22_latex_tables.py` | Exports count-test results for one or more diseases as LaTeX table fragments and a standalone LaTeX document. |
| `23_yacht_containment.py` | Extracts detected taxa's matching/exclusive k-mer counts and containment values from YACHT workbooks. |
| `24_taxa_abundance.py` | Calculates per-sample taxon relative-abundance proxies from normalized containment and weighted detection thresholds. |
| `plot_group_density.py` | Draws group density curves for counts or sequencing depth, with optional log transformation and coverage panels. |

The analysis groups are:

| Group | Definition |
| --- | --- |
| `disease` | Samples labeled as cases for the selected disease. |
| `control` | Samples labeled as controls for the selected disease, excluding its cases. |
| `other_control` | Controls for another disease, with no case label for any disease and no label for the selected disease. |

## Before running

The current paths and resource settings target the original Linux research server. Adapt them before using another environment:

- `sql/build_disease_db.sql`: source database, WGS accession manifest, temporary directory, and external `bmi_accessions.txt` input. Run this SQL from the repository root so its relative `.read` path resolves.
- `nextflow/nextflow.config`: KO/YACHT references, output and inflight directories, and resource settings.
- `disease_analysis/eda_scripts/common.py`: sketches, KO profiles, YACHT results, GRCh38 reference sketch, and metadata database.
- `disease_analysis/eda_scripts/00_labels.sh`: database and SQL paths; the DuckDB executable can be overridden with `DUCKDB`.
- `disease_analysis/eda_scripts/Makefile`: Python executable, overridable with `PY`.

The pipeline requires Nextflow with its Java runtime, Bash, AWS CLI, zstd, sourmash, `funcprofiler`, YACHT, and GNU/Linux command-line utilities. Analysis scripts use Python with NumPy, pandas, SciPy, Matplotlib, sourmash, scikit-bio, DuckDB, and an Excel reader such as openpyxl. Markdown table export also needs tabulate. Use GNU Make with grouped-target support for the analysis Makefile. This repository does not include a pinned dependency environment.

Comments reference external helpers such as `scripts/make_todo.py` and `scripts/run_disease.sh`; those files are not present here. In particular, the systemd hard resource cap described in the Nextflow configuration depends on that external launcher; it is not established by the included configuration alone.

## Running the pipeline

Provide a tab-separated todo file with the following header:

```text
accession	fetch_unitigs	fetch_contigs	sketch_unitigs	ko_unitigs	yacht_unitigs	sketch_contigs	ko_contigs	existing_unitigs_sig
```

Stage/download flags use `1` to request work and `0` otherwise. Leave `existing_unitigs_sig` empty unless supplying an existing unitig sketch for YACHT. Set fetch flags consistently with sketch/KO work that requires an assembly.

After configuring paths and resources, run from the repository root:

```bash
nextflow run nextflow --todo /path/to/todo.tsv
```

## Running analyses

After preparing the metadata database and processing outputs, run from the analysis directory:

```bash
cd disease_analysis/eda_scripts
make list-diseases PY=python
make sample-sketch-similar-viz DISEASE="Colorectal cancer" PY=python PROCS=8
```

Choose other Make targets according to the desired analysis:

| Target | Analysis |
| --- | --- |
| `sample-sketch-similar-viz` | Sketch similarity heatmap, summaries, ordination/tests, and hash-count plot. |
| `taxa-pres-abs` | Taxon prevalence, distinct counts, and coverage-specific ordination/tests. |
| `funprofile-pres-abs` | KO prevalence, distinct counts, and count tests. |
| `funprofile-abundance` | KO CLR transformation, Aitchison distances, heatmap, and ordination/tests. |
| `seq-depth` | Sequencing-depth tables and plots for taxon and KO sample groups. |
| `metadata-eda` | Metadata tables, figures, and study/geography confounding summaries. |
| `stratified-permanova` | Disease/control distance comparisons accounting for BioProject. |
| `count-tests` | Group count tests and within-BioProject comparisons. |
| `taxa-abundance` | Containment- and detection-score-based taxon abundance proxies. |

Each disease target accepts `DISEASE="..."` and `PY=python`. Additional Make variables include `PROCS`, `PERMS`, `MIN_PREV`, `N_SUB`, `MIN_PER_GROUP`, and `DETECT_WEIGHTS`; see the Makefile for defaults.

Generated analysis files go under `disease_analysis/eda_scripts/results/`. Shared label and YACHT tables live directly there; disease-specific outputs use `<disease>/sketch/`, `<disease>/taxa/`, `<disease>/func/`, and `<disease>/metadata/`, with some combined tables at the disease-directory level. Disease names use underscores instead of spaces and omit apostrophes.

## License

MIT; see `LICENSE`.
