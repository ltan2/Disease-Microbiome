#!/bin/bash
set -euo pipefail

wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz

sourmash sketch dna \
  -p k=31,scaled=1000 \
  --name GRCh38 \
  -o GRCh38.k31.zip \
  GCF_000001405.26_GRCh38_genomic.fna.gz
