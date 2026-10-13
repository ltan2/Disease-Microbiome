#!/bin/bash
# Write (acc, disease, status = case|control) for every disease to $1 by running
# sql/disease_groups.sql in memory against the read-only metadata DB.
set -euo pipefail
duckdb=${DUCKDB:-/home/grads/lbt5343/miniconda3/envs/logan/bin/duckdb}
db=/scratch/lbt5343/disease-sra/metadata/disease_sra_metadata.duckdb
sql=/home/grads/lbt5343/Disease_Inv/sql/disease_groups.sql
mkdir -p "$(dirname "$1")"
{
  echo "ATTACH '$db' AS m (READ_ONLY); USE m;"
  # the SQL creates a regular table; make it TEMP so nothing is written to the read-only DB
  sed 's/CREATE OR REPLACE TABLE disease_labels/CREATE OR REPLACE TEMP TABLE disease_labels/' "$sql"
  echo "COPY (SELECT acc, disease, status FROM disease_labels ORDER BY ALL) TO '$1' (HEADER, DELIMITER '\t');"
} | "$duckdb"
