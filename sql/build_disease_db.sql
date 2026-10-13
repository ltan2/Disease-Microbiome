-- Builds /scratch/lbt5343/disease-sra/metadata/disease_sra_metadata.duckdb: SRA metadata for every
-- WGS-metagenome accession (Koslicki_R1 manifest) that carries host disease/health-state or BMI
-- information, plus everything listed in Disease_Inv/bmi_accessions.txt.
-- Run from Disease_Inv/:  duckdb <db> < sql/build_disease_db.sql
SET threads = 16;
SET memory_limit = '64GB';
SET temp_directory = '/scratch/lbt5343/disease-sra/metadata/tmp';

ATTACH '/home/grads/lbt5343/Koslicki_R1/sra_metadata.db' AS src (READ_ONLY);

CREATE OR REPLACE TABLE wgs_acc AS
SELECT column0 AS acc FROM read_csv('/home/grads/lbt5343/Koslicki_R1/wgs_metagenome_accessions.txt',
                                     header = false, columns = {'column0': 'VARCHAR'});

-- long-form BioSample/run attributes of every WGS accession (input to the selection query)
CREATE OR REPLACE TABLE kv AS
SELECT m.acc, lower(trim(m.a.k)) AS k, m.a.v
FROM (SELECT acc, unnest(attributes) AS a FROM src.sra_metadata WHERE acc IN (SELECT acc FROM wgs_acc)) m;

.read sql/select_disease_accessions.sql

CREATE OR REPLACE TABLE bmi_file AS
SELECT DISTINCT trim(column0) AS acc
FROM read_csv('bmi_accessions.txt', header = false, columns = {'column0': 'VARCHAR'})
WHERE trim(column0) <> '';

-- final selection: metadata-derived set UNION the user-provided BMI list
CREATE OR REPLACE TABLE selection AS
SELECT a.acc,
       coalesce(d.has_disease_health, false) AS has_disease_health,
       coalesce(d.has_bmi, false)            AS has_bmi,
       a.acc IN (SELECT acc FROM bmi_file)   AS in_bmi_file,
       a.acc IN (SELECT acc FROM wgs_acc)    AS in_wgs_manifest,
       d.matched_keys
FROM (SELECT acc FROM disease_accessions UNION SELECT acc FROM bmi_file) a
LEFT JOIN disease_accessions d USING (acc);

-- full SRA row for every selected accession
CREATE OR REPLACE TABLE sra_metadata AS
SELECT m.*, m.organism ILIKE 'human%' AS human_host
FROM src.sra_metadata m WHERE m.acc IN (SELECT acc FROM selection);

-- every attribute of the selected accessions (includes env_broad_scale/env_local_scale/env_medium,
-- geo, age, sex, diet, ... for downstream covariates)
CREATE OR REPLACE TABLE attributes AS
SELECT m.acc, lower(trim(m.a.k)) AS k, m.a.v
FROM (SELECT acc, unnest(attributes) AS a FROM src.sra_metadata WHERE acc IN (SELECT acc FROM selection)) m;

-- only the disease/health/BMI fields that qualified an accession, labelled by category
CREATE OR REPLACE TABLE disease_bmi_fields AS
SELECT * FROM disease_kv WHERE acc IN (SELECT acc FROM selection);

-- convenience view of environment fields
CREATE OR REPLACE VIEW env_fields AS
SELECT * FROM attributes
WHERE regexp_matches(k, '^(env_broad_scale|env_local_scale|env_medium|env_biome|env_feature|env_material|host|host_body_site|body_site|isolation_source|geo_loc_name|lat_lon)');

-- one row per accession: selection flags + key SRA columns
CREATE OR REPLACE VIEW sample_summary AS
SELECT s.*, m.organism, m.human_host, m.bioproject, m.biosample, m.sra_study, m.mbases,
       m.geo_loc_name_country_calc, m.releasedate
FROM selection s LEFT JOIN sra_metadata m USING (acc);

DROP TABLE kv;
DROP TABLE disease_kv;
DROP TABLE disease_accessions;
DROP TABLE bmi_file;
DETACH src;

SELECT count(*) AS selected, count(*) FILTER (WHERE has_disease_health) AS disease_health,
       count(*) FILTER (WHERE has_bmi) AS bmi, count(*) FILTER (WHERE in_bmi_file) AS in_bmi_file,
       count(*) FILTER (WHERE in_wgs_manifest) AS in_wgs_manifest
FROM selection;
