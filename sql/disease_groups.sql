-- Normalises the messy disease/health metadata into (accession, disease, status = case | control)
-- and summarises case/control counts per disease. Run against the metadata DB:
--   duckdb /scratch/lbt5343/disease-sra/metadata/disease_sra_metadata.duckdb < sql/disease_groups.sql
--
-- Three kinds of evidence:
--  1. Flag/questionnaire fields named after one condition (American Gut `ibd_sam`, `diabetes_sam`,
--     `sars_cov2_status`, ...): yes/true/1/positive/"diagnosed by a medical professional" = case,
--     no/false/0/negative/"I do not have this condition" = control *for that condition*.
--     Self-diagnosed / alternative-practitioner answers are left out.
--  2. Free-text fields (`host_disease`, `diagnosis`, `phenotype`, ...): the value is matched against
--     the dictionary below for a named disease, or against the control vocabulary.
--  3. Generic values with no disease named (control, healthy, case, diseased, positive, ...): the
--     disease comes from the study (BioProject) -- a generic control is counted as a control for every
--     disease that has named cases in the same BioProject; a generic case is resolved only if its
--     BioProject has exactly one named disease, else 'Unspecified disease'.
-- An accession that is a case is never also counted as a generic control.

CREATE OR REPLACE TEMP TABLE dict(label VARCHAR, rx VARCHAR, prio INT);
INSERT INTO dict VALUES
  ('Crohn''s disease',          'crohn|^_?cd_?$|ileum|cd_disease|colon alone|ileocolon', 10),
  ('Ulcerative colitis',        'ulcerative colitis|^_?uc_?$|pancolit|left side colit|proctitis|pouchitis', 10),
  ('IBD (unspecified)',         '(^|[^a-z])ibd($|[^a-z])|inflammatory bowel|indeterminate colitis|ibd-u|^colitis$|colitis_sam|crohns_disease_sam', 20),
  ('C. difficile infection',    'c\.? ?difficile|^cdi$', 10),
  ('Colorectal cancer',         'colorectal (cancer|carcinoma)|^crc|colon cancer|rectal cancer', 10),
  ('Colorectal adenoma',        'adenoma', 15),
  ('Irritable bowel syndrome',  'irritable bowel|^ibs|ibsd|piibs', 10),
  ('Parkinson''s disease',      'parkinson|^pd$', 10),
  ('Type 1 diabetes',           'type (1|i) diabetes|^t1d|t1d \+|^ot1d', 10),
  ('Type 2 diabetes',           'type (2|ii) diabetes|^t2d|type2', 10),
  ('Diabetes (type unspecified)','diabet', 30),
  ('Lung cancer',               'nsclc|lung (cancer|adenocarcinoma)|bronchus and lung|lung neoplasm', 10),
  ('Breast cancer',             'breast', 10),
  ('Prostate cancer',           'prostate', 10),
  ('Pancreatic cancer',         'pdac|pancreatic (cancer|neoplasm|adenocarcinoma)', 10),
  ('Blood cancer',              'leukemia|lymphoma|myeloma|^aml$|^all$|dlbcl', 10),
  ('Cancer (other/unspecified)','cancer|carcinoma|tumou?r|neoplasm|malignan', 40),
  ('Liver disease',             'cirrhosis|aclf|hepatitis|nash|nafld|fatty liver|cholangitis|liver', 20),
  ('Kidney disease',            'kidney|renal|ckd|dialysis', 20),
  ('COVID-19',                  'covid|sars-?cov|sars_cov2', 10),
  ('HIV',                       '(^|[^a-z])hiv', 10),
  ('Pneumonia',                 'pneumonia', 10),
  ('Periodontal disease',       'periodont|gingivit|peri-implantitis|caries', 10),
  ('Enteric infection/diarrhea','diarrh|cholera|etec|(campylobacter|salmonella|clostridium|shigella):\+|salmonellosis|campylobacteriosis|shigellosis|enteric infection', 15),
  ('Alcohol use disorder',      'alcohol use|alcoholic|_aud$|^aud$|alcohol dependence', 10),
  ('Schizophrenia',             'schizophren', 10),
  ('Depression',                'depress', 10),
  ('Bipolar disorder',          'bipolar', 10),
  ('Eating disorder',           'anorexia|bulimia', 10),
  ('Substance abuse',           'substance', 10),
  ('Mental illness (any)',      'mental_illness_sam', 50),
  ('Autism',                    'autis|asd', 10),
  ('Multiple sclerosis',        'multiple sclerosis|^rrms|^ms$', 10),
  ('Rheumatoid arthritis',      'rheumatoid', 10),
  ('Lupus (SLE)',               'lupus|^sle$', 10),
  ('Psoriasis',                 'psoria', 10),
  ('Atopic dermatitis',         'atopic dermatitis|eczema', 10),
  ('Asthma',                    'asthma', 10),
  ('COPD',                      'copd|chronic obstructive', 10),
  ('Lung disease (any)',        'lung_disease', 50),
  ('Cystic fibrosis',           'cystic fibrosis', 10),
  ('Allergy',                   'allerg', 15),
  ('Cardiovascular disease',    'cardiovascular|stroke|hypertens|heart failure|coronary|atheroscl|myocard', 15),
  ('Epilepsy',                  'epilep|seizure', 10),
  ('Celiac disease',            'celiac|coeliac', 10),
  ('Meningitis/encephalitis',   'meningitis|encephalitis', 10),
  ('Urinary tract infection',   'urinary tract infection', 10),
  ('H. pylori infection',       'pylori', 10),
  ('Fibromyalgia',              'fibromyalgia|fibromylagia', 10),
  ('Obesity/metabolic syndrome','metabolic syndrome|^obese$|^obesity$', 20),
  ('Chronic GI disease',        'chronic_gi_disease|gastrointestinal_tract_disorder|intestinal_disease', 50);

CREATE OR REPLACE TEMP TABLE r AS
SELECT f.acc, f.k, lower(trim(f.v)) AS v, m.bioproject, m.human_host
FROM disease_bmi_fields f JOIN sra_metadata m USING (acc)
WHERE f.category = 'disease_health'
  AND f.k NOT IN ('clinical_condition_sam', 'mental_illness_type_unspecified_sam')
  -- American Gut subset_* flags duplicate its per-condition answers (subset_healthy=false is a
  -- catch-all "not healthy" with no disease), so they are not used
  AND NOT starts_with(f.k, 'subset_');

-- field named after one condition -> its label
CREATE OR REPLACE TEMP TABLE flag_keys AS
SELECT DISTINCT r.k, arg_min(d.label, d.prio) OVER (PARTITION BY r.k) AS label
FROM (SELECT DISTINCT k FROM r) r JOIN dict d ON regexp_matches(r.k, d.rx)
WHERE regexp_matches(r.k, '^(subset_)?(ibd|diabetes|cardiovascular_disease|lung_disease|cancer|liver_disease|kidney_disease|mental_illness|mental_illness_type_[a-z_]+|epilepsy_or_seizure_disorder|sars_cov2_status|chronic_gi_disease|celiac_disease|crohns_disease|colitis|intestinal_disease|gastrointestinal_tract_disorder|autism|asthma|hiv_status|is_tumor)(_sam.*|_exp|_run)?$');

CREATE OR REPLACE TEMP TABLE yes_rx AS SELECT '^(yes|y|true|1|positive|diagnosed by a medical professional.*)$' AS rx;
CREATE OR REPLACE TEMP TABLE no_rx  AS SELECT '^(no|n|false|0|negative|i do not have this condition)$' AS rx;
-- generic control vocabulary (free-text fields)
CREATE OR REPLACE TEMP TABLE ctl_rx AS SELECT
 '^(_?healthy_?|healthy[ _]control.*|.*healthy control.*|control.*|.*control group.*|normal|hc|ctr|ctrl|nc|hv|non[-_ ]?ibd|unaffected|apparently healthy|health|healthy mucosa|no ?conditions|non-diabetic|nt1d|non-t1d.*|non-alcoholic control|bv negative|healthy volunteer.*|no disease|disease[- ]free|negative|no|none|[a-z. ]+:-(;[a-z. ]+:-)*)$' AS rx;  -- last: all pathogen tests negative
CREATE OR REPLACE TEMP TABLE case_rx AS SELECT
 '^(case|cases|diseased|disease|patient|positive|infected|yes|y|true|1|sick|ill|affected)$' AS rx;

-- 1. flag fields
CREATE OR REPLACE TEMP TABLE ev_flag AS
SELECT r.acc, fk.label AS disease,
       CASE WHEN regexp_matches(r.v, (SELECT rx FROM yes_rx)) THEN 'case'
            WHEN regexp_matches(r.v, (SELECT rx FROM no_rx))  THEN 'control' END AS status,
       r.k, r.v, r.bioproject, r.human_host
FROM r JOIN flag_keys fk USING (k);

-- 2+3. free-text fields (everything that isn't a flag field)
CREATE OR REPLACE TEMP TABLE ft AS SELECT * FROM r WHERE k NOT IN (SELECT k FROM flag_keys);
CREATE OR REPLACE TEMP TABLE ft_named AS
SELECT DISTINCT ON (ft.acc, ft.k) ft.acc, d.label AS disease, 'case' AS status, ft.k, ft.v, ft.bioproject, ft.human_host
FROM ft JOIN dict d ON regexp_matches(ft.v, d.rx)
WHERE NOT regexp_matches(ft.v, (SELECT rx FROM ctl_rx)) AND NOT regexp_matches(ft.v, '^(non|no)[-_ ]')
ORDER BY ft.acc, ft.k, d.prio;
-- `healthy`-named keys invert: subset_healthy=true / is_healthy=yes -> control
-- Generic values. If the *field name* names a disease (t1d_diagnosed, hpylori_diagnosis, ...),
-- yes/no on it is case/control for that disease. `healthy`-named fields invert:
-- is_healthy=yes -> control.
CREATE OR REPLACE TEMP TABLE ft_generic AS
WITH g AS (
    SELECT ft.*, (SELECT d.label FROM dict d WHERE regexp_matches(ft.k, d.rx) ORDER BY d.prio LIMIT 1) AS key_disease
    FROM ft WHERE NOT EXISTS (SELECT 1 FROM ft_named n WHERE n.acc = ft.acc AND n.k = ft.k)
)
SELECT acc, key_disease AS disease,
       CASE WHEN regexp_matches(k, 'healthy') AND regexp_matches(v, (SELECT rx FROM yes_rx)) THEN 'control'
            WHEN regexp_matches(k, 'healthy') AND regexp_matches(v, (SELECT rx FROM no_rx))  THEN 'case'
            WHEN key_disease IS NOT NULL AND regexp_matches(v, (SELECT rx FROM no_rx))  THEN 'control'
            WHEN regexp_matches(v, (SELECT rx FROM ctl_rx)) THEN 'control'
            WHEN regexp_matches(v, (SELECT rx FROM case_rx)) THEN 'case' END AS status,
       k, v, bioproject, human_host
FROM g;

CREATE OR REPLACE TEMP TABLE ev AS
SELECT * FROM ev_flag WHERE status IS NOT NULL
UNION ALL SELECT * FROM ft_named
UNION ALL SELECT * FROM ft_generic WHERE status IS NOT NULL;

-- named diseases per BioProject (cases only), to give generic rows a disease
CREATE OR REPLACE TEMP TABLE proj_dis AS
SELECT bioproject, disease FROM ev WHERE status = 'case' AND disease IS NOT NULL GROUP BY ALL;
CREATE OR REPLACE TEMP TABLE proj_one AS
SELECT bioproject, any_value(disease) AS disease FROM proj_dis GROUP BY bioproject HAVING count(*) = 1;

CREATE OR REPLACE TABLE disease_labels AS
WITH resolved AS (
    SELECT acc, disease, status, k, v, bioproject, human_host, 'field' AS how FROM ev WHERE disease IS NOT NULL
    UNION ALL   -- generic case: the study's single disease, else unspecified
    SELECT e.acc, coalesce(p.disease, 'Unspecified disease'), 'case', e.k, e.v, e.bioproject, e.human_host,
           CASE WHEN p.disease IS NULL THEN 'unresolved' ELSE 'study' END
    FROM ev e LEFT JOIN proj_one p USING (bioproject) WHERE e.disease IS NULL AND e.status = 'case'
    UNION ALL   -- generic control: control for every named disease in its study
    SELECT e.acc, coalesce(p.disease, 'Healthy (no case group in study)'), 'control', e.k, e.v, e.bioproject, e.human_host, 'study'
    FROM ev e LEFT JOIN proj_dis p USING (bioproject) WHERE e.disease IS NULL AND e.status = 'control'
), dedup AS (
    SELECT acc, disease, status, any_value(bioproject) bioproject, any_value(human_host) human_host,
           string_agg(DISTINCT k || '=' || v, '; ') evidence, min(how) how
    FROM resolved GROUP BY acc, disease, status
)
SELECT * FROM dedup d
WHERE NOT (d.status = 'control' AND (
        -- a case for this disease is not also its control
        EXISTS (SELECT 1 FROM dedup x WHERE x.acc = d.acc AND x.disease = d.disease AND x.status = 'case')
        -- a free-text/study case of anything is not a generic healthy control
     OR (d.how = 'study' AND EXISTS (SELECT 1 FROM dedup x WHERE x.acc = d.acc AND x.status = 'case'))));
