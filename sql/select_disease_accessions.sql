-- Selects WGS-metagenome accessions (Koslicki_R1 manifest) whose BioSample attributes carry
-- host health / disease-state / BMI information. Run inside a DuckDB that already has the
-- `kv` table (acc, k, v) built by build_disease_db.sql.

CREATE OR REPLACE TABLE disease_kv AS
SELECT acc, k, trim(v) AS v,
  CASE
    WHEN regexp_matches(k, 'bmi|body_mass_index|weight_cat|obes') THEN 'bmi'
    ELSE 'disease_health'
  END AS category
FROM kv
WHERE
  (
    -- disease / diagnosis / health-state fields
    regexp_matches(k, '(^|_)(host_)?(disease|diseases|diagnosis|diagnosed|health_state|health_status|disease_stat|disease_status|clinical_condition|condition)(_|$)')
    OR regexp_matches(k, '(^|_)(ibd|crohn|crohns|colitis|celiac|diabetes|cancer|tumor|is_tumor|mental_illness|epilepsy|infected|healthy|healthy_uc_cd|case_or_control|case_status|case_control|host_phenotype|phenotype|sars_cov2_status|covid_status|dysbiosis|autism|parkinson|cirrhosis|hiv_status)(_|$)')
    -- BMI / weight category fields
    OR regexp_matches(k, '(^|_)(host_)?(bmi|body_mass_index|bmi_cat|weight_cat|obesity|obese)(_|$)')
  )
  -- exclude keys that describe something other than the sampled host's state
  AND NOT regexp_matches(k, 'mom_|mother|maternal|prepreg|father|donor_|recipient_|cleaning|film|soil|extraction|storage|growth|pcr|fermentation|treatment|seniority|onset|medication|_med_|duration|icd|past_3|family|relative|roommate|pet|dog|cat_|plant|leaf|root|sample_type|submitter|insdc|ena_|ngdc|reactor|experimental_condition|sample_id|_date|age_at|body_condition|initial_body|skin_condition|isolation|tumor_(her2|er|hr|pr|groups)')
  AND NOT regexp_matches(k, '^covid_(symptoms|patient|chronic)')
  -- drop missing-value placeholders
  AND v IS NOT NULL
  AND NOT regexp_matches(lower(trim(v)),
      '^(|na|n/a|nan|none|null|-|--|\.|\?|missing.*|not[ _]applicable.*|not[ _]provided.*|not[ _]collected.*|not[ _]available.*|not[ _]determined|unknown|unspecified|restricted access|labcontrol.*|lab control|blank|no data|nd|not reported|not recorded|undetermined|other)$');

CREATE OR REPLACE TABLE disease_accessions AS
SELECT acc,
       bool_or(category = 'disease_health') AS has_disease_health,
       bool_or(category = 'bmi')            AS has_bmi,
       string_agg(DISTINCT k, ';')          AS matched_keys
FROM disease_kv GROUP BY acc;
