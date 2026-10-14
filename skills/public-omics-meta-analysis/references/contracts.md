# Data contracts

## Cohort manifest

`cohort_id, accession, repository, publication, country, ancestry, tissue, disease_definition, comparator, assay, platform, n_donors, n_samples, paired, role, overlap_group, inclusion_status, exclusion_reason`

## Sample metadata

`cohort_id, sample_id, donor_id, tissue, region, condition, stage, treatment, timepoint, sex, age, ancestry, batch, platform, original_label, qc_status`

## Effect table

`cohort_id, feature_id, contrast, effect, se, p, n_case, n_control, direction, model, covariates, qc_status`

## Outputs

- search and screening ledger
- cohort and sample manifests
- per-cohort QC and effect tables
- meta-analysis estimates, heterogeneity, and prediction intervals
- leave-one-study-out and influence diagnostics
- discovery/replication summary
- reproducibility manifest

## Target handoff

Freeze candidate-selection rules before drug ranking. Export `target, contrast, tissue, disease_stage, meta_effect, meta_se, meta_p, q_value, heterogeneity_i2, prediction_interval, loso_stable, replication_status, orthogonal_validation, direction, evidence_limitations`. Expression evidence alone is contextual, not causal. Define minimum replication and orthogonal-validation criteria prospectively; downgrade rather than discard candidates that fail them.

## Operational ledger

For batch retrieval record API/query, release, cache key, checksum, pagination, rate-limit response, retry count, and terminal status. Use bounded parallelism, resumable downloads, immutable cached inputs, and per-cohort completion markers.
