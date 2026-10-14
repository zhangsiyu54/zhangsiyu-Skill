---
name: public-omics-meta-analysis
description: Discover, screen, harmonize, reanalyze, and meta-analyze public transcriptomic and multi-omics cohorts from GEO, SRA, ArrayExpress, TCGA, GTEx, and related repositories. Use for systematic public-dataset searches, eligibility criteria, metadata curation, batch-aware cohort analysis, effect-size synthesis, heterogeneity, leave-one-study-out validation, independent discovery/validation splits, or large multi-cohort RNA/protein/metabolite projects.
---

# Public omics meta-analysis

Separate dataset discovery, within-cohort analysis, and cross-cohort synthesis. Never merge matrices from unrelated studies as a substitute for meta-analysis.

## Build the protocol

1. Define disease, tissue, specimen, assay, species, treatment, time point, comparator, and minimum independent donor count.
2. Predefine exclusions: cell lines versus tissue, pooled samples, missing group labels, incompatible endpoints, duplicated cohorts, unavailable counts, or severe confounding.
3. Record search strings, repositories, dates, accession IDs, and reasons for exclusion.
4. Assign cohorts to discovery, replication, or contextual evidence before inspecting target results.

## Curate metadata

Create controlled fields for study, donor, sample, tissue, region, disease, stage, treatment, platform, library, sex, age, ancestry, batch, and pairing. Detect duplicate samples and publications reusing the same cohort. Keep original and harmonized labels.

## Analyze cohorts independently

1. Use raw counts when valid and available; otherwise document the processed scale and transformations.
2. Apply assay-specific QC and normalization within each cohort.
3. Fit a design reflecting pairing, batches, covariates, and repeated measures.
4. Export effect estimates, SE, sample sizes, direction, and QC flags rather than only significant gene lists.
5. Delegate assay-specific work to `bulkrna-analysis`, `proteomics-metabolomics-analysis`, `public-db-survival`, or another matching Skill.

## Synthesize evidence

1. Harmonize feature identifiers and effect directions.
2. Combine comparable effects using a justified fixed- or random-effects model.
3. Report heterogeneity, prediction intervals where useful, multiple-testing correction, and leave-one-study-out stability.
4. Keep vote counting or rank aggregation secondary when comparable effect sizes are unavailable.
5. Validate key findings in an independent cohort and orthogonal modality.

Read [references/protocol.md](references/protocol.md) for screening and duplication controls. Read [references/contracts.md](references/contracts.md) for tables, outputs, and the formal handoff to target/drug analysis.

## Executable meta-analysis

Use the dependency-light [scripts/meta_analyze.py](scripts/meta_analyze.py) when cohort-level effects and standard errors are comparable. It strictly validates every cohort/effect field in [references/contracts.md](references/contracts.md), numeric ranges, duplicate rows, cohort membership, cohort inclusion status, and effect QC. Only `included` cohorts with `pass`, `included`, or `ok` effect rows may enter synthesis. It then writes fixed-effect and DerSimonian-Laird random-effects estimates, confidence intervals, Q, I2, tau2, prediction intervals (for at least three cohorts), and leave-one-cohort-out results.

```bash
python3 scripts/meta_analyze.py --cohorts cohort_manifest.csv --effects effects.csv --out-prefix results/meta
```

The script intentionally does not harmonize biological contrasts or feature identifiers automatically. Confirm phenotype, tissue, assay, unit, and effect direction compatibility first. Treat DerSimonian-Laird as a transparent baseline; use a prospectively justified estimator for small or highly heterogeneous collections.

## Completion gate

Require a search ledger, inclusion flow, cohort manifest, metadata dictionary, per-cohort QC, effect-size tables, heterogeneity results, sensitivity analyses, independent-validation status, code/environment versions, and null findings.
