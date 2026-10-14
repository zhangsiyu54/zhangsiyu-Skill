---
name: gwas-causal-target-analysis
description: Analyze GWAS summary statistics and molecular QTL data with phenotype harmonization, allele QC, LD-aware locus definition, fine-mapping, Mendelian randomization, colocalization, TWAS/PWAS triangulation, cross-cohort replication, and target prioritization. Use for UKB, IEU OpenGWAS, FinnGen, GWAS Catalog, eQTL, pQTL, mQTL, drug-target MR, causal variants, causal genes, genetic correlation, or variant-to-gene analysis.
---

# GWAS causal target analysis

Build an auditable chain from phenotype to locus to molecular trait to target. Never equate association, MR, colocalization, or nearest-gene mapping with proof of causality.

## Start with a written design

1. Define exposure, outcome, ancestry, genome build, unit, case definition, covariates, discovery cohort, and independent replication cohort.
2. Create a cohort-overlap table. Treat IEU studies containing UKB and FinnGen-UKB meta-analyses as non-independent from UKB.
3. Pre-specify primary tests, sensitivity tests, significance thresholds, and direction of the proposed intervention.
4. Record provenance, release, study ID, sample size, ancestry, build, allele coding, and access date for every dataset.

## Run the workflow

1. Audit the local runtime before analysis:
   `python3 "${CODEX_HOME:-$HOME/.codex}/skills/gwas-causal-target-analysis/scripts/audit_environment.py" --output environment_audit.json`.
   Treat missing optional tools as capability gaps, not proof that the data are invalid.
2. Validate each raw summary-statistic file before harmonization:
   `python3 "${CODEX_HOME:-$HOME/.codex}/skills/gwas-causal-target-analysis/scripts/validate_sumstats.py" INPUT --report qc_report.json --normalized normalized.tsv`.
   Review inferred column mappings and all errors. The normalized table is a schema/QC staging artifact, not allele harmonization.
   Use the same validator for molecular QTL summary statistics. Manually verify QTL platform, molecular trait, tissue/cell context, cis/trans definition, genome build, and study-specific fields because the generic schema validator cannot establish them.
   Accept biallelic SNPs and sequence-resolved indels containing only A/C/G/T. Reject symbolic, missing, or unsplit multi-allelic allele fields; split multi-allelic records upstream before continuing.
3. Validate the study manifest:
   `python3 "${CODEX_HOME:-$HOME/.codex}/skills/gwas-causal-target-analysis/scripts/validate_manifest.py" study_manifest.csv --report manifest_report.json`.
   Resolve errors and document overlap warnings before causal analyses.
4. Harmonize variants. Resolve build, chromosome notation, rsID, effect allele, other allele, beta/OR, SE, P, EAF, N, INFO, and phenotype units.
5. Remove or flag impossible alleles, duplicate variants, low-quality imputation, ambiguous palindromic variants without reliable frequency, inconsistent effect directions, and extreme allele-frequency mismatches.
6. Define loci using ancestry-matched LD. Keep LD reference, clumping thresholds, windows, and software versions in the manifest.
7. Replicate lead variants or locus-level signals in an independent cohort before target nomination where possible.
8. Fine-map loci with a method compatible with the inputs. Diagnose LD mismatch; do not interpret credible sets when summary statistics and LD reference are incompatible.
9. Run MR only when the instrument and estimand are defensible. Report strength, heterogeneity, horizontal pleiotropy, directionality, outliers, leave-one-out results, and sample overlap.
10. Run colocalization using full regional data and trait-appropriate priors. Prefer conditional or SuSiE-aware methods for multiple signals. Report prior sensitivity, H3/H4 evidence, regional coverage, and variant density.
11. Triangulate with eQTL/pQTL/mQTL, TWAS/PWAS, chromatin links, tissue relevance, single-cell expression, and spatial localization.
12. Score targets from independent evidence categories; never count correlated datasets as separate confirmations.

The bundled validators use only the Python standard library. Exit code `0` means valid, `1` means validation errors, and `2` means invocation or file-format failure. JSON reports are the audit record; preserve them with commands and logs.

Read [references/workflow.md](references/workflow.md) for method selection and mandatory sensitivity analyses. Read [references/contracts.md](references/contracts.md) before creating input tables or final deliverables.

## Hand off

- Send cell-type and tissue hypotheses to `singlecell-analysis`, `spatial-transcriptomics-analysis`, or `large-scale-singlecell-engineering`.
- Send molecular target tables to `drug-target-repurposing`.
- Send pathway-level results to `enrichment-analysis` or `network-pathway-analysis`.

## Completion gate

Do not call the analysis complete without phenotype harmonization, ancestry/build checks, overlap assessment, independent replication status, multiple-testing control, sensitivity analyses, null/failed-result reporting, and a reproducible manifest.
