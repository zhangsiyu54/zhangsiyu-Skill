# Method selection and sensitivity checks

## Contents

- Summary-statistic QC
- Replication and meta-analysis
- MR
- Fine-mapping and colocalization
- Target triangulation

## Summary-statistic QC

Require chromosome, position or rsID, alleles, effect estimate, SE or convertible statistic, P value, and sample size. Prefer EAF and INFO. Confirm GRCh37 versus GRCh38 and lift over only with documented chain files and post-liftover validation. Represent SNPs and sequence-resolved indels with explicit A/C/G/T alleles. Split multi-allelic records into auditable biallelic rows upstream. Align all effects to the same effect allele.

## Replication and meta-analysis

Match phenotype definition, ancestry, effect scale, and allele direction. Pre-register genome-wide discovery, replication, multiple-testing, direction-consistency, and exact/proxy/locus-level criteria. Report those replication levels separately. Assess winner's curse and do not silently promote failed-replication loci. Use fixed or random effects according to the design and report heterogeneity. Do not label cohorts independent when participants overlap.

## MR

Report F statistics and instrument count. Use IVW as primary only when appropriate; add weighted median, MR-Egger, heterogeneity, intercept, leave-one-out, Steiger directionality, and outlier-robust analyses. For cis-drug-target MR, assess LD between instruments and use generalized models when variants are correlated. Interpret reverse MR as a distinct analysis.

## Fine-mapping and colocalization

Use ancestry-matched LD. Check allele alignment, regional completeness, effective sample size, and LD mismatch. For multiple signals, condition or use SuSiE-aware approaches. Explore plausible prior grids. Distinguish shared-causal-variant evidence from two associations in the same region.

## Target triangulation

Prefer cis-QTL for target inference and analyze cis/trans evidence separately. Check pQTL platform, assay range, sample size, protein-altering variants, and epitope-binding artifacts. Treat blood QTL as proxy evidence unless blood is the relevant tissue. Use separate evidence classes: human genetics, molecular QTL, functional genomics, cell/tissue specificity, perturbation, disease expression, pharmacology, and safety. Pre-register evidence weights, minimum advancement rules, and direction-conflict penalties. Preserve negative evidence.
