# Proteomics Metabolomics Advanced Modules

## QC

- Report missing values by sample and feature.
- Check sample clustering/PCA and batch/group separation.
- For proteomics, distinguish LFQ, TMT, DIA, phosphoproteomics and peptide/protein level.
- For metabolomics, distinguish targeted/untargeted and identified/unknown features.

## Missing Values

- State whether missingness is MCAR/MAR/MNAR-like when inferable.
- Use conservative imputation; do not impute without recording method.
- For heavy missingness, report filtered feature counts.

## Differential Analysis

- Proteomics: DEP/MSstats/limma-style models depending on design.
- Metabolomics: fold change, moderated testing, multiple testing; optional PLS-DA only as exploratory and with caution.
- Deliver full feature table, significant features, volcano, heatmap and PCA.

## Pathway and Multi-omics

- Proteins can hand off to `enrichment-analysis` or `network-pathway-analysis`.
- Metabolites need ID mapping and pathway database notes.
- Multi-omics integration can include correlation heatmaps, feature networks or mixOmics/DIABLO-style models.
