# Paper Practice

Follow MSstats/DEP/MetaboAnalyst/mixOmics-style omics papers: missing values, normalization, and platform details must be visible in the workflow.

## Reference Sources

- MSstats GitHub: https://github.com/Vitek-Lab/MSstats
- MSstats Bioconductor: https://bioconductor.org/packages/MSstats/
- DEP Bioconductor: https://bioconductor.org/packages/DEP/
- MetaboAnalyst: https://www.metaboanalyst.ca/
- mixOmics: https://mixomics.org/

## Paper-style Flow

1. Confirm feature matrix, metadata, platform, normalization status, missing value rule, and contrast.
2. Report missingness and sample QC before differential testing.
3. Apply appropriate normalization/imputation or state that input is already normalized.
4. Generate PCA, differential features, volcano, top feature heatmap, and pathway/multi-omics handoff.
5. Document platform, preprocessing, imputation, statistical model, thresholds, and databases.

## Visualization Pattern

- PCA should emphasize groups; avoid labeling every sample unless count is small.
- Volcano labels should be top features only and staggered.
- Heatmaps should cap feature count and use a centered color scale.
