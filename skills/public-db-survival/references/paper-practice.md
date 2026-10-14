# Paper Practice

Follow TCGA/cBioPortal/GEPIA/Xena-style public cohort papers: cohort definition and cutoff rule are as important as the curve.

## Reference Sources

- cBioPortal: https://www.cbioportal.org/
- GEPIA2: http://gepia2.cancer-pku.cn/
- TCGAbiolinks: https://bioconductor.org/packages/TCGAbiolinks/
- UCSCXenaTools: https://cran.r-project.org/package=UCSCXenaTools

## Paper-style Flow

1. Define cohort, cancer type, data source, normalization, survival endpoint, event coding, and inclusion criteria.
2. Confirm cutoff rule: median, tertile, optimal cutoff, or user-provided grouping.
3. Generate expression comparison, KM/log-rank, and Cox model when covariates exist.
4. Add validation cohort or stratified analysis when data allow.
5. Report sample counts, censoring, HR/CI/p values, and limitations.

## Visualization Pattern

- KM plots must show group sizes, p value, endpoint, and cutoff rule.
- Forest plots should show HR, 95% CI, p value, and covariate names.
- Expression boxplots should include sample count per group.
