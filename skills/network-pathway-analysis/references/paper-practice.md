# Paper Practice

Follow WGCNA/GSVA/SCENIC-style studies: network or pathway results should be tied to sample traits and then narrowed to interpretable modules or pathways.

## Reference Sources

- WGCNA CRAN: https://cran.r-project.org/package=WGCNA
- GSVA Bioconductor: https://bioconductor.org/packages/GSVA/
- pySCENIC GitHub: https://github.com/aertslab/pySCENIC
- SCENIC Nature Methods: https://www.nature.com/articles/nmeth.4463

## Paper-style Flow

1. Confirm matrix normalization, sample count, phenotype columns, species, and gene IDs.
2. For WGCNA, filter genes, inspect sample outliers, choose soft threshold, define modules, and relate modules to traits.
3. For GSVA/ssGSEA, compute pathway scores and compare scores by group.
4. For regulatory networks, report database and pruning method.
5. Provide module/pathway tables, hub genes, enrichment, and methods.

## Visualization Pattern

- Module-trait heatmap: blue-white-red centered at zero, with correlation and p value when available.
- Pathway heatmap: cluster rows, keep labels readable, and cap pathway count.
- Network plots should show top edges/nodes only.
