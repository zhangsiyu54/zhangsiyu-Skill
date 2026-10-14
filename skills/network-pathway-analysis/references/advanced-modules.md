# Network Pathway Advanced Modules

## WGCNA

- Recommended for moderate/large sample cohorts, not tiny two-group studies.
- Steps: filter low-expression genes, check samples/outliers, pick soft-threshold power, build adjacency/TOM, detect modules, merge similar modules.
- Deliver module eigengenes, module-trait correlation, hub genes and module enrichment.
- State correlation method, sample count, soft power and module cut height.

## GSVA and ssGSEA

- Use normalized/log expression or appropriate omics matrix.
- Gene sets: Hallmark, KEGG, Reactome, GO, user-provided GMT.
- Deliver pathway score matrix, group comparison table, heatmap and pathway activity boxplots.

## Hub Genes

- Hub gene criteria can include module membership, gene significance and network connectivity.
- Do not call a gene mechanistic driver without external validation.

## Integration with Other Skills

- For DEG-to-pathway interpretation, combine with `enrichment-analysis`.
- For survival validation of modules or genes, combine with `public-db-survival`.
- For multi-omics network interpretation, combine with `proteomics-metabolomics-analysis`.
