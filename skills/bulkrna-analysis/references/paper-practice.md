# Paper Practice

Follow DESeq2/limma-voom style bulk RNA-seq papers rather than a plot-only workflow.

## Reference Sources

- DESeq2 Genome Biology: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
- DESeq2 Bioconductor: https://bioconductor.org/packages/DESeq2/
- limma Bioconductor: https://bioconductor.org/packages/limma/
- fgsea: https://bioconductor.org/packages/fgsea/
- clusterProfiler: https://bioconductor.org/packages/clusterProfiler/

## Paper-style Flow

1. Confirm raw counts, design, organism, gene ID type, and contrast.
2. Match metadata to matrix columns and report mismatches before modeling.
3. For formal analysis, run DESeq2 or limma-voom; include normalization, dispersion/model design, thresholds, and multiple testing correction.
4. Include sample-level evidence: library size/QC, PCA, sample correlation or clustering when data allow.
5. Report DE table, volcano/MA, top DEG heatmap, and enrichment/GSEA handoff.
6. Write methods with package versions, design formula, thresholds, and limitations.

## Visualization Pattern

- Volcano labels should be sparse and staggered; usually top 2 up and top 2 down at most.
- Heatmaps should show top genes only, grouped samples, and clear color scale.
- Long gene names should be truncated or moved to a table.
