# CRISPR Screen Input Contract

Accepted inputs:
- gene summary table with `gene`, `logFC`, `pvalue`, `fdr` or `padj`.
- MAGeCK-style gene summary.
- sgRNA count table plus sample design for formal analysis.

Required parameters:
- contrast and selection direction.
- FDR threshold.
- gene identifier style.

Stop and ask if contrast or screen design is missing.
