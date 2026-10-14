# R Package Requirements

Formal route (`scripts/run_deseq2.R`, `scripts/run_limma_voom.R`) needs:

```r
BiocManager::install(c("DESeq2","edgeR","limma","ggplot2","pheatmap","apeglm","ashr"))
```

Check what's actually installed before choosing the formal-R vs Python-quicklook route:

```bash
python3 scripts/check_environment.py --out output --r-packages DESeq2,edgeR,limma,ggplot2,pheatmap,apeglm,ashr
```

If any package is missing, `run_deseq2.R`/`run_limma_voom.R` will fail fast with a clear
"Missing R packages: ..." message rather than a partial/incorrect run — fall back to
`scripts/quickplot_bulk.py` (Python-only, no R packages needed) in that case.

Verified (2026-07): a real installed environment ran `run_deseq2.R` end to end (real DESeq2 +
apeglm + ggplot2 + pheatmap) with correct PCA/volcano/heatmap output. `vst()` requires the
`nsub` argument to be capped at the gene count for inputs with fewer than 1000 genes — this is
already handled in the script.
