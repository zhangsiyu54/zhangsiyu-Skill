# R Package Requirements

Formal route (`scripts/run_atac_templates.R`) needs:

```r
BiocManager::install(c("DESeq2","apeglm","ashr","ggplot2","pheatmap"))
```

```bash
python3 scripts/check_environment.py --out output --r-packages DESeq2,apeglm,ashr,ggplot2,pheatmap
```

Missing packages make the script fail fast with a clear message; fall back to
`scripts/quickplot_atac.py` (Python-only) in that case.

Verified (2026-07): ran end to end with real DESeq2 + apeglm + ggplot2 + pheatmap, including
PCA/volcano/peak-heatmap output.
