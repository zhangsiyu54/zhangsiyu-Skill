# R Package Requirements

Formal route (`scripts/run_omics_templates.R`) needs:

```r
BiocManager::install(c("limma","ggplot2","pheatmap"))
```

```bash
python3 scripts/check_environment.py --out output --r-packages limma,ggplot2,pheatmap
```

Missing packages make the script fail fast with a clear message; fall back to
`scripts/quickplot_omics.py` (Python-only) in that case.

Verified (2026-07): ran end to end with real limma `lmFit`/`eBayes` moderated t-test on a
synthetic 80-feature/10-sample matrix.
