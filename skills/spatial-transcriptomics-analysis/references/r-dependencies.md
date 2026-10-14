# R Package Requirements

Formal route (`scripts/run_spatial_seurat.R`) needs:

```r
BiocManager::install(c("Seurat","ggplot2"))
```

```bash
python3 scripts/check_environment.py --out output --r-packages Seurat,ggplot2
```

Missing packages make the script fail fast with a clear message; fall back to
`scripts/quickplot_spatial.py` (Python-only) in that case.

Verified (2026-07): Seurat itself was confirmed working (see `singlecell-analysis`'s
`run_seurat.R` verification); this script shares the same Seurat spatial API surface and was
reviewed but not separately re-executed.
