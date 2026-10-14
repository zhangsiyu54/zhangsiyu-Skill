# R Package Requirements

Formal route (`scripts/run_wgcna_gsva.R`) needs:

```r
BiocManager::install(c("WGCNA","GSVA","GSEABase","ggplot2","pheatmap"))
```

```bash
python3 scripts/check_environment.py --out output --r-packages WGCNA,GSVA,GSEABase,ggplot2,pheatmap
```

Missing packages make the script fail fast with a clear message; fall back to
`scripts/quickplot_network_pathway.py` (Python-only) in that case.

Verified (2026-07): ran end to end with real WGCNA + GSVA on both a weak-signal dataset (0
modules detected -- correctly reported, not crashed) and a strong-signal dataset (4 real modules
detected, correct hub genes). The all-grey (0-module) degenerate case previously crashed the
script; this is now fixed.
