# R Package Requirements

Formal route (`scripts/run_enrichment.R`, `scripts/run_gsea_fgsea.R`) needs:

```r
BiocManager::install(c("clusterProfiler","fgsea","GOSemSim","ggplot2","org.Hs.eg.db"))
```

`org.Hs.eg.db` is only needed for human SYMBOL/ENTREZ ID mapping (KEGG route); swap for the
relevant `org.<Xx>.eg.db` package for other organisms.

```bash
python3 scripts/check_environment.py --out output --r-packages clusterProfiler,fgsea,GOSemSim,ggplot2,org.Hs.eg.db
```

Missing packages make the script fail fast with a clear message; fall back to
`scripts/quickplot_enrichment.py`/`scripts/quickplot_gsea.py` (Python-only) in that case.

Verified (2026-07): `run_enrichment.R` ran end to end with real clusterProfiler `enrichGO` on a
real cell-cycle gene panel and correctly recovered "cell cycle process" as the top enriched term
(padj=4.99e-07). `run_gsea_fgsea.R` ran end to end with real fgsea and correctly recovered a
seeded enriched pathway.
