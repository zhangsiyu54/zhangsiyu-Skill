# R Package Requirements

Formal route (`scripts/run_phyloseq_ancombc.R`) needs:

```r
BiocManager::install(c("phyloseq","ANCOMBC","mia","microbiome","vegan","ggplot2"))
```

`microbiome` is easy to miss — `ANCOMBC::ancombc2()` fails at runtime (not at the
`requireNamespace` check, since it's a transitive dependency) without it, with the error
"The 'microbiome' package is needed to process the imported data but is not installed."

```bash
python3 scripts/check_environment.py --out output --r-packages phyloseq,ANCOMBC,mia,microbiome,vegan,ggplot2
```

Missing packages make the script fail fast with a clear message (or, for the transitive
`microbiome` dependency, ANCOM-BC2 fails and the script skips just that section with a clear
note); fall back to `scripts/quickplot_microbiome.py` (Python-only) for anything else.

Verified (2026-07): ran end to end with real phyloseq + ANCOMBC2 differential abundance
analysis (12-sample, 2-group synthetic dataset).
