# R Package Requirements

Formal route (`scripts/run_maftools_signature.R`) needs:

```r
BiocManager::install(c("maftools","NMF"))
```

De novo mutational signature extraction additionally needs a BSgenome reference data package
(large download, ~800MB-1GB) matching your genome build:

```r
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")   # or hg38
```

Without the BSgenome package, the script still runs the oncoplot/MAF-summary/TMB parts and
cleanly skips signature extraction with an explicit message — it does not crash.

```bash
python3 scripts/check_environment.py --out output --r-packages maftools,NMF
```

Fall back to `scripts/quickplot_cancer_genomics.py` (Python-only) if maftools/NMF are
unavailable.

Verified (2026-07): ran end to end with real maftools oncoplot/MAF-summary/`tmb()` on a
synthetic 12-sample MAF; signature extraction (`trinucleotideMatrix`/`estimateSignatures`)
verified working after installing `BSgenome.Hsapiens.UCSC.hg19`.
