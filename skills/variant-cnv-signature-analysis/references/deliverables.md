# Variant CNV Signature Deliverables

Minimum output structure (Python quickplot route, `scripts/quickplot_cancer_genomics.py`):
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/mutation_frequency.tsv`
- `output/results/cnv_matrix.tsv` when available
- `output/results/signature_exposure.tsv` when available
- `output/results/tmb_summary.tsv` when a `--tmb` table is supplied
- `output/figures/mutation_frequency_barplot.png` and `.pdf`
- `output/figures/cnv_heatmap.png` and `.pdf`
- `output/figures/signature_exposure.png` and `.pdf`
- `output/figures/tmb_summary_boxplot.png` and `.pdf` when a `--tmb` table is supplied
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/qc.json`

Additional outputs when the formal R/maftools route is used (`scripts/run_maftools_signature.R`, requires maftools + NMF + a BSgenome data package; see script header for exact requirements per step):
- `output/results/sample_summary.tsv`, `output/results/gene_summary.tsv`
- `output/results/tmb_summary.tsv` (from `maftools::tmb()`, with the assay-specific `capture_size_mb` denominator stated in `output/report.md`)
- `output/results/trinucleotide_matrix.csv`, `output/results/apobec_scores.tsv` when a BSgenome package is available
- `output/results/signature_exposure.csv`, `output/results/signature_cosmic_similarity.csv` when the cohort is large enough for de novo extraction (>=20 samples); otherwise `output/results/signature_analysis_skipped.txt` explaining why, per `references/pitfalls.md` #2/#3
- `output/figures/oncoplot.png`/`.pdf`, `output/figures/maf_summary_dashboard.png`/`.pdf`, `output/figures/tmb_summary_boxplot.png`/`.pdf`, `output/figures/signature_cophenetic.png`/`.pdf` and `output/figures/signature_exposure.png`/`.pdf` when produced
- `output/report.md` (includes an explicit limitations section: no purity/ploidy-aware CNV, tumor-only germline caveat, de novo vs refit distinction)
