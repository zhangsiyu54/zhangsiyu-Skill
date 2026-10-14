# CRISPR Screen Deliverables

Minimum output structure (always produced by `scripts/quickplot_crispr.py` from a gene summary):
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md` — a results-and-conclusions summary, not a file listing: gene/hit counts, named top positive/negative hits with logFC/FDR and (when `--sgrna-counts` was supplied) a per-hit sgRNA-consistency note, only the `references/pitfalls.md` caveats actually triggered by this run's data (single-guide-driven hits, uneven library representation, low guide multiplicity, missing non-targeting controls, weak replicate correlation, uncalibrated essentiality, orthogonal-validation reminder), and a plain-language conclusion.
- `output/methods.md`
- `output/summary.tsv` (gene/hit counts, thresholds used, and sgRNA QC stats when available)
- `output/results/gene_summary.tsv` (real tab-separated copy of the input gene summary)
- `output/results/top_hits.tsv`
- `output/figures/crispr_volcano.png` and `.pdf`
- `output/figures/gene_rank_plot.png` and `.pdf`
- `output/figures/top_hits_barplot.png` and `.pdf` (top labeled positive/negative hits, horizontal bars with FDR annotations)
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/qc.json`

Additional outputs produced only when `--sgrna-counts` (an sgRNA-level raw count table) is supplied:
- `output/figures/qc_distribution.png` and `.pdf` — 4-panel QC: read count distribution per sample, sgRNAs-per-gene distribution, control replicate correlation (Pearson r), and Gini index per sample (library evenness).
- `output/figures/sgrna_consistency.png` and `.pdf` — per-guide log2FC for the top hit genes, to check that a gene-level call is not driven by a single sgRNA (see `references/pitfalls.md` #6).
- `output/results/sgrna_summary.tsv` — per-gene guide count, mean/sd/min/max per-guide log2FC, `n_concordant_guides` (guides matching the gene's direction with non-trivial magnitude), and `single_guide_driven` (True when a >=2-guide gene has <=1 concordant guide; see `references/pitfalls.md` #6).
- `output/run/qc.json` additionally includes a `sgrna_counts` block (control/treatment column split, Gini index per sample, mean guides per gene, replicate Pearson r).

`scripts/make_demo_data.py` writes both a gene summary (`gene_summary.csv`) and an sgRNA-level count table (`sgrna_counts.csv`, control/treatment x2 replicates) so the full figure set above can be exercised end-to-end; pass `--no-sgrna` to skip the count table.
