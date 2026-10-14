---
name: variant-cnv-signature-analysis
description: Portable cancer genomics variant, CNV and mutational signature workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides MAF-like mutation tables, gene mutation frequencies, CNV segments or gene-level CNV calls, mutational signature exposures, clinical groups, or asks for oncoplot-style summaries, mutation frequency barplots, lollipop-style summaries, CNV heatmaps, signature exposure plots, tumor mutation burden summaries, or cancer genomics deliverables.
---

# Variant CNV Signature Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: MAF-like mutation table, gene frequency table, CNV matrix, signature exposure table, per-sample TMB table, or existing cancer genomics summary.
2. Confirm genome/build, sample IDs, cohort, gene symbols and grouping columns.
3. Use `scripts/check_environment.py` first (see `references/r-dependencies.md` for the exact R package list and install command). Use Python quickplot when R/maftools (+ NMF + a BSgenome data package) is not available; use `scripts/run_maftools_signature.R` when it is.
4. Use `scripts/quickplot_cancer_genomics.py` for mutation frequency, CNV heatmap, signature exposure and (with `--tmb`) TMB-by-assay figures. Use `scripts/run_maftools_signature.R` for a formal MAF-based oncoplot, MAF summary dashboard, `tmb()`-based TMB, and trinucleotide-context de novo signature extraction with COSMIC comparison.
5. Read `references/visual-style.md` before plotting. Export only PNG/PDF unless explicitly requested otherwise.
6. After every run, write `output/methods.md` with `scripts/write_methods_doc.py`.
7. Do not infer driver status from mutation frequency alone.
8. Before finalizing conclusions, check `references/pitfalls.md` (TMB panel normalization, signature overfitting at low mutation count, refit-vs-de-novo confusion, purity/ploidy-uncorrected CNV, germline contamination, and related failure modes).
9. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name variant-cnv-signature-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Scripts

- `scripts/make_demo_data.py`: create demo mutation frequency, CNV matrix, signature exposure, per-sample TMB (`tmb_summary.csv`), and MAF-route demo (`demo.maf` + `demo_clinical.tsv`) tables.
- `scripts/quickplot_cancer_genomics.py`: generate mutation frequency barplot, oncoplot-ordered CNV heatmap, signature exposure stacked bar, and (with `--tmb`) a TMB-by-assay boxplot, as PNG/PDF plus result tables.
- `scripts/run_maftools_signature.R`: formal maftools route from a MAF file -- oncoplot, MAF summary dashboard, `tmb()` TMB (assay-specific `capture_size_mb` denominator), and trinucleotide-matrix-based de novo signature extraction (`estimateSignatures`/`extractSignatures`/`compareSignatures`/`plotSignatures`) with a cohort-size gate before attempting de novo extraction. Requires maftools + NMF + a BSgenome data package (not installed/executable in the authoring sandbox; manually reviewed against the maftools API, not run -- validate once in a provisioned R environment before production use). Does not perform purity/ploidy-aware CNV calling; see the script header and `references/pitfalls.md` #4/#12.
- `scripts/check_environment.py`: inspect R/Python availability.
- `scripts/write_methods_doc.py`: write methods details.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
