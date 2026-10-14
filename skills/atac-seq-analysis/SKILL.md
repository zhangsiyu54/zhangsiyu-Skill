---
name: atac-seq-analysis
description: Portable ATAC-seq and scATAC-seq workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides bulk ATAC peak count matrix, peak annotation, sample metadata, differential peak table, motif enrichment table, fragment files, Signac/ArchR object, scATAC clusters, gene activity scores, peak-to-gene links, or asks for ATAC QC, FRiP/TSS enrichment review, differential accessibility, peak annotation, motif enrichment, footprint-style interpretation, scATAC LSI/UMAP/clustering, gene score maps, peak heatmaps, or ATAC deliverables.
---

# ATAC-seq Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: bulk peak count matrix, differential peak table, motif table, fragment files, Signac/ArchR object, or exported scATAC tables.
2. Confirm genome build, sample groups, peak format, blacklist/filtering status, and whether counts are raw.
3. Read `references/advanced-modules.md` for QC metrics, bulk ATAC differential accessibility, motif enrichment, or scATAC analysis.
4. Use `scripts/quickplot_atac.py` for differential peak, motif and accessibility heatmap figures from tables.
5. Use `scripts/run_atac_templates.R` as a formal Signac/ArchR/DESeq2 template when dependencies and objects are available.
6. Export every figure as PNG and PDF.
7. Do not annotate peaks to genes or infer TF activity without stating database, genome and method.
8. Before finalizing conclusions, check `references/pitfalls.md` for known failure modes (peak-width/GC bias, QC-metric blind spots, promoter-vs-distal interpretation, motif background bias, peak-merging inflation, scATAC depth artifacts, batch over-correction, footprint depth limits, gene-activity-as-expression conflation, version drift, blacklist filtering, and more).

## Visual Standard

Read `references/visual-style.md` before plotting. Peak and motif figures must keep chromosome/peak labels readable.

## Scripts

- `scripts/make_demo_data.py`: create demo differential peak, motif and peak count tables.
- `scripts/quickplot_atac.py`: generate differential peak volcano, motif enrichment barplot, peak heatmap and sample QC overview (library depth, peaks detected) PNG/PDF from tables. Significance/display thresholds are CLI-configurable: `--padj-threshold` (default 0.05), `--log2fc-threshold` (default 1.0), `--top-n-heatmap` (default 30), `--top-n-motifs` (default 12); thresholds propagate into classification, volcano lines/title, and `report.md`.
- `scripts/run_atac_templates.R`: complete, runnable bulk ATAC DESeq2 differential-accessibility pipeline (peak count matrix + metadata in, differential peaks with apeglm/ashr shrinkage, PCA, sample correlation heatmap, top-peak heatmap, optional motif barplot, `report.md`/`summary.tsv` out); scATAC (Signac/ArchR) is out of scope and must be run interactively in its native object/session. Same thresholds are configurable via optional `--padj-threshold=`, `--log2fc-threshold=`, `--top-n-heatmap=`, `--top-n-motifs=` flags appended after the existing positional arguments (backward compatible; defaults 0.05/1/50/12).
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name atac-seq-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
