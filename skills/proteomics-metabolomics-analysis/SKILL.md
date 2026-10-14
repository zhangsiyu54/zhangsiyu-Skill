---
name: proteomics-metabolomics-analysis
description: Portable proteomics, metabolomics and multi-omics workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides protein intensity matrix, LFQ/TMT table, metabolite abundance matrix, feature annotation, sample metadata, missing value table, differential feature table, pathway annotation, or asks for proteomics QC, metabolomics QC, normalization, missing value review, PCA, differential proteins/metabolites, volcano plot, heatmap, pathway enrichment, correlation network, multi-omics integration, or publication-ready omics deliverables.
---

# Proteomics Metabolomics Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: protein matrix, metabolite matrix, LFQ/TMT table, annotation table, metadata, or existing differential feature table.
2. Confirm sample IDs, group column, platform, normalization status, missing value rule, annotation IDs and contrast.
3. Read `references/advanced-modules.md` for missingness, imputation, DEP/MSstats/MetaboAnalyst-style workflows, or multi-omics integration.
4. Use `scripts/quickplot_omics.py` for PCA, volcano and heatmap from matrix/metadata/differential table.
5. Use `scripts/run_omics_templates.R` as a formal DEP/MSstats/MetaboAnalyst/mixOmics template guide.
6. Export every figure as PNG and PDF.
7. Do not hide missing-value or batch limitations; report them explicitly.
8. Before finalizing conclusions, check `references/pitfalls.md` for common methodological failure modes (MNAR imputation, TMT ratio compression, ID confidence, multi-omics linearity assumptions, etc.).

## Visual Standard

Read `references/visual-style.md` before plotting. Use consistent feature labels and avoid crowded heatmaps.

## Scripts

- `scripts/make_demo_data.py`: create demo omics matrix, metadata and differential table.
- `scripts/quickplot_omics.py`: generate PCA, volcano and heatmap PNG/PDF. Differential significance (`--padj-threshold`, default 0.05; `--log2fc-threshold`, default 1.0) and the group-aware missingness-completeness check (`--min-group-completeness`, default 0.5) are CLI-configurable; the chosen values propagate into the volcano plot lines/title and `report.md` text so a per-project threshold change is always reflected in the deliverables.
- `scripts/run_omics_templates.R`: real limma-based differential proteomics/metabolomics pipeline (group-aware missingness filtering, `lmFit`/`eBayes` moderated t-test, PCA, missingness bar, sample correlation heatmap, volcano, top-feature heatmap). Same `--padj-threshold`, `--log2fc-threshold`, and `--min-group-completeness` flags as `quickplot_omics.py` (append after the positional args), with matching defaults — the completeness filter here actually drops features before testing, not just informationally. Treats the input matrix as already-normalized log2-scale intensity per `references/input-contract.md`. DEP/MSstats/MetaboAnalyst/mixOmics routes for TMT ratio correction and multi-omics integration remain documented in `references/advanced-modules.md` for when those specific packages/objects are available.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name proteomics-metabolomics-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
