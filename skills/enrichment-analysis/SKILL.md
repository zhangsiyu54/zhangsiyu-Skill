---
name: enrichment-analysis
description: Portable deep enrichment and pathway analysis workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides DEG tables, up/down gene lists, ranked gene tables, gene set GMT files, enrichment result tables, organism and gene ID information, and asks for GO, KEGG, Reactome, MSigDB, Hallmark, GSEA, ORA, fgsea, clusterProfiler, ID mapping, dotplot, barplot, bubble plot, ridgeplot, enrichment map, leading-edge genes, pathway summary, or enrichment deliverables.
---

# Enrichment Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: DEG table, gene list, ranked list, or existing enrichment result table.
2. Confirm organism, gene ID type, database scope, and whether up/down genes should be analyzed separately.
3. If gene IDs and organism are clear, use `scripts/run_enrichment.R` for clusterProfiler-style ORA (GO/KEGG) with an explicit background/universe, and `scripts/run_gsea_fgsea.R` for fgsea-based GSEA on a full ranked list.
4. Read `references/advanced-modules.md` when the user requests GSEA, MSigDB, leading-edge genes, multiple databases, custom GMT, or enrichment-map style interpretation.
5. If an enrichment result table is already provided, use `scripts/quickplot_enrichment.py` for clean bubble and bar plots.
6. Never map IDs by guess. Ask for organism and ID type when missing.
7. Before finalizing conclusions, check `references/pitfalls.md` and rule out each applicable pitfall (background universe, ranking metric, GO redundancy, ID mapping loss, multiple-testing scope, leading-edge size, species/ortholog mismatch, circularity, etc.).
8. Write outputs under `output/` and state whether enrichment was newly computed or only visualized from provided results.

## Visual Standard

Read `references/visual-style.md` before plotting. Sort terms by adjusted p value, wrap long pathway names, and keep database categories visually distinct.

## Scripts

- `scripts/make_demo_data.py`: create demo DEG, ORA, ranked-list, and gene-set tables (enrichment_results.csv includes a `geneID` column for gene-concept network demos).
- `scripts/quickplot_enrichment.py`: enrichment bubble/bar PNG/PDF generator; also draws a gene-concept network (`enrichment_cnetplot`) when a `geneID` column is present and enough terms are significant. `--padj-threshold` (default 0.05) controls which terms are starred/counted as significant in the bubble/bar/cnetplot and report.md; `--top-n` (default 20) caps the bubble plot.
- `scripts/quickplot_gsea.py`: GSEA-style running score and NES summary dotplot PNG/PDF generator (long pathway names wrapped). `--nes-threshold` (default 1.0, descriptive only — no p-value in this quicklook) controls which gene sets are starred/counted as notable in the summary dotplot and report.md; `--top-n` (default 20) caps the dotplot.
- `scripts/run_enrichment.R`: real ORA pipeline (`clusterProfiler::enrichGO`/`enrichKEGG`) — requires an explicit background/universe (derived from a full DE table or supplied as a gene list; never silently defaults to whole-genome background), writes the full tested-term table plus dotplot and barplot, and optionally reduces redundant GO terms via `simplify()`. Positional arg 7 (`padj_cutoff`) sets which DE-table genes form the foreground; the separate positional arg 9 (`enrichment_padj_threshold`, default 0.05) sets which *tested terms* are flagged significant in the output table, plot titles, and report.md.
- `scripts/run_gsea_fgsea.R`: ranked-list GSEA via `fgsea::fgsea()` for custom GMT/MSigDB-style gene sets — writes the full pathway table, a separate leading-edge gene table, an NES summary dotplot, and a running-score plot for the top pathway. Positional arg 6 (`padj_threshold`, default 0.05) controls which pathways are flagged significant in the dotplot title and report_gsea.md.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name enrichment-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
