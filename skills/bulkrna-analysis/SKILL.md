---
name: bulkrna-analysis
description: Portable deep bulk RNA-seq analysis workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides raw count/count matrix, normalized expression matrix, sample metadata, groups, batches, covariates or contrasts, and asks for QC, PCA, sample correlation, differential expression with DESeq2/edgeR/limma-voom, multi-factor design, paired design, volcano plot, MA plot, heatmap, normalized expression, gene set enrichment handoff, GSVA/WGCNA handoff, or a reproducible bulkRNA deliverable package.
---

# bulkRNA Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Read the request, matrix, metadata, and contrast before running analysis.
2. Validate that sample IDs match exactly between matrix columns and metadata rows.
3. Confirm group column, case/control labels, organism, gene ID type, and whether the matrix is raw counts.
4. If raw counts and R/Bioconductor are available, use `scripts/run_deseq2.R` for DESeq2-style analysis.
5. Use `scripts/run_limma_voom.R` when the user has multi-factor, batch, paired, or continuous covariate designs and the inputs are raw counts.
5a. If the input is NOT raw counts (already normalized/log-scale expression -- TPM/FPKM/CPM/logCPM, a microarray-style matrix, or a matrix exported from an Excel workbook/public dataset/collaborator's pipeline) and R/Bioconductor is available, use `scripts/run_limma_trend.R` for a formal limma `eBayes(trend=TRUE)` route. Do NOT write a new ad hoc analysis script for this scenario or hand-roll a t-test -- this input pattern (non-integer, already-normalized values, sometimes from a non-CSV source) is common enough that it has its own packaged, fully-reported script; if the raw file is Excel/other, extract/convert it to the `gene` + sample-column CSV/TSV shape this script expects first, then invoke it, so the delivered report/caveats/thresholds/pptx pipeline stays consistent instead of being reinvented per request.
6. Read `references/advanced-modules.md` when the request includes batch effects, multi-group comparisons, sample QC, pathway scores, WGCNA, or downstream interpretation.
7. If dependencies are missing or the user only needs a quick preview, use `scripts/quickplot_bulk.py` and clearly label it as exploratory. For non-integer normalized/log-like expression matrices, keep `--expression-scale auto` or pass `--expression-scale log-expression` so the script does not apply a second `log2(x+1)` transform.
8. Write outputs under `output/` with `report.md`, `summary.tsv`, `results/`, `figures/`, `qc/`, and `run/`.
9. Read `references/pitfalls.md` before finalizing conclusions or writing the methods section, and check the delivered results against it.
10. Do not invent contrasts, species, ID mappings, or missing sample labels. Write `output/questions.md` when required information is absent.

## Visual Standard

Read `references/visual-style.md` before generating final plots. Keep all figures on a clean white theme, consistent blue/red DEG colors, and export every figure as both PNG and PDF with publication-friendly labels.

## Scripts

- `scripts/make_demo_data.py`: create a small demo count matrix and metadata.
- `scripts/quickplot_bulk.py`: exploratory DEG summary plus volcano, top-DEG heatmap, PCA, sample correlation heatmap, and library-size/detected-gene QC bar chart, all exported as PNG/PDF; significance/effect-size/heatmap-size thresholds are configurable via `--padj-threshold`, `--log2fc-threshold`, `--top-n-heatmap` (defaults 0.05/1.0/30). It also supports `--expression-scale auto|raw-count|log-expression|linear-expression`; `auto` uses `log2(count+1)` for integer raw counts, avoids a second log transform for small non-integer log-like expression matrices, and uses `log2(value+1)` for large continuous abundance/intensity matrices.
- `scripts/run_deseq2.R`: formal R/Bioconductor workflow when DESeq2 is available; same thresholds configurable via optional `--padj-threshold=`, `--log2fc-threshold=`, `--top-n-heatmap=` flags (defaults 0.05/1.0/50), backward compatible with the existing positional-arg invocation.
- `scripts/run_limma_voom.R`: formal multi-factor/batch-aware limma-voom workflow (raw counts only) -- `edgeR::filterByExpr`+`calcNormFactors`+`voom`+`lmFit`+`eBayes`, PCA, sample correlation heatmap, volcano, heatmap, and a full results-and-conclusions `report_limma_voom.md` (named top genes, QC verdict, data-triggered `pitfalls.md` caveats, plain-language conclusion). Supports the same optional `--padj-threshold=`/`--log2fc-threshold=`/`--top-n-heatmap=` flags as `run_deseq2.R`.
- `scripts/run_limma_trend.R`: formal limma workflow for ALREADY-NORMALIZED/log-scale expression input (not raw counts) -- `lmFit`+`contrasts.fit`+`eBayes(trend=TRUE)` (no voom, since there is no raw read-depth to weight by), same PCA/sample-correlation/volcano/heatmap/QC-summary figure set and the same report/caveats/conclusion structure as `run_deseq2.R`. Auto-detects and flags in the report when the input actually looks like raw counts and should use `run_deseq2.R`/`run_limma_voom.R` instead. Same optional `--padj-threshold=`/`--log2fc-threshold=`/`--top-n-heatmap=` flags.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name bulkrna-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md` for accepted file names, columns, and defaults. Load it whenever file names or column names are unclear.

## Deliverable Checklist

See `references/deliverables.md` before final handoff. Every referenced file must exist and be openable.
