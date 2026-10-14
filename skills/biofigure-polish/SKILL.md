---
name: biofigure-polish
description: Portable bioinformatics figure polishing workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user has generated result tables, PNG/PDF plots, manuscript figures, or dashboard-like scientific panels and asks for unified color palette, typography, axis labels, panel layout, legend cleanup, journal-style aesthetics, figure manifest, or publication-ready bioinformatics visualization packaging.
---

# Biofigure Polish

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Inventory all source figures and result tables before editing.
2. Read `references/visual-style.md` and select one consistent palette for the full figure set.
3. Use `scripts/style_manifest.py` to create a figure manifest and naming plan.
4. Use `scripts/polish_figures.py` for lightweight PNG/PDF figure framing when source PNG files are available.
5. Preserve scientific meaning: do not change data values, thresholds, p values, or group labels.
6. Before finalizing conclusions, check the figure set against `references/pitfalls.md` (colorblind-unsafe palettes, rainbow colormaps, inconsistent colors/labels, font size, resolution/format, overplotting, truncated axes, and related failure modes).
7. Export polished files under `output/figures/` and write `output/figure_manifest.tsv`, `output/figure_qc_report.tsv`, and `output/report.md`. Read `output/report.md` and relay its per-figure flags and plain-language conclusion to the user — do not just report "polished N files".

## Visual Standard

Read `references/visual-style.md` before any figure changes. Avoid crowded legends, inconsistent fonts, random color reuse, low-contrast labels, and unaligned panels.

## Scripts

- `scripts/style_manifest.py`: scan figure files (keyed by path relative to the figure folder, so same-named figures in different subfolders are never merged) and create a naming/style manifest.
- `scripts/polish_figures.py`: add consistent panel labels and export polished PNG (raster preview) + PDF (vector, byte-identical to source when the source is already a PDF) figure files; loads colors from `assets/palette.json` (single source of truth); when run with `--manifest`, validates panel-label consistency, skips (rather than crashing on) missing/corrupt/unsupported source files while reporting each one clearly, and writes both `output/figure_qc_report.tsv` (machine-readable flags for low print resolution, suspected rainbow/jet colormaps, and suspected red-green diverging color use — heuristic, pixel-based only, see the script's module docstring for what it cannot detect or fix) and `output/report.md` (plain-language batch summary: counts processed/flagged/failed, each flagged figure named explicitly with its recommended fix, and a one-line conclusion). Always read `output/report.md` before summarizing results to the user. Run `--check-palette` to validate `assets/palette.json` for colorblind-safety before final delivery.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` with relevant R/Python package requirements. Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name biofigure-polish` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
