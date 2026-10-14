---
name: public-db-survival
description: Portable public database expression and survival analysis workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides TCGA, GEO, public cohort, expression-clinical table, gene expression columns, phenotype groups, survival time/status, clinical covariates, validation cohorts, or risk scores and asks for expression comparison, Kaplan-Meier curve, log-rank test, univariate or multivariate Cox model, forest plot, time-dependent ROC-style exploratory plot, stratified survival, public database validation, or figure deliverables.
---

# Public Database Survival

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Validate sample IDs, expression columns, clinical grouping, survival time, and event status.
2. Confirm whether expression should be split by median, tertile, user cutoff, or a provided group.
3. Use `scripts/run_survival.R` for formal survival and Cox outputs when R packages are available.
4. Read `references/advanced-modules.md` when covariates, validation cohorts, stratification, Cox forest plots, or risk scores are requested.
5. Use `scripts/quickplot_survival.py` for KM and expression PNG/PDF previews.
6. State data source, cohort, normalization status, censoring definition, and sample count in every report.
7. Do not claim clinical significance from tiny or simulated cohorts. Mark exploratory results clearly.
8. Before finalizing conclusions, read `references/pitfalls.md` and check the result against each pitfall (cutoff selection bias, unchecked proportional hazards, immortal time bias, competing risks, lead-time bias, batch/platform cutoff transfer, multiple testing, underpowered subgroups, confounding, endpoint conflation, and more).

## Visual Standard

Read `references/visual-style.md` before final plots. KM curves must show group sizes, p value, and cutoff rule.

## Scripts

- `scripts/make_demo_data.py`: create demo expression-clinical table.
- `scripts/quickplot_survival.py`: dependency-light (matplotlib + numpy only) KM curve (with risk table, censoring ticks, HR and log-rank p), expression boxplot, and Cox forest plot (from-scratch Newton-Raphson Cox fit, Breslow tie approximation) PNG/PDF generator. Supports `--cutoff-method median|tertile|quartile` and `--covariates` for a multivariate forest plot.
- `scripts/run_survival.R`: KM curve + log-rank test + Cox proportional-hazards check (`survival::cox.zph`) workflow. Only `survival` is required; supports `median`/`tertile`/`quartile`/`maxstat` cutoff rules (maxstat requires survminer and is always labeled exploratory). Uses `survminer`/`ggplot2` for a polished `ggsurvplot` when available, otherwise falls back to base-R graphics (`plot.survfit` + a manual risk-table panel) so PNG/PDF figures are still produced with a minimal R install.
- `scripts/run_cox_forest.R`: univariate (per covariate) + multivariate Cox model, proportional-hazards check for every model, and an events-per-variable (EPV) check, when covariates are available. Uses `ggplot2` for a polished forest plot when available, otherwise falls back to a manual base-R forest plot (`plot()`/`segments()`/`points()`).
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name public-db-survival` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
