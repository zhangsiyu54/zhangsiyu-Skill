# Public Database Survival Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/expression_group.tsv`
- `output/results/survival_summary.tsv` (KM/log-rank/Cox HR summary; `scripts/run_survival.R` also writes `survival_group_summary.tsv` with per-group n/events/median survival)
- `output/results/cox_summary.tsv` (univariate Cox on the target gene, written by `scripts/run_survival.R`)
- `output/results/cox_univariate.tsv` and `output/results/cox_forest.tsv` when covariates are provided (per-covariate univariate table and the joint multivariate table, written by `scripts/run_cox_forest.R`; `scripts/quickplot_survival.py` writes a combined univariate+multivariate `cox_forest.tsv` in the Python route)
- `output/results/ph_test.tsv` / `output/results/ph_test_cox_forest.tsv`: `survival::cox.zph` proportional-hazards check for every Cox model fit (see references/pitfalls.md #2); state explicitly whether the assumption was violated
- `output/figures/km_curve.png` and `output/figures/km_curve.pdf` (KM curve with risk table, censoring ticks, group sizes, cutoff rule, HR/CI and log-rank p in the title/subtitle)
- `output/figures/expression_boxplot.png` and `output/figures/expression_boxplot.pdf`
- `output/figures/cox_forest.png` and `.pdf` when Cox model is run (univariate + multivariate HRs shown together when covariates are provided)
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Notes on execution route:
- Both `scripts/run_survival.R` and `scripts/run_cox_forest.R` only require the base `survival` package. When `survminer`+`ggplot2` (and, for a nicer manual forest plot, nothing extra beyond `survival`) are unavailable, both scripts automatically fall back to base-R graphics (`plot.survfit` + a manual risk-table panel; `plot()`/`segments()`/`points()` forest plot) and still produce PNG/PDF output.
- `scripts/quickplot_survival.py` is fully dependency-light (matplotlib + numpy only) and now also produces a Cox forest plot using a from-scratch Newton-Raphson Cox fit (Breslow tie approximation); label it exploratory and prefer the R route's Efron-corrected model for a formal report.
- Every Cox model (R or Python) must be reported together with its proportional-hazards check result; never report a hazard ratio without it (references/pitfalls.md #2).

Report sections:
- cohort/source summary
- cutoff rule (median/tertile/quartile/maxstat; label maxstat exploratory per references/pitfalls.md #1)
- sample and event counts
- KM/log-rank result
- Cox result when available, including the proportional-hazards check outcome and, for multivariate models, the events-per-variable (EPV) check
- limitations
