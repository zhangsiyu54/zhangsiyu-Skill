# CNS Method Guide

Use this guide before running public cohort expression, clinical association, and survival analysis.

## Minimal Inputs

- Require cohort/source, expression matrix, clinical table, sample IDs, survival time, event status, endpoint definition, target gene/risk score, and covariates when adjusted models are requested.
- Confirm units for survival time and event coding. Convert consistently before modeling.
- Keep primary tumor, normal, metastatic, recurrence, or treatment subsets explicit.

## Cohort and Expression QC

- Report cohort name, data portal, download date, normalization type, gene ID type, and sample count after each filter.
- Remove duplicated samples only by transparent rule, such as primary tumor preference or highest-quality aliquot.
- Check missing clinical covariates and censoring distribution.
- For public validation, keep discovery and validation cohorts separate.

## Survival Modeling

- Kaplan-Meier plots require a pre-specified or transparently chosen cutoff. Prefer median/tertile or biologically pre-defined cutoff for simple delivery.
- Maximally selected cutoffs must be labeled exploratory and adjusted/validated due to optimistic bias.
- Use log-rank test for KM group comparison.
- Use Cox proportional hazards for effect estimates; include univariate and multivariate models when covariates are present.
- Check proportional hazards assumptions when possible and state if not assessed.

## Expression and Clinical Association

- Use Wilcoxon/t-test or ANOVA/Kruskal-Wallis based on groups and distribution; report sample counts per group.
- Correct multiple testing when testing many genes/features.
- Avoid causal claims from public retrospective cohorts.

## Required Figures

- Expression box/violin plot with sample counts.
- KM curve with risk table, hazard ratio when available, confidence interval, log-rank P value, and censoring ticks.
- Cox forest plot for covariates.
- Time-dependent ROC/calibration only when model and follow-up structure support it; label exploratory if using simple scripts.

## Methods Text Must Include

- Database/cohort, download date, inclusion/exclusion filters, normalization, endpoint/event definition, cutoff rule, tests, Cox covariates, software packages, and limitations.

## Algorithm and Design Decision Tree

Read `references/pitfalls.md` alongside this table before finalizing conclusions; each row below maps to a pitfall in that file.

| Scenario | Tool/method/parameter | Rationale |
|---|---|---|
| Choosing an expression cutoff for KM | Default to a pre-specified rule: median split, tertile, or a literature-established threshold. Only if the user explicitly wants a data-driven cutoff, use maximally selected rank statistics (`maxstat`/`survminer::surv_cutpoint`, Hothorn & Lausen 2003) and report the correction-adjusted p value, not the raw minimum p. | A cutoff chosen to maximize KM separation is an implicit multiple-testing search; the raw log-rank p at the winning cutoff is optimistically biased unless labeled exploratory and validated independently. |
| Non-proportional hazards detected | Run `survival::cox.zph` (Schoenfeld residuals, Grambsch & Therneau 1994) first. If the global or covariate-specific test is significant, switch to a time-varying coefficient Cox (`tt()` term), a piecewise/interval-split Cox, or report restricted mean survival time (RMST) instead of a single HR. | A constant-HR summary is misleading when the true effect changes sign or magnitude over follow-up; time-varying models or RMST describe the effect honestly. |
| Competing risks present (e.g., death from other causes in disease-specific survival) | Fit a cause-specific Cox model for the etiologic question and a Fine-Gray subdistribution hazard model (Fine & Gray 1999; `cmprsk`/`tidycmprsk` in R) for the cumulative-incidence/prediction question. Do not use naive KM/Cox that censors competing events as if non-informative. | Cause-specific and subdistribution hazards answer different clinical questions and can disagree in direction; naive censoring of competing events overestimates the cumulative incidence of the event of interest. |
| Small number of events for multivariate Cox | Apply the events-per-variable (EPV) heuristic (~10 events per covariate, Peduzzi et al. 1996; can be relaxed toward ~5-10 per Vittinghoff & McCulloch 2007 depending on context). If EPV is too low, reduce covariates to clinically pre-specified ones, or use penalized/ridge Cox (`glmnet` Cox path) instead of unrestricted multivariate Cox. | Below roughly 10 events per variable, Cox coefficients become biased, variance is misestimated, and confidence intervals lose nominal coverage. |
| Cross-cohort validation of a prognostic signature | Freeze the cutoff rule and/or model coefficients from the discovery cohort; apply the identical rule (not the identical numeric value) to the validation cohort after re-deriving distributional summaries (e.g., "median of this cohort") within that cohort. Never re-tune on the validation set. | Absolute cutoffs and coefficients do not transfer across platforms/normalizations (see pitfalls #1 and #6); tuning on the validation set defeats the purpose of validation. |
| Time-dependent covariate/exposure (e.g., treatment started after diagnosis, biomarker response status) | Use a counting-process (start, stop, event) Cox model with the covariate coded as time-dependent, or a landmark analysis anchored at a pre-specified time point. Never assign delayed-availability information back to time zero. | Static coding of a covariate that is only knowable after a delay produces immortal time bias (Suissa 2008), inflating the apparent benefit of the exposure. |
| Testing many genes/features for prognostic significance | Apply Benjamini-Hochberg FDR or permutation-based null correction across the full number of features screened; report q values next to raw p values; require independent-cohort replication before calling a multi-gene panel "prognostic." | With thousands of genes tested, a nominal p < 0.05 threshold guarantees numerous false positives by chance, a classic driver of irreproducible expression signatures. |
| Batch/platform difference between discovery and validation cohort (e.g., microarray vs. RNA-seq) | Re-derive the cutoff within each cohort using the same rule rather than transferring the absolute numeric value; or rank-transform/z-score standardize expression within each platform before applying any threshold. | Absolute expression values are not comparable across platforms or normalization pipelines; naive value transfer produces grossly unbalanced or non-replicating groups. |
| Endpoint choice ambiguous (OS vs. DSS vs. PFS available) | Pre-specify one primary endpoint matching the clinical question before running any model; label every figure/table with the exact endpoint; report other endpoints as clearly secondary. | Switching endpoints mid-analysis (endpoint-shopping) can surface a favorable p value that does not reflect a pre-specified hypothesis, and conflates distinct event definitions. |
| Known strong clinical confounders available (age, stage, treatment) | Always pair a univariate Cox/log-rank result with a multivariate Cox model adjusting for stage, age, and treatment before using language like "independent prognostic factor." | Expression is frequently correlated with stage; a univariate-only association can be entirely explained by confounding and disappear after adjustment. |
