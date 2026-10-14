# Pitfalls

Read this file before finalizing any conclusion, hazard ratio, or "prognostic" claim. Each pitfall lists the trigger that should make you suspicious, the statistical mechanism, the symptom that reveals it in the output, and the fix to apply or state as a limitation.

## 1. Post-hoc "optimal" expression cutoff (maximally selected rank statistic / X-tile-style scanning)

- **Trigger**: user or script scans many candidate cutpoints (e.g., every percentile of expression) and reports only the cutoff giving the smallest log-rank p value.
- **Mechanism**: the maximally selected rank statistic (Hothorn & Lausen, 2003) is itself a maximum over many correlated test statistics; the nominal p value from a single log-rank test at the winning cutoff does not account for the search, so it is a multiple-testing problem disguised as a single test. This is the same failure mode formalized by Lausen & Schumacher (1992) and packaged in tools like Cutoff Finder (Budczies et al., 2012).
- **Symptom**: a suspiciously small p value (e.g., p < 0.001) for a cutoff that sits far from the median/tertile, especially in a small cohort; the KM curves separate dramatically but the cutoff value looks "tuned."
- **Fix**: prefer a pre-specified cutoff (median, tertile, or literature-established threshold) for the primary result. If a data-driven optimal cutoff is used, label it explicitly as exploratory, report the corrected p value (e.g., Hothorn-Lausen or permutation-adjusted p, not the raw minimum p), and require confirmation in an independent validation cohort before any prognostic claim.

## 2. Proportional hazards assumption never checked

- **Trigger**: a single Cox hazard ratio is reported for a covariate across the entire follow-up without any diagnostic.
- **Mechanism**: the Cox model assumes the hazard ratio is constant over time. Many oncology effects are not constant (e.g., an early surgical or treatment effect that fades, or immune-related effects that strengthen over time). A single averaged HR can be an artifact that hides a crossing or attenuating effect (Grambsch & Therneau, 1994).
- **Symptom**: KM curves cross or converge over time; Schoenfeld residual plot (`cox.zph` in R) shows a non-flat trend; global test p value for `cox.zph` is significant.
- **Fix**: always run `cox.zph`/Schoenfeld residual testing when a formal Cox report is requested. If violated, report a time-varying coefficient model, split follow-up into intervals with interaction terms, use restricted mean survival time (RMST) as a robustness check, or explicitly state "the proportional hazards assumption was not met; the reported HR is a time-averaged approximation."

## 3. Immortal time bias from delayed exposure/treatment definition

- **Trigger**: treatment or biomarker-guided therapy status is defined using information that could only be known after some delay from time zero (diagnosis or surgery), then patients are classified as "treated" vs "untreated" from time zero.
- **Mechanism**: patients must survive long enough to receive the treatment, so the "treated" group accrues a period of follow-up during which, by construction, they could not have died (Suissa, 2008). This inflates the apparent survival benefit of treatment or of any exposure ascertained after baseline.
- **Symptom**: an implausibly large protective effect for a treatment or exposure, especially when the exposure requires weeks-to-months of survival to be classified as "exposed" (e.g., "received adjuvant therapy," "achieved molecular response").
- **Fix**: use a time-dependent covariate Cox model (counting-process/start-stop format) or a landmark analysis anchored at a pre-specified time point, and state the landmark time explicitly. Never assign exposure status retrospectively to time zero.

## 4. Competing risks ignored (treating competing events as simple censoring)

- **Trigger**: the endpoint is disease-specific or cause-specific (e.g., cancer death, relapse) in an elderly or comorbid cohort where death from other causes is common, but the analysis censors patients who died of other causes as if they were still at risk.
- **Mechanism**: standard KM/Cox treats competing-event patients as non-informatively censored, which overestimates the cumulative incidence of the event of interest, because those patients can no longer experience it. Cause-specific Cox and Fine-Gray subdistribution Cox (Fine & Gray, 1999) answer different questions and can point in different directions for the same covariate.
- **Symptom**: cumulative incidence from naive 1-KM exceeds what a competing-risks cumulative incidence function (CIF) would show; effect sizes for age or comorbidity-associated covariates look implausibly large on the event of interest.
- **Fix**: when a competing event is biologically plausible and non-negligible (typically >5-10% of censoring), fit both a cause-specific hazard model (etiologic question: what drives the biological process) and a Fine-Gray subdistribution hazard model (clinical/prediction question: what drives observed cumulative incidence), and report which question each answers.

## 5. Lead-time bias when comparing cohorts with different screening/diagnosis timing

- **Trigger**: comparing survival time from diagnosis across cohorts, eras, or countries where screening intensity, imaging frequency, or diagnostic criteria differ (e.g., a modern cohort diagnosed via routine screening vs. a historical cohort diagnosed at symptom onset).
- **Mechanism**: earlier detection moves the diagnosis date earlier without necessarily changing the date of death, mechanically lengthening measured survival time even if the natural history of disease and true benefit are unchanged.
- **Symptom**: a large survival improvement across eras/cohorts that is not accompanied by a corresponding improvement in stage-specific outcomes or mortality rates; cross-cohort comparisons where the "improved" cohort has a much higher proportion of early-stage or incidentally detected cases.
- **Fix**: avoid casual "survival improved over time/across cohorts" claims when screening or diagnostic practice differs. Prefer stage-matched or stage-stratified comparisons, and prefer mortality-rate-based comparisons (population-level cause-specific mortality) over survival-time comparisons when screening intensity differs.

## 6. Batch/platform differences invalidating cutoff transfer from discovery to "validation" cohort

- **Trigger**: a cutoff (absolute expression value, e.g., FPKM = 5.2) is derived in one cohort (e.g., microarray, RMA-normalized) and applied unchanged to a second cohort measured on a different platform (e.g., RNA-seq, TPM/FPKM) or with different normalization.
- **Mechanism**: microarray intensities, RNA-seq counts, and different normalization pipelines are not on the same numeric scale; absolute cutoffs are platform- and pipeline-specific and do not transfer, even when both cohorts are biologically comparable.
- **Symptom**: the validation cohort shows a wildly unbalanced high/low group split (e.g., 95%/5% instead of roughly 50/50), or the "validation" completely fails to replicate the discovery effect, and the code is silently applying a fixed numeric cutoff from the first cohort.
- **Fix**: re-derive the cutoff within each cohort using the same rule (e.g., "median of this cohort's expression distribution"), not the same absolute value; or use rank-based/quantile transformation before applying a cutoff, or z-score standardize within each cohort/platform before comparison. State normalization/platform explicitly for both cohorts in the methods text.

## 7. Multiple gene/biomarker testing for "prognostic significance" without correction

- **Trigger**: screening tens to thousands of genes for association with survival (e.g., a "prognostic signature" study scanning a full expression matrix) and reporting only genes below p < 0.05.
- **Mechanism**: with thousands of genes tested independently, a nominal 0.05 threshold guarantees many false positives by chance alone (classic multiple-testing inflation); this is a well-documented driver of irreproducible "prognostic gene signatures" in the oncology literature.
- **Symptom**: a large list of "significant" genes with no shared pathway/biological coherence; the signature fails to replicate in an independent cohort; p values cluster just under 0.05 rather than showing a strong separation from the null distribution.
- **Fix**: apply FDR correction (Benjamini-Hochberg) or permutation-based null distributions across the full number of genes/features screened, report q values alongside raw p values, and require independent-cohort validation of any multi-gene signature before calling it prognostic. State the total number of genes/features tested in the methods text.

## 8. Small subgroup survival analyses with very few events reported with unstable, over-precise HR/CI

- **Trigger**: a subgroup (e.g., stage IV, a rare histology, an uncommon mutation) has fewer than roughly 10-20 events, but a Cox HR with a narrow-looking or point-estimate-heavy interpretation is reported anyway.
- **Mechanism**: with very few events, the Cox partial likelihood is poorly identified; small changes in one or two patients' outcomes swing the HR substantially, and asymptotic confidence intervals are unreliable (related to the events-per-variable problem, Peduzzi et al., 1996).
- **Symptom**: an extremely wide 95% CI (e.g., HR 0.05-15) reported alongside a p value treated as if it were precise; or, conversely, a suspiciously tight CI given the tiny event count that likely reflects a software/quasi-separation artifact.
- **Fix**: report the exact event count (not just sample size) for every subgroup Cox model; flag any subgroup with fewer than ~10 events per estimated parameter as underpowered/exploratory; prefer Firth's penalized likelihood or exact methods for very small event counts; avoid multivariate adjustment within tiny subgroups (violates events-per-variable rule, see below).

## 9. Univariate-only "independent prognostic factor" claims that ignore known clinical confounders

- **Trigger**: an expression-survival association is reported as "an independent prognostic factor" based on a univariate Cox/log-rank result alone, with no adjustment for stage, age, sex, or treatment.
- **Mechanism**: strong clinical variables (tumor stage in particular) are usually the dominant survival predictors and are frequently correlated with expression of many genes (e.g., higher-stage tumors have systematically different transcriptomes). A univariate association can be entirely explained by stage confounding and vanish after adjustment.
- **Symptom**: a gene's univariate HR is large and significant, but the gene is also strongly stage-associated in the boxplot/ANOVA; the word "independent" appears in the report without a corresponding multivariate Cox table.
- **Fix**: never use the word "independent" without a multivariate Cox model adjusting at minimum for stage, age, and treatment when those variables are available; if covariates are unavailable, explicitly state "univariate only; independence from clinical stage/age/treatment not assessed."

## 10. Conflating overall survival (OS), disease-specific survival (DSS), and progression-free survival (PFS) within one report

- **Trigger**: different figures/tables in the same report use different endpoints (e.g., KM curve built on OS, Cox table built on PFS) without the endpoint being clearly labeled on each output.
- **Mechanism**: OS, DSS, and PFS have different event definitions, different competing-risk structure, and different censoring patterns; results are not interchangeable, and switching endpoints mid-analysis (endpoint-shopping) can be used, intentionally or not, to present the most favorable p value.
- **Symptom**: figure captions or the summary table do not state which endpoint was used; a "significant" PFS result is described in prose as if it were an overall survival benefit.
- **Fix**: pre-specify one primary endpoint matching the clinical question before analysis; label every figure, table, and sentence with the exact endpoint used; if multiple endpoints are reported, present all of them (not just the significant one) and state which was primary vs. secondary/exploratory.

## 11. Informative censoring treated as non-informative

- **Trigger**: patients are lost to follow-up, withdraw, or are censored administratively at very different rates across comparison groups (e.g., the "high expression" group has systematically shorter follow-up or higher dropout).
- **Mechanism**: KM and Cox both assume censoring is non-informative (independent of the risk of the event given covariates). If sicker patients drop out earlier for reasons related to their prognosis (e.g., transferred to hospice, moved to unlinked registry), censoring is informative and the KM estimator is biased.
- **Symptom**: markedly different censoring rates or median follow-up time between groups; a KM curve for one group that ends abruptly with many early censoring ticks clustered together.
- **Fix**: report censoring distribution per group (not just overall), compare follow-up time between groups, and flag informative censoring as a limitation when censoring patterns differ substantially by group; consider inverse-probability-of-censoring weighting if the imbalance is severe and the data support it.

## 12. Multi-gene risk score built and "validated" on overlapping or leaked data

- **Trigger**: a LASSO/Cox-based multi-gene risk score is trained (feature selection + coefficient fitting) on the full dataset, and then the same data (or a random split of it) is used to report the "validation" performance (C-index, time-dependent AUC).
- **Mechanism**: feature selection performed on the full cohort before any train/test split leaks information into the "test" set (selection bias / optimism), inflating apparent discrimination; this is distinct from, but related to, standard overfitting.
- **Symptom**: near-perfect discrimination (C-index > 0.85 in a small cohort) that is unusually strong compared to typical published prognostic signatures (usually C-index 0.6-0.75 in solid tumors); no genuinely independent (different cohort, different platform) validation performed.
- **Fix**: perform feature selection only within training folds (nested cross-validation) or, preferably, freeze the signature (genes + coefficients) from a discovery cohort and evaluate it unchanged in a separate, independent public cohort. Report both internal (cross-validated) and external validation performance separately.

## 13. Tumor purity/stromal content confounding expression-survival correlation

- **Trigger**: an immune- or stroma-related gene (e.g., an immune checkpoint, a fibroblast marker, an extracellular matrix gene) shows a strong survival association in bulk RNA-seq/microarray data.
- **Mechanism**: bulk tumor expression is a mixture of malignant, immune, and stromal cell signals; low tumor purity (more stroma/immune infiltrate) is itself often prognostic (usually favorably, sometimes unfavorably depending on tumor type), so a gene's apparent prognostic value may simply be a proxy for tumor purity rather than a tumor-cell-intrinsic effect.
- **Symptom**: the "prognostic" gene correlates strongly with established immune/stromal signature scores or with histological purity estimates (e.g., ESTIMATE, ABSOLUTE) rather than showing a purity-independent effect.
- **Fix**: when biologically plausible, check correlation with a tumor purity estimate before claiming a cell-intrinsic prognostic effect; consider adjusting for purity or explicitly note that the observed association may reflect microenvironment composition rather than tumor-cell biology.

## 14. Duplicate or overlapping samples across "independent" cohorts

- **Trigger**: combining or comparing multiple public datasets (e.g., a GEO series and a TCGA cohort, or two GEO series from the same consortium) without checking for shared patients.
- **Mechanism**: some GEO series reuse samples from the same underlying patient collection (re-deposited, reprocessed, or partially overlapping cohorts), or a "validation" cohort is actually a subset of the discovery cohort. Treating them as independent double-counts patients and produces spuriously concordant "replication."
- **Symptom**: unexpectedly high overlap in extreme expression values or clinical annotations between "independent" cohorts; identical rare covariate combinations (e.g., same age + stage + exact survival time) appearing in both.
- **Fix**: check sample/patient identifiers, submission dates, and platform/GPL accessions across cohorts before treating them as independent; when in doubt, note the possibility of cohort overlap as a limitation.

## 15. Median survival or HR reported for a curve where the median is not reached

- **Trigger**: a KM curve for a good-prognosis group never crosses the 50% survival probability line within the observed follow-up, but a "median survival" number or a precise-looking HR is reported anyway.
- **Mechanism**: when the median is not reached, any point estimate is an extrapolation beyond the data and is highly sensitive to the tail behavior of a handful of long-follow-up patients; presenting it as a firm number overstates precision.
- **Symptom**: a report states "median survival = NA" in the underlying summary table but the prose still asserts a specific number, or the KM curve plateau extends flat with wide/undefined confidence bounds at the point quoted.
- **Fix**: report "not reached" explicitly when applicable, give the follow-up duration at last observation, and prefer reporting a fixed-time survival probability (e.g., 3-year or 5-year survival with CI) instead of an unreached median when comparing groups.
