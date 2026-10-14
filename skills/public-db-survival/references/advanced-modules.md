# Public Database Survival Advanced Modules

Load this file for clinical/public cohort work beyond a single KM plot.

## Cohort and Preprocessing

- Record cohort/source, tumor type, sample filter, expression unit, log transform and normalization.
- Remove duplicated samples intentionally and document rules.
- Harmonize survival time units and event coding before analysis.

## Cutoff Rules

- Median cutoff is default exploratory.
- Tertile/quartile cutoffs can exaggerate separation; state excluded middle group if used.
- User-provided cutoffs must be recorded exactly.
- Do not search many cutoffs and present only the best one without correction.

## Cox Modeling

- Univariate Cox: gene/risk score only.
- Multivariate Cox: include age, stage, sex, batch or other provided covariates when appropriate.
- Deliver HR, 95% CI, p value, sample count and event count.
- Check proportional hazards if the user asks for a formal clinical-style report.

## Validation

- Keep discovery and validation cohorts separate.
- Do not tune cutoff on the validation cohort.
- Report concordant and discordant results across cohorts.

## Figures

- KM curve with group sizes, cutoff rule and p value.
- Expression box/violin across tumor-normal or clinical groups.
- Cox forest plot when covariates exist.
- Optional ROC-style exploratory plot only when time/event definitions are clear.
