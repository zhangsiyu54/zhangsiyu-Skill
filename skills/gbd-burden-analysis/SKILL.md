---
name: gbd-burden-analysis
description: Validate IHME GBD burden exports, estimate stratified EAPC trends, and decompose changes in counts into population, age structure, and age-specific rates. Guide GBD tables, maps, risk summaries, SDI, APC and BAPC analyses with explicit input and uncertainty requirements.
---

# GBD burden analysis

Use the executable base-R tools for validation, EAPC and count decomposition. Read [module-index.md](references/module-index.md) for other analyses: maps, SDI, attributable risks, APC and BAPC require project-specific implementation and validation. Do not describe those modules as tested executable pipelines.

## Establish the analysis

Record GBD release, cause, measure, metric, geography, sex, age range, years, risk definition and population standard. Resolve ambiguous choices before computing results. Keep each release and estimand distinguishable. Read [input-schema.md](references/input-schema.md) for supported schemas and command examples.

Resolve scripts relative to this Skill's directory; input/output paths belong to the user's project. Never source course notebooks automatically. The optional local archive is not needed by the maintained tools and is excluded from the public distribution; see [source-audit.md](references/source-audit.md).

## Run the maintained tools

```sh
Rscript scripts/check_dependencies.R core
Rscript scripts/validate_gbd_input.R /path/to/gbd.csv
Rscript scripts/run_analysis.R eapc /path/to/rates.csv /path/to/new-eapc-results
Rscript scripts/run_analysis.R decompose /path/to/age-input.csv /path/to/new-decomposition-results 100000
Rscript tests/test_core.R
```

The first three analysis helpers need only base R. Output directories must be new. Results include `results.csv`, `manifest.tsv` and `sessionInfo.txt`. Check EAPC `status` for every stratum; zero rates and short series produce explicit failed-stratum records rather than silently disappearing. Inspect year gaps, log-linear fit and residual autocorrelation before inference.

## Interpret correctly

- EAPC is `100 * (exp(beta_year) - 1)` from log-rate regression. Its t-based 95% CI is conditional on the regression assumptions; it does not propagate GBD uncertainty draws. Keep age-specific, crude and age-standardized estimates separate.
- Decomposition averages all six replacement orders (Shapley). Use matching, mutually exclusive age bands covering the specified population and age-specific rates with an explicit scale. Verify reconstructed counts against the source. The rate component is a mathematical remainder for rate change, not proof of a causal risk factor.
- Do not sum marginal GBD lower/upper bounds to obtain a combined uncertainty interval. Use aligned draws where available; otherwise report point estimates and the missing uncertainty propagation.
- Do not sum overlapping risk-attributable burdens or interpret SDI associations as individual-level causal effects.
- APC needs documented identifiability constraints and compatible age/period intervals. BAPC needs matched population forecasts, priors, holdout validation and predictive intervals. Estimated GBD counts are not raw observed event counts; document any Poisson likelihood approximation. Do not apply a count model directly to rates or DALYs.

Record exclusions, model assumptions, input provenance, package versions and unresolved limitations. For figures use vector PDF and PNG when requested. Keep generated datasets, fitted models and user data outside the installed Skill.
