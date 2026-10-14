# Module status and selection

| Question | Implementation | Required checks |
|---|---|---|
| Is this export valid? | `scripts/validate_gbd_input.R` | Dimensions, finite values, bounds, unique estimates |
| How fast did the rate change? | `scripts/run_analysis.R eapc` | Rate selection, log-linearity, zeros, gaps, residual correlation |
| What accounts for the count change? | `scripts/run_analysis.R decompose` | Same population coverage, matching age bands, rate scale, count reconciliation |
| Tables and age/sex profiles | Agent-generated project R code | Number vs Rate, age-standard definition, interval provenance |
| Maps | Project code using sf and a documented geography crosswalk | Unmatched locations, boundaries and missing values |
| SDI association | Project-specific ecological model | Repeated country-years, time effects, weighting; correlation is not an inequality index |
| Risk-attributable burden | Summarize supplied GBD risk estimates | Risk hierarchy, overlap; distinguish supplied PAF from newly estimated PAF |
| APC | Project-specific model | Identifiability, age/period width ratio, reference choices, estimable contrasts |
| BAPC | Project-specific BAPC/INLA model | Counts/exposure alignment, priors, backtesting, posterior vs predictive uncertainty |

Optional course archives on the author's machine are educational references, excluded from public distribution. None of the maintained helpers depends on them. For historical provenance see `source-audit.md`.
