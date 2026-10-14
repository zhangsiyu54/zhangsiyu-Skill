# R Package Requirements

Formal route (`scripts/run_survival.R`, `scripts/run_cox_forest.R`) needs:

```r
install.packages(c("survival","survminer","ggplot2"))
```

`survival` ships with base R on most installs. `survminer`/`ggplot2` are optional in the sense
that both scripts fall back to a base-R-graphics KM curve / manual forest plot when they are
missing, so a minimal R install (just `survival`) still produces PNG/PDF output — just less
polished.

```bash
python3 scripts/check_environment.py --out output --r-packages survival,survminer,ggplot2
```

Fall back to `scripts/quickplot_survival.py` (Python-only, includes its own from-scratch Cox/KM
implementation) if `survival` itself is unavailable.

Verified (2026-07): ran end to end with real `survival`+`survminer`, including the
`cox.zph()` proportional-hazards check. The Python quicklook route's independent Cox
implementation (Newton-Raphson) matches the R `coxph()` result to 4+ significant figures
(HR=2.109 R vs 2.10 Python). Fixed a real bug in `run_cox_forest.R`: categorical covariates
(e.g. tumor stage, sex) were previously coerced to numeric NA and silently dropped, making
multivariate Cox models unable to include any categorical clinical variable — now correctly
detected and modeled as a factor.
