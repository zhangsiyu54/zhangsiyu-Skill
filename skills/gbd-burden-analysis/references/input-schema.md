# Inputs and commands

## GBD estimate CSV (validation and EAPC)

One row per estimate. Require `year,val,lower,upper` and either `_id` or `_name` for each of `measure,location,sex,age,cause,metric`. Values and interval bounds must be finite and nonnegative; year must be an integer. Names are trimmed. Empty tables, missing dimensions, duplicated columns/estimate keys, inconsistent ID/name mappings and inverted intervals fail validation.

All `*_id`, `*_name` columns and optional `release`, `gbd_release`, `scenario` fields define strata. Remove administrative file/row IDs before analysis: they are not scientific strata. Add a release column when combining releases. Additional dimensions not following this convention must be renamed or analyzed separately. The validator cannot infer omitted scientific dimensions, source licenses, age-band overlap or release compatibility.

EAPC requires `metric_name=Rate` (case-insensitive) or, in ID-only exports, the standard GBD `metric_id=3`. Confirm ID dictionaries for the supplied release. The tool fits each stratum separately, retains short/zero-rate series with a status, and reports year gaps. Uncertainty bounds validate the source table but are not used as regression weights. At least three years are required for residual degrees of freedom; adequate scientific power may require more.

Example (synthetic; never use as observed data):

```csv
measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,lower,upper
Deaths,Synthetic,Both,Age-standardized,Synthetic disease,Rate,2000,100,90,110
Deaths,Synthetic,Both,Age-standardized,Synthetic disease,Rate,2001,102,91.8,112.2
Deaths,Synthetic,Both,Age-standardized,Synthetic disease,Rate,2002,104.04,93.636,114.444
```

## Decomposition CSV

Analyze exactly one location/sex/cause/measure pair of years per file. Columns: `age,population_start,population_end,rate_start,rate_end`. Unique age labels must refer to the same nonoverlapping age intervals at both endpoints. Supply all relevant age bands, including zero-population bands if appropriate. Never mix all-ages totals or age-standardized rates with age-specific rows.

Rates default to per 100,000. Set the final CLI argument to `1` for per-person rates. Population is the count of persons, not a rate. The decomposition reconstructs `sum(population * rate / scale)` and does not invent endpoint uncertainty. Contributions can be negative; a zero total change may still have nonzero cancelling contributions, so percentage shares are not produced.

```csv
age,population_start,population_end,rate_start,rate_end
0-49,800,960,100,100
50+,200,240,1000,1000
```

Expected count increase: 0.56, entirely from population growth. Counts here are synthetic expected values and need not be integers.

## Limitations

These helpers read the CSV into memory; they are intended for aggregate GBD exports. They do not read Excel directly, estimate joint-risk PAF, fit Joinpoint/APC/BAPC, or compute uncertainty draws. Convert Excel explicitly and validate the converted table. Package checks are module-specific; see [dependencies.md](dependencies.md).
