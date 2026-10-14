# Source code audit

## Public release status (2026-10-14)

The public release contains original base-R helpers and documentation. Course-derived `references/source-code/` is an optional LOCAL archive, excluded from publication because no redistribution license was supplied. It is not required for validation, EAPC, decomposition, or tests. Preserve any local archive and consult it only when the user requests work on those materials. Do not package it or derive a license for it from this repository's documentation.

The following audit describes the historical local archive, not current executable dependencies. Advanced course methods were not validated by the core helper tests. Do not claim a full GBD workflow has passed merely because those tests pass.

## What was found

The source directory contained no standard R package source: there was no `DESCRIPTION`, `NAMESPACE`, or package-level `R/` directory. It contained R scripts, R Markdown notebooks, example datasets, rendered outputs, PDFs, videos, and duplicated BAPC materials. The code was separated by analytical purpose and preserved under `references/source-code/`.

## Required modernization

1. Replace hard-coded Windows paths and `setwd()` with explicit input/output arguments and project-relative paths.
2. Remove `install.packages()` calls from analysis code. Manage dependencies before execution and record versions.
3. Replace retired `rgdal`, `maptools`, and related legacy spatial workflows with `sf`, `geodata`, `rnaturalearth`, or another maintained approach after checking geographic definitions.
4. Treat `BAPC`, `INLA`, `ggfan`, and R-Forge installation instructions as version-sensitive. Validate repository URLs and compatible versions before running.
5. Do not reuse GBD 2019 population inputs or standard weights automatically for GBD 2021 or later exports.
6. Replace embedded disease-specific filenames, measures, sexes, locations, and output names with configuration parameters.
7. Remove runtime package installation and implicit object reuse from `.RData`/`.Rhistory`.
8. Confirm that uncertainty limits remain paired with each estimate after joins and aggregation.
9. Test duplicated BAPC templates against a small fixture before using them on a full dataset.
10. Review the 2,055-line APC helper before publication use. It is preserved for provenance but has not been independently validated against an authoritative current APC implementation.

## Excluded from the Skill

- Course videos and PDFs.
- Rendered figures and intermediate `.RData`/`.Rhistory` files.
- Large disease-specific GBD and population CSV files.
- Duplicate BAPC archive data.

These exclusions keep the Skill portable. Analysts should supply current authorized GBD exports and population inputs for each project.
