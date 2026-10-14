# Public Database Survival Input Contract

Required table columns:
- `sample`
- survival time, e.g. `time`, `survival_month`, `OS.time`
- event status, coded as 0/1 or alive/dead
- one or more gene expression columns

Optional columns:
- `group`, `type`, `stage`, `risk`, `age`, `sex`
- batch, dataset, cohort

Required parameters:
- target gene.
- cutoff rule: median, tertile, upper/lower quartile, or user-provided cutoff.
- event definition.

Stop and ask if time/status columns or target gene are missing.
