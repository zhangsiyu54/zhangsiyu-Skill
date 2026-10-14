# Dependencies

Run `Rscript scripts/check_dependencies.R <module>`; it reports versions and returns exit code 2 for missing selected dependencies. It installs nothing.

| Module | Packages |
|---|---|
| core (default) | Base R: stats, utils, tools; tested locally with R 4.4.1 |
| plots | ggplot2, dplyr, tidyr |
| maps | sf, ggplot2 |
| excel | readxl |
| bapc | BAPC, INLA |
| all | Union of optional packages above |

The public distribution's tests and CLI tools need only core. Optional module checks indicate software availability, not model validity or an implemented pipeline. Pin versions within a project using renv when optional software is installed. Verify BAPC/INLA repository instructions and platform compatibility at analysis time. Do not run package-install commands embedded in old course examples. No retired spatial packages are required by the maintained scripts.
