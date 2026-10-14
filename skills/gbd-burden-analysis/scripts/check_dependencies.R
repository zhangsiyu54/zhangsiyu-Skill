#!/usr/bin/env Rscript
groups <- list(core = character(), plots = c("ggplot2", "dplyr", "tidyr"),
  maps = c("sf", "ggplot2"), excel = c("readxl"), bapc = c("BAPC", "INLA"))
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "core"
if (length(args) > 1 || !mode %in% c(names(groups), "all")) stop("Choose core, plots, maps, excel, bapc or all")
packages <- if (mode == "all") unique(unlist(groups)) else groups[[mode]]
installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
cat(R.version.string, "\nModule: ", mode, "\n", sep = "")
if (!length(packages)) cat("Base R is sufficient for validation, EAPC and decomposition.\n")
for (i in seq_along(packages)) cat(packages[i], if (installed[i]) as.character(utils::packageVersion(packages[i])) else "MISSING", "\n")
quit(status = if (all(installed)) 0L else 2L)
