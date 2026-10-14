#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "gbd_core.R"))
args <- commandArgs(TRUE)
if (length(args) != 1) stop("Usage: Rscript validate_gbd_input.R input.csv", call. = FALSE)
x <- gbd_read(args[1])
cat(sprintf("PASS: %d rows; years %d-%d\n", nrow(x), min(x$year), max(x$year)))
