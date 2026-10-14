#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "gbd_core.R"))
args <- commandArgs(TRUE)
if (length(args) < 3 || length(args) > 4 || !args[1] %in% c("eapc", "decompose"))
  stop("Usage: Rscript run_analysis.R eapc|decompose input.csv new-output-dir [rate-scale=100000]", call. = FALSE)
if (args[1] == "eapc" && length(args) == 4) stop("rate-scale applies only to decompose")
if (file.exists(args[3])) stop("Output path already exists; choose a new directory")
if (args[1] == "eapc") result <- gbd_eapc(gbd_read(args[2])) else {
  x <- utils::read.csv(args[2], check.names = FALSE, stringsAsFactors = FALSE)
  result <- gbd_decompose(x, if (length(args) == 4) as.numeric(args[4]) else 100000)
}
if (!dir.create(args[3], recursive = TRUE)) stop("Could not create output directory")
utils::write.csv(result, file.path(args[3], "results.csv"), row.names = FALSE)
utils::write.table(data.frame(input = normalizePath(args[2]), md5 = unname(tools::md5sum(args[2])),
  method = args[1], rate_scale = if (args[1] == "decompose") if (length(args) == 4) args[4] else "100000" else NA,
  script_md5 = unname(tools::md5sum(file.path(dirname(script), "gbd_core.R")))),
  file.path(args[3], "manifest.tsv"), sep = "\t", row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[3], "sessionInfo.txt"))
if (args[1] == "eapc") cat("EAPC output includes status for every stratum; inspect non-ok rows.\n")
cat("Results written to ", normalizePath(args[3]), "\n", sep = "")
