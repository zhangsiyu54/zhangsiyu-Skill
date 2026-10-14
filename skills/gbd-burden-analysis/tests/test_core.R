script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "scripts", "gbd_core.R"))
checks <- 0L
check <- function(condition) { stopifnot(condition); checks <<- checks + 1L }
rejects <- function(expr) check(inherits(tryCatch(force(expr), error = identity), "error"))
d <- data.frame(measure_name = "Deaths", location_name = "Synthetic", sex_name = "Both",
  age_name = "Age-standardized", cause_name = "Synthetic disease", metric_name = "Rate",
  year = 2000:2009, val = 100 * 1.02^(0:9))
d$lower <- d$val * 0.9; d$upper <- d$val * 1.1
check(nrow(gbd_validate(d)) == 10)
rejects(gbd_validate(d[FALSE, ]))
for (bad in c(NA_real_, NaN, Inf, -Inf, -1)) {
  b <- d; b$val[1] <- bad; rejects(gbd_validate(b))
}
b <- d; b$lower[1] <- 200; rejects(gbd_validate(b))
b <- d; b$year[1] <- 2000.5; rejects(gbd_validate(b))
b <- d; b$sex_name[1] <- " "; rejects(gbd_validate(b))
b <- d; names(b)[2] <- names(b)[1]; rejects(gbd_validate(b))
rejects(gbd_validate(rbind(d, d[1, ])))
b <- d; b$location_id <- "1"; b$location_name[1] <- "Other"; rejects(gbd_validate(b))
b <- d; b$metric_id <- "3"; b$metric_name <- NULL
check(abs(gbd_eapc(b)$eapc - 2) < 1e-7)
r <- gbd_eapc(d)
check(abs(r$eapc - 2) < 1e-7 && r$status == "ok")
check(r$ci_lower <= r$eapc && r$eapc <= r$ci_upper)
check(abs(gbd_eapc(d[10:1, ])$eapc - r$eapc) < 1e-10)
check(gbd_eapc(d[-4, ])$missing_years == 1)
check(gbd_eapc(d[1:2, ])$status == "insufficient_years")
b <- d; b$val[1] <- b$lower[1] <- 0; check(gbd_eapc(b)$status == "nonpositive_rate")
b <- d; b$metric_name <- "Number"; rejects(gbd_eapc(b))
b <- d; b$location_name <- "Other"; b$val <- 50 * 0.97^(0:9); b$lower <- b$val * .9; b$upper <- b$val * 1.1
r <- gbd_eapc(rbind(d, b)); check(nrow(r) == 2 && max(abs(sort(r$eapc) - c(-3, 2))) < 1e-7)
x <- data.frame(age = c("0-49", "50+"), population_start = c(800,200), population_end = c(960,240), rate_start = c(100,1000), rate_end = c(100,1000))
r <- gbd_decompose(x)
check(abs(r$contribution[1] - .56) < 1e-10 && all(abs(r$contribution[-1]) < 1e-10))
x$population_end <- c(700,300); x$rate_end <- c(80,700)
r <- gbd_decompose(x); check(abs(sum(r$contribution) - (2.66-2.8)) < 1e-10)
reverse <- x; reverse[c("population_start","population_end","rate_start","rate_end")] <- x[c("population_end","population_start","rate_end","rate_start")]
check(max(abs(gbd_decompose(reverse)$contribution + r$contribution)) < 1e-10)
check(max(abs(gbd_decompose(x[2:1, ])$contribution - r$contribution)) < 1e-10)
b <- x; b$population_start <- 0; rejects(gbd_decompose(b))
b <- x; b$age[2] <- b$age[1]; rejects(gbd_decompose(b))
rejects(gbd_decompose(x, 0))
# End-to-end CLIs: execute from a different working directory and protect output.
fixture <- tempfile(fileext = ".csv"); utils::write.csv(d, fixture, row.names = FALSE)
out <- tempfile("gbd-cli-")
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(path, args) {
  log <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(path), shQuote(args)), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status"); if (is.null(status)) 0L else status
}
oldwd <- getwd(); setwd(tempdir())
check(run(file.path(root, "scripts/validate_gbd_input.R"), fixture) == 0)
check(run(file.path(root, "scripts/run_analysis.R"), c("eapc", fixture, out)) == 0)
check(file.exists(file.path(out, "manifest.tsv")))
check(run(file.path(root, "scripts/run_analysis.R"), c("eapc", fixture, out)) != 0)
check(run(file.path(root, "scripts/check_dependencies.R"), "core") == 0)
setwd(oldwd)
cat(sprintf("PASS: %d behavioral checks\n", checks))
