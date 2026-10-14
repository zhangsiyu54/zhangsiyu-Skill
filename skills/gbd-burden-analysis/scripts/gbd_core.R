# Base-R helpers. Sourcing this file performs no I/O or installation.
gbd_fail <- function(...) stop(..., call. = FALSE)

gbd_validate <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) gbd_fail("Input must contain estimate rows")
  names(x) <- sub("^\ufeff", "", trimws(names(x)))
  if (anyDuplicated(names(x)) || any(!nzchar(names(x)))) gbd_fail("Duplicate or empty column names")
  dimensions <- c("measure", "location", "sex", "age", "cause", "metric")
  keys <- character()
  for (d in dimensions) {
    available <- intersect(paste0(d, c("_id", "_name")), names(x))
    if (!length(available)) gbd_fail("Missing dimension: ", d, " (_id or _name required)")
    keys <- c(keys, available)
  }
  extra <- grep("(_id|_name)$", names(x), value = TRUE)
  keys <- unique(c(keys, extra, intersect(c("release", "gbd_release", "scenario"), names(x))))
  for (nm in keys) {
    z <- trimws(as.character(x[[nm]]))
    if (anyNA(z) || any(!nzchar(z))) gbd_fail("Missing dimension value: ", nm)
    x[[nm]] <- z
  }
  for (nm in c("year", "val", "lower", "upper")) {
    if (!nm %in% names(x)) gbd_fail("Missing column: ", nm)
    z <- suppressWarnings(as.numeric(as.character(x[[nm]])))
    if (any(!is.finite(z))) gbd_fail("Missing or non-finite numeric value: ", nm)
    x[[nm]] <- z
  }
  if (any(x$year != floor(x$year))) gbd_fail("Year must be an integer")
  if (any(x[c("val", "lower", "upper")] < 0)) gbd_fail("Negative burden estimate or bound")
  if (any(x$lower > x$val | x$val > x$upper)) gbd_fail("Expected lower <= val <= upper")
  releases <- intersect(c("release", "gbd_release"), names(x))
  for (d in unique(sub("_(id|name)$", "", extra))) {
    pair <- paste0(d, c("_id", "_name"))
    if (all(pair %in% names(x))) {
      mapping <- unique(x[c(releases, pair)])
      if (anyDuplicated(mapping[c(releases, pair[1])]) ||
          anyDuplicated(mapping[c(releases, pair[2])])) gbd_fail("Inconsistent ID/name mapping: ", d)
    }
  }
  keys <- c(keys, "year")
  if (anyDuplicated(x[keys])) gbd_fail("Duplicate estimate dimensions; split files or resolve duplicates")
  attr(x, "gbd_keys") <- keys
  x
}

gbd_read <- function(path) {
  if (!file.exists(path)) gbd_fail("Input not found: ", path)
  gbd_validate(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
}

gbd_eapc <- function(x) {
  x <- gbd_validate(x)
  is_rate <- if ("metric_name" %in% names(x)) tolower(x$metric_name) == "rate" else x$metric_id == "3"
  if (!all(is_rate)) gbd_fail("EAPC input must contain Rate only; explicitly filter Number/Percent")
  keys <- setdiff(attr(x, "gbd_keys"), "year")
  ordering <- do.call(order, x[keys])
  sorted_keys <- x[ordering, keys, drop = FALSE]
  changed <- rep(FALSE, nrow(x)); changed[1] <- TRUE
  if (nrow(x) > 1) for (key in keys)
    changed[-1] <- changed[-1] | sorted_keys[[key]][-1] != sorted_keys[[key]][-nrow(x)]
  groups <- integer(nrow(x)); groups[ordering] <- cumsum(changed)
  out <- lapply(split(seq_len(nrow(x)), groups), function(idx) {
    d <- x[idx, , drop = FALSE]; d <- d[order(d$year), , drop = FALSE]
    row <- d[1, keys, drop = FALSE]
    row$year_start <- min(d$year); row$year_end <- max(d$year); row$n_years <- nrow(d)
    row$missing_years <- max(d$year) - min(d$year) + 1 - nrow(d)
    row$eapc <- row$ci_lower <- row$ci_upper <- row$r_squared <- NA_real_
    row$status <- "ok"
    if (nrow(d) < 3) row$status <- "insufficient_years"
    else if (any(d$val <= 0)) row$status <- "nonpositive_rate"
    else {
      fit <- stats::lm(log(val) ~ year, data = d)
      ci <- suppressWarnings(stats::confint(fit, "year", level = 0.95))
      row$eapc <- 100 * expm1(stats::coef(fit)[["year"]])
      row$ci_lower <- 100 * expm1(ci[1]); row$ci_upper <- 100 * expm1(ci[2])
      row$r_squared <- suppressWarnings(summary(fit)$r.squared)
      if (row$missing_years > 0) row$status <- "ok_with_year_gaps"
    }
    row
  })
  result <- do.call(rbind, out); rownames(result) <- NULL
  result
}

# Shapley decomposition of count = population * sum(age_share * age_rate).
# Average six replacement orders; point estimates only.
gbd_decompose <- function(x, rate_scale = 100000) {
  required <- c("age", "population_start", "population_end", "rate_start", "rate_end")
  if (!is.data.frame(x) || !nrow(x) || !all(required %in% names(x))) gbd_fail("Decomposition requires: ", paste(required, collapse = ", "))
  if (anyNA(x$age) || any(!nzchar(trimws(as.character(x$age)))) || anyDuplicated(x$age)) gbd_fail("Age bands must be nonempty and unique")
  for (nm in required[-1]) {
    x[[nm]] <- suppressWarnings(as.numeric(as.character(x[[nm]])))
    if (any(!is.finite(x[[nm]]) | x[[nm]] < 0)) gbd_fail("Invalid nonnegative numeric column: ", nm)
  }
  if (length(rate_scale) != 1 || !is.finite(rate_scale) || rate_scale <= 0) gbd_fail("rate_scale must be positive")
  n0 <- sum(x$population_start); n1 <- sum(x$population_end)
  if (n0 <= 0 || n1 <= 0) gbd_fail("Population totals must be positive")
  a <- list(n0, x$population_start / n0, x$rate_start / rate_scale)
  b <- list(n1, x$population_end / n1, x$rate_end / rate_scale)
  count <- function(v) v[[1]] * sum(v[[2]] * v[[3]])
  orders <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  contributions <- numeric(3)
  for (order in orders) {
    current <- a
    for (j in order) {
      before <- count(current); current[[j]] <- b[[j]]
      contributions[j] <- contributions[j] + (count(current) - before) / 6
    }
  }
  delta <- count(b) - count(a)
  if (abs(sum(contributions) - delta) > 1e-8 * max(1, abs(delta))) gbd_fail("Decomposition reconciliation failed")
  data.frame(component = c("population_growth", "age_structure", "age_specific_rate"),
             contribution = contributions, count_start = count(a), count_end = count(b),
             total_change = delta, residual = delta - sum(contributions))
}
