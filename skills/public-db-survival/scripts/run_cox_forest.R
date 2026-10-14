#!/usr/bin/env Rscript
## Public database survival workflow: univariate + multivariate Cox model,
## proportional-hazards check (survival::cox.zph), events-per-variable (EPV)
## check, and a forest plot.
##
## Usage:
##   Rscript run_cox_forest.R expression_clinical.csv time_col status_col output_dir variable1 [variable2 ...]
##
## Only `survival` is strictly required. When `ggplot2` (and, if present,
## `survminer`) is installed, a polished ggplot2 forest plot is produced;
## otherwise the script falls back to a manual base-R forest plot
## (plot() + segments() + points()) so PNG/PDF figures are still produced
## with a minimal R install.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_cox_forest.R expression_clinical.csv time_col status_col output_dir variable1 [variable2 ...]")
}
table_file <- args[[1]]
time_col   <- args[[2]]
status_col <- args[[3]]
out_dir    <- args[[4]]
variables  <- args[5:length(args)]

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Missing required R package: survival. Install it with install.packages('survival').")
}
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
if (!has_ggplot2) {
  message("Note: ggplot2 not available. Falling back to a manual base-R forest plot ",
          "(plot() + segments() + points()). Install ggplot2 for the polished version.")
}

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(table_file, check.names = FALSE)
missing_cols <- setdiff(c(time_col, status_col, variables), colnames(dat))
if (length(missing_cols)) {
  stop(
    "Column(s) not found in '", table_file, "': ", paste(missing_cols, collapse = ", "), ". ",
    "Available columns: ", paste(colnames(dat), collapse = ", ")
  )
}

## Coerce time/covariates to numeric BEFORE complete.cases(), so a
## non-numeric value (e.g. "unknown") is caught and counted as dropped,
## instead of silently becoming NA after the filter and corrupting
## downstream Surv()/coxph() calculations without changing the reported
## "dropped for missing values" count.
n_raw <- nrow(dat)
if (!is.numeric(dat[[time_col]])) suppressWarnings(dat[[time_col]] <- as.numeric(as.character(dat[[time_col]])))
## Covariates: only coerce to numeric when every non-missing value genuinely parses as a
## number. A real clinical covariate (stage, sex, treatment arm, histology) is often
## categorical text -- forcing as.numeric() on it would turn the whole column into NA and
## complete.cases() would then silently drop every row, which previously made any
## categorical covariate unusable (this is a real bug fixed after live testing, not a
## hypothetical). Genuinely categorical columns are converted to a factor instead, which
## coxph()/model.matrix() already handle natively via dummy coding.
for (v in variables) {
  if (is.numeric(dat[[v]])) next
  raw <- dat[[v]]
  non_missing <- !is.na(raw) & nzchar(trimws(as.character(raw)))
  as_num <- suppressWarnings(as.numeric(as.character(raw)))
  if (sum(non_missing) > 0 && sum(!is.na(as_num) & non_missing) == sum(non_missing)) {
    dat[[v]] <- as_num  # every non-missing value parsed cleanly as numeric: genuinely continuous
  } else {
    dat[[v]] <- factor(raw)  # not cleanly numeric: treat as a categorical covariate
    message("Note: '", v, "' is not purely numeric; treating it as a categorical covariate (factor) with ",
            "levels: ", paste(levels(dat[[v]]), collapse = ", "))
  }
}
keep_cols <- c(time_col, status_col, variables)
dat <- dat[stats::complete.cases(dat[, keep_cols]), ]
n_dropped_na <- n_raw - nrow(dat)
if (n_dropped_na > 0) {
  message("Note: dropped ", n_dropped_na, " row(s) with missing/non-numeric time/status/covariate values.")
}
if (!nrow(dat)) stop("No complete rows remain after removing missing/non-numeric values in time/status/covariate columns.")

status_raw <- dat[[status_col]]
if (is.character(status_raw) || is.factor(status_raw)) {
  status_chr <- tolower(trimws(as.character(status_raw)))
  event <- ifelse(status_chr %in% c("1", "dead", "deceased", "event", "yes", "true"), 1,
            ifelse(status_chr %in% c("0", "alive", "censored", "no", "false"), 0, NA))
  if (any(is.na(event))) stop("Could not map all values of '", status_col, "' to 0/1.")
  dat$.event <- event
} else {
  dat$.event <- as.numeric(status_raw)
}

n_events <- sum(dat$.event == 1)
if (n_events == 0) {
  stop(
    "No events observed in this cohort: all ", nrow(dat), " samples have status=0 (censored/alive). ",
    "Cox hazard ratios and proportional-hazards checks are undefined for all-censored data ",
    "(survival::cox.zph itself errors on a zero-event Cox model). Check the '", status_col, "' column encoding, ",
    "or confirm this cohort truly has zero events before re-running."
  )
}
n_covariates <- length(variables)
epv <- n_events / n_covariates
epv_note <- if (epv < 10) {
  sprintf(
    "EPV (events-per-variable) = %.1f (%d events / %d covariates) is BELOW the ~10 EPV heuristic (Peduzzi et al. 1996). ",
    epv, n_events, n_covariates
  )
} else {
  sprintf("EPV (events-per-variable) = %.1f (%d events / %d covariates), meets the ~10 EPV heuristic.", epv, n_events, n_covariates)
}
if (epv < 10) {
  epv_note <- paste0(
    epv_note,
    "Consider reducing covariates to clinically pre-specified ones or using penalized/ridge Cox (glmnet) ",
    "instead of an unrestricted multivariate model (see references/pitfalls.md #8 and cns-method-guide.md)."
  )
}
message(epv_note)

fit_cox <- function(vars, dat) {
  rhs <- paste(sprintf("`%s`", vars), collapse = " + ")
  formula <- stats::as.formula(paste0("survival::Surv(`", time_col, "`, .event) ~ ", rhs))
  survival::coxph(formula, data = dat)
}

as_tab <- function(fit, model_label) {
  ci <- as.data.frame(summary(fit)$conf.int)
  coef <- as.data.frame(summary(fit)$coefficients)
  data.frame(
    model = model_label,
    variable = rownames(coef),
    HR = signif(ci$`exp(coef)`, 4),
    lower95 = signif(ci$`lower .95`, 4),
    upper95 = signif(ci$`upper .95`, 4),
    pvalue = signif(coef$`Pr(>|z|)`, 4),
    row.names = NULL
  )
}

## --- univariate Cox: one model per covariate --------------------------------
uni_fits <- lapply(variables, function(v) fit_cox(v, dat))
names(uni_fits) <- variables
uni_tab <- do.call(rbind, lapply(seq_along(variables), function(i) as_tab(uni_fits[[i]], "univariate")))

## --- multivariate Cox: all covariates together -------------------------------
multi_fit <- fit_cox(variables, dat)
multi_tab <- as_tab(multi_fit, "multivariate")

cat("\n=== Univariate Cox models (one covariate at a time) ===\n")
for (v in variables) {
  cat("\n---", v, "---\n")
  print(summary(uni_fits[[v]]))
}

cat("\n=== Multivariate Cox model (all covariates jointly) ===\n")
print(summary(multi_fit))

## --- proportional-hazards check for every Cox model fit ----------------------
## survival::cox.zph can itself error out (e.g. quasi-complete separation with
## very few events); catch that instead of letting the whole script crash and
## report "not assessed" for the PH check in that case.
safe_zph <- function(cox_fit, label) {
  tryCatch(survival::cox.zph(cox_fit), error = function(e) {
    message("Note: cox.zph(", label, ") failed (", conditionMessage(e), "); PH check not assessed for this model (likely too few events or quasi-complete separation).")
    structure(list(table = matrix(NA_real_, nrow = 1, ncol = 3, dimnames = list("GLOBAL", c("chisq", "df", "p")))), class = "cox.zph.unavailable")
  })
}

cat("\n=== Proportional-hazards check: cox.zph(multivariate model) ===\n")
multi_zph <- safe_zph(multi_fit, "multivariate model")
if (inherits(multi_zph, "cox.zph")) print(multi_zph) else cat("(not assessed; see note above)\n")

uni_zph <- lapply(variables, function(v) safe_zph(uni_fits[[v]], paste0("univariate: ", v)))
names(uni_zph) <- variables
for (v in variables) {
  cat("\n=== Proportional-hazards check: cox.zph(univariate model:", v, ") ===\n")
  if (inherits(uni_zph[[v]], "cox.zph")) print(uni_zph[[v]]) else cat("(not assessed; see note above)\n")
}

ph_status <- function(zph, label) {
  gp <- zph$table["GLOBAL", "p"]
  if (is.null(gp) || is.na(gp)) return(list(p = NA_real_, status = "not assessed", violated = FALSE, violated_terms = character(0)))
  terms <- setdiff(rownames(zph$table), "GLOBAL")
  term_p <- zph$table[terms, "p"]
  violated_terms <- terms[!is.na(term_p) & term_p < 0.05]
  if (gp < 0.05) {
    detail <- if (length(violated_terms)) {
      sprintf("VIOLATED for covariate(s) %s (global p=%.3g < 0.05)", paste(violated_terms, collapse = ", "), gp)
    } else {
      sprintf("VIOLATED (global p=%.3g < 0.05; no single covariate individually significant but the overall model trends over time)", gp)
    }
    status <- paste0(detail, ": HR(s) should be interpreted with caution as time-averaged approximations; consider a time-varying Cox term, interval-split Cox, or RMST (see references/pitfalls.md #2).")
  } else {
    status <- sprintf("satisfied for %s (global p=%.3g >= 0.05)", label, gp)
  }
  list(p = gp, status = status, violated = gp < 0.05, violated_terms = violated_terms)
}
multi_ph <- ph_status(multi_zph, "the multivariate model")
uni_ph <- lapply(variables, function(v) ph_status(uni_zph[[v]], v))
names(uni_ph) <- variables
cat("\nPH check summary: multivariate model ->", multi_ph$status, "\n")
for (v in variables) cat("PH check summary:", v, "(univariate) ->", uni_ph[[v]]$status, "\n")

ph_table <- rbind(
  data.frame(model = "multivariate", term = rownames(multi_zph$table),
             chisq = signif(multi_zph$table[, "chisq"], 4), df = multi_zph$table[, "df"],
             p = signif(multi_zph$table[, "p"], 4)),
  do.call(rbind, lapply(variables, function(v) {
    zt <- uni_zph[[v]]$table
    data.frame(model = paste0("univariate: ", v), term = rownames(zt),
               chisq = signif(zt[, "chisq"], 4), df = zt[, "df"], p = signif(zt[, "p"], 4))
  }))
)

## --- write result tables -----------------------------------------------------
write.table(uni_tab, file.path(out_dir, "results", "cox_univariate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(multi_tab, file.path(out_dir, "results", "cox_forest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ph_table, file.path(out_dir, "results", "ph_test_cox_forest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## --- forest plot (combined univariate + multivariate) ------------------------
combined <- rbind(uni_tab, multi_tab)
combined$variable <- factor(combined$variable, levels = rev(unique(combined$variable)))
combined$model <- factor(combined$model, levels = c("univariate", "multivariate"))

if (has_ggplot2) {
  plot_width <- 9.6
  ## Wrap the (often long) EPV note so it doesn't run off the right edge of
  ## the figure; extra wrapped lines get extra plot height.
  epv_wrapped <- paste(strwrap(epv_note, width = 95), collapse = "\n")
  n_subtitle_lines <- lengths(regmatches(epv_wrapped, gregexpr("\n", epv_wrapped))) + 1
  plot_height <- max(4, 1 + 0.7 * length(variables)) + 0.28 * n_subtitle_lines

  p <- ggplot2::ggplot(combined, ggplot2::aes(HR, variable, color = model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#A7C6DD") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower95, xmax = upper95), height = 0.18,
                             position = ggplot2::position_dodge(width = 0.5), linewidth = 0.9) +
    ggplot2::geom_point(size = 3.0, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_color_manual(values = c(univariate = "#667085", multivariate = "#D84A4A")) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Cox model forest plot",
      subtitle = epv_wrapped,
      x = "Hazard ratio (log scale)", y = NULL, color = NULL
    ) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#173B5F"),
      legend.position = "top"
    )
  ggplot2::ggsave(file.path(out_dir, "figures", "cox_forest.png"), p, width = plot_width, height = plot_height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", "cox_forest.pdf"), p, width = plot_width, height = plot_height, bg = "white")
} else {
  ## Manual base-R forest plot: two rows per variable (univariate + multivariate),
  ## log-scale HR axis, dashed reference line at HR=1, error bars via segments(),
  ## point estimates via points(), HR (95% CI) / p-value text annotated at right.
  plot_forest_base <- function(out_base, width = 9.4) {
    vars <- levels(combined$variable)
    n <- length(vars)
    epv_wrapped <- strwrap(epv_note, width = 118)
    mar_bottom <- 5.2 + 1.5 * length(epv_wrapped)
    height <- max(5.2, 1.3 + 0.85 * n + 0.5 * mar_bottom)
    finite_bounds <- combined[is.finite(combined$lower95) & is.finite(combined$upper95), ]
    xr <- range(c(finite_bounds$lower95, finite_bounds$upper95, 1))
    xr <- c(xr[1] / 1.3, xr[2] * 6)
    model_offset <- c(univariate = -0.16, multivariate = 0.16)
    model_color  <- c(univariate = "#667085", multivariate = "#D84A4A")

    draw <- function() {
      graphics::par(mar = c(mar_bottom, 11.5, 3.6, 11.5), xpd = FALSE)
      y_base <- setNames(n:1, vars)
      plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), log = "x", xaxt = "n", yaxt = "n",
           xlab = "", ylab = "",
           main = "Cox model forest plot", font.main = 2, cex.lab = 1.15)
      graphics::mtext("Hazard ratio (log scale)", side = 1, line = 2.3, cex = 1.0)
      graphics::axis(1, at = pretty(xr, n = 5))
      graphics::axis(2, at = y_base, labels = vars, las = 1, cex.axis = 1.0, tick = FALSE)
      graphics::abline(v = 1, lty = 2, col = "#A7C6DD")
      for (i in seq_len(nrow(combined))) {
        row <- combined[i, ]
        y <- y_base[[as.character(row$variable)]] + model_offset[[as.character(row$model)]]
        col <- model_color[[as.character(row$model)]]
        graphics::segments(row$lower95, y, row$upper95, y, col = col, lwd = 2)
        graphics::points(row$HR, y, pch = 18, cex = 1.5, col = col)
        graphics::text(x = xr[2], y = y, pos = 2, xpd = NA, cex = 0.72, col = col,
                        labels = sprintf("%.2f (%.2f-%.2f), p=%.3g", row$HR, row$lower95, row$upper95, row$pvalue))
      }
      graphics::legend("topright", inset = c(-0.32, 0), xpd = NA, bty = "n", cex = 0.95,
                        legend = c("Univariate", "Multivariate"), col = c("#667085", "#D84A4A"), pch = 18)
      for (i in seq_along(epv_wrapped)) {
        graphics::mtext(epv_wrapped[i], side = 1, line = 3.6 + i, cex = 0.62, col = "#173B5F", adj = 0)
      }
    }
    grDevices::png(paste0(out_base, ".png"), width = width, height = height, units = "in", res = 320, bg = "white")
    draw()
    grDevices::dev.off()
    grDevices::pdf(paste0(out_base, ".pdf"), width = width, height = height, bg = "white")
    draw()
    grDevices::dev.off()
  }
  plot_forest_base(file.path(out_dir, "figures", "cox_forest"))
}

## --- plain-language conclusion, keyed off the first covariate as the primary variable ---
## A factor/categorical covariate (e.g. `type` with levels Normal/Tumor) is expanded by
## coxph()/model.matrix() into per-level dummy coefficients named "<variable><level>" (e.g.
## "typeTumor"), NOT a row literally named "type". An exact match() against the raw variable
## name silently returns NA for every categorical covariate (previously misreported as a
## "degenerate fit" and HR=NA in the report even though the model fit fine) -- match by name
## prefix instead, and report every expanded level for multi-level factors.
rows_for_variable <- function(tab, v) {
  exact <- tab[tab$variable == v, , drop = FALSE]
  if (nrow(exact)) return(exact)
  tab[startsWith(tab$variable, v), , drop = FALSE]
}

primary_var <- variables[1]
primary_multi <- rows_for_variable(multi_tab, primary_var)
primary_uni <- rows_for_variable(uni_tab, primary_var)
if (nrow(primary_multi) > 1) {
  # multi-level factor: report the largest-|log(HR)| level as "primary" for the headline
  # conclusion sentence, since there is no single scalar HR for the variable as a whole
  primary_multi <- primary_multi[order(-abs(log(primary_multi$HR))), ][1, , drop = FALSE]
}
if (nrow(primary_uni) > 1) {
  primary_uni <- primary_uni[order(-abs(log(primary_uni$HR))), ][1, , drop = FALSE]
}

hr_lines <- unlist(lapply(variables, function(v) {
  uni_rows <- rows_for_variable(uni_tab, v)
  multi_rows <- rows_for_variable(multi_tab, v)
  levels_n <- max(nrow(uni_rows), nrow(multi_rows), 1)
  vapply(seq_len(levels_n), function(i) {
    u <- if (i <= nrow(uni_rows)) uni_rows[i, ] else data.frame(variable = v, HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
    m <- if (i <= nrow(multi_rows)) multi_rows[i, ] else data.frame(variable = v, HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
    label <- if (levels_n > 1) u$variable else v  # use the expanded "<var><level>" name when there's more than one dummy
    sprintf(
      "%s: univariate HR=%.3g (95%% CI %.3g-%.3g, p=%.3g); multivariate HR=%.3g (95%% CI %.3g-%.3g, p=%.3g)",
      label, u$HR, u$lower95, u$upper95, u$pvalue, m$HR, m$lower95, m$upper95, m$pvalue
    )
  }, character(1))
}))

caveats <- character(0)
if (epv < 10) {
  caveats <- c(caveats, sprintf("EPV=%.1f is below the ~10-events-per-variable heuristic, so multivariate coefficients/CIs may be unstable", epv))
}
if (isTRUE(multi_ph$violated)) {
  caveats <- c(caveats, "the multivariate proportional-hazards assumption was violated (see PH check above), so the multivariate HR(s) should be read as time-averaged approximations")
}
caveat_txt <- if (length(caveats)) paste0(" Caveats: ", paste(caveats, collapse = "; "), ".") else ""

oxford_join <- function(x) {
  if (length(x) <= 1) return(paste(x, collapse = ""))
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}

if (!nrow(primary_multi) || is.na(primary_multi$HR[1]) || is.na(primary_multi$pvalue[1])) {
  ## Degenerate fit (e.g. a constant/zero-variance covariate, or complete
  ## separation): coxph() returns NA coefficients rather than erroring, so
  ## guard here instead of crashing on `if (NA > 1)` below.
  conclusion <- paste0(
    "The Cox model for ", primary_var, " could not produce a usable hazard ratio (coefficient/CI came back NA), ",
    "which typically means ", primary_var, " (or another covariate) is constant, collinear, or perfectly separates ",
    "the outcome in this cohort. Check the input data (variance, duplicate columns, extreme class imbalance) before ",
    "drawing any conclusion; see the univariate/multivariate tables above for details."
  )
} else {
  # "per unit increase" only makes sense for a continuous covariate; a factor/categorical
  # covariate's coefficient is a level-vs-reference contrast (e.g. "typeTumor" = Tumor vs
  # the reference level Normal), so phrase it accordingly instead of implying a dose-response.
  is_primary_factor <- is.factor(dat[[primary_var]])
  unit_txt <- if (is_primary_factor) paste0("for ", primary_multi$variable[1], " vs. the reference level") else "per unit increase"
  direction <- if (primary_multi$HR[1] > 1) "higher hazard (worse survival)" else "lower hazard (better survival)"
  sig_txt <- if (primary_multi$pvalue[1] < 0.05) {
    sprintf("remained statistically significant after adjustment (p=%.3g)", primary_multi$pvalue[1])
  } else {
    sprintf("was not statistically significant after adjustment (p=%.3g)", primary_multi$pvalue[1])
  }
  indep_txt <- if (length(variables) > 1 && !is.na(primary_uni$pvalue[1]) && primary_uni$pvalue[1] < 0.05 && primary_multi$pvalue[1] >= 0.05) {
    " Note: this variable was significant univariately but lost significance after adjusting for the other covariates, i.e. its apparent effect may be confounded rather than independent (see references/pitfalls.md #9)."
  } else if (length(variables) > 1 && primary_multi$pvalue[1] < 0.05) {
    " This adjustment supports (but does not by itself definitively prove) independence from the other covariates in the model."
  } else {
    ""
  }

  conclusion <- paste0(
    "In the multivariate model, ", primary_var, " was associated with a ", direction,
    " ", unit_txt, " (HR=", signif(primary_multi$HR[1], 3), ", 95% CI ", signif(primary_multi$lower95[1], 3),
    "-", signif(primary_multi$upper95[1], 3), ") and ", sig_txt, " when adjusting for ",
    oxford_join(setdiff(variables, primary_var)), ".", indep_txt, caveat_txt,
    " This is a retrospective, single-cohort association and does not by itself establish causality (see references/pitfalls.md)."
  )
  if (length(variables) == 1) {
    sig_txt_single <- if (primary_multi$pvalue[1] < 0.05) {
      sprintf("statistically significant (p=%.3g)", primary_multi$pvalue[1])
    } else {
      sprintf("not statistically significant (p=%.3g)", primary_multi$pvalue[1])
    }
    conclusion <- paste0(
      "The univariate Cox model for ", primary_var, " gave HR=", signif(primary_multi$HR[1], 3), " (95% CI ",
      signif(primary_multi$lower95[1], 3), "-", signif(primary_multi$upper95[1], 3), "), which was ", sig_txt_single, ".",
      " No other covariates were included, so this result cannot be described as an independent prognostic factor ",
      "(see references/pitfalls.md #9).", caveat_txt
    )
  }
}

writeLines(c(
  "# Cox forest report",
  "",
  paste0("Covariates: ", paste(variables, collapse = ", ")),
  paste0("Samples used: ", nrow(dat), " (dropped for missing values: ", n_dropped_na, ")"),
  paste0("Events: ", n_events),
  epv_note,
  "",
  "## Hazard ratios (univariate vs. multivariate)",
  "",
  paste0("- ", hr_lines),
  "",
  paste0("Proportional-hazards check (multivariate model): PH assumption ", multi_ph$status),
  paste0("Proportional-hazards check (univariate models): ",
         paste(sprintf("%s: PH assumption %s", variables, vapply(uni_ph, `[[`, character(1), "status")), collapse = "; ")),
  "",
  "## Conclusion",
  "",
  conclusion,
  "",
  "## Limitations",
  "",
  "- Never call an expression covariate an \"independent prognostic factor\" from a univariate result alone; use the multivariate table above.",
  "- Read references/pitfalls.md before reporting any hazard ratio (proportional hazards, EPV, competing risks, and more)."
), file.path(out_dir, "report_cox_forest.md"))
