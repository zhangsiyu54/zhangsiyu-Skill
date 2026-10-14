#!/usr/bin/env Rscript
## Public database survival workflow: expression cutoff + Kaplan-Meier + log-rank
## + Cox proportional-hazards check (survival::cox.zph).
##
## Usage:
##   Rscript run_survival.R expression_clinical.csv gene time_col status_col output_dir [cutoff_method]
##   cutoff_method: median (default) | tertile | quartile | maxstat
##
## Only `survival` (base R recommended package) is strictly required. When
## `survminer`+`ggplot2` are installed, nicer ggsurvplot output (with a risk
## table panel) is produced; otherwise the script falls back to base-R
## graphics (plot.survfit + a manual risk-table panel) so PNG/PDF figures are
## still produced with a minimal R install.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript run_survival.R expression_clinical.csv gene time_col status_col output_dir [cutoff_method]\n",
    "  cutoff_method: median (default) | tertile | quartile | maxstat"
  )
}
table_file    <- args[[1]]
gene          <- args[[2]]
time_col      <- args[[3]]
status_col    <- args[[4]]
out_dir       <- args[[5]]
cutoff_method <- if (length(args) >= 6) args[[6]] else "median"

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Missing required R package: survival. Install it with install.packages('survival').")
}

has_survminer <- requireNamespace("survminer", quietly = TRUE)
has_ggplot2   <- requireNamespace("ggplot2", quietly = TRUE)
use_survminer <- has_survminer && has_ggplot2
if (!use_survminer) {
  message(
    "Note: survminer/ggplot2 not both available (survminer=", has_survminer,
    ", ggplot2=", has_ggplot2, "). Falling back to base-R graphics for the KM plot ",
    "(plot.survfit + manual risk table). Install survminer+ggplot2 for the polished ggsurvplot output."
  )
}

dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

## --- load + basic QC --------------------------------------------------------
dat <- read.csv(table_file, check.names = FALSE)
missing_cols <- setdiff(c(gene, time_col, status_col), colnames(dat))
if (length(missing_cols)) {
  stop(
    "Column(s) not found in '", table_file, "': ", paste(missing_cols, collapse = ", "), ". ",
    "Available columns: ", paste(colnames(dat), collapse = ", ")
  )
}

## Coerce gene/time to numeric BEFORE the missing-value filter, so a
## non-numeric value (e.g. "unknown" in the time column) is caught and
## counted as dropped, instead of silently becoming NA after the filter and
## corrupting downstream Surv()/coxph() calculations without changing the
## reported "dropped for missing values" count.
n_raw <- nrow(dat)
if (!is.numeric(dat[[gene]])) {
  suppressWarnings(dat[[gene]] <- as.numeric(as.character(dat[[gene]])))
}
if (!is.numeric(dat[[time_col]])) {
  suppressWarnings(dat[[time_col]] <- as.numeric(as.character(dat[[time_col]])))
}
dat <- dat[!is.na(dat[[gene]]) & !is.na(dat[[time_col]]) & !is.na(dat[[status_col]]) & dat[[status_col]] != "", ]
n_dropped_na <- n_raw - nrow(dat)
if (n_dropped_na > 0) {
  message(
    "Note: dropped ", n_dropped_na, " row(s) with missing/non-numeric '", gene, "', '", time_col,
    "', or empty '", status_col, "' values."
  )
}
if (!nrow(dat)) stop("No complete rows remain after removing missing/non-numeric gene/time/status values.")

status_raw <- dat[[status_col]]
if (is.character(status_raw) || is.factor(status_raw)) {
  status_chr <- tolower(trimws(as.character(status_raw)))
  event <- ifelse(status_chr %in% c("1", "dead", "deceased", "event", "yes", "true"), 1,
            ifelse(status_chr %in% c("0", "alive", "censored", "no", "false"), 0, NA))
  if (any(is.na(event))) {
    stop("Could not map all values of '", status_col, "' to 0/1. Unique values found: ",
         paste(unique(status_chr), collapse = ", "))
  }
  dat$.event <- event
} else {
  dat$.event <- as.numeric(status_raw)
}
if (!all(dat$.event %in% c(0, 1))) stop("Event status must code to 0 (censored/alive) / 1 (event/dead).")
if (any(dat[[time_col]] <= 0, na.rm = TRUE)) {
  warning("Non-positive survival time values found; check units/censoring definition.")
}

## --- cutoff rule (see references/pitfalls.md #1 and cns-method-guide.md) ---
apply_cutoff <- function(dat, gene, method) {
  x <- dat[[gene]]
  method <- tolower(method)
  n_before <- nrow(dat)
  if (method == "median") {
    cutoff <- stats::median(x, na.rm = TRUE)
    dat$expression_group <- ifelse(x >= cutoff, "High", "Low")
    note <- sprintf(
      "Median split at %.4g (pre-specified rule; primary/default, low risk of cutoff-selection bias).",
      cutoff
    )
  } else if (method %in% c("tertile", "quartile")) {
    probs <- if (method == "tertile") c(1 / 3, 2 / 3) else c(0.25, 0.75)
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    grp <- ifelse(x <= q[1], "Low", ifelse(x >= q[2], "High", "Middle"))
    dat$expression_group <- grp
    dat <- dat[dat$expression_group != "Middle", ]
    cutoff <- q
    label <- if (method == "tertile") c("P33", "P67") else c("P25", "P75")
    note <- sprintf(
      "%s split (<=%s Low, >=%s High); middle group (n=%d) excluded from KM/log-rank/Cox.",
      tools::toTitleCase(method), label[1], label[2], n_before - nrow(dat)
    )
  } else if (method == "maxstat") {
    if (!has_survminer) {
      stop(
        "cutoff_method='maxstat' requires the survminer package (surv_cutpoint). ",
        "Install survminer, or use median/tertile/quartile instead. ",
        "See references/pitfalls.md #1: a raw maxstat log-rank p is optimistically biased ",
        "and must not be reported as a primary result without correction/validation."
      )
    }
    sc <- survminer::surv_cutpoint(dat, time = time_col, event = ".event", variables = gene)
    cutoff <- sc$cutpoint$cutpoint
    catdat <- survminer::surv_categorize(sc)
    dat$expression_group <- ifelse(catdat[[gene]] == "high", "High", "Low")
    note <- paste0(
      "EXPLORATORY maxstat cutoff at ", signif(cutoff, 4), " (survminer::surv_cutpoint, ",
      "maximally selected rank statistic). This maximizes group separation, so the log-rank p ",
      "below is OPTIMISTICALLY BIASED (Hothorn & Lausen 2003) and must be labeled exploratory ",
      "and confirmed in an independent validation cohort before any prognostic claim ",
      "(see references/pitfalls.md #1)."
    )
  } else {
    stop("Unknown cutoff_method: '", method, "'. Use one of: median, tertile, quartile, maxstat.")
  }
  list(dat = dat, cutoff = cutoff, note = note, method = method)
}

cr <- apply_cutoff(dat, gene, cutoff_method)
dat <- cr$dat
dat$expression_group <- factor(dat$expression_group, levels = c("Low", "High"))
dat <- droplevels(dat)
if (nlevels(dat$expression_group) < 2 || any(table(dat$expression_group) < 2)) {
  stop("Cutoff rule '", cutoff_method, "' did not produce two groups with >=2 samples each; ",
       "check the gene column, cutoff method, or cohort size.")
}

group_n      <- table(dat$expression_group)
group_events <- tapply(dat$.event, dat$expression_group, sum)

n_total_events <- sum(dat$.event)
if (n_total_events == 0) {
  stop(
    "No events observed in this cohort: all ", nrow(dat), " samples have status=0 (censored/alive). ",
    "Kaplan-Meier median survival, log-rank tests, and Cox hazard ratios are all undefined for all-censored data ",
    "(survival::cox.zph itself errors on a zero-event Cox model). Check the '", status_col, "' column encoding, ",
    "or confirm this cohort truly has zero events before re-running."
  )
}

## --- KM + log-rank + Cox check ----------------------------------------------
surv_obj <- survival::Surv(dat[[time_col]], dat$.event)
fit      <- survival::survfit(surv_obj ~ expression_group, data = dat)
logrank  <- survival::survdiff(surv_obj ~ expression_group, data = dat)
logrank_p <- stats::pchisq(logrank$chisq, length(logrank$n) - 1, lower.tail = FALSE)

cox_group      <- survival::coxph(surv_obj ~ expression_group, data = dat)
cox_cont       <- survival::coxph(survival::Surv(dat[[time_col]], dat$.event) ~ dat[[gene]])

## survival::cox.zph can itself error out (e.g. quasi-complete separation with
## very few events); catch that instead of letting the whole script crash and
## report "not assessed" for the PH check in that case.
safe_zph <- function(cox_fit, label) {
  tryCatch(survival::cox.zph(cox_fit), error = function(e) {
    message("Note: cox.zph(", label, ") failed (", conditionMessage(e), "); PH check not assessed for this model (likely too few events or quasi-complete separation).")
    structure(list(table = matrix(NA_real_, nrow = 1, ncol = 3, dimnames = list("GLOBAL", c("chisq", "df", "p")))), class = "cox.zph.unavailable")
  })
}
cox_group_zph  <- safe_zph(cox_group, "expression_group model")
cox_cont_zph   <- safe_zph(cox_cont, paste0(gene, " continuous model"))

cat("\n=== Log-rank test (expression_group: High vs Low) ===\n")
print(logrank)
cat(sprintf("Log-rank chi-sq = %.4f, df = %d, p = %.4g\n", logrank$chisq, length(logrank$n) - 1, logrank_p))

cat("\n=== Cox model: Surv(time, status) ~ expression_group (categorical HR) ===\n")
print(summary(cox_group))

cat("\n=== Proportional-hazards check: cox.zph(cox_group) ===\n")
if (inherits(cox_group_zph, "cox.zph")) print(cox_group_zph) else cat("(not assessed; see note above)\n")

cat("\n=== Cox model: Surv(time, status) ~ ", gene, " (continuous, per-unit HR) ===\n", sep = "")
print(summary(cox_cont))

cat("\n=== Proportional-hazards check: cox.zph(cox_continuous) ===\n")
if (inherits(cox_cont_zph, "cox.zph")) print(cox_cont_zph) else cat("(not assessed; see note above)\n")

ph_status <- function(zph, term_label) {
  gp <- zph$table["GLOBAL", "p"]
  if (is.na(gp)) return(list(p = NA_real_, status = "not assessed", violated = FALSE))
  if (gp < 0.05) {
    list(
      p = gp, violated = TRUE,
      status = sprintf(
        "VIOLATED for %s (global p=%.3g < 0.05): the reported HR should be interpreted with caution as a time-averaged approximation; consider a time-varying Cox term, interval-split Cox, or RMST (see references/pitfalls.md #2).",
        term_label, gp
      )
    )
  } else {
    list(p = gp, violated = FALSE, status = sprintf("satisfied for %s (global p=%.3g >= 0.05)", term_label, gp))
  }
}
ph_group <- ph_status(cox_group_zph, "expression_group")
ph_cont  <- ph_status(cox_cont_zph, gene)
cat("\nPH check summary: group model ->", ph_group$status, "\n")
cat("PH check summary: continuous model ->", ph_cont$status, "\n\n")

## --- figures ----------------------------------------------------------------
group_colors <- c(Low = "#4EA3F1", High = "#E45B5B")
hr_ci <- summary(cox_group)$conf.int
hr_label <- sprintf(
  "HR(High vs Low) = %.2f (95%% CI %.2f-%.2f), log-rank p = %.3g",
  hr_ci[1, "exp(coef)"], hr_ci[1, "lower .95"], hr_ci[1, "upper .95"], logrank_p
)

if (use_survminer) {
  save_gg <- function(plot, name, width = 8.4, height = 7.2) {
    ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
    ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
  }
  gg <- survminer::ggsurvplot(
    fit, data = dat, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    censor = TRUE, palette = c(group_colors[["Low"]], group_colors[["High"]]),
    legend.labs = paste0(names(group_n), " (n=", group_n, ")"),
    title = "Kaplan-Meier curve",
    xlab = "Time", ylab = "Survival probability",
    risk.table.title = "Number at risk",
    ggtheme = ggplot2::theme_classic(base_size = 15)
  )
  gg$plot <- gg$plot +
    ggplot2::labs(subtitle = paste0(cr$note, "\n", hr_label)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#173B5F"),
      axis.title = ggplot2::element_text(size = 17, color = "#173B5F"),
      axis.text = ggplot2::element_text(size = 13, color = "#173B5F"),
      legend.position = "top"
    )
  combined <- survminer::arrange_ggsurvplots(list(gg), print = FALSE, ncol = 1, nrow = 1)
  ggplot2::ggsave(file.path(out_dir, "figures", "km_curve.png"), combined, width = 8.4, height = 7.4, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", "km_curve.pdf"), combined, width = 8.4, height = 7.4, bg = "white")
} else {
  ## Base-R fallback: KM curve (with censoring ticks) + manual risk-table panel.
  plot_km_base <- function(out_base, width = 8.2, height = 7.4) {
    time_breaks <- pretty(c(0, max(fit$time)), n = 6)
    time_breaks <- sort(unique(pmax(0, pmin(time_breaks, max(fit$time)))))
    strata_names <- sub(".*=", "", names(fit$strata))
    smy <- summary(fit, times = time_breaks, extend = TRUE)

    draw <- function() {
      graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1.2))
      graphics::par(mar = c(2.3, 5.2, 4.2, 2))
      plot(fit, col = group_colors[strata_names], lwd = 2.6, mark.time = TRUE, pch = 3,
           xaxt = "n", xlab = "", ylab = "Survival probability", ylim = c(0, 1.05),
           main = "Kaplan-Meier curve", font.main = 2, cex.lab = 1.25, cex.axis = 1.05)
      graphics::axis(1, at = time_breaks)
      graphics::legend("bottomleft", legend = paste0(strata_names, " (n=", group_n[strata_names], ")"),
                        col = group_colors[strata_names], lwd = 2.6, bty = "n", cex = 1.05)
      graphics::mtext(cr$note, side = 3, line = 0.9, cex = 0.72, col = "#173B5F")
      graphics::mtext(hr_label, side = 3, line = 0.1, cex = 0.78, col = "#173B5F")

      graphics::par(mar = c(3.2, 5.2, 1.3, 2))
      plot(NA, xlim = range(time_breaks), ylim = c(0.5, length(strata_names) + 0.5),
           xaxt = "n", yaxt = "n", xlab = "Time", ylab = "", bty = "n")
      graphics::axis(1, at = time_breaks)
      graphics::axis(2, at = length(strata_names):1, labels = strata_names, las = 1, tick = FALSE, cex.axis = 1.05)
      for (i in seq_along(strata_names)) {
        y <- length(strata_names) - i + 1
        idx <- smy$strata == names(fit$strata)[i]
        graphics::text(x = time_breaks, y = y, labels = smy$n.risk[idx], cex = 0.95, col = group_colors[strata_names[i]])
      }
      graphics::mtext("Number at risk", side = 3, line = -0.1, adj = 0, cex = 0.85, col = "#173B5F")
    }

    grDevices::png(paste0(out_base, ".png"), width = width, height = height, units = "in", res = 320, bg = "white")
    draw()
    grDevices::dev.off()
    grDevices::pdf(paste0(out_base, ".pdf"), width = width, height = height, bg = "white")
    draw()
    grDevices::dev.off()
  }
  plot_km_base(file.path(out_dir, "figures", "km_curve"))
}

## --- result tables ------------------------------------------------------------
km_table <- summary(fit)$table
median_rows <- lapply(rownames(km_table), function(nm) {
  grp <- sub(".*=", "", nm)
  med <- km_table[nm, "median"]
  last_time <- max(dat[[time_col]][dat$expression_group == grp])
  data.frame(
    group = grp,
    n = group_n[[grp]],
    events = group_events[[grp]],
    median_survival = if (is.na(med)) "not reached" else signif(med, 4),
    last_followup = signif(last_time, 4)
  )
})
group_summary <- do.call(rbind, median_rows)

survival_summary <- data.frame(
  gene = gene,
  cutoff_method = cr$method,
  cutoff_note = cr$note,
  samples_used = nrow(dat),
  samples_dropped_missing = n_dropped_na,
  logrank_chisq = signif(logrank$chisq, 4),
  logrank_df = length(logrank$n) - 1,
  logrank_p = signif(logrank_p, 4),
  HR_high_vs_low = signif(hr_ci[1, "exp(coef)"], 4),
  HR_lower95 = signif(hr_ci[1, "lower .95"], 4),
  HR_upper95 = signif(hr_ci[1, "upper .95"], 4),
  HR_pvalue = signif(summary(cox_group)$coefficients[1, "Pr(>|z|)"], 4),
  cox_zph_group_global_p = signif(ph_group$p, 4),
  ph_status_group = ph_group$status,
  cox_zph_continuous_global_p = signif(ph_cont$p, 4),
  ph_status_continuous = ph_cont$status
)

write.table(group_summary, file.path(out_dir, "results", "survival_group_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(survival_summary, file.path(out_dir, "results", "survival_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dat[, setdiff(colnames(dat), ".event")], file.path(out_dir, "results", "expression_group.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cox_cont_ci <- summary(cox_cont)$conf.int
cox_summary_tab <- data.frame(
  term = gene,
  coef = signif(summary(cox_cont)$coefficients[1, "coef"], 4),
  HR = signif(cox_cont_ci[1, "exp(coef)"], 4),
  lower95 = signif(cox_cont_ci[1, "lower .95"], 4),
  upper95 = signif(cox_cont_ci[1, "upper .95"], 4),
  pvalue = signif(summary(cox_cont)$coefficients[1, "Pr(>|z|)"], 4)
)
write.table(cox_summary_tab, file.path(out_dir, "results", "cox_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

ph_table <- rbind(
  data.frame(model = "expression_group (categorical)", term = rownames(cox_group_zph$table),
             chisq = signif(cox_group_zph$table[, "chisq"], 4), df = cox_group_zph$table[, "df"],
             p = signif(cox_group_zph$table[, "p"], 4)),
  data.frame(model = paste0(gene, " (continuous)"), term = rownames(cox_cont_zph$table),
             chisq = signif(cox_cont_zph$table[, "chisq"], 4), df = cox_cont_zph$table[, "df"],
             p = signif(cox_cont_zph$table[, "p"], 4))
)
write.table(ph_table, file.path(out_dir, "results", "ph_test.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## --- report + top-level summary ----------------------------------------------
median_line <- paste0(
  "Median survival: ",
  paste(sprintf("%s = %s (last follow-up %s)", group_summary$group, group_summary$median_survival, group_summary$last_followup), collapse = "; ")
)

hr_val <- hr_ci[1, "exp(coef)"]
direction <- if (hr_val > 1) "worse" else "better"
sig_txt <- if (logrank_p < 0.05) {
  sprintf("a statistically significant difference (log-rank p=%.3g)", logrank_p)
} else {
  sprintf("no statistically significant difference (log-rank p=%.3g)", logrank_p)
}
caveats <- character(0)
if (isTRUE(ph_group$violated) || isTRUE(ph_cont$violated)) {
  caveats <- c(caveats, "the proportional-hazards assumption was violated for at least one model above, so this HR should be read as a time-averaged approximation rather than a constant effect over follow-up")
}
if (identical(cr$method, "maxstat")) {
  caveats <- c(caveats, "the expression cutoff was chosen by a data-driven (maxstat) search, which optimistically biases the log-rank p value at the winning cutoff; treat this result as exploratory until confirmed with a corrected p value in an independent validation cohort")
}
if (any(grepl("not reached", group_summary$median_survival))) {
  caveats <- c(caveats, "median survival was not reached in at least one group within the observed follow-up, so any such median is an extrapolation beyond the data; the fixed-time group sizes/events above are more reliable than an unreached median")
}
if (!is.finite(hr_ci[1, "upper .95"]) || hr_ci[1, "lower .95"] == 0 || hr_val > 1000 || hr_val < 1e-3) {
  caveats <- c(caveats, sprintf(
    "the 95%% CI is extremely wide or unbounded (Low events=%d, High events=%d), consistent with quasi-complete separation from very few events; treat this HR as unstable/exploratory rather than a precise estimate (see references/pitfalls.md #8)",
    group_events[["Low"]], group_events[["High"]]
  ))
} else if (min(group_events) < 10) {
  caveats <- c(caveats, sprintf(
    "at least one group has fewer than 10 events (Low=%d, High=%d), so this HR/CI is relatively imprecise (see references/pitfalls.md #8)",
    group_events[["Low"]], group_events[["High"]]
  ))
}
caveat_txt <- if (length(caveats)) paste0(" Caveats: ", paste(caveats, collapse = "; "), ".") else ""

conclusion <- paste0(
  "In this cohort, patients with High ", gene, " expression showed ", direction,
  " survival than patients with Low expression, with ", sig_txt,
  " and a Cox hazard ratio of ", signif(hr_val, 3), " (95% CI ", signif(hr_ci[1, "lower .95"], 3),
  "-", signif(hr_ci[1, "upper .95"], 3), ") for High vs Low.", caveat_txt,
  " This is a single-cohort, retrospective association and does not by itself establish causality, clinical utility, or independence from clinical confounders (see references/pitfalls.md)."
)

writeLines(c(
  "# Public database survival report",
  "",
  paste0("Target gene: ", gene),
  paste0("Cutoff rule: ", cr$note),
  paste0("Samples used: ", nrow(dat), " (dropped for missing values: ", n_dropped_na, ")"),
  paste0("Group sizes: ", paste(sprintf("%s n=%d (events=%d)", names(group_n), group_n, group_events[names(group_n)]), collapse = "; ")),
  median_line,
  paste0("Log-rank test: chi-sq=", signif(logrank$chisq, 4), ", df=", length(logrank$n) - 1, ", p=", signif(logrank_p, 4)),
  paste0("Cox HR (High vs Low): ", hr_label),
  paste0("Proportional-hazards check (expression_group model): PH assumption ", ph_group$status),
  paste0("Proportional-hazards check (continuous ", gene, " model): PH assumption ", ph_cont$status),
  "",
  "## Conclusion",
  "",
  conclusion,
  "",
  "## Limitations",
  "",
  "- Mark small or simulated cohorts as exploratory; this KM/log-rank/Cox result does not by itself establish clinical utility.",
  "- Read references/pitfalls.md before reporting any hazard ratio or prognostic claim (cutoff-selection bias, unchecked proportional hazards, immortal time bias, competing risks, and more)."
), file.path(out_dir, "report.md"))

write.table(
  data.frame(
    item = c("samples", "events_high", "events_low", "cutoff_method", "cutoff_note", "logrank_p", "HR_high_vs_low", "HR_p", "ph_violated_group", "ph_violated_continuous"),
    value = c(
      nrow(dat),
      group_events[["High"]], group_events[["Low"]],
      cr$method, cr$note, signif(logrank_p, 4),
      signif(hr_ci[1, "exp(coef)"], 4), signif(summary(cox_group)$coefficients[1, "Pr(>|z|)"], 4),
      grepl("VIOLATED", ph_group$status), grepl("VIOLATED", ph_cont$status)
    )
  ),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Wrote survival outputs to", out_dir, "\n")
