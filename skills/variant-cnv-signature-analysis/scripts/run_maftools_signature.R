#!/usr/bin/env Rscript
# Formal R/maftools route for variant-cnv-signature-analysis.
#
# NOTE ON VERIFICATION: this script was written and manually reviewed against
# the maftools documentation/vignette (read.maf, oncoplot, plotmafSummary,
# tmb, trinucleotideMatrix, estimateSignatures, extractSignatures,
# compareSignatures, plotSignatures) but could NOT be executed in this
# sandbox: R 4.6.1 is present but maftools, NMF and BSgenome reference data
# packages are not installed here. Before production use, run it once in an
# environment with those packages installed and confirm the exact argument
# and return-value field names against the maftools version you have
# (function names/arguments are stable across recent releases, but this was
# not confirmed by execution). Several steps defensively save the raw R
# object (.rds) alongside a best-effort CSV/TSV extraction so nothing is
# silently lost if a list field name differs from what is assumed here.
#
# Usage:
#   Rscript run_maftools_signature.R maf_file [output_dir] [ref_genome] [capture_size_mb] [clinical_file] [n_signatures]
#
# Arguments:
#   maf_file        Path to a MAF-like tab-delimited file (see
#                    references/input-contract.md; scripts/make_demo_data.py
#                    writes a compatible demo.maf).
#   output_dir      Default "output".
#   ref_genome      "hg19" or "hg38". Default "hg19". Selects the BSgenome
#                    package used for trinucleotide-context extraction. If the
#                    matching BSgenome package is not installed, the signature
#                    section is skipped (oncoplot/TMB/summary still run).
#   capture_size_mb Numeric callable footprint in Mb used as the TMB
#                    denominator (see references/pitfalls.md #1: never use a
#                    generic "30 Mb exome" constant for panel data). Default 50.
#   clinical_file   Optional tab-delimited file with Tumor_Sample_Barcode plus
#                    a grouping column, used to annotate/sort the oncoplot.
#   n_signatures    Optional fixed signature count for extractSignatures().
#                    estimateSignatures()/plotCophenetic() do not return a
#                    single programmatic "best rank" in maftools -- rank
#                    selection is normally done by visually inspecting the
#                    cophenetic-correlation elbow plot. Default 3 (matches the
#                    maftools vignette's own low-mutation-count example);
#                    override after reviewing figures/signature_cophenetic.pdf
#                    from a first run.
#
# This script does NOT perform purity/ploidy-aware CNV calling. CNV in this
# skill is handled by the gene/segment matrix + Python quickplot route
# (scripts/quickplot_cancer_genomics.py). For purity/ploidy-aware absolute
# copy number, run a dedicated caller upstream (ASCAT, FACETS, PureCN,
# Sequenza) and see references/pitfalls.md #4 and #12 before interpreting any
# CNV output alongside this script's mutation/signature results.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_maftools_signature.R maf_file [output_dir] [ref_genome] [capture_size_mb] [clinical_file] [n_signatures]")
}
maf_file <- args[[1]]
out_dir <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "output"
ref_genome_key <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "hg19"
capture_size_mb <- if (length(args) >= 4 && nzchar(args[[4]])) as.numeric(args[[4]]) else 50
clinical_file <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else NA
n_signatures_override <- if (length(args) >= 6 && nzchar(args[[6]])) as.numeric(args[[6]]) else NA

if (!file.exists(maf_file)) stop("MAF file not found: ", maf_file)
if (is.na(capture_size_mb) || capture_size_mb <= 0) {
  stop("capture_size_mb must be a positive number (assay-specific callable Mb; see references/pitfalls.md #1)")
}

required <- c("maftools")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall via BiocManager::install('maftools') before running the formal R route.")
}

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

log_lines <- character(0)
log <- function(...) {
  msg <- paste0(...)
  message(msg)
  log_lines <<- c(log_lines, msg)
}

# Minimal JSON writer so the script does not require the 'jsonlite' package
# just to log run parameters (keeps the dependency list to maftools + NMF +
# a BSgenome data package). Defined before first use (top-level R scripts
# execute sequentially; a function must exist before it is called).
manual_json <- function(x) {
  esc <- function(v) gsub('"', '\\\\"', as.character(v))
  kv <- vapply(names(x), function(k) paste0('"', k, '": "', esc(x[[k]]), '"'), character(1))
  paste0("{\n  ", paste(kv, collapse = ",\n  "), "\n}")
}

# --- helper: save a base-graphics plotting call as PNG + PDF ----------------
# maftools plotting functions (oncoplot, plotmafSummary, tmb, plotSignatures,
# plotCophenetic) draw directly to the open graphics device as a side effect
# and do not return a re-playable plot object (unlike ggplot2, whose objects
# run_deseq2.R saves via ggplot2::ggsave()). `expr` is captured UNEVALUATED
# via substitute() and re-evaluated once per device with eval(); simply
# forcing a promise twice would only draw once; R evaluates a lazy argument
# promise the first time it is accessed and caches that result, so passing
# an already-evaluated value here would silently produce a blank second file.
save_base_plot <- function(name, expr, width = 10, height = 8, res = 320) {
  expr <- substitute(expr)
  env <- parent.frame()
  png_path <- file.path(out_dir, "figures", paste0(name, ".png"))
  pdf_path <- file.path(out_dir, "figures", paste0(name, ".pdf"))
  grDevices::png(png_path, width = width, height = height, units = "in", res = res, bg = "white")
  ok <- tryCatch({ eval(expr, envir = env); TRUE },
                 error = function(e) { log("Plot failed (", name, "): ", conditionMessage(e)); FALSE })
  grDevices::dev.off()
  if (ok) {
    grDevices::pdf(pdf_path, width = width, height = height, bg = "white")
    tryCatch(eval(expr, envir = env),
             error = function(e) log("PDF re-render failed (", name, "): ", conditionMessage(e)))
    grDevices::dev.off()
  }
  invisible(ok)
}

# --- helper: run a compute-only call without leaking a stray default plot --
# Some maftools/NMF functions may draw to the current device even when only
# their return value is wanted (e.g. tmb() plots as well as returning a
# data.table). Calling any plotting-capable function with no device open in
# Rscript triggers R's default device, which would otherwise leave a stray
# "Rplots.pdf" in the working directory. pdf(NULL) opens a device that
# accepts drawing commands but writes no file.
with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expr
}

# --- 1. Read MAF (+ optional clinical annotation) ---------------------------
clinical_data <- NULL
if (!is.na(clinical_file)) {
  if (!file.exists(clinical_file)) stop("clinical_file not found: ", clinical_file)
  clinical_data <- clinical_file
}

maf <- maftools::read.maf(maf = maf_file, clinicalData = clinical_data, verbose = TRUE)
log("Read MAF: ", maf_file)

sample_summary <- maftools::getSampleSummary(maf)
gene_summary <- maftools::getGeneSummary(maf)
n_samples <- nrow(sample_summary)
log("Samples: ", n_samples, " | Genes: ", nrow(gene_summary))

utils::write.table(as.data.frame(sample_summary), file.path(out_dir, "results", "sample_summary.tsv"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(as.data.frame(gene_summary), file.path(out_dir, "results", "gene_summary.tsv"),
                    sep = "\t", row.names = FALSE, quote = FALSE)

# --- 2. Oncoplot (top mutated genes) ----------------------------------------
top_n <- min(20, nrow(gene_summary))
onco_args <- list(maf = maf, top = top_n)
if (!is.null(clinical_data)) {
  clin_df <- tryCatch(as.data.frame(maftools::getClinicalData(maf)), error = function(e) NULL)
  group_cols <- if (!is.null(clin_df)) setdiff(colnames(clin_df), "Tumor_Sample_Barcode") else character(0)
  if (length(group_cols) >= 1) {
    onco_args$clinicalFeatures <- group_cols[[1]]
    onco_args$sortByAnnotation <- TRUE
  }
}
save_base_plot("oncoplot", do.call(maftools::oncoplot, onco_args), width = 11, height = 8.5)

# --- 3. MAF summary dashboard (variant classification, TiTv, burden panel) --
save_base_plot("maf_summary_dashboard",
               maftools::plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE),
               width = 11, height = 8.5)

# --- 4. TMB (mutations / Mb), explicit assay-specific denominator ----------
# captureSize is the assay's own callable footprint in Mb (pitfalls.md #1):
# do not use a generic "~30-38 Mb exome" constant for targeted-panel data.
tmb_res <- tryCatch(
  with_null_device(maftools::tmb(maf = maf, captureSize = capture_size_mb, logScale = TRUE)),
  error = function(e) { log("tmb() failed: ", conditionMessage(e)); NULL }
)
if (!is.null(tmb_res)) {
  save_base_plot("tmb_summary_boxplot",
                 maftools::tmb(maf = maf, captureSize = capture_size_mb, logScale = TRUE),
                 width = 8, height = 7)
  utils::write.table(as.data.frame(tmb_res), file.path(out_dir, "results", "tmb_summary.tsv"),
                      sep = "\t", row.names = FALSE, quote = FALSE)
  log("TMB computed with captureSize = ", capture_size_mb, " Mb (state this denominator and the assay it corresponds to in methods)")
} else {
  log("TMB step skipped due to error above.")
}

# --- 5. Mutational signatures: trinucleotide context -> de novo NMF --------
# Route selection follows references/cns-method-guide.md and pitfalls.md #2/#3:
#   - de novo extraction only when the cohort is large enough for a stable
#     NMF factorization (commonly >=20 samples with adequate mutations);
#   - otherwise skip de novo and report only the mutation-context summary,
#     flagging the cohort as too small for reliable single-run extraction and
#     pointing to reference-based refitting tools (deconstructSigs,
#     SigProfilerAssignment) as the appropriate alternative for small cohorts.
MIN_SAMPLES_DE_NOVO <- 20
bsgenome_pkg <- switch(ref_genome_key,
  hg19 = "BSgenome.Hsapiens.UCSC.hg19",
  hg38 = "BSgenome.Hsapiens.UCSC.hg38",
  stop("Unsupported ref_genome '", ref_genome_key, "'; use 'hg19' or 'hg38'")
)

have_bsgenome <- requireNamespace(bsgenome_pkg, quietly = TRUE)
have_nmf <- requireNamespace("NMF", quietly = TRUE)

if (!have_bsgenome) {
  log("Skipping signature extraction: ", bsgenome_pkg, " is not installed. ",
      "Install it (BiocManager::install('", bsgenome_pkg, "')) to enable trinucleotideMatrix().")
} else if (!have_nmf) {
  log("Skipping signature extraction: package 'NMF' is not installed (required by estimateSignatures/extractSignatures).")
} else {
  tnm <- tryCatch(
    maftools::trinucleotideMatrix(maf = maf, ref_genome = bsgenome_pkg, prefix = NULL, add = TRUE, useSyn = TRUE),
    error = function(e) { log("trinucleotideMatrix() failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(tnm)) {
    tryCatch(utils::write.csv(tnm$nmf_matrix, file.path(out_dir, "results", "trinucleotide_matrix.csv")),
             error = function(e) log("Could not write trinucleotide_matrix.csv: ", conditionMessage(e)))
    if (!is.null(tnm$APOBEC_scores)) {
      utils::write.table(as.data.frame(tnm$APOBEC_scores), file.path(out_dir, "results", "apobec_scores.tsv"),
                          sep = "\t", row.names = FALSE, quote = FALSE)
    }

    if (n_samples < MIN_SAMPLES_DE_NOVO) {
      log("Cohort has ", n_samples, " samples (< ", MIN_SAMPLES_DE_NOVO, "). ",
          "Skipping de novo NMF signature extraction as unstable for this cohort size (pitfalls.md #2/#3). ",
          "Wrote trinucleotide context matrix only; use a reference-based refitting tool ",
          "(e.g. deconstructSigs, SigProfilerAssignment) against COSMIC SBS signatures instead of de novo extraction.")
      writeLines(
        c("# Signature extraction skipped",
          paste0("Cohort size (", n_samples, ") is below the de novo NMF stability threshold (", MIN_SAMPLES_DE_NOVO, " samples)."),
          "De novo extraction was not run to avoid an unstable/overfit factorization (see references/pitfalls.md #2, #3).",
          "The 96-context trinucleotide matrix was written to results/trinucleotide_matrix.csv.",
          "Recommended next step: refit against a fixed COSMIC SBS reference set (e.g. deconstructSigs, SigProfilerAssignment, or mSigAct) rather than de novo extraction."),
        file.path(out_dir, "results", "signature_analysis_skipped.txt")
      )
    } else {
      est <- tryCatch(
        with_null_device(maftools::estimateSignatures(mat = tnm, nTry = 6, pConstant = 0.1, plotBestFitRes = FALSE, parallel = 1)),
        error = function(e) { log("estimateSignatures() failed: ", conditionMessage(e)); NULL }
      )
      if (!is.null(est)) saveRDS(est, file.path(out_dir, "results", "signature_estimateSignatures_result.rds"))

      # maftools does not return a single programmatic "best rank" from
      # estimateSignatures(); rank selection is normally done by visually
      # inspecting the cophenetic-correlation elbow plot below. Use the
      # user-supplied override if given, else fall back to n=3 (the
      # maftools vignette's own low-mutation-count default) and say so.
      n_sig <- if (!is.na(n_signatures_override)) n_signatures_override else 3
      log("Using n = ", n_sig, " signatures for extractSignatures(). ",
          if (is.na(n_signatures_override)) "This is a default, NOT an automatically estimated rank -- inspect figures/signature_cophenetic.pdf and rerun with an explicit n_signatures argument if a different rank fits better." else "(user-specified via n_signatures argument)")
      if (!is.null(est)) {
        save_base_plot("signature_cophenetic", maftools::plotCophenetic(res = est), width = 8, height = 6)
      }

      sig <- tryCatch(
        maftools::extractSignatures(mat = tnm, n = n_sig, pConstant = 0.1, parallel = 1),
        error = function(e) { log("extractSignatures() failed: ", conditionMessage(e)); NULL }
      )
      if (!is.null(sig)) {
        saveRDS(sig, file.path(out_dir, "results", "signature_extractSignatures_result.rds"))
        if (!is.null(sig$contributions)) {
          utils::write.csv(sig$contributions, file.path(out_dir, "results", "signature_exposure.csv"))
        } else {
          log("extractSignatures() result had no 'contributions' field under that name; ",
              "full object saved to results/signature_extractSignatures_result.rds for manual inspection.")
        }

        cmp <- tryCatch(
          maftools::compareSignatures(nmfRes = sig, sig_db = "SBS_v34", verbose = TRUE),
          error = function(e) { log("compareSignatures() failed: ", conditionMessage(e)); NULL }
        )
        if (!is.null(cmp)) {
          saveRDS(cmp, file.path(out_dir, "results", "signature_compareSignatures_result.rds"))
          cos_field <- intersect(c("cosine_similarities", "cosine_similarity", "cosine_sim"), names(cmp))
          if (length(cos_field)) {
            utils::write.csv(as.data.frame(cmp[[cos_field[1]]]), file.path(out_dir, "results", "signature_cosmic_similarity.csv"))
          } else {
            log("compareSignatures() result field name for cosine similarity not found among expected names; ",
                "full object saved to results/signature_compareSignatures_result.rds for manual inspection.")
          }
        }

        save_base_plot("signature_exposure",
                        maftools::plotSignatures(nmfRes = sig, contributions = TRUE, sig_db = "SBS_v34", show_barcodes = TRUE),
                        width = 10, height = 7)
        log("Signature extraction: de novo NMF (n=", n_sig, ") matched to COSMIC SBS_v34 by cosine similarity. ",
            "This is a de novo + reference-comparison workflow, NOT a constrained refit (pitfalls.md #3): ",
            "state this distinction explicitly in methods/report text.")
      }
    }
  }
}

# --- 6. Report / summary / methods notes ------------------------------------
# Build a genuine results-and-conclusions summary (not just a file listing),
# using defensive field-name lookups since exact maftools return-object column
# names were not confirmed by execution in this sandbox (see script header).

LOW_MUTATION_COUNT_THRESHOLD <- 50  # pitfalls.md #2
SMALL_COHORT_THRESHOLD <- 50  # pitfalls.md #8

top_genes_lines <- tryCatch({
  gs <- as.data.frame(gene_summary)
  count_col <- intersect(c("MutatedSamples", "AlteredSamples", "total"), names(gs))
  if (length(count_col) && "Hugo_Symbol" %in% names(gs)) {
    gs_sorted <- gs[order(-gs[[count_col[1]]]), , drop = FALSE]
    top_n_report <- min(5, nrow(gs_sorted))
    precision <- if (n_samples < SMALL_COHORT_THRESHOLD) 0 else 1
    pct <- round(100 * gs_sorted[[count_col[1]]][seq_len(top_n_report)] / n_samples, precision)
    sprintf("- **%s**: %d sample(s) mutated (%s%%)",
            gs_sorted$Hugo_Symbol[seq_len(top_n_report)],
            gs_sorted[[count_col[1]]][seq_len(top_n_report)], pct)
  } else {
    "Gene summary field names did not match expected columns (MutatedSamples/AlteredSamples/total, Hugo_Symbol) for this maftools version; see results/gene_summary.tsv for raw counts."
  }
}, error = function(e) paste0("Could not summarize top mutated genes: ", conditionMessage(e)))

low_mut_note <- tryCatch({
  ss <- as.data.frame(sample_summary)
  total_col <- intersect(c("total", "Total"), names(ss))
  barcode_col <- intersect(c("Tumor_Sample_Barcode"), names(ss))
  if (length(total_col) && length(barcode_col)) {
    low <- ss[ss[[total_col[1]]] < LOW_MUTATION_COUNT_THRESHOLD, , drop = FALSE]
    if (nrow(low) > 0) {
      paste0(nrow(low), "/", nrow(ss), " sample(s) have fewer than ", LOW_MUTATION_COUNT_THRESHOLD,
             " total mutations (", paste(low[[barcode_col[1]]], collapse = ", "),
             "); signature exposure attributable to these samples should be treated as lower-confidence even under cohort-level de novo extraction (pitfalls.md #2).")
    } else NULL
  } else NULL
}, error = function(e) NULL)

tmb_summary_line <- if (!is.null(tmb_res)) {
  tryCatch({
    td <- as.data.frame(tmb_res)
    permb_col <- intersect(c("total_perMB", "total_perMB_log"), names(td))
    if (length(permb_col)) {
      vals <- td[[permb_col[1]]]
      paste0("Median TMB (", permb_col[1], "): ", sprintf("%.2f", stats::median(vals, na.rm = TRUE)),
             " across ", nrow(td), " sample(s); captureSize denominator = ", capture_size_mb, " Mb.")
    } else {
      "TMB result field names did not match expected columns for this maftools version; see results/tmb_summary.tsv for raw values."
    }
  }, error = function(e) paste0("Could not summarize TMB: ", conditionMessage(e)))
} else {
  "TMB step was skipped or failed; see Log section below."
}

sig_summary_line <- if (exists("sig") && !is.null(sig) && !is.null(sig$contributions)) {
  tryCatch({
    contrib <- sig$contributions
    dominant_idx <- apply(contrib, 2, which.max)
    dominant_sig <- rownames(contrib)[dominant_idx]
    tab <- sort(table(dominant_sig), decreasing = TRUE)
    paste0("Cohort-wide dominant de novo signature: **", names(tab)[1], "** is the top contributor in ",
           tab[1], "/", ncol(contrib), " samples (matched to COSMIC SBS_v34 by cosine similarity; see results/signature_cosmic_similarity.csv).")
  }, error = function(e) paste0("Could not summarize signature exposure: ", conditionMessage(e)))
} else if (n_samples < MIN_SAMPLES_DE_NOVO) {
  paste0("Signature de novo extraction was skipped (cohort n=", n_samples, " < ", MIN_SAMPLES_DE_NOVO,
         "); see results/signature_analysis_skipped.txt and use a reference-based refitting tool instead (pitfalls.md #2, #3).")
} else {
  "Signature extraction did not produce a contributions matrix; see Log section below."
}

conclusion_line <- {
  parts <- character(0)
  parts <- c(parts, paste0("Across ", n_samples, " sample(s) and ", nrow(gene_summary), " gene(s) summarized, the mutation and TMB results above describe this cohort's variant landscape."))
  if (!is.null(low_mut_note)) parts <- c(parts, "Some samples have low mutation counts and should be treated cautiously for signature interpretation.")
  parts <- c(parts, "CNV in this cohort (if analyzed separately via the Python quickplot route) is not purity/ploidy-corrected, and signature results are de novo extraction compared to a reference set, not a constrained refit.")
  parts <- c(parts, "Treat all findings here as descriptive and hypothesis-generating; confirm driver/actionability claims against references/pitfalls.md before finalizing conclusions.")
  paste(parts, collapse = " ")
}

writeLines(c(
  "# variant-cnv-signature-analysis: maftools report",
  "",
  paste0("MAF file: ", maf_file),
  paste0("Samples: ", n_samples, " | Genes summarized: ", nrow(gene_summary)),
  paste0("TMB captureSize used (Mb): ", capture_size_mb, " -- state the assay this footprint corresponds to in methods (pitfalls.md #1)"),
  paste0("Reference genome / BSgenome package: ", ref_genome_key, " / ", bsgenome_pkg,
         if (!have_bsgenome) " (NOT installed; signature analysis skipped)" else ""),
  "",
  "## Top mutated genes",
  "",
  top_genes_lines,
  "",
  "## Tumor mutation burden (TMB) summary",
  "",
  tmb_summary_line,
  "",
  "## Mutational signature summary",
  "",
  sig_summary_line,
  if (!is.null(low_mut_note)) c("", paste0("Caveat: ", low_mut_note)) else NULL,
  "",
  "## Limitations",
  "- This script does not perform purity/ploidy-aware CNV calling; use ASCAT/FACETS/PureCN/Sequenza upstream for absolute copy number (pitfalls.md #4, #12).",
  "- Tumor-only vs tumor-normal calling status of the input MAF was not verified here; confirm germline filtering before treating calls as somatic (pitfalls.md #5).",
  "- Signature results (if produced) reflect de novo NMF extraction compared to COSMIC SBS_v34 by cosine similarity, not a constrained reference-set refit (pitfalls.md #3).",
  "- This script was manually reviewed but not executed against real maftools/NMF/BSgenome installs; verify output once in a fully provisioned environment.",
  "",
  "## Conclusion",
  "",
  conclusion_line,
  "",
  "## Log",
  log_lines
), file.path(out_dir, "report.md"))

utils::write.table(
  data.frame(item = c("samples", "genes", "capture_size_mb"),
             value = c(n_samples, nrow(gene_summary), capture_size_mb)),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

run_info <- list(
  maf_file = maf_file,
  ref_genome = ref_genome_key,
  bsgenome_package = bsgenome_pkg,
  capture_size_mb = capture_size_mb,
  n_samples = n_samples,
  n_genes = nrow(gene_summary),
  maftools_version = as.character(utils::packageVersion("maftools"))
)
writeLines(manual_json(run_info), file.path(out_dir, "run", "methods.json"))
