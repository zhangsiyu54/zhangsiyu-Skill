# Common Pitfalls and Failure Modes

Read this file before finalizing conclusions or writing the methods section. Each entry gives the input/request that triggers the pitfall, why it happens, what it looks like when it happens, and the correct fix. These are recurring, concrete failure modes seen in real bulk RNA-seq DE work — not generic caveats.

## 1. Reference level flips the sign (and meaning) of log2FC

- **Trigger:** Group column is a plain character vector or an unordered factor (e.g. `condition = c("treated","control",...)`) passed to DESeq2/edgeR/limma without an explicit reference level.
- **Mechanism:** R's default factor level ordering is alphabetical. DESeq2/limma pick the first level as the reference (denominator) unless told otherwise, so `"control"` (alphabetically after `"Control"`/`"Treated"` depending on capitalization) can silently become the numerator.
- **Symptom:** log2FC direction is reversed relative to what the user intended (upregulated genes reported as downregulated); "up in treatment" volcano plot actually shows genes up in control.
- **Fix:** Always set `relevel(factor(condition), ref = "control")` (or `factor(condition, levels = c("control","treated"))`) explicitly before building the DESeqDataSet/design, and state the reference level in the methods text and report.

## 2. `results()` called with no `name`/`contrast` returns the wrong coefficient

- **Trigger:** Multi-factor design (`~ batch + condition`) or any design with more than one non-intercept coefficient; the user asks for "the treatment effect" without specifying which comparison.
- **Mechanism:** `DESeq2::results()` with no arguments returns the **last variable in the design formula's last coefficient** by default, not necessarily the biologically intended contrast. In edgeR/limma, `topTable()`/`glmLRT()` without an explicit `coef=`/contrast likewise defaults to the last model term.
- **Symptom:** Reported DEGs correspond to a batch effect, an unintended factor level comparison, or the wrong pairwise contrast — numbers look "plausible" but do not answer the stated question.
- **Fix:** Always call `results(dds, contrast = c("condition","treated","control"))` or `results(dds, name = "condition_treated_vs_control")` explicitly; run `resultsNames(dds)` first and confirm the coefficient name matches the intended comparison before reporting.

## 3. Mixing shrunken log2FC with unshrunken test statistics/p-values

- **Trigger:** User wants both a clean volcano/ranking (shrunken LFC) and the significance calls from the original `results()` object.
- **Mechanism:** `lfcShrink()` (apeglm/ashr/normal) only replaces the effect-size estimate and its SE; by design it does **not** recompute the Wald statistic or p-value from the original `results()` call unless `lfcShrink(..., res=dds_res)` is used consistently downstream. Treating the shrunken LFC as if it generated the reported p-value is an internal inconsistency.
- **Symptom:** Genes with strong shrunken LFC but the p-value/padj table pulled from an unshrunken `results()` call (or vice versa), producing a report where the effect size and the significance test appear to come from different models.
- **Fix:** Compute `res <- results(dds, contrast=...)`, then `res_shrunk <- lfcShrink(dds, res=res, type="apeglm")`, and report shrunken LFC alongside the p-value/padj **from the same `res_shrunk` object** (padj is unchanged by shrinkage, only LFC/SE are). Explicitly say "log2FC values are apeglm-shrunken; p-values/padj are from the Wald test" in methods.

## 4. Independent filtering / Cook's distance sets padj = NA for genes that look biologically real

- **Trigger:** A gene of known biological interest shows `padj = NA` in the DESeq2 results table.
- **Mechanism:** Two distinct DESeq2 mechanisms cause this: (a) independent filtering by mean normalized count sets padj to NA below an automatically chosen `baseMean` threshold to improve power (Bourgon et al. 2010); (b) Cook's distance outlier detection sets padj to NA when a single sample dominates a gene's within-group variance (a potential outlier), especially with small n and >3 samples per group.
- **Symptom:** A user-nominated candidate gene is missing from the significant list with `padj = NA` instead of a numeric value, and the agent risks either silently dropping it or misreporting it as "not significant" (NA is not the same as non-significant).
- **Fix:** Distinguish the two causes: check `results(dds)@metadata$filterThreshold` and the gene's baseMean (low baseMean -> independent filtering) versus per-sample Cook's distances via `assays(dds)[["cooks"]]` (large outlier value in one sample -> outlier flag). Report NA genes separately with the reason; do not silently drop them from delivered tables. Consider `cooksCutoff=FALSE` or `independentFiltering=FALSE` only with explicit justification, and note the power/type-I-error tradeoff.

## 5. LRT reports an omnibus test but its logged LFC is only the last model coefficient

- **Trigger:** Multi-level factor (e.g. time course with 4+ levels, or a 3+ level treatment) and the user asks "which genes change at all across groups."
- **Mechanism:** `DESeq2` likelihood-ratio test (`test="LRT"`, `reduced = ~1` or a reduced model) correctly gives an omnibus p-value for "does this gene's expression depend on the factor at all," but the `log2FoldChange` column in that same results table is only the Wald estimate for one specific coefficient (typically the last factor level vs. reference) — it is not an average or maximum effect across levels.
- **Symptom:** A gene flagged as highly significant by the LRT p-value shows a deceptively small (or even near-zero) log2FoldChange, because the reported LFC column reflects a single pairwise contrast that doesn't capture where the real variation lies (e.g. change is between two middle time points).
- **Fix:** Use the LRT p-value/padj strictly for significance ("does this gene respond to the factor"), and separately report individual pairwise Wald contrasts (or a per-level LFC matrix) to characterize the direction/magnitude of the effect. Never present the LRT-table LFC column as "the effect size" without this caveat.

## 6. Using tximport/salmon-derived offsets incorrectly, or reverting to a plain count matrix

- **Trigger:** Input is transcript-level quantification (Salmon/kallisto), and the pipeline either skips tximport or uses `countsFromAbundance="scaledTPM"`/`"lengthScaledTPM"` when raw counts + offset were expected.
- **Mechanism:** Gene-level counts naively summed from transcript counts ignore that average transcript length can change between conditions (differential isoform usage), biasing counts independent of true expression change. `tximport()` with `countsFromAbundance="no"` produces a counts matrix plus a matched average-transcript-length offset matrix that DESeq2/edgeR use via `DESeqDataSetFromTximport()`/`normMat`, correcting for this. Using `scaledTPM`/`lengthScaledTPM` instead produces counts already length-normalized, which is appropriate for edgeR/limma-voom without the offset but must not be treated as raw counts elsewhere.
- **Symptom:** Spurious DE calls for genes whose isoform usage (not overall transcription) differs between conditions; or a "raw counts" label applied to a matrix that was already length/TPM-normalized upstream.
- **Fix:** Confirm which quantifier and tximport mode produced the matrix; for DESeq2 use `tximport(..., countsFromAbundance="no")` + `DESeqDataSetFromTximport()` so the length offset is used; for edgeR/limma use `tximport(..., type="salmon")` with `edgeR::catchSalmon()`/`tximport`'s recommended offset workflow. State the quantifier, tximport mode, and offset handling explicitly in methods.

## 7. Median-of-ratios / TMM normalization breaks down when most of the transcriptome shifts

- **Trigger:** Conditions expected to alter a very large fraction of expressed genes: viral infection with host shutoff, heat-shock/global stress response, strong prokaryotic stringent response, or a global perturbation like a general transcription/translation inhibitor.
- **Mechanism:** DESeq2's median-of-ratios and edgeR's TMM both assume most genes are not differentially expressed, so the per-sample size factor/normalization factor is estimated from a "stable majority." When >~50% of genes truly shift (especially in one direction), this assumption is violated and the estimated normalization factor is biased, compressing or inflating fold changes and can flip the direction of "control" genes.
- **Symptom:** Housekeeping genes or spike-in controls show implausible fold changes; MA plots show off-center bulk shift rather than a symmetric cloud around log2FC=0; global "downregulation" or "upregulation" bias across nearly all genes.
- **Fix:** Flag the biological scenario explicitly and consider alternative normalization: spike-in-based normalization (ERCC), total-RNA/rRNA-depletion-based scaling, or a user-supplied normalization factor from an independent assay (e.g. RT-qPCR of a validated stable reference). Report the assumption and why standard normalization may be unreliable rather than silently trusting default size factors.

## 8. Batch correction included for visualization but not in the DE model (or vice versa)

- **Trigger:** Data has a known batch/lab/sequencing-run variable and the user wants "batch-corrected" plots and DE results.
- **Mechanism:** Two distinct, commonly confused workflows: (a) including `batch` as a covariate in the DE design formula (`~ batch + condition`) correctly accounts for batch when testing `condition`, preserving proper variance/uncertainty; (b) applying `ComBat`/`limma::removeBatchEffect()` to produce a "batch-removed" expression matrix and then running DE testing (or even just computing p-values) directly on that adjusted matrix artificially removes variance and inflates significance — Nygaard, Rødland & Hovig (2016) demonstrated this produces exaggerated confidence when the group of interest is not perfectly balanced across batches.
- **Symptom:** Suspiciously clean PCA/heatmap after ComBat but DE results that show far more significant genes than the batch-in-design-formula approach on the same data; or, conversely, a design formula with batch included but PCA/heatmap plots still showing raw batch-driven clustering that confuses interpretation.
- **Fix:** For hypothesis testing, always include batch as a covariate in the design formula and let DESeq2/edgeR/limma model it jointly with condition — never test on ComBat-adjusted values. `ComBat`/`removeBatchEffect()` output is for **visualization only** (PCA, heatmaps, clustering); label these plots explicitly as "batch-adjusted for visualization only, not used for statistical testing."

## 9. Paired/repeated-measures samples analyzed as if independent

- **Trigger:** Design includes matched pairs (pre/post treatment in the same patient, tumor/normal from the same individual, longitudinal timepoints per subject) but the design formula only includes `condition`.
- **Mechanism:** Ignoring the subject/pair structure treats within-subject correlation as if it were between-subject variance, which both wastes statistical power (failing to remove consistent inter-subject baseline differences) and can produce invalid standard errors.
- **Symptom:** Fewer DEGs than expected given clear paired trends in individual gene plots; PCA/heatmap shows samples clustering by subject rather than by condition, but the model doesn't account for it.
- **Fix:** Add the subject/pair identifier to the design formula (`~ patient + condition` in DESeq2/limma) so each subject acts as its own blocking factor; for limma with more complex random-effect-like structure, use `duplicateCorrelation()`. Confirm the pairing column exists and is complete (no missing pairs) before modeling.

## 10. Underpowered designs (n=1, n=2 per group) reported as if inferentially valid

- **Trigger:** Only 1-2 replicates per group, or the user asks for "significant genes" from a pilot/exploratory dataset.
- **Mechanism:** DESeq2/edgeR need enough replicates per group to estimate gene-wise dispersion reliably; with n=2, dispersion estimates are almost entirely dependent on the shared trend/empirical Bayes shrinkage across genes rather than actual within-group replication, making individual gene calls fragile. With n=1 per group, no within-group variance exists at all and DESeq2 will refuse or require `--no-replicates`-style workarounds that assume Poisson-only technical variance (unrealistic for biological samples).
- **Symptom:** Tool runs and produces a results table (possibly with a warning), giving a false impression of statistical rigor; large numbers of "significant" genes driven entirely by shrinkage assumptions rather than genuine replication.
- **Fix:** Explicitly flag n<3/group as underpowered for confident biological inference; report fold-change/ranking as exploratory only, and state in methods and limitations that dispersion could not be reliably estimated from biological replicates. Do not present n=1-per-group comparisons as formal differential expression.

## 11. Feeding already-normalized/microarray-style continuous values into count-based tools

- **Trigger:** Input matrix contains TPM/FPKM/RPKM, log-transformed values, or microarray intensities, but the user (or a previous pipeline step) requests DESeq2/edgeR.
- **Mechanism:** DESeq2 and edgeR model raw integer counts with a negative binomial distribution; non-integer, already-normalized, or log-scaled values violate this likelihood model entirely (the tools may error, silently round, or produce meaningless dispersion estimates).
- **Symptom:** DESeq2 throws an error about non-integer values, or (worse) coerces/rounds values and produces results that look plausible but rest on an invalid statistical model.
- **Fix:** Detect non-integer/negative/already-scaled inputs before modeling; for TPM/FPKM/microarray-style continuous data, use `limma` with `lmFit()`+`eBayes()` (moderated t-test on log-scale values) instead, and state clearly that DESeq2/edgeR were not used because the input was not raw counts.

## 12. Reporting padj-only "significance" without an effect-size floor

- **Trigger:** Very large sample size per group (common in public cohort re-analysis, e.g. n>50), or extremely deep sequencing.
- **Mechanism:** Statistical power to detect small effect sizes as significant scales with sample size; with large n, genes with biologically trivial log2FC (e.g. 0.05) can reach padj<0.05 purely from precision, not biological relevance.
- **Symptom:** Very long "significant" gene lists dominated by tiny fold changes; downstream enrichment analysis diluted by noise-level effects.
- **Fix:** Always pair the adjusted p-value threshold with an explicit effect-size threshold (e.g. |log2FC|>=1) as already required by this skill's default thresholds, and call this out explicitly when n is large, since padj alone becomes a weak filter at scale.

## 13. Confounded design / rank-deficient model matrix

- **Trigger:** A covariate (e.g. batch) is perfectly or nearly perfectly correlated with the condition of interest (e.g. all "treated" samples were run in batch 2, all "control" in batch 1).
- **Mechanism:** If batch and condition are fully confounded, the design matrix is not full rank (or nearly so), so the model cannot separate the batch effect from the condition effect — DESeq2/limma will either error ("model matrix not full rank") or return coefficients with enormous, meaningless standard errors.
- **Symptom:** Tool errors on `DESeq()`/`lmFit()`, or runs but produces absurdly large/unstable log2FC and SE values for the condition coefficient.
- **Fix:** Check the cross-tabulation of every covariate against the condition of interest before modeling (`table(batch, condition)`); if fully confounded, state explicitly in the report that batch and condition cannot be statistically disentangled with this design, and do not attempt to "adjust" for a fully confounded batch term.

## 14. Multiple contrasts from one model reusing padj across comparisons

- **Trigger:** A single multi-factor model is used to extract more than one pairwise contrast (e.g. `treated_A_vs_control`, `treated_B_vs_control`) and the user wants one combined "significant gene" table.
- **Mechanism:** Each call to `results(dds, contrast=...)` performs its own independent Benjamini-Hochberg correction across all genes for that specific contrast. Padj values from different `results()` calls are not on a jointly corrected scale, so combining/comparing them as if drawn from one FDR-controlled family is invalid.
- **Symptom:** A gene appears "significant" in contrast A's table and "not significant" in contrast B's table, and the two padj columns get merged into a single "combined significance" column without accounting for the fact that each was corrected separately.
- **Fix:** Report and threshold each contrast's padj independently; when a combined multi-contrast summary is needed, either use a joint LRT for the overall question or clearly label each contrast's FDR as computed within that contrast only, not across contrasts.
