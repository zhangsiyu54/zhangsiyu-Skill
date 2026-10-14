# Common Pitfalls

Check this list before finalizing any proteomics/metabolomics conclusion. Each entry: Trigger (when it shows up) / Mechanism (why it happens) / Symptom (how it looks in the data or output) / Fix (what to do instead).

## 1. Treating MNAR missingness as MCAR/MAR

- Trigger: proteomics or metabolomics matrix with missingness concentrated in low-abundance features or in one condition where a protein/metabolite drops below the detection limit.
- Mechanism: low-abundance ions are stochastically not sampled/fragmented (left-censored, "missing not at random"), which is fundamentally different from missingness that is unrelated to abundance (MCAR) or explained by an observed covariate (MAR). Mean, median, or KNN imputation assumes the missing value would look like the observed distribution, which is false for MNAR.
- Symptom: imputed values sit near the center of the observed distribution instead of the tail, fold changes for features that are "all-or-nothing" between groups are compressed toward zero, and differential feature counts differ sharply depending on which imputation method was used on the same input.
- Fix: inspect the missingness-vs-intensity relationship first (e.g., proportion missing per abundance bin); when missingness rises as abundance falls, treat as MNAR and use left-censored methods (MinProb, QRILC, deterministic minimum/half-minimum) instead of mean/median/KNN. When a mixture is plausible, consider hybrid strategies (e.g., MAR-appropriate method for features missing at random within a condition, MNAR-appropriate method for features missing in an entire condition), and always report which method was used and why (Lazar et al. 2016).

## 2. TMT/iTRAQ ratio compression from reporter-ion interference

- Trigger: isobaric multiplexed (TMT/iTRAQ) quantification performed at MS2 without co-isolation/interference filtering, especially in complex, deep-fractionated samples.
- Mechanism: co-eluting/co-isolated near-isobaric precursor ions contribute background signal to the reporter ion channels, diluting the true reporter ratio toward 1:1. This "ratio compression" is systematic, not random noise, and increases with sample complexity and isolation window width.
- Symptom: true fold changes are underestimated (often several-fold compression for large true ratios), differential features look more modest than the biology suggests, and effect sizes do not reproduce well between MS2-only and MS3/complementary-ion runs of the same sample.
- Fix: prefer MS3/SPS-MS3 or complementary-ion (TMTproC) acquisition when available; if only MS2 data exists, report this as a limitation and avoid stating exact fold-change magnitudes with confidence — favor rank/direction claims. Apply isolation-interference or co-isolation filtering (e.g., excluding PSMs with high precursor co-isolation) before quantification, and never compare TMT fold changes directly to LFQ fold changes from a different platform without noting the compression risk (Ting et al. 2011).

## 3. Conflating LFQ intensity, iBAQ, and spectral counts

- Trigger: a table labeled generically "protein abundance" or "intensity" without stating whether values are raw LFQ intensity, MaxLFQ-normalized intensity, iBAQ, or spectral counts.
- Mechanism: these are different quantification bases. Raw/MaxLFQ intensity is a relative, per-protein comparison across samples; iBAQ divides intensity by the number of theoretical peptides to approximate within-sample absolute abundance/stoichiometry; spectral counts are a much coarser, count-based proxy that saturates for abundant proteins and is unreliable at low abundance.
- Symptom: cross-sample fold-change claims computed on iBAQ values, or abundance-ranking/stoichiometry claims computed on LFQ intensity, or spectral counts used for small fold-change detection — all give misleading conclusions because the metric does not support the claim being made.
- Fix: confirm the exact quantification basis from the source software (MaxQuant, Spectronaut, Proteome Discoverer, FragPipe) before analysis; use MaxLFQ/normalized intensity for between-sample differential comparisons, iBAQ (or riBAQ) only for within-sample relative abundance/stoichiometry statements, and treat spectral counts as semi-quantitative — flag any spectral-count-based fold change as low confidence, especially near the count floor (Cox et al. 2014; Schwanhäusser et al. 2011).

## 4. Uncorrected batch/run effects across MS acquisitions

- Trigger: samples processed across multiple MS runs, injection batches, TMT plexes, or instrument time windows, especially when group and batch are partially or fully confounded.
- Mechanism: instrument sensitivity drift, column aging, and reagent lot differences introduce systematic shifts unrelated to biology; if batch correlates with the group of interest, this shift is indistinguishable from a true biological effect using naive statistics.
- Symptom: PCA/clustering separates samples by run date or plex rather than (or in addition to) biological group; differential features track the batch variable when checked; a technical replicate pair looks less similar than two different biological groups within the same batch.
- Fix: always record and inspect the batch/run column before differential testing. If group and batch are not fully confounded, include batch as a covariate in the linear model (limma/MSstats design) or apply ComBat-style empirical Bayes batch correction; for TMT with multiple plexes, use a bridging/reference channel plus internal reference scaling (IRS) across plexes. If batch and group are fully confounded, say so explicitly and do not claim the group effect is batch-free (Johnson, Li & Rabinovic 2007).

## 5. Peptide-to-protein rollup ambiguity (shared peptides, razor assignment)

- Trigger: protein groups with high sequence homology (paralogs, isoforms, splice variants) where some identified peptides are shared across multiple protein entries.
- Mechanism: a peptide observed in the MS run can belong to more than one protein sequence in the database; naive summation assigns the full peptide intensity to every matching protein, inflating apparent abundance and creating false-positive "detection" of low-abundance paralogs that were never actually distinguished from their more abundant relative.
- Symptom: two paralogous proteins show near-identical intensity profiles and correlated differential results purely because they share most of their peptides; a protein reported as differentially abundant is actually driven by a peptide belonging to a different, unreported isoform.
- Fix: use razor-peptide logic (each shared peptide assigned to the protein group with the most evidence, e.g., MaxQuant "razor + unique peptides") or unique-peptides-only quantification when isoform resolution matters; always report results at the protein-group level as delivered by the search engine rather than re-splitting groups ad hoc, and flag any biological claim that rests on a paralog-specific peptide count that is low.

## 6. Equal-confidence treatment of metabolite identification levels

- Trigger: an untargeted metabolomics feature table where some rows are confirmed against an authentic standard and others are annotated only by accurate mass, MS/MS spectral library match, or in-silico prediction.
- Mechanism: the Metabolomics Standards Initiative (MSI) defines four identification confidence levels (1: confirmed by standard under identical conditions; 2: putatively annotated, e.g. library MS/MS match without standard; 3: putative class; 4: unknown). Downstream pathway or biological interpretation is only as strong as the weakest identification in the chain.
- Symptom: a pathway enrichment or mechanistic narrative is built on a level-3/4 feature as if it were a confirmed metabolite; two "same-named" metabolites from different studies are compared without checking they were identified with comparable confidence.
- Fix: carry an identification-confidence (MSI level) column through the whole pipeline; report it alongside every named metabolite in results tables and figures, and explicitly downgrade the strength of biological claims (hedge language) for level 2–4 features. Confirm level-1 identity for any metabolite that is central to the paper's conclusion (Sumner et al. 2007).

## 7. Normalization method chosen without regard to matrix-specific dilution effects

- Trigger: urine or highly variable-dilution biofluid metabolomics/proteomics data normalized the same way as plasma/serum or cell-lysate data, or a normalization method picked by default rather than by matrix behavior.
- Mechanism: urine concentration varies with hydration and renal function independent of biology, so a single global scaling factor (e.g., total ion count) is easily skewed by a few very abundant, biologically irrelevant compounds (e.g., excreted protein or a dominant osmolyte). Plasma/serum and cell lysates have comparatively stable total protein/metabolite content, so median or quantile normalization is usually adequate there.
- Symptom: after normalization, samples still cluster by apparent dilution (e.g., creatinine or total-intensity proxy) rather than by group; normalized values correlate strongly with a "cheap" proxy like total intensity or specific gravity that should have been removed.
- Fix: for urine (and other high-dilution-variability matrices), prefer probabilistic quotient normalization (PQN) over simple total-intensity/TIC scaling, and check whether a small subset of very abundant features (e.g., urinary protein signal) is driving the PQN coefficient — if so, consider excluding those regions before computing the reference spectrum/quotient. For LFQ proteomics, median or variance-stabilizing normalization (VSN) is generally preferred over quantile normalization when the assumption "most features unchanged" may be violated; validate the choice by checking whether replicate CVs improve after normalization, not just by default habit (Dieterle et al. 2006; Välikangas, Suomi & Elo 2018).

## 8. Multiple-testing burden from many features and few samples

- Trigger: proteomics (thousands of proteins) or untargeted metabolomics (thousands of features) datasets with typical pilot-scale n (3–6 per group).
- Mechanism: testing thousands of features with small per-group sample sizes gives each individual test very low power, while the number of tests still demands FDR correction; the combination produces both a high false-negative rate for true effects and possible false-positive clusters if p-values are not moderated.
- Symptom: raw (uncorrected) p < 0.05 lists that shrink to almost nothing after Benjamini-Hochberg correction; volcano plots that look dramatic on p-value alone but are dominated by features with only 1–2 non-missing values per group; results that fail to reproduce in a validation cohort.
- Fix: use moderated/empirical-Bayes variance shrinkage (limma-style eBayes, as implemented in DEP/MSstats) rather than per-feature t-tests, since borrowing variance information across features partially offsets small-n instability; always report BH-FDR-adjusted q-values as the primary threshold, keep the full tested-feature table (not just "significant" hits) so effect sizes can be re-examined, and state the achieved power/limitations explicitly rather than implying the study is adequately powered (Kammers et al. 2015).

## 9. Assuming linear/direct correspondence in multi-omics integration

- Trigger: proteomics compared directly to transcriptomics (or proteomics to metabolomics) with an implicit expectation that fold changes should agree in direction and magnitude, feature by feature.
- Mechanism: post-transcriptional and post-translational regulation (translation efficiency, protein half-life, degradation, feedback control) decouples mRNA abundance from protein abundance; genome-wide mRNA-protein correlations are typically only moderate, and for any individual gene the two layers can move independently or even in opposite directions. Proteomics-metabolomics relationships are similarly indirect because a single metabolite level reflects the net balance of multiple enzymes, transporters, and consumption.
- Symptom: a gene highlighted as "significant" at the transcript level is absent, unchanged, or moving oppositely at the protein level (or vice versa), and this is reported as a data-quality problem rather than a biological one; a correlation network treats every cross-omics edge as evidence of the same regulatory mechanism.
- Fix: report cross-omics agreement as a rate/fraction, not an expectation of 1:1 concordance, and explicitly discuss discordant genes as potential post-transcriptional/post-translational regulation candidates rather than errors. Prefer pathway- or module-level integration (e.g., MOFA2, DIABLO/mixOmics) over feature-by-feature correlation when sample size allows, and always describe correlational multi-omics findings as associative, not causal, absent independent validation (Schwanhäusser et al. 2011).

## 10. Metabolite-to-pathway ID mapping mismatches

- Trigger: pathway enrichment for metabolomics features run through a gene-centric or loosely mapped database, or feature IDs (KEGG compound C-numbers, HMDB IDs, PubChem CIDs, InChIKeys) mixed without a verified crosswalk.
- Mechanism: pathway enrichment tools (e.g., KEGG-based MSEA in MetaboAnalyst) require exact compound ID matches to their internal pathway compound lists; the same metabolite can have different identifiers across databases, and many-to-one or one-to-many mappings (e.g., stereoisomers collapsing to one KEGG ID, or one HMDB ID splitting across several KEGG entries) silently drop or duplicate features during mapping.
- Symptom: a metabolite known to participate in a pathway is absent from the enrichment result because its ID didn't map; enrichment background list uses "all KEGG compounds" instead of "all compounds actually measured/detectable on this platform," inflating apparent significance; isotopologues/adducts/in-source fragments of the same underlying metabolite are counted as separate hits, artificially boosting a pathway's hit count.
- Fix: verify and record the ID crosswalk (name -> InChIKey -> KEGG/HMDB) explicitly before enrichment, report the mapping/coverage rate (features mapped / total features), use the measured-feature list (not the full pathway database) as the statistical background, and collapse known isotopologue/adduct duplicates of the same compound before counting pathway hits.

## 11. Ignoring qualitative "detected-in-one-group-only" features in fold-change ranking

- Trigger: a feature (protein or metabolite) has essentially all valid values in one condition and all missing values in the other.
- Mechanism: computing log2FC after imputing the missing side with a small constant produces an arbitrarily large but not truly quantitative fold change — the imputed value sets the ceiling, so the reported effect size reflects the imputation choice more than biology.
- Symptom: the top of the volcano plot / ranked list is dominated by "all-or-nothing" features with implausibly large fold changes and very small p-values driven by low within-group variance of the imputed values.
- Fix: flag presence/absence features separately (e.g., a qualitative "on/off" table) from the quantitative differential table; when included in ranked/quantitative lists, cap or annotate the fold change as a lower bound and disclose the imputation-dependent nature of the number rather than presenting it as a precisely measured effect size.

## 12. Untargeted metabolomics feature redundancy (isotopologues, adducts, in-source fragments)

- Trigger: LC-MS untargeted metabolomics feature tables from peak-picking software (e.g., XCMS, MZmine) without adduct/isotope annotation collapsing.
- Mechanism: a single metabolite typically generates multiple detected features — the monoisotopic peak plus 13C isotopologues, multiple adducts ([M+H]+, [M+Na]+, [M+NH4]+), and in-source fragments — each of which is picked up as an independent "feature" by untargeted software.
- Symptom: feature counts vastly exceed the plausible number of distinct metabolites in the sample; several "different" significant features in the differential table turn out to be the same underlying compound, inflating both the apparent hit count and downstream pathway enrichment significance (pseudo-replication of one biological signal across several rows).
- Fix: run adduct/isotope grouping (e.g., CAMERA, xMSannotator, or platform-native grouping) before statistics if identity claims will be made, or at minimum flag correlated same-retention-time feature clusters and report them as one biological entity in interpretation, even if they remain separate rows in the raw statistical table.

## 13. Using default search-engine FDR without accounting for protein-level inference

- Trigger: reporting protein identification/quantification confidence using only the PSM- or peptide-level 1% FDR reported by the search engine.
- Mechanism: PSM-level FDR controls the false-discovery rate among peptide-spectrum matches, but proteins are inferred from potentially many peptides (and shared peptides across homologs, see pitfall 5); the false-discovery rate at the protein-group level is not the same quantity and is typically higher unless separately controlled (protein-level FDR / picked-protein FDR).
- Symptom: a protein list described as "1% FDR" that, on closer inspection, was only peptide/PSM-level filtered, with no separate protein-level FDR control step, especially in single-peptide protein identifications.
- Fix: confirm which level (PSM, peptide, or protein-group) the reported FDR threshold applies to; require a separate protein-level FDR control (e.g., picked-protein FDR, or search-engine protein-level q-value) before treating single- or two-peptide proteins as confidently identified, and consider excluding or flagging single-peptide-supported proteins in the final report.
