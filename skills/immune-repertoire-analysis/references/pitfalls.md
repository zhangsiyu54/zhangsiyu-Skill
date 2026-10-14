# Immune Repertoire Analysis Pitfalls

Read this before finalizing any diversity, expansion, overlap, or clonal lineage conclusion. Each pitfall lists the trigger condition, the underlying mechanism, how it typically shows up in output, and the fix.

## 1. Diversity indices confounded by sequencing depth / cell number

- **Trigger**: Comparing Shannon entropy, Simpson/inverse Simpson, Gini coefficient, clonality (1 - normalized Shannon entropy), or Chao1 richness across samples with unequal read depth, unequal number of recovered cells, or unequal number of unique templates/UMIs.
- **Mechanism**: All of these indices are estimated from a finite sample of an unknown, far larger underlying repertoire. Richness-sensitive metrics (Chao1, observed species, D50) rise monotonically with sampling effort because rarer clones keep being discovered; entropy-based metrics (Shannon, Simpson) are less depth-sensitive at high N but still biased at low N, and clonality inflates artificially in undersampled repertoires because a few high-frequency clones dominate a small sample even when the true repertoire is broad.
- **Symptom**: A sample with 2,000 cells recovered looks "less diverse"/"more clonal" than a sample with 20,000 cells even when the two are biologically identical; diversity rank order flips after downsampling.
- **Fix**: Report cell/read count and unique template count alongside every diversity number. Use rarefaction curves (interpolation) and, where extrapolation is needed, depth-aware estimators such as Chao1/Chao2 with bias correction, Recon (maximum-likelihood diversity estimation with error bars, Chibnik/Kaplinsky & Arnaout), or DivE (Laydon et al., model-averaged extrapolation). At minimum, downsample all samples in a comparison to the same number of cells/templates (e.g., minimum common depth) before computing entropy/richness, and repeat with multiple random subsamples to get a distribution rather than a point estimate.

## 2. Clonotype definition ambiguity changes every downstream number

- **Trigger**: Pipeline or user does not state whether a "clonotype" is defined by CDR3 nucleotide identity, CDR3 amino acid identity, full V(D)J nucleotide sequence, or paired-chain (alpha+beta / heavy+light) identity.
- **Mechanism**: These definitions nest at different resolutions. CDR3-nt is the most granular (separates PCR/somatic silent variants and different V/J calls sharing the same CDR3aa); CDR3-aa collapses convergent nucleotide solutions into one clonotype; paired-chain definitions additionally require both chains to co-occur in the same cell barcode. Switching definitions changes clone counts, singleton fraction, expansion index, and diversity simultaneously — often without any change in the underlying biology.
- **Symptom**: Two analyses of the same dataset report different "top clone frequency" or "number of expanded clones" and neither is wrong, they used different clonotype grain.
- **Fix**: State the clonotype definition explicitly in the methods text before reporting any number (this skill's cns-method-guide.md and paper-practice.md already require this). When comparing to a prior publication or a second dataset, verify both used the same definition before treating the comparison as apples-to-apples. Prefer CDR3-aa + V/J gene (not raw CDR3-nt) as the default "clonotype" for cross-sample biological comparison, since nt-level identity is sensitive to PCR/sequencing noise unless UMI error-corrected.

## 3. Multiplet/doublet contamination fabricates dual-TCR "clonotypes" in single-cell data

- **Trigger**: scTCR/scBCR (10x Genomics or similar droplet-based) data shows cells with 3+ productive contigs, or an unexpectedly high fraction of "dual-clonotype" or multi-chain cells.
- **Mechanism**: Droplet-based capture occasionally co-encapsulates two cells (a true doublet) or captures ambient mRNA/contigs leaked into a droplet, producing a barcode with chains from more than one lymphocyte. This is distinct from genuine biological dual-TCR expression from incomplete allelic exclusion, which does occur in ~10-20% of normal T cells (thymus/blood) but should show one dominant alpha+beta pair, not arbitrary chain combinations.
- **Symptom**: Clonotype table contains many low-frequency 3- or 4-chain "clonotypes" that inflate apparent clonal diversity and richness, or a large expanded clonotype turns out to be an artifact of ambient contamination rather than true clonal expansion.
- **Fix**: Apply the assay pipeline's built-in doublet/multiplet filtering (e.g., Cell Ranger VDJ's doublet test, which removes exact subclonotypes consistent with multiplets) before downstream analysis. Report the fraction of cells excluded. Do not silently keep all raw contigs; state whether multi-chain cells were collapsed to the dominant alpha-beta (or heavy-light) pair by UMI count, or excluded entirely, and be consistent about this choice across all samples being compared.

## 4. Cross-sample comparison without normalizing for cell number / library size

- **Trigger**: Comparing raw clone counts, raw clonotype richness, or overlap counts between samples that differ substantially in cells recovered, sequencing depth, or input material (e.g., tumor biopsy vs. peripheral blood).
- **Mechanism**: Almost every repertoire summary statistic (richness, overlap size, number of expanded clones) scales with sampling effort. A sample with 10x more cells will show more unique clonotypes and more pairwise overlap by chance alone, independent of any biological difference.
- **Symptom**: "Group A has significantly higher clonal diversity than Group B" turns out to be fully explained by Group A having 5x the cell recovery.
- **Fix**: Always normalize before comparing: convert counts to frequencies, downsample to a common depth/cell number, or use overlap/diversity metrics designed to be sample-size robust (e.g., Morisita-Horn index for overlap, which is less sensitive to sample size than raw shared-clone counts; Chao-Jost Hill-number framework for diversity at a common order q). Report both raw N and normalized metric side by side so a reviewer can audit the correction.

## 5. Convergent recombination distorts expansion/diversity depending on clonotype grain

- **Trigger**: Reporting "clonal expansion" or "diversity" using a single clonotype definition without checking whether identical CDR3 amino acid sequences arise from multiple distinct nucleotide rearrangements or V gene calls.
- **Mechanism**: V(D)J recombination is biased — the same CDR3aa sequence can be produced independently by different N-addition/deletion events, sometimes using different V or J genes ("convergent recombination", first characterized by Venturi et al. 2006, PNAS). If clonotyping is done at the amino-acid level only, independently-generated (non-clonal) sequences get merged into one "clone," inflating apparent expansion and deflating apparent diversity. If clonotyping is done at the nucleotide level only, truly clonal descendants that accumulated silent mutations or were called against different reference alleles can be split into separate "clones," deflating apparent expansion.
- **Symptom**: A "highly expanded" clone by CDR3aa criteria shows 2+ distinct V genes or nucleotide sequences on inspection, meaning it is not actually a single expanded clone; or, conversely, plainly identical biological clones are reported as separate singletons at the nucleotide level.
- **Fix**: Report clonotype definition and, when expansion claims matter (e.g., putative antigen-driven clones), verify that top "expanded" aa-level clonotypes share consistent V/J usage and low nucleotide diversity before calling them a true clone. Flag CDR3aa sequences known to be recombination-biased/public (short CDR3, few N-additions, germline-encoded like "public" motifs) as candidates for convergent recombination rather than genuine expansion.

## 6. Public/private overlap interpreted without a convergent-recombination null model

- **Trigger**: Reporting that clonotypes are "shared" (public) between unrelated individuals/samples and interpreting this as evidence of a common antigen-driven response.
- **Mechanism**: Some CDR3 sequences are generated at high probability by the recombination machinery regardless of antigen exposure (short CDR3, minimal N-nucleotide additions, frequently-used V/J pairs). These "public" sequences will recur across unrelated repertoires purely by chance at a rate that depends on repertoire size and generation probability, not shared immunological history.
- **Symptom**: Overlap between two completely unrelated cohorts (e.g., different diseases, different species-matched controls) is nonzero and gets over-interpreted as biologically meaningful sharing.
- **Fix**: Before interpreting overlap as biologically meaningful, compare observed overlap to a background/null expectation (e.g., permutation of samples, or generation-probability models such as OLGA/IGoR-based pgen estimates) and account for repertoire size in the comparison (larger sampled repertoires have higher expected chance overlap). State explicitly whether the overlap metric used (Jaccard, Morisita-Horn, raw shared-clone count) is sequence-level or clonotype-ID-level, since ID-level overlap from different pipelines is not directly comparable.

## 7. BCR clonal lineages require lineage reconstruction, not exact-match clonotyping

- **Trigger**: Applying a TCR-style exact CDR3 (nt or aa) clonotype match to BCR/Ig data to call "the same B cell clone" across sequences or timepoints.
- **Mechanism**: Unlike T cells, B cells undergo somatic hypermutation (SHM) and affinity maturation after initial V(D)J recombination, so descendants of a single naive B cell clone diverge in CDR3 (and framework) sequence over time. Exact-match clonotyping will therefore split a single clonal lineage into many "different" clonotypes, systematically undercounting true clonal expansion and overstating diversity for BCR data.
- **Symptom**: BCR diversity looks implausibly high, and no expanded clones are detected even in a known germinal-center or vaccine-response sample.
- **Fix**: Define B cell clonal families by lineage-aware criteria (same V and J gene, same CDR3 length, and sequence similarity above a threshold, e.g., hierarchical clustering with a chosen distance cutoff) using tools built for this purpose (Change-O `DefineClones`, partis, SCOPer), then reconstruct lineage trees (IgPhyML, dowser) with germline-corrected reference alleles (TIgGER) before making any statement about clonal expansion, selection pressure, or affinity maturation.

## 8. V/J gene usage bias from library preparation method, not biology

- **Trigger**: Comparing V/J gene usage frequencies across samples or studies processed with different library prep protocols (multiplex PCR with V-gene-specific primers vs. 5' RACE/template-switching vs. capture-based enrichment).
- **Mechanism**: Multiplex PCR primer sets have unequal amplification efficiency across V/J gene primers, systematically over- or under-representing specific genes independent of the true repertoire. 5' RACE and UMI-based methods reduce but do not eliminate this bias (ligation/reverse-transcription efficiency still varies by transcript). Comparing raw usage frequencies across differently-prepared datasets attributes a technical artifact to biology.
- **Symptom**: "Differential V gene usage" between two cohorts turns out to correlate with which sequencing center/kit processed each cohort, not with disease status.
- **Fix**: Only compare V/J usage across samples processed with the same library prep chemistry and primer set. If cross-protocol comparison is unavoidable, note it as a caveat, restrict comparison to genes known to amplify comparably, or apply a documented amplification-bias correction fit from spike-in/synthetic repertoire controls.

## 9. Small cell numbers per sample produce unstable metrics reported with false precision

- **Trigger**: Diversity, clonality, or overlap statistics computed and reported to 2-3 decimal places for samples with only tens to a few hundred cells recovered (common in FACS-sorted rare populations or low-input scVDJ).
- **Mechanism**: With few observations, richness and evenness estimators have wide sampling variance; a single additional expanded clone observed or missed can move the diversity index substantially. Point estimates without confidence intervals give a false impression of precision.
- **Symptom**: A table reports "Shannon diversity = 4.312 vs. 4.298" for two 80-cell samples and treats the difference as meaningful.
- **Fix**: State the cell/template count next to every diversity number. For samples below a few hundred unique templates, compute bootstrap confidence intervals (resample with replacement) or jackknife estimates rather than reporting a single point value, and avoid ranking or statistically testing samples whose CIs overlap substantially. Prefer categorical description ("low cell number, diversity metrics unstable") over spurious rank ordering when N is very small.

## 10. Alluvial/Sankey clone-tracking plots implying unsupported lineage causality across timepoints

- **Trigger**: Building an alluvial/Sankey/circos plot that connects the "same" clonotype across timepoints or tissues and describing it in prose as a clone "migrating," "expanding," or "differentiating" between states.
- **Mechanism**: Alluvial plots visualize co-occurrence of a clonotype label across categories/timepoints; they do not, by themselves, establish a causal or even temporally ordered biological relationship. Independent draws from a large, overlapping public/convergent clonotype pool can produce the same visual pattern without any true lineage relationship, especially for CDR3aa-level clonotype definitions (see pitfall 5).
- **Symptom**: A results narrative states "clone X migrated from blood to tumor and differentiated into an exhausted phenotype" based solely on the clonotype ID reappearing in both compartments, with no statistical test of enrichment above chance and no single-cell trajectory evidence.
- **Fix**: Limit alluvial/tracking plots to top clonotypes with adequate supporting cell numbers per node, state the clonotype definition used, and restrict causal/directional language to cases supported by paired longitudinal sampling design, statistical enrichment tests (e.g., permutation test for shared-clone counts exceeding chance), or independent lineage evidence (BCR: shared lineage tree; TCR: paired scRNA trajectory).

## 11. UMI/PCR duplicate handling errors distort clone size estimates

- **Trigger**: Amplicon-based repertoire sequencing without unique molecular identifiers (UMIs), or UMI-based data where UMI collapsing/error-correction is skipped or misconfigured.
- **Mechanism**: PCR amplification is exponential and stochastic ("PCR jackpotting"); without UMIs to trace reads back to original template molecules, read counts per clonotype reflect amplification efficiency and cycle number as much as the true number of input molecules/cells. Even with UMIs, incorrect deduplication (e.g., not collapsing UMIs with sequencing errors) can artificially split or merge molecule counts.
- **Symptom**: Clone size rank order or frequency shifts substantially between technical replicates of the same sample, or clone frequency does not track cell number in matched single-cell data.
- **Fix**: Prefer UMI-based library prep for frequency-sensitive analyses (clonality, expansion ranking) and confirm the pipeline performs UMI error correction (e.g., directional/network-based UMI collapsing, not naive exact-match). When only read-count (non-UMI) data is available, treat frequency estimates as approximate and disclose this in methods; do not present read-count-derived frequencies with the same confidence as UMI- or cell-count-derived frequencies.

## 12. Clonally related cells treated as statistically independent observations

- **Trigger**: Running differential abundance, differential expression, or other per-cell statistical tests across groups/conditions on scRNA+VDJ data without accounting for clonal membership.
- **Mechanism**: Cells belonging to the same expanded clone are not independent draws — they share a common TCR/BCR and often correlated transcriptional states (having descended from one activated/proliferating parent cell). Standard per-cell tests (e.g., naive t-test/Wilcoxon on a marker gene across all cells in a condition) implicitly assume independence, so a single large expanded clone can dominate the test statistic and produce an apparently significant result driven by one clonal family, not by the condition.
- **Symptom**: A "significant" transcriptional or compositional difference between groups disappears or reverses when the analysis is re-run at clonotype level (one representative/aggregate value per clone) rather than cell level.
- **Fix**: Report clone size distribution alongside any per-cell statistical comparison; consider clone-aware or pseudobulk-by-clone aggregation, mixed-effects models with clonotype as a random effect, or explicitly note when a result is driven by one or a few dominant clones and verify robustness by excluding them.
