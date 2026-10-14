# CRISPR Screen Pitfalls

Check this list before finalizing any hit call, ranking, or written conclusion. Each entry: Trigger (when it applies), Mechanism (why it happens), Symptom (how it shows up in the data), Fix (what to do).

## 1. Library representation collapses before selection even starts
- Trigger: Transduction/puromycin selection performed at low cells-per-sgRNA, or early passage bottleneck before T0 sampling.
- Mechanism: Stochastic loss of low-abundance sgRNAs during a narrow cell-number bottleneck is indistinguishable from true biological dropout; sub-~200-500x coverage per sgRNA (library-size dependent, often cited as >=300-500x for genome-wide libraries and often higher for larger libraries) sharply increases this risk.
- Symptom: High Gini index or skewed lognormal count distribution already at T0/plasmid vs. cell-pellet comparison, before any treatment; many sgRNAs with count 0 or near-0 in every replicate including control.
- Fix: Always inspect T0 (or plasmid pool) representation independent of the treatment/control contrast; if T0 coverage is degraded, flag the whole screen as underpowered rather than trusting downstream depletion calls. Do not silently proceed to hit calling.

## 2. Positive- and negative-selection hits need different statistical handling
- Trigger: Screen has both an enrichment arm (dropout survivors / resistance) and a depletion arm (essentiality / sensitization) in the same contrast.
- Mechanism: MAGeCK RRA aggregates ranks separately for each direction and is comparatively conservative for enrichment (positive selection) because a handful of strongly outlying sgRNAs can dominate rank aggregation, while MLE's beta score and per-gene NB GLM can behave differently for genes with only 1-2 strongly responding guides. BAGEL/BAGEL2 is built around a depletion-only reference-set model and is not appropriate for scoring positive selection at all.
- Symptom: A resistance-gene hit list that looks implausibly short compared to the essentiality hit list from the same run, or a "hit" that appears only in one algorithm's default direction settings.
- Fix: State which direction was tested and run the tool in the matching mode (`--gene-lfc-method`/`neg`/`pos` selection in `mageck test`, or a BAGEL2 depletion run plus a separate RRA/MLE positive-selection pass). Never let a negative-selection-tuned tool report positive-selection conclusions without checking its assumptions.

## 3. Essential-gene calls without a real non-targeting/safe-harbor baseline
- Trigger: Screen lacks non-targeting controls (NTCs), safe-harbor-targeting guides, or an established nonessential/essential reference gene set (e.g., Hart/DepMap common-essential and nonessential lists).
- Mechanism: Without a null distribution anchored to guides/genes known a priori to have no fitness effect, any depletion threshold is arbitrary; cutting toxicity itself (Cas9 double-strand-break response) can create a floor effect that mimics "essentiality" for many genes at once.
- Symptom: An unusually large fraction of the genome (well above the ~10% commonly reported as core-essential in DepMap-style pan-cancer screens) called "essential"; z-score or beta-score distributions that are not clearly bimodal.
- Fix: Require NTC or reference essential/nonessential gene sets before making any essentiality claim; if absent, explicitly label results as an uncalibrated rank list rather than a validated essential-gene set, and say so in the methods/limitations text.

## 4. Low-count sgRNAs are ranked as if they were reliably measured
- Trigger: Raw sgRNA counts contain guides with very low reads (e.g., single digits) in the control/T0 sample.
- Mechanism: A guide already near the Poisson/negative-binomial noise floor in the reference sample will show a large, unstable fold-change with a small absolute change in reads; this is a library-complexity/sampling artifact, not biological depletion.
- Symptom: Top-ranked "hits" driven by fold-changes computed from denominators of <20-30 reads; volatile ranking when a pseudocount or normalization method is changed.
- Fix: Apply an explicit low-count filter (e.g., minimum control-sample read count, or minimum total reads across replicates) before hit calling, and report the threshold used. Confirm the tool's built-in filtering (e.g., MAGeCK's default handling) was not silently overridden.

## 5. T0 vs. late-timepoint contrasts conflate selection with generic passage/fitness drift
- Trigger: Comparison is framed as "treatment vs. control" but is actually "early timepoint vs. late timepoint," with or without a drug/treatment arm.
- Mechanism: Extended passaging imposes a general proliferation/fitness selection (culture adaptation, media stress, contact inhibition) on top of any treatment-specific effect; genes affecting general fitness will appear as hits in every arm, treatment or not.
- Symptom: Large overlap between "treatment hit list" and "vehicle/DMSO-only hit list" when both are compared back to T0; a long list of common-essential genes dominating a supposedly drug-specific hit list.
- Fix: Always run the matched vehicle/untreated-vs-T0 contrast alongside the treatment-vs-T0 contrast and subtract or model out the shared fitness signal (e.g., treatment-vs-vehicle contrast, or MLE with a design matrix including a time/passage covariate) before claiming treatment specificity.

## 6. Gene-level "hit" driven by a single sgRNA out of several
- Trigger: A gene has multiple independent sgRNAs (commonly 4-10 per gene) but only one or two show the phenotype.
- Mechanism: Individual sgRNA efficiency varies substantially (on-target cutting efficiency, chromatin accessibility, off-target activity); a single highly active or off-target guide can produce a strong signal that RRA's rank aggregation or an unweighted mean can still surface as gene-significant, especially with permissive parameters.
- Symptom: Gene-level FDR is significant but the sgRNA-level fold-change plot shows 1 of 4-6 guides carrying essentially all the signal, with the rest near zero.
- Fix: Always pull and plot the sgRNA-level summary (not just the gene summary) for candidate hits; require concordant direction in at least half the guides (or the tool's own guide-consistency diagnostic, e.g., MAGeCK's `--gene-lfc-method` or CRISPRcleanR's normLRT) before treating a gene as high-confidence.

## 7. FDR corrected within each contrast but not across the full set of comparisons reported
- Trigger: Multiple contrasts are tested (multiple drugs, multiple timepoints, multiple cell lines) and each is FDR-corrected independently before being combined into one results table or headline claim.
- Mechanism: Per-comparison FDR controls the false discovery rate within that single test set; it does not control the family-wise or overall false-discovery rate across N independently-corrected comparisons, so the combined table is more permissive than its stated threshold implies.
- Symptom: A "hits across all conditions" summary table where the total false-positive burden is understated relative to the per-column FDR values shown.
- Fix: State clearly that FDR was calculated per contrast; if a combined or "hit in >=1 condition" claim is made, note the inflation risk explicitly (or apply an additional correction, e.g., Bonferroni across the number of independent contrasts) rather than presenting per-contrast FDRs as if they applied to the pooled claim.

## 8. Genome-scale screen hits treated as validated "essential/target genes" without orthogonal follow-up
- Trigger: Any mechanistic or biological claim ("Gene X is essential for Y", "Gene X is a drug target") built solely from pooled-screen ranking output.
- Mechanism: Pooled screens have false-positive sources (off-target guides, copy-number artifacts, passenger mutations in the cell line, batch effects) that a single genome-wide pass cannot itself rule out; screen output is a prioritized hypothesis list, not a validated causal claim.
- Symptom: Report language stating a gene "is required for" or "drives" a phenotype based only on screen rank/FDR, with no arrayed knockout, rescue, or secondary-screen citation.
- Fix: Explicitly frame screen output as hypothesis-generating; require arrayed validation (individual knockout/knockdown + rescue, competition assay, or an independent orthogonal screen/library) before mechanistic wording, and say so in the deliverable's limitations section (see `references/deliverables.md`).

## 9. CRISPRi/CRISPRa screens analyzed with knockout-screen assumptions
- Trigger: dCas9-KRAB (CRISPRi) or dCas9-VPR/SAM (CRISPRa) screen data run through a pipeline/parameter set tuned for Cas9 knockout screens.
- Mechanism: CRISPRi/CRISPRa efficacy depends heavily on guide position relative to the TSS (typically most effective within roughly -50 to +300 bp of the TSS for CRISPRi, similar or slightly different optimal windows for CRISPRa) rather than exon-disruption logic used for knockout guide design; per-guide efficiency variance is often larger, and there is no simple "frameshift = null allele" assumption. Off-target/spurious binding at unintended promoters is also a distinct failure mode from Cas9 indel off-targeting.
- Symptom: Guide-level inconsistency far exceeds what is typical for knockout screens of the same gene; strong hits at genes whose guides cluster far from the annotated TSS.
- Fix: Confirm the library was designed for the correct modality (TSS-window guide design, not exon-tiling), check guide-to-TSS distance for top hits, and prefer per-guide-weighted models (e.g., MAGeCK MLE with guide efficiency terms, or tools designed for CRISPRi/a) rather than assuming knockout-screen guide behavior transfers directly.

## 10. Copy-number amplified regions create false-positive "essential gene" signal in knockout screens
- Trigger: Genome-wide Cas9 knockout screen in an aneuploid cancer cell line with known amplified regions.
- Mechanism: Cas9-induced double-strand breaks trigger a DNA-damage/cell-cycle-arrest response whose severity scales with the number of cut sites, i.e., with local copy number, independent of whether the targeted gene has any true fitness effect; genes in highly amplified loci can therefore appear strongly "essential" purely from cutting toxicity.
- Symptom: Clusters of "hit" genes that map to the same chromosomal amplicon, including neighboring genes with no plausible shared function; correlation between hit strength and known copy-number segments from orthogonal CNV data.
- Fix: Apply a copy-number correction step (CERES, CRISPRcleanR, or Chronos) whenever copy-number data for the cell line are available, or at minimum flag amplicon-clustered hits as likely CN artifacts and cross-check against public copy-number/dependency resources (e.g., DepMap) before including them as biological essential-gene calls.

## 11. Wrong or inconsistent normalization method changes which genes look significant
- Trigger: Switching between median-ratio, total-count, or control-sgRNA (non-targeting-guide) normalization without re-validating downstream thresholds.
- Mechanism: Each normalization assumes a different "most sgRNAs/genes are unchanged" baseline; in screens with a genuinely large fraction of the library depleted (very strong dropout / high-lethality treatment), total-count or median normalization can be biased by the shifting bulk of the distribution, inflating or deflating fold-changes for unaffected genes.
- Symptom: Hit list composition changes substantially when re-normalizing the same raw counts with a different method; nonessential control genes drift away from zero fold-change.
- Fix: Prefer control-sgRNA or nonessential-gene-anchored normalization when a nonessential reference set exists; otherwise sanity-check that nonessential control genes cluster near zero fold-change after normalization, and report the exact normalization method used.

## 12. High multiplicity of infection (MOI) confounds sgRNA-to-phenotype mapping
- Trigger: Transduction performed at MOI above the low range conventionally targeted for pooled screens (commonly ~0.3-0.5 for standard single-guide, single-integration designs).
- Mechanism: At higher MOI, a meaningful fraction of cells receive more than one integrated sgRNA, so an observed phenotype in a "sgRNA X" read cannot be unambiguously attributed to gene X; this dilutes true signal and can create spurious co-occurrence artifacts.
- Symptom: Weaker-than-expected separation between essential and nonessential control distributions given adequate library coverage and read depth; unexplained genotype-phenotype inconsistency across replicates.
- Fix: Confirm and report the transduction MOI used; if MOI is unknown or was not controlled, flag this as a screen-design limitation that reduces confidence in single-gene attribution, especially for weaker/borderline hits.

## 13. Off-target and multi-mapping sgRNAs are scored as if they were gene-specific
- Trigger: Library annotation includes guides that map to paralogous sequences, repetitive elements, or multiple genomic loci (common for older or non-curated libraries).
- Mechanism: A guide that cuts at several loci produces a fold-change reflecting the combined fitness effect of all targeted sites, which can be misattributed entirely to the annotated gene, especially problematic if one off-target site sits in an essential paralog or amplified region.
- Symptom: A "hit" gene supported by guides that, on inspection of the library annotation/BLAST-style mapping, have known secondary target sites; unusually strong effect size for a gene with no prior essentiality evidence.
- Fix: Use a current, curated library annotation with documented off-target/multi-mapping guides removed or flagged (e.g., BAGEL2's multi-target correction), and check top-hit guide sequences against a recent genome build before finalizing gene-level claims.

## 14. Pathway/GO enrichment run against the whole genome instead of the screened background
- Trigger: Enrichment analysis (GO, KEGG, Reactome, etc.) performed on the hit list using a default whole-genome/whole-transcriptome background.
- Mechanism: A pooled library only covers a defined gene set (e.g., genome-wide library, druggable genome, or a focused sub-library); using the whole genome as background inflates apparent enrichment for any pathway that happens to be over-represented in the library design itself, independent of the screen's biology.
- Symptom: Enrichment results dominated by pathways known to be library-design artifacts (e.g., over-represented gene families in a focused library) rather than screen-specific biology.
- Fix: Always set the enrichment background to the full set of genes/sgRNAs actually present in the screened library, consistent with `references/cns-method-guide.md`'s guidance on background universe.
