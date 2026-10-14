# CNS Method Guide

Use this guide before cancer genomics mutation, CNV, tumor mutation burden, or mutational signature analysis.

## Minimal Inputs

- Require MAF-like mutation table or variant table, sample IDs, genome build, cohort, variant classification, gene symbols, and clinical/group metadata if stratified.
- For CNV, require segment-level or gene-level calls, genome build, calling method, and amplitude/status definition.
- For signatures, require trinucleotide mutation counts or exposure matrix plus signature reference version.

## QC and Harmonization

- Harmonize sample IDs across mutation, CNV, signature, and clinical tables.
- Filter artifacts, synonymous variants for oncoplots/TMB when appropriate, and low-confidence calls using explicit rules.
- Confirm tumor-only versus tumor-normal calling; report limitations for tumor-only data.
- Keep genome build and gene annotation consistent.

## Mutation Analysis

- Report mutation frequency, variant classes, top genes, co-occurrence/mutual exclusivity only with suitable sample size.
- TMB requires clear denominator, variant inclusion rule, and panel/exome size; avoid comparing TMB across platforms without normalization.
- Lollipop/domain plots require transcript/protein coordinate consistency.

## CNV and Signature Analysis

- CNV heatmaps should separate amplification, gain, loss, and deletion definitions.
- Do not interpret gene-level CNV as expression activation without expression evidence.
- Mutational signatures require enough mutations per sample; low-burden samples should be excluded or labeled unreliable.
- Report signature reference version such as COSMIC SBS version and fitting method.

## Required Figures

- Oncoplot, mutation frequency bar, variant classification summary, TMB boxplot, CNV heatmap, signature exposure stacked bar/heatmap, group comparison plots when metadata exist.

## Methods Text Must Include

- Variant source/caller, genome build, filtering rules, MAF columns used, TMB denominator, CNV threshold/caller, signature reference/fitting method, statistical tests, thresholds, and package versions.

## Algorithm and Design Decision Tree

Use this table to choose the right method/parameter for a given scenario rather than defaulting to one workflow for every input. See `references/pitfalls.md` for the failure modes each row is guarding against.

| Scenario | Recommended tool/method | Key parameter | Rationale |
|---|---|---|---|
| TMB from a targeted panel vs WES/WGS | Compute mutations/Mb using the assay's own callable/interrogated footprint (panel BED or exome capture manifest), not a generic "30 Mb exome" constant | Denominator = callable Mb at adequate depth for that specific assay; keep panel- and WES/WGS-derived TMB in separate reported groups | Panels enrich for hotspot/driver-dense regions, so raw mutations/Mb is not on the same scale as WES/WGS; mixing them without harmonization biases any TMB-based grouping or threshold comparison |
| Signature analysis with low mutation count (roughly <50-100 SNVs/sample) | Cohort-level de novo extraction or cohort-aggregated refitting rather than single-sample fitting; label low-count samples explicitly | Minimum mutation-count threshold (commonly >=50, ideally >=100 SNVs) before trusting a single-sample exposure; report bootstrap/cosine-similarity confidence | Below this range the 96-channel spectrum is dominated by sampling noise, and single-sample NMF/regression exposures become unstable and can assign activity to biologically implausible signatures |
| De novo vs refit signature extraction | De novo (e.g. SigProfilerExtractor-style NMF) when cohort is reasonably sized (commonly >=20-50 samples with adequate mutations) and novel processes are plausible; refitting (e.g. deconstructSigs/SigProfilerAssignment-style) against a COSMIC reference set when cohort/sample size is small or the question is "which known processes are present" | Rank/signature-number selection via stability (e.g. silhouette/reproducibility across NMF restarts) for de novo; fixed reference signature set + regression method for refitting | Refitting is bounded by the reference basis and cannot discover new biology, but is appropriate and more stable for small data; de novo can reveal cohort-specific novel signatures but needs enough samples/mutations for a stable factorization |
| Purity/ploidy-aware CNV interpretation | Purity/ploidy-aware caller (e.g. ASCAT, ABSOLUTE, FACETS, PureCN, Sequenza) rather than raw log2-ratio thresholding, especially when purity is roughly <30-40% or whole-genome doubling is suspected | Report estimated purity and ploidy per sample; use purity-adjusted integer copy number and per-sample-adjusted amplitude thresholds instead of one fixed log2 cutoff for the whole cohort | Raw copy ratio is diluted by normal-cell admixture; a fixed threshold systematically under-calls focal events (e.g. homozygous deletions) in low-purity samples and misclassifies neutral copy number in genome-doubled tumors |
| Driver gene / mutual exclusivity testing | Restrict pairwise testing to genes above a minimum frequency/count (e.g. >=3-5% or >=5-10 mutated samples); Fisher's exact test per pair with Benjamini-Hochberg FDR correction across all pairs tested; MutSigCV-style background-mutation-rate correction for driver significance rather than raw frequency ranking | FDR threshold (e.g. q<0.1-0.2 exploratory, q<0.05 confirmatory) and minimum cell-count per 2x2 table before reporting a pair | Low-frequency genes produce sparse, unstable contingency tables; testing many pairs without correction yields chance "hits"; raw mutation frequency ranking does not correct for gene length/background rate, so large non-driver genes (e.g. TTN) can outrank true drivers |
| Tumor-only vs tumor-normal variant calling | Tumor-normal calling when a matched normal is available; for tumor-only, apply population allele-frequency filtering (e.g. gnomAD) plus a panel-of-normals before treating calls as "somatic" | Population AF threshold (commonly <0.0001-0.001, adjusted for cohort ancestry) and PoN recurrence filter | Without a matched normal, germline heterozygous/homozygous variants are indistinguishable from somatic calls by allele count alone; population/PoN filtering removes most but not all germline contamination and must be disclosed as a limitation |
| Focal vs arm-level/broad CNV events | GISTIC2-style focal vs broad separation (or an explicit segment-length threshold) rather than one amplitude cutoff applied uniformly | Segment length threshold distinguishing focal (gene/small-region) from arm-scale calls; report focal peak q-value/amplitude separately from broad event calls | Collapsing focal and arm-level events into one "gain/loss" category conflates passenger arm-level aneuploidy with driver-level focal amplification/deletion targeting a specific gene |
| Clonal vs subclonal mutation burden / signature computation | Report VAF/CCF-based clonal-only subset alongside the all-mutations result when purity/multi-region data allow; state the VAF/depth inclusion filter explicitly | VAF or CCF threshold used to call "clonal" (cohort- and purity-dependent, no universal cutoff) | Low-VAF subclonal calls sit closer to the noise floor and are handled inconsistently across callers/filters; TMB and signature spectra computed on an inconsistent clonal/subclonal mix are not comparable across samples of different purity/depth |
| FFPE-derived or low-input library signature analysis | Screen for and flag known deamination/oxidative artifact signatures before biological interpretation; prefer a COSMIC reference set version that includes recognized artifact signatures | Specimen-type flag (FFPE vs fresh-frozen/blood) carried through to the signature step | FFPE fixation and oxidative damage during library prep generate C>T- and G>T-dominated artifacts that mimic real mutational signatures; treating them as biological without screening produces false signature/etiology claims |
| Cross-cohort or cross-build mutation/CNV comparison | Confirm genome build (hg19/hg38) for every input before joining; perform an explicit versioned liftover (e.g. UCSC liftOver/CrossMap) if builds differ, rather than assuming coordinate compatibility | Documented build per source table; liftover tool + version recorded in methods | Coordinate-based joins (mutation-to-CNV-segment, gene-to-cytoband) silently produce wrong matches for a subset of loci when builds differ, with no runtime error to catch it |
