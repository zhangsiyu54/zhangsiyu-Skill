# Paper Practice

Follow MAGeCK/CRISPRAnalyzeR/BAGEL/DepMap-style screen papers: screen QC and positive/negative selection should precede hit claims.

## Reference Sources

- MAGeCK GitHub: https://github.com/liulab-dfci/MAGeCK
- MAGeCK Genome Biology: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4
- CRISPRAnalyzeR GitHub: https://github.com/boutroslab/CRISPRAnalyzeR
- BAGEL2 GitHub: https://github.com/hart-lab/bagel
- DepMap Chronos: https://github.com/broadinstitute/chronos

## Paper-style Flow

1. Confirm sgRNA count matrix, treatment/control samples, library, contrast, and selection direction.
2. Review read count distribution, sgRNA representation, replicate consistency, and essential/nonessential control behavior when data allow.
3. Run MAGeCK/BAGEL/drugZ/Chronos-style scoring when dependencies and raw counts exist; otherwise visualize existing gene summary.
4. Report positive and negative selection hits, FDR thresholds, rank table, and pathway handoff.
5. Do not call validation without orthogonal evidence.

## Visualization Pattern

- Volcano and gene rank plots should show FDR threshold and sparse staggered labels.
- QC plots should show read depth, sgRNA distribution, and replicate correlation when available.
- Split positive and negative selection interpretation.
