# Paper Practice

Follow ArchR/Signac/chromVAR-style ATAC papers: QC metrics and genome build are mandatory before accessibility interpretation.

## Reference Sources

- ArchR GitHub: https://github.com/GreenleafLab/ArchR
- ArchR Nature Genetics: https://www.nature.com/articles/s41588-021-00790-6
- Signac GitHub: https://github.com/stuart-lab/signac
- Signac documentation: https://stuartlab.org/signac/
- chromVAR: https://greenleaflab.github.io/chromVAR/

## Paper-style Flow

1. Confirm genome build, peak coordinates, blacklist/filtering, sample groups, and fragment/peak count inputs.
2. Review FRiP, TSS enrichment, fragment size, peak count distribution, and sample/cluster QC.
3. For bulk ATAC, run differential accessibility and peak annotation.
4. For scATAC, run LSI/UMAP/clustering, gene activity, motif deviation, and peak-to-gene links when possible.
5. Report motif database, peak annotation source, genome build, and thresholds.

## Visualization Pattern

- Differential peak volcano: sparse labels and threshold lines.
- Motif enrichment: horizontal bars sorted by adjusted p value.
- Peak heatmap: blue-white-red centered at zero; keep peak labels short or omit them.
