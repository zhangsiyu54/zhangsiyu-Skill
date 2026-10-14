# CNS Method Guide

Use this guide before polishing bioinformatics figures for manuscript-style delivery.

## Minimal Inputs

- Require source figures or result tables, intended journal/page format, panel labels, group labels, color meaning, significance annotations, and final output format.
- If only raster images are provided, preserve scientific content and document that data-level replotting was not possible.
- If source tables are provided, prefer replotting over pixel-level editing.

## Scientific Integrity

- Never alter data values, group labels, thresholds, P values, axes, or significance markers to improve appearance.
- Do not crop out outliers, failed samples, or unfavorable groups unless the user explicitly requests a scientifically justified exclusion and it is documented.
- Keep figure legends and methods consistent with plotted statistics.

## Visual Standards

- Use consistent typography, panel label style, line widths, point sizes, axis formatting, color palette, and legend placement across the full figure set.
- Use colorblind-safe palettes and avoid red/green-only contrasts.
- Make raster exports high resolution and keep vector/PDF outputs when possible.
- Align panels to a grid with equal margins; keep legends outside data-dense panels unless space requires otherwise.
- Remove decorative effects that do not encode data.

## Design Decision Tree

Use this table to pick a palette/approach before touching any script. When a figure matches more than one row (e.g., a diverging heatmap that is also a dense scatter), apply every matching row's guidance rather than picking one.

| Scenario | Recommended palette/approach | Rationale |
| --- | --- | --- |
| Diverging data (up/down regulation, log2FC, correlation coefficient, z-score) | Blue-red diverging scale (e.g., `#3B4CC0`-`#B40426`, ColorBrewer `RdBu`) or an Okabe-Ito-derived blue/vermillion pair; center the scale at zero with symmetric limits | Red-green is unreadable to ~8% of men with red-green color blindness; blue-red diverging pairs remain distinguishable under protanopia/deuteranopia and in grayscale print |
| Continuous sequential data (expression level, abundance, density, intensity) | Perceptually uniform sequential colormap: viridis, magma, or cividis (cividis is itself CVD-optimized) | Equal steps in the colormap correspond to equal steps in perceived brightness/color, so gradients are not visually distorted; avoids the false banding of rainbow/jet |
| Categorical group data (cell type, cluster identity, treatment arm) with <=8-10 groups | Okabe-Ito 8-color qualitative palette, or ColorBrewer `Set2`/`Dark2` | These qualitative sets are explicitly validated for simultaneous discriminability under common color vision deficiencies and hold up in grayscale/print |
| Colorblind-accessible requirement stated explicitly (by user or target journal) | Okabe-Ito for categorical, viridis/cividis for sequential, blue-red for diverging; verify with a CVD simulator (e.g., Coblis, Color Oracle) before final export | Matches the palette families most journals (e.g., Nature-family) now recommend or require, and gives an objective pass/fail check rather than a subjective one |
| Multi-panel manuscript figure consistency (same category must look identical everywhere) | Define one category-to-hex lookup table (shared config/theme) and reuse it in every script that renders that category, across main and supplementary figures | Prevents the same "treatment" group from being red in Figure 1 and blue in Figure 3, which reviewers flag as sloppy or inconsistent analysis |
| Vector vs. raster export choice | Vector PDF/EPS/SVG for line art, bar/box/line plots, and schematics; raster TIFF/PNG (>=300 dpi, often 600 dpi) only for micrographs, gel scans, or extremely dense point/heatmap rasterizations | Journals require editable vector text/lines for crisp reproduction at any zoom; raster is appropriate only where the content is inherently pixel-based |
| Dense scatter/volcano overplotting (thousands+ of points) | Alpha blending (~0.1-0.3 opacity), 2D density contours, or hexbin/heatmap representation instead of fully opaque markers | Fully opaque overlapping points saturate visually and hide true density differences; density-aware encodings restore that information |
| Heatmap combining magnitude and direction (e.g., row-scaled z-score expression across conditions) | Diverging colormap centered at zero with symmetric color limits, plus row/column dendrograms or annotation bars for group structure | An asymmetric or sequential-only scale on signed data visually exaggerates one direction over the other |
| Large qualitative category count (>8-10 clusters, e.g., scRNA-seq UMAP with dozens of clusters) | Combine color with a redundant channel: direct on-plot text/number labels, shape, or faceting into small multiples, rather than adding more hues | Human simultaneous color discrimination degrades sharply beyond roughly 8-10 hues regardless of color vision status; more hues just add ambiguity |
| Survival/time-to-event or other line-based group comparison for print reproduction | Consistent categorical palette plus redundant linetype (solid/dashed/dotted) or marker per group | Ensures the figure remains interpretable in grayscale print or when colors shift slightly in production |

See `references/pitfalls.md` for the failure modes each of these choices is designed to avoid.

## Figure-Specific Rules

- Volcano plots: show thresholds, label only selected genes, avoid label overlap, and keep effect direction clear.
- Heatmaps: include row/column annotations, readable color scale, and avoid overcrowding with too many features.
- UMAP/spatial plots: use small points, stable colors, and legends that do not cover data.
- Survival plots: include risk table, censoring, confidence interval when appropriate, and exact group counts.
- Enrichment plots: encode adjusted P value, gene ratio/count, and term names without truncating meaning.

## Deliverables

- Provide polished PNG and PDF, source or style manifest, figure inventory, and a short change log of visual-only edits.
- For journal-like panels, export assembled figure plus individual panels.

## Methods Text Must Include

- Whether figures were replotted or polished from raster, source files used, software/packages, palette, export resolution, and any visual-only transformations.
