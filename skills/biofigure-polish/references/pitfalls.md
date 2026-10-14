# Biofigure Polish Pitfalls

Review this checklist before finalizing any figure. Each entry: Trigger (when this applies) / Mechanism (why it goes wrong) / Symptom (what a reviewer or reader sees) / Fix (what to do instead).

## 1. Red-green diverging scheme for up/down-regulation

- **Trigger**: Volcano plots, MA plots, differential-expression heatmaps, or any panel encoding "up vs down" with red = increased and green = decreased.
- **Mechanism**: Red-green color-blindness (protanopia/deuteranopia) affects roughly 8% of men and ~0.5% of women; red and green channels collapse to similar brownish/yellowish hues, so direction of change becomes indistinguishable.
- **Symptom**: A meaningful fraction of readers (and possibly reviewers) cannot tell which genes are up- versus down-regulated from color alone; grayscale print reproduction shows almost no contrast between the two hues.
- **Fix**: Replace with a colorblind-safe diverging pair such as blue-red (e.g., `#3B4CC0`/`#B40426`, or ColorBrewer `RdBu`/`PuOr`), or an Okabe-Ito-derived pair (blue `#0072B2` vs vermillion `#D55E00`). Keep the scale symmetric around zero and add a redundant encoding (triangle/arrow direction, up/down text label) when the figure is central to the paper's claim.

## 2. Rainbow/jet colormap for continuous data

- **Trigger**: Heatmaps, spatial expression maps, or correlation matrices using the classic "jet"/rainbow gradient (blue-cyan-green-yellow-red).
- **Mechanism**: Jet is not perceptually uniform — equal steps in data value do not map to equal steps in perceived color difference (lightness oscillates non-monotonically), and it creates false banding/edges at yellow-green and cyan transitions.
- **Symptom**: Readers perceive artificial "boundaries" or hot spots that do not correspond to real thresholds in the data; the same colormap looks different in grayscale (some jet mid-tones become indistinguishable when converted to grayscale).
- **Fix**: Use a perceptually uniform sequential colormap (viridis, magma, cividis, or the scientific-colour-map families such as `batlow`) for magnitude-only data. Reserve diverging maps for signed data (see #1). Never use rainbow for anything encoding a continuous ordered variable.

## 3. Inconsistent color-to-category mapping across panels/figures

- **Trigger**: A multi-panel figure or full manuscript where "treatment" is red in Figure 1 but blue in Figure 3, or a cluster labeled cluster 3 changes color between a UMAP and its corresponding bar chart.
- **Mechanism**: Each panel or figure is generated independently (different scripts, different default palettes) without a shared color-to-label dictionary.
- **Symptom**: Readers must re-learn the color key for every panel; cross-referencing findings across figures becomes error-prone, and reviewers flag this as sloppy or (worse) as evidence the panels were not built from a coherent analysis.
- **Fix**: Define one explicit mapping (category -> hex code) as a single source of truth (a shared config, theme, or lookup table) and reuse it for every script/panel that plots that category, for the entire manuscript, including supplementary figures.

## 4. Font size too small at final print size

- **Trigger**: Multi-panel figures built at a large on-screen canvas size but that must be shrunk to fit a single (~89 mm) or double (~183 mm) journal column.
- **Mechanism**: Font sizes chosen for a 15-20 cm on-screen preview do not scale correctly once the figure is compressed to fit the column width required at submission; axis text can drop below the ~5-7 pt effective print size most journals require.
- **Symptom**: Axis labels, tick numbers, or legend text become illegible or blur at final print/PDF size, even though they looked fine on screen.
- **Fix**: Design and export at the true final physical dimensions (mm) and verify legibility by rendering/printing at 100% scale, not by eyeballing an on-screen zoom. Set explicit point sizes (not "auto-fit") and check the smallest text element against the target journal's minimum (commonly ~5-7 pt for Nature-family journals, but confirm against the specific journal's guide since minimums vary; hedge toward larger sizes when uncertain).

## 5. Raster PNG where vector is required

- **Trigger**: Line art, bar/box/line plots, schematics, or any figure with vector-native elements (text, straight lines, simple shapes) exported only as a low-DPI PNG.
- **Mechanism**: Many journals (e.g., Nature-family) require editable vector formats (PDF/EPS/SVG/AI) for line art so fonts remain sharp and editable at any zoom; a rasterized PNG bakes in a fixed pixel grid.
- **Symptom**: Text and lines look jagged or blurry when the reviewer/printer zooms in; production may reject the file at proofing stage; text becomes non-selectable/non-editable, blocking journal copyediting.
- **Fix**: Export line-art panels as vector PDF/EPS/SVG with fonts embedded or outlined only if required. Reserve raster (TIFF/PNG, typically >=300 dpi, often 600 dpi for combination art) for genuinely pixel-based content such as micrographs, gel scans, or very dense scatter/heatmap rasterizations.

## 6. Overplotting hides true density in scatter/volcano plots

- **Trigger**: Scatter or volcano plots with thousands to millions of points (e.g., single-cell UMAP, genome-wide SNP plots, large volcano plots) rendered as fully opaque markers.
- **Mechanism**: Points stack on top of each other; once a region is saturated, additional points add no visible information, so the eye cannot distinguish "1000 overlapping points" from "10 overlapping points."
- **Symptom**: Dense regions look like a single flat blob of solid color; apparent cluster size or effect-size clustering is misleading, and outliers in sparse regions get disproportionate visual weight.
- **Fix**: Use alpha blending (partial transparency, e.g., 0.1-0.3), point downsampling with a stated random seed, 2D density contours, or a hexbin/heatmap representation for the densest regions. State the chosen method in the methods text so density visualization choices are reproducible.

## 7. Truncated/non-zero-origin bar chart y-axis

- **Trigger**: Bar charts (and sometimes line charts) comparing group means or counts where the y-axis starts near the data minimum rather than at zero.
- **Mechanism**: Bar length is visually interpreted as proportional to magnitude; truncating the axis inflates the apparent ratio between bars even when the underlying difference is small.
- **Symptom**: A 5% real difference visually reads as a 3x or 5x difference; readers and reviewers accuse the figure of being misleading or "cherry-picked."
- **Fix**: Start bar-chart y-axes at zero. If the meaningful comparison genuinely requires a zoomed range (e.g., small differences on a large baseline), switch to a dot/point-range or line plot (which do not carry the "length = magnitude" convention) instead of truncating a bar chart, and mark any axis break explicitly and prominently if unavoidable.

## 8. Excessive precision in axis ticks or legends

- **Trigger**: Axis tick labels or legend entries showing values like `0.0000123456` or p-values like `p = 0.0483271`.
- **Mechanism**: Default plotting-library number formatting propagates full floating-point precision instead of a human-scaled rounding appropriate to the data's actual precision/uncertainty.
- **Symptom**: Cluttered, hard-to-scan axes; false impression of precision beyond what the assay/statistics actually support; inconsistent decimal places between adjacent tick labels.
- **Fix**: Round to a consistent, meaningful number of significant figures (commonly 2-3 sig figs for continuous axes; standard p-value conventions such as `p < 0.001` or 2-3 significant digits) and apply the same rounding rule to every panel in the figure set.

## 9. Missing, inconsistent, or misplaced panel labels

- **Trigger**: Multi-panel assembled figures where some panels have bold "A/B/C" labels in the top-left and others do not, or label font/size/case differs between panels.
- **Mechanism**: Panels are drawn by different scripts/tools (e.g., one from ggplot2, one from a manual PowerPoint pasted image) without a shared labeling pass.
- **Symptom**: Readers and reviewers cannot map in-text callouts ("as shown in Figure 2C") to the correct panel; the figure looks visibly stitched-together rather than a coherent assembled composite.
- **Fix**: Apply panel labels in a single final assembly pass (not per-source-script) using one consistent style (typically bold, sans-serif, same size, same corner position, same case) across the entire figure — this is exactly what `scripts/polish_figures.py` should standardize.

## 10. Legend placement obscuring data

- **Trigger**: Default plotting-library legend placement ("best" or "auto") on scatter plots, heatmaps, or UMAPs with data spread across the whole plotting area.
- **Mechanism**: Automatic legend positioning optimizes for "least apparent overlap at render time" but does not know which data points are scientifically important, and dense datasets can occupy every corner of the plot.
- **Symptom**: Legend box sits directly on top of key data points, a cluster label, or a significance bracket; reviewer cannot see the covered data.
- **Fix**: Move the legend fully outside the plotting area (right margin or below), or replace it with direct labels/callouts placed in empty space near the relevant series when the legend would otherwise have to overlap dense data.

## 11. Inconsistent significant figures or unit labeling across panels

- **Trigger**: One panel reports "Expression (TPM)" and another reports "log2(TPM+1)" without saying so, or one panel rounds to whole numbers while an adjacent panel of the same metric shows three decimal places.
- **Mechanism**: Panels are generated at different analysis stages (raw counts vs. normalized vs. log-transformed) and merged into one figure without reconciling axis titles/units.
- **Symptom**: Readers misinterpret magnitude or misunderstand what transformation was applied; apparent inconsistency undermines confidence in data handling even when the underlying analysis is correct.
- **Fix**: Explicitly state units and any transformation in every axis title (e.g., `log2(TPM + 1)`, not just "Expression"), and standardize decimal precision/units for the same measured quantity across every panel that reports it.

## 12. Single hue distinguishing more categories than the eye can resolve

- **Trigger**: Categorical palettes (cluster identity, cell type, sample group) with more than ~8-10 discrete colors, common in single-cell UMAPs or multi-sample cohort plots.
- **Mechanism**: Human color discrimination for small marks on a crowded plot degrades sharply beyond roughly 8-10 simultaneously distinguishable hues, regardless of color vision status; extended qualitative palettes often reuse similar hues at different lightness that photograph or print indistinguishably.
- **Symptom**: Two "different-colored" clusters look identical at normal viewing distance or after print/projector reproduction; legend entries become impossible to match back to marks.
- **Fix**: Cap direct color-only encoding at ~8-10 categories. Beyond that, combine color with a second channel (shape, direct on-plot text labels, small multiples/faceting, or grouping into a coarser category with color while retaining fine labels as text) rather than relying on hue alone.

## 13. Redundant or missing error-bar/uncertainty definition

- **Trigger**: Bar or point plots with error bars where the caption/legend does not state whether the whiskers represent SD, SEM, or a confidence interval.
- **Mechanism**: SD, SEM, and CI have very different widths and implications (SEM shrinks with sample size and is often mistaken for a measure of population spread), and journals differ in default convention.
- **Symptom**: Readers cannot correctly judge variability or statistical significance from the figure alone; apparent "tight" error bars may simply reflect large n rather than low biological variability.
- **Fix**: State explicitly in the axis/caption/legend what the error bars represent (SD, SEM, or XX% CI) and the sample size n per group; prefer showing individual data points (dot/strip plot or box+jitter) over bare bars with error bars whenever n is small enough to plot individually (see also `references/paper-practice.md` on data presentation).
