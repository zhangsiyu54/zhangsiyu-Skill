# Biofigure Polish Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/figure_manifest.tsv`
- `output/figure_qc_report.tsv` (written by `scripts/polish_figures.py --manifest`)
- `output/report.md` (written by `scripts/polish_figures.py --manifest`; plain-language batch summary -- see below)
- `output/methods.md`
- `output/figures/` polished figure files (PNG raster preview + PDF vector, one pair per input)
- `output/style_notes.md`
- `output/run/environment.json`
- `output/run/methods.json`

Manifest columns (`output/figure_manifest.tsv`, written by `scripts/style_manifest.py`):
- source_file
- output_stem
- panel_label
- figure_type
- changes
- status

QC report columns (`output/figure_qc_report.tsv`, written by `scripts/polish_figures.py --manifest`):
- source_file
- output_stem
- panel_label
- width_px / height_px (`n/a_vector_source` for a PDF whose preview could not be rasterized)
- effective_dpi_at_89mm — heuristic effective print resolution at a standard single journal column width; used as a proxy for "font size too small at final print size" (references/pitfalls.md #4) since exact font metrics cannot be extracted from raster pixels without OCR
- low_resolution_risk — true when effective_dpi_at_89mm is below 300
- suspected_rainbow_or_jet_colormap — heuristic hue-band scan (references/pitfalls.md #2)
- suspected_red_green_diverging — heuristic hue-band scan (references/pitfalls.md #1)
- notes — free-text warnings, or `none`

`output/report.md` is the human-facing counterpart to `figure_qc_report.tsv`: a Markdown summary stating how many
figures were processed vs. flagged vs. failed, one section per flagged figure naming it explicitly with what was
flagged and the recommended fix (quoted from the heuristic warning text / `references/pitfalls.md`), a one-line
plain-language conclusion (e.g. "2/2 figures polished; 0 flagged for colormap issues; 1 figure has a low estimated
print DPI (~150) at column width -- consider requesting the original vector source rather than a rasterized PNG."),
and any processing errors (missing files, corrupt images, unsupported formats) that caused a figure to be skipped.
Always read this file (not just the raw TSV) before summarizing results to the user.

**Known limitation (state this in `output/style_notes.md` and `output/methods.md` whenever `scripts/polish_figures.py` is used on raster/vector inputs):** this script re-frames existing image files; it cannot fix overplotting, a rainbow/jet colormap, a red-green diverging scheme, a truncated axis, or undefined error bars that are already baked into the pixels or vector drawing commands of the source figure. `figure_qc_report.tsv` only flags a subset of these heuristically. When a flag appears (or overplotting is visually obvious), replotting the panel from its source data with `references/cns-method-guide.md`'s Design Decision Tree is the correct fix, not further polishing of the existing image.

Style notes:
- palette (reference `assets/palette.json` as the single source of truth; note the run of `scripts/polish_figures.py --check-palette`)
- typography
- label cleanup
- files changed
- unresolved issues (include any `figure_qc_report.tsv` flags not resolved by replotting, and the raster/vector limitation above when applicable)
