# Biofigure Polish Visual Style

- Before final delivery, inspect PDF/PNG for overlapping labels; if labels collide, reduce label count, stagger annotations, wrap long labels, or move legends outside the panel.
- Export every figure as exactly two formats by default: PNG (raster preview) and PDF (vector/print) only. `scripts/polish_figures.py` always writes both for every input; for a source PDF, the original vector file is preserved byte-for-byte and only a PNG preview is regenerated (see script docstring).
- Colors come from `assets/palette.json` — the single source of truth (references/pitfalls.md #3). Do not hardcode alternate hex values in scripts or docs; if you change a color, change it in `assets/palette.json` first.
  - Categorical/UI: blue `#3B82C4`, red `#D84A4A`, green `#3FA37A`, amber `#D99A2B`, purple `#7C65C1`.
  - Neutral text: `#1F2937`; secondary/muted text: `#667085`; grid/axis edge: `#E5E7EB`.
  - Diverging (signed data: log2FC, z-score, correlation): blue-red pair `diverging_low` `#3B4CC0` / `diverging_high` `#B40426`. Never use a red-green pair for direction (references/pitfalls.md #1).
  - Continuous/sequential (expression, abundance, density): `sequential_colormap` = matplotlib `viridis` (or `cividis`/`magma`). Never use rainbow/jet (references/pitfalls.md #2).
  - Categorical color-only encoding is capped at `max_categorical_hues` (10); beyond that, add a redundant channel (shape, direct labels, or facets) instead of more hues (references/pitfalls.md #12).
  - Run `python3 scripts/polish_figures.py --check-palette` to validate `assets/palette.json` against a simplified colorblind-safety heuristic before final delivery; it is a pre-check, not a substitute for a real CVD simulator (Coblis/Color Oracle).
- Use one font family across all panels.
- Use centered bold titles, axis title size around 21 pt, tick text around 16 pt, and keep x/y axis titles centered on their axes.
- Panel labels: bold uppercase A, B, C in top-left.
- Align plot widths, legend positions, and axis title spacing.
- Prefer direct labels when legends become too large.
- Keep backgrounds white for journal compatibility.
