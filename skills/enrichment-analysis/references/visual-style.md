# Enrichment Visual Style

- Before final delivery, inspect PDF/PNG for overlapping labels; if labels collide, reduce label count, stagger annotations, wrap long labels, or move legends outside the panel.
- Export every figure as exactly two formats by default: PNG and PDF only.
- Use centered bold titles, axis title size around 21 pt, tick text around 16 pt, and keep x/y axis titles centered on their axes.
- Sort terms by `padj` ascending.
- Show at most 20 terms per figure unless requested.
- Wrap long term names.
- Bubble plot: x axis gene ratio, size count, color `-log10(padj)`.
- Bar plot: bar length count or gene ratio, color by category.
- Use distinct but restrained category colors.
- Always state whether enrichment was computed or only visualized from an existing table.
