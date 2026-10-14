#!/usr/bin/env python3
import argparse
import csv
import json
import statistics
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except Exception as exc:
    raise SystemExit(
        "quickplot_singlecell.py requires matplotlib to export PNG/PDF. "
        "Install it with: python3 -m pip install matplotlib"
    ) from exc


# Categorical palette and text color match assets/palette.json ("clusters" / "text") so that
# Python quicklook figures and the R (run_seurat.R etc.) figures use the same visual identity.
PALETTE = ["#3B82C4", "#D84A4A", "#3FA37A", "#D99A2B", "#7C65C1", "#4DA6A6", "#F08AB4"]
TEXT = "#1F2937"
GRID = "#DCEAF5"
# Continuous scale for FeaturePlot / module score maps: light grey -> saturated blue,
# matching assets/palette.json's feature_low/feature_high and references/visual-style.md.
FEATURE_LOW = "#E5E7EB"
FEATURE_HIGH = "#2563A8"
FEATURE_CMAP = LinearSegmentedColormap.from_list("feature_scale", [FEATURE_LOW, FEATURE_HIGH])

# QC metrics recognized in an exported cell metadata table, with display labels.
QC_METRIC_CANDIDATES = [
    ("nCount", "nCount"),
    ("nCount_RNA", "nCount"),
    ("nFeature", "nFeature"),
    ("nFeature_RNA", "nFeature"),
    ("percent_mt", "percent.mt"),
    ("percent.mt", "percent.mt"),
    ("pct_counts_mt", "percent.mt"),
]


def safe_float(value):
    """Parse a CSV cell to float; return None (not an exception) for blank/non-numeric/NaN values.

    Exported tables commonly contain blank cells for missing values, and "NaN"/"NA" text for
    all-missing columns; treating either as a hard crash makes the whole pipeline fail on one bad
    cell. Callers filter out None before aggregating so an all-missing column degrades to "no data
    for this metric" instead of a traceback.
    """
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.upper() in ("NA", "NAN", "NULL", "NONE"):
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if parsed != parsed:  # NaN
        return None
    return parsed


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 16,
            "axes.titlesize": 24,
            "axes.labelsize": 21,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 15,
            "figure.titlesize": 22,
            "axes.edgecolor": "#B9D7EE",
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "text.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_plot(fig, path_base):
    path_base = Path(path_base)
    path_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_base.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(path_base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def read_csv(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def detect_qc_metrics(meta):
    """Return a list of (column_name, display_label) present in the metadata table."""
    if not meta:
        return []
    columns = meta[0].keys()
    found = []
    seen_labels = set()
    for col, label in QC_METRIC_CANDIDATES:
        if col in columns and label not in seen_labels:
            found.append((col, label))
            seen_labels.add(label)
    return found


def write_qc_violin(meta, out_base):
    """QC violin + jitter plot for nCount/nFeature/percent.mt style columns (required figure)."""
    metrics = detect_qc_metrics(meta)
    if not metrics:
        return False
    apply_style()
    usable = []
    for col, label in metrics:
        vals = [v for v in (safe_float(r.get(col)) for r in meta) if v is not None]
        if len(vals) >= 2:
            usable.append((col, label, vals))
    if not usable:
        return False
    metrics = [(col, label) for col, label, _ in usable]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.6 * len(metrics), 6.0))
    if len(metrics) == 1:
        axes = [axes]
    for ax, (col, label, vals) in zip(axes, usable):
        parts = ax.violinplot(vals, showmeans=False, showmedians=True, showextrema=True)
        for body in parts["bodies"]:
            body.set_facecolor(PALETTE[0])
            body.set_edgecolor(PALETTE[0])
            body.set_alpha(0.55)
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            if key in parts:
                parts[key].set_color(TEXT)
                parts[key].set_linewidth(1.1)
        jitter_x = [1 + (hash((label, i)) % 1000 / 1000 - 0.5) * 0.32 for i in range(len(vals))]
        ax.scatter(jitter_x, vals, s=8, color=PALETTE[1], alpha=0.35, linewidths=0, zorder=3)
        ax.set_title(label, pad=14, loc="center", fontweight="bold")
        ax.set_xticks([1])
        ax.set_xticklabels([label])
        ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    axes[0].set_ylabel("Value", labelpad=12)
    fig.suptitle("QC metrics", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_plot(fig, out_base)
    return True


def write_umap(meta, color_col, x_col, y_col, out_base, title=None):
    apply_style()
    labels = sorted({r[color_col] for r in meta})
    colors = {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(8.4, 6.8))
    ax.set_facecolor("#FBFDFF")
    n_dropped = 0
    for label in labels:
        subset = [r for r in meta if r[color_col] == label]
        pts = [(safe_float(r[x_col]), safe_float(r[y_col])) for r in subset]
        valid = [(x, y) for x, y in pts if x is not None and y is not None]
        n_dropped += len(pts) - len(valid)
        if not valid:
            continue
        ax.scatter(
            [x for x, _ in valid],
            [y for _, y in valid],
            s=18,
            c=colors[label],
            alpha=0.78,
            linewidths=0,
            label=label,
        )
    if n_dropped:
        print(f"note: dropped {n_dropped} cell(s) with non-numeric {x_col}/{y_col} from {out_base}")
    ax.set_title(title or f"UMAP by {color_col}", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("UMAP 1", labelpad=12)
    ax.set_ylabel("UMAP 2", labelpad=12)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    save_plot(fig, out_base)


def write_featureplot(meta, expr, gene, x_col, y_col, out_base):
    """Continuous single-gene expression UMAP (FeaturePlot-style), light grey -> saturated blue.

    Returns False without writing a file if no cell has a valid (numeric x, y, expression) triple
    (e.g. an all-blank/NaN expression column, or no metadata/expression cell IDs overlap) so the
    caller can skip listing this figure as generated instead of shipping a misleading empty plot.
    """
    apply_style()
    expr_by_cell = {r["cell"]: r for r in expr}
    cells = [r for r in meta if r["cell"] in expr_by_cell]
    triples = []
    for r in cells:
        x = safe_float(r[x_col])
        y = safe_float(r[y_col])
        v = safe_float(expr_by_cell[r["cell"]].get(gene))
        if x is not None and y is not None and v is not None:
            triples.append((x, y, v))
    if not triples:
        print(f"note: skipped {out_base} - no cells with numeric {x_col}/{y_col}/{gene}")
        return False
    xs, ys, vals = zip(*triples)
    fig, ax = plt.subplots(figsize=(7.8, 6.6))
    ax.set_facecolor("#FBFDFF")
    sc = ax.scatter(xs, ys, c=vals, cmap=FEATURE_CMAP, s=20, alpha=0.85, linewidths=0)
    ax.set_title(gene, pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("UMAP 1", labelpad=12)
    ax.set_ylabel("UMAP 2", labelpad=12)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Expression", fontsize=16, labelpad=10)
    fig.tight_layout()
    save_plot(fig, out_base)
    return True


def write_composition(meta, celltype_col, group_col, out_base):
    apply_style()
    groups = sorted({r[group_col] for r in meta})
    celltypes = sorted({r[celltype_col] for r in meta})
    counts = {(g, ct): 0 for g in groups for ct in celltypes}
    totals = {g: 0 for g in groups}
    for r in meta:
        counts[(r[group_col], r[celltype_col])] += 1
        totals[r[group_col]] += 1

    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    bottom = [0.0] * len(groups)
    for ci, ct in enumerate(celltypes):
        vals = [counts[(g, ct)] / max(1, totals[g]) * 100 for g in groups]
        ax.bar(groups, vals, bottom=bottom, color=PALETTE[ci % len(PALETTE)], width=0.62, label=ct)
        bottom = [bottom[i] + vals[i] for i in range(len(groups))]
    ax.set_title("Cell type proportion", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("Group", labelpad=12)
    ax.set_ylabel("Proportion (%)", labelpad=12)
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    save_plot(fig, out_base)


def write_dotplot(meta, expr, celltype_col, out_base):
    apply_style()
    expr_by_cell = {r["cell"]: r for r in expr}
    genes = [g for g in expr[0].keys() if g != "cell"]
    celltypes = sorted({r[celltype_col] for r in meta})
    xs, ys, sizes, colors = [], [], [], []
    for i, ct in enumerate(celltypes):
        cells = [r["cell"] for r in meta if r[celltype_col] == ct]
        for j, gene in enumerate(genes):
            vals = [v for v in (safe_float(expr_by_cell[c].get(gene)) for c in cells if c in expr_by_cell) if v is not None]
            avg = statistics.mean(vals) if vals else 0
            pct = sum(v > 0.2 for v in vals) / max(1, len(vals))
            xs.append(j)
            ys.append(i)
            sizes.append(70 + pct * 620)
            colors.append(avg)

    fig_w = max(12.4, 0.9 * len(genes) + 7.0)
    fig_h = max(5.8, 0.52 * len(celltypes) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="Blues", edgecolors="#4EA3F1", linewidths=1.0, alpha=0.9)
    ax.set_title("Marker dotplot", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("Marker genes", labelpad=12)
    ax.set_ylabel("Cell type", labelpad=12)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=35, ha="right")
    ax.set_yticks(range(len(celltypes)))
    ax.set_yticklabels(celltypes)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(-0.6, len(celltypes) - 0.4)
    ax.grid(True, color=GRID, linewidth=0.8)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.025)
    cbar.set_label("Average expression", fontsize=16, labelpad=10)
    for pct, size in [(0.25, 70 + 0.25 * 620), (0.50, 70 + 0.50 * 620), (0.75, 70 + 0.75 * 620)]:
        ax.scatter([], [], s=size, c="#D7ECFF", edgecolors="#4EA3F1", label=f"{int(pct * 100)}%")
    ax.legend(title="Expressed cells", title_fontsize=16, frameon=False, loc="center left", bbox_to_anchor=(1.34, 0.5))
    fig.tight_layout(rect=(0, 0, 0.74, 1))
    save_plot(fig, out_base)


def write_marker_heatmap(meta, expr, celltype_col, out_base):
    """Top marker heatmap: average expression of each marker gene across cell types/clusters."""
    apply_style()
    expr_by_cell = {r["cell"]: r for r in expr}
    genes = [g for g in expr[0].keys() if g != "cell"]
    celltypes = sorted({r[celltype_col] for r in meta})
    matrix = []
    for ct in celltypes:
        cells = [r["cell"] for r in meta if r[celltype_col] == ct]
        row = []
        for gene in genes:
            vals = [v for v in (safe_float(expr_by_cell[c].get(gene)) for c in cells if c in expr_by_cell) if v is not None]
            row.append(statistics.mean(vals) if vals else 0.0)
        matrix.append(row)

    fig_w = max(8.0, 0.9 * len(genes) + 3.5)
    fig_h = max(4.8, 0.62 * len(celltypes) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, cmap=FEATURE_CMAP, aspect="auto")
    ax.set_title("Top marker heatmap", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("Marker genes", labelpad=12)
    ax.set_ylabel("Cell type", labelpad=12)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=35, ha="right")
    ax.set_yticks(range(len(celltypes)))
    ax.set_yticklabels(celltypes)
    vmax = max((max(row) for row in matrix), default=1) or 1
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            text_color = "white" if val > 0.6 * vmax else TEXT
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=12, color=text_color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Average expression", fontsize=16, labelpad=10)
    fig.tight_layout()
    save_plot(fig, out_base)


def compute_module_scores(expr, genes):
    """Per-cell module score as the mean min-max scaled expression across the gene list.

    This is a lightweight, dependency-free approximation for quicklook use only; it is not
    a substitute for Seurat::AddModuleScore / Scanpy score_genes control-set-based scoring
    (see references/pitfalls.md #14 on module score comparability).
    """
    requested = [g for g in genes if g in expr[0]]
    present = []
    ranges = {}
    for g in requested:
        vals = [v for v in (safe_float(r.get(g)) for r in expr) if v is not None]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        ranges[g] = (lo, hi if hi > lo else lo + 1e-9)
        present.append(g)
    if not present:
        return {}, []
    scores = {}
    for r in expr:
        scaled = []
        for g in present:
            v = safe_float(r.get(g))
            if v is None:
                continue
            scaled.append((v - ranges[g][0]) / (ranges[g][1] - ranges[g][0]))
        if scaled:
            scores[r["cell"]] = statistics.mean(scaled)
    return scores, present


def write_module_score_umap(meta, scores, x_col, y_col, out_base, label):
    apply_style()
    triples = []
    for r in meta:
        if r["cell"] not in scores:
            continue
        x = safe_float(r[x_col])
        y = safe_float(r[y_col])
        if x is not None and y is not None:
            triples.append((x, y, scores[r["cell"]]))
    if not triples:
        print(f"note: skipped {out_base} - no cells with numeric {x_col}/{y_col} and a module score")
        return False
    xs, ys, vals = zip(*triples)
    fig, ax = plt.subplots(figsize=(7.8, 6.6))
    ax.set_facecolor("#FBFDFF")
    sc = ax.scatter(xs, ys, c=vals, cmap=FEATURE_CMAP, s=20, alpha=0.85, linewidths=0)
    ax.set_title(f"Module score: {label}", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("UMAP 1", labelpad=12)
    ax.set_ylabel("UMAP 2", labelpad=12)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Module score", fontsize=16, labelpad=10)
    fig.tight_layout()
    save_plot(fig, out_base)
    return True


def compute_top_markers(meta, expr, celltype_col, top_n=3):
    """Return {group: [(gene, mean_in_group, specificity), ...]}, sorted by specificity desc.

    Specificity is mean expression inside the group minus the mean across the other groups: a
    lightweight, dependency-free stand-in for a real marker/DE test (Wilcoxon/DESeq2) so the
    report can name candidate marker genes per cluster/cell type instead of only pointing at the
    heatmap/dotplot image. Not a substitute for `FindAllMarkers`/DESeq2 statistical testing.
    """
    expr_by_cell = {r["cell"]: r for r in expr}
    genes = [g for g in expr[0].keys() if g != "cell"]
    groups = sorted({r[celltype_col] for r in meta})
    group_cells = {g: [r["cell"] for r in meta if r[celltype_col] == g and r["cell"] in expr_by_cell] for g in groups}
    gene_means = {}
    for g, cells in group_cells.items():
        for gene in genes:
            vals = [v for v in (safe_float(expr_by_cell[c].get(gene)) for c in cells) if v is not None]
            gene_means[(g, gene)] = statistics.mean(vals) if vals else None
    result = {}
    for g in groups:
        rows = []
        for gene in genes:
            own = gene_means[(g, gene)]
            if own is None:
                continue
            others = [gene_means[(other, gene)] for other in groups if other != g and gene_means[(other, gene)] is not None]
            other_mean = statistics.mean(others) if others else 0.0
            rows.append((gene, own, own - other_mean))
        rows.sort(key=lambda x: x[2], reverse=True)
        result[g] = rows[:top_n]
    return result


def group_counts(meta, col):
    counts = {}
    for r in meta:
        counts[r[col]] = counts.get(r[col], 0) + 1
    return counts


def small_groups(counts, threshold=20):
    """Return {label: n} for groups under threshold cells (references/pitfalls.md #11)."""
    return {k: v for k, v in counts.items() if v < threshold}


def write_module_score_violin(meta, scores, celltype_col, out_base, label):
    apply_style()
    celltypes = sorted({r[celltype_col] for r in meta if r["cell"] in scores})
    data = [[scores[r["cell"]] for r in meta if r[celltype_col] == ct and r["cell"] in scores] for ct in celltypes]
    if not celltypes or not any(data):
        print(f"note: skipped {out_base} - no cells with a module score")
        return False
    fig, ax = plt.subplots(figsize=(max(6.4, 1.1 * len(celltypes) + 3.0), 6.0))
    parts = ax.violinplot(data, showmedians=True, showextrema=True)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(PALETTE[i % len(PALETTE)])
        body.set_edgecolor(PALETTE[i % len(PALETTE)])
        body.set_alpha(0.6)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        if key in parts:
            parts[key].set_color(TEXT)
    ax.set_xticks(range(1, len(celltypes) + 1))
    ax.set_xticklabels(celltypes, rotation=20, ha="right")
    ax.set_title(f"Module score by cell type: {label}", pad=18, loc="center", fontweight="bold")
    ax.set_ylabel("Module score", labelpad=12)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    fig.tight_layout()
    save_plot(fig, out_base)
    return True


def main():
    parser = argparse.ArgumentParser(description="Single-cell quick plots from exported tables with PNG/PDF export.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--group-col", default="condition")
    parser.add_argument("--cluster-col", default="cluster", help="Column with unannotated cluster IDs, if present")
    parser.add_argument("--sample-col", default="sample", help="Column with biological sample IDs, if present")
    parser.add_argument("--x-col", default="UMAP_1")
    parser.add_argument("--y-col", default="UMAP_2")
    parser.add_argument("--feature-genes", default="", help="Comma-separated genes for FeaturePlot; default is all genes in --expression")
    parser.add_argument("--module-genes", default="", help="Comma-separated gene list to compute a quicklook module score")
    parser.add_argument("--module-name", default="ModuleScore", help="Label for the module score plots/tables")
    parser.add_argument("--out", default="output")
    args = parser.parse_args()

    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    meta = read_csv(args.metadata)
    expr = read_csv(args.expression)
    if not meta:
        raise SystemExit(
            f"No rows found in metadata table: {args.metadata}. "
            "Check the file has a header row plus at least one data row."
        )
    if not expr:
        raise SystemExit(
            f"No rows found in expression table: {args.expression}. "
            "Check the file has a header row plus at least one data row."
        )
    if "cell" not in meta[0]:
        raise SystemExit(f"Metadata table {args.metadata} has no 'cell' column: found {list(meta[0].keys())}")
    if "cell" not in expr[0]:
        raise SystemExit(f"Expression table {args.expression} has no 'cell' column: found {list(expr[0].keys())}")
    for col_name, col_val in (("--celltype-col", args.celltype_col), ("--x-col", args.x_col), ("--y-col", args.y_col)):
        if col_val not in meta[0]:
            raise SystemExit(
                f"Column '{col_val}' (from {col_name}) not found in metadata table. "
                f"Available columns: {list(meta[0].keys())}"
            )
    meta_ids = {r["cell"] for r in meta}
    expr_ids = {r["cell"] for r in expr}
    overlap = meta_ids & expr_ids
    if not overlap:
        raise SystemExit(
            "No overlapping cell IDs between metadata and expression tables — cannot join them. "
            f"Example metadata cell IDs: {sorted(meta_ids)[:3]}; example expression cell IDs: {sorted(expr_ids)[:3]}. "
            "Check for barcode suffix/prefix mismatches (e.g. '-1', sample prefixes) or wrong input files."
        )
    if len(overlap) < len(meta_ids):
        print(
            f"note: {len(meta_ids) - len(overlap)} of {len(meta_ids)} metadata cells have no matching row "
            "in the expression table and will be skipped in expression-based plots (dotplot/heatmap/featureplot)."
        )

    # write_composition/celltype_proportion.tsv require args.group_col to exist; when the
    # metadata has no such column (e.g. a single-condition dataset), fall back to a constant
    # "all" group instead of crashing with a KeyError, mirroring run_seurat.R's fallback.
    group_col_present = args.group_col in meta[0]
    if not group_col_present:
        for r in meta:
            r[args.group_col] = "all"
        print(f"note: group column '{args.group_col}' not found in metadata; treating all cells as one group ('all').")

    written_figures = []

    # Required: UMAP by cell type, and by cluster/sample/condition when those columns exist.
    write_umap(meta, args.celltype_col, args.x_col, args.y_col, out / "figures" / "umap_celltype")
    written_figures.append("figures/umap_celltype.png")
    if args.cluster_col in meta[0]:
        write_umap(
            meta, args.cluster_col, args.x_col, args.y_col, out / "figures" / "umap_clusters", title="UMAP by cluster"
        )
        written_figures.append("figures/umap_clusters.png")
    if args.sample_col in meta[0]:
        write_umap(
            meta, args.sample_col, args.x_col, args.y_col, out / "figures" / "umap_sample", title="UMAP by sample"
        )
        written_figures.append("figures/umap_sample.png")
    if group_col_present:
        write_umap(
            meta, args.group_col, args.x_col, args.y_col, out / "figures" / "umap_condition", title="UMAP by condition"
        )
        written_figures.append("figures/umap_condition.png")

    # Required: QC violin/histogram from nCount/nFeature/percent.mt style columns.
    if write_qc_violin(meta, out / "figures" / "qc_violin"):
        written_figures.append("figures/qc_violin.png")

    write_composition(meta, args.celltype_col, args.group_col, out / "figures" / "celltype_proportion")
    written_figures.append("figures/celltype_proportion.png")
    write_dotplot(meta, expr, args.celltype_col, out / "figures" / "dotplot_markers")
    written_figures.append("figures/dotplot_markers.png")

    # Required: top marker heatmap.
    write_marker_heatmap(meta, expr, args.celltype_col, out / "figures" / "heatmap_markers")
    written_figures.append("figures/heatmap_markers.png")

    # Required when target genes exist: FeaturePlot (continuous expression on UMAP).
    all_genes = [g for g in expr[0].keys() if g != "cell"]
    feature_genes = [g.strip() for g in args.feature_genes.split(",") if g.strip()] or all_genes
    feature_genes = [g for g in feature_genes if g in all_genes]
    first_feature_written = False
    for gene in feature_genes:
        if write_featureplot(meta, expr, gene, args.x_col, args.y_col, out / "figures" / f"featureplot_{gene}"):
            written_figures.append(f"figures/featureplot_{gene}.png")
            if not first_feature_written:
                # Also write the canonical singular filename expected by references/deliverables.md.
                write_featureplot(meta, expr, gene, args.x_col, args.y_col, out / "figures" / "featureplot")
                written_figures.append("figures/featureplot.png")
                first_feature_written = True

    # Optional: module score plot when a gene list is provided.
    module_genes = [g.strip() for g in args.module_genes.split(",") if g.strip()]
    module_summary_path = None
    if module_genes:
        scores, used_genes = compute_module_scores(expr, module_genes)
        if not scores:
            print(
                f"note: module score '{args.module_name}' skipped - none of the requested genes "
                f"({', '.join(module_genes)}) had usable numeric expression values"
            )
        if scores:
            if write_module_score_umap(
                meta, scores, args.x_col, args.y_col, out / "figures" / "module_score_umap", args.module_name
            ):
                written_figures.append("figures/module_score_umap.png")
            if write_module_score_violin(
                meta, scores, args.celltype_col, out / "figures" / "module_score_violin", args.module_name
            ):
                written_figures.append("figures/module_score_violin.png")
            with (out / "results" / "module_score_percell.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["cell", args.celltype_col, "module_score"])
                for r in meta:
                    if r["cell"] in scores:
                        writer.writerow([r["cell"], r[args.celltype_col], f"{scores[r['cell']]:.4f}"])
            by_ct = {}
            for r in meta:
                if r["cell"] in scores:
                    by_ct.setdefault(r[args.celltype_col], []).append(scores[r["cell"]])
            with (out / "results" / "module_score_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow([args.celltype_col, "n_cells", "mean_score"])
                for ct, vals in sorted(by_ct.items()):
                    writer.writerow([ct, len(vals), f"{statistics.mean(vals):.4f}"])
            module_summary_path = "results/module_score_summary.tsv"

    counts = {}
    for row in meta:
        key = (row[args.group_col], row[args.celltype_col])
        counts[key] = counts.get(key, 0) + 1
    with (out / "results" / "celltype_proportion.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "celltype", "count"])
        for (group, celltype), count in sorted(counts.items()):
            writer.writerow([group, celltype, count])

    qc_metrics = detect_qc_metrics(meta)
    qc_summary = {}
    for col, label in qc_metrics:
        vals = [v for v in (safe_float(r.get(col)) for r in meta) if v is not None]
        if vals:
            qc_summary[label] = {
                "median": round(statistics.median(vals), 3),
                "mean": round(statistics.mean(vals), 3),
                "min": round(min(vals), 3),
                "max": round(max(vals), 3),
            }

    # Generic QC pass/fail estimate using the same default thresholds as run_seurat.R
    # (nFeature >= 200, percent.mt <= 15). These are defaults, not dataset-tuned cutoffs
    # (references/pitfalls.md #5) -- reported transparently as an estimate, not a filtering action.
    qc_col_by_label = {label: col for col, label in qc_metrics}
    nfeat_col = qc_col_by_label.get("nFeature")
    mt_col = qc_col_by_label.get("percent.mt")
    qc_pass = None
    if nfeat_col or mt_col:
        n_pass = 0
        n_checked = 0
        for r in meta:
            nfeat = safe_float(r.get(nfeat_col)) if nfeat_col else None
            mt = safe_float(r.get(mt_col)) if mt_col else None
            if nfeat is None and mt is None:
                continue
            n_checked += 1
            ok = (nfeat is None or nfeat >= 200) and (mt is None or mt <= 15)
            n_pass += 1 if ok else 0
        if n_checked:
            criteria = ", ".join(
                x
                for x in [f"{nfeat_col} >= 200" if nfeat_col else None, f"{mt_col} <= 15" if mt_col else None]
                if x
            )
            qc_pass = {"checked": n_checked, "pass": n_pass, "fail": n_checked - n_pass, "criteria": criteria}

    # Cell type composition breakdown (overall, across all groups).
    celltype_counts = group_counts(meta, args.celltype_col)
    total_cells = len(meta)
    top_celltypes = sorted(celltype_counts.items(), key=lambda kv: kv[1], reverse=True)
    composition_report_lines = [f"  - {ct}: {n} cells ({100 * n / total_cells:.1f}%)" for ct, n in top_celltypes[:6]]

    # Low-n group caveats (references/pitfalls.md #11).
    small_celltype_groups = small_groups(celltype_counts)
    small_cluster_groups = small_groups(group_counts(meta, args.cluster_col)) if args.cluster_col in meta[0] else {}

    # Marker findings: top genes per cell type by a simple mean-expression specificity score.
    top_markers_by_group = compute_top_markers(meta, expr, args.celltype_col, top_n=3)
    marker_lines = [
        f"  - {ct}: " + ", ".join(f"{gene} (avg={mean:.2f})" for gene, mean, _spec in rows)
        for ct, rows in top_markers_by_group.items()
        if rows
    ]

    (out / "summary.tsv").write_text(
        "item\tvalue\n"
        f"cells\t{total_cells}\n"
        f"celltypes\t{len(celltype_counts)}\n"
        + (
            f"qc_cells_checked\t{qc_pass['checked']}\nqc_cells_pass\t{qc_pass['pass']}\nqc_cells_fail\t{qc_pass['fail']}\n"
            if qc_pass
            else ""
        ),
        encoding="utf-8",
    )
    (out / "run" / "qc.json").write_text(
        json.dumps(
            {"cells": total_cells, "quicklook": True, "qc_metrics": qc_summary, "qc_pass_fail_estimate": qc_pass},
            indent=2,
        ),
        encoding="utf-8",
    )
    (out / "run" / "artifacts.json").write_text(
        json.dumps({"figures": sorted(set(written_figures)), "module_score_table": module_summary_path}, indent=2),
        encoding="utf-8",
    )

    # ---- Caveats tied to references/pitfalls.md; only emitted when they actually apply ----
    caveats = []
    if small_celltype_groups:
        detail = ", ".join(f"{k} (n={v})" for k, v in sorted(small_celltype_groups.items(), key=lambda kv: kv[1]))
        caveats.append(
            f"Low cell count in `{args.celltype_col}`: {detail}. Annotation/statistics for these groups "
            "should be treated as unstable, not confirmatory (references/pitfalls.md #11)."
        )
    if small_cluster_groups:
        detail = ", ".join(f"{k} (n={v})" for k, v in sorted(small_cluster_groups.items(), key=lambda kv: kv[1]))
        caveats.append(
            f"Low cell count in `{args.cluster_col}`: {detail}. Treat these clusters as exploratory "
            "(references/pitfalls.md #11)."
        )
    if module_summary_path and len(module_genes) == 1:
        caveats.append(
            f"Module score '{args.module_name}' is based on a single gene ('{module_genes[0]}'); a single "
            "marker cannot confirm cell identity or pathway activity on its own — treat as exploratory "
            "(references/pitfalls.md #6 and #14)."
        )
    if group_col_present:
        sample_col_present = args.sample_col in meta[0]
        replicate_caveat_needed = True
        if sample_col_present:
            group_sample_counts = {}
            for r in meta:
                group_sample_counts.setdefault(r[args.group_col], set()).add(r[args.sample_col])
            min_samples = min((len(s) for s in group_sample_counts.values()), default=0)
            replicate_caveat_needed = min_samples < 2
        if replicate_caveat_needed:
            reason = "fewer than 2 samples in at least one group" if sample_col_present else "no biological sample column was found"
            caveats.append(
                f"Cell type proportions are pooled across cells per `{args.group_col}`, and {reason} — a "
                "statistically valid composition comparison needs per-sample proportions and a sample-aware "
                "test (references/pitfalls.md #12), not a pooled cell-level count."
            )
    if len(overlap) < len(meta_ids):
        caveats.append(
            f"{len(meta_ids) - len(overlap)} of {len(meta_ids)} metadata cells had no matching row in the "
            "expression table and were excluded from marker/dotplot/heatmap/featureplot statistics."
        )

    # ---- Plain-language conclusion synthesizing composition + markers + caveats ----
    dominant = ", ".join(f"{ct} ({100 * n / total_cells:.0f}%)" for ct, n in top_celltypes[:3])
    conclusion = (
        f"This quicklook covers {total_cells} cells across {len(celltype_counts)} group(s) in "
        f"`{args.celltype_col}`, dominated by {dominant}."
    )
    first_marked = next(((ct, rows) for ct, rows in top_markers_by_group.items() if rows), None)
    if first_marked:
        first_ct, first_rows = first_marked
        conclusion += (
            f" The clearest marker signal in this panel is {first_rows[0][0]} in {first_ct}; "
            "see the per-group marker list below for the rest."
        )
    if qc_pass:
        conclusion += (
            f" Of {qc_pass['checked']} cells with recorded QC metrics, {qc_pass['pass']} "
            f"({100 * qc_pass['pass'] / qc_pass['checked']:.1f}%) pass generic thresholds ({qc_pass['criteria']}); "
            "confirm against this dataset's own QC distributions before filtering (references/pitfalls.md #5)."
        )
    if caveats:
        conclusion += " See caveats below before treating these results as confirmatory."
    else:
        conclusion += (
            " No automatic caveats were triggered for this input, but this quicklook route remains "
            "exploratory relative to a full Seurat/Scanpy pipeline."
        )

    report_lines = [
        "# Single-cell quicklook report",
        "",
        "Generated from exported tables (quicklook route, not a full Seurat/Scanpy pipeline).",
        "",
        "## Results",
        "",
        f"- Cells: {total_cells}",
        f"- Cell types/groups (`{args.celltype_col}`): {len(celltype_counts)}",
        f"- Figures written: {len(set(written_figures))}",
    ]
    if qc_summary:
        report_lines.append("- QC metrics summarized: " + ", ".join(qc_summary.keys()))
    if qc_pass:
        report_lines.append(
            f"- QC pass/fail estimate ({qc_pass['criteria']}): {qc_pass['pass']} pass / {qc_pass['fail']} fail "
            f"out of {qc_pass['checked']} cells with recorded QC metrics"
        )
    if module_summary_path:
        report_lines.append(
            f"- Module score ({args.module_name}) computed for {len(used_genes)}/{len(module_genes)} requested genes"
        )

    report_lines += ["", "### Cell type composition (top groups)", ""] + (
        composition_report_lines or ["  - (no cell type data)"]
    )

    if marker_lines:
        report_lines += [
            "",
            "### Top marker genes per group (simple specificity score, not a statistical test)",
            "",
        ] + marker_lines

    if caveats:
        report_lines += ["", "### Caveats", ""] + [f"- {c}" for c in caveats]

    report_lines += [
        "",
        "### Conclusion",
        "",
        conclusion,
        "",
        "Module scores here are a lightweight min-max-scaled mean, not Seurat::AddModuleScore's "
        "control-gene-bin method; treat as exploratory only (see references/pitfalls.md #14). Marker genes "
        "above are ranked by a simple mean-expression specificity score computed only from the exported "
        "marker panel, not a full differential expression test.",
    ]
    (out / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"wrote single-cell quicklook outputs to {out}")


if __name__ == "__main__":
    main()
