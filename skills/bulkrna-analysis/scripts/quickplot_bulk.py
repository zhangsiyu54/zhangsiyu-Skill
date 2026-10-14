#!/usr/bin/env python3
import argparse
import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except Exception as exc:
    raise SystemExit(
        "quickplot_bulk.py requires matplotlib and numpy to export PNG/PDF and compute PCA/correlation. "
        "Install them with: python3 -m pip install matplotlib numpy"
    ) from exc


UP = "#E45B5B"
DOWN = "#4EA3F1"
NS = "#B8C2CC"
TEXT = "#173B5F"
GRID = "#DCEAF5"

OUTLIER_CORR_GAP = 0.1  # heuristic: flag a sample if its mean same-group correlation
# trails the group average by more than this many Pearson-r units.
SMALL_N_THRESHOLD = 3  # pitfalls.md #10
LARGE_N_THRESHOLD = 30  # pitfalls.md #12
LOG_LIKE_MAX_HEURISTIC = 50.0


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


def read_table(path):
    """Read a CSV/TSV into (rows, fieldnames), with a clear error on unreadable input."""
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError as exc:
        raise SystemExit(f"Input file '{path}' does not exist. Check the path passed to --counts/--metadata.") from exc
    if not text.strip():
        raise SystemExit(f"Input file '{path}' is empty (no header or rows).")
    try:
        dialect = csv.Sniffer().sniff(text[:2048], delimiters=",\t")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    rows = list(reader)
    fieldnames = reader.fieldnames or []
    return rows, fieldnames


def bh_adjust(pairs):
    ordered = sorted(pairs, key=lambda x: x[1])
    m = len(ordered)
    adjusted = {}
    prev = 1.0
    for rank, (gene, p) in enumerate(reversed(ordered), start=1):
        idx = m - rank + 1
        val = min(prev, p * m / idx)
        adjusted[gene] = min(1.0, val)
        prev = val
    return adjusted


def normal_p_from_effect(effect, spread, n):
    if spread <= 0:
        return 1e-6 if abs(effect) > 0 else 1.0
    z = abs(effect) / (spread / math.sqrt(max(1, n)))
    return max(1e-12, math.erfc(z / math.sqrt(2)))


def infer_expression_scale(raw_matrix, has_non_integer, requested_scale):
    """Return how values should be transformed for quicklook statistics/plots.

    raw-count: integer count-like data; use log2(x+1).
    log-expression: already log/normalized continuous expression; use values as-is.
    linear-expression: non-count abundance/intensity; use log2(x+1).

    The auto rule is intentionally conservative for the recurring failure mode:
    small non-integer expression values (TPM/FPKM already log transformed,
    voom/logCPM, microarray/log intensity) must not be logged again.
    """
    if requested_scale != "auto":
        return requested_scale
    values = [v for vals in raw_matrix.values() for v in vals.values()]
    max_value = max(values) if values else 0.0
    if has_non_integer and max_value <= LOG_LIKE_MAX_HEURISTIC:
        return "log-expression"
    return "raw-count" if not has_non_integer else "linear-expression"


def transform_value(value, expression_scale):
    if expression_scale == "log-expression":
        return value
    return math.log2(value + 1)


def transform_label(expression_scale):
    if expression_scale == "log-expression":
        return "input expression value"
    if expression_scale == "linear-expression":
        return "log2(input value + 1)"
    return "log2(count + 1)"


def write_volcano(rows, out_base, padj_threshold=0.05, log2fc_threshold=1.0):
    apply_style()
    xs = [r["log2FC"] for r in rows]
    ys = [-math.log10(max(r["padj"], 1e-12)) for r in rows]
    xlim = max(log2fc_threshold * 1.6, max(abs(x) for x in xs) * 1.16)
    ylim = max(-math.log10(padj_threshold) * 1.15, max(ys) * 1.18, 2.2)
    groups = {
        "Up": ([r["log2FC"] for r in rows if r["change"] == "up"], [-math.log10(max(r["padj"], 1e-12)) for r in rows if r["change"] == "up"], UP),
        "Down": ([r["log2FC"] for r in rows if r["change"] == "down"], [-math.log10(max(r["padj"], 1e-12)) for r in rows if r["change"] == "down"], DOWN),
        "Not significant": ([r["log2FC"] for r in rows if r["change"] == "ns"], [-math.log10(max(r["padj"], 1e-12)) for r in rows if r["change"] == "ns"], NS),
    }

    fig, ax = plt.subplots(figsize=(9.8, 6.8))
    ax.set_facecolor("#FBFDFF")
    for label, (gx, gy, color) in groups.items():
        if gx:
            ax.scatter(gx, gy, s=34 if label != "Not significant" else 22, c=color, alpha=0.82, linewidths=0, label=label)
    ax.axvline(-log2fc_threshold, color="#A7C6DD", linestyle="--", linewidth=1.2)
    ax.axvline(log2fc_threshold, color="#A7C6DD", linestyle="--", linewidth=1.2)
    ax.axhline(-math.log10(padj_threshold), color="#A7C6DD", linestyle="--", linewidth=1.2)
    label_rows = [row for row in rows if row["padj"] < padj_threshold and abs(row["log2FC"]) >= log2fc_threshold][:4]
    side_counts = {"left": 0, "right": 0}
    for row in label_rows:
        y = -math.log10(max(row["padj"], 1e-12))
        side = "left" if row["log2FC"] < 0 else "right"
        rank_on_side = side_counts[side]
        side_counts[side] += 1
        x_offset = -24 if side == "left" else 24
        y_offset = 16 + rank_on_side * 28
        ax.annotate(
            row["gene"],
            xy=(row["log2FC"], y),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            ha="right" if row["log2FC"] < 0 else "left",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
            arrowprops={"arrowstyle": "-", "color": "#A7C6DD", "linewidth": 0.8},
        )
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(0, ylim)
    ax.set_title(
        f"Volcano plot\npadj < {padj_threshold:g}, |log2FC| >= {log2fc_threshold:g}",
        pad=18,
        loc="center",
        fontweight="bold",
    )
    ax.set_xlabel("log2 fold change", labelpad=12)
    ax.set_ylabel("-log10 adjusted p value", labelpad=12)
    ax.grid(True, color=GRID, linewidth=0.9)
    # Labeled points can sit at the xlim edge with a +/-24pt outward text offset
    # (see label loop above); push the legend further out than the default 1.02
    # anchor so it never overlaps the rightmost gene-name annotation.
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.28, 1.0), borderaxespad=0)
    fig.tight_layout()
    save_plot(fig, out_base)


def write_heatmap(matrix, samples, rows, out_base, top_n=30, expression_scale="raw-count"):
    apply_style()
    top_rows = rows[: min(top_n, len(rows))]
    values = []
    genes = []
    for row in top_rows:
        vals = [transform_value(matrix[row["gene"]][s], expression_scale) for s in samples]
        mean = statistics.mean(vals)
        sd = statistics.pstdev(vals) or 1
        values.append([(v - mean) / sd for v in vals])
        genes.append(row["gene"])
    fig_w = max(7.2, 1.02 * len(samples) + 3.2)
    fig_h = max(6.4, 0.26 * len(genes) + 2.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
    ax.set_title(f"Top {len(genes)} DEG heatmap", pad=18, loc="center", fontweight="bold")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=38, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=12)
    ax.set_xlabel("Samples", labelpad=12)
    ax.set_ylabel("Genes", labelpad=12)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label(f"Row-scaled {transform_label(expression_scale)}", fontsize=16, labelpad=10)
    fig.tight_layout()
    save_plot(fig, out_base)


def assign_colors(samples, groups, case_label, control_label):
    """Map every group label seen in the metadata to a stable color.

    case/control keep the standard up/down palette; any additional group
    labels (e.g. extra cohorts present in the count matrix but not part of
    the current contrast) get an overflow palette so PCA/QC plots never
    silently drop samples.
    """
    extra_palette = ["#6FCF97", "#F2C94C", "#9B7ED9", "#3FBFB0", "#E08E45"]
    color_map = {case_label: UP, control_label: DOWN}
    extra_i = 0
    for s in samples:
        g = groups[s]
        if g not in color_map:
            color_map[g] = extra_palette[extra_i % len(extra_palette)]
            extra_i += 1
    return color_map


def write_pca(matrix, samples, groups, color_map, out_base, expression_scale="raw-count"):
    apply_style()
    genes = list(matrix.keys())
    # samples x genes matrix, centered per gene (column) so
    # PCA reflects deviation from each gene's mean across samples.
    data = np.array([[transform_value(matrix[g][s], expression_scale) for g in genes] for s in samples], dtype=float)
    data = data - data.mean(axis=0, keepdims=True)
    u, s_vals, _vt = np.linalg.svd(data, full_matrices=False)
    scores = u * s_vals
    total_var = float(np.sum(s_vals**2))
    var_pct = (s_vals**2) / total_var * 100 if total_var > 0 else np.zeros_like(s_vals)
    n_comp = scores.shape[1]
    pc1 = scores[:, 0] if n_comp > 0 else np.zeros(len(samples))
    pc2 = scores[:, 1] if n_comp > 1 else np.zeros(len(samples))
    pc1_pct = float(var_pct[0]) if len(var_pct) > 0 else 0.0
    pc2_pct = float(var_pct[1]) if len(var_pct) > 1 else 0.0

    fig, ax = plt.subplots(figsize=(8.6, 6.8))
    ax.set_facecolor("#FBFDFF")
    seen = set()
    for i, sample in enumerate(samples):
        g = groups[sample]
        color = color_map[g]
        label = g if g not in seen else None
        seen.add(g)
        ax.scatter(pc1[i], pc2[i], s=100, c=color, edgecolors="white", linewidths=0.9, label=label, zorder=3)
        ax.annotate(sample, (pc1[i], pc2[i]), textcoords="offset points", xytext=(7, 7), fontsize=11, color=TEXT)
    ax.axhline(0, color=GRID, linewidth=0.9, zorder=1)
    ax.axvline(0, color=GRID, linewidth=0.9, zorder=1)
    ax.set_title("PCA of samples", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% variance)", labelpad=12)
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% variance)", labelpad=12)
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    fig.tight_layout()
    save_plot(fig, out_base)
    return {"pc1": dict(zip(samples, pc1.tolist())), "pc1_pct": pc1_pct, "pc2_pct": pc2_pct}


def write_sample_correlation(matrix, samples, out_base, expression_scale="raw-count"):
    apply_style()
    corr = sample_correlation_matrix(matrix, samples, expression_scale=expression_scale)
    n = len(samples)
    fig_size = max(6.8, 0.62 * n + 3.2)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_title("Sample correlation heatmap", pad=18, loc="center", fontweight="bold")
    ax.set_xticks(range(n))
    ax.set_xticklabels(samples, rotation=38, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(samples)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(n):
        for j in range(n):
            ax.text(
                j,
                i,
                f"{corr[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=10,
                color="white" if abs(corr[i, j]) > 0.6 else TEXT,
            )
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Pearson r", fontsize=16, labelpad=10)
    fig.tight_layout()
    save_plot(fig, out_base)
    return corr


def sample_correlation_matrix(matrix, samples, expression_scale="raw-count"):
    genes = list(matrix.keys())
    data = np.array([[transform_value(matrix[g][s], expression_scale) for s in samples] for g in genes], dtype=float)
    if len(samples) < 2:
        return np.array([[1.0]])
    return np.corrcoef(data, rowvar=False)


def write_qc_summary(samples, groups, lib_sizes, detected, color_map, out_base):
    apply_style()
    colors = [color_map[groups[s]] for s in samples]
    fig_w = max(10.4, 1.05 * len(samples) + 6.0)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, 5.8))

    ax = axes[0]
    ax.bar(range(len(samples)), [lib_sizes[s] for s in samples], color=colors, edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=38, ha="right")
    ax.set_title("Library size", pad=14, loc="center", fontweight="bold")
    ax.set_ylabel("Total counts", labelpad=10)
    ax.grid(True, axis="y", color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    ax = axes[1]
    handles = []
    seen = set()
    for s in samples:
        g = groups[s]
        if g not in seen:
            seen.add(g)
            handles.append(plt.Rectangle((0, 0), 1, 1, color=color_map[g], label=g))
    ax.bar(range(len(samples)), [detected[s] for s in samples], color=colors, edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=38, ha="right")
    ax.set_title("Detected genes", pad=14, loc="center", fontweight="bold")
    ax.set_ylabel("Genes with count > 0", labelpad=10)
    ax.grid(True, axis="y", color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)

    fig.suptitle("Sample QC summary", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_plot(fig, out_base)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def load_and_validate_inputs(counts_path, metadata_path, group_col, case_label, control_label):
    """Load counts/metadata and fail fast with actionable messages on malformed input.

    Covers: missing required columns, duplicate sample/gene IDs, sample ID
    mismatches between the two files, non-numeric/negative count values, and
    empty case/control groups. See references/pitfalls.md and
    references/input-contract.md for the underlying rules.
    """
    meta_rows, meta_fields = read_table(metadata_path)
    if "sample" not in meta_fields:
        raise SystemExit(
            f"metadata file '{metadata_path}' is missing the required 'sample' column. "
            f"Found columns: {', '.join(meta_fields) or '(none)'}. Add a 'sample' column with IDs matching "
            "the count matrix header (see references/input-contract.md), or write output/questions.md if unclear."
        )
    if group_col not in meta_fields:
        other_cols = [c for c in meta_fields if c != "sample"]
        raise SystemExit(
            f"metadata file '{metadata_path}' has no column named '{group_col}' (--group-col). "
            f"Available non-sample columns: {', '.join(other_cols) or '(none)'}. "
            "Pass the correct --group-col, or write output/questions.md if the grouping is unclear."
        )
    if not meta_rows:
        raise SystemExit(f"metadata file '{metadata_path}' has a header but zero sample rows.")

    meta_sample_ids = [row["sample"] for row in meta_rows]
    dup_meta = sorted(s for s, n in Counter(meta_sample_ids).items() if n > 1)
    if dup_meta:
        raise SystemExit(
            f"metadata file '{metadata_path}' has duplicate 'sample' values: {', '.join(dup_meta)}. "
            "Each sample must appear exactly once; deduplicate or rename before re-running."
        )

    groups = {row["sample"]: row[group_col] for row in meta_rows}
    extra_covariate_cols = [c for c in meta_fields if c not in ("sample", group_col)]

    count_rows, count_fields = read_table(counts_path)
    if "gene" not in count_fields:
        raise SystemExit(
            f"count matrix '{counts_path}' is missing the required 'gene' column (first column). "
            f"Found columns: {', '.join(count_fields) or '(none)'}."
        )
    if not count_rows:
        raise SystemExit(f"count matrix '{counts_path}' has a header but zero gene rows.")
    samples = [k for k in count_fields if k != "gene"]
    if not samples:
        raise SystemExit(f"count matrix '{counts_path}' has a 'gene' column but no sample columns after it.")

    gene_ids = [row["gene"] for row in count_rows]
    dup_genes = sorted(g for g, n in Counter(gene_ids).items() if n > 1)
    if dup_genes:
        preview = ", ".join(dup_genes[:10]) + (", ..." if len(dup_genes) > 10 else "")
        raise SystemExit(
            f"count matrix '{counts_path}' has duplicate gene IDs: {preview}. "
            "Aggregate duplicates (e.g. sum counts per gene ID) before re-running so each gene has exactly one row."
        )

    matrix_only = sorted(set(samples) - set(groups))
    metadata_only = sorted(set(groups) - set(samples))
    if matrix_only or metadata_only:
        parts = []
        if matrix_only:
            parts.append(f"in count matrix but not metadata: {', '.join(matrix_only)}")
        if metadata_only:
            parts.append(f"in metadata but not count matrix: {', '.join(metadata_only)}")
        raise SystemExit(
            "Sample ID mismatch between counts and metadata (" + "; ".join(parts) + "). "
            "Sample IDs must match exactly, including case, between the count matrix header and metadata['sample']."
        )

    case_samples = [s for s, g in groups.items() if g == case_label]
    control_samples = [s for s, g in groups.items() if g == control_label]
    if not case_samples:
        raise SystemExit(
            f"No samples have {group_col} == '{case_label}' (--case). "
            f"Observed {group_col} values in metadata: {sorted(set(groups.values()))}."
        )
    if not control_samples:
        raise SystemExit(
            f"No samples have {group_col} == '{control_label}' (--control). "
            f"Observed {group_col} values in metadata: {sorted(set(groups.values()))}."
        )

    raw_matrix = {}
    for row in count_rows:
        gene = row["gene"]
        vals = {}
        for s in samples:
            raw = row[s]
            try:
                vals[s] = float(raw)
            except (TypeError, ValueError) as exc:
                raise SystemExit(
                    f"count matrix '{counts_path}' has a non-numeric value ('{raw}') for gene '{gene}', "
                    f"sample '{s}'. Fix or remove that cell before re-running."
                ) from exc
        raw_matrix[gene] = vals

    negative = [(g, s, v) for g, vals in raw_matrix.items() for s, v in vals.items() if v < 0]
    if negative:
        g, s, v = negative[0]
        raise SystemExit(
            f"count matrix '{counts_path}' contains negative values (e.g. gene '{g}', sample '{s}' = {v}). "
            "Counts must be non-negative; check upstream quantification/export before re-running "
            "(see references/input-contract.md)."
        )

    has_non_integer = any(abs(v - round(v)) > 1e-6 for vals in raw_matrix.values() for v in vals.values())

    return {
        "gene_order": gene_ids,
        "raw_matrix": raw_matrix,
        "samples": samples,
        "groups": groups,
        "case_samples": case_samples,
        "control_samples": control_samples,
        "extra_covariate_cols": extra_covariate_cols,
        "meta_rows": meta_rows,
        "has_non_integer": has_non_integer,
    }


def filter_detected_genes(gene_order, raw_matrix):
    """Drop genes with zero counts across every sample (minimal detection filter,
    see references/cns-method-guide.md 'Filter low-count genes before modeling')."""
    kept = [g for g in gene_order if any(v > 0 for v in raw_matrix[g].values())]
    if not kept:
        raise SystemExit(
            "Every gene has zero counts across all samples after filtering; there is nothing to test. "
            "Check that the count matrix was loaded/exported correctly."
        )
    matrix = {g: raw_matrix[g] for g in kept}
    return matrix, len(gene_order) - len(kept)


# ---------------------------------------------------------------------------
# QC verdict + pitfalls-driven caveats
# ---------------------------------------------------------------------------


def compute_qc_verdict(matrix, samples, groups, extra_covariate_cols, meta_rows, expression_scale="raw-count"):
    """Heuristic QC read on the same correlation matrix rendered in the sample
    correlation heatmap: (1) flag samples whose mean same-group correlation
    trails their group's average (possible outlier/mislabel/batch effect);
    (2) flag any extra metadata covariate that is perfectly confounded with
    the intended group (contingency check, pitfalls.md #13)."""
    corr = sample_correlation_matrix(matrix, samples, expression_scale=expression_scale)
    idx = {s: i for i, s in enumerate(samples)}

    group_members = {}
    for s in samples:
        group_members.setdefault(groups[s], []).append(s)

    outlier_notes = []
    for group_label, members in group_members.items():
        if len(members) < 2:
            outlier_notes.append(
                f"Group '{group_label}' has only {len(members)} sample(s); within-group correlation consistency "
                "cannot be assessed."
            )
            continue
        own_group_corr = {
            s: float(np.mean([corr[idx[s], idx[o]] for o in members if o != s])) for s in members
        }
        group_avg = float(np.mean(list(own_group_corr.values())))
        for s in members:
            gap = group_avg - own_group_corr[s]
            if gap > OUTLIER_CORR_GAP:
                outlier_notes.append(
                    f"Sample '{s}' (group '{group_label}') correlates with its own group at r={own_group_corr[s]:.2f}, "
                    f"{gap:.2f} below the group average (r={group_avg:.2f}) — possible outlier, mislabel, or batch effect."
                )

    confounded_covariates = []
    if extra_covariate_cols:
        cov_lookup = {row["sample"]: row for row in meta_rows if row.get("sample") in groups}
        for col in extra_covariate_cols:
            by_cov = {}
            for s in samples:
                v = cov_lookup.get(s, {}).get(col)
                by_cov.setdefault(v, set()).add(groups[s])
            if len(by_cov) < 2:
                continue  # constant covariate (zero variance) carries no separating information
            if all(len(gset) == 1 for gset in by_cov.values()):
                confounded_covariates.append(col)

    return {"outlier_notes": outlier_notes, "confounded_covariates": confounded_covariates}


def compute_caveats(n_case, n_control, has_non_integer, confounded_covariates, group_col, case_label, control_label,
                    padj_threshold=0.05, log2fc_threshold=1.0, expression_scale="raw-count"):
    """Check the run's actual data against references/pitfalls.md and return only
    the caveats that are triggered — never a static list of all pitfalls.

    Each caveat is a (title, detail) pair so callers can render a full bullet
    and also build a short conclusion clause without re-parsing prose."""
    caveats = []
    if n_case < SMALL_N_THRESHOLD or n_control < SMALL_N_THRESHOLD:
        caveats.append((
            "underpowered design",
            f"**Underpowered design** (pitfalls.md #10): n={n_case} ({case_label}) vs n={n_control} ({control_label}) "
            f"is below the recommended minimum of {SMALL_N_THRESHOLD} biological replicates per group. Dispersion/variance "
            "estimates are unreliable at this replication level; treat fold-changes and rankings here as exploratory "
            "hypothesis-generation, not confirmed differential expression.",
        ))
    if n_case > LARGE_N_THRESHOLD and n_control > LARGE_N_THRESHOLD:
        caveats.append((
            "effect-size floor at large n",
            f"**Effect-size floor at large n** (pitfalls.md #12): with n={n_case}/{n_control} per group, statistical "
            f"significance can be reached for biologically trivial fold-changes. The abs(log2FC) >= {log2fc_threshold:g} filter is kept "
            f"as a hard floor alongside padj < {padj_threshold:g}; consider ranking by effect size rather than padj alone.",
        ))
    for col in confounded_covariates:
        caveats.append((
            f"'{col}' confounded with '{group_col}'",
            f"**Confounded design** (pitfalls.md #13): metadata column '{col}' is perfectly confounded with "
            f"'{group_col}' — every '{col}' value maps to exactly one group. The effect of '{col}' and '{group_col}' "
            "cannot be statistically disentangled with this design; do not attribute observed differences to "
            f"'{group_col}' alone.",
        ))
    if has_non_integer:
        if expression_scale == "log-expression":
            scale_note = (
                "The quicklook route detected small continuous non-integer values and used them as already "
                "log/normalized expression, avoiding a second `log2(x+1)` transform."
            )
        else:
            scale_note = (
                "The quicklook route treated these as linear continuous abundance/intensity values and used "
                "`log2(x+1)` for quicklook statistics/plots."
            )
        caveats.append((
            "non-integer count values",
            "**Non-integer count values** (pitfalls.md #11): the input matrix contains non-integer values, "
            "suggesting it may already be normalized (TPM/FPKM/CPM/logCPM/log expression) rather than raw counts. "
            f"{scale_note} Do not run `run_deseq2.R`/edgeR on this matrix — those require raw integer counts. "
            "Use limma on log-scale values for a formal model instead.",
        ))
    return caveats


def build_report(ctx):
    """Assemble a results-and-conclusions report.md: input summary, key findings
    (named top genes, not just counts), a computed QC verdict, only the
    pitfalls.md caveats actually triggered by this run's data, and a
    plain-language conclusion — not a file/parameter listing."""
    lines = ["# bulkRNA quicklook report", ""]
    lines.append(
        f"Exploratory analysis for **{ctx['case']}** vs **{ctx['control']}** ({ctx['control']} = reference/control). "
        f"This is the Python quicklook route (`quickplot_bulk.py`), using **{ctx['expression_scale']}** scale "
        f"({ctx['transform_description']}) for quicklook statistics and plots: approximate mean-difference effect sizes and "
        "a normal-approximation p value, not a formal negative-binomial (DESeq2/edgeR) or linear (limma) model. "
        "Use `run_deseq2.R` (or `run_limma_voom.R` for multi-factor/batch/paired designs) for confirmatory statistics "
        "when raw counts and Bioconductor are available."
    )
    lines.append("")

    lines.append("## Input summary")
    lines.append("")
    lines.append(f"- Samples: {len(ctx['case_samples'])} {ctx['case']} vs {len(ctx['control_samples'])} {ctx['control']}")
    filt = ctx["filtered_out"]
    lines.append(
        f"- Genes in input matrix: {ctx['total_genes_input']}; genes tested after removing {filt} "
        f"undetected (all-zero) gene(s): {ctx['genes_tested']}"
    )
    lines.append(
        f"- Contrast: {ctx['case']} vs {ctx['control']}; thresholds: padj < {ctx['padj_threshold']:g} and "
        f"abs(log2FC) >= {ctx['log2fc_threshold']:g}"
    )
    lines.append(f"- Expression scale handling: {ctx['expression_scale']} ({ctx['transform_description']})")
    if ctx["extra_covariate_cols"]:
        lines.append(f"- Additional metadata columns present (checked for confounding): {', '.join(ctx['extra_covariate_cols'])}")
    lines.append("")

    lines.append("## Key findings")
    lines.append("")
    lines.append(f"- Significant up in {ctx['case']}: {ctx['up']} gene(s)")
    lines.append(f"- Significant down in {ctx['case']}: {ctx['down']} gene(s)")
    # "Top genes" means significant (up/down) genes, not just small-padj genes -- a gene can have
    # a tiny padj but sit below the log2FC floor (change == "ns") and would otherwise misleadingly
    # rank #1 in a "top gene" list despite not clearing the effect-size threshold.
    sig_results = [r for r in ctx["results"] if r["change"] != "ns"]
    top_n = min(5, len(sig_results))
    if top_n:
        lines.append(f"- Top {top_n} gene(s) by adjusted p value:")
        for rank, r in enumerate(sig_results[:top_n], start=1):
            lines.append(f"  {rank}. `{r['gene']}`: log2FC = {r['log2FC']:.2f}, padj = {r['padj']:.3g} ({r['change']})")
    lines.append("")

    lines.append("## QC verdict")
    lines.append("")
    verdict = ctx["qc_verdict"]
    if verdict["outlier_notes"]:
        for note in verdict["outlier_notes"]:
            lines.append(f"- {note}")
    else:
        lines.append("- No sample fell outside the same-group correlation heuristic; no outlier flagged.")
    if verdict["confounded_covariates"]:
        for col in verdict["confounded_covariates"]:
            lines.append(f"- Covariate '{col}' is fully confounded with the group of interest (see caveats below).")
    else:
        if ctx["extra_covariate_cols"]:
            lines.append("- No additional metadata covariate is fully confounded with the group of interest.")
    lines.append("")

    lines.append("## Caveats (automatically checked against references/pitfalls.md)")
    lines.append("")
    if ctx["caveats"]:
        for _title, detail in ctx["caveats"]:
            lines.append(f"- {detail}")
    else:
        lines.append(
            "- No data-triggered caveats from the automatic checks (design >= "
            f"{SMALL_N_THRESHOLD}/group, no full covariate confounding detected, counts are integer-valued). "
            "This does not replace manual review of `references/pitfalls.md`."
        )
    lines.append("")

    lines.append("## Conclusion")
    lines.append("")
    n_case, n_control = len(ctx["case_samples"]), len(ctx["control_samples"])
    qc_sentence = (
        "PCA and sample correlation show no evidence of outlier samples or covariate confounding"
        if not verdict["outlier_notes"] and not verdict["confounded_covariates"]
        else "PCA/sample correlation flagged one or more issues noted above"
    )
    robustness = (
        "exploratory: this quicklook route uses approximate statistics and should be confirmed with the formal "
        "DESeq2/limma route before drawing biological conclusions"
    )
    caveat_clause = (
        f", with the following flagged concerns: {', '.join(title for title, _detail in ctx['caveats'])}"
        if ctx["caveats"]
        else ""
    )
    lines.append(
        f"{ctx['up']} gene(s) were significantly upregulated and {ctx['down']} downregulated in {ctx['case']} vs "
        f"{ctx['control']} (n={n_case}/{n_control}) at padj<{ctx['padj_threshold']:g} and "
        f"abs(log2FC)>={ctx['log2fc_threshold']:g}. {qc_sentence}. "
        f"Results should be treated as {robustness}{caveat_clause}."
    )
    lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append(
        "Figures: `volcano`, `heatmap_top_degs`, `pca`, `sample_correlation_heatmap`, `qc_summary` "
        "(library size and detected genes), each exported as PNG and PDF. Full gene table: "
        "`results/differential_expression.tsv`. Per-sample QC: `qc/sample_qc.tsv`."
    )
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="bulk RNA-seq quicklook plots with PNG/PDF export.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--group-col", default="group")
    parser.add_argument("--case", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--out", default="output")
    parser.add_argument(
        "--padj-threshold", type=float, default=0.05,
        help="Adjusted p-value (BH FDR) cutoff for calling a gene significant (default 0.05).",
    )
    parser.add_argument(
        "--log2fc-threshold", type=float, default=1.0,
        help="Absolute log2 fold-change cutoff for calling a gene significant (default 1.0).",
    )
    parser.add_argument(
        "--top-n-heatmap", type=int, default=30,
        help="Number of top genes (by adjusted p-value) shown in the top-DEG heatmap (default 30).",
    )
    parser.add_argument(
        "--expression-scale",
        choices=["auto", "raw-count", "log-expression", "linear-expression"],
        default="auto",
        help=(
            "How to transform input values for quicklook statistics/plots. auto: integer matrices use log2(count+1); "
            "small non-integer matrices are treated as already log/normalized expression and are not logged again; "
            "large non-integer abundance/intensity matrices use log2(value+1)."
        ),
    )
    args = parser.parse_args()

    padj_threshold = args.padj_threshold
    log2fc_threshold = args.log2fc_threshold

    out = Path(args.out)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "qc").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    loaded = load_and_validate_inputs(args.counts, args.metadata, args.group_col, args.case, args.control)
    samples = loaded["samples"]
    groups = loaded["groups"]
    case_samples = loaded["case_samples"]
    control_samples = loaded["control_samples"]
    total_genes_input = len(loaded["gene_order"])

    matrix, filtered_out = filter_detected_genes(loaded["gene_order"], loaded["raw_matrix"])
    expression_scale = infer_expression_scale(loaded["raw_matrix"], loaded["has_non_integer"], args.expression_scale)
    transform_description = transform_label(expression_scale)

    results = []
    p_pairs = []
    for gene in matrix:
        cvals = [transform_value(matrix[gene][s], expression_scale) for s in control_samples]
        tvals = [transform_value(matrix[gene][s], expression_scale) for s in case_samples]
        log2fc = statistics.mean(tvals) - statistics.mean(cvals)
        spread = statistics.pstdev(tvals + cvals) or 1
        p = normal_p_from_effect(log2fc, spread, len(tvals) + len(cvals))
        results.append({"gene": gene, "log2FC": log2fc, "pvalue": p})
        p_pairs.append((gene, p))

    padj = bh_adjust(p_pairs)
    for row in results:
        row["padj"] = padj[row["gene"]]
        row["change"] = (
            "up" if row["log2FC"] >= log2fc_threshold and row["padj"] < padj_threshold
            else "down" if row["log2FC"] <= -log2fc_threshold and row["padj"] < padj_threshold
            else "ns"
        )
    results.sort(key=lambda r: (r["padj"], -abs(r["log2FC"])))

    with (out / "results" / "differential_expression.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "log2FC", "pvalue", "padj", "change"], delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow({k: f"{v:.6g}" if isinstance(v, float) else v for k, v in row.items()})

    write_volcano(results, out / "figures" / "volcano", padj_threshold=padj_threshold, log2fc_threshold=log2fc_threshold)
    write_heatmap(
        matrix,
        samples,
        results,
        out / "figures" / "heatmap_top_degs",
        top_n=args.top_n_heatmap,
        expression_scale=expression_scale,
    )

    # QC summary: library size, detected genes, PCA, and sample correlation
    # (see references/cns-method-guide.md Required Figures / advanced-modules.md).
    lib_sizes = {s: sum(matrix[gene][s] for gene in matrix) for s in samples}
    detected = {s: sum(1 for gene in matrix if matrix[gene][s] > 0) for s in samples}
    color_map = assign_colors(samples, groups, args.case, args.control)
    write_pca(matrix, samples, groups, color_map, out / "figures" / "pca", expression_scale=expression_scale)
    write_sample_correlation(
        matrix,
        samples,
        out / "figures" / "sample_correlation_heatmap",
        expression_scale=expression_scale,
    )
    write_qc_summary(samples, groups, lib_sizes, detected, color_map, out / "figures" / "qc_summary")

    qc_verdict = compute_qc_verdict(
        matrix,
        samples,
        groups,
        loaded["extra_covariate_cols"],
        loaded["meta_rows"],
        expression_scale=expression_scale,
    )
    caveats = compute_caveats(
        len(case_samples),
        len(control_samples),
        loaded["has_non_integer"],
        qc_verdict["confounded_covariates"],
        args.group_col,
        args.case,
        args.control,
        padj_threshold=padj_threshold,
        log2fc_threshold=log2fc_threshold,
        expression_scale=expression_scale,
    )

    up = sum(1 for r in results if r["change"] == "up")
    down = sum(1 for r in results if r["change"] == "down")
    (out / "summary.tsv").write_text(
        "item\tvalue\n"
        f"case_samples\t{len(case_samples)}\n"
        f"control_samples\t{len(control_samples)}\n"
        f"genes_input\t{total_genes_input}\n"
        f"genes_filtered_out\t{filtered_out}\n"
        f"genes_tested\t{len(results)}\n"
        f"up_genes\t{up}\n"
        f"down_genes\t{down}\n"
        f"padj_threshold\t{padj_threshold:g}\n"
        f"log2fc_threshold\t{log2fc_threshold:g}\n"
        f"top_n_heatmap\t{args.top_n_heatmap}\n"
        f"expression_scale\t{expression_scale}\n"
        f"transform\t{transform_description}\n",
        encoding="utf-8",
    )
    (out / "qc" / "sample_qc.tsv").write_text(
        "sample\tgroup\tlibrary_size\tdetected_genes\n"
        + "\n".join(f"{s}\t{groups[s]}\t{int(lib_sizes[s])}\t{detected[s]}" for s in samples)
        + "\n",
        encoding="utf-8",
    )
    artifacts = [
        "summary.tsv",
        "results/differential_expression.tsv",
        "qc/sample_qc.tsv",
        "figures/volcano.png",
        "figures/volcano.pdf",
        "figures/heatmap_top_degs.png",
        "figures/heatmap_top_degs.pdf",
        "figures/pca.png",
        "figures/pca.pdf",
        "figures/sample_correlation_heatmap.png",
        "figures/sample_correlation_heatmap.pdf",
        "figures/qc_summary.png",
        "figures/qc_summary.pdf",
    ]
    (out / "run" / "artifacts.json").write_text(json.dumps({"artifacts": [{"path": p, "exists": (out / p).exists()} for p in artifacts]}, indent=2), encoding="utf-8")
    (out / "run" / "qc.json").write_text(
        json.dumps(
            {
                "genes_input": total_genes_input,
                "genes_filtered_out": filtered_out,
                "genes_tested": len(results),
                "samples": len(samples),
                "case": args.case,
                "control": args.control,
                "exploratory": True,
                "outlier_notes": qc_verdict["outlier_notes"],
                "confounded_covariates": qc_verdict["confounded_covariates"],
                "caveats_triggered": len(caveats),
                "padj_threshold": padj_threshold,
                "log2fc_threshold": log2fc_threshold,
                "top_n_heatmap": args.top_n_heatmap,
                "expression_scale": expression_scale,
                "transform": transform_description,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    report_ctx = {
        "case": args.case,
        "control": args.control,
        "case_samples": case_samples,
        "control_samples": control_samples,
        "total_genes_input": total_genes_input,
        "filtered_out": filtered_out,
        "genes_tested": len(results),
        "results": results,
        "up": up,
        "down": down,
        "qc_verdict": qc_verdict,
        "caveats": caveats,
        "extra_covariate_cols": loaded["extra_covariate_cols"],
        "padj_threshold": padj_threshold,
        "log2fc_threshold": log2fc_threshold,
        "expression_scale": expression_scale,
        "transform_description": transform_description,
    }
    (out / "report.md").write_text(build_report(report_ctx), encoding="utf-8")
    print(f"wrote quicklook outputs to {out}")


if __name__ == "__main__":
    main()
