#!/usr/bin/env python3
"""CRISPR screen quick-analysis plots and result tables.

Primary analysis route when no R/MAGeCK installation is available (see
references/cns-method-guide.md and SKILL.md). Always required:
  --gene-summary  a MAGeCK-style or generic gene summary (gene, logFC/beta,
                  pvalue, fdr/padj).

Optional, unlocks QC and sgRNA-consistency figures:
  --sgrna-counts  an sgRNA-level raw count table (sgRNA, gene, >=2 numeric
                  sample columns; column names containing control/ctrl/t0/
                  plasmid/dmso/vehicle are treated as the control group, the
                  rest as treatment).

Figures follow references/visual-style.md (PNG+PDF only, centered bold
titles, sparse staggered labels, POS/NEG/NS palette).
"""
import argparse
import csv
import json
import math
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
except Exception as exc:
    raise SystemExit("Install matplotlib: python3 -m pip install matplotlib") from exc

try:
    import numpy as np
except Exception as exc:
    raise SystemExit("Install numpy: python3 -m pip install numpy") from exc

TEXT = "#173B5F"
GRID = "#DCEAF5"
POS = "#E45B5B"
NEG = "#4EA3F1"
NS = "#B8C2CC"
ACCENT = "#A7C6DD"

SCORE_CANDIDATES = ["logfc", "lfc", "beta", "score", "log2fc"]
FDR_CANDIDATES = ["fdr", "padj", "adj.pvalue", "adj_pvalue", "qvalue", "q.value"]
GENE_CANDIDATES = ["gene", "id", "symbol", "genes", "target_gene"]
CONTROL_HINTS = ("control", "ctrl", "t0", "plasmid", "dmso", "vehicle", "baseline")


def read(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def pick_column(fieldnames, candidates, label, required=True):
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c in lower:
            return lower[c]
    if required:
        raise SystemExit(f"gene summary is missing a {label} column (looked for {candidates}); found {list(fieldnames)}")
    return None


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 16,
            "axes.titlesize": 24,
            "axes.labelsize": 21,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 13,
            "axes.edgecolor": "#B9D7EE",
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "text.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def classify(score, fdr, lfc_threshold, fdr_threshold):
    if score > lfc_threshold and fdr < fdr_threshold:
        return "pos"
    if score < -lfc_threshold and fdr < fdr_threshold:
        return "neg"
    return "ns"


def _boxes_overlap(a, b, pad=3.0):
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    return not (ax1 + pad < bx0 or bx1 + pad < ax0 or ay1 + pad < by0 or by1 + pad < ay0)


def add_staggered_labels(ax, points, fontsize=11):
    """Label points without an external adjustText dependency.

    For each point (in the order given), try a ring of candidate offsets
    (near -> far, alternating above/below and left/right/centered) and place
    the label at the first candidate whose estimated bounding box (in pixel
    space) does not overlap any label already placed. This keeps dense hit
    clusters legible instead of stacking labels on top of each other
    (references/visual-style.md: stagger annotations so labels never collide).
    """
    if not points:
        return
    fig = ax.figure
    fig.canvas.draw()
    dpi = fig.dpi
    pt_to_px = dpi / 72.0

    candidates = []
    for k in range(6):
        base = 13 + 12 * k
        candidates.extend([(0, base), (0, -base), (18, base), (18, -base), (-18, base), (-18, -base)])

    placed_boxes = []
    for x, y, text, color in points:
        cx0, cy0 = ax.transData.transform((x, y))
        w_px = len(text) * fontsize * 0.62 * pt_to_px
        h_px = fontsize * 1.4 * pt_to_px
        chosen = candidates[-1]
        chosen_box = None
        for dx, dy in candidates:
            cx = cx0 + dx * pt_to_px
            cy = cy0 + dy * pt_to_px
            bx0, bx1 = (cx - w_px / 2, cx + w_px / 2) if dx == 0 else ((cx, cx + w_px) if dx > 0 else (cx - w_px, cx))
            by0, by1 = (cy, cy + h_px) if dy >= 0 else (cy - h_px, cy)
            box = (bx0, by0, bx1, by1)
            if not any(_boxes_overlap(box, pb) for pb in placed_boxes):
                chosen, chosen_box = (dx, dy), box
                break
        if chosen_box is None:
            dx, dy = chosen
            cx, cy = cx0 + dx * pt_to_px, cy0 + dy * pt_to_px
            chosen_box = (cx - w_px / 2, cy, cx + w_px / 2, cy + h_px) if dy >= 0 else (cx - w_px / 2, cy - h_px, cx + w_px / 2, cy)
        placed_boxes.append(chosen_box)
        dx, dy = chosen
        ha = "center" if dx == 0 else ("left" if dx > 0 else "right")
        va = "bottom" if dy >= 0 else "top"
        ax.annotate(
            text,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=fontsize,
            color=color,
            arrowprops=dict(arrowstyle="-", color=color, lw=0.7, alpha=0.65),
        )


def top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, direction):
    sig = [r for r in rows if classify(float(r[score_key]), float(r[fdr_key]), lfc_threshold, fdr_threshold) == direction]
    reverse = direction == "pos"
    sig.sort(key=lambda r: float(r[score_key]), reverse=reverse)
    return sig[:top_n]


def gini(values):
    values = np.sort(np.asarray([v for v in values if v >= 0], dtype=float))
    n = values.size
    total = values.sum()
    if n == 0 or total == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return float((2 * np.sum(index * values)) / (n * total) - (n + 1) / n)


def split_groups(sample_cols):
    control = [c for c in sample_cols if any(k in c.lower() for k in CONTROL_HINTS)]
    treatment = [c for c in sample_cols if c not in control]
    if not control and len(sample_cols) >= 2:
        half = len(sample_cols) // 2
        control, treatment = sample_cols[:half], sample_cols[half:]
    return control, treatment


def plot_volcano(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, out):
    xs = np.array([float(r[score_key]) for r in rows])
    ys = np.array([-math.log10(max(float(r[fdr_key]), 1e-12)) for r in rows])
    cats = [classify(x, float(r[fdr_key]), lfc_threshold, fdr_threshold) for x, r in zip(xs, rows)]
    colors = [{"pos": POS, "neg": NEG, "ns": NS}[c] for c in cats]

    fig, ax = plt.subplots(figsize=(9.2, 6.8))
    ax.scatter(xs, ys, c=colors, s=42, alpha=0.86, linewidths=0, zorder=2)
    ax.axvline(0, color=ACCENT, linewidth=1.1)
    ax.axvline(lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.axvline(-lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.axhline(-math.log10(fdr_threshold), color=ACCENT, linestyle="--", linewidth=1.1)
    ax.set_title("CRISPR screen volcano", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("log fold change (beta / logFC)", labelpad=12)
    ax.set_ylabel("-log10 FDR", labelpad=12)
    ax.grid(True, color=GRID, zorder=0)

    pos_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "pos")
    neg_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "neg")
    points = [(float(r[score_key]), -math.log10(max(float(r[fdr_key]), 1e-12)), r[gene_key], POS) for r in pos_hits]
    points += [(float(r[score_key]), -math.log10(max(float(r[fdr_key]), 1e-12)), r[gene_key], NEG) for r in neg_hits]
    add_staggered_labels(ax, points)

    n_pos, n_neg, n_ns = cats.count("pos"), cats.count("neg"), cats.count("ns")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=POS, markersize=9, label=f"positive selection (n={n_pos})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NEG, markersize=9, label=f"negative selection (n={n_neg})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NS, markersize=9, label=f"not significant (n={n_ns})"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save(fig, out / "figures" / "crispr_volcano")
    return n_pos, n_neg, n_ns


def plot_gene_rank(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, out):
    ranked = sorted(rows, key=lambda r: float(r[score_key]))
    rank_x = np.arange(len(ranked))
    rank_y = np.array([float(r[score_key]) for r in ranked])
    cats = [classify(y, float(r[fdr_key]), lfc_threshold, fdr_threshold) for y, r in zip(rank_y, ranked)]
    colors = [{"pos": POS, "neg": NEG, "ns": NS}[c] for c in cats]

    fig, ax = plt.subplots(figsize=(9.6, 6.6))
    ax.plot(rank_x, rank_y, color="#C9D6E3", lw=1.6, zorder=1)
    ax.scatter(rank_x, rank_y, c=colors, s=26, zorder=2, linewidths=0)
    ax.axhline(0, color=ACCENT, linewidth=1.0)
    ax.axhline(lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.axhline(-lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.set_title("Gene rank plot", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Gene rank (by score, low to high)", labelpad=12)
    ax.set_ylabel("log fold change (beta / logFC)", labelpad=12)
    ax.grid(True, color=GRID, zorder=0)

    pos_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "pos")
    neg_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "neg")
    order_index = {id(r): i for i, r in enumerate(ranked)}
    points = []
    for r in pos_hits:
        points.append((order_index[id(r)], float(r[score_key]), r[gene_key], POS))
    for r in neg_hits:
        points.append((order_index[id(r)], float(r[score_key]), r[gene_key], NEG))
    add_staggered_labels(ax, points)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=POS, markersize=9, label="positive selection"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NEG, markersize=9, label="negative selection"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NS, markersize=9, label="not significant"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    save(fig, out / "figures" / "gene_rank_plot")


def plot_top_hits_bar(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, out):
    pos_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "pos")
    neg_hits = top_hit_rows(rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "neg")
    if not pos_hits and not neg_hits:
        return None
    combined = list(reversed(neg_hits)) + list(reversed(pos_hits))
    labels = [r[gene_key] for r in combined]
    values = [float(r[score_key]) for r in combined]
    colors = [NEG if v < 0 else POS for v in values]

    fig, ax = plt.subplots(figsize=(8.6, max(4.5, 0.42 * len(combined) + 1.6)))
    y_pos = np.arange(len(combined))
    ax.barh(y_pos, values, color=colors, height=0.62, zorder=2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.axvline(0, color=ACCENT, linewidth=1.1)
    ax.set_title("Top positive/negative selection hits", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("log fold change (beta / logFC)", labelpad=12)
    ax.grid(True, axis="x", color=GRID, zorder=0)
    ax.margins(x=0.32)  # leave room for the FDR annotation beyond each bar tip so it never collides with the y-axis gene labels
    for y, v, r in zip(y_pos, values, combined):
        fdr = float(r[fdr_key])
        ax.annotate(
            f"FDR={fdr:.1e}",
            xy=(v, y),
            xytext=(6 if v >= 0 else -6, 0),
            textcoords="offset points",
            ha="left" if v >= 0 else "right",
            va="center",
            fontsize=9.5,
            color=TEXT,
        )
    save(fig, out / "figures" / "top_hits_barplot")
    return combined


def read_sgrna_counts(path):
    rows = read(path)
    if not rows:
        raise SystemExit(f"{path} has no data rows")
    fieldnames = list(rows[0].keys())
    lower = {f.lower(): f for f in fieldnames}
    sgrna_col = lower.get("sgrna") or lower.get("guide") or fieldnames[0]
    gene_col = lower.get("gene") or lower.get("target") or lower.get("target_gene") or fieldnames[1]
    sample_cols = [f for f in fieldnames if f not in (sgrna_col, gene_col)]
    return rows, sgrna_col, gene_col, sample_cols


def plot_qc(sg_rows, sgrna_col, gene_col, sample_cols, out):
    control_cols, treatment_cols = split_groups(sample_cols)
    samples = {c: np.array([float(r[c]) for r in sg_rows]) for c in sample_cols}
    gini_by_sample = {c: gini(v) for c, v in samples.items()}
    is_ntc = np.array([r[gene_col].upper() == "NTC" for r in sg_rows])
    guides_per_gene = {}
    for r in sg_rows:
        if r[gene_col].upper() == "NTC":
            continue
        guides_per_gene[r[gene_col]] = guides_per_gene.get(r[gene_col], 0) + 1

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))

    # (a) read count distribution per sample
    ax = axes[0, 0]
    data = [np.log10(samples[c] + 1) for c in sample_cols]
    try:
        bp = ax.boxplot(data, tick_labels=sample_cols, patch_artist=True, widths=0.55)
    except TypeError:  # matplotlib < 3.9 does not know tick_labels
        bp = ax.boxplot(data, labels=sample_cols, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], sample_cols):
        patch.set_facecolor(NEG if c in control_cols else POS)
        patch.set_alpha(0.55)
    ax.set_title("Read count distribution", loc="center", pad=14, color=TEXT, fontweight="bold")
    ax.set_ylabel("log10(count + 1)", labelpad=10)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(True, axis="y", color=GRID, zorder=0)

    # (b) sgRNAs per gene: a discrete bar chart (not a histogram) since guide
    # counts per gene are small integers; avoids fractional bin-edge ticks
    # when most/all genes share the same guide count.
    ax = axes[0, 1]
    counts = list(guides_per_gene.values())
    if counts:
        dist = {}
        for c in counts:
            dist[c] = dist.get(c, 0) + 1
        xs_sorted = sorted(dist)
        ax.bar(xs_sorted, [dist[k] for k in xs_sorted], color=ACCENT, edgecolor="white", width=0.6, zorder=2)
        ax.set_xticks(xs_sorted)
        # With a single guide-count category (common when every gene has the same
        # number of guides), matplotlib's default autoscale margin hugs the one bar
        # so tightly it renders as a full-width color block; force a fixed pad so
        # a lone bar still reads as a bar chart, not a solid fill.
        if len(xs_sorted) == 1:
            ax.set_xlim(xs_sorted[0] - 1.5, xs_sorted[0] + 1.5)
    ax.set_title("sgRNAs per gene", loc="center", pad=14, color=TEXT, fontweight="bold")
    ax.set_xlabel("sgRNAs per gene", labelpad=10)
    ax.set_ylabel("number of genes", labelpad=10)
    ax.grid(True, axis="y", color=GRID, zorder=0)

    # (c) replicate correlation (control group)
    ax = axes[1, 0]
    if len(control_cols) >= 2:
        a, b = np.log10(samples[control_cols[0]] + 1), np.log10(samples[control_cols[1]] + 1)
        r = float(np.corrcoef(a, b)[0, 1])
        ax.scatter(a, b, s=16, color=NEG, alpha=0.5, linewidths=0)
        lim = [min(a.min(), b.min()), max(a.max(), b.max())]
        ax.plot(lim, lim, color=ACCENT, linestyle="--", linewidth=1.1)
        ax.set_xlabel(f"log10({control_cols[0]} + 1)", labelpad=10)
        ax.set_ylabel(f"log10({control_cols[1]} + 1)", labelpad=10)
        ax.annotate(f"Pearson r = {r:.3f}", xy=(0.05, 0.92), xycoords="axes fraction", fontsize=13, color=TEXT)
    else:
        ax.text(0.5, 0.5, "Need >=2 control replicates\nfor correlation", ha="center", va="center", fontsize=13, color=TEXT)
        ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Control replicate correlation", loc="center", pad=14, color=TEXT, fontweight="bold")
    ax.grid(True, color=GRID, zorder=0)

    # (d) Gini index per sample (library evenness)
    ax = axes[1, 1]
    labels = sample_cols
    vals = [gini_by_sample[c] for c in labels]
    colors = [NEG if c in control_cols else POS for c in labels]
    ax.bar(np.arange(len(labels)), vals, color=colors, alpha=0.75)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.axhline(0.2, color=ACCENT, linestyle="--", linewidth=1.0)
    ax.set_title("Gini index (library evenness)", loc="center", pad=14, color=TEXT, fontweight="bold")
    ax.set_ylabel("Gini index", labelpad=10)
    ax.grid(True, axis="y", color=GRID, zorder=0)

    fig.suptitle("sgRNA count QC", fontsize=22, fontweight="bold", color=TEXT, y=1.02)
    fig.tight_layout()
    save(fig, out / "figures" / "qc_distribution")

    return {
        "control_columns": control_cols,
        "treatment_columns": treatment_cols,
        "gini_by_sample": gini_by_sample,
        "mean_guides_per_gene": float(np.mean(list(guides_per_gene.values()))) if guides_per_gene else None,
        "n_genes_with_guides": len(guides_per_gene),
        "n_sgrnas": len(sg_rows),
        "n_ntc_guides": int(is_ntc.sum()),
        "replicate_pearson_r_control": (
            float(np.corrcoef(np.log10(samples[control_cols[0]] + 1), np.log10(samples[control_cols[1]] + 1))[0, 1])
            if len(control_cols) >= 2
            else None
        ),
    }


def plot_sgrna_consistency(sg_rows, sgrna_col, gene_col, sample_cols, gene_rows, gene_key, score_key, fdr_key,
                            lfc_threshold, fdr_threshold, top_n, out):
    control_cols, treatment_cols = split_groups(sample_cols)
    if not control_cols or not treatment_cols:
        return None
    pseudo = 5.0
    per_guide = {}
    for r in sg_rows:
        gene = r[gene_col]
        if gene.upper() == "NTC":
            continue
        c_mean = np.mean([float(r[c]) for c in control_cols])
        t_mean = np.mean([float(r[c]) for c in treatment_cols])
        lfc = math.log2((t_mean + pseudo) / (c_mean + pseudo))
        per_guide.setdefault(gene, []).append((r[sgrna_col], lfc))

    pos_hits = top_hit_rows(gene_rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "pos")
    neg_hits = top_hit_rows(gene_rows, score_key, fdr_key, gene_key, lfc_threshold, fdr_threshold, top_n, "neg")
    hit_genes = [(r[gene_key], POS) for r in pos_hits] + [(r[gene_key], NEG) for r in neg_hits]
    hit_genes = [(g, c) for g, c in hit_genes if g in per_guide]
    if not hit_genes:
        return None

    fig, ax = plt.subplots(figsize=(10.2, max(4.2, 0.5 * len(hit_genes) + 1.6)))
    y_labels = []
    for i, (gene, color) in enumerate(hit_genes):
        guides = per_guide[gene]
        xs = [lfc for _, lfc in guides]
        ys = [i + np.random.uniform(-0.12, 0.12) for _ in guides]
        ax.scatter(xs, ys, color=color, s=48, alpha=0.75, zorder=3, edgecolors="white", linewidths=0.6)
        ax.scatter([np.mean(xs)], [i], color=TEXT, marker="D", s=42, zorder=4)
        y_labels.append(gene)
    ax.set_yticks(range(len(hit_genes)))
    ax.set_yticklabels(y_labels)
    ax.axvline(0, color=ACCENT, linewidth=1.0)
    ax.axvline(lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.axvline(-lfc_threshold, color=ACCENT, linestyle=":", linewidth=1.0)
    ax.set_title("sgRNA-level consistency for top hits", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("per-guide log2 fold change (treatment vs control)", labelpad=12)
    ax.grid(True, axis="x", color=GRID, zorder=0)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NS, markersize=8, label="individual sgRNA"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=TEXT, markersize=7, label="gene mean"),
    ]
    # Placed outside the axes (not "lower right"/"best"): with hit genes spanning
    # both directions, an in-plot legend corner can sit directly on top of the
    # bottom- or top-most gene's data points (e.g. the last plotted gene at y=0).
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save(fig, out / "figures" / "sgrna_consistency")

    summary_rows = []
    for gene, guides in per_guide.items():
        xs = [lfc for _, lfc in guides]
        n_guides = len(guides)
        mean_lfc = float(np.mean(xs))
        sign_mean = 1 if mean_lfc >= 0 else -1
        active_threshold = lfc_threshold / 2.0
        # A guide "concordant" with the gene-level call: same sign as the gene's
        # mean per-guide log2FC and with a non-trivial magnitude (not noise-floor).
        # See references/pitfalls.md #6: a gene-level hit driven by 1 of N guides
        # (n_concordant_guides <= 1 with n_guides >= 2) is weaker evidence.
        n_concordant = sum(1 for x in xs if ((1 if x >= 0 else -1) == sign_mean) and abs(x) >= active_threshold)
        summary_rows.append(
            {
                "gene": gene,
                "n_guides": n_guides,
                "mean_log2fc": round(mean_lfc, 4),
                "sd_log2fc": round(float(np.std(xs)), 4),
                "min_log2fc": round(float(np.min(xs)), 4),
                "max_log2fc": round(float(np.max(xs)), 4),
                "n_concordant_guides": n_concordant,
                "single_guide_driven": n_guides >= 2 and n_concordant <= 1,
            }
        )
    return summary_rows


def get_consistency_note(gene, has_sgrna, consistency_lookup):
    """One short clause appended after a named top hit describing whether its
    sgRNA-level signal is guide-consistent or single-guide-driven (see
    references/pitfalls.md #6). Returns '' when no sgRNA data was supplied."""
    if not has_sgrna:
        return ""
    info = consistency_lookup.get(gene)
    if not info:
        return " (no sgRNA-level rows found for this gene in the count table)"
    if info["single_guide_driven"]:
        return f" — sgRNA consistency: {info['n_concordant_guides']}/{info['n_guides']} guides concordant (single-guide-driven, weaker evidence)"
    return f" — sgRNA consistency: {info['n_concordant_guides']}/{info['n_guides']} guides concordant"


def compute_caveats(n_pos, n_neg, has_sgrna, consistency_lookup, pos_hits, neg_hits, gene_key, qc_stats):
    """Check this run's actual data against references/pitfalls.md and return
    only the caveats that are triggered (never a static list of all pitfalls).
    Each caveat is a (title, detail) pair: title feeds the short conclusion
    clause, detail is the full bullet text."""
    caveats = []

    if has_sgrna:
        driven = []
        for r in pos_hits + neg_hits:
            gene = r[gene_key]
            info = consistency_lookup.get(gene)
            if info and info["single_guide_driven"]:
                driven.append(f"`{gene}` ({info['n_concordant_guides']}/{info['n_guides']} guides concordant)")
        if driven:
            caveats.append((
                "single-guide-driven hit(s)",
                "**Single-guide-driven hit(s)** (pitfalls.md #6): " + "; ".join(driven) + " reach gene-level "
                "significance with essentially one sgRNA carrying the signal while the remaining guides show "
                "near-zero effect. Treat these as weaker evidence than hits with concordant multi-guide support; "
                "inspect `figures/sgrna_consistency.png` and rule out off-target or guide-efficiency explanations "
                "before treating them as high-confidence.",
            ))
    elif n_pos or n_neg:
        caveats.append((
            "no sgRNA-level data supplied",
            "**No sgRNA-level data supplied** (pitfalls.md #6): only the gene-level summary was provided, so this "
            "run cannot check whether any top hit is driven by a single sgRNA out of several. Re-run with "
            "`--sgrna-counts` to generate the guide-consistency check before treating any hit as high-confidence.",
        ))

    if has_sgrna and qc_stats:
        high_gini = [f"{s} (Gini={g:.2f})" for s, g in qc_stats["gini_by_sample"].items() if g > 0.2]
        if high_gini:
            caveats.append((
                "uneven library representation",
                "**Uneven library representation** (pitfalls.md #1): Gini index exceeds 0.2 for " +
                ", ".join(high_gini) + ". If this includes a T0/plasmid-pool sample, treat the whole screen as "
                "potentially underpowered rather than trusting downstream depletion calls; verify cells-per-sgRNA "
                "coverage at transduction before finalizing hit calls.",
            ))
        mgpg = qc_stats.get("mean_guides_per_gene")
        if mgpg is not None and mgpg < 3:
            caveats.append((
                "low guide multiplicity",
                f"**Low guide multiplicity** (pitfalls.md #6): the count table averages {mgpg:.1f} sgRNAs per "
                "gene. Fewer independent guides per gene weakens confidence in any single gene-level call; require "
                "concordant direction in at least half the guides before treating a gene as high-confidence.",
            ))
        if qc_stats.get("n_ntc_guides", 0) == 0:
            caveats.append((
                "no non-targeting controls",
                "**No non-targeting controls detected** (pitfalls.md #3): the count table has no guides labelled "
                "`NTC`. Without a null/reference distribution, any depletion/enrichment threshold is arbitrary; "
                "label results as an uncalibrated rank list rather than validated essential/resistance gene calls.",
            ))
        r = qc_stats.get("replicate_pearson_r_control")
        if r is not None and r < 0.8:
            caveats.append((
                "weak control replicate correlation",
                f"**Weak control replicate correlation** (Pearson r={r:.2f}): control replicates correlate more "
                "weakly than the r>=0.8 rule of thumb for pooled-screen QC, which reduces confidence in any "
                "downstream fold-change estimate. Inspect `figures/qc_distribution.png` before trusting hit calls.",
            ))

    if n_neg:
        caveats.append((
            "essentiality calls are uncalibrated",
            "**Essentiality calls are uncalibrated** (pitfalls.md #3): this run has no essential/nonessential "
            "reference gene set or non-targeting-control baseline as input, so negative-selection hits are a "
            "rank/FDR-based list, not a validated essential-gene call. Cross-check against a reference (e.g. "
            "Hart/DepMap common-essential and nonessential lists) before using essentiality language.",
        ))

    if n_pos or n_neg:
        caveats.append((
            "orthogonal validation required",
            "**Orthogonal validation required** (pitfalls.md #8): pooled-screen ranking alone cannot rule out "
            "off-target guides, copy-number artifacts, or batch effects. Do not call a hit \"validated\" or "
            "\"essential\"/a confirmed drug target from this output alone; arrayed knockout/knockdown + rescue, a "
            "competition assay, or an independent orthogonal screen is required before mechanistic claims.",
        ))

    return caveats


def build_report(ctx):
    """Assemble a results-and-conclusions report.md: named top hits (not just
    counts) with sgRNA-consistency notes, a QC read on the sgRNA data when
    supplied, only the pitfalls.md caveats actually triggered by this run's
    data, and a plain-language conclusion — not a file/parameter listing."""
    gene_key, score_key, fdr_key = ctx["gene_key"], ctx["score_key"], ctx["fdr_key"]
    lines = ["# CRISPR screen quicklook report", ""]
    lines.append(
        "Exploratory hit ranking and QC figures generated from the provided gene summary"
        + (" and sgRNA-level count table" if ctx["has_sgrna"] else "")
        + " (`quickplot_crispr.py`). This is the Python quicklook route: it visualizes an existing MAGeCK-style "
        "gene summary rather than re-running MAGeCK/BAGEL2/drugZ statistics. Use the formal R/MAGeCK route "
        "(see `references/cns-method-guide.md`) when raw counts and a design matrix are available and "
        "confirmatory statistics are required."
    )
    lines.append("")

    lines.append("## Input summary")
    lines.append("")
    lines.append(f"- Genes in gene summary: {ctx['n_genes']}")
    lines.append(f"- Thresholds: FDR < {ctx['fdr_threshold']}, abs(logFC/beta) >= {ctx['lfc_threshold']}")
    if ctx["has_sgrna"]:
        qs = ctx["qc_stats"]
        mgpg = f"{qs['mean_guides_per_gene']:.1f}" if qs.get("mean_guides_per_gene") is not None else "n/a"
        lines.append(
            f"- sgRNA-level count table: {qs['n_sgrnas']} guides across {qs['n_genes_with_guides']} genes "
            f"(mean {mgpg} guides/gene) plus {qs['n_ntc_guides']} non-targeting guide(s); "
            f"control columns: {', '.join(qs['control_columns']) or 'none detected'}; "
            f"treatment columns: {', '.join(qs['treatment_columns']) or 'none detected'}"
        )
    else:
        lines.append("- sgRNA-level count table: not provided (gene-level summary only, no guide-consistency check)")
    lines.append("")

    lines.append("## Key findings")
    lines.append("")
    lines.append(f"- Positive selection (enriched) hits: {ctx['n_pos']}")
    lines.append(f"- Negative selection (depleted) hits: {ctx['n_neg']}")
    lines.append(f"- Not significant: {ctx['n_ns']}")
    lines.append("")
    if ctx["pos_hits"]:
        lines.append("Top positive selection hits:")
        for r in ctx["pos_hits"]:
            gene = r[gene_key]
            note = get_consistency_note(gene, ctx["has_sgrna"], ctx["consistency_lookup"])
            lines.append(f"  1. `{gene}`: logFC = {float(r[score_key]):.2f}, FDR = {float(r[fdr_key]):.2e}{note}")
    else:
        lines.append("- No positive selection hits above threshold.")
    lines.append("")
    if ctx["neg_hits"]:
        lines.append("Top negative selection hits:")
        for r in ctx["neg_hits"]:
            gene = r[gene_key]
            note = get_consistency_note(gene, ctx["has_sgrna"], ctx["consistency_lookup"])
            lines.append(f"  1. `{gene}`: logFC = {float(r[score_key]):.2f}, FDR = {float(r[fdr_key]):.2e}{note}")
    else:
        lines.append("- No negative selection hits above threshold.")
    lines.append("")

    lines.append("## Caveats (automatically checked against references/pitfalls.md)")
    lines.append("")
    if ctx["caveats"]:
        for _title, detail in ctx["caveats"]:
            lines.append(f"- {detail}")
    else:
        lines.append(
            "- No data-triggered caveats from the automatic checks (no single-guide-driven top hit, library "
            "representation/replicate QC within range, non-targeting controls present, no significant hits). "
            "This does not replace manual review of `references/pitfalls.md`."
        )
    lines.append("")

    lines.append("## Conclusion")
    lines.append("")
    robustness = (
        "hypothesis-generating: this quicklook route ranks hits from the supplied summary and should be confirmed "
        "with a formal MAGeCK/BAGEL2/drugZ run and orthogonal validation before any biological or mechanistic claim"
    )
    caveat_clause = (
        ", with the following flagged concerns: " + ", ".join(title for title, _detail in ctx["caveats"])
        if ctx["caveats"]
        else ""
    )
    lines.append(
        f"{ctx['n_pos']} gene(s) showed positive selection and {ctx['n_neg']} showed negative selection out of "
        f"{ctx['n_genes']} tested, at FDR<{ctx['fdr_threshold']} and abs(logFC)>={ctx['lfc_threshold']}. "
        f"Results should be treated as {robustness}{caveat_clause}."
    )
    lines.append("")

    lines.append("## Files")
    lines.append("")
    fig_list = "`crispr_volcano`, `gene_rank_plot`"
    if ctx["top_hits_bar_drawn"]:
        fig_list += ", `top_hits_barplot`"
    if ctx["has_sgrna"]:
        fig_list += ", `qc_distribution`"
        if ctx["consistency_generated"]:
            fig_list += ", `sgrna_consistency`"
    lines.append(f"Figures: {fig_list}, each exported as PNG and PDF. Full gene table: `results/gene_summary.tsv`. "
                  f"Top hits: `results/top_hits.tsv`.")
    if ctx["consistency_generated"]:
        lines.append("Per-gene guide consistency: `results/sgrna_summary.tsv`.")

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="CRISPR screen quick plots and result tables from a gene summary (and optional sgRNA counts).")
    parser.add_argument("--gene-summary", required=True)
    parser.add_argument("--sgrna-counts", default=None, help="optional sgRNA-level raw count table for QC/consistency plots")
    parser.add_argument("--out", default="output")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--lfc-threshold", type=float, default=1.0)
    parser.add_argument("--top-n-labels", type=int, default=8, help="max labeled hit genes per direction on volcano/rank/bar plots")
    parser.add_argument("--top-n-consistency", type=int, default=6, help="max hit genes per direction on the sgRNA consistency plot")
    args = parser.parse_args()
    apply_style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    rows = read(args.gene_summary)
    if not rows:
        raise SystemExit(f"{args.gene_summary} has no data rows")
    fieldnames = list(rows[0].keys())
    score_key = pick_column(fieldnames, SCORE_CANDIDATES, "score (logFC/beta)")
    fdr_key = pick_column(fieldnames, FDR_CANDIDATES, "FDR/padj")
    gene_key = pick_column(fieldnames, GENE_CANDIDATES, "gene identifier")

    n_pos, n_neg, n_ns = plot_volcano(rows, score_key, fdr_key, gene_key, args.lfc_threshold, args.fdr_threshold, args.top_n_labels, out)
    plot_gene_rank(rows, score_key, fdr_key, gene_key, args.lfc_threshold, args.fdr_threshold, args.top_n_labels, out)
    top_hits_bar_drawn = plot_top_hits_bar(rows, score_key, fdr_key, gene_key, args.lfc_threshold, args.fdr_threshold, args.top_n_labels, out) is not None

    qc_stats = None
    consistency_rows = None
    if args.sgrna_counts:
        sg_rows, sgrna_col, gene_col, sample_cols = read_sgrna_counts(args.sgrna_counts)
        qc_stats = plot_qc(sg_rows, sgrna_col, gene_col, sample_cols, out)
        consistency_rows = plot_sgrna_consistency(
            sg_rows, sgrna_col, gene_col, sample_cols, rows, gene_key, score_key, fdr_key,
            args.lfc_threshold, args.fdr_threshold, args.top_n_consistency, out,
        )
        if consistency_rows:
            with (out / "results" / "sgrna_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(consistency_rows[0].keys()), delimiter="\t")
                writer.writeheader()
                writer.writerows(consistency_rows)

    # results/gene_summary.tsv: real tab-separated conversion (previous version copied the
    # CSV bytes verbatim with a .tsv extension, which is not actually tab-separated).
    with (out / "results" / "gene_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    top = sorted(rows, key=lambda r: float(r[fdr_key]))[:20]
    with (out / "results" / "top_hits.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(top)

    summary_lines = [
        "item\tvalue",
        f"genes\t{len(rows)}",
        f"positive_selection_hits\t{n_pos}",
        f"negative_selection_hits\t{n_neg}",
        f"not_significant\t{n_ns}",
        f"fdr_threshold\t{args.fdr_threshold}",
        f"lfc_threshold\t{args.lfc_threshold}",
        f"top_hits_reported\t{len(top)}",
    ]
    if qc_stats:
        summary_lines.append(f"n_sgrnas\t{qc_stats['n_sgrnas']}")
        summary_lines.append(f"n_ntc_guides\t{qc_stats['n_ntc_guides']}")
        summary_lines.append(f"mean_guides_per_gene\t{qc_stats['mean_guides_per_gene']}")
        if qc_stats["replicate_pearson_r_control"] is not None:
            summary_lines.append(f"replicate_pearson_r_control\t{qc_stats['replicate_pearson_r_control']:.4f}")
    (out / "summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    qc_payload = {
        "genes": len(rows),
        "positive_selection_hits": n_pos,
        "negative_selection_hits": n_neg,
        "not_significant": n_ns,
        "fdr_threshold": args.fdr_threshold,
        "lfc_threshold": args.lfc_threshold,
    }
    if qc_stats:
        qc_payload["sgrna_counts"] = qc_stats
    (out / "run" / "qc.json").write_text(json.dumps(qc_payload, indent=2), encoding="utf-8")

    top_n_report = 5
    pos_hits_report = top_hit_rows(rows, score_key, fdr_key, gene_key, args.lfc_threshold, args.fdr_threshold, top_n_report, "pos")
    neg_hits_report = top_hit_rows(rows, score_key, fdr_key, gene_key, args.lfc_threshold, args.fdr_threshold, top_n_report, "neg")
    consistency_lookup = {r["gene"]: r for r in consistency_rows} if consistency_rows else {}
    caveats = compute_caveats(
        n_pos, n_neg, bool(args.sgrna_counts), consistency_lookup, pos_hits_report, neg_hits_report, gene_key, qc_stats,
    )
    report_ctx = {
        "gene_key": gene_key,
        "score_key": score_key,
        "fdr_key": fdr_key,
        "n_genes": len(rows),
        "fdr_threshold": args.fdr_threshold,
        "lfc_threshold": args.lfc_threshold,
        "has_sgrna": bool(args.sgrna_counts),
        "qc_stats": qc_stats,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_ns": n_ns,
        "pos_hits": pos_hits_report,
        "neg_hits": neg_hits_report,
        "consistency_lookup": consistency_lookup,
        "caveats": caveats,
        "top_hits_bar_drawn": top_hits_bar_drawn,
        "consistency_generated": bool(consistency_rows),
    }
    (out / "report.md").write_text(build_report(report_ctx), encoding="utf-8")
    print(f"wrote CRISPR outputs to {out}")


if __name__ == "__main__":
    main()
