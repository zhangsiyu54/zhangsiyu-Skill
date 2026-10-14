#!/usr/bin/env python3
import argparse
import csv
import json
import math
import re
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit("Install matplotlib: python3 -m pip install matplotlib") from exc


TEXT = "#173B5F"
GRID = "#DCEAF5"
OPEN = "#E45B5B"
CLOSE = "#4EA3F1"
NS = "#B8C2CC"

PEAK_ID_RE = re.compile(r"^[^:]+:([0-9]+)-([0-9]+)$")


def read_csv(path):
    text = Path(path).read_text(encoding="utf-8-sig")
    rows = list(csv.DictReader(text.splitlines()))
    return rows


def require_table(rows, required_cols, table_name, path):
    """Fail with a clear, actionable message instead of an opaque
    IndexError/KeyError on empty or malformed tables."""
    if not rows:
        raise SystemExit(
            f"{table_name} table '{path}' has no data rows. Provide at least one row "
            f"with columns {required_cols}."
        )
    available = list(rows[0].keys())
    missing = [c for c in required_cols if c not in available]
    if missing:
        raise SystemExit(
            f"{table_name} table '{path}' is missing required column(s) {missing}. "
            f"Found columns: {available}. See references/input-contract.md."
        )


def resolve_column(rows, options, table_name, path):
    available = list(rows[0].keys())
    for opt in options:
        if opt in available:
            return opt
    raise SystemExit(
        f"{table_name} table '{path}' must contain one of {options}. "
        f"Found columns: {available}. See references/input-contract.md."
    )


def to_float(value, field, row_id, table_name, path):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise SystemExit(
            f"{table_name} table '{path}': could not parse {field}='{value}' as a number "
            f"for row '{row_id}'."
        )


def classify(log2fc, padj, padj_threshold, log2fc_threshold):
    if padj < padj_threshold and log2fc >= log2fc_threshold:
        return "opening"
    if padj < padj_threshold and log2fc <= -log2fc_threshold:
        return "closing"
    return "ns"


def peak_width(peak_id):
    m = PEAK_ID_RE.match(peak_id or "")
    if not m:
        return None
    start, end = int(m.group(1)), int(m.group(2))
    width = end - start
    return width if width > 0 else None


def spearman(xs, ys):
    """Spearman rank correlation with tie-averaging, no scipy dependency."""
    n = len(xs)
    if n < 2:
        return None

    def rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[order[k]] = avg_rank
            i = j + 1
        return ranks

    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    if vx == 0 or vy == 0:
        return None
    return cov / math.sqrt(vx * vy)


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 16, "axes.titlesize": 24, "axes.labelsize": 21, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15, "axes.edgecolor": "#B9D7EE", "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT, "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42})


def fmt_num(value, digits=3):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}g}"


def build_report(peaks, motifs_full, top_motifs, motif_sort_col, motif_label_col,
                  motif_has_matched_bg, all_samples, width_corr, width_note,
                  padj_threshold, log2fc_threshold):
    n_tested = len(peaks)
    n_up = sum(1 for p in peaks if p["change"] == "opening")
    n_down = sum(1 for p in peaks if p["change"] == "closing")
    n_ns = n_tested - n_up - n_down
    n_sig = n_up + n_down
    sig_frac = n_sig / n_tested if n_tested else 0.0

    top_peaks = sorted(peaks, key=lambda p: p["padj"])[:5]
    top_motif_rows = top_motifs[:5]

    lines = []
    lines.append("# ATAC quicklook: results and conclusions")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append(f"- Peaks tested: {n_tested}")
    lines.append(f"- Samples in peak count matrix: {len(all_samples)} ({', '.join(all_samples)})")
    lines.append(f"- Motifs tested: {len(motifs_full)}")
    lines.append(
        f"- Significant differential peaks (padj < {padj_threshold}, |log2FC| >= {log2fc_threshold}): "
        f"{n_sig} total -- {n_up} opening (more accessible), {n_down} closing (less accessible), "
        f"{n_ns} not significant."
    )
    lines.append("")

    lines.append("## Top differential peaks")
    lines.append("")
    if n_tested == 0:
        lines.append("No peaks were tested.")
    else:
        for i, p in enumerate(top_peaks, 1):
            direction = {"opening": "opening", "closing": "closing", "ns": "not significant"}[p["change"]]
            lines.append(
                f"{i}. `{p['peak']}` -- log2FC = {fmt_num(p['log2fc'])}, padj = {fmt_num(p['padj'], 3)} ({direction})"
            )
    lines.append("")

    lines.append("## Top enriched motifs")
    lines.append("")
    if not top_motif_rows:
        lines.append("No motifs were tested.")
    else:
        for i, m in enumerate(top_motif_rows, 1):
            extra = f", enrichment = {m['enrichment']}" if m.get("enrichment") not in (None, "") else ""
            lines.append(
                f"{i}. **{m['label']}** -- {motif_sort_col} = {fmt_num(m['sort_value'], 3)}{extra}"
            )
    lines.append("")

    caveats = []
    if width_corr is not None and abs(width_corr) > 0.3:
        caveats.append(
            f"**Peak-width bias (pitfalls.md #1)**: peak width correlates with log2FC "
            f"(Spearman r = {width_corr:.3f}, |r| > 0.3). Differential peaks may skew toward "
            f"wider peaks rather than truly more accessible regions -- consider fixed-width "
            f"summit-centered peaks or a peak-width covariate before trusting this list."
        )
    elif width_corr is not None:
        caveats.append(
            f"Peak-width vs. log2FC check: no strong trend detected (Spearman r = {width_corr:.3f}, "
            f"|r| <= 0.3), so wide-peak bias (pitfalls.md #1) is unlikely to be driving the results."
        )
    else:
        caveats.append(width_note)

    if n_tested > 0 and n_sig == 0:
        caveats.append(
            "**No significant differential peaks were detected.** This may reflect a true null "
            "result, an underpowered design (too few replicates or peaks), or overly strict "
            "thresholds -- do not report a positive differential accessibility finding from this run."
        )
    elif n_tested >= 20 and 0 < n_sig < 3:
        caveats.append(
            f"**Possible underpowered result**: only {n_sig}/{n_tested} peaks reached significance. "
            f"Treat the differential peak list as preliminary and check replicate number/sequencing depth."
        )
    if n_tested >= 20 and sig_frac > 0.5:
        caveats.append(
            f"**Implausibly large differential peak fraction (pitfalls.md #1)**: {n_sig}/{n_tested} "
            f"peaks ({sig_frac:.0%}) are called significant. This can indicate a normalization or "
            f"majority-shift problem (e.g., global accessibility shift, size-factor distortion) rather "
            f"than genuine focal differential accessibility -- inspect normalization and diagnostic plots "
            f"before trusting this list."
        )

    if motifs_full and not motif_has_matched_bg:
        caveats.append(
            "**Motif background not confirmed (pitfalls.md #4)**: this table does not state that a "
            "GC/length-matched background peak set was used. If the motif test was run against a "
            "generic genomic background, GC-rich motifs may appear enriched purely from base "
            "composition -- confirm the background strategy before interpreting the ranked list."
        )

    lines.append("## Caveats")
    lines.append("")
    if caveats:
        for c in caveats:
            lines.append(f"- {c}")
    else:
        lines.append("- No automatic caveats were triggered for this run.")
    lines.append("")

    lines.append("## Conclusion")
    lines.append("")
    if n_tested == 0:
        conclusion = "No peaks were available to test; no conclusion can be drawn from this run."
    elif n_sig == 0:
        conclusion = (
            f"Across {n_tested} tested peaks, none met the significance thresholds "
            f"(padj < {padj_threshold}, |log2FC| >= {log2fc_threshold}), so this run does not support a differential "
            f"accessibility conclusion between the compared conditions. "
        )
    else:
        lead_peak = top_peaks[0]
        conclusion = (
            f"Of {n_tested} tested peaks, {n_sig} ({sig_frac:.0%}) show significant differential "
            f"accessibility ({n_up} more accessible, {n_down} less accessible) between the compared "
            f"conditions, led by `{lead_peak['peak']}` (log2FC = {fmt_num(lead_peak['log2fc'])}, "
            f"padj = {fmt_num(lead_peak['padj'], 3)}). "
        )
    if top_motif_rows:
        top_tf_names = ", ".join(m["label"] for m in top_motif_rows[:3])
        conclusion += (
            f"Motif enrichment highlights {top_tf_names} as the top candidate transcription-factor "
            f"motifs overlapping the tested peak set; this indicates motif accessibility, not "
            f"confirmed TF binding or activity (pitfalls.md #4, #14), and any nearest-gene or "
            f"promoter-accessibility interpretation should be reported as a hypothesis, not proven "
            f"expression change (pitfalls.md #3), unless corroborated by independent data."
        )
    else:
        conclusion += "No motif enrichment results were available to interpret."
    lines.append(conclusion)
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="ATAC quick plots from differential peak and motif tables.")
    parser.add_argument("--differential-peaks", required=True)
    parser.add_argument("--motifs", required=True)
    parser.add_argument("--peak-counts", required=True)
    parser.add_argument("--out", default="output")
    parser.add_argument(
        "--padj-threshold", type=float, default=0.05,
        help="Adjusted p-value cutoff for calling a peak significant (default: 0.05).",
    )
    parser.add_argument(
        "--log2fc-threshold", type=float, default=1.0,
        help="Minimum |log2 fold change| for calling a peak significant (default: 1.0).",
    )
    parser.add_argument(
        "--top-n-heatmap", type=int, default=30,
        help="Number of peaks (by input table order) shown in the peak accessibility heatmap (default: 30).",
    )
    parser.add_argument(
        "--top-n-motifs", type=int, default=12,
        help="Number of top motifs shown in the motif enrichment barplot (default: 12).",
    )
    args = parser.parse_args()
    style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    # ---- Load and validate differential peaks table -----------------------
    peaks_raw = read_csv(args.differential_peaks)
    require_table(peaks_raw, ["peak", "log2FC", "padj"], "Differential-peaks", args.differential_peaks)
    peaks = []
    for r in peaks_raw:
        log2fc = to_float(r["log2FC"], "log2FC", r.get("peak", "?"), "Differential-peaks", args.differential_peaks)
        padj = to_float(r["padj"], "padj", r.get("peak", "?"), "Differential-peaks", args.differential_peaks)
        padj_clamped = min(max(padj, 0.0), 1.0)
        peaks.append({
            "peak": r["peak"],
            "log2fc": log2fc,
            "padj": padj_clamped,
            "change": classify(log2fc, padj_clamped, args.padj_threshold, args.log2fc_threshold),
        })

    xs = [p["log2fc"] for p in peaks]
    ys = [-math.log10(max(p["padj"], 1e-12)) for p in peaks]
    colors = [OPEN if p["change"] == "opening" else CLOSE if p["change"] == "closing" else NS for p in peaks]
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    ax.scatter(xs, ys, c=colors, s=34, alpha=0.84, linewidths=0)
    ax.axvline(-args.log2fc_threshold, color="#A7C6DD", linestyle="--")
    ax.axvline(args.log2fc_threshold, color="#A7C6DD", linestyle="--")
    ax.axhline(-math.log10(args.padj_threshold), color="#A7C6DD", linestyle="--")
    ax.set_title(
        f"Differential accessibility (padj<{args.padj_threshold}, |log2FC|>={args.log2fc_threshold})",
        loc="center", pad=18, fontweight="bold", color=TEXT,
    )
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p value")
    ax.grid(True, color=GRID)
    save(fig, out / "figures" / "differential_peak_volcano")

    # ---- Load and validate motif table -------------------------------------
    motifs_full = read_csv(args.motifs)
    top_motifs = []
    motif_sort_col = None
    motif_label_col = None
    motif_has_matched_bg = False
    if motifs_full:
        motif_label_col = "tf" if "tf" in motifs_full[0] else resolve_column(motifs_full, ["tf", "motif"], "Motif", args.motifs)
        motif_sort_col = resolve_column(motifs_full, ["padj", "pvalue"], "Motif", args.motifs)
        motif_has_matched_bg = any(
            str(r.get("background", r.get("background_matched", ""))).strip().lower() in ("1", "true", "yes", "matched")
            for r in motifs_full
        )
        for r in motifs_full:
            sort_value = to_float(r[motif_sort_col], motif_sort_col, r.get(motif_label_col, "?"), "Motif", args.motifs)
            top_motifs.append({
                "label": r[motif_label_col],
                "sort_value": min(max(sort_value, 0.0), 1.0),
                "enrichment": r.get("enrichment"),
            })
        top_motifs.sort(key=lambda m: m["sort_value"])
        top_motifs = top_motifs[: args.top_n_motifs]

        fig, ax = plt.subplots(figsize=(8.2, max(4.8, len(top_motifs) * 0.42 + 2)))
        labels = [m["label"] for m in reversed(top_motifs)]
        vals = [-math.log10(max(m["sort_value"], 1e-12)) for m in reversed(top_motifs)]
        ax.barh(labels, vals, color="#4EA3F1", alpha=0.86)
        ax.set_title("Motif enrichment", loc="center", pad=18, fontweight="bold", color=TEXT)
        ax.set_xlabel(f"-log10 adjusted {motif_sort_col}" if motif_sort_col == "padj" else f"-log10 {motif_sort_col}")
        ax.grid(True, axis="x", color=GRID)
        save(fig, out / "figures" / "motif_enrichment_barplot")
    else:
        print("warning: motif table has no rows; skipping motif_enrichment_barplot.png/.pdf")

    # ---- Load and validate peak count matrix -------------------------------
    peak_counts_full = read_csv(args.peak_counts)
    require_table(peak_counts_full, ["peak"], "Peak-counts", args.peak_counts)
    all_samples = [k for k in peak_counts_full[0].keys() if k != "peak"]
    if not all_samples:
        raise SystemExit(
            f"Peak-counts table '{args.peak_counts}' has a 'peak' column but no sample columns. "
            "Add at least one sample count column (see references/input-contract.md)."
        )
    lib_depth = {s: 0.0 for s in all_samples}
    peaks_detected = {s: 0 for s in all_samples}
    for r in peak_counts_full:
        for s in all_samples:
            v = to_float(r[s], s, r.get("peak", "?"), "Peak-counts", args.peak_counts)
            lib_depth[s] += v
            if v > 0:
                peaks_detected[s] += 1

    fig, axes = plt.subplots(1, 2, figsize=(max(9.6, len(all_samples) * 0.95), 5.8))
    axes[0].bar(all_samples, [lib_depth[s] for s in all_samples], color=CLOSE, alpha=0.86)
    axes[0].set_title("Library depth", loc="center", pad=14, fontweight="bold", color=TEXT)
    axes[0].set_ylabel("Total peak counts")
    axes[1].bar(all_samples, [peaks_detected[s] for s in all_samples], color=OPEN, alpha=0.86)
    axes[1].set_title("Peaks detected", loc="center", pad=14, fontweight="bold", color=TEXT)
    axes[1].set_ylabel("Peaks with count > 0")
    for ax_ in axes:
        ax_.grid(True, axis="y", color=GRID)
        ax_.tick_params(axis="x", rotation=35)
        for label in ax_.get_xticklabels():
            label.set_ha("right")
    fig.suptitle("Sample QC overview", fontweight="bold", color=TEXT, fontsize=24, y=1.08)
    fig.tight_layout()
    save(fig, out / "figures" / "qc_overview")

    counts = peak_counts_full[: min(args.top_n_heatmap, len(peak_counts_full))]
    samples = [k for k in counts[0].keys() if k != "peak"]
    mat = []
    names = []
    for r in counts:
        vals = [math.log2(to_float(r[s], s, r.get("peak", "?"), "Peak-counts", args.peak_counts) + 1) for s in samples]
        mean = sum(vals) / len(vals)
        sd = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5 or 1
        mat.append([(v - mean) / sd for v in vals])
        names.append(r["peak"])
    fig, ax = plt.subplots(figsize=(max(7.2, len(samples) * 0.55 + 3), max(6, len(mat) * 0.18 + 2.8)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
    ax.set_title("Peak accessibility heatmap", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=35, ha="right")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Row-scaled accessibility", fontsize=16, labelpad=10)
    save(fig, out / "figures" / "peak_heatmap")

    # ---- Peak-width vs. log2FC diagnostic (pitfalls.md #1) -----------------
    widths = [peak_width(p["peak"]) for p in peaks]
    parseable = [(w, p["log2fc"]) for w, p in zip(widths, peaks) if w is not None]
    width_corr = None
    if len(parseable) < 10:
        width_note = (
            "Peak-width vs. log2FC check skipped: fewer than 10 peak IDs were parseable as "
            "`chr:start-end`, so peak width could not be computed."
        )
    else:
        width_corr = spearman([w for w, _ in parseable], [fc for _, fc in parseable])
        if width_corr is None:
            width_note = "Peak-width vs. log2FC check skipped: peak widths or log2FC values had no variance."

    # ---- Write tables and JSON ----------------------------------------------
    (out / "results" / "differential_peaks.tsv").write_text(
        "\n".join(["\t".join(peaks_raw[0].keys())] + ["\t".join(r.values()) for r in peaks_raw]) + "\n",
        encoding="utf-8",
    )
    (out / "results" / "motif_enrichment.tsv").write_text(Path(args.motifs).read_text(encoding="utf-8"), encoding="utf-8")

    n_up = sum(1 for p in peaks if p["change"] == "opening")
    n_down = sum(1 for p in peaks if p["change"] == "closing")
    (out / "summary.tsv").write_text(
        "item\tvalue\n"
        f"peaks_tested\t{len(peaks)}\n"
        f"opening_peaks\t{n_up}\n"
        f"closing_peaks\t{n_down}\n"
        f"samples\t{len(all_samples)}\n"
        f"motifs_tested\t{len(motifs_full)}\n",
        encoding="utf-8",
    )
    (out / "run" / "qc.json").write_text(json.dumps({
        "peaks": len(peaks),
        "opening_peaks": n_up,
        "closing_peaks": n_down,
        "motifs": len(motifs_full),
        "samples": len(all_samples),
        "library_depth": {s: lib_depth[s] for s in all_samples},
        "peaks_detected": {s: peaks_detected[s] for s in all_samples},
        "peak_width_log2fc_spearman_r": width_corr,
        "padj_threshold": args.padj_threshold,
        "log2fc_threshold": args.log2fc_threshold,
        "top_n_heatmap": args.top_n_heatmap,
        "top_n_motifs": args.top_n_motifs,
    }, indent=2), encoding="utf-8")

    report = build_report(
        peaks=peaks,
        motifs_full=motifs_full,
        top_motifs=top_motifs,
        motif_sort_col=motif_sort_col or "padj",
        motif_label_col=motif_label_col,
        motif_has_matched_bg=motif_has_matched_bg,
        all_samples=all_samples,
        width_corr=width_corr,
        width_note=width_note,
        padj_threshold=args.padj_threshold,
        log2fc_threshold=args.log2fc_threshold,
    )
    (out / "report.md").write_text(report, encoding="utf-8")
    print(f"wrote ATAC outputs to {out}")


if __name__ == "__main__":
    main()
