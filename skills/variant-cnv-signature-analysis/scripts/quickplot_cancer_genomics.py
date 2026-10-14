#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit("Install matplotlib and numpy: python3 -m pip install matplotlib numpy") from exc

TEXT = "#173B5F"
GRID = "#DCEAF5"
RED = "#E45B5B"
BLUE = "#4EA3F1"
PALETTE = [RED, BLUE, "#4CB984", "#F2B84B", "#9B7EDE"]

LOW_MUTATION_COUNT_THRESHOLD = 50  # references/pitfalls.md #2
SMALL_COHORT_THRESHOLD = 50  # references/pitfalls.md #8


def to_float(value, context):
    """Parse a numeric cell, failing with a message that names the offending
    row/column instead of a bare ValueError/TypeError traceback."""
    try:
        return float(value)
    except (TypeError, ValueError):
        raise SystemExit(f"Could not parse '{value}' as a number ({context}). Check the input file for non-numeric or missing values.")


def read(path):
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError:
        raise SystemExit(f"Input file not found: {path}")
    return list(csv.DictReader(text.splitlines()))


def require_columns(rows, required_cols, path, kind):
    if not rows:
        return
    missing = [c for c in required_cols if c not in rows[0]]
    if missing:
        raise SystemExit(
            f"{kind} table at {path} is missing required column(s): {missing}. "
            f"Expected at least: {required_cols}. Found columns: {list(rows[0].keys())}."
        )


def require_nonempty(rows, path, kind):
    if not rows:
        raise SystemExit(
            f"{kind} table at {path} has no data rows (file is empty or header-only). "
            f"Provide at least one row before running quickplot."
        )


def require_sample_columns(rows, id_col, path, kind):
    samples = [k for k in rows[0] if k != id_col]
    if not samples:
        raise SystemExit(
            f"{kind} table at {path} has a '{id_col}' column but no sample columns. "
            f"Add at least one sample column."
        )
    return samples


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


def order_oncoplot_style(genes, samples, mat):
    """Order genes by alteration frequency (desc) and samples by alteration
    burden (desc), matching the oncoplot convention of "genes by frequency,
    samples by event" from references/paper-practice.md and visual-style.md."""
    gene_freq = [sum(1 for v in row if v != 0) for row in mat]
    gene_order = sorted(range(len(genes)), key=lambda i: gene_freq[i], reverse=True)
    sample_burden = [sum(abs(mat[i][j]) for i in range(len(genes))) for j in range(len(samples))]
    sample_order = sorted(range(len(samples)), key=lambda j: sample_burden[j], reverse=True)
    genes_sorted = [genes[i] for i in gene_order]
    samples_sorted = [samples[j] for j in sample_order]
    mat_sorted = [[mat[i][j] for j in sample_order] for i in gene_order]
    return genes_sorted, samples_sorted, mat_sorted


def plot_mutation_frequency(mut, out, path="mutation-frequency table"):
    mut_sorted = sorted(mut, key=lambda r: to_float(r["frequency"], f"frequency for gene {r.get('gene', '?')} in {path}"))
    genes = [r["gene"] for r in mut_sorted]
    freqs = [to_float(r["frequency"], f"frequency for gene {r.get('gene', '?')} in {path}") * 100 for r in mut_sorted]
    fig, ax = plt.subplots(figsize=(8.4, max(5.8, len(mut_sorted) * 0.42 + 2.0)))
    bars = ax.barh(genes, freqs, color=RED, alpha=0.88)
    max_freq = max(freqs) if freqs else 0.0
    xlim_upper = max_freq * 1.2 if max_freq > 0 else 1.0
    for bar, r in zip(bars, mut_sorted):
        # report count alongside percentage (pitfalls.md #8: avoid false precision)
        ax.text(
            bar.get_width() + xlim_upper * 0.015, bar.get_y() + bar.get_height() / 2,
            f"n={r['count']}", va="center", ha="left", fontsize=12, color=TEXT,
        )
    if max_freq == 0 and mut_sorted:
        ax.text(0.5, 0.5, "No mutations observed (all frequencies = 0)", transform=ax.transAxes,
                 ha="center", va="center", fontsize=14, color=TEXT, style="italic")
    ax.set_title("Mutation frequency", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Frequency (%)", labelpad=12)
    ax.set_ylabel("Gene", labelpad=12)
    ax.set_xlim(0, xlim_upper)
    ax.grid(True, axis="x", color=GRID)
    save(fig, out / "figures" / "mutation_frequency_barplot")


def plot_cnv_heatmap(cnv, out, path="CNV matrix"):
    samples = require_sample_columns(cnv, "gene", path, "CNV")
    genes = [r["gene"] for r in cnv]
    mat = [[to_float(r[s], f"CNV value for gene {r.get('gene', '?')} sample {s} in {path}") for s in samples] for r in cnv]
    genes, samples, mat = order_oncoplot_style(genes, samples, mat)
    fig, ax = plt.subplots(figsize=(max(7.4, len(samples) * 0.5 + 3), max(6.2, len(genes) * 0.34 + 2.0)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_title("CNV heatmap", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Sample (ordered by CNV burden)", labelpad=12)
    ax.set_ylabel("Gene (ordered by alteration frequency)", labelpad=12)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=35, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Copy-number score", fontsize=16, labelpad=10)
    save(fig, out / "figures" / "cnv_heatmap")
    return genes, samples


def plot_signature_exposure(sig, out, path="signature exposure table"):
    samples = require_sample_columns(sig, "signature", path, "Signature exposure")
    bottom = [0] * len(samples)
    fig, ax = plt.subplots(figsize=(max(6.5, len(samples) * 0.6 + 3.5), 6.0))
    for i, row in enumerate(sig):
        vals = [to_float(row[s], f"exposure for signature {row.get('signature', '?')} sample {s} in {path}") for s in samples]
        ax.bar(samples, vals, bottom=bottom, label=row["signature"], color=PALETTE[i % len(PALETTE)], alpha=0.9)
        bottom = [bottom[j] + vals[j] for j in range(len(vals))]
    ax.set_title("Mutational signature exposure", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Sample", labelpad=12)
    ax.set_ylabel("Exposure", labelpad=12)
    ax.tick_params(axis="x", rotation=35)
    if len(samples) > 1:
        for label in ax.get_xticklabels():
            label.set_ha("right")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left", title="Signature")
    save(fig, out / "figures" / "signature_exposure")


def read_tmb(path):
    rows = read(path)
    require_nonempty(rows, path, "TMB")
    require_columns(rows, ["sample", "assay", "panel_size_mb", "mutation_count"], path, "TMB")
    for r in rows:
        r["mutation_count"] = to_float(r["mutation_count"], f"mutation_count for sample {r.get('sample', '?')} in {path}")
        r["panel_size_mb"] = to_float(r["panel_size_mb"], f"panel_size_mb for sample {r.get('sample', '?')} in {path}")
        r["assay"] = r.get("assay") or "assay"
        r["tmb"] = r["mutation_count"] / r["panel_size_mb"] if r["panel_size_mb"] else float("nan")
    return rows


def plot_tmb_summary(tmb_rows, out):
    assays = sorted({r["assay"] for r in tmb_rows})
    data = [[r["tmb"] for r in tmb_rows if r["assay"] == a] for a in assays]
    fig, ax = plt.subplots(figsize=(max(6.5, len(assays) * 1.8 + 3), 6.2))
    try:
        bp = ax.boxplot(data, tick_labels=assays, patch_artist=True, widths=0.45, showfliers=False)
    except TypeError:
        bp = ax.boxplot(data, labels=assays, patch_artist=True, widths=0.45, showfliers=False)
    for patch, color in zip(bp["boxes"], PALETTE):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
        patch.set_edgecolor(TEXT)
        patch.set_linewidth(1.3)
    for key in ("whiskers", "caps", "medians"):
        for artist in bp[key]:
            artist.set_color(TEXT)
            artist.set_linewidth(1.3)
    rng = np.random.default_rng(7)
    for i, a in enumerate(assays):
        vals = [r["tmb"] for r in tmb_rows if r["assay"] == a]
        xs = rng.normal(i + 1, 0.05, size=len(vals))
        ax.scatter(xs, vals, color=TEXT, s=26, alpha=0.75, zorder=3, edgecolor="white", linewidth=0.4)
    ax.set_title("Tumor mutation burden by assay", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Assay", labelpad=12)
    ax.set_ylabel("TMB (mutations / Mb)", labelpad=12)
    ax.grid(True, axis="y", color=GRID)
    # Denominator caveat as a separate footnote below the tick labels/xlabel
    # (not appended to the xlabel string) so it cannot collide with long
    # assay tick labels such as "targeted_panel" (visual-style.md: inspect
    # for overlapping labels before delivery).
    ax.annotate(
        "Denominator = callable panel/exome Mb, per assay (pitfalls.md #1)",
        xy=(0.5, 0), xycoords="axes fraction", xytext=(0, -95), textcoords="offset points",
        ha="center", va="top", fontsize=13, color=TEXT, style="italic",
    )
    save(fig, out / "figures" / "tmb_summary_boxplot")


def check_sample_id_overlap(cnv_samples, sig_samples, tmb_rows):
    """Cross-check sample IDs across tables (cns-method-guide.md QC section:
    'Harmonize sample IDs across mutation, CNV, signature, and clinical
    tables'). Returns a list of note strings; does not raise, since a
    mismatch may be intentional (e.g. signature run on a subset)."""
    notes = []
    cnv_set, sig_set = set(cnv_samples), set(sig_samples)
    if cnv_set and sig_set:
        overlap = cnv_set & sig_set
        if not overlap:
            notes.append(
                f"Sample ID mismatch: the CNV matrix ({len(cnv_set)} samples) and the signature exposure "
                f"table ({len(sig_set)} samples) share zero sample IDs. Verify these tables describe the "
                f"same cohort before comparing CNV and signature results per-sample."
            )
        elif overlap != cnv_set or overlap != sig_set:
            notes.append(
                f"Sample ID overlap is partial between the CNV matrix and signature exposure table "
                f"({len(overlap)} shared / {len(cnv_set)} CNV / {len(sig_set)} signature samples). "
                f"Confirm both tables represent the same cohort before drawing cross-table conclusions."
            )
    if tmb_rows is not None:
        tmb_set = {r["sample"] for r in tmb_rows}
        if cnv_set and tmb_set and not (tmb_set & cnv_set):
            notes.append(
                f"Sample ID mismatch: the TMB table ({len(tmb_set)} samples) shares zero sample IDs with "
                f"the CNV matrix ({len(cnv_set)} samples). Verify these tables describe the same cohort."
            )
    return notes


def build_report(mut, cnv, sig, tmb_rows, cnv_samples, sig_samples, id_notes, out):
    """Write output/report.md: a results-and-conclusions summary, not just a
    file listing. See references/deliverables.md and references/pitfalls.md."""
    lines = ["# Variant / CNV / Signature Analysis Report", ""]

    n_genes = len(mut)
    n_cnv_samples = len(cnv_samples)
    n_sig_samples = len(sig_samples)
    lines.append("## Cohort overview")
    lines.append("")
    lines.append(f"- Genes in mutation frequency table: {n_genes}")
    lines.append(f"- Samples in CNV matrix: {n_cnv_samples}")
    lines.append(f"- Samples in signature exposure table: {n_sig_samples}")
    if tmb_rows is not None:
        lines.append(f"- Samples in TMB table: {len(tmb_rows)}")
    lines.append("")

    # --- Mutation frequency summary ----------------------------------------
    lines.append("## Top mutated genes")
    lines.append("")
    mutated = [r for r in mut if to_float(r.get("frequency", 0) or 0, "frequency") > 0 or to_float(r.get("count", 0) or 0, "count") > 0]
    if not mut:
        lines.append("No mutation frequency data provided.")
    elif not mutated:
        lines.append(
            "No mutated genes: all genes in the mutation frequency table have count = 0 / frequency = 0%. "
            "Confirm this reflects a true hypomutated cohort rather than an upstream filtering or join error."
        )
    else:
        top = sorted(mutated, key=lambda r: to_float(r["frequency"], "frequency"), reverse=True)[:5]
        precision = 0 if n_cnv_samples and n_cnv_samples < SMALL_COHORT_THRESHOLD else 1
        for r in top:
            pct = to_float(r["frequency"], "frequency") * 100
            lines.append(f"- **{r['gene']}**: {r['count']} sample(s) mutated ({pct:.{precision}f}%)")
        if n_cnv_samples and n_cnv_samples < SMALL_COHORT_THRESHOLD:
            lines.append("")
            lines.append(
                f"  Note: cohort size is small (n={n_cnv_samples}); percentages are rounded to whole numbers "
                f"and should be read alongside the raw counts above, not compared to published cohorts with false precision (pitfalls.md #8)."
            )
    lines.append("")

    # --- CNV summary ---------------------------------------------------------
    lines.append("## CNV summary")
    lines.append("")
    if not cnv:
        lines.append("No CNV data provided.")
    else:
        burden = []
        for r in cnv:
            vals = [to_float(r[s], f"CNV value for gene {r.get('gene', '?')} sample {s}") for s in cnv_samples]
            amp = sum(1 for v in vals if v > 0)
            deletion = sum(1 for v in vals if v < 0)
            burden.append({"gene": r["gene"], "amp": amp, "del": deletion})
        top_amp = sorted([b for b in burden if b["amp"] > 0], key=lambda b: b["amp"], reverse=True)[:3]
        top_del = sorted([b for b in burden if b["del"] > 0], key=lambda b: b["del"], reverse=True)[:3]
        if not top_amp and not top_del:
            lines.append(
                "All genes are copy-number neutral across all samples in this matrix (no amplification or "
                "deletion calls). Confirm this reflects the true calling result rather than a placeholder/"
                "all-zero matrix or an upstream join error."
            )
        else:
            if top_amp:
                lines.append("Highest amplification burden: " + ", ".join(f"{b['gene']} ({b['amp']}/{n_cnv_samples} samples)" for b in top_amp))
            else:
                lines.append("Highest amplification burden: none (no gene shows amplification in any sample).")
            if top_del:
                lines.append("Highest deletion burden: " + ", ".join(f"{b['gene']} ({b['del']}/{n_cnv_samples} samples)" for b in top_del))
            else:
                lines.append("Highest deletion burden: none (no gene shows deletion in any sample).")
        lines.append("")
        lines.append(
            "Caveat: these calls are raw gene-level copy-number scores, not purity/ploidy-corrected integer "
            "copy number. Low-purity or genome-doubled samples can under-call true focal events with a fixed "
            "threshold (pitfalls.md #4); treat amplification/deletion status as approximate without a "
            "dedicated purity/ploidy-aware caller (ASCAT/ABSOLUTE/FACETS/PureCN/Sequenza)."
        )
    lines.append("")

    # --- Signature summary ----------------------------------------------------
    lines.append("## Mutational signature exposure summary")
    lines.append("")
    if not sig:
        lines.append("No signature exposure data provided.")
    else:
        dominant = {}
        for s in sig_samples:
            best_sig, best_val = None, 0.0
            for row in sig:
                v = to_float(row[s], f"exposure for signature {row.get('signature', '?')} sample {s}")
                if v > best_val:
                    best_val, best_sig = v, row["signature"]
            dominant[s] = best_sig
        undetermined = [s for s, sg in dominant.items() if sg is None]
        determined = {s: sg for s, sg in dominant.items() if sg is not None}
        if not determined:
            lines.append("No sample has a positive exposure value for any signature; dominant signature is undetermined for all samples.")
        elif n_sig_samples == 1:
            s = sig_samples[0]
            lines.append(f"Sample **{s}**: dominant signature is **{determined.get(s, 'undetermined')}**.")
        else:
            counts = Counter(determined.values())
            top_sig, top_n = counts.most_common(1)[0]
            lines.append(
                f"Cohort-wide dominant signature: **{top_sig}** is the top signature in {top_n}/{n_sig_samples} samples."
            )
            if len(counts) > 1:
                breakdown = ", ".join(f"{sg} in {c} sample(s)" for sg, c in counts.most_common())
                lines.append(f"Per-signature dominance breakdown: {breakdown}.")
        if undetermined:
            lines.append(f"Dominant signature undetermined (all exposures <= 0) for: {', '.join(undetermined)}.")
        if tmb_rows is not None:
            overlap = set(sig_samples) & {r["sample"] for r in tmb_rows}
            low = [r for r in tmb_rows if r["sample"] in overlap and r["mutation_count"] < LOW_MUTATION_COUNT_THRESHOLD]
            if low:
                lines.append("")
                lines.append(
                    f"Caveat: {len(low)} sample(s) with signature exposure have fewer than "
                    f"{LOW_MUTATION_COUNT_THRESHOLD} mutations per the TMB table "
                    f"({', '.join(r['sample'] for r in low)}) -- single-sample signature deconvolution is "
                    f"unstable below this range (pitfalls.md #2); treat their exposures as low-confidence "
                    f"or prefer cohort-level fitting for these samples."
                )
        else:
            lines.append("")
            lines.append(
                "Note: no per-sample mutation counts were supplied (no --tmb table), so signature reliability "
                "at low mutation count could not be checked automatically; confirm each sample has enough "
                "mutations (commonly >=50-100 SNVs) for stable single-sample fitting (pitfalls.md #2)."
            )
    lines.append("")

    # --- TMB summary -----------------------------------------------------------
    if tmb_rows is not None:
        lines.append("## Tumor mutation burden (TMB) summary")
        lines.append("")
        assays = sorted({r["assay"] for r in tmb_rows})
        for a in assays:
            vals = sorted(r["tmb"] for r in tmb_rows if r["assay"] == a)
            n = len(vals)
            median = vals[n // 2] if n % 2 == 1 else (vals[n // 2 - 1] + vals[n // 2]) / 2
            lines.append(f"- **{a}** (n={n}): median TMB = {median:.2f} mut/Mb")
        if len(assays) > 1:
            lines.append("")
            lines.append(
                "TMB was computed and is reported separately per assay (own callable-Mb denominator per "
                "assay). Do not pool these distributions or compare medians across assays without explicit "
                "cross-platform harmonization (pitfalls.md #1) -- the per-assay boxplot in "
                "`figures/tmb_summary_boxplot.png` keeps assays visually separated for this reason."
            )
        lines.append("")

    # --- Sample ID harmonization notes -----------------------------------------
    if id_notes:
        lines.append("## Data quality flags")
        lines.append("")
        for note in id_notes:
            lines.append(f"- {note}")
        lines.append("")

    # --- Conclusion --------------------------------------------------------
    lines.append("## Conclusion")
    lines.append("")
    conclusion = []
    if mutated:
        top1 = sorted(mutated, key=lambda r: to_float(r["frequency"], "frequency"), reverse=True)[0]
        conclusion.append(
            f"Across {n_cnv_samples or n_genes} sample(s), {top1['gene']} is the most frequently mutated gene "
            f"({top1['count']} sample(s))."
        )
    else:
        conclusion.append("No mutated genes were observed in the supplied mutation frequency table.")
    if cnv:
        if top_amp or top_del:
            focus = top_amp[0]["gene"] if top_amp else top_del[0]["gene"]
            conclusion.append(f"CNV burden is concentrated around {focus}, but calls are not purity/ploidy-corrected.")
        else:
            conclusion.append("No CNV alterations were detected in the supplied matrix.")
    if sig and determined:
        conclusion.append(f"The dominant mutational process across the cohort is consistent with {top_sig if n_sig_samples > 1 else list(determined.values())[0]}.")
    if tmb_rows is not None and len(assays) > 1:
        conclusion.append("TMB was computed on separate assay-specific denominators and must be interpreted per-assay, not pooled or compared directly.")
    conclusion.append(
        "These are descriptive, hypothesis-generating summaries; treat driver, actionability, and "
        "cross-cohort comparisons cautiously pending the caveats above and a review of "
        "references/pitfalls.md before finalizing conclusions."
    )
    lines.append(" ".join(conclusion))
    lines.append("")

    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Cancer genomics quick plots from mutation, CNV, signature and TMB tables.")
    parser.add_argument("--mutation-frequency", required=True)
    parser.add_argument("--cnv", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--tmb", default=None, help="Optional per-sample TMB table: sample,assay,panel_size_mb,mutation_count")
    parser.add_argument("--out", default="output")
    args = parser.parse_args()
    apply_style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    mut = read(args.mutation_frequency)
    require_nonempty(mut, args.mutation_frequency, "Mutation frequency")
    require_columns(mut, ["gene", "count", "frequency"], args.mutation_frequency, "Mutation frequency")
    plot_mutation_frequency(mut, out, args.mutation_frequency)

    cnv = read(args.cnv)
    require_nonempty(cnv, args.cnv, "CNV")
    require_columns(cnv, ["gene"], args.cnv, "CNV")
    cnv_samples = require_sample_columns(cnv, "gene", args.cnv, "CNV")
    plot_cnv_heatmap(cnv, out, args.cnv)

    sig = read(args.signature)
    require_nonempty(sig, args.signature, "Signature exposure")
    require_columns(sig, ["signature"], args.signature, "Signature exposure")
    sig_samples = require_sample_columns(sig, "signature", args.signature, "Signature exposure")
    plot_signature_exposure(sig, out, args.signature)

    tmb_rows = None
    if args.tmb:
        tmb_rows = read_tmb(args.tmb)
        plot_tmb_summary(tmb_rows, out)
        with (out / "results" / "tmb_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["sample", "assay", "panel_size_mb", "mutation_count", "tmb_mut_per_mb"])
            for r in tmb_rows:
                writer.writerow([r["sample"], r["assay"], r["panel_size_mb"], int(r["mutation_count"]), round(r["tmb"], 3)])

    (out / "results" / "mutation_frequency.tsv").write_text(Path(args.mutation_frequency).read_text(encoding="utf-8"), encoding="utf-8")
    (out / "results" / "cnv_matrix.tsv").write_text(Path(args.cnv).read_text(encoding="utf-8"), encoding="utf-8")
    (out / "results" / "signature_exposure.tsv").write_text(Path(args.signature).read_text(encoding="utf-8"), encoding="utf-8")

    id_notes = check_sample_id_overlap(cnv_samples, sig_samples, tmb_rows)
    for note in id_notes:
        print(f"WARNING: {note}", file=sys.stderr)

    build_report(mut, cnv, sig, tmb_rows, cnv_samples, sig_samples, id_notes, out)

    summary_lines = ["item\tvalue", f"genes\t{len(mut)}", f"samples\t{len(cnv_samples)}"]
    qc = {"genes": len(mut), "samples": len(cnv_samples), "sample_id_warnings": id_notes}
    if tmb_rows:
        assays = sorted({r["assay"] for r in tmb_rows})
        for a in assays:
            vals = sorted(r["tmb"] for r in tmb_rows if r["assay"] == a)
            median = vals[len(vals) // 2]
            summary_lines.append(f"tmb_median_{a}\t{median:.3f}")
            qc[f"tmb_median_{a}"] = median
        if len(assays) > 1:
            summary_lines.append("tmb_note\tTMB computed per-assay denominator; do not pool across assays without harmonization")
    (out / "summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    (out / "run" / "qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(f"wrote cancer genomics outputs to {out}")


if __name__ == "__main__":
    main()
