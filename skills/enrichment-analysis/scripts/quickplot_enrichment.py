#!/usr/bin/env python3
import argparse
import csv
import json
import math
import textwrap
from pathlib import Path

try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "quickplot_enrichment.py requires matplotlib and numpy to export PNG/PDF. "
        "Install them with: python3 -m pip install matplotlib numpy"
    ) from exc


TEXT = "#173B5F"
GRID = "#DCEAF5"
COLORS = {"KEGG": "#4EA3F1", "GO_BP": "#4CB984", "GO_CC": "#E2A13B", "GO_MF": "#8B72D6", "Reactome": "#E45B5B"}


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
    text = Path(path).read_text(encoding="utf-8-sig")
    if not text.strip():
        raise SystemExit(f"Enrichment table '{path}' is empty — nothing to plot or summarize.")
    dialect = csv.Sniffer().sniff(text[:2048], delimiters=",\t")
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    return rows


REQUIRED_COLUMNS = ["term", "category", "gene_ratio", "padj", "count"]


def validate_rows(rows, path):
    """Fail fast with a clear message instead of a raw KeyError/ValueError deep in a plot function."""
    if not rows:
        raise SystemExit(
            f"Enrichment table '{path}' has a header but no data rows — nothing to plot or summarize."
        )
    header = set(rows[0].keys())
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise SystemExit(
            f"Enrichment table '{path}' is missing required column(s): {', '.join(missing)}. "
            f"Found columns: {', '.join(sorted(h for h in header if h))}. "
            "Expected 'term', 'category', 'gene_ratio', 'padj', 'count' — see references/input-contract.md."
        )
    for i, row in enumerate(rows, start=1):
        for col in ("padj", "gene_ratio", "count"):
            raw = row.get(col, "")
            try:
                float(raw)
            except (TypeError, ValueError):
                raise SystemExit(
                    f"Enrichment table '{path}' row {i} (term={row.get('term', '?')!r}) has a "
                    f"non-numeric '{col}' value ({raw!r}). Fix the input table and retry."
                )


def wrap(text, width=28, max_lines=2):
    """Genuinely wrap long term names onto multiple lines (references/visual-style.md:
    'Wrap long term names') instead of truncating with an ellipsis. Capped at
    `max_lines` so labels don't outgrow their allotted per-row plot height."""
    lines = textwrap.wrap(text, width=width) or [text]
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip() + "…"
    return "\n".join(lines)


def write_bubble(rows, out_base, padj_threshold=0.05, top_n=20):
    apply_style()
    rows = sorted(rows, key=lambda r: float(r["padj"]))[:top_n]
    rows = list(reversed(rows))
    # Trailing "*" marks terms passing the significance threshold, so the highlight
    # is visible on the axis itself, not just in the title/report text.
    terms = [wrap(r["term"]) + (" *" if float(r["padj"]) < padj_threshold else "") for r in rows]
    ratios = [float(r["gene_ratio"]) for r in rows]
    counts = [float(r["count"]) for r in rows]
    scores = [-math.log10(max(float(r["padj"]), 1e-12)) for r in rows]
    fig_h = max(7.0, 0.52 * len(rows) + 2.2)  # extra room per row for 2-line wrapped labels
    fig, ax = plt.subplots(figsize=(9.8, fig_h))
    sc = ax.scatter(ratios, range(len(rows)), s=[80 + c / max(counts) * 680 for c in counts], c=scores, cmap="Blues", alpha=0.9, edgecolors="#347FC2", linewidths=0.9)
    ax.set_title(f"Enrichment bubble plot (* padj < {padj_threshold:g})", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("Gene ratio", labelpad=12)
    ax.set_ylabel("Enriched term", labelpad=12)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(terms)
    ax.grid(True, axis="x", color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("-log10 adjusted p value", fontsize=16, labelpad=10)
    for c, size in [(min(counts), 80 + min(counts) / max(counts) * 680), (max(counts), 760)]:
        ax.scatter([], [], s=size, c="#D7ECFF", edgecolors="#347FC2", label=f"Count {int(c)}")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_plot(fig, out_base)


def write_bar(rows, out_base, padj_threshold=0.05, top_n=15):
    apply_style()
    rows = sorted(rows, key=lambda r: float(r["padj"]))[:top_n]
    rows = list(reversed(rows))
    terms = [wrap(r["term"]) + (" *" if float(r["padj"]) < padj_threshold else "") for r in rows]
    counts = [float(r["count"]) for r in rows]
    colors = [COLORS.get(r["category"], "#6B7C8D") for r in rows]
    fig_h = max(6.4, 0.56 * len(rows) + 2.0)  # extra room per row for 2-line wrapped labels
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    ax.barh(range(len(rows)), counts, color=colors, alpha=0.88, height=0.68)
    ax.set_title(f"Top enriched terms (* padj < {padj_threshold:g})", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel("Gene count", labelpad=12)
    ax.set_ylabel("Enriched term", labelpad=12)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(terms)
    ax.grid(True, axis="x", color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    for i, val in enumerate(counts):
        ax.text(val + max(counts) * 0.018, i, str(int(val)), va="center", fontsize=14)
    handles = []
    for category in sorted({r["category"] for r in rows}):
        handles.append(plt.Line2D([0], [0], marker="s", linestyle="", markersize=9, color=COLORS.get(category, "#6B7C8D"), label=category))
    ax.legend(handles=handles, frameon=False, loc="lower right")
    fig.tight_layout()
    save_plot(fig, out_base)


def _split_genes(raw):
    raw = raw.replace(";", "/")
    return [g.strip() for g in raw.split("/") if g.strip()]


def write_cnetplot(rows, out_base, padj_threshold=0.05, max_terms=8, max_genes=70, min_significant=3):
    """Gene-concept network: top enriched terms linked to their member genes.

    Dependency-light (matplotlib + numpy only, no networkx): terms are placed
    on an inner ring, genes on an outer ring, and edges are drawn term->gene.
    Per cns-method-guide.md this figure is only drawn "when enough significant
    terms exist"; returns False (and draws nothing) otherwise.
    """
    if not rows or "geneID" not in rows[0]:
        return False
    sig_rows = [r for r in rows if float(r["padj"]) < padj_threshold]
    if len(sig_rows) < min_significant:
        return False
    apply_style()
    rows_sorted = sorted(sig_rows, key=lambda r: float(r["padj"]))[:max_terms]
    term_genes = {r["term"]: _split_genes(r["geneID"]) for r in rows_sorted}
    all_genes = sorted({g for genes in term_genes.values() for g in genes})
    n_genes_total = len(all_genes)
    truncated = n_genes_total > max_genes
    if truncated:
        all_genes = all_genes[:max_genes]

    n_terms = len(term_genes)
    n_genes = len(all_genes)
    term_angles = np.linspace(0, 2 * math.pi, n_terms, endpoint=False)
    gene_angles = np.linspace(0, 2 * math.pi, n_genes, endpoint=False) if n_genes else np.array([])
    term_pos = {t: (0.38 * math.cos(a), 0.38 * math.sin(a)) for t, a in zip(term_genes, term_angles)}
    # Stagger alternating genes at two slightly different radii so densely packed
    # neighbors (large gene counts) don't sit on top of each other near the ring.
    gene_radius = [1.0 if i % 2 == 0 else 1.12 for i in range(n_genes)]
    gene_pos = {g: (r * math.cos(a), r * math.sin(a)) for g, a, r in zip(all_genes, gene_angles, gene_radius)}

    fig_size = 10.5 if n_genes <= 50 else 13.0
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    for term, genes in term_genes.items():
        x0, y0 = term_pos[term]
        for g in genes:
            if g not in gene_pos:
                continue
            x1, y1 = gene_pos[g]
            ax.plot([x0, x1], [y0, y1], color="#C9D6E3", linewidth=0.8, alpha=0.75, zorder=1)

    if all_genes:
        gx = [gene_pos[g][0] for g in all_genes]
        gy = [gene_pos[g][1] for g in all_genes]
        ax.scatter(gx, gy, s=60, color="#DCEAF5", edgecolors="#4EA3F1", linewidths=0.9, zorder=2)
        # Rotate each gene label to radiate outward along its own angle (instead of
        # a fixed horizontal orientation) so densely packed labels near the top/
        # bottom of the ring stack along their short axis instead of overlapping
        # edge-to-edge — the same trick circular gene-network plots use.
        label_gap = 0.07
        for g, radius in zip(all_genes, gene_radius):
            x, y = gene_pos[g]
            angle_deg = math.degrees(math.atan2(y, x))
            on_right = math.cos(math.radians(angle_deg)) >= 0
            ha = "left" if on_right else "right"
            rot = angle_deg if on_right else angle_deg + 180
            lx = math.cos(math.radians(angle_deg)) * (radius + label_gap)
            ly = math.sin(math.radians(angle_deg)) * (radius + label_gap)
            ax.text(lx, ly, g, fontsize=7.5, ha=ha, va="center", rotation=rot, rotation_mode="anchor", color=TEXT)

    term_sizes = [220 + len(genes) * 14 for genes in term_genes.values()]
    tx = [term_pos[t][0] for t in term_genes]
    ty = [term_pos[t][1] for t in term_genes]
    term_colors = {t: COLORS.get(next((r["category"] for r in rows_sorted if r["term"] == t), ""), "#6B7C8D") for t in term_genes}
    ax.scatter(tx, ty, s=term_sizes, color=[term_colors[t] for t in term_genes], edgecolors="white", linewidths=1.4, zorder=3)
    # Term labels sit outside their (small) node marker in a bordered callout box,
    # rather than inside the circle, so long term names never get clipped.
    label_radius = 0.68
    for t, a in zip(term_genes, term_angles):
        lx, ly = label_radius * math.cos(a), label_radius * math.sin(a)
        ax.plot([term_pos[t][0], lx], [term_pos[t][1], ly], color=term_colors[t], linewidth=1.1, alpha=0.8, zorder=3)
        ax.text(
            lx, ly, wrap(t, 20, max_lines=3), fontsize=11, ha="center", va="center", fontweight="bold", color=term_colors[t],
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=term_colors[t], linewidth=1.3, alpha=0.95),
            zorder=5,
        )

    ax.set_title("Gene-concept network (top enriched terms)", pad=18, loc="center", fontweight="bold")
    lim = 1.3
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    if truncated:
        # Figure-level caption (not axes data coordinates) so it can never collide
        # with the outermost gene labels regardless of ring density.
        fig.subplots_adjust(bottom=0.08)
        fig.text(
            0.5, 0.02, f"Showing {n_genes} of {n_genes_total} member genes (truncated for readability).",
            fontsize=11, ha="center", va="bottom", color=TEXT, style="italic",
        )
    save_plot(fig, out_base)
    return True


def _jaccard_redundancy_flags(rows, padj_threshold=0.05, jaccard_threshold=0.5, top_n=10):
    """Simple heuristic for pitfalls.md #3 (redundant GO/KEGG terms): pairwise Jaccard
    overlap of member-gene sets among the top significant terms. Only meaningful when a
    `geneID` column is present; returns [] otherwise."""
    if not rows or "geneID" not in rows[0]:
        return []
    sig = [r for r in rows if float(r["padj"]) < padj_threshold and r.get("geneID")]
    sig = sorted(sig, key=lambda r: float(r["padj"]))[:top_n]
    flags = []
    for i in range(len(sig)):
        gi = set(_split_genes(sig[i]["geneID"]))
        if not gi:
            continue
        for j in range(i + 1, len(sig)):
            gj = set(_split_genes(sig[j]["geneID"]))
            if not gj:
                continue
            union = gi | gj
            jac = len(gi & gj) / len(union) if union else 0.0
            if jac >= jaccard_threshold:
                flags.append((sig[i]["term"], sig[j]["term"], jac, len(gi & gj)))
    flags.sort(key=lambda f: f[2], reverse=True)
    return flags


def _small_overlap_flags(rows, padj_threshold=0.05, min_count=3):
    """pitfalls.md #6: a 'significant' term whose gene overlap/count is tiny (<3) is
    statistically nominal, not a coherent biological signal — flag it explicitly."""
    sig = [r for r in rows if float(r["padj"]) < padj_threshold]
    return [r for r in sig if float(r["count"]) < min_count]


def build_report(rows, cnet_drawn, universe_size, mapped_genes, input_genes, padj_threshold=0.05, top_n=5):
    n_terms = len(rows)
    categories = sorted({r.get("category", "") for r in rows if r.get("category")})
    by_padj = sorted(rows, key=lambda r: float(r["padj"]))
    sig = [r for r in by_padj if float(r["padj"]) < padj_threshold]
    n_sig = len(sig)

    lines = ["# Enrichment quicklook report", ""]
    lines += [
        "## Input",
        "",
        "- Route: visualization of an already-computed enrichment result table (no new enrichment test was run here).",
        f"- Terms in table: {n_terms}",
        f"- Significance threshold used for this report/plots: padj < {padj_threshold:g}",
        f"- Database categories present: {', '.join(categories) if categories else 'not labeled'}",
    ]
    if universe_size is not None:
        lines.append(f"- Background/universe size used upstream (as supplied): {universe_size}")
    else:
        lines.append(
            "- Background/universe size: not provided to this quicklook script. Before trusting these "
            "results, confirm the upstream analysis used an explicit tested-gene universe rather than a "
            "whole-genome default (references/pitfalls.md #1)."
        )
    lines += ["", "## Results", ""]
    lines.append(f"- Significant terms (padj < {padj_threshold:g}): {n_sig} of {n_terms}")

    top_source = sig if n_sig > 0 else by_padj
    top_label = "Top terms by adjusted p value" if n_sig > 0 else (
        f"No term reached padj < {padj_threshold:g}; top terms by padj are listed for reference only (not claimed as significant)"
    )
    lines.append(f"- {top_label}:")
    for r in top_source[:top_n]:
        lines.append(
            f"  - {r['term']} ({r.get('category', 'NA')}) — padj={float(r['padj']):.3g}, "
            f"gene overlap count={int(float(r['count']))}, gene ratio={r.get('gene_ratio', 'NA')}"
        )
    lines.append(
        f"- Gene-concept network (cnetplot): {'drawn' if cnet_drawn else 'skipped'}"
        + ("" if cnet_drawn else f" (needs a 'geneID' column and >=3 terms with padj < {padj_threshold:g}).")
    )

    caveats = []
    if mapped_genes is not None and input_genes:
        rate = mapped_genes / input_genes
        if rate < 0.9:
            caveats.append(
                f"ID mapping rate is low: {mapped_genes}/{input_genes} genes mapped ({rate:.1%}). "
                "A low mapping rate can silently shrink both the foreground and background and bias "
                "results (references/pitfalls.md #4); check `unmapped_genes.tsv` if available."
            )
    small = _small_overlap_flags(rows, padj_threshold=padj_threshold)
    if small:
        names = ", ".join(f"{r['term']} (n={int(float(r['count']))})" for r in small[:5])
        caveats.append(
            f"{len(small)} significant term(s) have a gene overlap count below 3, which is statistically "
            f"nominal rather than a coherent biological signal (references/pitfalls.md #6): {names}."
        )
    redundant = _jaccard_redundancy_flags(rows, padj_threshold=padj_threshold)
    if redundant:
        examples = "; ".join(
            f"'{a}' vs '{b}' share {n} genes (Jaccard={j:.2f})" for a, b, j, n in redundant[:3]
        )
        caveats.append(
            f"{len(redundant)} pair(s) of top significant terms show high gene-set overlap, suggesting "
            f"redundant (parent/child or near-duplicate) biology rather than independent findings "
            f"(references/pitfalls.md #3): {examples}. Consider clusterProfiler::simplify()/REVIGO-style "
            "reduction before reporting each term as an independent hit."
        )
    if caveats:
        lines += ["", "## Caveats", ""]
        lines += [f"- {c}" for c in caveats]

    lines += ["", "## Conclusion", ""]
    if n_sig == 0:
        lines.append(
            f"None of the {n_terms} tested terms reached padj < {padj_threshold:g} in this table; the terms listed above "
            "are the closest to significance but should not be reported as enriched findings without "
            "independent support."
        )
    else:
        top_name = top_source[0]["term"]
        top_padj = float(top_source[0]["padj"])
        top_count = int(float(top_source[0]["count"]))
        lines.append(
            f"{n_sig} of {n_terms} tested terms were significant at padj < {padj_threshold:g}, led by '{top_name}' "
            f"(padj={top_padj:.3g}, {top_count} overlapping genes). "
            + ("No redundancy or small-overlap issues were flagged above this threshold. " if not caveats else "See caveats above before treating every listed term as an independent biological finding. ")
            + "This is a visualization-only quicklook of a result table that was computed elsewhere; it does "
            "not re-run the statistical test or ID mapping (references/pitfalls.md #9 on not conflating "
            "ORA/GSEA interpretation)."
        )
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Enrichment plots from result table with PNG/PDF export.")
    parser.add_argument("--table", required=True)
    parser.add_argument("--out", default="output")
    parser.add_argument(
        "--padj-threshold", type=float, default=0.05,
        help="Adjusted p-value cutoff for which terms are highlighted/counted as significant in the "
        "barplot/bubble/cnetplot and in report.md (default: 0.05). Independent of any DE-table padj "
        "cutoff used upstream to build the input table.",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Max terms shown in the bubble plot (default: 20).")
    parser.add_argument("--universe-size", type=int, default=None, help="Background/universe size used upstream, for report transparency (pitfalls.md #1).")
    parser.add_argument("--mapped-genes", type=int, default=None, help="Optional: number of genes successfully ID-mapped upstream, for the mapping-rate caveat.")
    parser.add_argument("--input-genes", type=int, default=None, help="Optional: number of genes submitted to ID mapping upstream (denominator for --mapped-genes).")
    args = parser.parse_args()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    rows = read_table(args.table)
    validate_rows(rows, args.table)
    write_bubble(rows, out / "figures" / "enrichment_bubble", padj_threshold=args.padj_threshold, top_n=args.top_n)
    write_bar(rows, out / "figures" / "enrichment_barplot", padj_threshold=args.padj_threshold)
    cnet_drawn = write_cnetplot(rows, out / "figures" / "enrichment_cnetplot", padj_threshold=args.padj_threshold)
    if not cnet_drawn:
        print(f"skipped enrichment_cnetplot: needs a 'geneID' column and >=3 terms with padj < {args.padj_threshold:g}")
    with (out / "results" / "enrichment_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    n_sig = sum(1 for r in rows if float(r["padj"]) < args.padj_threshold)
    (out / "summary.tsv").write_text(
        f"item\tvalue\nterms\t{len(rows)}\nterms_significant\t{n_sig}\npadj_threshold\t{args.padj_threshold:g}\n", encoding="utf-8"
    )
    (out / "run" / "qc.json").write_text(
        json.dumps(
            {
                "terms": len(rows),
                "terms_significant": n_sig,
                "padj_threshold": args.padj_threshold,
                "visualized_existing_table": True,
                "cnetplot_drawn": cnet_drawn,
                "universe_size": args.universe_size,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    report = build_report(rows, cnet_drawn, args.universe_size, args.mapped_genes, args.input_genes, padj_threshold=args.padj_threshold)
    (out / "report.md").write_text(report, encoding="utf-8")
    print(f"wrote enrichment quicklook outputs to {out}")


if __name__ == "__main__":
    main()
