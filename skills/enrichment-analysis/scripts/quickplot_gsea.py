#!/usr/bin/env python3
import argparse
import csv
import json
import math
import textwrap
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "quickplot_gsea.py requires matplotlib to export PNG/PDF. "
        "Install it with: python3 -m pip install matplotlib"
    ) from exc


TEXT = "#173B5F"
GRID = "#DCEAF5"
RED = "#E45B5B"
BLUE = "#4EA3F1"


def read_csv(path):
    text = Path(path).read_text(encoding="utf-8-sig")
    if not text.strip():
        raise SystemExit(f"'{path}' is empty — no data rows to read.")
    return list(csv.DictReader(text.splitlines()))


def save_plot(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


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


def wrap(text, width=28, max_lines=2):
    """Genuinely wrap long pathway names onto multiple lines (references/visual-style.md:
    'Wrap long term names') instead of truncating with an ellipsis."""
    lines = textwrap.wrap(text, width=width) or [text]
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip() + "…"
    return "\n".join(lines)


def load_ranked(path):
    """Load a ranked gene list, validating columns and resolving duplicate gene IDs.

    Duplicate gene IDs are resolved by keeping the row with the largest absolute
    score (references/pitfalls.md #4's "explicit resolution rule" for many-to-one
    collisions), and the number of rows dropped this way is returned for reporting.
    """
    rows = read_csv(path)
    if not rows:
        raise SystemExit(f"Ranked gene list '{path}' has a header but no data rows.")
    header = set(rows[0].keys())
    missing = [c for c in ("gene", "score") if c not in header]
    if missing:
        raise SystemExit(
            f"Ranked gene list '{path}' is missing required column(s): {', '.join(missing)}. "
            f"Found columns: {', '.join(sorted(h for h in header if h))}."
        )
    parsed = []
    for i, row in enumerate(rows, start=1):
        gene = (row.get("gene") or "").strip()
        if not gene:
            continue
        try:
            score = float(row["score"])
        except (TypeError, ValueError):
            raise SystemExit(f"Ranked gene list '{path}' row {i} (gene={gene!r}) has a non-numeric score ({row['score']!r}).")
        if not math.isfinite(score):
            continue
        parsed.append({"gene": gene, "score": score})
    if not parsed:
        raise SystemExit(f"Ranked gene list '{path}' has no usable rows (all missing/non-numeric/non-finite scores).")
    best_by_gene = {}
    for r in parsed:
        prev = best_by_gene.get(r["gene"])
        if prev is None or abs(r["score"]) > abs(prev["score"]):
            best_by_gene[r["gene"]] = r
    n_dropped_dup = len(parsed) - len(best_by_gene)
    ranked = sorted(best_by_gene.values(), key=lambda r: r["score"], reverse=True)
    return ranked, n_dropped_dup


def load_sets(path):
    rows = read_csv(path)
    if not rows:
        raise SystemExit(f"Gene sets file '{path}' has a header but no data rows.")
    header = set(rows[0].keys())
    missing = [c for c in ("set", "gene") if c not in header]
    if missing:
        raise SystemExit(
            f"Gene sets file '{path}' is missing required column(s): {', '.join(missing)}. "
            f"Found columns: {', '.join(sorted(h for h in header if h))}."
        )
    sets = {}
    for row in rows:
        set_name = (row.get("set") or "").strip()
        gene = (row.get("gene") or "").strip()
        if not set_name or not gene:
            continue
        sets.setdefault(set_name, set()).add(gene)
    if not sets:
        raise SystemExit(f"Gene sets file '{path}' has no valid (non-empty) 'set'/'gene' rows.")
    return sets


def running_score(ranked, geneset):
    genes = [r["gene"] for r in ranked]
    scores = [abs(float(r["score"])) for r in ranked]
    hit_sum = sum(s for g, s in zip(genes, scores) if g in geneset) or 1
    miss_penalty = 1 / max(1, len(genes) - len(geneset))
    current = 0
    trace = []
    hit_positions = []
    for i, (gene, score) in enumerate(zip(genes, scores), start=1):
        if gene in geneset:
            current += score / hit_sum
            hit_positions.append(i)
        else:
            current -= miss_penalty
        trace.append(current)
    es = max(trace, key=abs)
    return trace, hit_positions, es


def leading_edge_count(trace, hit_positions, es):
    """Approximate leading-edge/core gene count: hits on the enriched side of the
    running-score peak (before the peak for a positive ES, after it for a negative ES).
    This mirrors the fgsea/GSEA-desktop leading-edge definition but is computed from
    the quicklook's own running score, not a permutation-tested statistic."""
    if not trace or not hit_positions:
        return 0
    peak_pos = max(range(len(trace)), key=lambda i: abs(trace[i])) + 1  # 1-based rank
    if es >= 0:
        return sum(1 for p in hit_positions if p <= peak_pos)
    return sum(1 for p in hit_positions if p >= peak_pos)


def draw_running(ranked, set_name, geneset, out_base):
    trace, hits, es = running_score(ranked, geneset)
    xs = list(range(1, len(trace) + 1))
    fig, axes = plt.subplots(3, 1, figsize=(9.4, 7.2), gridspec_kw={"height_ratios": [3.0, 0.7, 1.2]}, sharex=True)
    axes[0].plot(xs, trace, color=RED if es >= 0 else BLUE, linewidth=2.6)
    axes[0].axhline(0, color="#A7C6DD", linewidth=1)
    axes[0].set_ylabel("Running ES")
    axes[0].set_title(f"GSEA running score: {wrap(set_name, 60)}", loc="center", pad=18, fontweight="bold", color=TEXT)
    axes[0].grid(True, color=GRID)
    for h in hits:
        axes[1].axvline(h, color=TEXT, linewidth=0.75, alpha=0.82)
    axes[1].set_yticks([])
    axes[1].set_ylabel("Hits")
    scores = [float(r["score"]) for r in ranked]
    axes[2].fill_between(xs, scores, 0, color="#B9D7EE", alpha=0.85)
    axes[2].axhline(0, color="#A7C6DD", linewidth=1)
    axes[2].set_ylabel("Rank metric")
    axes[2].set_xlabel("Rank in ordered gene list")
    fig.tight_layout()
    save_plot(fig, out_base)


def draw_summary(ranked, sets, out_base, nes_threshold=1.0, max_sets=20):
    rows = []
    for name, genes in sets.items():
        trace, hits, es = running_score(ranked, genes)
        overlap = sum(1 for r in ranked if r["gene"] in genes)
        rows.append({
            "set": name,
            "ES": es,
            "NES": es * math.sqrt(max(overlap, 1)),
            "overlap": overlap,
            "leading_edge": leading_edge_count(trace, hits, es),
        })
    rows.sort(key=lambda r: abs(r["NES"]), reverse=True)
    rows = rows[:max_sets]
    plot_rows = list(reversed(rows))  # largest |NES| plotted at the top
    fig, ax = plt.subplots(figsize=(8.6, max(4.8, len(plot_rows) * 0.62 + 2)))  # extra room per row for 2-line wrapped labels
    colors = [RED if r["NES"] >= 0 else BLUE for r in plot_rows]
    ax.scatter(
        [r["NES"] for r in plot_rows], range(len(plot_rows)),
        s=[90 + r["overlap"] * 18 for r in plot_rows], color=colors, alpha=0.88, edgecolors="#1F2937", linewidths=0.5,
    )
    ax.axvline(0, color="#A7C6DD", linewidth=1)
    ax.axvline(nes_threshold, color="#8A99A8", linewidth=1, linestyle="--")
    ax.axvline(-nes_threshold, color="#8A99A8", linewidth=1, linestyle="--")
    # Trailing "*" marks gene sets passing the notable-NES threshold, so the
    # highlight is visible on the axis itself, not just in the title/report text.
    labels = [wrap(r["set"]) + (" *" if abs(r["NES"]) >= nes_threshold else "") for r in plot_rows]
    ax.set_yticks(range(len(plot_rows)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("NES-like score")
    ax.set_ylabel("Gene set")
    ax.set_title(f"GSEA summary (* |NES| >= {nes_threshold:g})", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.grid(True, axis="x", color=GRID)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_plot(fig, out_base)
    return rows


NES_NOTABLE_DEFAULT = 1.0  # descriptive-only effect-size heuristic; NOT a significance threshold (no p-value is computed)


def _gene_set_redundancy_flags(sets, notable_names, jaccard_threshold=0.5):
    """pitfalls.md #3-style heuristic applied to gene-set membership (e.g. overlapping
    MSigDB/Hallmark sets): pairwise Jaccard overlap among the notable sets' full gene
    membership (not just the ranked-list overlap)."""
    names = [n for n in notable_names if n in sets][:10]
    flags = []
    for i in range(len(names)):
        gi = sets[names[i]]
        for j in range(i + 1, len(names)):
            gj = sets[names[j]]
            union = gi | gj
            jac = len(gi & gj) / len(union) if union else 0.0
            if jac >= jaccard_threshold:
                flags.append((names[i], names[j], jac, len(gi & gj)))
    flags.sort(key=lambda f: f[2], reverse=True)
    return flags


def build_report(rows, sets, n_genes, n_dropped_dup, n_dropped_nonfinite, highlighted_set, nes_threshold=NES_NOTABLE_DEFAULT, top_n=5):
    up = [r for r in rows if r["NES"] >= nes_threshold]
    down = [r for r in rows if r["NES"] <= -nes_threshold]
    by_abs_nes = sorted(rows, key=lambda r: abs(r["NES"]), reverse=True)

    lines = ["# GSEA quicklook report", ""]
    lines += [
        "## Input",
        "",
        "- Route: dependency-light running-score quicklook (no permutation test; NOT a substitute for fgsea/clusterProfiler formal GSEA).",
        f"- Ranked genes used: {n_genes}",
        f"- Gene sets tested: {len(sets)}",
    ]
    if n_dropped_dup or n_dropped_nonfinite:
        lines.append(
            f"- Genes dropped before ranking: {n_dropped_dup} duplicate gene ID(s) (kept the occurrence with "
            f"the largest absolute score, per references/pitfalls.md #4), {n_dropped_nonfinite} non-finite score(s)."
        )
    lines += ["", "## Results", ""]
    lines.append(
        f"- Gene sets with a notable NES-like score (|NES| >= {nes_threshold:g}, descriptive only, no p-value): "
        f"{len(up)} skewed up, {len(down)} skewed down, out of {len(rows)} tested."
    )
    lines.append("- Top gene sets by |NES-like score|:")
    for r in by_abs_nes[:top_n]:
        lines.append(
            f"  - {r['set']} — NES-like={r['NES']:.2f}, ES={r['ES']:.2f}, ranked-list overlap={r['overlap']} genes, "
            f"leading-edge (approx.)={r['leading_edge']} genes"
        )
    if highlighted_set is not None:
        lines.append(f"- Running-score plot drawn for: {highlighted_set}")

    caveats = [
        "No permutation or adaptive-multilevel p value is computed by this quicklook (references/pitfalls.md #13); "
        "treat the NES-like scores above as descriptive ranking only, not statistical significance. Use "
        "scripts/run_gsea_fgsea.R for a formally testable padj."
    ]
    small_edge = [r for r in by_abs_nes[:top_n] if abs(r["NES"]) >= nes_threshold and r["leading_edge"] < 3]
    if small_edge:
        names = ", ".join(f"{r['set']} (leading-edge={r['leading_edge']})" for r in small_edge)
        caveats.append(
            f"Leading-edge gene count is below 3 for: {names}. A handful of genes driving the whole running-score "
            f"peak is statistically nominal, not a coordinated pathway signal (references/pitfalls.md #6)."
        )
    notable_names = [r["set"] for r in up] + [r["set"] for r in down]
    redundant = _gene_set_redundancy_flags(sets, notable_names)
    if redundant:
        examples = "; ".join(f"'{a}' vs '{b}' share {n} genes (Jaccard={j:.2f})" for a, b, j, n in redundant[:3])
        caveats.append(
            f"{len(redundant)} pair(s) of notable gene sets share most of their member genes, suggesting "
            f"overlapping/redundant gene-set definitions rather than independent pathway hits "
            f"(references/pitfalls.md #3): {examples}."
        )
    lines += ["", "## Caveats", ""]
    lines += [f"- {c}" for c in caveats]

    lines += ["", "## Conclusion", ""]
    if not up and not down:
        lines.append(
            f"None of the {len(rows)} tested gene sets showed a notable NES-like skew (|NES| >= {nes_threshold:g}) "
            "in this ranked list under this quicklook heuristic."
        )
    else:
        top = by_abs_nes[0]
        direction = "up" if top["NES"] >= 0 else "down"
        lines.append(
            f"'{top['set']}' shows the strongest skew ({direction}, NES-like={top['NES']:.2f}) with "
            f"{top['leading_edge']} approximate leading-edge genes out of {top['overlap']} ranked-list overlap genes. "
            f"Overall, {len(up)} gene set(s) skew toward the up-ranked end and {len(down)} toward the down-ranked "
            "end of the list. These are descriptive quicklook results without a formal p value; confirm any "
            "candidate pathway with scripts/run_gsea_fgsea.R before reporting it as significant."
        )
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="GSEA-style running score plots from ranked genes and gene sets.")
    parser.add_argument("--ranked", required=True)
    parser.add_argument("--gene-sets", required=True)
    parser.add_argument("--set-name", default="")
    parser.add_argument("--out", default="output")
    parser.add_argument(
        "--nes-threshold", type=float, default=NES_NOTABLE_DEFAULT,
        help="|NES-like score| cutoff for which gene sets are highlighted/counted as notable in the "
        f"summary dotplot and report.md (default: {NES_NOTABLE_DEFAULT:g}). Descriptive only — this "
        "quicklook computes no p-value; use run_gsea_fgsea.R's --padj-threshold for a formal cutoff.",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Max gene sets shown in the summary dotplot (default: 20).")
    args = parser.parse_args()
    apply_style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    ranked, n_dropped_dup = load_ranked(args.ranked)
    n_dropped_nonfinite = 0  # already filtered inside load_ranked; kept as a named field for report symmetry with run_gsea_fgsea.R
    sets = load_sets(args.gene_sets)
    set_name = args.set_name or next(iter(sets))
    if set_name not in sets:
        raise SystemExit(f"gene set not found: {set_name!r}. Available sets: {', '.join(sorted(sets))}")
    draw_running(ranked, set_name, sets[set_name], out / "figures" / "gsea_running_score")
    rows = draw_summary(ranked, sets, out / "figures" / "gsea_summary_dotplot", nes_threshold=args.nes_threshold, max_sets=args.top_n)
    with (out / "results" / "gsea_quick_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["set", "ES", "NES", "overlap", "leading_edge"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (out / "summary.tsv").write_text(
        f"item\tvalue\ngenes\t{len(ranked)}\ngene_sets\t{len(sets)}\nnes_threshold\t{args.nes_threshold:g}\n", encoding="utf-8"
    )
    (out / "run" / "qc.json").write_text(
        json.dumps(
            {
                "genes": len(ranked),
                "gene_sets": len(sets),
                "quicklook": True,
                "genes_dropped_duplicate": n_dropped_dup,
                "nes_threshold": args.nes_threshold,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (out / "methods.md").write_text(
        "# Methods\n\nGSEA-style quicklook was generated from a ranked gene list and gene-set membership table. "
        "The running enrichment score follows the ranked-list hit/miss intuition of GSEA and is intended for visualization; "
        "formal inference should use fgsea/clusterProfiler with permutation or adaptive multilevel statistics. "
        f"Duplicate gene IDs in the ranked list ({n_dropped_dup} found) were resolved by keeping the occurrence "
        f"with the largest absolute score (references/pitfalls.md #4). Gene sets with |NES-like score| >= "
        f"{args.nes_threshold:g} were treated as notable in this run (descriptive only, no p-value).\n",
        encoding="utf-8",
    )
    report = build_report(rows, sets, len(ranked), n_dropped_dup, n_dropped_nonfinite, set_name, nes_threshold=args.nes_threshold)
    # Named report_gsea.md, not report.md: quickplot_enrichment.py (ORA) also writes report.md, and
    # running both routes into the same --out directory (a very normal thing to do -- "give me both
    # GO/KEGG ORA and GSEA") would otherwise have GSEA silently overwrite the ORA report. Matches the
    # contract already documented in references/deliverables.md for the R route (run_gsea_fgsea.R).
    (out / "report_gsea.md").write_text(report, encoding="utf-8")
    print(f"wrote GSEA outputs to {out}")


if __name__ == "__main__":
    main()
