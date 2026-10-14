#!/usr/bin/env python3
"""Dependency-light spatial transcriptomics quicklook: feature/domain/QC maps and result tables
from exported spot and expression tables (see references/input-contract.md)."""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except Exception as exc:
    raise SystemExit("Install matplotlib: python3 -m pip install matplotlib") from exc


PALETTE = ["#4EA3F1", "#E45B5B", "#4CB984", "#E2A13B", "#8B72D6", "#56BFC0"]
TEXT = "#173B5F"
GRID = "#DCEAF5"
BLUE = "#4EA3F1"
RED = "#E45B5B"
# references/visual-style.md: "Feature maps: light grey/blue low values to saturated red high values."
FEATURE_CMAP = LinearSegmentedColormap.from_list("spatial_feature", [GRID, BLUE, RED])


def read_csv(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def pearson_corr(xs, ys):
    """Pearson correlation over paired values, skipping any pair with a missing entry.
    Returns None if fewer than 3 usable pairs or either series is constant (undefined)."""
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return None
    mx = statistics.mean(p[0] for p in pairs)
    my = statistics.mean(p[1] for p in pairs)
    sx = sum((p[0] - mx) ** 2 for p in pairs)
    sy = sum((p[1] - my) ** 2 for p in pairs)
    if sx == 0 or sy == 0:
        return None
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pairs)
    return sxy / math.sqrt(sx * sy)


def morans_i(xs, ys, values, k=6, max_n=1500):
    """Lightweight k-nearest-neighbor Moran's I as a quicklook spatial-autocorrelation proxy.
    Not a substitute for a formal SVG method (SPARK-X/SpatialDE/Moran's I in squidpy on
    depth-corrected data) -- skipped above max_n spots to keep the quicklook fast, and
    returns None when expression is constant (denominator zero) or too few valid points."""
    triplets = [(x, y, v) for x, y, v in zip(xs, ys, values) if x is not None and y is not None and v is not None]
    n = len(triplets)
    if n < max(10, k + 1) or n > max_n:
        return None
    px = [t[0] for t in triplets]
    py = [t[1] for t in triplets]
    mean_v = statistics.mean(t[2] for t in triplets)
    z = [t[2] - mean_v for t in triplets]
    denom = sum(v * v for v in z)
    if denom == 0:
        return None
    numerator = 0.0
    for i in range(n):
        dists = sorted(range(n), key=lambda j: (px[i] - px[j]) ** 2 + (py[i] - py[j]) ** 2 if j != i else float("inf"))
        neighbors = dists[:k]
        w = 1.0 / len(neighbors)
        numerator += sum(w * z[i] * z[j] for j in neighbors)
    return numerator / denom


def describe_hotspot(rows, values, frac=0.1, min_n=3):
    """Describe the approximate tissue quadrant where the highest values concentrate,
    relative to the median x/y of all spots. Returns None if uninformative (all-zero/constant)."""
    pts = [(to_float(r["x"]), to_float(r["y"]), to_float(v)) for r, v in zip(rows, values)]
    pts = [p for p in pts if p[0] is not None and p[1] is not None and p[2] is not None]
    if len(pts) < min_n:
        return None
    vals = [p[2] for p in pts]
    if max(vals) == min(vals):
        return None
    med_x = statistics.median(p[0] for p in pts)
    med_y = statistics.median(p[1] for p in pts)
    top_n = max(min_n, int(len(pts) * frac))
    top = sorted(pts, key=lambda p: p[2], reverse=True)[:top_n]
    top_x = statistics.mean(p[0] for p in top)
    top_y = statistics.mean(p[1] for p in top)
    ns = "upper" if top_y > med_y else "lower"
    ew = "right" if top_x > med_x else "left"
    return f"{ns}-{ew} region of the plotted coordinate range (top {top_n} spots by value)"


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 16, "axes.titlesize": 24, "axes.labelsize": 21,
        "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15,
        "axes.edgecolor": "#B9D7EE", "axes.labelcolor": TEXT, "xtick.color": TEXT,
        "ytick.color": TEXT, "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def spatial_scatter(rows, color_values, title, color_label, out, categorical=False, flip_y=False):
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    xs = [float(r["x"]) for r in rows]
    ys = [float(r["y"]) for r in rows]
    if categorical:
        cats = sorted(set(color_values))
        cmap = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(cats)}
        for c in cats:
            idx = [i for i, v in enumerate(color_values) if v == c]
            ax.scatter([xs[i] for i in idx], [ys[i] for i in idx], s=44, color=cmap[c], label=c, alpha=0.88, linewidths=0)
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), title=color_label)
    else:
        numeric = [to_float(v, float("nan")) for v in color_values]
        sc = ax.scatter(xs, ys, c=numeric, s=44, cmap=FEATURE_CMAP, alpha=0.92, linewidths=0)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.03)
        cbar.set_label(color_label, fontsize=16)
    ax.set_title(title, loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_aspect("equal", adjustable="box")
    # references/visual-style.md + SKILL.md: preserve coordinate orientation unless a flip is requested.
    if flip_y:
        ax.invert_yaxis()
    save(fig, out)


def main():
    parser = argparse.ArgumentParser(description="Spatial transcriptomics quick plots and result tables from spot/expression tables.")
    parser.add_argument("--spots", required=True, help="spots.csv with spot,x,y[,domain,nCount,nFeature,sample]")
    parser.add_argument("--expression", required=True, help="expression.csv with spot + gene columns")
    parser.add_argument("--gene", default="GeneA", help="comma-separated target gene(s); first gene is written as the canonical spatial_feature_map")
    parser.add_argument("--domain-col", default="domain")
    parser.add_argument("--deconvolution", default=None, help="optional spot-by-celltype proportion CSV (spot + celltype columns)")
    parser.add_argument("--flip-y", action="store_true", help="invert the y axis (only when y follows image-pixel row convention with origin top-left)")
    parser.add_argument("--out", default="output")
    args = parser.parse_args()
    style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    spots = read_csv(args.spots)
    if not spots:
        raise SystemExit(f"No rows found in {args.spots}")
    required_spot_cols = [c for c in ("spot", "x", "y") if c not in spots[0]]
    if required_spot_cols:
        raise SystemExit(
            f"Spots table `{args.spots}` is missing required column(s): {', '.join(required_spot_cols)}; "
            "see references/input-contract.md (spots table needs spot, x, y)."
        )
    expr_records = read_csv(args.expression)
    if not expr_records:
        raise SystemExit(f"No rows found in {args.expression}")
    if "spot" not in expr_records[0]:
        raise SystemExit(f"Expression table `{args.expression}` is missing required column: spot.")
    expr = {r["spot"]: r for r in expr_records}
    matched = [r for r in spots if r["spot"] in expr]
    if not matched:
        raise SystemExit("No overlapping spot IDs between spots table and expression table; check references/input-contract.md.")
    n_matched_ids = len(matched)
    rows = [r for r in matched if to_float(r["x"]) is not None and to_float(r["y"]) is not None]
    n_bad_coord = n_matched_ids - len(rows)
    if not rows:
        raise SystemExit("No spots with both matching expression data and numeric x/y coordinates; check references/input-contract.md.")

    genes = [g.strip() for g in args.gene.split(",") if g.strip()]
    expr_columns = set(expr_records[0].keys()) if expr_records else set()
    missing_genes = [g for g in genes if g not in expr_columns]
    if missing_genes:
        raise SystemExit(f"Gene column(s) not found in expression table: {', '.join(missing_genes)}")

    has_domain = args.domain_col in rows[0]
    qc_metric = "nCount" if "nCount" in rows[0] else ("nFeature" if "nFeature" in rows[0] else None)

    # ---- spatial feature maps (one per requested gene; first gene keeps the canonical name) ----
    xs_all = [to_float(r["x"]) for r in rows]
    ys_all = [to_float(r["y"]) for r in rows]
    qc_numeric = [to_float(r.get(qc_metric)) for r in rows] if qc_metric else None

    gene_values = {}
    gene_stats = {}
    for i, gene in enumerate(genes):
        values = [expr[r["spot"]].get(gene, "") for r in rows]
        gene_values[gene] = values
        stem = out / "figures" / "spatial_feature_map" if i == 0 else out / "figures" / f"spatial_feature_map_{gene}"
        spatial_scatter(rows, values, f"Spatial feature: {gene}", gene, stem, flip_y=args.flip_y)

        numeric = [to_float(v) for v in values]
        valid = [v for v in numeric if v is not None]
        stats = {
            "n": len(valid),
            "min": min(valid) if valid else None,
            "max": max(valid) if valid else None,
            "mean": round(statistics.mean(valid), 4) if valid else None,
            "all_zero": bool(valid) and max(valid) == 0 == min(valid),
        }
        stats["hotspot"] = None if stats["all_zero"] or not valid else describe_hotspot(rows, values)
        stats["moran_i"] = None if stats["all_zero"] else morans_i(xs_all, ys_all, numeric)
        stats["depth_corr"] = pearson_corr(numeric, qc_numeric) if qc_numeric else None
        gene_stats[gene] = stats

    # ---- spatial domain map ----
    if has_domain:
        domain_values = [r[args.domain_col] for r in rows]
        spatial_scatter(rows, domain_values, "Spatial domains", args.domain_col, out / "figures" / "spatial_domain_map", categorical=True, flip_y=args.flip_y)

    # ---- spot QC map ----
    if qc_metric:
        qc_values = [r[qc_metric] for r in rows]
        qc_title = "Spot count QC" if qc_metric == "nCount" else "Spot feature QC"
        spatial_scatter(rows, qc_values, qc_title, qc_metric, out / "figures" / "spot_qc_map", flip_y=args.flip_y)

    # ---- optional deconvolution: dominant-type map + proportions table ----
    deconv_rows = None
    if args.deconvolution:
        deconv_records = read_csv(args.deconvolution)
        deconv_table = {r["spot"]: r for r in deconv_records}
        celltype_cols = [c for c in deconv_records[0].keys() if c != "spot"] if deconv_records else []
        deconv_rows = [r for r in rows if r["spot"] in deconv_table]
        if deconv_rows and celltype_cols:
            dominant = []
            for r in deconv_rows:
                proportions = deconv_table[r["spot"]]
                dominant.append(max(celltype_cols, key=lambda c: to_float(proportions.get(c), 0.0)))
            spatial_scatter(
                deconv_rows, dominant, "Dominant estimated cell type (deconvolution)", "dominant_type",
                out / "figures" / "deconvolution_composition_map", categorical=True, flip_y=args.flip_y,
            )
            with (out / "results" / "deconvolution_proportions.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["spot", "x", "y"] + celltype_cols + ["dominant_type"])
                for r, dom in zip(deconv_rows, dominant):
                    proportions = deconv_table[r["spot"]]
                    writer.writerow([r["spot"], r["x"], r["y"]] + [proportions.get(c, "") for c in celltype_cols] + [dom])
        else:
            deconv_rows = None

    # ---- result tables ----
    with (out / "results" / "spatial_feature_values.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["spot", "x", "y"] + genes)
        for idx, r in enumerate(rows):
            writer.writerow([r["spot"], r["x"], r["y"]] + [gene_values[g][idx] for g in genes])

    qc_fields = [f for f in ("nCount", "nFeature") if f in rows[0]]
    with (out / "results" / "spot_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["spot", "x", "y"] + qc_fields + ([args.domain_col] if has_domain else []))
        for r in rows:
            writer.writerow([r["spot"], r["x"], r["y"]] + [r.get(f, "") for f in qc_fields] + ([r[args.domain_col]] if has_domain else []))

    domains_sorted = sorted(set(r[args.domain_col] for r in rows)) if has_domain else []
    domain_sizes = {}
    small_domain_threshold = max(5, round(0.03 * len(rows)))
    if has_domain:
        with (out / "results" / "domain_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            header = ["domain", "n_spots"] + [f"mean_{g}" for g in genes] + ([f"mean_{qc_metric}"] if qc_metric else [])
            writer.writerow(header)
            for d in domains_sorted:
                d_rows = [r for r in rows if r[args.domain_col] == d]
                domain_sizes[d] = len(d_rows)
                row_out = [d, len(d_rows)]
                for g in genes:
                    vals = [v for v in (to_float(expr[r["spot"]].get(g)) for r in d_rows) if v is not None]
                    row_out.append(round(statistics.mean(vals), 4) if vals else "")
                if qc_metric:
                    vals = [v for v in (to_float(r.get(qc_metric)) for r in d_rows) if v is not None]
                    row_out.append(round(statistics.mean(vals), 2) if vals else "")
                writer.writerow(row_out)
    small_domains = [d for d in domains_sorted if domain_sizes.get(d, 0) < small_domain_threshold]

    # ---- summary.tsv, run/qc.json, report.md ----
    mean_qc = None
    if qc_metric:
        vals = [v for v in (to_float(r.get(qc_metric)) for r in rows) if v is not None]
        mean_qc = round(statistics.mean(vals), 2) if vals else None

    summary_lines = ["item\tvalue", f"spots\t{len(rows)}", f"genes\t{','.join(genes)}"]
    if has_domain:
        summary_lines.append(f"domains\t{len(domains_sorted)}")
    if mean_qc is not None:
        summary_lines.append(f"mean_{qc_metric}\t{mean_qc}")
    if deconv_rows:
        summary_lines.append(f"deconvolution_spots\t{len(deconv_rows)}")
    (out / "summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    qc_payload = {"spots": len(rows), "genes": genes, "flip_y": args.flip_y}
    if has_domain:
        qc_payload["domains"] = domains_sorted
    if mean_qc is not None:
        qc_payload[f"mean_{qc_metric}"] = mean_qc
    (out / "run" / "qc.json").write_text(json.dumps(qc_payload, indent=2), encoding="utf-8")

    dropped_ids = len(spots) - n_matched_ids

    def moran_strength(value):
        mag = abs(value)
        if mag >= 0.3:
            return "moderate-to-strong positive spatial pattern" if value > 0 else "moderate-to-strong negative (checkerboard) spatial pattern"
        if mag >= 0.1:
            return "weak-to-moderate spatial pattern"
        return "little/no detectable spatial pattern at this neighbor scale"

    # ---- per-gene feature summary lines + ranking for "top" spatially patterned genes ----
    feature_lines = []
    ranked_by_moran = []
    for gene in genes:
        st = gene_stats[gene]
        if st["all_zero"]:
            feature_lines.append(f"- **{gene}**: all-zero expression across all {st['n']} matched spots; feature map has no signal to interpret.")
            continue
        if st["n"] == 0:
            feature_lines.append(f"- **{gene}**: no parseable numeric values found in the expression table.")
            continue
        line = f"- **{gene}**: n={st['n']}, range [{st['min']:.3f}, {st['max']:.3f}], mean={st['mean']}"
        if st["hotspot"]:
            line += f"; highest values concentrate in the {st['hotspot']}"
        if st["moran_i"] is not None:
            line += f"; quicklook spatial-autocorrelation proxy (k-NN Moran's I, k=6) = {st['moran_i']:.3f} ({moran_strength(st['moran_i'])})"
            ranked_by_moran.append((gene, st["moran_i"]))
        elif len(rows) > 1500:
            line += "; Moran's I proxy skipped (dataset exceeds the quicklook size cap of 1500 spots)"
        feature_lines.append(line)
    ranked_by_moran.sort(key=lambda t: abs(t[1]), reverse=True)

    # ---- automatic caveats tied to references/pitfalls.md, only printed when they apply ----
    caveats = []
    if dropped_ids:
        caveats.append(f"- {dropped_ids} spot ID(s) in `{args.spots}` had no matching row in `{args.expression}` (or vice versa) and were excluded; check ID formatting if this count is unexpectedly high (references/input-contract.md).")
    if n_bad_coord:
        caveats.append(f"- {n_bad_coord} matched spot(s) had non-numeric x/y coordinates and were dropped before plotting.")
    if genes and any(not gene_stats[g]["all_zero"] and gene_stats[g]["n"] for g in genes):
        caveats.append(
            "- Spatially-patterned signal above (k-NN Moran's I) is a dependency-light quicklook proxy computed on the raw/as-provided expression values, "
            "not depth-corrected, variance-stabilized, or FDR-controlled across genes; a gene correlated with sequencing depth can appear spatially patterned "
            "artifactually (references/pitfalls.md #4). Use `scripts/run_spatial_seurat.R` (SCTransform + Moran's I) or a dedicated SVG method "
            "(SPARK-X/SpatialDE) for a statistically supported spatially-variable-gene call."
        )
    if qc_metric:
        for gene in genes:
            dc = gene_stats[gene]["depth_corr"]
            if dc is not None and abs(dc) >= 0.4:
                caveats.append(f"- **{gene}** expression correlates with `{qc_metric}` across matched spots (r={dc:.2f}); its spatial pattern may partly reflect a sequencing-depth/QC gradient rather than true biology (references/pitfalls.md #4) — cross-check against `spot_qc_map`.")
    if has_domain and small_domains:
        detail = ", ".join(f"{d} (n={domain_sizes[d]})" for d in small_domains)
        caveats.append(f"- Domain(s) with very few spots (< {small_domain_threshold}): {detail}; per-domain means and spatial patterns for these domains are unstable and low-confidence (references/pitfalls.md #2-3).")
    if has_domain and len(domains_sorted) == 1:
        caveats.append("- Only one domain/cluster label is present in the input; the domain map and any domain-level comparison are uninformative for spatial heterogeneity.")
    if deconv_rows:
        caveats.append("- Deconvolution proportions were plotted/tabulated from the supplied table as-is; this script does not confirm that the single-cell reference used to generate them matches this tissue's type, disease state, or developmental stage. Treat proportions as approximate/qualitative until reference-tissue match is confirmed (references/pitfalls.md #9), and treat any dominant-type call as a per-spot majority estimate, not a single-cell identity (references/pitfalls.md #1).")

    # ---- plain-language conclusion ----
    conclusion_parts = [
        f"Of {len(spots)} spots in the input table, {len(rows)} had matching expression data and valid coordinates and were used in this analysis"
        + (f" across {len(domains_sorted)} domain(s)/cluster(s)" if has_domain else "") + "."
    ]
    if ranked_by_moran:
        top_gene, top_i = ranked_by_moran[0]
        conclusion_parts.append(
            f"Among the requested gene(s), **{top_gene}** showed the strongest quicklook spatial-autocorrelation signal (Moran's I={top_i:.3f})"
            + (f", concentrated in the {gene_stats[top_gene]['hotspot']}" if gene_stats[top_gene]["hotspot"] else "")
            + "."
        )
    if deconv_rows:
        conclusion_parts.append(f"A dominant estimated cell type was assigned for {len(deconv_rows)} spots from the supplied deconvolution table.")
    conclusion_parts.append(
        "These are exploratory, descriptive results from a dependency-light quicklook pass, not statistically validated spatial-domain or "
        "spatially-variable-gene calls; address the caveats below (and re-run the formal Seurat/Squidpy route where relevant) before drawing biological conclusions."
    )

    report_lines = [
        "# Spatial transcriptomics quicklook",
        "",
        "## Input summary",
        f"- Spots table: `{args.spots}` ({len(spots)} rows)",
        f"- Expression table: `{args.expression}` ({len(expr_records)} rows)",
        f"- Spots with matched expression and valid coordinates: {len(rows)} (of {len(spots)} total; {dropped_ids} unmatched by ID, {n_bad_coord} dropped for non-numeric coordinates)",
        "",
        "## QC observations",
    ]
    if mean_qc is not None:
        report_lines.append(f"- Mean {qc_metric} across the {len(rows)} matched spots: {mean_qc}")
    else:
        report_lines.append("- No `nCount`/`nFeature` column found; spot QC map/table skipped.")
    report_lines += ["", "## Spatial feature maps"] + feature_lines
    if ranked_by_moran:
        top_list = ", ".join(f"{g} (I={v:.3f})" for g, v in ranked_by_moran[:3])
        report_lines.append(f"- Top candidate spatially patterned gene(s) by quicklook Moran's I: {top_list}.")
    report_lines += ["", "## Spatial domains"]
    if has_domain:
        size_list = ", ".join(f"{d}: {domain_sizes[d]}" for d in domains_sorted)
        report_lines.append(f"- Domains in `{args.domain_col}` column ({len(domains_sorted)} total): {size_list} (labels taken as-is, not re-clustered here).")
    else:
        report_lines.append(f"- No `{args.domain_col}` column found; spatial domain map skipped.")
    if deconv_rows:
        report_lines += [
            "", "## Deconvolution",
            f"- Dominant estimated cell type mapped for {len(deconv_rows)} spots from `{args.deconvolution}`; see results/deconvolution_proportions.tsv.",
        ]
    report_lines += ["", "## Conclusion", " ".join(conclusion_parts)]
    report_lines += ["", "## Caveats specific to this run"]
    report_lines += caveats if caveats else ["- No additional run-specific caveats triggered beyond the general limitations below."]
    report_lines += [
        "",
        "## Limitations",
        "- This is a lightweight, dependency-light exploratory pass on exported tables; it does not perform image alignment, SCTransform normalization, or model-based spatially variable gene testing (use scripts/run_spatial_seurat.R for the formal Seurat route when dependencies are available).",
        "- Domain labels are taken as-is from the input table, not derived from spatial-aware clustering; treat as descriptive, not inferential (references/pitfalls.md #2).",
        "- Coordinate orientation is preserved by default (pass --flip-y only if y follows an image-pixel row convention); verify against any available histology image before interpreting spatial position (references/pitfalls.md #6).",
        "- Any dominant-type call from --deconvolution is a per-spot majority estimate, not a single-cell identity (references/pitfalls.md #1).",
    ]
    (out / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"wrote spatial outputs to {out}")


if __name__ == "__main__":
    main()
