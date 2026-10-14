#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit("Install matplotlib: python3 -m pip install matplotlib") from exc

try:
    from scipy import stats as scipy_stats
except Exception:
    scipy_stats = None


TEXT = "#173B5F"
GRID = "#DCEAF5"
GROUP_COLORS = ["#4EA3F1", "#E45B5B", "#4CB984", "#E2A13B", "#8B72D6"]
NA_TOKENS = {"", "na", "n/a", "nan", "null", "none", "-"}
SAMPLE_SIZE_CAUTION = 20
SAMPLE_SIZE_EXPLORATORY = 15
DEFAULT_PADJ_THRESHOLD = 0.05


def read_csv(path):
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError:
        raise SystemExit(f"Input file not found: {path}")
    lines = text.splitlines()
    # Auto-detect tab vs comma so .tsv outputs from run_wgcna_gsva.R (e.g. hub_genes.tsv)
    # can be read with the same helper as the comma-separated quickplot inputs.
    delimiter = "\t" if lines and lines[0].count("\t") > lines[0].count(",") else ","
    rows = list(csv.DictReader(lines, delimiter=delimiter))
    if not rows:
        raise SystemExit(f"{path} has no data rows (header only or empty file).")
    return rows


def require_columns(rows, required_any_of, path, label):
    """Validate that each column group in `required_any_of` has at least one match
    in the table header, raising a clear message (not a bare KeyError) otherwise.

    `required_any_of` is a list where each item is either a column name (required)
    or a tuple of alternative column names (at least one required).
    """
    found = set(rows[0].keys())
    missing = []
    for spec in required_any_of:
        alts = spec if isinstance(spec, tuple) else (spec,)
        if not any(a in found for a in alts):
            missing.append(" or ".join(alts))
    if missing:
        raise SystemExit(
            f"{label} ({path}) is missing required column(s): {', '.join(missing)}. "
            f"Found columns: {', '.join(found)}."
        )


def to_float(value, context=""):
    """Parse a numeric cell, treating blank/NA-like tokens as missing (None) instead
    of crashing, and raising a clear SystemExit (not a bare ValueError/traceback) for
    genuinely unparseable values.
    """
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in NA_TOKENS:
        return None
    try:
        return float(text)
    except ValueError:
        raise SystemExit(f"Could not parse numeric value {value!r}{(' in ' + context) if context else ''}; expected a number or a blank/NA marker.")


def bh_adjust(pvalues):
    """Benjamini-Hochberg FDR correction (pure Python, no extra dependency).

    Applied across the full set of p-values passed in (e.g. the whole module x trait
    matrix), matching references/pitfalls.md #3: correcting per-trait only understates
    the true false-discovery rate.
    """
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = pvalues[i] * n / rank
        running_min = min(running_min, val)
        adjusted[i] = min(running_min, 1.0)
    return adjusted


def significance_stars(p, threshold=DEFAULT_PADJ_THRESHOLD):
    """Star tiers scale off the single significance `threshold` (default 0.05) so the
    conventional ***/**/* nesting stays consistent for any custom cutoff: the loosest
    tier is the threshold itself, with the tighter tiers at threshold/5 and
    threshold/50 (matching the default 0.05/0.01/0.001 convention exactly)."""
    if p is None:
        return ""
    if p < threshold / 50:
        return "***"
    if p < threshold / 5:
        return "**"
    if p < threshold:
        return "*"
    return ""


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 16, "axes.titlesize": 24, "axes.labelsize": 21, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15, "axes.edgecolor": "#B9D7EE", "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT, "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42})


def prepare_module_trait_significance(rows):
    """Fill in a `padj` column (BH-FDR across the whole module x trait matrix) when the
    input table only carries raw `pvalue`, and leave a supplied `padj` untouched.

    Missing/NA p-values are carried through as None rather than crashing; BH-FDR is
    computed only over the rows with a usable p-value (references/pitfalls.md #3).
    Mutates and returns `rows` so both the plotted stars and the exported TSV reflect
    the corrected values.
    """
    has_padj_col = rows and "padj" in rows[0]
    parsed_padj = [to_float(r.get("padj"), context="'padj' column") for r in rows] if has_padj_col else [None] * len(rows)
    has_padj = has_padj_col and any(p is not None for p in parsed_padj)
    if not has_padj and "pvalue" in rows[0]:
        parsed_pvalue = [to_float(r.get("pvalue"), context="'pvalue' column") for r in rows]
        valid_idx = [i for i, p in enumerate(parsed_pvalue) if p is not None]
        if valid_idx:
            adjusted = bh_adjust([parsed_pvalue[i] for i in valid_idx])
            for i, p in zip(valid_idx, adjusted):
                parsed_padj[i] = p
            has_padj = True
    for r, p in zip(rows, parsed_padj):
        r["padj"] = "" if p is None else f"{p:.6g}"
    return rows, has_padj


def plot_module_trait(rows, out, padj_threshold=DEFAULT_PADJ_THRESHOLD):
    style()
    require_columns(rows, [("module",), ("trait",), ("correlation",)], "--module-trait input", "module-trait table")
    rows, has_padj = prepare_module_trait_significance(rows)
    parsed_corr = [to_float(r.get("correlation"), context="'correlation' column") for r in rows]
    if not any(v is not None for v in parsed_corr):
        raise SystemExit("module-trait table contains no valid numeric values in the 'correlation' column (all blank/NA); nothing to plot.")
    for r, v in zip(rows, parsed_corr):
        r["correlation"] = "" if v is None else str(v)
    modules = sorted({r["module"] for r in rows})
    traits = sorted({r["trait"] for r in rows})
    corr_lookup = dict(zip(((r["module"], r["trait"]) for r in rows), parsed_corr))
    padj_lookup = {(r["module"], r["trait"]): to_float(r["padj"]) for r in rows} if has_padj else {}
    mat = [[corr_lookup.get((m, t)) if corr_lookup.get((m, t)) is not None else float("nan") for t in traits] for m in modules]
    n_na_cells = sum(1 for row in mat for v in row if v != v)  # NaN != NaN
    fig, ax = plt.subplots(figsize=(max(6, len(traits) * 1.4 + 3), max(4.8, len(modules) * 0.48 + 2.2)))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#DADFE3")
    im = ax.imshow(mat, cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    ax.set_title("Module-trait correlation", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xticks(range(len(traits)))
    ax.set_xticklabels(traits, rotation=35, ha="right")
    ax.set_yticks(range(len(modules)))
    ax.set_yticklabels(modules)
    for i, module in enumerate(modules):
        for j, trait in enumerate(traits):
            corr = corr_lookup.get((module, trait))
            if corr is None:
                ax.text(j, i, "NA", ha="center", va="center", fontsize=12, color=TEXT)
                continue
            label = f"{corr:.2f}"
            if has_padj:
                label += significance_stars(padj_lookup.get((module, trait)), padj_threshold)
            ax.text(j, i, label, ha="center", va="center", fontsize=12, color=TEXT)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("Correlation", fontsize=16)
    if has_padj:
        # Anchor to the colorbar (top-right, away from the rotated x tick labels) to
        # avoid the label-collision failure mode called out in visual-style.md.
        legend_text = (f"*p<{padj_threshold:g}\n**p<{padj_threshold / 5:g}\n"
                        f"***p<{padj_threshold / 50:g}\n(BH-FDR)")
        cbar.ax.text(0.5, 1.03, legend_text, transform=cbar.ax.transAxes,
                     ha="center", va="bottom", fontsize=10, color=TEXT, linespacing=1.3)
    save(fig, out / "module_trait_heatmap")
    meta = {
        "has_padj": has_padj,
        "n_na_cells": n_na_cells,
        "n_total_cells": len(rows),
        "n_modules": len(modules),
        "n_traits": len(traits),
    }
    return rows, meta


def plot_pathways(rows, out):
    style()
    require_columns(rows, [("pathway",)], "--pathway-scores input", "pathway score table")
    samples = [k for k in rows[0].keys() if k != "pathway"]
    if not samples:
        raise SystemExit("pathway score table has a 'pathway' column but no sample columns; nothing to plot.")
    pathways = [r["pathway"] for r in rows]
    mat = [[to_float(r[s], context=f"pathway '{r['pathway']}' / sample '{s}'") for s in samples] for r in rows]
    if any(v is None for row in mat for v in row):
        raise SystemExit("pathway score table contains blank/non-numeric values; fill or drop them before plotting (every pathway x sample cell must be numeric).")
    lim = max((abs(v) for row in mat for v in row), default=1.0) or 1.0
    fig, ax = plt.subplots(figsize=(max(8.2, len(samples) * 0.5 + 3), max(5.4, len(pathways) * 0.44 + 2.4)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_title("Pathway activity heatmap", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=35, ha="right")
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(pathways)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Pathway score", fontsize=16)
    save(fig, out / "pathway_activity_heatmap")

    means = [sum(float(r[s]) for s in samples) / len(samples) for r in rows]
    fig, ax = plt.subplots(figsize=(8.4, max(4.8, len(pathways) * 0.42 + 2)))
    colors = ["#E45B5B" if v > 0 else "#4EA3F1" for v in means]
    ax.barh(pathways, means, color=colors, alpha=0.86)
    ax.axvline(0, color="#A7C6DD", linewidth=1)
    ax.set_title("Mean pathway activity", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xlabel("Mean score")
    ax.grid(True, axis="x", color=GRID)
    save(fig, out / "pathway_activity_barplot")


def plot_pathway_boxplot(pathway_rows, sample_groups, out):
    """Grouped boxplot of per-sample pathway scores (Required Figures: pathway boxplot).

    `sample_groups` maps sample id -> group label, taken from an optional metadata
    table so the lightweight route can also show a group comparison, not just the
    heatmap/barplot summaries.
    """
    style()
    samples = [k for k in pathway_rows[0].keys() if k != "pathway"]
    known_samples = [s for s in samples if s in sample_groups]
    if not known_samples:
        return None
    groups = sorted({sample_groups[s] for s in known_samples})
    color_map = {g: GROUP_COLORS[i % len(GROUP_COLORS)] for i, g in enumerate(groups)}
    pathways = [r["pathway"] for r in pathway_rows]
    ncols = min(3, len(pathways))
    nrows = -(-len(pathways) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.6 * nrows), squeeze=False)
    for idx, row in enumerate(pathway_rows):
        ax = axes[idx // ncols][idx % ncols]
        data = [[float(row[s]) for s in known_samples if sample_groups[s] == g] for g in groups]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=True)
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups)
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(color_map[g])
            patch.set_alpha(0.75)
        ax.set_title(row["pathway"], fontsize=14, color=TEXT, wrap=True)
        ax.tick_params(axis="x", labelsize=12)
        ax.grid(True, axis="y", color=GRID)
    for idx in range(len(pathways), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig.suptitle("Pathway activity by group", fontsize=22, fontweight="bold", color=TEXT, y=1.02)
    fig.tight_layout()
    save(fig, out / "pathway_activity_boxplot")
    return groups


def pathway_group_summary(pathway_rows, sample_groups):
    """Per-pathway group means + a rank-based group-difference test (Mann-Whitney for
    2 groups, Kruskal-Wallis for >2), BH-adjusted across pathways. Requires scipy.
    """
    samples = [k for k in pathway_rows[0].keys() if k != "pathway"]
    known_samples = [s for s in samples if s in sample_groups]
    if not known_samples or scipy_stats is None:
        return None
    groups = sorted({sample_groups[s] for s in known_samples})
    rows_out = []
    pvals = []
    for row in pathway_rows:
        by_group = {g: [float(row[s]) for s in known_samples if sample_groups[s] == g] for g in groups}
        entry = {"pathway": row["pathway"]}
        for g in groups:
            vals = by_group[g]
            entry[f"mean_{g}"] = sum(vals) / len(vals) if vals else float("nan")
        if len(groups) == 2 and all(len(by_group[g]) >= 1 for g in groups):
            stat, p = scipy_stats.mannwhitneyu(by_group[groups[0]], by_group[groups[1]], alternative="two-sided")
            entry["test"] = "mann_whitney_u"
        elif len(groups) > 2 and all(len(by_group[g]) >= 1 for g in groups):
            stat, p = scipy_stats.kruskal(*[by_group[g] for g in groups])
            entry["test"] = "kruskal_wallis"
        else:
            stat, p = float("nan"), float("nan")
            entry["test"] = "not_applicable"
        entry["statistic"] = stat
        entry["pvalue"] = p
        rows_out.append(entry)
        pvals.append(p if p == p else 1.0)  # NaN-safe for BH ranking
    padj = bh_adjust(pvals)
    for entry, p in zip(rows_out, padj):
        entry["padj"] = p
    return rows_out


def top_correlations(module_rows, n=3):
    """Return up to `n` (module, trait, correlation, padj) tuples ranked by |correlation|,
    skipping rows with no usable correlation value."""
    parsed = []
    for r in module_rows:
        c = to_float(r.get("correlation"))
        if c is None:
            continue
        padj_raw = r.get("padj")
        p = to_float(padj_raw) if padj_raw not in (None, "") else None
        parsed.append((r["module"], r["trait"], c, p))
    parsed.sort(key=lambda x: abs(x[2]), reverse=True)
    return parsed[:n]


def fdr_survival_count(module_rows, padj_threshold=DEFAULT_PADJ_THRESHOLD):
    """Count how many module-trait cells have a usable padj, and how many survive
    < padj_threshold."""
    total, survived = 0, 0
    for r in module_rows:
        padj_raw = r.get("padj")
        p = to_float(padj_raw) if padj_raw not in (None, "") else None
        if p is None:
            continue
        total += 1
        if p < padj_threshold:
            survived += 1
    return survived, total


def top_pathway_diffs(group_summary_rows, n=3):
    """Return up to `n` group-summary rows ranked by padj (falls back to raw pvalue),
    skipping rows where no test could be run."""
    if not group_summary_rows:
        return []
    valid = [r for r in group_summary_rows if r.get("test") != "not_applicable" and r.get("pvalue") == r.get("pvalue")]
    valid.sort(key=lambda r: r.get("padj") if r.get("padj") == r.get("padj") else 1.0)
    return valid[:n]


def pathway_fdr_survival_count(group_summary_rows, padj_threshold=DEFAULT_PADJ_THRESHOLD):
    """Count how many pathway group-comparison rows have a usable padj, and how many
    survive < padj_threshold (mirrors fdr_survival_count for the module-trait matrix)."""
    if not group_summary_rows:
        return 0, 0
    total, survived = 0, 0
    for r in group_summary_rows:
        p = r.get("padj")
        if p is None or p != p:  # None or NaN
            continue
        total += 1
        if p < padj_threshold:
            survived += 1
    return survived, total


def infer_sample_size(module_rows, pathway_rows):
    """Best-effort sample count for the small-n WGCNA-reliability caveat: prefer an
    explicit n_samples column (e.g. written by run_wgcna_gsva.R), else fall back to the
    number of sample columns in the pathway score table."""
    if module_rows and "n_samples" in module_rows[0]:
        n = to_float(module_rows[0]["n_samples"])
        if n is not None:
            return int(n)
    if pathway_rows:
        return len([k for k in pathway_rows[0].keys() if k != "pathway"])
    return None


def top_hub_genes(hub_path, top_module, n=5):
    """Read an optional hub_genes.tsv (e.g. from run_wgcna_gsva.R) and return up to `n`
    gene names for `top_module`, preferring rows flagged is_hub, else ranked by
    kME_own_module/GS. Returns [] if the file/module/columns are not usable."""
    if not hub_path or not top_module:
        return []
    try:
        rows = read_csv(hub_path)
    except SystemExit:
        return []
    if "gene" not in rows[0] or "module" not in rows[0]:
        return []
    module_key = top_module[2:] if top_module.startswith("ME") else top_module
    candidates = [r for r in rows if r["module"] in (top_module, module_key)]
    if not candidates:
        return []
    if "is_hub" in candidates[0]:
        flagged = [r for r in candidates if str(r.get("is_hub")).strip().lower() == "true"]
        if flagged:
            candidates = flagged
    rank_col = "kME_own_module" if "kME_own_module" in candidates[0] else ("GS" if "GS" in candidates[0] else None)
    if rank_col:
        candidates.sort(key=lambda r: abs(to_float(r.get(rank_col)) or 0), reverse=True)
    return [r["gene"] for r in candidates[:n]]


def build_report(module_rows, pathway_rows, module_meta, groups_used, group_summary_rows, hub_genes, top_module_name, padj_threshold=DEFAULT_PADJ_THRESHOLD):
    n_modules = len(set(r["module"] for r in module_rows))
    n_traits = len(set(r["trait"] for r in module_rows))
    n_pathways = len(pathway_rows)
    n_samples = infer_sample_size(module_rows, pathway_rows)
    top_corrs = top_correlations(module_rows, n=3)
    survived, tested = fdr_survival_count(module_rows, padj_threshold)
    top_diffs = top_pathway_diffs(group_summary_rows, n=3) if group_summary_rows else []

    lines = ["# Network pathway quicklook", ""]
    lines.append("Generated module-trait and pathway activity figures from existing tables (lightweight/no-WGCNA route).")
    lines.append(f"Input: {n_modules} module(s) x {n_traits} trait(s); {n_pathways} pathway(s) scored"
                 + (f" across {n_samples} sample(s)." if n_samples else "."))
    if module_meta.get("n_na_cells"):
        lines.append(f"Note: {module_meta['n_na_cells']} of {module_meta['n_total_cells']} module x trait cell(s) had no usable correlation value and are shown as blank/NA.")

    lines.append("")
    lines.append("## Module-trait results")
    if top_corrs:
        for module, trait, corr, padj in top_corrs:
            sig = significance_stars(padj, padj_threshold) if padj is not None else ""
            p_txt = f", padj={padj:.3g}{sig}" if padj is not None else " (no p-value supplied)"
            lines.append(f"- **{module}** vs **{trait}**: r={corr:.2f}{p_txt}")
    else:
        lines.append("- No usable correlation values found in the module-trait table.")

    lines.append("")
    lines.append("## Pathway activity results")
    if top_diffs:
        lines.append(f"Group comparison across: {', '.join(groups_used)}.")
        for r in top_diffs:
            lines.append(f"- **{r['pathway']}**: {r['test']} statistic={r['statistic']:.3g}, p={r['pvalue']:.3g}, padj={r['padj']:.3g}{significance_stars(r['padj'], padj_threshold)}")
    elif groups_used:
        lines.append("Group summary computed but no pathway had a valid test result.")
    elif module_meta.get("metadata_no_overlap"):
        lines.append("Metadata provided but no sample overlap found with pathway score columns; boxplot/group summary skipped.")
    else:
        lines.append("No group metadata supplied; reporting overall pathway activity only (see barplot/heatmap).")

    if hub_genes:
        lines.append("")
        lines.append("## Hub genes")
        lines.append(f"Top module **{top_module_name}** hub gene candidate(s) (topological centrality, from supplied hub gene table): {', '.join(hub_genes)}.")

    lines.append("")
    lines.append("## Caveats")
    caveats = []
    if n_samples is not None and n_samples < SAMPLE_SIZE_CAUTION:
        tier = "strongly exploratory" if n_samples < SAMPLE_SIZE_EXPLORATORY else "exploratory"
        caveats.append(f"Sample size (n={n_samples}) is below the ~15-20 samples recommended for stable WGCNA module detection; treat module-trait associations as {tier} and avoid strong claims if these modules came from a formal WGCNA run (pitfalls.md #2).")
    if tested > 0:
        caveats.append(f"BH-FDR multiple-testing correction was applied across the full module x trait matrix ({tested} cell(s) with a usable p-value); {survived}/{tested} remained significant at padj<{padj_threshold:g} (pitfalls.md #3).")
    else:
        caveats.append("No usable p-values were supplied for the module-trait table; correlations are reported without a significance test and should be treated as descriptive only.")
    pw_survived, pw_tested = pathway_fdr_survival_count(group_summary_rows, padj_threshold)
    if pw_tested > 0:
        caveats.append(f"Pathway group-comparison p-values were BH-FDR corrected across all scored pathways ({pw_tested} tested); {pw_survived}/{pw_tested} remained significant at padj<{padj_threshold:g}.")
    caveats.append("Pathway activity scores are only comparable within this single scoring run; do not directly compare these values to scores computed in a separate GSVA/ssGSEA invocation on a different sample set (pitfalls.md #5).")
    if hub_genes:
        caveats.append("Hub gene status reflects intramodular connectivity/correlation strength only, not evidence of a causal or regulatory role (pitfalls.md #4).")
    lines.extend(f"- {c}" for c in caveats)

    lines.append("")
    lines.append("## Conclusion")
    conclusion_parts = []
    if top_corrs:
        m, t, c, p = top_corrs[0]
        direction = "positively" if c > 0 else "negatively"
        conclusion_parts.append(f"Module {m} is most strongly associated with {t} (r={c:.2f}, {direction} correlated" + (f", padj={p:.3g}" if p is not None else "") + ")")
    if top_diffs:
        d = top_diffs[0]
        conclusion_parts.append(f"pathway activity differs most between groups for {d['pathway']} ({d['test']}, p={d['pvalue']:.3g})")
    if hub_genes:
        conclusion_parts.append(f"candidate hub gene(s) for that module include {', '.join(hub_genes[:3])}")
    if conclusion_parts:
        lines.append("; ".join(conclusion_parts) + ".")
    else:
        lines.append("No module-trait or pathway-group signal could be summarized from the supplied tables.")
    lines.append("These results are correlational and exploratory: they describe co-expression and pathway-activity patterns in this dataset, not causal or mechanistic relationships. Validate hub genes and pathway differences in an independent cohort before drawing biological conclusions (pitfalls.md #4, #10).")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Plot module-trait and pathway activity tables.")
    parser.add_argument("--module-trait", required=True)
    parser.add_argument("--pathway-scores", required=True)
    parser.add_argument("--metadata", help="optional sample,group[,...] table; enables pathway boxplot + group summary")
    parser.add_argument("--group-col", default="group", help="column in --metadata holding the group label (default: group)")
    parser.add_argument("--hub-genes", help="optional hub_genes.tsv (e.g. from run_wgcna_gsva.R) to name hub genes for the top module in report.md")
    parser.add_argument("--padj-threshold", type=float, default=DEFAULT_PADJ_THRESHOLD,
                         help=f"BH-FDR adjusted-p significance cutoff for heatmap stars, pathway group-comparison flagging, and report.md text (default: {DEFAULT_PADJ_THRESHOLD})")
    parser.add_argument("--out", default="output")
    args = parser.parse_args()
    if not (0 < args.padj_threshold < 1):
        raise SystemExit(f"--padj-threshold must be between 0 and 1 (got {args.padj_threshold}).")
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    module_rows = read_csv(args.module_trait)
    pathway_rows = read_csv(args.pathway_scores)
    module_rows, module_meta = plot_module_trait(module_rows, out / "figures", args.padj_threshold)
    plot_pathways(pathway_rows, out / "figures")

    group_summary_rows = None
    groups_used = None
    if args.metadata:
        metadata_rows = read_csv(args.metadata)
        require_columns(metadata_rows, [("sample",)], args.metadata, "metadata table")
        if not all(args.group_col in r for r in metadata_rows):
            raise SystemExit(f"--group-col '{args.group_col}' not found in {args.metadata}")
        sample_groups = {r["sample"]: r[args.group_col] for r in metadata_rows}
        groups_used = plot_pathway_boxplot(pathway_rows, sample_groups, out / "figures")
        group_summary_rows = pathway_group_summary(pathway_rows, sample_groups)
        if group_summary_rows:
            fields = list(group_summary_rows[0].keys())
            lines = ["\t".join(fields)] + ["\t".join(str(r[f]) for f in fields) for r in group_summary_rows]
            (out / "results" / "pathway_group_summary.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
        else:
            module_meta["metadata_no_overlap"] = True

    top_module_name = top_correlations(module_rows, n=1)[0][0] if top_correlations(module_rows, n=1) else None
    hub_genes = top_hub_genes(args.hub_genes, top_module_name)

    module_fields = list(module_rows[0].keys())
    module_lines = ["\t".join(module_fields)] + ["\t".join(str(r.get(f, "")) for f in module_fields) for r in module_rows]
    (out / "results" / "module_trait_correlation.tsv").write_text("\n".join(module_lines) + "\n", encoding="utf-8")
    (out / "results" / "pathway_scores.tsv").write_text(Path(args.pathway_scores).read_text(encoding="utf-8"), encoding="utf-8")
    (out / "summary.tsv").write_text(f"item\tvalue\nmodules\t{len(set(r['module'] for r in module_rows))}\npathways\t{len(pathway_rows)}\n", encoding="utf-8")
    (out / "run" / "qc.json").write_text(json.dumps({
        "module_trait_rows": len(module_rows),
        "pathways": len(pathway_rows),
        "padj_method": "benjamini_hochberg",
        "padj_threshold": args.padj_threshold,
        "groups_compared": groups_used,
        "module_trait_na_cells": module_meta.get("n_na_cells", 0),
    }, indent=2), encoding="utf-8")

    report_text = build_report(module_rows, pathway_rows, module_meta, groups_used, group_summary_rows, hub_genes, top_module_name, args.padj_threshold)
    (out / "report.md").write_text(report_text, encoding="utf-8")
    print(f"wrote network/pathway outputs to {out}")


if __name__ == "__main__":
    main()
