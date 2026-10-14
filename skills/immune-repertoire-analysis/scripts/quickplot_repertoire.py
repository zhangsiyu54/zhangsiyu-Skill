#!/usr/bin/env python3
"""Quick-plot immune repertoire figures/tables from a clonotype table.

This is the primary analysis script for this skill when no R environment is
available (see scripts/check_environment.py). It implements, with only
matplotlib/numpy-grade dependencies:

  - clonal expansion by clone-size class (not a per-clone label dump)
  - top-N clonotype frequency / simple cross-sample clone tracking
  - diversity metrics (Shannon, Gini-Simpson, inverse Simpson, Chao1,
    clonality) with an explicit, non-silent depth-normalization caveat and an
    optional rarefaction (subsampling) estimate at a common depth
  - V and J gene usage (horizontal bars, sorted by count)
  - CDR3 length distribution (if a cdr3 column is present)
  - pairwise clonotype overlap between samples (Morisita-Horn + Jaccard)

See references/cns-method-guide.md (Required Figures, Algorithm Decision
Tree) and references/pitfalls.md (especially pitfalls 1, 4, 5, 6, 10) for the
rationale behind each design choice below.
"""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
except Exception as exc:
    raise SystemExit("Install matplotlib/numpy: python3 -m pip install matplotlib numpy") from exc

TEXT = "#173B5F"
GRID = "#DCEAF5"
BLUE = "#4EA3F1"
RED = "#E45B5B"
GREEN = "#4CB984"
AMBER = "#F0A202"
LIGHT_BLUE = "#B9D7EE"
OTHER_GRAY = "#C9D6E3"

TREAT_LIKE_GROUPS = {"treat", "treatment", "case", "tumor", "disease", "expanded", "post"}

# Clone-size-class thresholds as a fraction of each sample's total count.
# This is the immunarch-style "clonal space homeostasis" convention. The
# thresholds are dataset-scale dependent (they assume the whole repertoire,
# not just a handful of sorted cells, was sampled) -- adjust for repertoires
# with very different sequencing depth or cell recovery.
CLONE_SIZE_BINS = [
    (0.0, 0.01, "Rare"),
    (0.01, 0.05, "Small"),
    (0.05, 0.10, "Medium"),
    (0.10, 0.20, "Large"),
    (0.20, 1.000001, "Hyperexpanded"),
]
CLASS_COLORS = {
    "Rare": LIGHT_BLUE,
    "Small": BLUE,
    "Medium": GREEN,
    "Large": AMBER,
    "Hyperexpanded": RED,
}

DEFAULT_TOP_N_CLONES = 8
DEFAULT_RAREFACTION_ITERS = 200

# Trigger threshold for pitfalls.md #9 ("few hundred cells or fewer" produce
# unstable, falsely-precise diversity/overlap metrics).
SMALL_SAMPLE_DEPTH_THRESHOLD = 300
# Trigger threshold for flagging thin group-level / overlap comparisons.
MIN_SAMPLES_FOR_ROBUST_COMPARISON = 3

DEPTH_CAVEAT = (
    "Diversity indices (Shannon, Simpson, Chao1, clonality) are sensitive to "
    "sequencing depth / recovered cell number (references/pitfalls.md #1, #4): "
    "richness-sensitive metrics rise with sampling effort and clonality "
    "inflates in undersampled repertoires. Raw per-sample values below are "
    "NOT directly comparable across samples with different total counts. "
    "A rarefied Shannon estimate at the common minimum depth is reported "
    "alongside the raw value for a fairer comparison; treat samples whose "
    "rarefaction could not be computed (too few reads/cells) as descriptive "
    "only, and prefer the rarefied column when ranking samples."
)

REQUIRED_COLUMNS = ("sample", "clonotype")
COUNT_COLUMN_ALIASES = ("count", "frequency", "freq", "reads", "umis")


def read(path):
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"Clonotype table not found: {path}")
    if p.is_dir():
        raise SystemExit(f"Expected a clonotype table file but got a directory: {path}")
    try:
        text = p.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as exc:
        raise SystemExit(f"Could not read {path} as UTF-8 text: {exc}") from exc
    lines = [ln for ln in text.splitlines() if ln.strip() != ""]
    if not lines:
        return []
    rows = list(csv.DictReader(lines))
    # csv.DictReader silently drops rows that are pure whitespace/short; make
    # sure header parsing itself did not fail (e.g. no header row at all).
    if rows and rows[0] is None:
        raise SystemExit(f"Could not parse a header row from {path}.")
    return rows


def validate_and_normalize_columns(rows, path):
    """Fail fast with a clear message instead of a bare KeyError deep in the
    pipeline when required columns are missing (see references/input-contract.md).
    Also accepts 'frequency'/'freq'/'reads'/'umis' as an alias for 'count' and
    rescales fractional (0-1) frequency values to integer pseudo-counts so
    downstream integer-based logic (rarefaction, Chao1 singleton/doubleton
    counts) still works. Returns a list of human-readable notes to surface in
    report.md.
    """
    if not rows:
        raise SystemExit(f"No data rows found in {path}.")
    fieldnames = set(rows[0].keys())
    notes = []

    missing_required = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
    if missing_required:
        raise SystemExit(
            f"Missing required column(s) {missing_required} in {path}. "
            f"Expected at least: sample, clonotype, count (see references/input-contract.md). "
            f"Found columns: {sorted(fieldnames)}"
        )

    if "count" not in fieldnames:
        alias = next((c for c in COUNT_COLUMN_ALIASES if c in fieldnames and c != "count"), None)
        if alias is None:
            raise SystemExit(
                f"Missing a clonotype count/frequency column in {path}. "
                f"Expected one of {COUNT_COLUMN_ALIASES} (see references/input-contract.md). "
                f"Found columns: {sorted(fieldnames)}"
            )
        try:
            values = [float(row[alias]) for row in rows if row.get(alias, "") != ""]
        except ValueError as exc:
            raise SystemExit(f"Column '{alias}' in {path} contains a non-numeric value: {exc}") from exc
        fractional = bool(values) and max(values) <= 1.0001
        if fractional:
            scale = 1_000_000
            for row in rows:
                raw = row.get(alias, "")
                row["count"] = str(max(1, round(float(raw) * scale)) if raw not in ("", None) and float(raw) > 0 else 0)
            notes.append(
                f"No 'count' column was found; the '{alias}' column looked like a fractional frequency "
                f"(max value <= 1) and was rescaled by x{scale:,} to integer pseudo-counts for internal "
                "rarefaction/Chao1 computation. Treat absolute count-based numbers (e.g. top_clone_count) "
                "as relative only; frequencies/proportions remain valid."
            )
        else:
            for row in rows:
                row["count"] = row.get(alias, "0")
            notes.append(f"No 'count' column was found; used '{alias}' as the clonotype count column.")
    return notes


def to_count(value, context=""):
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        raise SystemExit(
            f"Could not parse clonotype count value {value!r}{context} as a number. "
            "Check the 'count' (or count-alias) column for blanks/non-numeric entries."
        )


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


def sample_group_map(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["sample"], row.get("group", ""))
    return groups


def sample_color(sample, group):
    if group and group.strip().lower() in TREAT_LIKE_GROUPS:
        return RED
    if not group and sample.upper().startswith("T"):
        return RED
    return BLUE


# ---------------------------------------------------------------------------
# Diversity metrics
# ---------------------------------------------------------------------------


def shannon(vals):
    """Shannon entropy on clonotype counts.

    NOTE: sensitive to sequencing depth / cell number (pitfalls.md #1) --
    do not compare Shannon values across samples with very different total
    counts without also consulting the rarefied estimate below.
    """
    total = sum(vals)
    ps = [v / total for v in vals if total and v > 0]
    return -sum(p * math.log(p) for p in ps)


def gini_simpson(vals):
    """1 - sum(p_i^2); higher = more diverse. Less depth-sensitive than
    richness-based metrics but still biased low-N (pitfalls.md #1, #9)."""
    total = sum(vals)
    if not total:
        return 0.0
    return 1 - sum((v / total) ** 2 for v in vals)


def inverse_simpson(vals):
    total = sum(vals)
    if not total:
        return 0.0
    denom = sum((v / total) ** 2 for v in vals)
    return (1 / denom) if denom else 0.0


def chao1(vals):
    """Bias-corrected Chao1 richness estimator from clonotype abundance
    counts (observed richness + correction from singleton/doubleton ratio)."""
    s_obs = sum(1 for v in vals if v > 0)
    f1 = sum(1 for v in vals if v == 1)
    f2 = sum(1 for v in vals if v == 2)
    if f2 > 0:
        return s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    if f1 > 1:
        return s_obs + (f1 * (f1 - 1)) / 2
    return float(s_obs)


def clonality(vals):
    """1 - normalized Shannon entropy (0 = maximally diverse, 1 = monoclonal)."""
    s_obs = sum(1 for v in vals if v > 0)
    if s_obs <= 1:
        return 1.0
    return 1 - shannon(vals) / math.log(s_obs)


def rarefied_shannon(vals, depth, iterations=DEFAULT_RAREFACTION_ITERS, seed=13):
    """Subsample without replacement to a common depth and average Shannon
    entropy over repeated draws (multivariate hypergeometric resampling).
    This is the "at minimum" fix from pitfalls.md #1 for comparing diversity
    across samples of unequal depth. Returns None if the sample has fewer
    than `depth` counted units.

    Uses numpy's multivariate_hypergeometric sampler rather than building an
    explicit per-unit pool: a naive pool of `sum(vals)` Python objects would
    scale with raw read/UMI count and can become extremely slow or exhaust
    memory for bulk repertoire data with counts in the hundreds of thousands
    to millions (or after alias-column rescaling of fractional frequencies).
    """
    if depth <= 0:
        return None
    colors = np.array([v for v in vals if v > 0], dtype=np.int64)
    if colors.size == 0 or depth > colors.sum():
        return None
    rng = np.random.default_rng(seed)
    entropies = []
    for _ in range(iterations):
        draw = rng.multivariate_hypergeometric(colors, int(depth))
        nz = draw[draw > 0]
        if nz.size == 0:
            entropies.append(0.0)
            continue
        p = nz / nz.sum()
        entropies.append(float(-(p * np.log(p)).sum()))
    return sum(entropies) / len(entropies)


def size_class(fraction):
    for lo, hi, label in CLONE_SIZE_BINS:
        if lo <= fraction < hi:
            return label
    return CLONE_SIZE_BINS[-1][2]


# ---------------------------------------------------------------------------
# Overlap metrics
# ---------------------------------------------------------------------------


def morisita_horn(x, y):
    """Morisita-Horn overlap index between two aligned clonotype count
    vectors (0-filled for clonotypes absent in one sample). Less sensitive
    to sample-size differences than a raw shared-clone count (pitfalls.md
    #4, #6); 0 = no overlap, 1 = identical clonotype-frequency structure."""
    x_total, y_total = sum(x), sum(y)
    if x_total == 0 or y_total == 0:
        return 0.0
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi ** 2 for xi in x)
    sum_y2 = sum(yi ** 2 for yi in y)
    denom = (sum_x2 / x_total ** 2 + sum_y2 / y_total ** 2) * x_total * y_total
    return (2 * sum_xy / denom) if denom else 0.0


def jaccard(set_a, set_b):
    """Presence/absence clonotype-ID-level overlap (not sequence-level)."""
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def compute_overlap_matrices(per_sample_counts, samples):
    all_clonotypes = sorted({c for counts in per_sample_counts.values() for c in counts})
    mh = [[0.0] * len(samples) for _ in samples]
    jac = [[0.0] * len(samples) for _ in samples]
    vectors = {s: [per_sample_counts[s].get(c, 0) for c in all_clonotypes] for s in samples}
    id_sets = {s: set(per_sample_counts[s]) for s in samples}
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i == j:
                mh[i][j] = 1.0
                jac[i][j] = 1.0
            else:
                mh[i][j] = morisita_horn(vectors[s1], vectors[s2])
                jac[i][j] = jaccard(id_sets[s1], id_sets[s2])
    return mh, jac


def write_matrix_tsv(matrix, samples, path):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + samples)
        for sample, row in zip(samples, matrix):
            writer.writerow([sample] + [round(v, 4) for v in row])


def plot_overlap_heatmap(matrix, samples, out_path, title):
    cmap = LinearSegmentedColormap.from_list("overlap", ["#FFFFFF", BLUE])
    fig, ax = plt.subplots(figsize=(max(6.5, len(samples) * 0.9), max(6.0, len(samples) * 0.9)))
    im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(samples)))
    ax.set_yticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_yticklabels(samples)
    for i in range(len(samples)):
        for j in range(len(samples)):
            val = matrix[i][j]
            color = "white" if val > 0.55 else TEXT
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=11)
    ax.set_title(title, loc="center", pad=18, color=TEXT, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=13)
    fig.tight_layout()
    save(fig, out_path)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_clone_size_classes(by_sample, samples, colors, out_path):
    """Clonal expansion by size class -- explicitly avoids labeling/plotting
    every individual clone (references/visual-style.md: "clonal expansion
    should not label all clones")."""
    class_labels = [c[2] for c in CLONE_SIZE_BINS]
    fractions = {label: [] for label in class_labels}
    for s in samples:
        vals = by_sample[s]
        total = sum(vals)
        class_totals = {label: 0 for label in class_labels}
        for v in vals:
            frac = (v / total) if total else 0
            class_totals[size_class(frac)] += v
        for label in class_labels:
            fractions[label].append((class_totals[label] / total) if total else 0)
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    bottom = [0.0] * len(samples)
    for label in class_labels:
        ax.bar(samples, fractions[label], bottom=bottom, label=label, color=CLASS_COLORS[label], alpha=0.92)
        bottom = [b + f for b, f in zip(bottom, fractions[label])]
    ax.set_title("Clonal expansion by clone-size class", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Sample", labelpad=12)
    ax.set_ylabel("Fraction of sample total count", labelpad=12)
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylim(0, 1.02)
    ax.grid(True, axis="y", color=GRID)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, title="Clone size")
    save(fig, out_path)
    del colors  # colors reserved for future per-group faceting; unused here on purpose


def compute_clonotype_totals(rows):
    """Single source of truth for per-clonotype totals used by both the
    top-clonotype plot and report.md, so the two never disagree."""
    total_by_clonotype = defaultdict(int)
    per_sample_clonotype = defaultdict(lambda: defaultdict(int))
    sample_total = defaultdict(int)
    for row in rows:
        c = to_count(row.get("count"), context=f" (sample={row.get('sample')!r}, clonotype={row.get('clonotype')!r})")
        total_by_clonotype[row["clonotype"]] += c
        per_sample_clonotype[row["sample"]][row["clonotype"]] += c
        sample_total[row["sample"]] += c
    return total_by_clonotype, per_sample_clonotype, sample_total


def plot_top_clonotype_frequency(
    total_by_clonotype, per_sample_clonotype, sample_total, samples, out_path, top_n=DEFAULT_TOP_N_CLONES
):
    """Top-N clonotype frequency, stacked with an explicit "Other" bucket.

    Only the top N clonotypes across the whole dataset are drawn/labeled
    (never all clones), and tracking the same colored band across samples
    gives a lightweight, dependency-free stand-in for alluvial/Sankey clone
    tracking (references/cns-method-guide.md Required Figures; pitfalls.md
    #10 on not over-interpreting such tracking as causal).
    """
    top_ids = [cid for cid, _ in sorted(total_by_clonotype.items(), key=lambda kv: -kv[1])[:top_n]]
    cmap = plt.get_cmap("plasma")
    colors = [cmap(i / max(1, len(top_ids) - 1)) for i in range(len(top_ids))] if top_ids else []
    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    bottom = [0.0] * len(samples)
    for cid, color in zip(top_ids, colors):
        fracs = [
            (per_sample_clonotype[s].get(cid, 0) / sample_total[s]) if sample_total[s] else 0.0
            for s in samples
        ]
        ax.bar(samples, fracs, bottom=bottom, label=cid, color=color)
        bottom = [b + f for b, f in zip(bottom, fracs)]
    other_frac = [max(0.0, 1.0 - b) for b in bottom]
    ax.bar(samples, other_frac, bottom=bottom, label="Other clones", color=OTHER_GRAY)
    ax.set_title(f"Top {len(top_ids)} clonotype frequency (rest = Other)", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Sample", labelpad=12)
    ax.set_ylabel("Fraction of sample total", labelpad=12)
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylim(0, 1.02)
    ax.grid(True, axis="y", color=GRID)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=11, title="Clonotype (top N)")
    save(fig, out_path)
    return top_ids


def plot_diversity(samples, colors, raw_shannon, rarefied, totals, common_depth, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(16.2, 6.4))
    ax = axes[0]
    width = 0.38
    x = list(range(len(samples)))
    ax.bar([i - width / 2 for i in x], raw_shannon, width=width, color=colors, alpha=0.92, label="Shannon (raw)")
    rarefied_plot = [v if v is not None else 0.0 for v in rarefied]
    ax.bar(
        [i + width / 2 for i in x],
        rarefied_plot,
        width=width,
        color=colors,
        alpha=0.45,
        hatch="//",
        label=f"Shannon (rarefied, n={common_depth})" if common_depth else "Shannon (rarefied, n/a)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=35)
    ax.set_title("Clonal diversity: raw vs. depth-rarefied", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_ylabel("Shannon diversity", labelpad=12)
    ax.grid(True, axis="y", color=GRID)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    # Shannon diversity is non-negative; avoid matplotlib's symmetric
    # +/-0.05 autoscale fallback when every sample is monoclonal (all 0s).
    top = ax.get_ylim()[1]
    ax.set_ylim(0, max(top, 0.1))

    ax2 = axes[1]
    ax2.bar(samples, totals, color=colors, alpha=0.86)
    if common_depth:
        ax2.axhline(common_depth, color=TEXT, linestyle="--", linewidth=1.3)
        ax2.text(
            len(samples) - 0.5,
            common_depth,
            f" common depth = {common_depth}",
            va="bottom",
            ha="right",
            color=TEXT,
            fontsize=12,
        )
    ax2.set_title("Sample depth (total count)", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax2.set_xlabel("Sample", labelpad=12)
    ax2.set_ylabel("Total clonotype count", labelpad=12)
    ax2.tick_params(axis="x", rotation=35)
    ax2.grid(True, axis="y", color=GRID)

    fig.text(
        0.5,
        -0.04,
        "Diversity is depth-sensitive: compare the rarefied (right-hatched) bars, not raw Shannon,\n"
        "when samples differ in total count. See references/pitfalls.md #1.",
        ha="center",
        va="top",
        fontsize=12,
        color=TEXT,
    )
    fig.tight_layout()
    save(fig, out_path)


def gene_usage_plot(rows, gene_col, out_path, title):
    usage = defaultdict(int)
    for row in rows:
        gene = row.get(gene_col)
        if not gene:
            continue
        usage[gene] += to_count(row.get("count"), context=f" (sample={row.get('sample')!r}, gene={gene!r})")
    if not usage:
        return None
    items = sorted(usage.items(), key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(8.6, max(5.8, len(items) * 0.44 + 2.2)))
    ax.barh([x[0] for x in items], [x[1] for x in items], color=GREEN, alpha=0.88)
    ax.set_title(title, loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("Count", labelpad=12)
    ax.set_ylabel(gene_col.replace("_", " "), labelpad=12)
    ax.grid(True, axis="x", color=GRID)
    save(fig, out_path)
    return dict(items)


def plot_cdr3_length(rows, out_path):
    if not rows or "cdr3" not in rows[0] or not any(row.get("cdr3") for row in rows):
        return False
    lengths_by_group = defaultdict(list)
    for row in rows:
        cdr3 = row.get("cdr3")
        if not cdr3:
            continue
        length = len(cdr3)
        weight = to_count(row.get("count"), context=f" (sample={row.get('sample')!r}, cdr3={cdr3!r})")
        group = row.get("group") or "all"
        lengths_by_group[group].extend([length] * weight)
    all_lengths = [l for vals in lengths_by_group.values() for l in vals]
    if not all_lengths:
        return False
    bins = range(min(all_lengths), max(all_lengths) + 2)
    fig, ax = plt.subplots(figsize=(8.6, 6.0))
    for group, lengths in lengths_by_group.items():
        color = RED if group.strip().lower() in TREAT_LIKE_GROUPS else BLUE
        ax.hist(lengths, bins=bins, alpha=0.55, label=group, color=color, edgecolor="white")
    ax.set_title("CDR3 length distribution", loc="center", pad=18, color=TEXT, fontweight="bold")
    ax.set_xlabel("CDR3 length (aa)", labelpad=12)
    ax.set_ylabel("Count-weighted frequency", labelpad=12)
    ax.grid(True, axis="y", color=GRID)
    if len(lengths_by_group) > 1:
        ax.legend(frameon=False)
    save(fig, out_path)
    return True


# ---------------------------------------------------------------------------
# Caveat assessment (pitfalls.md-driven, data-checked -- only fires when the
# data actually shows the trigger condition)
# ---------------------------------------------------------------------------


def assess_clonotype_definition(rows):
    """Data-driven check for pitfalls.md #2/#5: does each clonotype ID
    correspond to exactly one CDR3(+V/J) identity, or is there evidence of
    merging (multiple distinct CDR3 sequences sharing one clonotype ID --
    amino-acid-level/ID-level grouping that convergent recombination can
    inflate into a false "clone") or splitting (one CDR3+V/J identity spread
    across multiple clonotype IDs -- possible nucleotide-level
    over-splitting of a truly clonal lineage)?

    Returns None if there is no cdr3 column to check against (cannot
    assess), otherwise a list of caveat strings (always non-empty: at
    minimum a reminder to state the clonotype definition explicitly).
    """
    if not any(row.get("cdr3") for row in rows):
        return None
    cdr3_by_clonotype = defaultdict(set)
    clonotype_by_identity = defaultdict(set)
    vj_by_cdr3 = defaultdict(set)
    for row in rows:
        cdr3 = row.get("cdr3")
        clonotype = row.get("clonotype")
        if not cdr3 or not clonotype:
            continue
        cdr3_by_clonotype[clonotype].add(cdr3)
        v_gene, j_gene = row.get("v_gene") or "", row.get("j_gene") or ""
        clonotype_by_identity[(cdr3, v_gene, j_gene)].add(clonotype)
        if v_gene or j_gene:
            vj_by_cdr3[cdr3].add((v_gene, j_gene))
    merged = {k: v for k, v in cdr3_by_clonotype.items() if len(v) > 1}
    split = {k: v for k, v in clonotype_by_identity.items() if len(v) > 1}
    convergent = {k: v for k, v in vj_by_cdr3.items() if len(v) > 1}
    notes = []
    if convergent:
        example_cdr3 = sorted(convergent)[0]
        example_vj = sorted(convergent[example_cdr3])
        notes.append(
            f"{len(convergent)} CDR3 sequence(s) (e.g. '{example_cdr3}') are found paired with more than one "
            f"distinct V/J gene combination (e.g. {example_vj}). This is the classic convergent-recombination "
            "signature: the same amino-acid CDR3 can arise from independent rearrangement events. Do not treat "
            "these as a single expanded clone based on CDR3 identity alone -- verify nucleotide-level identity "
            "or accept them as biologically distinct clones with a shared amino-acid solution "
            "(references/pitfalls.md #5, #6)."
        )
    if merged:
        example = sorted(merged)[0]
        notes.append(
            f"{len(merged)} clonotype ID(s) (e.g. '{example}') map to more than one distinct CDR3 "
            "sequence. If clonotyping was done at the amino-acid level, independently-generated "
            "convergent-recombination sequences may have been merged into one apparent 'clone', "
            "inflating expansion and deflating diversity (references/pitfalls.md #2, #5). Verify "
            "shared V/J usage and low nucleotide diversity before calling merged IDs a single true clone."
        )
    if split:
        example = sorted(split, key=lambda k: k[0])[0]
        notes.append(
            f"{len(split)} identical CDR3(+V/J) identity/identities (e.g. CDR3 '{example[0]}') appear "
            "under more than one distinct clonotype ID. If clonotyping was done at the nucleotide level, "
            "silent-mutation or allele-calling differences may have split one true clonal lineage into "
            "separate artificial singletons, deflating apparent expansion (references/pitfalls.md #2, #5)."
        )
    if not notes:
        notes.append(
            "No evidence of clonotype-ID merging/splitting was detected against CDR3(+V/J) identity in "
            "this table, but the clonotype definition (CDR3 nucleotide vs. amino acid vs. paired-chain) "
            "should still be confirmed and stated before treating cross-sample or cross-study comparisons "
            "as apples-to-apples (references/pitfalls.md #2)."
        )
    return notes


def assess_small_sample_sizes(samples, totals, groups):
    """pitfalls.md #9: flag samples with very low depth/cell recovery, and
    thin group/overlap comparisons where conclusions would be over-precise."""
    notes = []
    thin = [(s, t) for s, t in zip(samples, totals) if t < SMALL_SAMPLE_DEPTH_THRESHOLD]
    if thin:
        listing = ", ".join(f"{s} (n={t})" for s, t in thin)
        notes.append(
            f"{len(thin)} of {len(samples)} sample(s) have fewer than {SMALL_SAMPLE_DEPTH_THRESHOLD} total "
            f"counted units: {listing}. Diversity/clonality/overlap point estimates for these samples carry "
            "high sampling variance and should not be ranked or statistically compared without bootstrap or "
            "jackknife confidence intervals (references/pitfalls.md #9). Treat their values as descriptive only."
        )
    if len(samples) < MIN_SAMPLES_FOR_ROBUST_COMPARISON:
        notes.append(
            f"Only {len(samples)} sample(s) are present. Group-level diversity/overlap comparisons and any "
            "'group A vs. group B' language are not statistically supported with this few replicates; "
            "report per-sample values descriptively rather than as a tested group difference "
            "(references/pitfalls.md #4, #9)."
        )
    if groups:
        group_counts = Counter(g for g in groups.values() if g)
        thin_groups = [g for g, n in group_counts.items() if n < 2]
        if thin_groups and len(group_counts) > 1:
            notes.append(
                f"Group(s) {sorted(thin_groups)} have fewer than 2 samples. A group with a single replicate "
                "cannot support a variance estimate; do not report group-level significance for it."
            )
    return notes


# ---------------------------------------------------------------------------
# report.md
# ---------------------------------------------------------------------------


def write_report_md(
    out,
    input_path,
    rows,
    samples,
    groups,
    totals,
    n_clonotypes_obs,
    unique_clonotypes,
    common_depth,
    raw_shannon,
    rarefied,
    gini,
    inv_simpson,
    chao,
    clon,
    total_by_clonotype,
    per_sample_clonotype,
    sample_total,
    top_n_clones,
    v_usage,
    j_usage,
    mh_matrix,
    jac_matrix,
    column_notes,
    definition_notes,
    small_sample_notes,
):
    lines = []
    lines.append("# Immune Repertoire Analysis Report")
    lines.append("")
    lines.append(f"Input: `{input_path}`")
    lines.append("")

    # --- Overview -----------------------------------------------------
    lines.append("## Overview")
    lines.append("")
    group_list = sorted({g for g in groups.values() if g}) or ["(none provided)"]
    lines.append(f"- Samples: {len(samples)} ({', '.join(samples)})")
    lines.append(f"- Groups: {', '.join(group_list)}")
    lines.append(f"- Clonotype rows: {len(rows)}; unique clonotype IDs: {unique_clonotypes}")
    lines.append(
        f"- V gene column present: {'yes' if v_usage is not None else 'no'}; "
        f"J gene column present: {'yes' if j_usage is not None else 'no'}; "
        f"CDR3 column present: {'yes' if any(row.get('cdr3') for row in rows) else 'no'}"
    )
    lines.append("")

    # --- QC / per-sample summary ---------------------------------------
    lines.append("## QC / Per-Sample Summary")
    lines.append("")
    lines.append("| Sample | Group | Total count (depth) | Clonotypes observed | Top clone count |")
    lines.append("|---|---|---|---|---|")
    for i, s in enumerate(samples):
        top_count = max(
            (per_sample_clonotype[s].values() if per_sample_clonotype[s] else [0]),
            default=0,
        )
        lines.append(f"| {s} | {groups.get(s, '')} | {totals[i]} | {n_clonotypes_obs[i]} | {top_count} |")
    lines.append("")

    # --- Clonal expansion / top clones ---------------------------------
    lines.append("## Clonal Expansion: Top Clonotypes")
    lines.append("")
    grand_total = sum(totals) or 1
    top_sorted = sorted(total_by_clonotype.items(), key=lambda kv: -kv[1])[:top_n_clones]
    if top_sorted:
        lines.append("| Rank | Clonotype | Overall frequency | Peak sample | Peak sample frequency |")
        lines.append("|---|---|---|---|---|")
        for rank, (cid, total_c) in enumerate(top_sorted, start=1):
            best_sample, best_frac = "", 0.0
            for s in samples:
                c = per_sample_clonotype[s].get(cid, 0)
                frac = (c / sample_total[s]) if sample_total[s] else 0.0
                if frac > best_frac:
                    best_frac, best_sample = frac, s
            lines.append(
                f"| {rank} | `{cid}` | {total_c / grand_total:.2%} | {best_sample} | {best_frac:.2%} |"
            )
    else:
        lines.append("No clonotypes found.")
    lines.append("")
    lines.append("See `figures/clonotype_expansion.png` (clone-size-class composition per sample) and "
                 "`figures/top_clonotype_frequency.png` (named top clones tracked across samples).")
    lines.append("")

    # --- Diversity -------------------------------------------------------
    lines.append("## Diversity Comparison (Raw vs. Rarefied)")
    lines.append("")
    lines.append(
        "| Sample | Group | Depth | Shannon (raw) | Shannon (rarefied, n={}) | Gini-Simpson | Inverse Simpson | Chao1 | Clonality |".format(
            common_depth if common_depth else "n/a"
        )
    )
    lines.append("|---|---|---|---|---|---|---|---|---|")

    def z(v):
        """Avoid printing '-0.000' from floating-point negative zero."""
        return 0.0 if v == 0 else v

    for i, s in enumerate(samples):
        rare = f"{z(rarefied[i]):.3f}" if rarefied[i] is not None else "NA (below common depth)"
        lines.append(
            f"| {s} | {groups.get(s, '')} | {totals[i]} | {z(raw_shannon[i]):.3f} | {rare} | "
            f"{z(gini[i]):.3f} | {z(inv_simpson[i]):.3f} | {z(chao[i]):.2f} | {z(clon[i]):.3f} |"
        )
    lines.append("")
    if group_list != ["(none provided)"] and len(group_list) > 1:
        lines.append("Group-level mean rarefied Shannon (descriptive; see sample-size caveat below if applicable):")
        lines.append("")
        for g in group_list:
            vals = [rarefied[i] for i, s in enumerate(samples) if groups.get(s) == g and rarefied[i] is not None]
            if vals:
                lines.append(f"- {g}: mean rarefied Shannon = {sum(vals) / len(vals):.3f} (n={len(vals)} samples)")
            else:
                lines.append(f"- {g}: rarefied Shannon not available for any sample in this group")
        lines.append("")

    # --- V/J usage ---------------------------------------------------------
    if v_usage or j_usage:
        lines.append("## V/J Gene Usage Highlights")
        lines.append("")
        for label, usage in (("V gene", v_usage), ("J gene", j_usage)):
            if not usage:
                continue
            usage_total = sum(usage.values()) or 1
            top_genes = sorted(usage.items(), key=lambda kv: -kv[1])[:3]
            formatted = ", ".join(f"{g} ({c / usage_total:.1%})" for g, c in top_genes)
            lines.append(f"- Top {label}s by count: {formatted}")
        lines.append("")

    # --- Overlap -------------------------------------------------------
    if mh_matrix is not None:
        lines.append("## Clonotype Overlap")
        lines.append("")
        pairs = []
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                pairs.append((samples[i], samples[j], mh_matrix[i][j], jac_matrix[i][j]))
        pairs.sort(key=lambda x: -x[2])
        lines.append("| Sample A | Sample B | Morisita-Horn | Jaccard (ID-level) |")
        lines.append("|---|---|---|---|")
        for a, b, mh, jac in pairs[:10]:
            lines.append(f"| {a} | {b} | {mh:.3f} | {jac:.3f} |")
        if len(pairs) > 10:
            lines.append(f"\n_({len(pairs) - 10} additional pair(s) omitted; see `results/overlap_matrix_*.tsv` for the full matrix.)_")
        lines.append("")
        lines.append(
            "Shared clonotypes between samples can arise from convergent recombination rather than a "
            "shared antigen-driven response; treat overlap as biologically meaningful only after comparison "
            "to a permutation or generation-probability null model (references/pitfalls.md #6)."
        )
        lines.append("")
    else:
        lines.append("## Clonotype Overlap")
        lines.append("")
        lines.append("Skipped: fewer than 2 samples are present, so pairwise overlap is not defined.")
        lines.append("")

    # --- Caveats ---------------------------------------------------------
    lines.append("## Caveats")
    lines.append("")
    lines.append(f"- **Sequencing depth / cell number sensitivity** (references/pitfalls.md #1, #4): {DEPTH_CAVEAT}")
    for note in small_sample_notes:
        lines.append(f"- **Small sample size** (references/pitfalls.md #9): {note}")
    if definition_notes:
        for note in definition_notes:
            lines.append(f"- **Clonotype definition / convergent recombination** (references/pitfalls.md #2, #5): {note}")
    else:
        lines.append(
            "- **Clonotype definition**: no `cdr3` column was available, so merging/splitting could not be "
            "checked against sequence identity; confirm the clonotype definition (nucleotide vs. amino acid "
            "vs. paired-chain) with the data provider before reporting expansion/diversity numbers "
            "(references/pitfalls.md #2)."
        )
    for note in column_notes:
        lines.append(f"- **Input normalization**: {note}")
    lines.append("")

    # --- Conclusion --------------------------------------------------------
    lines.append("## Conclusion")
    lines.append("")
    if top_sorted:
        top_cid, top_total = top_sorted[0]
        top_frac = top_total / grand_total
    else:
        top_cid, top_frac = "n/a", 0.0
    mean_clonality = sum(clon) / len(clon) if clon else 0.0
    conclusion = (
        f"Across {len(samples)} sample(s) and {unique_clonotypes} unique clonotypes, the most expanded "
        f"clonotype overall is `{top_cid}` at {top_frac:.1%} of total counted units, and mean clonality "
        f"across samples is {mean_clonality:.2f} (0 = maximally diverse, 1 = monoclonal). "
    )
    if common_depth and any(v is not None for v in rarefied):
        conclusion += (
            "Diversity comparisons above use the depth-rarefied Shannon value where available; raw values "
            "are shown alongside for transparency but should not drive cross-sample ranking on their own. "
        )
    if mh_matrix is not None and pairs:
        top_pair = pairs[0]
        conclusion += (
            f"The most similar sample pair by clonotype-frequency structure is {top_pair[0]}/{top_pair[1]} "
            f"(Morisita-Horn = {top_pair[2]:.2f}); this reflects shared clonotype usage and does not by "
            "itself establish a shared antigen-specific response. "
        )
    conclusion += (
        "All expansion, diversity, and overlap statements above are descriptive summaries of the clonotype "
        "table as provided; treat them as hypothesis-generating rather than confirmatory until the caveats "
        "above are addressed (definition, depth normalization, sample size) and, where relevant, validated "
        "against antigen-specificity or lineage evidence."
    )
    lines.append(conclusion)
    lines.append("")

    lines.append("## File Index")
    lines.append("")
    lines.append("- `results/clonotype_summary.tsv`, `results/diversity_summary.tsv`, `results/diversity_caveat.txt`")
    if mh_matrix is not None:
        lines.append("- `results/overlap_matrix_morisita_horn.tsv`, `results/overlap_matrix_jaccard.tsv`")
    lines.append(
        "- `figures/clonotype_expansion.{png,pdf}`, `figures/top_clonotype_frequency.{png,pdf}`, "
        "`figures/diversity_barplot.{png,pdf}`"
    )
    if v_usage:
        lines.append("- `figures/v_gene_usage.{png,pdf}`")
    if j_usage:
        lines.append("- `figures/j_gene_usage.{png,pdf}`")
    if any(row.get("cdr3") for row in rows):
        lines.append("- `figures/cdr3_length_distribution.{png,pdf}`")
    if mh_matrix is not None:
        lines.append("- `figures/clonotype_overlap_heatmap.{png,pdf}`")
    lines.append("- `summary.tsv`, `run/qc.json`, `methods.md` (written by `scripts/write_methods_doc.py`)")
    lines.append("")

    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Immune repertoire quick plots from a clonotype table.")
    parser.add_argument("--clonotypes", required=True)
    parser.add_argument("--out", default="output")
    parser.add_argument("--top-n-clones", type=int, default=DEFAULT_TOP_N_CLONES)
    parser.add_argument("--rarefaction-iters", type=int, default=DEFAULT_RAREFACTION_ITERS)
    args = parser.parse_args()
    apply_style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)

    rows = read(args.clonotypes)
    if not rows:
        raise SystemExit(f"No rows found in clonotype table: {args.clonotypes}")
    column_notes = validate_and_normalize_columns(rows, args.clonotypes)

    groups = sample_group_map(rows)
    by_sample = defaultdict(list)
    per_sample_counts = defaultdict(lambda: defaultdict(int))
    v_use_present = any(row.get("v_gene") for row in rows)
    j_use_present = any(row.get("j_gene") for row in rows)
    cdr3_present = any(row.get("cdr3") for row in rows)

    for row in rows:
        c = to_count(row.get("count"), context=f" (sample={row.get('sample')!r}, clonotype={row.get('clonotype')!r})")
        by_sample[row["sample"]].append(c)
        per_sample_counts[row["sample"]][row["clonotype"]] += c

    samples = sorted(by_sample)
    colors = [sample_color(s, groups.get(s, "")) for s in samples]

    # 1. Clonal expansion by clone-size class (no per-clone labels).
    plot_clone_size_classes(by_sample, samples, colors, out / "figures" / "clonotype_expansion")

    # 2. Top-N clonotype frequency / lightweight clone tracking across samples.
    total_by_clonotype, per_sample_clonotype, sample_total = compute_clonotype_totals(rows)
    top_ids = plot_top_clonotype_frequency(
        total_by_clonotype,
        per_sample_clonotype,
        sample_total,
        samples,
        out / "figures" / "top_clonotype_frequency",
        top_n=args.top_n_clones,
    )

    # 3. Diversity: raw + rarefied Shannon, Gini-Simpson, inverse Simpson,
    #    Chao1, clonality -- always reported next to total count (depth).
    raw_shannon = [shannon(by_sample[s]) for s in samples]
    totals = [sum(by_sample[s]) for s in samples]
    common_depth = min(totals) if totals and min(totals) > 0 else 0
    rarefied = [
        rarefied_shannon(by_sample[s], common_depth, iterations=args.rarefaction_iters) if common_depth else None
        for s in samples
    ]
    gini = [gini_simpson(by_sample[s]) for s in samples]
    inv_simpson = [inverse_simpson(by_sample[s]) for s in samples]
    chao = [chao1(by_sample[s]) for s in samples]
    clon = [clonality(by_sample[s]) for s in samples]
    n_clonotypes_obs = [sum(1 for v in by_sample[s] if v > 0) for s in samples]

    plot_diversity(samples, colors, raw_shannon, rarefied, totals, common_depth, out / "figures" / "diversity_barplot")

    print("CAVEAT:", DEPTH_CAVEAT)
    (out / "results" / "diversity_caveat.txt").write_text(DEPTH_CAVEAT + "\n", encoding="utf-8")

    with (out / "results" / "diversity_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "group",
                "total_count_depth",
                "n_clonotypes_observed",
                "top_clone_count",
                "shannon_raw",
                "shannon_rarefied_common_depth",
                "common_depth_used",
                "gini_simpson",
                "inverse_simpson",
                "chao1",
                "clonality",
            ]
        )
        for i, s in enumerate(samples):
            writer.writerow(
                [
                    s,
                    groups.get(s, ""),
                    totals[i],
                    n_clonotypes_obs[i],
                    max(by_sample[s]) if by_sample[s] else 0,
                    round(raw_shannon[i], 4),
                    round(rarefied[i], 4) if rarefied[i] is not None else "NA",
                    common_depth,
                    round(gini[i], 4),
                    round(inv_simpson[i], 4),
                    round(chao[i], 4),
                    round(clon[i], 4),
                ]
            )

    # 4. V/J gene usage (horizontal bars, sorted by count).
    v_usage = None
    j_usage = None
    if v_use_present:
        v_usage = gene_usage_plot(rows, "v_gene", out / "figures" / "v_gene_usage", "V gene usage")
    else:
        print("Note: no v_gene column found; skipped V gene usage plot.")
    if j_use_present:
        j_usage = gene_usage_plot(rows, "j_gene", out / "figures" / "j_gene_usage", "J gene usage")
    else:
        print("Note: no j_gene column found; skipped J gene usage plot.")

    # 5. CDR3 length distribution.
    if cdr3_present:
        plot_cdr3_length(rows, out / "figures" / "cdr3_length_distribution")
    else:
        print("Note: no cdr3 column found; skipped CDR3 length distribution plot.")

    # 6. Pairwise clonotype overlap (only meaningful with 2+ samples).
    mh_matrix = None
    jac_matrix = None
    if len(samples) >= 2:
        mh_matrix, jac_matrix = compute_overlap_matrices(per_sample_counts, samples)
        plot_overlap_heatmap(
            mh_matrix, samples, out / "figures" / "clonotype_overlap_heatmap", "Clonotype overlap (Morisita-Horn)"
        )
        write_matrix_tsv(mh_matrix, samples, out / "results" / "overlap_matrix_morisita_horn.tsv")
        write_matrix_tsv(jac_matrix, samples, out / "results" / "overlap_matrix_jaccard.tsv")
    else:
        print("Note: fewer than 2 samples; skipped overlap heatmap.")

    # Raw table copy + top-level summary.
    (out / "results" / "clonotype_summary.tsv").write_text(
        Path(args.clonotypes).read_text(encoding="utf-8"), encoding="utf-8"
    )
    unique_clonotypes = len({row["clonotype"] for row in rows})
    (out / "summary.tsv").write_text(
        "item\tvalue\n"
        f"samples\t{len(samples)}\n"
        f"clonotype_rows\t{len(rows)}\n"
        f"unique_clonotypes\t{unique_clonotypes}\n"
        f"common_depth_for_rarefaction\t{common_depth}\n",
        encoding="utf-8",
    )
    definition_notes = assess_clonotype_definition(rows)
    small_sample_notes = assess_small_sample_sizes(samples, totals, groups)
    for note in small_sample_notes:
        print("CAVEAT (small sample size):", note)
    if definition_notes:
        for note in definition_notes:
            print("CAVEAT (clonotype definition):", note)

    (out / "run" / "qc.json").write_text(
        json.dumps(
            {
                "samples": len(samples),
                "rows": len(rows),
                "unique_clonotypes": unique_clonotypes,
                "has_v_gene": v_use_present,
                "has_j_gene": j_use_present,
                "has_cdr3": cdr3_present,
                "common_depth_for_rarefaction": common_depth,
                "top_clonotypes_tracked": top_ids,
                "depth_caveat": DEPTH_CAVEAT,
                "small_sample_caveats": small_sample_notes,
                "clonotype_definition_caveats": definition_notes,
                "input_column_notes": column_notes,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    write_report_md(
        out,
        args.clonotypes,
        rows,
        samples,
        groups,
        totals,
        n_clonotypes_obs,
        unique_clonotypes,
        common_depth,
        raw_shannon,
        rarefied,
        gini,
        inv_simpson,
        chao,
        clon,
        total_by_clonotype,
        per_sample_clonotype,
        sample_total,
        args.top_n_clones,
        v_usage,
        j_usage,
        mh_matrix,
        jac_matrix,
        column_notes,
        definition_notes,
        small_sample_notes,
    )
    print(f"wrote repertoire outputs to {out}")
    print(f"wrote {out / 'report.md'}")


if __name__ == "__main__":
    main()
