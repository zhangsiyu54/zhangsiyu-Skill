#!/usr/bin/env python3
import argparse
import csv
import json
import math
import sys
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except Exception as exc:
    raise SystemExit("Install matplotlib and numpy: python3 -m pip install matplotlib numpy") from exc


PALETTE = ["#4EA3F1", "#E45B5B", "#4CB984", "#E2A13B", "#8B72D6", "#56BFC0", "#F08AB4", "#8AA1B4", "#B8C2CC", "#D9C25C"]
TEXT = "#173B5F"
GRID = "#DCEAF5"


def read_csv(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 16, "axes.titlesize": 24, "axes.labelsize": 21, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15, "axes.edgecolor": "#B9D7EE", "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT, "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42})


def shannon(vals):
    total = sum(vals)
    ps = [v / total for v in vals if total and v > 0]
    return -sum(p * math.log(p) for p in ps)


def simpson(vals):
    total = sum(vals)
    if not total:
        return 0.0
    ps = [v / total for v in vals if v > 0]
    return 1.0 - sum(p * p for p in ps)


def bray_curtis_distance(matrix):
    """Pairwise Bray-Curtis distance between samples (columns) of a taxa x sample count matrix."""
    n = matrix.shape[1]
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(matrix[:, i] - matrix[:, j]).sum()
            denom = (matrix[:, i] + matrix[:, j]).sum()
            d = num / denom if denom > 0 else 0.0
            dist[i, j] = dist[j, i] = d
    return dist


def classical_mds(dist, k=2):
    """Classical (Torgerson) MDS / PCoA via double-centering, equivalent to R's stats::cmdscale."""
    n = dist.shape[0]
    d2 = dist ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j.dot(d2).dot(j)
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = np.clip(eigvals[:k], 0, None)
    coords = eigvecs[:, :k] * np.sqrt(pos)
    denom = eigvals[eigvals > 0].sum()
    var_explained = (eigvals[:k] / denom * 100) if denom > 0 else np.zeros(k)
    return coords, var_explained


def _ss_within(d2, groups):
    total = 0.0
    for g in set(groups):
        idx = np.where(groups == g)[0]
        n_g = len(idx)
        if n_g <= 1:
            continue
        sub = d2[np.ix_(idx, idx)]
        total += sub[np.triu_indices(n_g, k=1)].sum() / n_g
    return total


def permanova_test(dist, groups, n_perm=999, seed=0):
    """Distance-based one-way PERMANOVA (McArdle & Anderson 2001), equivalent in spirit to vegan::adonis2."""
    groups = np.asarray(groups)
    n = len(groups)
    a = len(set(groups))
    if n <= a:
        return {"F": float("nan"), "R2": float("nan"), "p_value": float("nan"), "permutations": 0}
    d2 = dist ** 2
    ss_t = d2[np.triu_indices(n, k=1)].sum() / n
    ss_w_obs = _ss_within(d2, groups)
    ss_a_obs = ss_t - ss_w_obs
    f_obs = (ss_a_obs / (a - 1)) / (ss_w_obs / (n - a)) if ss_w_obs > 0 else float("inf")
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        ss_w = _ss_within(d2, perm)
        ss_a = ss_t - ss_w
        f_perm = (ss_a / (a - 1)) / (ss_w / (n - a)) if ss_w > 0 else float("inf")
        if f_perm >= f_obs:
            count += 1
    p_value = (count + 1) / (n_perm + 1)
    r2 = ss_a_obs / ss_t if ss_t > 0 else float("nan")
    return {"F": f_obs, "R2": r2, "p_value": p_value, "permutations": n_perm}


def betadisper_test(coords, groups, n_perm=999, seed=1):
    """Approximate PERMDISP2/betadisper dispersion check using distance-to-centroid in truncated PCoA space."""
    groups = np.asarray(groups)
    uniq = sorted(set(groups))

    def dist_to_centroid(grp):
        vals = np.zeros(len(grp))
        for g in uniq:
            idx = np.where(grp == g)[0]
            if len(idx) == 0:
                continue
            centroid = coords[idx].mean(axis=0)
            vals[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)
        return vals

    def anova_f(vals, grp):
        grand_mean = vals.mean()
        ss_b = sum(len(vals[grp == g]) * (vals[grp == g].mean() - grand_mean) ** 2 for g in uniq if len(vals[grp == g]))
        ss_w = sum(((vals[grp == g] - vals[grp == g].mean()) ** 2).sum() for g in uniq if len(vals[grp == g]))
        dfb = len(uniq) - 1
        dfw = len(vals) - len(uniq)
        if dfb <= 0 or dfw <= 0 or ss_w <= 0:
            return float("nan")
        return (ss_b / dfb) / (ss_w / dfw)

    obs_vals = dist_to_centroid(groups)
    f_obs = anova_f(obs_vals, groups)
    if math.isnan(f_obs):
        return {"F": float("nan"), "p_value": float("nan"), "permutations": 0}
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        f_perm = anova_f(obs_vals, perm)
        if not math.isnan(f_perm) and f_perm >= f_obs:
            count += 1
    p_value = (count + 1) / (n_perm + 1)
    return {"F": f_obs, "p_value": p_value, "permutations": n_perm}


def taxon_id_column(fieldnames):
    """Accept either 'taxon' or 'feature' as the first-column id per references/input-contract.md."""
    fieldnames = list(fieldnames)
    for cand in ("taxon", "feature"):
        if cand in fieldnames:
            return cand
    raise SystemExit(f"Abundance table must have a 'taxon' or 'feature' id column; found columns: {fieldnames}")


def rank_average(values):
    """Average (mid-)ranks with tie correction, for a rank-based (Kruskal-Wallis-style) test."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values))
    sorted_vals = values[order]
    i = 0
    n = len(values)
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    return ranks


def kruskal_h(values, groups):
    groups = np.asarray(groups)
    ranks = rank_average(values)
    n = len(values)
    h = 0.0
    for g in sorted(set(groups)):
        idx = np.where(groups == g)[0]
        n_g = len(idx)
        if n_g == 0:
            continue
        h += (ranks[idx].sum() ** 2) / n_g
    return (12.0 / (n * (n + 1))) * h - 3 * (n + 1)


def diversity_group_test(values, groups, n_perm=999, seed=0):
    """Permutation-based Kruskal-Wallis-style test for a group difference in an alpha diversity metric.

    Dependency-light equivalent of stats::kruskal.test / stats::wilcox.test used by the R route; uses an
    exact permutation p-value rather than a chi-square/normal approximation.
    """
    groups = np.asarray(groups)
    uniq = sorted(set(groups))
    counts = {g: int((groups == g).sum()) for g in uniq}
    if len(uniq) < 2 or any(c < 2 for c in counts.values()):
        return {"h": float("nan"), "p_value": float("nan"), "insufficient_replicates": True}
    values = np.asarray(values, dtype=float)
    h_obs = kruskal_h(values, groups)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        if kruskal_h(values, perm) >= h_obs:
            count += 1
    p_value = (count + 1) / (n_perm + 1)
    return {"h": h_obs, "p_value": p_value, "insufficient_replicates": False}


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment (numpy-only), used when >1 alpha diversity metric is tested."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    parser = argparse.ArgumentParser(description="Microbiome quick plots from abundance and metadata tables.")
    parser.add_argument("--abundance", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--group-col", default="group")
    parser.add_argument("--out", default="output")
    parser.add_argument("--top-taxa", type=int, default=9, help="Number of top taxa to show individually before collapsing the rest into 'Other' (max 9, per visual-style <=10 stacked segments).")
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.abundance)
    meta = read_csv(args.metadata)
    if not rows:
        raise SystemExit(f"Abundance table '{args.abundance}' has no data rows.")
    if not meta:
        raise SystemExit(f"Metadata table '{args.metadata}' has no data rows.")
    id_col = taxon_id_column(rows[0].keys())
    samples = [k for k in rows[0].keys() if k != id_col]
    if not samples:
        raise SystemExit("Abundance table has no sample columns besides the taxon/feature id column.")
    if "sample" not in meta[0] or args.group_col not in meta[0]:
        raise SystemExit(
            f"Metadata table must have a 'sample' column and a '{args.group_col}' group column; "
            f"found columns: {list(meta[0].keys())}"
        )
    groups = {r["sample"]: r[args.group_col] for r in meta}
    taxa = [r[id_col] for r in rows]
    try:
        matrix = np.array([[float(r[s]) for s in samples] for r in rows], dtype=float)
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"Abundance table contains a non-numeric or missing value ({exc}); check for blanks/text in sample columns.") from exc
    unmatched = [s for s in samples if s not in groups]
    matched = [s for s in samples if s in groups]
    if not matched:
        raise SystemExit(
            "No sample IDs in the abundance table match the metadata 'sample' column; check sample-ID spelling/casing."
        )
    if unmatched:
        print(
            f"WARNING: {len(unmatched)} sample(s) in the abundance table have no metadata match and are "
            f"grouped as 'NA': {', '.join(unmatched)}",
            file=sys.stderr,
        )
    rel = matrix / np.maximum(matrix.sum(axis=0, keepdims=True), 1)
    sample_group = np.array([groups.get(s, "NA") for s in samples])
    group_names = sorted({a for a in sample_group})

    ## ---- Read depth / QC -----------------------------------------------------
    library_size = matrix.sum(axis=0)
    zero_depth_samples = [samples[i] for i in range(len(samples)) if library_size[i] == 0]
    low_depth_samples = [samples[i] for i in range(len(samples)) if 0 < library_size[i] < 1000]
    nonzero_depth = library_size[library_size > 0]
    depth_ratio = float(nonzero_depth.max() / nonzero_depth.min()) if len(nonzero_depth) else float("nan")
    uneven_depth = bool(nonzero_depth.size) and depth_ratio >= 10
    fig, ax = plt.subplots(figsize=(max(8.5, len(samples) * 0.5 + 3), 5.8))
    color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(group_names)}
    order = np.argsort(-library_size)
    ax.bar([samples[i] for i in order], library_size[order], color=[color_map[sample_group[i]] for i in order], width=0.72)
    ax.set_title("Read depth per sample", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Library size (reads)")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(True, axis="y", color=GRID)
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[g]) for g in group_names]
    ax.legend(handles, group_names, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), title="Group")
    save(fig, out / "figures" / "read_depth")

    ## ---- Alpha diversity: Observed / Shannon / Simpson -----------------------
    alpha = [
        {
            "sample": s,
            "group": groups.get(s, "NA"),
            "observed_taxa": int((matrix[:, i] > 0).sum()),
            "shannon": shannon(matrix[:, i]),
            "simpson": simpson(matrix[:, i]),
        }
        for i, s in enumerate(samples)
    ]
    with (out / "results" / "alpha_diversity.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "group", "observed_taxa", "shannon", "simpson"], delimiter="\t")
        writer.writeheader()
        writer.writerows(alpha)

    metrics = [("observed_taxa", "Observed taxa"), ("shannon", "Shannon diversity"), ("simpson", "Simpson diversity")]

    ## ---- Alpha diversity group test (permutation Kruskal-Wallis-style, BH-adjusted) ----
    alpha_tests = {}
    raw_pvals, keys_for_p = [], []
    for seed_offset, (key, label) in enumerate(metrics):
        vals = np.array([a[key] for a in alpha], dtype=float)
        test = diversity_group_test(vals, sample_group, n_perm=args.permutations, seed=10 + seed_offset)
        group_means = {g: float(np.mean([a[key] for a in alpha if a["group"] == g])) for g in group_names}
        higher_group = max(group_means, key=group_means.get) if len(group_means) >= 2 else None
        alpha_tests[key] = {"label": label, "group_means": group_means, "higher_group": higher_group, "p_adj": None, **test}
        if not test["insufficient_replicates"] and not math.isnan(test["p_value"]):
            raw_pvals.append(test["p_value"])
            keys_for_p.append(key)
    if raw_pvals:
        for key, p_adj in zip(keys_for_p, bh_adjust(raw_pvals)):
            alpha_tests[key]["p_adj"] = float(p_adj)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.6))
    for ax, (key, label) in zip(axes, metrics):
        data = [[a[key] for a in alpha if a["group"] == g] for g in group_names]
        bp = ax.boxplot(data, patch_artist=True)
        for patch, color in zip(bp["boxes"], PALETTE):
            patch.set_facecolor(color)
            patch.set_alpha(0.24)
            patch.set_edgecolor(color)
        for i, vals in enumerate(data, start=1):
            jitter = np.random.default_rng(i).uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, color=PALETTE[(i - 1) % len(PALETTE)], s=22, alpha=0.85, zorder=3)
        ax.set_xticks(range(1, len(group_names) + 1))
        ax.set_xticklabels([f"{g}\nn={len(v)}" for g, v in zip(group_names, data)])
        ax.set_title(label, loc="center", pad=14, fontweight="bold", color=TEXT, fontsize=18)
        ax.set_xlabel("Group")
        ax.grid(True, axis="y", color=GRID)
    fig.suptitle("Alpha diversity", fontweight="bold", color=TEXT, fontsize=24, y=1.03)
    save(fig, out / "figures" / "alpha_diversity")

    ## ---- Taxa composition: top taxa + Other -----------------------------------
    top_n = max(1, min(args.top_taxa, 9, len(taxa)))
    order_idx = np.argsort(-rel.mean(axis=1))
    top = order_idx[:top_n]
    other = np.setdiff1d(np.arange(len(taxa)), top)
    comp = rel[top, :]
    if len(other):
        comp = np.vstack([comp, rel[other, :].sum(axis=0)])
        labels = [taxa[i] for i in top] + ["Other"]
    else:
        labels = [taxa[i] for i in top]
    sample_order = sorted(range(len(samples)), key=lambda i: (sample_group[i], samples[i]))
    ordered_samples = [samples[i] for i in sample_order]
    comp_ordered = comp[:, sample_order]
    fig, ax = plt.subplots(figsize=(max(8.5, len(samples) * 0.55 + 3), 6.2))
    bottom = np.zeros(len(ordered_samples))
    for i, label in enumerate(labels):
        ax.bar(ordered_samples, comp_ordered[i] * 100, bottom=bottom * 100, color=PALETTE[i % len(PALETTE)], width=0.72, label=label)
        bottom += comp_ordered[i]
    ax.set_title("Taxa composition", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Relative abundance (%)")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    save(fig, out / "figures" / "taxa_composition")

    ## ---- Beta diversity: Bray-Curtis + PCoA + PERMANOVA + betadisper ----------
    bray = bray_curtis_distance(matrix)
    coords, var_explained = classical_mds(bray, k=2)
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    for g in group_names:
        idx = [i for i, sample in enumerate(samples) if groups.get(sample, "NA") == g]
        ax.scatter(coords[idx, 0], coords[idx, 1], s=72, color=color_map[g], label=g, alpha=0.86)
    ax.set_title("PCoA (Bray-Curtis)", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xlabel(f"PCoA1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({var_explained[1]:.1f}%)")
    ax.grid(True, color=GRID)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    save(fig, out / "figures" / "ordination")

    if len(group_names) >= 2:
        permanova = permanova_test(bray, sample_group, n_perm=args.permutations, seed=0)
        dispersion = betadisper_test(coords, sample_group, n_perm=args.permutations, seed=1)
    else:
        permanova = {"F": float("nan"), "R2": float("nan"), "p_value": float("nan"), "permutations": 0}
        dispersion = {"F": float("nan"), "p_value": float("nan"), "permutations": 0}
    with (out / "results" / "permanova.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["term", "F", "R2", "p_value", "permutations"])
        writer.writerow(["group (PERMANOVA, Bray-Curtis)", f"{permanova['F']:.6g}", f"{permanova['R2']:.6g}", f"{permanova['p_value']:.6g}", permanova["permutations"]])
        writer.writerow(["betadisper (approx., PCoA1-2 based)", f"{dispersion['F']:.6g}", "NA", f"{dispersion['p_value']:.6g}", dispersion["permutations"]])

    ## ---- Result tables and run metadata ---------------------------------------
    with (out / "results" / "taxa_relative_abundance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", *samples])
        for taxon, vals in zip(taxa, rel):
            writer.writerow([taxon, *[f"{v:.6g}" for v in vals]])
    (out / "summary.tsv").write_text(
        "item\tvalue\n"
        f"samples\t{len(samples)}\n"
        f"taxa\t{len(taxa)}\n"
        f"permanova_R2\t{permanova['R2']:.4g}\n"
        f"permanova_p\t{permanova['p_value']:.4g}\n"
        f"betadisper_p\t{dispersion['p_value']:.4g}\n",
        encoding="utf-8",
    )
    (out / "run" / "qc.json").write_text(
        json.dumps(
            {
                "samples": len(samples),
                "taxa": len(taxa),
                "min_library_size": float(library_size.min()),
                "max_library_size": float(library_size.max()),
                "zero_depth_samples": zero_depth_samples,
                "low_depth_samples": low_depth_samples,
                "unmatched_samples": unmatched,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    artifacts = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
    (out / "run" / "artifacts.json").write_text(json.dumps(artifacts, indent=2), encoding="utf-8")

    ## ---- Exploratory (non-statistical) taxa shift, 2-group case only ----------
    top_shift = []
    if len(group_names) == 2:
        g0, g1 = group_names
        idx0 = [i for i, s in enumerate(samples) if sample_group[i] == g0]
        idx1 = [i for i, s in enumerate(samples) if sample_group[i] == g1]
        if idx0 and idx1:
            mean0 = rel[:, idx0].mean(axis=1)
            mean1 = rel[:, idx1].mean(axis=1)
            diff = mean1 - mean0
            for i in np.argsort(-np.abs(diff))[:5]:
                if abs(diff[i]) < 1e-9:
                    continue
                direction = g1 if diff[i] > 0 else g0
                top_shift.append((taxa[i], direction, abs(diff[i]) * 100))

    ## ---- report.md: results-and-conclusions summary ----------------------------
    group_counts = {g: int((sample_group == g).sum()) for g in group_names}
    lines = ["# Microbiome quicklook\n\n"]
    lines.append(f"- Samples: {len(samples)} ({', '.join(f'{g}: n={n}' for g, n in group_counts.items())})\n")
    lines.append(f"- Taxa (features): {len(taxa)}\n\n")

    lines.append("## Alpha diversity\n\n")
    for key, label in metrics:
        t = alpha_tests[key]
        means_str = ", ".join(f"{g}={t['group_means'][g]:.3g}" for g in group_names)
        if t["insufficient_replicates"]:
            lines.append(f"- {label}: mean by group ({means_str}); statistical test skipped (a group has fewer than 2 samples).\n")
        else:
            p_adj = t["p_adj"]
            sig = "significant" if (p_adj is not None and p_adj < 0.05) else "not significant"
            p_txt = f"p = {t['p_value']:.4g}, BH-adjusted p = {p_adj:.4g}" if p_adj is not None else f"p = {t['p_value']:.4g}"
            lines.append(f"- {label}: mean by group ({means_str}); highest in **{t['higher_group']}**; permutation test {p_txt} ({sig}).\n")
    lines.append("\n")

    lines.append("## Beta diversity (Bray-Curtis PERMANOVA)\n\n")
    if len(group_names) < 2:
        lines.append("- Only one group present in metadata; PERMANOVA/dispersion group comparison skipped.\n\n")
    else:
        perm_sig = "significant group separation" if permanova["p_value"] < 0.05 else "no significant group separation"
        lines.append(
            f"- PERMANOVA ({permanova['permutations']} permutations): R2 = {permanova['R2']:.4g}, "
            f"F = {permanova['F']:.4g}, p = {permanova['p_value']:.4g} -> {perm_sig}.\n"
        )
        if math.isnan(dispersion["p_value"]):
            lines.append("- betadisper dispersion check could not be computed for this data.\n\n")
        elif dispersion["p_value"] < 0.05:
            lines.append(
                f"- betadisper dispersion check (approximate, PCoA-axis 1-2 based): p = {dispersion['p_value']:.4g} "
                "(significant) -- within-group dispersion differs between groups, so part of the PERMANOVA result "
                "may reflect unequal spread rather than a pure centroid shift (references/pitfalls.md #5); interpret "
                "the PERMANOVA result with this caveat. For a rigorous full-dimensional betadisper, use "
                "scripts/run_phyloseq_ancombc.R.\n\n"
            )
        else:
            lines.append(
                f"- betadisper dispersion check (approximate, PCoA-axis 1-2 based): p = {dispersion['p_value']:.4g} "
                "(not significant) -- dispersion is comparable between groups, supporting a centroid-shift "
                "interpretation of the PERMANOVA result.\n\n"
            )

    lines.append("## Exploratory taxa shift (descriptive only, not a statistical test)\n\n")
    if top_shift:
        for taxon, direction, pct in top_shift:
            lines.append(f"- {taxon}: higher mean relative abundance in **{direction}** (delta ~= {pct:.2f} percentage points)\n")
        lines.append(
            "- These are descriptive relative-abundance mean differences ranked by magnitude, not a differential "
            "abundance test, and relative abundances are compositional (references/pitfalls.md #1, #2); use "
            "scripts/run_phyloseq_ancombc.R (ANCOM-BC2) for a compositionality-aware statistical test.\n\n"
        )
    else:
        lines.append(
            "- Not computed (this descriptive comparison requires exactly two groups). No differential abundance "
            "test is run by this script; use scripts/run_phyloseq_ancombc.R for ANCOM-BC2.\n\n"
        )

    caveats = []
    if zero_depth_samples:
        caveats.append(f"Sample(s) with zero total reads: {', '.join(zero_depth_samples)} -- excluded from meaningful diversity interpretation.")
    if low_depth_samples:
        caveats.append(
            f"Low read-depth sample(s) (<1000 reads): {', '.join(low_depth_samples)} -- review before treating "
            "diversity/composition as comparable across samples (references/pitfalls.md #3)."
        )
    if uneven_depth:
        caveats.append(
            f"Read depth varies about {depth_ratio:.0f}x between the deepest and shallowest sample -- plot "
            "diversity vs. depth before trusting alpha-diversity group differences (references/pitfalls.md #3)."
        )
    if unmatched:
        caveats.append(
            f"{len(unmatched)} sample ID(s) in the abundance table had no metadata match and were grouped as "
            f"'NA': {', '.join(unmatched)} -- verify sample IDs before trusting group comparisons."
        )
    if any(n == 1 for n in group_counts.values()):
        singleton_groups = [g for g, n in group_counts.items() if n == 1]
        caveats.append(f"Group(s) with only 1 sample: {', '.join(singleton_groups)} -- statistics for these groups have essentially no within-group variance estimate.")
    if len(group_names) >= 2 and not math.isnan(dispersion["p_value"]) and dispersion["p_value"] < 0.05:
        caveats.append("betadisper flagged a significant dispersion difference between groups (see Beta diversity section above).")
    if top_shift:
        caveats.append("Compositionality reminder: the exploratory taxa shift above uses raw relative-abundance means (not a compositional-aware method); see references/pitfalls.md #1, #2.")
    if caveats:
        lines.append("## Caveats\n\n")
        for c in caveats:
            lines.append(f"- {c}\n")
        lines.append("\n")

    lines.append("## Conclusion\n\n")
    conclusion = []
    if len(group_names) < 2:
        conclusion.append("Only one group is present in the metadata, so no between-group diversity or composition comparison could be made.")
    else:
        shannon_t = alpha_tests["shannon"]
        if shannon_t["insufficient_replicates"]:
            conclusion.append("Alpha diversity (Shannon) could not be tested statistically because at least one group has fewer than 2 samples.")
        elif shannon_t["p_adj"] is not None:
            if shannon_t["p_adj"] < 0.05:
                conclusion.append(f"Shannon diversity was significantly higher in {shannon_t['higher_group']} (BH-adjusted p = {shannon_t['p_adj']:.4g}).")
            else:
                conclusion.append(
                    f"Shannon diversity did not differ significantly between groups (BH-adjusted p = {shannon_t['p_adj']:.4g}), "
                    f"though the mean was nominally higher in {shannon_t['higher_group']}."
                )
        if not math.isnan(permanova["p_value"]):
            if permanova["p_value"] < 0.05:
                disp_note = (
                    ", though this may partly reflect a dispersion difference rather than a pure centroid shift"
                    if (not math.isnan(dispersion["p_value"]) and dispersion["p_value"] < 0.05)
                    else ""
                )
                conclusion.append(f"Community composition (Bray-Curtis PERMANOVA) differed significantly between groups (R2 = {permanova['R2']:.4g}, p = {permanova['p_value']:.4g}){disp_note}.")
            else:
                conclusion.append(f"Community composition did not differ significantly between groups by PERMANOVA (R2 = {permanova['R2']:.4g}, p = {permanova['p_value']:.4g}).")
        if top_shift:
            top_taxon, top_dir, top_pct = top_shift[0]
            conclusion.append(
                f"The largest descriptive taxa shift was {top_taxon}, higher in {top_dir} by about {top_pct:.2f} "
                "percentage points, though this is exploratory and not a formal differential abundance call."
            )
    lines.append(" ".join(conclusion) + "\n")

    (out / "report.md").write_text("".join(lines), encoding="utf-8")
    print(f"wrote microbiome outputs to {out}")


if __name__ == "__main__":
    main()
