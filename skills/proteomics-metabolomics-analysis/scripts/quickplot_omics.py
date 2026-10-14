#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except Exception as exc:
    raise SystemExit("Install matplotlib and numpy: python3 -m pip install matplotlib numpy") from exc


UP = "#E45B5B"
DOWN = "#4EA3F1"
NS = "#B8C2CC"
TEXT = "#173B5F"
GRID = "#DCEAF5"
AMBER = "#E2A13B"

HIGH_MISSING_FLAG = 0.3  # reference line only (visual flag), not a hard cutoff
DEFAULT_SIG_PADJ = 0.05
DEFAULT_SIG_LOG2FC = 1.0
DEFAULT_MIN_GROUP_COMPLETENESS = 0.5  # matches MIN_GROUP_COMPLETENESS in run_omics_templates.R
SMALL_GROUP_N = 3  # below this per-group n, flag pitfalls.md #8 (multiple-testing burden at small n)


def read_csv(path):
    rows = list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))
    return rows


def require_columns(rows, required, label):
    if not rows:
        raise SystemExit(f"{label} file has no data rows.")
    present = set(rows[0].keys())
    missing = [c for c in required if c not in present]
    if missing:
        raise SystemExit(f"{label} is missing required column(s): {', '.join(missing)} (found: {', '.join(present)})")


def to_float(value):
    """Parse a matrix cell to float, treating blank/NA/NaN tokens as missing (np.nan)."""
    text = "" if value is None else str(value).strip()
    if text == "" or text.lower() in ("na", "nan", "null"):
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 16, "axes.titlesize": 24, "axes.labelsize": 21, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15, "axes.edgecolor": "#B9D7EE", "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT, "text.color": TEXT, "pdf.fonttype": 42, "ps.fonttype": 42})


def main():
    parser = argparse.ArgumentParser(description="Omics quick plots from abundance matrix, metadata and differential table.")
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--differential", required=True)
    parser.add_argument("--group-col", default="group")
    parser.add_argument("--out", default="output")
    parser.add_argument("--padj-threshold", type=float, default=DEFAULT_SIG_PADJ,
                         help=f"Adjusted p-value cutoff for up/down classification (default {DEFAULT_SIG_PADJ:g}).")
    parser.add_argument("--log2fc-threshold", type=float, default=DEFAULT_SIG_LOG2FC,
                         help=f"Absolute log2 fold-change cutoff for up/down classification (default {DEFAULT_SIG_LOG2FC:g}).")
    parser.add_argument("--min-group-completeness", type=float, default=DEFAULT_MIN_GROUP_COMPLETENESS,
                         help="Fraction of non-missing values required WITHIN EACH group for a feature to pass the "
                              f"missingness-completeness filter (informational; default {DEFAULT_MIN_GROUP_COMPLETENESS:g}, "
                              "matches MIN_GROUP_COMPLETENESS in run_omics_templates.R).")
    args = parser.parse_args()
    if not (0 < args.padj_threshold <= 1):
        raise SystemExit(f"--padj-threshold must be in (0, 1]; got {args.padj_threshold}")
    if args.log2fc_threshold < 0:
        raise SystemExit(f"--log2fc-threshold must be >= 0; got {args.log2fc_threshold}")
    if not (0 <= args.min_group_completeness <= 1):
        raise SystemExit(f"--min-group-completeness must be in [0, 1]; got {args.min_group_completeness}")
    sig_padj = args.padj_threshold
    sig_log2fc = args.log2fc_threshold
    min_group_completeness = args.min_group_completeness
    style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.matrix)
    diff = read_csv(args.differential)
    meta = read_csv(args.metadata)

    require_columns(rows, ["feature"], "Matrix")
    require_columns(diff, ["feature", "log2FC", "padj"], "Differential table")
    require_columns(meta, ["sample", args.group_col], "Metadata")

    samples = [k for k in rows[0].keys() if k != "feature"]
    if not samples:
        raise SystemExit("Matrix file has a 'feature' column but no sample columns.")

    features = [r["feature"] for r in rows]
    dup_features = sorted({f for f in features if features.count(f) > 1})
    if dup_features:
        raise SystemExit(
            "Matrix has duplicate feature IDs (each feature must be unique): "
            + ", ".join(dup_features[:10]) + (" ..." if len(dup_features) > 10 else "")
        )

    meta_by_sample = {r["sample"]: r for r in meta}
    missing_meta = [s for s in samples if s not in meta_by_sample]
    if missing_meta:
        raise SystemExit(f"Sample(s) in matrix are missing from metadata: {', '.join(missing_meta)}")
    groups = {s: meta_by_sample[s][args.group_col] for s in samples}
    group_names = sorted({groups.get(s, "NA") for s in samples})
    if len(group_names) < 2:
        print(f"Warning: metadata group column '{args.group_col}' has only one level ({group_names}); "
              "PCA/heatmap will still render but there is no case-vs-control contrast to visualize by color.")
    group_counts = {g: sum(1 for s in samples if groups.get(s, "NA") == g) for g in group_names}
    small_groups = [g for g, n in group_counts.items() if n < SMALL_GROUP_N]

    mat = np.array([[to_float(r[s]) for s in samples] for r in rows], dtype=float)
    palette = {g: [DOWN, UP, "#4CB984", "#8B72D6"][i % 4] for i, g in enumerate(group_names)}

    # PCA: use only fully-observed features (no missing values) so the
    # projection is not distorted by an imputation choice (see
    # references/pitfalls.md #1). If too few complete features remain, skip.
    complete_mask = ~np.isnan(mat).any(axis=1)
    pca_mat = mat[complete_mask]
    if pca_mat.shape[0] < 3:
        print(f"Warning: only {pca_mat.shape[0]} fully-observed features; skipping PCA plot.")
    else:
        x = pca_mat.T - pca_mat.T.mean(axis=0, keepdims=True)
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        coords = u[:, :2] * s[:2]
        fig, ax = plt.subplots(figsize=(7.8, 6.2))
        for g in group_names:
            idx = [i for i, sample in enumerate(samples) if groups.get(sample, "NA") == g]
            ax.scatter(coords[idx, 0], coords[idx, 1], s=72, color=palette[g], label=g, alpha=0.88)
            for i in idx:
                ax.annotate(samples[i], (coords[i, 0], coords[i, 1]), fontsize=9, color=TEXT,
                            xytext=(3, 3), textcoords="offset points")
        ax.set_title("PCA", loc="center", pad=18, fontweight="bold", color=TEXT)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.grid(True, color=GRID)
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
        save(fig, out / "figures" / "pca")

    for row_num, r in enumerate(diff, start=2):
        try:
            float(r["log2FC"])
            float(r["padj"])
        except (TypeError, ValueError):
            raise SystemExit(
                f"Differential table row {row_num} (feature={r.get('feature')!r}) has a non-numeric "
                f"log2FC={r.get('log2FC')!r} or padj={r.get('padj')!r}."
            )

    xs = [float(r["log2FC"]) for r in diff]
    padjs = [float(r["padj"]) for r in diff]
    ys = [-math.log10(max(p, 1e-12)) for p in padjs]
    direction = [
        "up" if x >= sig_log2fc and p < sig_padj else "down" if x <= -sig_log2fc and p < sig_padj else "ns"
        for x, p in zip(xs, padjs)
    ]
    colors = [UP if d == "up" else DOWN if d == "down" else NS for d in direction]
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    ax.scatter(xs, ys, c=colors, s=34, alpha=0.84, linewidths=0)
    ax.axvline(-sig_log2fc, color="#A7C6DD", linestyle="--")
    ax.axvline(sig_log2fc, color="#A7C6DD", linestyle="--")
    ax.axhline(-math.log10(sig_padj), color="#A7C6DD", linestyle="--")
    ax.set_title(f"Differential features (padj<{sig_padj:g}, |log2FC|>={sig_log2fc:g})",
                 loc="center", pad=18, fontweight="bold", color=TEXT, fontsize=19)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p value")
    ax.grid(True, color=GRID)
    save(fig, out / "figures" / "volcano")

    n_up = direction.count("up")
    n_down = direction.count("down")

    diff_order = sorted(diff, key=lambda r: float(r["padj"]))[:30]
    top_names = [r["feature"] for r in diff_order]
    unmatched = [f for f in top_names if f not in features]
    if unmatched:
        print(f"Warning: {len(unmatched)} top differential feature ID(s) not found in the matrix "
              f"(e.g. {unmatched[:5]}); they are skipped in the heatmap.")
    idx_all = [features.index(f) for f in top_names if f in features]
    # Exclude features that are entirely missing across all samples: there is
    # no observed value to display or z-score, and a row-mean fill would be
    # fabricated data rather than a display convenience.
    idx = [i for i in idx_all if not np.isnan(mat[i, :]).all()]
    n_fully_missing_dropped = len(idx_all) - len(idx)
    if n_fully_missing_dropped:
        print(f"Note: {n_fully_missing_dropped} top differential feature(s) are missing in every sample "
              "and were excluded from the heatmap (no data to display).")
    if not idx:
        print("Warning: none of the top differential features have any observed data; skipping top-feature heatmap.")
    else:
        sub = mat[idx, :].copy()
        resid_na = int(np.isnan(sub).sum())
        if resid_na:
            # Visualization-only fill (never used for statistics): replace any
            # residual missing cell among the plotted top features with that
            # feature's own row mean so the heatmap renders cleanly.
            print(f"Note: {resid_na} residual missing cell(s) among top {len(idx)} heatmap features filled with the per-feature row mean for display only.")
            row_means = np.nanmean(sub, axis=1, keepdims=True)
            nan_pos = np.isnan(sub)
            sub[nan_pos] = np.broadcast_to(row_means, sub.shape)[nan_pos]
        z = (sub - sub.mean(axis=1, keepdims=True)) / np.maximum(sub.std(axis=1, keepdims=True), 1e-6)
        fig, ax = plt.subplots(figsize=(max(7.4, len(samples) * 0.5 + 3), max(6.2, len(idx) * 0.2 + 2.8)))
        im = ax.imshow(z, cmap="RdBu_r", vmin=-2.2, vmax=2.2, aspect="auto")
        ax.set_title("Top feature heatmap", loc="center", pad=18, fontweight="bold", color=TEXT)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=35, ha="right")
        ax.set_yticks(range(len(idx)))
        ax.set_yticklabels([features[i] for i in idx], fontsize=12)
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
        save(fig, out / "figures" / "top_feature_heatmap")

    # Missingness: per-sample and per-feature, written as one tidy table and
    # visualized as a bar chart (cns-method-guide.md lists "Missingness
    # heatmap/bar" as the first Required Figure). Derived from the already-
    # parsed numeric matrix (np.isnan) rather than re-parsing raw strings, so
    # this is consistent with to_float()'s NA-token handling and safe against
    # ragged/short CSV rows.
    n_features = len(rows)
    n_samples = len(samples)
    missing_mask = np.isnan(mat)
    sample_missing = {s: int(missing_mask[:, j].sum()) for j, s in enumerate(samples)}
    feature_missing = {f: int(missing_mask[i, :].sum()) for i, f in enumerate(features)}

    # Group-aware completeness filter (informational; see MIN_GROUP_COMPLETENESS
    # in run_omics_templates.R). This quicklook script does not drop matrix rows
    # (it only visualizes the matrix/differential table as supplied), but it
    # reports, for the --min-group-completeness threshold actually requested,
    # how many features would pass a rule requiring >= that fraction of
    # non-missing values WITHIN EACH group (case AND control), matching the
    # formal R pipeline's group-aware missingness filter (pitfalls.md #11).
    group_sample_idx = {g: [j for j, s in enumerate(samples) if groups.get(s, "NA") == g] for g in group_names}
    feature_pass_filter = {}
    for i, f in enumerate(features):
        ok = True
        for g in group_names:
            idxs = group_sample_idx[g]
            if not idxs:
                continue
            frac_non_missing = sum(1 for j in idxs if not np.isnan(mat[i, j])) / len(idxs)
            if frac_non_missing < min_group_completeness:
                ok = False
                break
        feature_pass_filter[f] = ok
    n_pass_filter = sum(1 for f in features if feature_pass_filter[f])
    n_fail_filter = n_features - n_pass_filter

    missingness_rows = []
    for f in features:
        n_miss = feature_missing[f]
        missingness_rows.append({"scope": "feature", "id": f, "group": "", "n_total": n_samples,
                                  "n_missing": n_miss, "frac_missing": round(n_miss / n_samples, 4),
                                  "pass_filter": feature_pass_filter[f]})
    for s in samples:
        n_miss = sample_missing[s]
        missingness_rows.append({"scope": "sample", "id": s, "group": groups.get(s, "NA"), "n_total": n_features,
                                  "n_missing": n_miss, "frac_missing": round(n_miss / n_features, 4),
                                  "pass_filter": ""})
    with (out / "results" / "missingness_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scope", "id", "group", "n_total", "n_missing", "frac_missing", "pass_filter"], delimiter="\t")
        writer.writeheader()
        writer.writerows(missingness_rows)

    sample_order = sorted(samples, key=lambda s: groups.get(s, "NA"))
    fracs = [sample_missing[s] / n_features for s in sample_order]
    bar_colors = [palette.get(groups.get(s, "NA"), NS) for s in sample_order]
    fig, ax = plt.subplots(figsize=(max(7.4, len(samples) * 0.42 + 3), 5.6))
    ax.bar(range(len(sample_order)), fracs, color=bar_colors, width=0.68, alpha=0.92)
    ax.axhline(HIGH_MISSING_FLAG, color=AMBER, linestyle="--", linewidth=1.2)
    ax.set_title("Per-sample missingness", loc="center", pad=18, fontweight="bold", color=TEXT)
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=40, ha="right")
    ax.set_ylabel("Fraction missing features")
    ax.set_ylim(0, max(0.35, max(fracs) * 1.15 if fracs else 0.35))
    ax.grid(True, axis="y", color=GRID)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[g]) for g in group_names]
    ax.legend(handles, group_names, frameon=False, loc="upper right")
    save(fig, out / "figures" / "missingness")

    (out / "results" / "differential_features.tsv").write_text(Path(args.differential).read_text(encoding="utf-8"), encoding="utf-8")
    (out / "summary.tsv").write_text(f"item\tvalue\nsamples\t{len(samples)}\nfeatures\t{len(features)}\n", encoding="utf-8")
    (out / "run" / "qc.json").write_text(json.dumps({
        "samples": len(samples),
        "features": len(features),
        "padj_threshold": sig_padj,
        "log2fc_threshold": sig_log2fc,
        "min_group_completeness": min_group_completeness,
        "features_passing_completeness_filter": n_pass_filter,
        "features_failing_completeness_filter": n_fail_filter,
    }, indent=2), encoding="utf-8")

    # ---- Results-and-conclusions report.md ----
    high_miss_features = [f for f in features if feature_missing[f] / n_samples > HIGH_MISSING_FLAG]
    high_miss_samples = [s for s in samples if sample_missing[s] / n_features > HIGH_MISSING_FLAG]

    top_report = sorted(
        [{"feature": r["feature"], "log2FC": float(r["log2FC"]), "padj": float(r["padj"]),
          "direction": d} for r, d in zip(diff, direction)],
        key=lambda r: r["padj"],
    )[:10]

    lines = ["# Omics quicklook report", ""]
    lines.append("## Input summary")
    lines.append(f"- Features in matrix: {n_features}")
    group_desc = ", ".join(f"{g}={group_counts[g]}" for g in group_names)
    lines.append(f"- Samples: {n_samples} ({group_desc})")
    lines.append(f"- Features tested in supplied differential table: {len(diff)}")
    lines.append("")

    lines.append("## Missingness overview")
    lines.append(
        f"- Features flagged high-missingness (>{int(HIGH_MISSING_FLAG * 100)}% missing across samples): "
        f"{len(high_miss_features)} / {n_features} ({100 * len(high_miss_features) / max(n_features, 1):.1f}%)"
    )
    samp_line = (
        f"- Samples flagged high-missingness (>{int(HIGH_MISSING_FLAG * 100)}% missing features): "
        f"{len(high_miss_samples)} / {n_samples}"
    )
    if high_miss_samples:
        samp_line += " (" + ", ".join(high_miss_samples) + ")"
    lines.append(samp_line)
    lines.append(
        f"- Missingness-completeness filter (informational; keep features with >= {min_group_completeness * 100:.0f}% "
        f"non-missing values in EACH group [{', '.join(group_names)}], matching MIN_GROUP_COMPLETENESS in "
        f"`run_omics_templates.R`): {n_pass_filter} / {n_features} feature(s) would pass, {n_fail_filter} would fail."
    )
    lines.append("- This quicklook script does not drop or impute matrix rows itself; the completeness check above is "
                  "reported for QC only (compare with the `pass_filter` column in `results/missingness_summary.tsv`); "
                  "run `run_omics_templates.R` to actually apply this filter before differential testing.")
    lines.append("")

    lines.append("## Differential results (from supplied differential table)")
    lines.append(f"- Significant up (log2FC >= {sig_log2fc:g}, padj < {sig_padj:g}): {n_up}")
    lines.append(f"- Significant down (log2FC <= -{sig_log2fc:g}, padj < {sig_padj:g}): {n_down}")
    lines.append(f"- Not significant / below threshold: {len(diff) - n_up - n_down}")
    lines.append("")
    lines.append("### Top differential features (by padj)")
    if top_report:
        for r in top_report:
            lines.append(f"- {r['feature']}: log2FC={r['log2FC']:.2f}, padj={r['padj']:.2e} ({r['direction']})")
    else:
        lines.append("- No rows in the supplied differential table.")
    lines.append("")

    lines.append("## Caveats")
    if high_miss_samples:
        lines.append(
            f"- {len(high_miss_samples)} sample(s) exceed {int(HIGH_MISSING_FLAG * 100)}% missing features "
            "(see amber line on `figures/missingness.png`); inspect these samples for run quality before "
            "trusting features that are largely observed only in the other samples (references/pitfalls.md #4, #11)."
        )
    if high_miss_features:
        lines.append(
            f"- {len(high_miss_features)} feature(s) exceed {int(HIGH_MISSING_FLAG * 100)}% missing values; "
            "check whether their missingness is concentrated in one group (possible MNAR / detection-limit "
            "signal, references/pitfalls.md #1) before ranking them by fold change."
        )
    lines.append(
        "- This quicklook script performs no imputation (MNAR-vs-MAR) — it only visualizes the matrix and the "
        "supplied differential table as-is. If imputation is required for downstream use, follow the decision "
        "tree in references/cns-method-guide.md and references/pitfalls.md #1 first."
    )
    lines.append(
        "- TMT/iTRAQ ratio-compression correction (references/pitfalls.md #2) and multi-omics integration "
        "(references/cns-method-guide.md) are out of scope for this quicklook and require dedicated workflows "
        "(see references/advanced-modules.md)."
    )
    if small_groups:
        lines.append(
            f"- Group(s) with n < {SMALL_GROUP_N} ({', '.join(f'{g}={group_counts[g]}' for g in small_groups)}): "
            "per-feature variance estimates are unstable at this sample size; prefer moderated/empirical-Bayes "
            "testing (e.g. limma, `run_omics_templates.R`) over plain t-tests (references/pitfalls.md #8)."
        )
    if len(group_names) < 2:
        lines.append(
            f"- Metadata group column '{args.group_col}' has only one level; no case-vs-control contrast could "
            "be visualized by color in PCA/missingness plots."
        )
    if n_fully_missing_dropped:
        lines.append(
            f"- {n_fully_missing_dropped} top-ranked differential feature(s) have zero observed values in this "
            "matrix and were omitted from the heatmap; verify feature ID matching between the differential "
            "table and the matrix."
        )
    lines.append("")

    lines.append("## Conclusion")
    if n_up + n_down == 0:
        concl = (
            f"Of {len(diff)} tested features across {n_samples} samples, none met the "
            f"log2FC >= {sig_log2fc:g} / padj < {sig_padj:g} threshold in either direction. "
            "This may reflect a true absence of a strong group effect, insufficient statistical power at the "
            "current sample size, or a threshold that is too strict for this dataset — inspect the volcano plot "
            "and consider whether the study is adequately powered before concluding there is no effect."
        )
    else:
        top_name = top_report[0]["feature"] if top_report else "the top feature"
        up_verb = "was" if n_up == 1 else "were"
        down_verb = "was" if n_down == 1 else "were"
        concl = (
            f"Of {len(diff)} tested features across {n_samples} samples, {n_up} {up_verb} up and {n_down} {down_verb} down "
            f"(log2FC >= {sig_log2fc:g}, padj < {sig_padj:g}); the strongest signal was {top_name} "
            f"(log2FC={top_report[0]['log2FC']:.2f}, padj={top_report[0]['padj']:.2e}). "
            "These results describe association/direction of change only — see the caveats above before treating "
            "any single feature's fold-change magnitude as a precisely measured, imputation-independent effect, "
            "and confirm the biological/annotation confidence of any feature central to the conclusion."
        )
    lines.append(concl)
    lines.append("")

    (out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote omics outputs to {out}")


if __name__ == "__main__":
    main()
