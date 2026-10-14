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
except Exception as exc:
    raise SystemExit(
        "quickplot_communication.py requires matplotlib to export PNG/PDF. "
        "Install it with: python3 -m pip install matplotlib"
    ) from exc


# Categorical palette and text color match assets/palette.json ("clusters" / "text") so that
# Python quicklook figures and the R (run_cellchat.R etc.) figures use the same visual identity.
PALETTE = ["#3B82C4", "#D84A4A", "#3FA37A", "#D99A2B", "#7C65C1", "#4DA6A6", "#F08AB4"]
TEXT = "#1F2937"
GRID = "#DCEAF5"


def read_csv(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def safe_float(value, default=0.0):
    """Parse a CSV cell to float; fall back to `default` (not an exception) for blank/NaN/bad values."""
    if value is None:
        return default
    text = str(value).strip()
    if text == "" or text.upper() in ("NA", "NAN", "NULL", "NONE"):
        return default
    try:
        parsed = float(text)
    except ValueError:
        return default
    return default if parsed != parsed else parsed


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


def draw_network(rows, out_base):
    celltypes = sorted({r["source"] for r in rows} | {r["target"] for r in rows})
    colors = {ct: PALETTE[i % len(PALETTE)] for i, ct in enumerate(celltypes)}
    n = len(celltypes)
    pos = {}
    for i, ct in enumerate(celltypes):
        angle = 2 * math.pi * i / n + math.pi / 8
        pos[ct] = (math.cos(angle), math.sin(angle))
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    ax.set_title("Cell-cell communication network", pad=18, loc="center", fontweight="bold", color=TEXT)
    for row in rows:
        sx, sy = pos[row["source"]]
        tx, ty = pos[row["target"]]
        prob = safe_float(row.get("probability"))
        ax.annotate(
            "",
            xy=(tx * 0.82, ty * 0.82),
            xytext=(sx * 0.82, sy * 0.82),
            arrowprops={
                "arrowstyle": "-|>",
                "lw": 1.2 + prob * 4.2,
                "color": colors[row["source"]],
                "alpha": 0.45 + prob * 0.45,
                "shrinkA": 18,
                "shrinkB": 18,
                "connectionstyle": "arc3,rad=0.18",
            },
        )
    for ct, (x, y) in pos.items():
        label = ct.replace("_", "\n")
        # A long cell type name on one line (e.g. "Myeloid", no underscore to split on) can be wider
        # than the fixed-size node circle and get clipped by the figure edge; shrink the font and
        # grow the marker for longer labels instead of letting text overrun the circle.
        longest_line = max(len(line) for line in label.split("\n"))
        fontsize = 14 if longest_line <= 5 else max(9, 14 - (longest_line - 5) * 1.1)
        marker_size = 1500 if longest_line <= 5 else 1500 + (longest_line - 5) * 220
        ax.scatter([x], [y], s=marker_size, color=colors[ct], alpha=0.92, edgecolor="white", linewidth=2.4, zorder=3)
        ax.text(x, y, label, ha="center", va="center", fontsize=fontsize, fontweight="bold", color="white", zorder=4)
    # Slightly wider than the unit circle layout to leave room for node labels/markers near the edge.
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")
    save_plot(fig, out_base)


def draw_dotplot(rows, out_base):
    pairs = [f"{r['ligand']}-{r['receptor']}" for r in rows]
    interactions = [f"{r['source']} -> {r['target']}" for r in rows]
    fig_h = max(5.4, 0.46 * len(rows) + 2.3)
    fig_w = max(11.8, len(rows) * 1.15 + 4.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    xs = list(range(len(rows)))
    ys = list(range(len(rows)))
    scores = [-math.log10(max(safe_float(r.get("pvalue"), 1.0), 1e-12)) for r in rows]
    sizes = [90 + safe_float(r.get("probability")) * 760 for r in rows]
    sc = ax.scatter(xs, ys, c=scores, s=sizes, cmap="Reds", edgecolors="#E45B5B", linewidths=0.9, alpha=0.88)
    ax.set_title("Ligand-receptor interactions", pad=18, loc="center", fontweight="bold", color=TEXT)
    ax.set_xlabel("Ligand-receptor pair", labelpad=12)
    ax.set_ylabel("Sender -> receiver", labelpad=12)
    ax.set_xticks(xs)
    ax.set_xticklabels(pairs, rotation=34, ha="right")
    ax.set_yticks(ys)
    ax.set_yticklabels(interactions)
    ax.grid(True, color=GRID, linewidth=0.8)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("-log10 p value", fontsize=16)
    fig.tight_layout()
    save_plot(fig, out_base)


def main():
    parser = argparse.ArgumentParser(description="Cell-cell communication quick plots from ligand-receptor edge table.")
    parser.add_argument("--edges", required=True)
    parser.add_argument("--out", default="output")
    args = parser.parse_args()
    apply_style()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.edges)
    if not rows:
        raise SystemExit(
            f"No rows found in edge table: {args.edges}. Check the file has a header row plus at least one data row."
        )
    required_cols = ["source", "target", "ligand", "receptor", "probability", "pvalue"]
    missing_cols = [c for c in required_cols if c not in rows[0]]
    if missing_cols:
        raise SystemExit(
            f"Edge table {args.edges} is missing required column(s): {', '.join(missing_cols)}. "
            f"Found columns: {list(rows[0].keys())}"
        )
    rows = sorted(rows, key=lambda r: safe_float(r.get("probability")), reverse=True)
    draw_network(rows, out / "figures" / "communication_network")
    draw_dotplot(rows, out / "figures" / "ligand_receptor_dotplot")
    with (out / "results" / "ligand_receptor_table.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    celltypes = sorted(set(r["source"] for r in rows) | set(r["target"] for r in rows))
    (out / "summary.tsv").write_text(
        f"item\tvalue\ninteractions\t{len(rows)}\ncelltypes\t{len(celltypes)}\n", encoding="utf-8"
    )
    (out / "run" / "qc.json").write_text(json.dumps({"interactions": len(rows), "quicklook": True}, indent=2), encoding="utf-8")
    (out / "methods.md").write_text(
        "# Methods\n\nCell-cell communication quicklook was generated from a ligand-receptor edge table. "
        "Edges were ranked by interaction probability; network edge width encodes probability and dot color encodes statistical significance. "
        "This visualization is exploratory and does not replace CellChat/NicheNet/LIANA model-based inference.\n",
        encoding="utf-8",
    )

    # ---- Results summary: top pairs, top senders/receivers, pathway breakdown when present ----
    def fmt_pair(r):
        p = safe_float(r.get("probability"))
        pv = safe_float(r.get("pvalue"), 1.0)
        return f"{r['source']} -> {r['target']}: {r['ligand']}-{r['receptor']} (prob={p:.2f}, p={pv:.2g})"

    top_pairs = [fmt_pair(r) for r in rows[:5]]

    sender_counts = {}
    receiver_counts = {}
    for r in rows:
        sender_counts[r["source"]] = sender_counts.get(r["source"], 0) + 1
        receiver_counts[r["target"]] = receiver_counts.get(r["target"], 0) + 1
    top_senders = sorted(sender_counts.items(), key=lambda kv: kv[1], reverse=True)[:3]
    top_receivers = sorted(receiver_counts.items(), key=lambda kv: kv[1], reverse=True)[:3]

    pathway_lines = []
    if "pathway" in rows[0]:
        pathway_counts = {}
        for r in rows:
            p = r.get("pathway") or "(unlabeled)"
            pathway_counts[p] = pathway_counts.get(p, 0) + 1
        top_pathways = sorted(pathway_counts.items(), key=lambda kv: kv[1], reverse=True)[:5]
        pathway_lines = [f"  - {p}: {n} interaction(s)" for p, n in top_pathways]

    # ---- Caveats (only emitted when they actually apply) ----
    caveats = [
        "This is a static edge-table quicklook, not CellChat/LIANA/CellPhoneDB's full statistical "
        "model (permutation testing, complex-aware ligand-receptor scoring); treat probability/p-value "
        "columns as whatever the upstream tool computed, not re-derived here.",
        "Cell-cell communication here is inferred from expression/interaction scores, not directly "
        "observed physical signaling; report it as inference, not proof (references/cns-method-guide.md).",
    ]
    n_low_p = sum(1 for r in rows if safe_float(r.get("pvalue"), 1.0) >= 0.05)
    if n_low_p:
        caveats.append(
            f"{n_low_p} of {len(rows)} edge(s) have p >= 0.05 in the input table; verify whether the "
            "upstream tool already filtered for significance before treating all rows as confirmed interactions."
        )

    conclusion = (
        f"This quicklook summarizes {len(rows)} ligand-receptor interaction(s) across {len(celltypes)} cell "
        f"group(s). The strongest interaction by probability is "
        f"{rows[0]['source']} -> {rows[0]['target']} via {rows[0]['ligand']}-{rows[0]['receptor']}"
        f" (prob={safe_float(rows[0].get('probability')):.2f}). "
        f"{top_senders[0][0]} is the top signal-sending cell type and {top_receivers[0][0]} is the top "
        f"signal-receiving cell type by interaction count. Treat these results as exploratory expression-based "
        "inference; see caveats below."
    )

    report_lines = [
        "# Cell-cell communication quicklook",
        "",
        "Generated from a ligand-receptor edge table (quicklook route, not a full CellChat/LIANA run).",
        "",
        "## Results",
        "",
        f"- Interactions: {len(rows)}",
        f"- Cell types/groups involved: {len(celltypes)}",
        "",
        "### Top interactions by probability",
        "",
    ] + [f"  - {line}" for line in top_pairs] + [
        "",
        "### Top signaling cell types",
        "",
        "  - Senders: " + ", ".join(f"{ct} ({n})" for ct, n in top_senders),
        "  - Receivers: " + ", ".join(f"{ct} ({n})" for ct, n in top_receivers),
    ]
    if pathway_lines:
        report_lines += ["", "### Top pathways", ""] + pathway_lines
    report_lines += ["", "### Caveats", ""] + [f"- {c}" for c in caveats]
    report_lines += ["", "### Conclusion", "", conclusion]
    (out / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"wrote communication outputs to {out}")


if __name__ == "__main__":
    main()
