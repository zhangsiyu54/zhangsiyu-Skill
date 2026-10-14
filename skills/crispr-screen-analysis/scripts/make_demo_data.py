#!/usr/bin/env python3
"""Create demo CRISPR screen data.

Writes two files that together cover both accepted input types from
`references/input-contract.md`:

- `gene_summary.csv`: MAGeCK-style gene summary (gene, logFC, pvalue, fdr).
- `sgrna_counts.csv`: sgRNA-level raw count table (sgRNA, gene, 2 control
  replicates, 2 treatment replicates), used by `quickplot_crispr.py` for the
  QC (read depth / sgRNA distribution / replicate correlation / Gini index)
  and sgRNA-consistency figures. Pass --no-sgrna to skip it.

Both files are derived from the same simulated per-guide counts (gene_summary's
logFC is the single most extreme guide per gene, mirroring how MAGeCK RRA's
rank aggregation can be swayed by one strongly-responding guide), so a gene's
gene-level call is always internally consistent with its guide-level pattern.
Gene009 is deliberately seeded with only 1 of 5 active guides: its gene-level
score is carried entirely by that guide, so it reliably reaches significance
in gene_summary.csv while the sgRNA consistency plot/report caveat correctly
flags it as single-guide-driven (references/pitfalls.md #6). With --no-sgrna,
gene_summary.csv is still written (from the same simulated guide data;
sgrna_counts.csv is simply not saved to disk).
"""
import argparse
import csv
import math
import random
from pathlib import Path

GENES = 100
POS_HITS = range(1, 9)  # Gene001-Gene008: positive selection (enriched in treatment)
NEG_HITS = range(9, 17)  # Gene009-Gene016: negative selection (depleted in treatment)
GUIDES_PER_GENE = 5
NTC_GUIDES = 60
PSEUDOCOUNT = 5.0  # matches quickplot_crispr.py's plot_sgrna_consistency pseudocount


def gene_effect(i):
    if i in POS_HITS:
        return 1.8
    if i in NEG_HITS:
        return -1.6
    return 0.0


def simulate_guides(rng):
    """Simulate per-guide control/treatment counts for every gene plus NTC
    guides, and the per-guide log2FC (same pseudocount formula the analysis
    script uses) needed to derive a consistent gene-level score."""
    sgrna_rows = []
    per_gene_guide_lfc = {}
    for i in range(1, GENES + 1):
        gene = f"Gene{i:03d}"
        effect = gene_effect(i)
        guide_lfcs = []
        for g in range(1, GUIDES_PER_GENE + 1):
            sg_id = f"{gene}_sg{g}"
            # Gene009: only guide 1 is active -> single-guide-driven hit demo.
            efficiency = 0.05 if (gene == "Gene009" and g > 1) else 1.0
            baseline = max(20.0, rng.lognormvariate(5.9, 0.35))
            c1 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
            c2 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
            fold = 2 ** (effect * efficiency)
            t_mean = baseline * fold
            t1 = max(0, round(rng.gauss(t_mean, max(t_mean * 0.15, 3))))
            t2 = max(0, round(rng.gauss(t_mean, max(t_mean * 0.15, 3))))
            sgrna_rows.append([sg_id, gene, c1, c2, t1, t2])
            c_mean, t_mean_obs = (c1 + c2) / 2.0, (t1 + t2) / 2.0
            guide_lfcs.append(math.log2((t_mean_obs + PSEUDOCOUNT) / (c_mean + PSEUDOCOUNT)))
        per_gene_guide_lfc[gene] = guide_lfcs
    for n in range(1, NTC_GUIDES + 1):
        sg_id = f"NTC_sg{n:03d}"
        baseline = max(20.0, rng.lognormvariate(5.9, 0.35))
        c1 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
        c2 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
        t1 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
        t2 = max(0, round(rng.gauss(baseline, baseline * 0.12)))
        sgrna_rows.append([sg_id, "NTC", c1, c2, t1, t2])
    return sgrna_rows, per_gene_guide_lfc


def write_gene_summary(out, per_gene_guide_lfc):
    path = out / "gene_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["gene", "logFC", "pvalue", "fdr"])
        for i in range(1, GENES + 1):
            gene = f"Gene{i:03d}"
            # Most extreme guide, not the mean: MAGeCK RRA-style rank aggregation
            # can be swayed by one strongly-responding guide (pitfalls.md #6),
            # so this keeps gene_summary.csv consistent with sgrna_counts.csv.
            fc = max(per_gene_guide_lfc[gene], key=abs)
            pval = min(1, 10 ** (-abs(fc) * 1.7))
            fdr = min(1, pval * 8)
            w.writerow([gene, round(fc, 3), round(pval, 6), round(fdr, 6)])
    return path


def write_sgrna_counts(out, sgrna_rows):
    path = out / "sgrna_counts.csv"
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["sgRNA", "gene", "control_rep1", "control_rep2", "treatment_rep1", "treatment_rep2"])
        w.writerows(sgrna_rows)
    return path


def main():
    p = argparse.ArgumentParser(description="Create demo CRISPR gene summary and sgRNA count table.")
    p.add_argument("--out", default="demo/crispr")
    p.add_argument("--no-sgrna", action="store_true", help="skip writing the sgRNA-level count table")
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(97)
    sgrna_rows, per_gene_guide_lfc = simulate_guides(rng)
    gs = write_gene_summary(out, per_gene_guide_lfc)
    print(f"wrote {gs}")
    if not a.no_sgrna:
        sg = write_sgrna_counts(out, sgrna_rows)
        print(f"wrote {sg}")
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
