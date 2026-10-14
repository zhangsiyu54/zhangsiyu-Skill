#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo DEG and enrichment tables.")
    parser.add_argument("--out", default="demo/enrichment")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "differential_expression.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "log2FC", "padj"])
        for i in range(1, 41):
            fc = 2.2 if i <= 12 else -1.7 if i <= 24 else 0.3
            padj = 0.001 * i if i <= 24 else 0.2 + i / 300
            writer.writerow([f"Gene{i:03d}", fc, round(padj, 5)])

    # Gene membership per term, sized to match each term's "count" below, with
    # deliberate overlap between related terms so the gene-concept network
    # (cnetplot) demo has shared-gene edges to draw.
    gene_pools = {
        "Cell cycle": [f"Gene{i:03d}" for i in list(range(1, 13)) + list(range(30, 36))],  # 18
        "DNA replication": [f"Gene{i:03d}" for i in list(range(30, 40)) + list(range(64, 69))],  # 15
        "Immune response": [f"Gene{i:03d}" for i in list(range(13, 25)) + [69]],  # 13
        "Apoptotic signaling": [f"Gene{i:03d}" for i in list(range(20, 29)) + [70, 71]],  # 11
        "Extracellular matrix": [f"Gene{i:03d}" for i in list(range(25, 33)) + [72]],  # 9
    }

    with (out / "enrichment_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["term", "category", "gene_ratio", "pvalue", "padj", "count", "geneID"])
        rows = [
            ("Cell cycle", "KEGG", 0.32, 0.0004, 0.004, 18),
            ("DNA replication", "KEGG", 0.27, 0.001, 0.009, 15),
            ("Immune response", "GO_BP", 0.24, 0.003, 0.018, 13),
            ("Apoptotic signaling", "GO_BP", 0.21, 0.006, 0.028, 11),
            ("Extracellular matrix", "GO_BP", 0.18, 0.011, 0.041, 9),
        ]
        rows = [row + ("/".join(gene_pools[row[0]]),) for row in rows]
        writer.writerows(rows)
    with (out / "ranked_genes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene", "score"])
        for i in range(1, 101):
            score = 3.0 - i * 0.055
            if 36 <= i <= 55:
                score -= 1.1
            writer.writerow([f"Gene{i:03d}", round(score, 4)])
    with (out / "gene_sets.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["set", "gene"])
        for gene in [f"Gene{i:03d}" for i in list(range(1, 15)) + list(range(64, 70))]:
            writer.writerow(["Interferon response", gene])
        for gene in [f"Gene{i:03d}" for i in list(range(30, 45)) + list(range(80, 86))]:
            writer.writerow(["Cell cycle checkpoint", gene])
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
