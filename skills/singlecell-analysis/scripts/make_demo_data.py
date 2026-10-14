#!/usr/bin/env python3
import argparse
import csv
import math
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo exported single-cell tables.")
    parser.add_argument("--out", default="demo/singlecell")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(11)

    celltypes = [("T_cell", -2.0, 0.6), ("B_cell", 1.6, 1.2), ("Myeloid", 0.2, -1.8)]
    samples = ["Sample1", "Sample2", "Sample3", "Sample4"]
    cells = []
    for idx in range(180):
        ct, cx, cy = celltypes[idx % len(celltypes)]
        cond = "Treat" if idx % 4 in (0, 1) else "Control"
        sample = samples[idx % len(samples)]
        x = random.gauss(cx, 0.42)
        y = random.gauss(cy, 0.38)
        # two sub-clusters per cell type so cluster labels are finer than celltype labels
        sub = idx % 2
        cluster = f"C{celltypes.index((ct, cx, cy)) * 2 + sub}"
        cell = f"Cell{idx+1:03d}"
        # simulate QC metrics with per-celltype variation and a couple of low-quality outliers
        ncount = max(300, int(random.gauss(6000, 1400)))
        nfeature = max(150, int(random.gauss(1800, 320)))
        mito_base = 3.0 if ct != "Myeloid" else 5.5
        percent_mt = max(0.0, round(random.gauss(mito_base, 1.6), 3))
        if idx % 37 == 0:
            # inject a rare low-quality / high-mito outlier cell for QC plot realism
            ncount = int(ncount * 0.15)
            nfeature = int(nfeature * 0.2)
            percent_mt = round(percent_mt + 18.0, 3)
        cells.append((cell, x, y, ct, cond, cluster, sample, ncount, nfeature, percent_mt))

    with (out / "cell_metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["cell", "UMAP_1", "UMAP_2", "celltype", "condition", "cluster", "sample", "nCount", "nFeature", "percent_mt"]
        )
        writer.writerows(cells)

    genes = ["CD3D", "MS4A1", "LYZ", "MKI67"]
    with (out / "marker_expression.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell", *genes])
        for cell, _, _, ct, *_ in cells:
            vals = []
            for gene in genes:
                base = 0.1
                if gene == "CD3D" and ct == "T_cell":
                    base = 2.4
                if gene == "MS4A1" and ct == "B_cell":
                    base = 2.2
                if gene == "LYZ" and ct == "Myeloid":
                    base = 2.5
                if gene == "MKI67":
                    base = 0.7 + 0.3 * math.sin(int(cell[-3:]) / 9)
                vals.append(round(max(0, random.gauss(base, 0.3)), 3))
            writer.writerow([cell, *vals])

    with (out / "communication_edges.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source", "target", "ligand", "receptor", "pathway", "probability", "pvalue"])
        rows = [
            ("T_cell", "Myeloid", "IFNG", "IFNGR1", "IFN-II", 0.82, 0.001),
            ("Myeloid", "T_cell", "IL15", "IL15RA", "Cytokine", 0.66, 0.006),
            ("B_cell", "T_cell", "CD40", "CD40LG", "Costimulation", 0.58, 0.014),
            ("Myeloid", "B_cell", "TNF", "TNFRSF1A", "TNF", 0.53, 0.022),
            ("T_cell", "B_cell", "IL21", "IL21R", "Cytokine", 0.47, 0.031),
            ("B_cell", "Myeloid", "CXCL12", "CXCR4", "Chemokine", 0.41, 0.045),
        ]
        writer.writerows(rows)

    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
