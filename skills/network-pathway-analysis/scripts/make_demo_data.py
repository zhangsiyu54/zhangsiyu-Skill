#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo module-trait and pathway score tables.")
    parser.add_argument("--out", default="demo/network-pathway")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(27)
    samples = [f"C{i}" for i in range(1, 7)] + [f"T{i}" for i in range(1, 7)]
    with (out / "metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "group", "severity"])
        for i, sample in enumerate(samples):
            treat = sample.startswith("T")
            writer.writerow([sample, "Treat" if treat else "Control", round((2.2 if treat else 0.7) + random.random(), 3)])
    modules = ["MEblue", "MEred", "MEgreen", "MEpurple"]
    traits = ["group_Treat", "severity"]
    with (out / "module_trait.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["module", "trait", "correlation", "pvalue"])
        for module in modules:
            for trait in traits:
                val = random.uniform(-0.75, 0.85)
                if module == "MEred" and trait == "severity":
                    val = 0.82
                writer.writerow([module, trait, round(val, 3), round(10 ** (-abs(val) * 3.2), 5)])
    pathways = ["Interferon response", "Oxidative phosphorylation", "Cell cycle", "ECM remodeling", "T cell activation"]
    with (out / "pathway_scores.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pathway", *samples])
        for pathway in pathways:
            row = [pathway]
            for sample in samples:
                shift = 0.8 if sample.startswith("T") and pathway in {"Interferon response", "ECM remodeling"} else 0
                row.append(round(random.gauss(0, 0.45) + shift, 3))
            writer.writerow(row)
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
