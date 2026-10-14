#!/usr/bin/env python3
import argparse, csv, random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"


def random_cdr3(rng, length):
    # crude but plausible TCR-beta-style CDR3aa: starts with C, ends with F
    middle = "".join(rng.choice(AA) for _ in range(max(1, length - 2)))
    return f"C{middle}F"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="demo/repertoire")
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(87)
    vgenes = ["TRBV5-1", "TRBV7-2", "TRBV20-1", "TRBV6-5", "TRBV12-3"]
    jgenes = ["TRBJ1-1", "TRBJ1-2", "TRBJ2-1", "TRBJ2-3", "TRBJ2-7"]
    n_clonotypes = 25

    # Fix one CDR3/V/J identity per clonotype index so the *same* clonotype ID
    # always represents the same receptor sequence whenever it recurs across
    # samples. This matters for overlap/tracking analyses to be meaningful
    # (per cns-method-guide.md: clonotype = CDR3aa + V/J gene).
    identity = {}
    for i in range(1, n_clonotypes + 1):
        identity[i] = {
            "cdr3": random_cdr3(rng, rng.randint(10, 20)),
            "v_gene": rng.choice(vgenes),
            "j_gene": rng.choice(jgenes),
        }

    with (out / "clonotypes.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["sample", "group", "clonotype", "count", "v_gene", "j_gene", "cdr3", "celltype"])
        for sample in [f"C{i}" for i in range(1, 5)] + [f"T{i}" for i in range(1, 5)]:
            group = "Treat" if sample.startswith("T") else "Control"
            for i in range(1, n_clonotypes + 1):
                rate = 0.18 if group == "Treat" and i < 6 else 0.35
                count = max(1, int(rng.expovariate(rate)))
                idn = identity[i]
                w.writerow(
                    [
                        sample,
                        group,
                        f"clono_{i:03d}",
                        count,
                        idn["v_gene"],
                        idn["j_gene"],
                        idn["cdr3"],
                        rng.choice(["CD4_T", "CD8_T", "Treg"]),
                    ]
                )
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
