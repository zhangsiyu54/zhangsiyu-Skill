#!/usr/bin/env python3
"""Create a transparent, conservative per-sample resource plan from a manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

GIB = 1024 ** 3


def delimiter(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def number(row: dict[str, str], key: str) -> int:
    value = (row.get(key) or "").strip()
    return int(value) if value else 0


def estimate(row: dict[str, str], headroom: float) -> dict:
    cells = number(row, "cells")
    nnz = number(row, "nonzero_entries")
    file_bytes = number(row, "file_size_bytes")
    # Three independent lower bounds: compressed input expansion, sparse working
    # matrices (48 bytes/nonzero), and per-cell metadata/embeddings (8 KiB/cell).
    working_gib = max(4.0, file_bytes * 3.0 / GIB, nnz * 48.0 / GIB, cells * 8192.0 / GIB)
    memory_gb = max(8, math.ceil((working_gib + 4.0) * headroom))
    cpus = 2 if cells < 50_000 else 4 if cells < 250_000 else 8
    hours = max(1, math.ceil((cells / 100_000 if cells else 1) + (file_bytes / (20 * GIB))))
    scratch_gb = max(10, math.ceil(max(file_bytes * 4.0 / GIB, nnz * 16.0 / GIB) * headroom))
    resource_class = "small" if memory_gb <= 32 else "medium" if memory_gb <= 128 else "large"
    return {
        "study_id": row.get("study_id", ""), "donor_id": row.get("donor_id", ""),
        "sample_id": row.get("sample_id", ""), "cells": cells, "nonzero_entries": nnz,
        "file_size_bytes": file_bytes, "cpus": cpus, "memory_gb": memory_gb,
        "wall_hours": hours, "scratch_gb": scratch_gb, "resource_class": resource_class,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", required=True, type=Path, help="JSON resource plan")
    parser.add_argument("--headroom", type=float, default=1.30)
    args = parser.parse_args()
    if args.headroom < 1.0:
        parser.error("--headroom must be at least 1.0")
    with args.manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter(args.manifest)))
    estimates = [estimate(row, args.headroom) for row in rows]
    plan = {
        "schema_version": 1,
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "manifest_row_count": len(rows),
        "estimation_status": "planning_only_not_a_benchmark",
        "headroom_multiplier": args.headroom,
        "assumptions": ["48 bytes per nonzero working-set lower bound", "8 KiB per cell lower bound", "3x compressed input expansion", "pilot benchmark required before production submission"],
        "samples": estimates,
        "array_envelope": {
            "cpus": max((x["cpus"] for x in estimates), default=2),
            "memory_gb": max((x["memory_gb"] for x in estimates), default=8),
            "wall_hours": max((x["wall_hours"] for x in estimates), default=1),
            "scratch_gb": max((x["scratch_gb"] for x in estimates), default=10),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(plan["array_envelope"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
