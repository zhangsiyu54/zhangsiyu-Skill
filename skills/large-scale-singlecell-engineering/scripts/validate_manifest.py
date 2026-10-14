#!/usr/bin/env python3
"""Validate a study/donor/sample manifest using only the Python standard library."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED = ("study_id", "donor_id", "sample_id", "input_path")
NONNEGATIVE_INTS = ("cells", "nonzero_entries", "file_size_bytes")


def dialect(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def clean(value: object) -> str:
    return str(value or "").strip()


def validate(path: Path, check_paths: bool = False) -> dict:
    errors: list[str] = []
    warnings: list[str] = []
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=dialect(path))
        fields = [clean(x) for x in (reader.fieldnames or [])]
        missing = [field for field in REQUIRED if field not in fields]
        if missing:
            return {"valid": False, "errors": [f"missing required columns: {', '.join(missing)}"], "warnings": [], "summary": {"rows": 0}}
        for line_number, raw in enumerate(reader, start=2):
            row = {clean(k): clean(v) for k, v in raw.items() if k is not None}
            rows.append(row)
            for field in REQUIRED:
                if not row.get(field):
                    errors.append(f"line {line_number}: empty {field}")
            for field in NONNEGATIVE_INTS:
                value = row.get(field, "")
                if value:
                    try:
                        if int(value) < 0:
                            raise ValueError
                    except ValueError:
                        errors.append(f"line {line_number}: {field} must be a non-negative integer")
            input_path = row.get("input_path", "")
            if input_path and check_paths and not Path(os.path.expandvars(os.path.expanduser(input_path))).exists():
                errors.append(f"line {line_number}: input_path does not exist: {input_path}")

    sample_counts = Counter(row.get("sample_id", "") for row in rows if row.get("sample_id"))
    for sample_id, count in sorted(sample_counts.items()):
        if count > 1:
            errors.append(f"sample_id is not unique: {sample_id} ({count} rows)")
    donor_studies: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row.get("donor_id") and row.get("study_id"):
            donor_studies[row["donor_id"]].add(row["study_id"])
    for donor_id, studies in sorted(donor_studies.items()):
        if len(studies) > 1:
            warnings.append(f"donor_id {donor_id} occurs in multiple studies; use globally unique donor IDs if these are different people")
    if not rows:
        errors.append("manifest has no data rows")
    if not any(row.get("cells") or row.get("nonzero_entries") or row.get("file_size_bytes") for row in rows):
        warnings.append("no sizing fields found; resource estimates will use conservative defaults")
    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "rows": len(rows),
            "studies": len({r.get("study_id") for r in rows if r.get("study_id")}),
            "donors": len({(r.get("study_id"), r.get("donor_id")) for r in rows if r.get("donor_id")}),
            "samples": len(sample_counts),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--check-paths", action="store_true", help="require every input_path to exist locally")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    try:
        report = validate(args.manifest, args.check_paths)
    except (OSError, csv.Error) as exc:
        report = {"valid": False, "errors": [str(exc)], "warnings": [], "summary": {"rows": 0}}
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered + "\n", encoding="utf-8")
    return 0 if report["valid"] else 2


if __name__ == "__main__":
    sys.exit(main())
