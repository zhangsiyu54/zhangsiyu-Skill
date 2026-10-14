#!/usr/bin/env python3
import argparse
import json
from datetime import datetime
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Write a reproducible methods document for a BioFlow analysis.")
    parser.add_argument("--out", default="output")
    parser.add_argument("--analysis-name", required=True)
    parser.add_argument("--route", default="python_quickplot", help="formal_r, python_quickplot, or setup_required")
    parser.add_argument("--inputs", default="", help="Comma-separated input files")
    parser.add_argument("--parameters", default="{}", help="JSON object or plain text parameter summary")
    parser.add_argument("--details", default="")
    parser.add_argument("--paper-style", default="", help="Paper/code practice basis, e.g. DESeq2/Seurat/CellChat style workflow")
    args = parser.parse_args()
    out = Path(args.out)
    (out / "run").mkdir(parents=True, exist_ok=True)
    try:
        params = json.loads(args.parameters)
    except Exception:
        params = {"summary": args.parameters}
    payload = {
        "analysis_name": args.analysis_name,
        "route": args.route,
        "inputs": [x.strip() for x in args.inputs.split(",") if x.strip()],
        "parameters": params,
        "details": args.details,
        "paper_style": args.paper_style,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (out / "run" / "methods.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    lines = [
        "# Methods and analysis details",
        "",
        f"## Analysis",
        "",
        args.analysis_name,
        "",
        "## Execution route",
        "",
        args.route,
        "",
        "## Input files",
        "",
    ]
    lines.extend([f"- `{x}`" for x in payload["inputs"]] or ["- Not recorded"])
    lines.extend(["", "## Parameters", ""])
    for key, value in payload["parameters"].items():
        lines.append(f"- `{key}`: {value}")
    if args.paper_style:
        lines.extend(["", "## Paper-style basis", "", args.paper_style])
    if args.details:
        lines.extend(["", "## Method details", "", args.details])
    lines.extend(
        [
            "",
            "## Reproducibility note",
            "",
            "This document records the route, inputs and parameters used for the analysis. "
            "For manuscript submission, add package versions, genome/database versions and exact statistical model terms from the final run logs.",
        ]
    )
    (out / "methods.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out / 'methods.md'}")


if __name__ == "__main__":
    main()
