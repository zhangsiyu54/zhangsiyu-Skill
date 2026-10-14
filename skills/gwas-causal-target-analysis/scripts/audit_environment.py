#!/usr/bin/env python3
"""Audit common GWAS command-line tools and Python/R packages."""
import argparse, importlib.util, json, platform, shutil, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

TOOLS = ["plink2", "bcftools", "tabix", "bgzip", "Rscript"]
PYTHON_PACKAGES = ["numpy", "pandas", "scipy", "statsmodels", "pysam"]
R_PACKAGES = ["data.table", "coloc", "susieR", "TwoSampleMR", "ieugwasr"]

def tool_version(path):
    for flag in ("--version", "-v"):
        try:
            p = subprocess.run([path, flag], text=True, capture_output=True, timeout=8)
            line = (p.stdout or p.stderr).strip().splitlines()
            if line:
                return line[0][:300]
        except (OSError, subprocess.TimeoutExpired):
            pass
    return None

def r_packages(rscript):
    if not rscript:
        return {x: False for x in R_PACKAGES}
    expr = "cat(paste(sapply(c(%s), requireNamespace, quietly=TRUE), collapse='\\n'))" % ",".join(repr(x) for x in R_PACKAGES)
    try:
        p = subprocess.run([rscript, "-e", expr], text=True, capture_output=True, timeout=30)
        vals = p.stdout.strip().splitlines()
        return dict(zip(R_PACKAGES, [x == "TRUE" for x in vals])) if len(vals) == len(R_PACKAGES) else {x: False for x in R_PACKAGES}
    except (OSError, subprocess.TimeoutExpired):
        return {x: False for x in R_PACKAGES}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    tools = {}
    for name in TOOLS:
        path = shutil.which(name)
        tools[name] = {"available": bool(path), "path": path, "version": tool_version(path) if path else None}
    report = {
        "schema_version": "1.0", "generated_at": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python": {"version": sys.version.split()[0], "executable": sys.executable},
        "tools": tools,
        "python_packages": {x: importlib.util.find_spec(x) is not None for x in PYTHON_PACKAGES},
        "r_packages": r_packages(tools["Rscript"]["path"]),
    }
    rp = report["r_packages"]
    report["capabilities"] = {
        "basic_validation": True,
        "ld_operations": tools["plink2"]["available"],
        "indexed_variant_io": tools["bcftools"]["available"] and tools["tabix"]["available"],
        "r_colocalization": tools["Rscript"]["available"] and rp["coloc"],
        "r_fine_mapping": tools["Rscript"]["available"] and rp["susieR"],
        "r_mr_core": tools["Rscript"]["available"] and rp["TwoSampleMR"],
        "r_opengwas_access": tools["Rscript"]["available"] and rp["ieugwasr"],
        "r_table_processing": tools["Rscript"]["available"] and rp["data.table"],
        "r_causal_workflow_complete": tools["Rscript"]["available"] and all(
            rp[x] for x in ("data.table", "coloc", "susieR", "TwoSampleMR", "ieugwasr")
        ),
    }
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0
if __name__ == "__main__": raise SystemExit(main())
