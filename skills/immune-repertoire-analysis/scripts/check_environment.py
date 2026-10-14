#!/usr/bin/env python3
import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path


def has_python_module(name):
    return importlib.util.find_spec(name) is not None


def r_packages(packages):
    rscript = shutil.which("Rscript")
    if not rscript or not packages:
        return {}
    expr = "cat(paste(%s, vapply(%s, requireNamespace, logical(1), quietly=TRUE), sep='=', collapse='\\n'))" % (
        "c(" + ",".join(repr(p) for p in packages) + ")",
        "c(" + ",".join(repr(p) for p in packages) + ")",
    )
    try:
        proc = subprocess.run([rscript, "-e", expr], text=True, capture_output=True, check=False, timeout=20)
        out = {}
        for line in proc.stdout.splitlines():
            if "=" in line:
                key, val = line.split("=", 1)
                out[key] = val == "TRUE"
        return out
    except Exception:
        return {p: False for p in packages}


def main():
    parser = argparse.ArgumentParser(description="Check R/Python availability for BioFlow skills.")
    parser.add_argument("--out", default="output")
    parser.add_argument("--r-packages", default="", help="Comma-separated R packages to check")
    parser.add_argument("--python-packages", default="matplotlib,numpy", help="Comma-separated Python packages to check")
    args = parser.parse_args()

    out = Path(args.out)
    (out / "run").mkdir(parents=True, exist_ok=True)
    py_packages = [p.strip() for p in args.python_packages.split(",") if p.strip()]
    r_pkgs = [p.strip() for p in args.r_packages.split(",") if p.strip()]
    status = {
        "python": sys.executable,
        "python_version": sys.version.split()[0],
        "r": shutil.which("R"),
        "rscript": shutil.which("Rscript"),
        "python_packages": {p: has_python_module(p) for p in py_packages},
        "r_packages": r_packages(r_pkgs),
        "r_package_check_supplied": bool(r_pkgs),
    }
    status["can_run_r_workflow"] = bool(status["rscript"]) and bool(r_pkgs) and all(status["r_packages"].values())
    status["can_run_python_quickplot"] = all(status["python_packages"].values())
    status["warnings"] = []
    if status["rscript"] and not r_pkgs:
        status["warnings"].append(
            "Rscript is available, but no analysis-specific R packages were checked; formal R workflow readiness is unknown. Re-run with --r-packages from references/r-dependencies.md."
        )
    if status["can_run_r_workflow"]:
        status["recommended_route"] = "formal_r"
    elif status["can_run_python_quickplot"]:
        status["recommended_route"] = "python_quickplot"
    else:
        status["recommended_route"] = "setup_required"

    (out / "run" / "environment.json").write_text(json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8")
    if status["recommended_route"] == "setup_required":
        (out / "setup_instructions.md").write_text(
            "# Environment setup required\n\n"
            "Install Python plotting dependencies:\n\n"
            "```bash\npython3 -m pip install matplotlib numpy\n```\n\n"
            "Install R with one of these routes:\n\n"
            "```bash\n# macOS Homebrew\nbrew install r\n\n"
            "# macOS CRAN package when you can use the installer UI\n"
            "open https://cran.r-project.org/bin/macosx/\n\n"
            "# macOS portable R when you do not have administrator permission\n"
            "mkdir -p ~/.local/share/portable-r ~/.local/bin\n"
            "# Download the arm64 or x86_64 portable-r archive, unpack it, then link R/Rscript into ~/.local/bin\n"
            "export PATH=\"$HOME/.local/bin:$PATH\"\n"
            "```\n\n"
            "Then install analysis-specific R packages listed in the skill references.\n",
            encoding="utf-8",
        )
    print(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    main()
