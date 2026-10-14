#!/usr/bin/env python3
"""Apply lightweight, consistent publication framing to already-rendered figure
files (PNG/PDF) and export a standardized PNG (raster preview) + PDF (vector/print)
pair for each input, plus a batch QC report.

SCOPE AND HARD LIMITATION (read before using this script on real submissions):
    This script re-frames EXISTING image files. It adds a panel label, standardizes
    surrounding typography/margins/background, and re-exports PNG+PDF. It does
    **not** re-render the underlying plot from source data. Any problem baked into
    the original pixels -- overplotting/overlapping points, a rainbow/jet colormap,
    a red-green diverging scheme, a truncated y-axis, undefined error bars, etc. --
    CANNOT be fixed by this script, because the pixels (or, for PDF, the vector
    drawing commands of the original plot) are copied/re-embedded as-is. Fixing
    those requires regenerating the figure from the source table/values with a
    colorblind-safe palette and density-aware rendering (see
    references/pitfalls.md and references/cns-method-guide.md).

    What this script CAN do about that gap: `analyze_source_image()` below runs a
    best-effort pixel heuristic scan for a subset of pitfalls (low print resolution,
    rainbow/jet-like hue distribution, red-green-dominant diverging color use) and
    reports them in `figure_qc_report.tsv` and on stdout, so a human/agent can
    decide whether a panel needs to be replotted from source instead of polished.
    This is a heuristic, not a guarantee -- it can both miss real issues and flag
    false positives, and it cannot check anything that requires the source data or
    caption text (error-bar definition, units, axis truncation, gene/label choice).
"""
import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
except Exception as exc:
    raise SystemExit(
        "polish_figures.py requires matplotlib to export PNG/PDF. "
        "Install it with: python3 -m pip install matplotlib"
    ) from exc


# ---------------------------------------------------------------------------
# Single source of truth for colors: load assets/palette.json so every script
# (and references/visual-style.md, which documents the same hex values) shares
# one category-to-color mapping. See references/pitfalls.md #3 (inconsistent
# color-to-category mapping across panels/figures) -- previously this script
# hardcoded its own text/edge colors that did not match assets/palette.json.
# ---------------------------------------------------------------------------
PALETTE_PATH = Path(__file__).resolve().parent.parent / "assets" / "palette.json"
DEFAULT_PALETTE = {
    "blue": "#3B82C4",
    "red": "#D84A4A",
    "green": "#3FA37A",
    "amber": "#D99A2B",
    "purple": "#7C65C1",
    "text": "#1F2937",
    "muted": "#667085",
    "grid": "#E5E7EB",
    "diverging_low": "#3B4CC0",
    "diverging_high": "#B40426",
    "sequential_colormap": "viridis",
    "max_categorical_hues": 10,
}


def load_palette():
    try:
        data = json.loads(PALETTE_PATH.read_text(encoding="utf-8"))
        merged = dict(DEFAULT_PALETTE)
        merged.update(data)
        return merged
    except Exception:
        return dict(DEFAULT_PALETTE)


PALETTE = load_palette()
TEXT = PALETTE["text"]
MUTED = PALETTE["muted"]
GRID = PALETTE["grid"]

MIN_PRINT_DPI = 300  # journal raster minimum at final physical size (references/cns-method-guide.md)
SINGLE_COLUMN_MM = 89.0  # standard single-column journal width, used only for the resolution heuristic


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 16,
            "axes.titlesize": 24,
            "axes.labelsize": 21,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 15,
            "axes.edgecolor": GRID,
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "text.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_png_pdf(fig, output_stem):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Pixel-heuristic QC scan (references/pitfalls.md #1, #2, #4)
# ---------------------------------------------------------------------------
def analyze_source_image(image, skip_resolution_check=False):
    """Best-effort heuristic scan of a source raster for a subset of documented
    pitfalls. Returns a dict of qc fields plus a list of human-readable warnings.
    See the module docstring for what this can and cannot detect.
    """
    height, width = image.shape[:2]
    warnings = []

    # --- Pitfall #4: font size too small at final print size -----------------
    # We cannot measure a baked-in font's point size from raster pixels without
    # OCR, so we use effective print resolution at a standard single journal
    # column width as a documented proxy: low effective DPI at final physical
    # size strongly correlates with illegible small text once the figure is
    # placed at its true submission size.
    if skip_resolution_check:
        effective_dpi = None
        low_resolution_risk = False
    else:
        effective_dpi = width / (SINGLE_COLUMN_MM / 25.4)
        low_resolution_risk = effective_dpi < MIN_PRINT_DPI
        if low_resolution_risk:
            warnings.append(
                f"low_resolution_risk: ~{effective_dpi:.0f} effective dpi at a "
                f"{SINGLE_COLUMN_MM:.0f}mm single column (<{MIN_PRINT_DPI} dpi target); "
                "baked-in text/lines may be illegible at final print size -- consider "
                "re-exporting or replotting from source at higher resolution."
            )

    rgb = image[..., :3].astype(float)
    if rgb.max() > 1.0:
        rgb = rgb / 255.0
    flat = rgb.reshape(-1, 3)
    luma = flat @ np.array([0.299, 0.587, 0.114])
    saturation = flat.max(axis=1) - flat.min(axis=1)
    # ignore near-white/near-black/low-saturation background pixels, which dominate
    # most figures (white canvas, black axis lines/text)
    mask = (luma > 0.05) & (luma < 0.97) & (saturation > 0.12)
    sample = flat[mask]
    if len(sample) > 4000:
        idx = np.random.RandomState(0).choice(len(sample), 4000, replace=False)
        sample = sample[idx]

    suspected_rainbow = False
    suspected_red_green = False
    if len(sample) > 200:
        hue = mcolors.rgb_to_hsv(sample)[:, 0] * 360.0
        red_mass = float(np.mean((hue < 15) | (hue > 345)))
        yellow_mass = float(np.mean((hue >= 45) & (hue < 75)))
        green_mass = float(np.mean((hue >= 90) & (hue < 160)))
        cyan_mass = float(np.mean((hue >= 170) & (hue < 200)))
        blue_mass = float(np.mean((hue >= 200) & (hue < 260)))

        # Pitfall #2: rainbow/jet spans blue->cyan->green->yellow->red with
        # roughly comparable mass in every band. A typical figure with a
        # handful of categorical/diverging colors will not populate all five
        # bands substantially at once.
        if min(red_mass, yellow_mass, green_mass, cyan_mass, blue_mass) > 0.05:
            suspected_rainbow = True
            warnings.append(
                "suspected_rainbow_or_jet_colormap: pixel hues populate red/yellow/"
                "green/cyan/blue bands roughly evenly, consistent with a jet/rainbow "
                "gradient. If this is a heatmap or other continuous-data panel, "
                f"replot from source with a perceptually uniform colormap "
                f"(matplotlib '{PALETTE.get('sequential_colormap', 'viridis')}') instead."
            )

        # Pitfall #1: red-green diverging encoding -- both red AND green present
        # with little/no blue is consistent with a red/green (rather than
        # blue/red) diverging scale. Thresholds are deliberately low: real
        # diverging data is often normally distributed, so only the extreme tail
        # actually reaches saturated red/green while most cells sit near the
        # (often yellow) midpoint -- calibrated against a synthetic RdYlGn
        # heatmap of N(0,1) data, which put red/green mass at ~0.03/0.05.
        if red_mass > 0.025 and green_mass > 0.025 and blue_mass < 0.08:
            suspected_red_green = True
            warnings.append(
                "suspected_red_green_diverging: substantial red and green pixel mass "
                "with little blue, consistent with a red-green up/down color scheme "
                "that is unreadable under red-green color blindness. If this encodes "
                "direction, replot with a blue-red diverging scale "
                f"({PALETTE.get('diverging_low')}/{PALETTE.get('diverging_high')})."
            )

    return {
        "width_px": width,
        "height_px": height,
        "effective_dpi_at_89mm": round(effective_dpi, 1) if effective_dpi is not None else "n/a_vector_source",
        "low_resolution_risk": low_resolution_risk,
        "suspected_rainbow_or_jet_colormap": suspected_rainbow,
        "suspected_red_green_diverging": suspected_red_green,
        "warnings": warnings,
    }


def _build_panel_figure(image, panel_label=""):
    apply_style()
    height, width = image.shape[:2]
    fig_w = max(5.0, min(12.0, width / 220))
    fig_h = max(4.0, min(10.0, height / 220))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(image)
    ax.axis("off")
    if panel_label:
        ax.text(
            0.015,
            0.98,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=20,
            fontweight="bold",
            color=TEXT,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        )
    return fig


def polish_png(source, output_stem, panel_label=""):
    source = Path(source)
    try:
        image = plt.imread(source)
    except Exception as exc:
        raise ValueError(
            f"could not decode {source} as an image -- the file may be corrupt, truncated, or not a "
            "real PNG despite its extension"
        ) from exc
    if image is None or getattr(image, "size", 0) == 0 or image.ndim < 2 or 0 in image.shape[:2]:
        raise ValueError(f"{source} decoded to an empty or invalid image (0-size array)")
    qc = analyze_source_image(image)
    fig = _build_panel_figure(image, panel_label)
    save_png_pdf(fig, output_stem)
    return qc


def render_pdf_page_to_png(source, target_png):
    target_png = Path(target_png)
    target_png.parent.mkdir(parents=True, exist_ok=True)
    try:
        import fitz

        document = fitz.open(str(source))
        page = document[0]
        pixmap = page.get_pixmap(matrix=fitz.Matrix(2.2, 2.2), alpha=False)
        pixmap.save(str(target_png))
        return True
    except Exception:
        pass
    commands = [
        ["pdftoppm", "-png", "-r", "220", "-singlefile", str(source), str(target_png.with_suffix(""))],
        ["magick", "-density", "220", f"{source}[0]", str(target_png)],
        ["sips", "-s", "format", "png", str(source), "--out", str(target_png)],
    ]
    for command in commands:
        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if target_png.exists():
                return True
        except Exception:
            continue
    return False


def write_pdf_preview_placeholder(source, target_png):
    apply_style()
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    ax.axis("off")
    ax.text(0.5, 0.58, Path(source).name, ha="center", va="center", fontsize=22, fontweight="bold", color=TEXT)
    ax.text(
        0.5,
        0.42,
        "PDF copied. Install PyMuPDF, Poppler or ImageMagick to rasterize PDF previews.",
        ha="center",
        va="center",
        fontsize=14,
        color=MUTED,
        wrap=True,
    )
    fig.savefig(target_png, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def polish_pdf(source, output_stem, panel_label=""):
    source = Path(source)
    output_stem = Path(output_stem)
    if source.stat().st_size == 0:
        raise ValueError(f"{source} is an empty (0-byte) file, not a valid PDF")
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    # Preserve the true vector PDF untouched -- this is the print/journal-submission
    # artifact (references/pitfalls.md #5: raster where vector is required). Only a
    # PNG preview is (re)generated below; it must never be allowed to overwrite this
    # copy with a rasterized re-save.
    shutil.copy2(source, output_stem.with_suffix(".pdf"))
    qc = {
        "width_px": "n/a_vector_source",
        "height_px": "n/a_vector_source",
        "effective_dpi_at_89mm": "n/a_vector_source",
        "low_resolution_risk": False,
        "suspected_rainbow_or_jet_colormap": False,
        "suspected_red_green_diverging": False,
        "warnings": [],
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        preview = Path(tmpdir) / "preview.png"
        if render_pdf_page_to_png(source, preview):
            image = plt.imread(preview)
            qc = analyze_source_image(image, skip_resolution_check=True)
            fig = _build_panel_figure(image, panel_label)
            # PNG preview only -- do NOT call save_png_pdf here, it would overwrite
            # the vector .pdf copied above with a rasterized re-save.
            fig.savefig(output_stem.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
            plt.close(fig)
        else:
            write_pdf_preview_placeholder(source, output_stem.with_suffix(".png"))
    return qc


def polish_file(source, output_stem, panel_label=""):
    source = Path(source)
    output_stem = Path(output_stem)
    if not source.exists():
        raise FileNotFoundError(
            f"source figure not found: {source} (check the source_file path in the manifest -- "
            "it is resolved relative to the current working directory)"
        )
    if source.is_dir():
        raise ValueError(f"{source} is a directory, not a figure file")
    if source.suffix.lower() == ".png":
        qc = polish_png(source, output_stem, panel_label)
        return ["png", "pdf"], qc
    if source.suffix.lower() == ".pdf":
        qc = polish_pdf(source, output_stem, panel_label)
        return ["png", "pdf"], qc
    raise ValueError(
        f"Unsupported source format: {source} (suffix '{source.suffix}'). Only .png and .pdf are supported."
    )


# ---------------------------------------------------------------------------
# Panel-label consistency check (references/pitfalls.md #9)
# ---------------------------------------------------------------------------
def validate_panel_labels(rows):
    labels = [row.get("panel_label", "").strip() for row in rows]
    warnings = []
    missing = [row.get("source_file", "?") for row, label in zip(rows, labels) if not label]
    if missing:
        warnings.append(f"panel_label missing for {len(missing)} row(s): {', '.join(missing)}")
    non_empty = [label for label in labels if label]
    if len(non_empty) != len(set(non_empty)):
        warnings.append("duplicate panel_label values found in manifest -- every panel must be uniquely labeled")
    expected = [chr(ord("A") + i) for i in range(len(non_empty))]
    if non_empty and non_empty != expected:
        warnings.append(
            f"panel_label sequence is {non_empty} but the standard consistent scheme is sequential "
            f"A, B, C... ({''.join(expected)}); confirm this is an intentional custom order "
            "(references/pitfalls.md #9)"
        )
    return warnings


# ---------------------------------------------------------------------------
# Palette colorblind-safety check (references/pitfalls.md #1, #3, #7)
# ---------------------------------------------------------------------------
def _hex_to_rgb(value):
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def _simulate_deuteranopia(rgb):
    """Rough linear approximation of deuteranopia (red-green color blindness) for a
    quick automated pre-check. This is a simplified heuristic matrix, not a
    clinically validated simulator (e.g. Coblis/Color Oracle) -- always confirm the
    final palette with a real CVD simulator as recommended in
    references/cns-method-guide.md before submission.
    """
    matrix = np.array(
        [
            [0.625, 0.375, 0.0],
            [0.700, 0.300, 0.0],
            [0.000, 0.300, 0.7],
        ]
    )
    return rgb @ matrix.T


DATA_ENCODING_KEYS = {"blue", "red", "green", "amber", "purple", "diverging_low", "diverging_high"}


def check_palette_colorblind_safe(palette=None):
    palette = palette or PALETTE
    # Only compare colors that are actually used to encode data categories/direction.
    # UI/typography colors (text, muted, grid) are never used as a category-vs-category
    # or direction encoding, so comparing them against data colors produces noise.
    named = [
        (k, v)
        for k, v in palette.items()
        if isinstance(v, str) and v.startswith("#") and k in DATA_ENCODING_KEYS
    ]
    problems = []
    if len(named) >= 2:
        rgb = np.array([_hex_to_rgb(v) for _, v in named])
        sim = _simulate_deuteranopia(rgb)
        for i in range(len(named)):
            for j in range(i + 1, len(named)):
                dist = float(np.linalg.norm(sim[i] - sim[j]))
                if dist < 0.12:
                    problems.append(
                        f"{named[i][0]} ({named[i][1]}) and {named[j][0]} ({named[j][1]}) look similar "
                        f"under a simulated deuteranopia pass (delta={dist:.3f}); avoid using both to "
                        "distinguish categories/direction in the same panel"
                    )
    diverging_low = palette.get("diverging_low")
    diverging_high = palette.get("diverging_high")
    if diverging_low and diverging_high:
        low_hue = mcolors.rgb_to_hsv(np.array([_hex_to_rgb(diverging_low)]))[0, 0] * 360
        high_hue = mcolors.rgb_to_hsv(np.array([_hex_to_rgb(diverging_high)]))[0, 0] * 360
        if not (180 <= low_hue <= 280):
            problems.append(f"diverging_low hue {low_hue:.0f} deg is not in the blue range (180-280 deg); pick a blue-family low color")
        if not (330 <= high_hue <= 360 or 0 <= high_hue <= 20):
            problems.append(f"diverging_high hue {high_hue:.0f} deg is not in the red range (330-360/0-20 deg); pick a red-family high color")
    max_hues = palette.get("max_categorical_hues")
    categorical_keys = {"blue", "red", "green", "amber", "purple"}
    categorical_count = len([k for k, v in named if k in categorical_keys])
    if max_hues and categorical_count > max_hues:
        problems.append(
            f"{categorical_count} categorical colors defined but max_categorical_hues is {max_hues} "
            "(references/pitfalls.md #12); add a redundant channel (shape/text/facet) instead of more hues"
        )
    return problems


def write_qc_report(qc_rows, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "source_file",
                "output_stem",
                "panel_label",
                "width_px",
                "height_px",
                "effective_dpi_at_89mm",
                "low_resolution_risk",
                "suspected_rainbow_or_jet_colormap",
                "suspected_red_green_diverging",
                "notes",
            ]
        )
        for row in qc_rows:
            writer.writerow(
                [
                    row["source_file"],
                    row["output_stem"],
                    row.get("panel_label", ""),
                    row["qc"]["width_px"],
                    row["qc"]["height_px"],
                    row["qc"]["effective_dpi_at_89mm"],
                    row["qc"]["low_resolution_risk"],
                    row["qc"]["suspected_rainbow_or_jet_colormap"],
                    row["qc"]["suspected_red_green_diverging"],
                    "; ".join(row["qc"]["warnings"]) or "none",
                ]
            )
    return out_path


# ---------------------------------------------------------------------------
# Human-readable batch report (references/deliverables.md: output/report.md)
# ---------------------------------------------------------------------------
def _error_qc(message):
    return {
        "width_px": "error",
        "height_px": "error",
        "effective_dpi_at_89mm": "error",
        "low_resolution_risk": False,
        "suspected_rainbow_or_jet_colormap": False,
        "suspected_red_green_diverging": False,
        "warnings": [f"PROCESSING_ERROR: {message}"],
        "processing_error": True,
        "error_message": message,
    }


def _flag_reasons(qc):
    reasons = []
    if qc.get("low_resolution_risk"):
        reasons.append("low_resolution_risk")
    if qc.get("suspected_rainbow_or_jet_colormap"):
        reasons.append("suspected_rainbow_or_jet_colormap")
    if qc.get("suspected_red_green_diverging"):
        reasons.append("suspected_red_green_diverging")
    return reasons


def _build_conclusion_sentence(total, processed, errored, n_low_res, n_rainbow, n_redgreen):
    parts = [f"{len(processed)}/{total} figures polished"]
    if errored:
        parts.append(f"{len(errored)} failed to process (see Processing errors)")
    colormap_flags = n_rainbow + n_redgreen
    parts.append(f"{colormap_flags} flagged for colormap issues (rainbow/jet or red-green diverging)")
    if n_low_res == 0:
        parts.append(f"0 figures below the {MIN_PRINT_DPI} dpi print-resolution target")
    elif n_low_res == 1:
        low_res_row = next(r for r in processed if r["qc"].get("low_resolution_risk"))
        dpi = low_res_row["qc"].get("effective_dpi_at_89mm")
        name = low_res_row["source_file"]
        parts.append(
            f"1 figure ({name}) has a low estimated print DPI (~{dpi}) at column width -- consider "
            "requesting the original vector source rather than a rasterized PNG"
        )
    else:
        parts.append(
            f"{n_low_res} figures have a low estimated print DPI at column width -- consider requesting "
            "original vector sources rather than rasterized PNGs"
        )
    return "; ".join(parts) + "."


def write_batch_report(qc_rows, out_path, manifest_path=None):
    """Write a plain-language Markdown summary of a --manifest batch run: how many
    figures were processed, how many were flagged by which heuristic, what to do
    about each flagged figure (naming it explicitly), and any processing errors.
    This is the human-facing counterpart to the machine-readable figure_qc_report.tsv.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total = len(qc_rows)
    errored = [r for r in qc_rows if r["qc"].get("processing_error")]
    processed = [r for r in qc_rows if not r["qc"].get("processing_error")]
    flagged = [r for r in processed if _flag_reasons(r["qc"])]
    n_low_res = sum(1 for r in processed if r["qc"].get("low_resolution_risk"))
    n_rainbow = sum(1 for r in processed if r["qc"].get("suspected_rainbow_or_jet_colormap"))
    n_redgreen = sum(1 for r in processed if r["qc"].get("suspected_red_green_diverging"))

    lines = ["# Figure Polish Report", ""]
    source_note = f" from `{manifest_path}`" if manifest_path else ""
    lines.append(f"Generated by `scripts/polish_figures.py --manifest`{source_note}.")
    lines.append("")

    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Figures in manifest: {total}")
    lines.append(f"- Successfully polished (PNG+PDF exported): {len(processed)}")
    if errored:
        lines.append(f"- Failed to process: {len(errored)} (see Processing errors below)")
    lines.append(f"- Flagged by QC heuristics: {len(flagged)} / {len(processed)}")
    lines.append(f"  - Suspected rainbow/jet colormap: {n_rainbow}")
    lines.append(f"  - Suspected red-green diverging use: {n_redgreen}")
    lines.append(
        f"  - Low estimated print resolution (<{MIN_PRINT_DPI} dpi at {SINGLE_COLUMN_MM:.0f}mm column): {n_low_res}"
    )
    lines.append("")

    lines.append("## Flagged figures (action recommended)")
    lines.append("")
    if flagged:
        for row in flagged:
            qc = row["qc"]
            panel = row.get("panel_label") or "?"
            lines.append(f"### `{row['source_file']}` (panel {panel})")
            lines.append("")
            lines.append(f"- Flagged for: {', '.join(_flag_reasons(qc))}")
            for warning in qc.get("warnings", []):
                lines.append(f"- {warning}")
            lines.append("")
    else:
        lines.append(
            "None of the processed figures triggered a QC heuristic flag. This does not certify the "
            "figures as publication-ready -- see Known limitation below; unflagged issues (overplotting, "
            "truncated axes, undefined error bars, unit/precision mismatches) require manual review against "
            "`references/pitfalls.md`."
        )
        lines.append("")

    if errored:
        lines.append("## Processing errors")
        lines.append("")
        lines.append("These figures were skipped and were NOT polished; fix the underlying file/path and re-run.")
        lines.append("")
        for row in errored:
            lines.append(f"- `{row['source_file']}`: {row['qc'].get('error_message', 'unknown error')}")
        lines.append("")

    lines.append("## Plain-language conclusion")
    lines.append("")
    lines.append(_build_conclusion_sentence(total, processed, errored, n_low_res, n_rainbow, n_redgreen))
    lines.append("")

    lines.append("## Known limitation")
    lines.append("")
    lines.append(
        "This QC scan is a heuristic, pixel-based check only (see the module docstring at the top of "
        "`scripts/polish_figures.py`). It cannot detect overplotting, truncated axes, undefined error bars, "
        "or unit/precision mismatches, and it can both miss real problems and raise false positives. A flag "
        "above means a figure is worth reviewing against `references/pitfalls.md`, not a certified pass/fail. "
        "Per `references/deliverables.md`, any flagged figure should be replotted from its source data/table "
        "using `references/cns-method-guide.md`'s Design Decision Tree -- this script cannot fix colors or "
        "resolution baked into an existing raster/vector file, only re-frame it."
    )
    lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Apply lightweight publication framing to PNG/PDF figures and export PNG/PDF.")
    parser.add_argument("--manifest", help="TSV manifest from style_manifest.py")
    parser.add_argument("--out", help="Directory for polished files when using --manifest")
    parser.add_argument("--input", help="Single PNG/PDF input")
    parser.add_argument("--output-stem", help="Single output path without extension")
    parser.add_argument("--panel-label", default="")
    parser.add_argument(
        "--check-palette",
        action="store_true",
        help="Validate assets/palette.json for colorblind-safe category/diverging colors and exit",
    )
    args = parser.parse_args()

    if args.check_palette:
        problems = check_palette_colorblind_safe()
        if problems:
            print(f"palette check: {len(problems)} issue(s) found in {PALETTE_PATH}")
            for problem in problems:
                print(f"  - {problem}")
            raise SystemExit(1)
        print(f"palette check: PASS ({PALETTE_PATH})")
        return

    if args.manifest:
        manifest_file = Path(args.manifest)
        if not manifest_file.exists():
            raise SystemExit(f"error: manifest not found: {manifest_file}")
        rows = list(csv.DictReader(manifest_file.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
        if not rows:
            raise SystemExit(f"error: manifest {manifest_file} has no data rows (run scripts/style_manifest.py first)")
        for warning in validate_panel_labels(rows):
            print(f"warning: {warning}")
        qc_rows = []
        n_errors = 0
        for row in rows:
            source_file = row.get("source_file", "")
            output_stem = Path(row.get("output_stem", source_file))
            if args.out:
                output_stem = Path(args.out) / output_stem.name
            try:
                _, qc = polish_file(source_file, output_stem, row.get("panel_label", ""))
            except Exception as exc:
                n_errors += 1
                print(f"error: {source_file}: {exc}")
                qc = _error_qc(str(exc))
            else:
                for warning in qc["warnings"]:
                    print(f"warning: {source_file}: {warning}")
            qc_rows.append(
                {
                    "source_file": source_file,
                    "output_stem": str(output_stem),
                    "panel_label": row.get("panel_label", ""),
                    "qc": qc,
                }
            )
        out_dir = Path(args.out) if args.out else Path("output")
        report_path = write_qc_report(qc_rows, out_dir / "figure_qc_report.tsv")
        summary_path = write_batch_report(qc_rows, out_dir / "report.md", manifest_path=args.manifest)
        print(f"polished {len(rows) - n_errors}/{len(rows)} figure files from {args.manifest}")
        if n_errors:
            print(f"{n_errors} figure(s) failed to process -- see {summary_path} for details")
        print(f"qc report written to {report_path}")
        print(f"batch report written to {summary_path}")
        if n_errors:
            raise SystemExit(1)
    elif args.input and args.output_stem:
        try:
            formats, qc = polish_file(args.input, args.output_stem, args.panel_label)
        except Exception as exc:
            raise SystemExit(f"error: {args.input}: {exc}") from exc
        for warning in qc["warnings"]:
            print(f"warning: {args.input}: {warning}")
        print(f"polished {args.input} -> {args.output_stem}.{{{','.join(formats)}}}")
    else:
        raise SystemExit("Use --manifest or --input/--output-stem, or --check-palette")


if __name__ == "__main__":
    main()
