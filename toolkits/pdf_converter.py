"""PDF-to-Markdown conversion with automatic fallback chain.

Priority: marker → nougat → pypdf (always available).
If a higher-priority converter fails or times out, the next one is tried.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger("optibench.pdf")


def _find_bin(name: str) -> str | None:
    """Locate an executable by *name*, also checking the active venv's bin/."""
    found = shutil.which(name)
    if found:
        return found
    venv_bin = Path(sys.prefix) / "bin" / name
    if venv_bin.is_file() and os.access(venv_bin, os.X_OK):
        return str(venv_bin)
    return None


def convert_pdf_to_markdown(pdf_path: str) -> str:
    """Convert a local PDF file to Markdown text.

    Tries converters in order (marker → nougat → pypdf).
    If one fails or times out, automatically falls back to the next.
    """
    pdf = Path(pdf_path)
    if not pdf.exists():
        return f"ERROR: file not found – {pdf_path}"

    backend = os.getenv("OPTIBENCH_PDF_BACKEND", "pypdf").strip().lower()

    if backend == "pypdf":
        return _convert_with_pypdf(pdf)
    if backend == "marker":
        return _convert_with_marker(pdf)
    if backend == "nougat":
        return _convert_with_nougat(pdf)

    if _find_bin("marker_single") or _find_bin("marker"):
        result = _convert_with_marker(pdf)
        if not result.startswith("ERROR"):
            return result
        logger.warning("marker failed (%s), falling back", result[:80])

    if _find_bin("nougat"):
        result = _convert_with_nougat(pdf)
        if not result.startswith("ERROR"):
            return result
        logger.warning("nougat failed (%s), falling back", result[:80])

    return _convert_with_pypdf(pdf)


# ------------------------------------------------------------------
# pypdf (always available)
# ------------------------------------------------------------------

def _convert_with_pypdf(pdf: Path) -> str:
    try:
        from pypdf import PdfReader
    except ImportError:
        return (
            "ERROR: no PDF converter available. Install one of:\n"
            "  pip install marker-pdf   (best quality)\n"
            "  pip install nougat-ocr\n"
            "  pip install pypdf        (basic fallback)"
        )
    try:
        reader = PdfReader(str(pdf))
        pages: list[str] = []
        for i, page in enumerate(reader.pages):
            text = page.extract_text() or ""
            if text.strip():
                pages.append(f"<!-- page {i + 1} -->\n{text}")
        if not pages:
            return "ERROR: pypdf extracted no text from the PDF."
        header = (
            "> **Note:** This text was extracted with pypdf (plain-text fallback). "
            "LaTeX formulas may be garbled. Install `marker-pdf` for better results.\n\n"
        )
        return header + "\n\n".join(pages)
    except Exception as exc:
        return f"ERROR: pypdf failed – {exc}"


# ------------------------------------------------------------------
# marker
# ------------------------------------------------------------------

def _convert_with_marker(pdf: Path) -> str:
    outdir = tempfile.mkdtemp(prefix="optibench_marker_")
    ms = _find_bin("marker_single")
    if ms:
        cmd = [ms, str(pdf), "--output_dir", outdir, "--output_format", "markdown"]
    else:
        cmd = [_find_bin("marker") or "marker", "convert", str(pdf),
               "--output_dir", outdir, "--output_format", "markdown"]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=int(os.getenv("OPTIBENCH_MARKER_TIMEOUT", "300")),
        )
        md_files = list(Path(outdir).rglob("*.md"))
        if not md_files:
            stderr = result.stderr.strip()
            return f"ERROR: marker produced no output. {stderr[:300]}"
        return md_files[0].read_text(encoding="utf-8")
    except FileNotFoundError:
        return "ERROR: marker command not found."
    except subprocess.TimeoutExpired:
        return "ERROR: marker timed out."
    finally:
        shutil.rmtree(outdir, ignore_errors=True)


# ------------------------------------------------------------------
# nougat
# ------------------------------------------------------------------

def _convert_with_nougat(pdf: Path) -> str:
    outdir = tempfile.mkdtemp(prefix="optibench_nougat_")
    try:
        subprocess.run(
            ["nougat", str(pdf), "-o", outdir],
            capture_output=True,
            text=True,
            timeout=int(os.getenv("OPTIBENCH_NOUGAT_TIMEOUT", "300")),
        )
        mmd_files = list(Path(outdir).rglob("*.mmd"))
        md_files = list(Path(outdir).rglob("*.md"))
        candidates = mmd_files or md_files
        if not candidates:
            return "ERROR: nougat produced no output."
        return candidates[0].read_text(encoding="utf-8")
    except FileNotFoundError:
        return "ERROR: nougat command not found."
    except subprocess.TimeoutExpired:
        return "ERROR: nougat timed out."
    finally:
        shutil.rmtree(outdir, ignore_errors=True)
