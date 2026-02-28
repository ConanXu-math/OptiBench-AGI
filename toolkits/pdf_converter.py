"""PDF-to-Markdown conversion using *marker* (preferred) or *nougat* CLI."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def convert_pdf_to_markdown(pdf_path: str) -> str:
    """Convert a local PDF file to Markdown text.

    Tries ``marker`` first; falls back to ``nougat`` if marker is unavailable.

    Args:
        pdf_path: Absolute or relative path to the PDF file.

    Returns:
        Markdown text extracted from the PDF, or an error message.
    """
    pdf = Path(pdf_path)
    if not pdf.exists():
        return f"ERROR: file not found – {pdf_path}"

    if shutil.which("marker"):
        return _convert_with_marker(pdf)
    if shutil.which("nougat"):
        return _convert_with_nougat(pdf)
    return (
        "ERROR: neither 'marker' nor 'nougat' CLI found on PATH. "
        "Install via: pip install marker-pdf  or  pip install nougat-ocr"
    )


def _convert_with_marker(pdf: Path) -> str:
    outdir = tempfile.mkdtemp(prefix="optibench_marker_")
    try:
        subprocess.run(
            ["marker_single", str(pdf), outdir],
            capture_output=True,
            text=True,
            timeout=int(os.getenv("OPTIBENCH_MARKER_TIMEOUT", "300")),
        )
        md_files = list(Path(outdir).rglob("*.md"))
        if not md_files:
            return "ERROR: marker produced no markdown output."
        return md_files[0].read_text(encoding="utf-8")
    except FileNotFoundError:
        return "ERROR: marker_single command not found."
    except subprocess.TimeoutExpired:
        return "ERROR: marker conversion timed out."
    finally:
        shutil.rmtree(outdir, ignore_errors=True)


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
        return "ERROR: nougat conversion timed out."
    finally:
        shutil.rmtree(outdir, ignore_errors=True)
