"""Extracts plain text from a protocol PDF for context-aware QC."""

from __future__ import annotations

from pathlib import Path


def parse_protocol(path: str | Path | None) -> str:
    """Return protocol text. Returns empty string if path is None or extraction fails."""
    if path is None:
        return ""
    path = Path(path)
    if not path.exists():
        return ""
    try:
        import pypdf
        reader = pypdf.PdfReader(str(path))
        return "\n".join(page.extract_text() or "" for page in reader.pages)
    except Exception:
        return ""
