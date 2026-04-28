"""Parses QuantStudio Ct export CSVs (header block + data table)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class CtExport:
    header: dict[str, str]
    data: pd.DataFrame

    # convenience accessors
    @property
    def plate_barcode(self) -> str:
        return self.header.get("Plate Barcode", "")

    @property
    def run_date(self) -> str:
        return self.header.get("Run Date", "")

    @property
    def operator(self) -> str:
        return self.header.get("Operator", "")

    @property
    def experiment_name(self) -> str:
        return self.header.get("Experiment Name", "")


def parse_ct_export(path: str | Path) -> CtExport:
    """Parse a QuantStudio Ct export CSV.

    The file starts with comment lines of the form  ``# Key,Value,...``
    followed by a blank line, then a standard CSV table.
    """
    path = Path(path)
    raw_lines = path.read_text(encoding="utf-8-sig").splitlines()

    header: dict[str, str] = {}
    data_start = 0

    for i, line in enumerate(raw_lines):
        stripped = line.strip()
        if stripped.startswith("#"):
            # '# Key,Value,...' — take first two comma-separated parts
            parts = stripped[1:].split(",", 1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().rstrip(",")
                if key:
                    header[key] = value
        elif header and stripped.replace(",", "") == "":
            # separator line (empty or all-comma) after header block
            for j in range(i + 1, len(raw_lines)):
                candidate = raw_lines[j].strip()
                if candidate and not candidate.startswith("#") and candidate.replace(",", "") != "":
                    data_start = j
                    break
            break

    # read data table
    from io import StringIO
    table_text = "\n".join(raw_lines[data_start:])
    data = pd.read_csv(StringIO(table_text))

    # normalise column names: strip whitespace
    data.columns = [c.strip() for c in data.columns]

    # coerce CT to numeric (Undetermined → NaN)
    if "CT" in data.columns:
        data["CT"] = pd.to_numeric(data["CT"], errors="coerce")

    return CtExport(header=header, data=data)
