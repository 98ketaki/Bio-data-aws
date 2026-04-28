"""Parses the plate samplesheet CSV and applies schema dtype coercion."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .schema import SchemaContract


def parse_samplesheet(path: str | Path, schema: SchemaContract) -> pd.DataFrame:
    df = pd.read_csv(Path(path))
    df.columns = [c.strip() for c in df.columns]

    for col_spec in schema.columns:
        if col_spec.name not in df.columns:
            continue
        col = col_spec.name
        dtype = col_spec.dtype

        if dtype in ("str",):
            df[col] = df[col].astype(str).replace("nan", pd.NA)
        elif dtype == "int":
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        elif dtype == "float":
            df[col] = pd.to_numeric(df[col], errors="coerce")
        elif dtype in ("category",):
            df[col] = df[col].astype(str).replace("nan", pd.NA).astype("category")
        elif dtype == "bool":
            df[col] = df[col].map({"True": True, "False": False, True: True, False: False})

    return df
