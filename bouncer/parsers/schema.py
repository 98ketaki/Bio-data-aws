"""Loads and exposes the YAML schema contract and QC contract."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class ColumnSpec:
    name: str
    dtype: str
    required: bool
    unique: bool = False
    ntc_exempt: bool = False
    allowed_values: list[str] = field(default_factory=list)
    pattern: str | None = None
    description: str = ""


@dataclass
class SchemaContract:
    version: str
    assay_type: str
    instrument_format: str
    data_stage: str
    index_column: str
    columns: list[ColumnSpec]
    output_features: list[Any]

    def required_columns(self) -> list[ColumnSpec]:
        return [c for c in self.columns if c.required]

    def column(self, name: str) -> ColumnSpec | None:
        return next((c for c in self.columns if c.name == name), None)


@dataclass
class QCContract:
    version: str
    assay_type: str
    ct_export_checks: list[dict]
    ntc_checks: list[dict]
    reference_gene_checks: list[dict]
    samplesheet_checks: list[dict]
    cross_file_checks: list[dict]
    design_checks: list[dict]


def load_schema(path: str | Path) -> SchemaContract:
    raw = yaml.safe_load(Path(path).read_text())
    columns = [
        ColumnSpec(
            name=col["name"],
            dtype=col.get("dtype", "str"),
            required=col.get("required", False),
            unique=col.get("unique", False),
            ntc_exempt=col.get("ntc_exempt", False),
            allowed_values=col.get("allowed_values", []),
            pattern=col.get("pattern"),
            description=col.get("description", ""),
        )
        for col in raw.get("metadata_columns", [])
    ]
    return SchemaContract(
        version=raw.get("version", ""),
        assay_type=raw.get("assay_type", ""),
        instrument_format=raw.get("instrument_format", ""),
        data_stage=raw.get("data_stage", ""),
        index_column=raw.get("index_column", "sample_id"),
        columns=columns,
        output_features=raw.get("output_features", []),
    )


def load_qc_contract(path: str | Path) -> QCContract:
    raw = yaml.safe_load(Path(path).read_text())
    return QCContract(
        version=raw.get("version", ""),
        assay_type=raw.get("assay_type", ""),
        ct_export_checks=raw.get("ct_export", []),
        ntc_checks=raw.get("ntc", []),
        reference_gene_checks=raw.get("reference_genes", []),
        samplesheet_checks=raw.get("samplesheet", []),
        cross_file_checks=raw.get("cross_file", []),
        design_checks=raw.get("design", []),
    )
