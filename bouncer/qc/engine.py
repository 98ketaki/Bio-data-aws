"""Orchestrates all QC check sections and returns a QCResult."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .finding import Finding, Severity
from .checks import ct_export as ct_checks
from .checks import ntc as ntc_checks
from .checks import reference_genes as ref_checks
from .checks import samplesheet_checks as ss_checks
from .checks import cross_file as cf_checks
from .checks import design as design_checks
from ..parsers import (
    CtExport, SchemaContract, QCContract,
    parse_ct_export, parse_samplesheet, load_schema, load_qc_contract, parse_protocol,
)


@dataclass
class QCResult:
    findings: list[Finding] = field(default_factory=list)
    ct: CtExport | None = None
    samplesheet: pd.DataFrame | None = None
    schema: SchemaContract | None = None
    qc_contract: QCContract | None = None
    protocol_text: str = ""

    @property
    def hard_failures(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == Severity.HARD]

    @property
    def soft_failures(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == Severity.SOFT]

    @property
    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == Severity.WARNING]

    @property
    def passed(self) -> bool:
        return len(self.hard_failures) == 0

    @property
    def assay_type(self) -> str:
        return self.schema.assay_type if self.schema else "unknown"


def run_qc(
    ct_export_path: str | Path,
    samplesheet_path: str | Path,
    schema_path: str | Path,
    qc_contract_path: str | Path,
    protocol_path: str | Path | None = None,
) -> QCResult:
    schema = load_schema(schema_path)
    qc_contract = load_qc_contract(qc_contract_path)
    ct = parse_ct_export(ct_export_path)
    samplesheet = parse_samplesheet(samplesheet_path, schema)
    protocol_text = parse_protocol(protocol_path)

    findings: list[Finding] = []

    findings += ct_checks.run(ct, qc_contract.ct_export_checks)
    findings += ntc_checks.run(ct, qc_contract.ntc_checks)
    findings += ref_checks.run(ct, samplesheet, qc_contract.reference_gene_checks)
    findings += ss_checks.run(samplesheet, schema, qc_contract.samplesheet_checks)
    findings += cf_checks.run(ct, samplesheet, qc_contract.cross_file_checks)
    findings += design_checks.run(samplesheet, qc_contract.design_checks)

    return QCResult(
        findings=findings,
        ct=ct,
        samplesheet=samplesheet,
        schema=schema,
        qc_contract=qc_contract,
        protocol_text=protocol_text,
    )
