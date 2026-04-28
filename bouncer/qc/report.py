"""Generates human-readable QC reports from a QCResult."""

from __future__ import annotations

from .engine import QCResult
from .finding import Severity

_ICONS = {Severity.HARD: "✖ HARD", Severity.SOFT: "⚠ SOFT", Severity.WARNING: "ℹ WARN"}
_SECTION_LABELS = {
    "ct_export": "CT Export Integrity",
    "ntc": "NTC Contamination",
    "reference_genes": "Reference Genes",
    "samplesheet": "Samplesheet Integrity",
    "cross_file": "Cross-File Concordance",
    "design": "Experimental Design",
}


def generate_report(result: QCResult, verbose: bool = False) -> str:
    lines: list[str] = []

    lines.append("=" * 70)
    lines.append("BOUNCER QC REPORT")
    lines.append("=" * 70)

    assay = result.assay_type.upper()
    schema_ver = result.schema.version if result.schema else "?"
    lines.append(f"Assay type : {assay}")
    lines.append(f"Schema     : v{schema_ver}")
    if result.ct:
        lines.append(f"Experiment : {result.ct.experiment_name}")
        lines.append(f"Plate      : {result.ct.plate_barcode}  Run date: {result.ct.run_date}")
        lines.append(f"Operator   : {result.ct.operator}")
    lines.append("")

    n_hard = len(result.hard_failures)
    n_soft = len(result.soft_failures)
    n_warn = len(result.warnings)

    status = "PASS ✓" if result.passed else "FAIL ✗"
    lines.append(f"Overall status : {status}")
    lines.append(f"  Hard failures : {n_hard}  (block registration)")
    lines.append(f"  Soft failures : {n_soft}  (logged as warnings)")
    lines.append(f"  Advisories    : {n_warn}  (design / metadata)")
    lines.append("")

    # group by section
    sections: dict[str, list] = {}
    for f in result.findings:
        sections.setdefault(f.section, []).append(f)

    for section_key, label in _SECTION_LABELS.items():
        section_findings = sections.get(section_key, [])
        if not section_findings:
            lines.append(f"[{label}]  — all checks passed")
            continue

        lines.append(f"[{label}]")
        for f in section_findings:
            icon = _ICONS.get(f.severity, "?")
            lines.append(f"  {icon}  {f.check}")
            lines.append(f"         {f.message}")
            if verbose and f.wells:
                lines.append(f"         Wells   : {', '.join(str(w) for w in f.wells[:20])}"
                             f"{'...' if len(f.wells) > 20 else ''}")
            if verbose and f.samples:
                lines.append(f"         Samples : {', '.join(str(s) for s in f.samples[:20])}"
                             f"{'...' if len(f.samples) > 20 else ''}")
        lines.append("")

    if result.passed:
        lines.append("Registration is allowed. Run bouncer.register() to push to feature store.")
    else:
        lines.append("Registration BLOCKED. Fix hard failures above before re-submitting.")
        lines.append("")
        lines.append("Fix guide:")
        for i, f in enumerate(result.hard_failures, 1):
            lines.append(f"  {i}. [{f.section}] {f.check}: {f.message}")

    lines.append("=" * 70)
    return "\n".join(lines)


def print_report(result: QCResult, verbose: bool = False) -> None:
    print(generate_report(result, verbose=verbose))
