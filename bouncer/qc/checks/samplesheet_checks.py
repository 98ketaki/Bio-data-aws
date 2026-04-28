"""Section 4: Samplesheet integrity checks."""

from __future__ import annotations

import re

import pandas as pd

from ..finding import Finding, Severity
from ...parsers import SchemaContract


def run(samplesheet: pd.DataFrame, schema: SchemaContract, checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = samplesheet
    check_set = {c["check"] for c in checks}

    is_ntc = df["target_type"].astype(str) == "ntc" if "target_type" in df.columns \
        else pd.Series([False] * len(df))
    non_ntc = df[~is_ntc]

    # sample_id_unique
    if "sample_id_unique" in check_set and "sample_id" in df.columns:
        dupes = df[df["sample_id"].duplicated(keep=False)]["sample_id"].dropna().unique()
        if len(dupes):
            findings.append(Finding(
                check="sample_id_unique",
                severity=Severity.HARD,
                section="samplesheet",
                message=f"Duplicate sample_id(s): {', '.join(str(x) for x in dupes[:10])}.",
                samples=list(dupes[:10]),
            ))

    # well_position_unique_per_plate
    if "well_position_unique_per_plate" in check_set and \
            "well_position" in df.columns and "plate_id" in df.columns:
        dupes = df[df.duplicated(["plate_id", "well_position"], keep=False)]
        if not dupes.empty:
            findings.append(Finding(
                check="well_position_unique_per_plate",
                severity=Severity.HARD,
                section="samplesheet",
                message=(f"{len(dupes)} rows share a plate_id+well_position combination "
                         f"(duplicate wells)."),
                wells=dupes["well_position"].tolist(),
            ))

    # no_empty_required_columns
    if "no_empty_required_columns" in check_set:
        for col_spec in schema.required_columns():
            col = col_spec.name
            if col not in df.columns:
                findings.append(Finding(
                    check="no_empty_required_columns",
                    severity=Severity.HARD,
                    section="samplesheet",
                    message=f"Required column '{col}' is missing from the samplesheet.",
                ))
                continue
            check_df = non_ntc if col_spec.ntc_exempt else df
            bad = check_df[check_df[col].isna() |
                           (check_df[col].astype(str).str.strip() == "") |
                           (check_df[col].astype(str).str.lower() == "nan")]
            if not bad.empty:
                findings.append(Finding(
                    check="no_empty_required_columns",
                    severity=Severity.HARD,
                    section="samplesheet",
                    message=(f"Required column '{col}' has {len(bad)} empty/null value(s). "
                             f"Affected sample_ids: "
                             f"{', '.join(bad['sample_id'].astype(str).tolist()[:5])}"
                             f"{'...' if len(bad) > 5 else ''}."),
                    samples=bad["sample_id"].astype(str).tolist()[:10],
                ))

    # target_type_values
    if "target_type_values" in check_set and "target_type" in df.columns:
        allowed = {"reference", "target", "ntc"}
        bad = df[~df["target_type"].astype(str).isin(allowed)]
        if not bad.empty:
            bad_vals = bad["target_type"].astype(str).unique()
            findings.append(Finding(
                check="target_type_values",
                severity=Severity.HARD,
                section="samplesheet",
                message=f"Invalid target_type value(s): {', '.join(bad_vals)}.",
                samples=bad["sample_id"].astype(str).tolist()[:10],
            ))

    # task_consistent_with_target_type
    if "task_consistent_with_target_type" in check_set and \
            "task" in df.columns and "target_type" in df.columns:
        ntc_wrong_task = df[
            (df["target_type"].astype(str) == "ntc") &
            (df["task"].astype(str).str.upper() != "NTC")
        ]
        non_ntc_wrong_task = df[
            (df["target_type"].astype(str) != "ntc") &
            (df["task"].astype(str).str.upper() == "NTC")
        ]
        for bad_df, msg in [
            (ntc_wrong_task, "NTC target_type rows with task ≠ NTC"),
            (non_ntc_wrong_task, "Non-NTC rows with task = NTC"),
        ]:
            if not bad_df.empty:
                findings.append(Finding(
                    check="task_consistent_with_target_type",
                    severity=Severity.SOFT,
                    section="samplesheet",
                    message=f"{len(bad_df)} {msg}.",
                    samples=bad_df["sample_id"].astype(str).tolist()[:10],
                ))

    # cdna_input_consistent_per_sample
    if "cdna_input_consistent_per_sample" in check_set and \
            "biological_sample_id" in df.columns or "biological_replicate" in df.columns:
        bio_col = "biological_sample_id" if "biological_sample_id" in df.columns \
            else None
        for c in ["cdna_input_ng", "cdna_concentration_ng_uL"]:
            if c not in non_ntc.columns or bio_col is None:
                continue
            inconsistent = non_ntc.groupby(bio_col)[c].nunique()
            bad_ids = inconsistent[inconsistent > 1].index.tolist()
            if bad_ids:
                findings.append(Finding(
                    check="cdna_input_consistent_per_sample",
                    severity=Severity.SOFT,
                    section="samplesheet",
                    message=(f"Column '{c}' inconsistent within biological_sample_id for: "
                             f"{', '.join(str(x) for x in bad_ids[:5])}."),
                ))

    # replicate_is_integer_non_ntc
    if "replicate_is_integer_non_ntc" in check_set:
        for rep_col in ["biological_replicate", "technical_replicate"]:
            if rep_col not in non_ntc.columns:
                continue
            bad = non_ntc[non_ntc[rep_col].isna() | (non_ntc[rep_col] <= 0)]
            if not bad.empty:
                findings.append(Finding(
                    check="replicate_is_integer_non_ntc",
                    severity=Severity.HARD,
                    section="samplesheet",
                    message=(f"'{rep_col}' has {len(bad)} non-positive/null value(s) "
                             f"in non-NTC rows."),
                    samples=bad["sample_id"].astype(str).tolist()[:10],
                ))

    # organism_in_allowed_values
    if "organism_in_allowed_values" in check_set and "organism" in df.columns:
        allowed = [c.allowed_values for c in schema.columns if c.name == "organism"]
        allowed = allowed[0] if allowed else []
        if allowed:
            bad = non_ntc[~non_ntc["organism"].astype(str).isin(allowed)]
            if not bad.empty:
                bad_vals = bad["organism"].astype(str).unique()
                findings.append(Finding(
                    check="organism_in_allowed_values",
                    severity=Severity.HARD,
                    section="samplesheet",
                    message=(f"Unknown organism value(s): {', '.join(bad_vals)}. "
                             f"Allowed: {', '.join(allowed)}."),
                    samples=bad["sample_id"].astype(str).tolist()[:10],
                ))

    # passage_number_range
    if "passage_number_range" in check_set and "passage_number" in df.columns:
        pn = df["passage_number"].dropna()
        bad = df[df["passage_number"].notna() &
                 ((df["passage_number"] <= 0) | (df["passage_number"] > 100))]
        if not bad.empty:
            findings.append(Finding(
                check="passage_number_range",
                severity=Severity.SOFT,
                section="samplesheet",
                message=f"passage_number out of [1, 100] for {len(bad)} row(s).",
                samples=bad["sample_id"].astype(str).tolist()[:10],
            ))

    return findings
