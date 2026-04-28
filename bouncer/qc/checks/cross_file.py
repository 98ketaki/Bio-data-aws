"""Section 5: Cross-file concordance checks."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..finding import Finding, Severity

if TYPE_CHECKING:
    from ...parsers import CtExport


def run(ct: "CtExport", samplesheet: pd.DataFrame, checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = ct.data
    ss = samplesheet
    check_set = {c["check"] for c in checks}

    ct_wells = set(df["Well Position"].dropna().astype(str)) if "Well Position" in df.columns else set()
    ss_wells = set(ss["well_position"].dropna().astype(str)) if "well_position" in ss.columns else set()

    # well_positions_match
    if "well_positions_match" in check_set:
        only_ct = ct_wells - ss_wells
        only_ss = ss_wells - ct_wells
        if only_ct:
            findings.append(Finding(
                check="well_positions_match",
                severity=Severity.HARD,
                section="cross_file",
                message=(f"{len(only_ct)} well(s) in Ct export not found in samplesheet: "
                         f"{', '.join(sorted(only_ct)[:10])}."),
                wells=sorted(only_ct),
            ))
        if only_ss:
            findings.append(Finding(
                check="well_positions_match",
                severity=Severity.HARD,
                section="cross_file",
                message=(f"{len(only_ss)} well(s) in samplesheet not found in Ct export: "
                         f"{', '.join(sorted(only_ss)[:10])}."),
                wells=sorted(only_ss),
            ))

    # sample_names_match / target_names_match
    if ({"sample_names_match", "target_names_match"} & check_set) and \
            "well_position" in ss.columns and "Well Position" in df.columns:
        merged = ss.merge(
            df[["Well Position", "Sample Name", "Target Name"]].drop_duplicates("Well Position"),
            left_on="well_position", right_on="Well Position", how="inner",
        )

        if "sample_names_match" in check_set and \
                "biological_sample_id" in merged.columns and "Sample Name" in merged.columns:
            # try matching on biological_sample_id
            mismatch = merged[merged["biological_sample_id"].astype(str) !=
                              merged["Sample Name"].astype(str)]
            if not mismatch.empty:
                findings.append(Finding(
                    check="sample_names_match",
                    severity=Severity.HARD,
                    section="cross_file",
                    message=(f"{len(mismatch)} well(s) have biological_sample_id ≠ "
                             f"Ct export Sample Name."),
                    wells=mismatch["well_position"].tolist()[:10],
                ))

        if "target_names_match" in check_set and \
                "target_name" in merged.columns and "Target Name" in merged.columns:
            mismatch = merged[merged["target_name"].astype(str) !=
                              merged["Target Name"].astype(str)]
            if not mismatch.empty:
                findings.append(Finding(
                    check="target_names_match",
                    severity=Severity.HARD,
                    section="cross_file",
                    message=(f"{len(mismatch)} well(s) have target_name ≠ "
                             f"Ct export Target Name."),
                    wells=mismatch["well_position"].tolist()[:10],
                ))

    # plate_barcode_matches_samplesheet
    if "plate_barcode_matches_samplesheet" in check_set and "plate_id" in ss.columns:
        plate_ids = set(ss["plate_id"].dropna().astype(str))
        if ct.plate_barcode and ct.plate_barcode not in plate_ids:
            findings.append(Finding(
                check="plate_barcode_matches_samplesheet",
                severity=Severity.SOFT,
                section="cross_file",
                message=(f"Ct export Plate Barcode '{ct.plate_barcode}' not found among "
                         f"samplesheet plate_id(s): {', '.join(sorted(plate_ids))}."),
            ))

    # ct_mean_consistent_with_raw_ct
    if "ct_mean_consistent_with_raw_ct" in check_set and \
            "Ct Mean" in df.columns and "CT" in df.columns:
        tol = 0.05
        non_ntc = df[df["Task"].str.upper() != "NTC"] if "Task" in df.columns else df
        grp = non_ntc.groupby(["Sample Name", "Target Name"])
        for (sn, tn), g in grp:
            raw_mean = g["CT"].mean()
            inst_mean = g["Ct Mean"].dropna()
            if inst_mean.empty:
                continue
            inst_val = inst_mean.iloc[0]
            if abs(raw_mean - inst_val) > tol:
                findings.append(Finding(
                    check="ct_mean_consistent_with_raw_ct",
                    severity=Severity.SOFT,
                    section="cross_file",
                    message=(f"Sample '{sn}' × target '{tn}': "
                             f"computed mean CT={raw_mean:.3f} vs instrument Ct Mean={inst_val:.3f} "
                             f"(diff={abs(raw_mean-inst_val):.3f} > {tol})."),
                ))

    return findings
