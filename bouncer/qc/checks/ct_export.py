"""Section 1: Ct export integrity checks."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd

from ..finding import Finding, Severity

if TYPE_CHECKING:
    from ...parsers import CtExport


def _non_ntc(df: pd.DataFrame) -> pd.DataFrame:
    if "Task" in df.columns:
        return df[df["Task"].str.upper() != "NTC"]
    return df


def run(ct: "CtExport", checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = ct.data
    non_ntc = _non_ntc(df)

    check_set = {c["check"] for c in checks}

    # ct_values_present_non_ntc
    if "ct_values_present_non_ntc" in check_set:
        missing = non_ntc[non_ntc["CT"].isna()]
        if not missing.empty:
            wells = missing["Well Position"].tolist()
            findings.append(Finding(
                check="ct_values_present_non_ntc",
                severity=Severity.HARD,
                section="ct_export",
                message=f"{len(wells)} non-NTC well(s) have missing/Undetermined CT values.",
                wells=wells,
            ))

    # ct_numeric_type — already coerced; non-numeric would be NaN
    # (covered by ct_values_present_non_ntc)

    # ct_range_non_ntc
    range_cfg = next((c for c in checks if c["check"] == "ct_range_non_ntc"), None)
    if range_cfg and not non_ntc.empty:
        ct_col = non_ntc["CT"].dropna()
        hard_min = range_cfg.get("hard_min", 5.0)
        hard_max = range_cfg.get("hard_max", 40.0)
        soft_min = range_cfg.get("soft_min", 10.0)
        soft_max = range_cfg.get("soft_max", 38.0)

        hard_fail = non_ntc[(non_ntc["CT"] < hard_min) | (non_ntc["CT"] > hard_max)]
        if not hard_fail.empty:
            findings.append(Finding(
                check="ct_range_non_ntc",
                severity=Severity.HARD,
                section="ct_export",
                message=(f"{len(hard_fail)} wells have CT outside hard bounds "
                         f"[{hard_min}, {hard_max}]."),
                wells=hard_fail["Well Position"].tolist(),
            ))

        soft_fail = non_ntc[
            ((non_ntc["CT"] < soft_min) | (non_ntc["CT"] > soft_max)) &
            ~((non_ntc["CT"] < hard_min) | (non_ntc["CT"] > hard_max))
        ]
        if not soft_fail.empty:
            findings.append(Finding(
                check="ct_range_non_ntc",
                severity=Severity.SOFT,
                section="ct_export",
                message=(f"{len(soft_fail)} wells have CT outside soft bounds "
                         f"[{soft_min}, {soft_max}] — investigate."),
                wells=soft_fail["Well Position"].tolist(),
            ))

    # ct_technical_replicate_sd
    sd_cfg = next((c for c in checks if c["check"] == "ct_technical_replicate_sd"), None)
    if sd_cfg and not non_ntc.empty:
        hard_max = sd_cfg.get("hard_max", 1.0)
        soft_max = sd_cfg.get("soft_max", 0.5)
        grp = non_ntc.groupby(["Sample Name", "Target Name"])["CT"].std()
        hard_fail = grp[grp > hard_max]
        soft_fail = grp[(grp > soft_max) & (grp <= hard_max)]
        if not hard_fail.empty:
            pairs = [f"{s}/{t}" for s, t in hard_fail.index]
            findings.append(Finding(
                check="ct_technical_replicate_sd",
                severity=Severity.HARD,
                section="ct_export",
                message=(f"CT SD > {hard_max} for {len(hard_fail)} sample×target pair(s): "
                         f"{', '.join(pairs[:5])}{'...' if len(pairs) > 5 else ''}"),
            ))
        if not soft_fail.empty:
            pairs = [f"{s}/{t}" for s, t in soft_fail.index]
            findings.append(Finding(
                check="ct_technical_replicate_sd",
                severity=Severity.SOFT,
                section="ct_export",
                message=(f"CT SD {soft_max}–{hard_max} for {len(soft_fail)} "
                         f"sample×target pair(s): {', '.join(pairs[:5])}"
                         f"{'...' if len(pairs) > 5 else ''}"),
            ))

    # ct_outlier_replicate_detection
    if "ct_outlier_replicate_detection" in check_set and not non_ntc.empty:
        flagged = []
        for (sn, tn), grp in non_ntc.groupby(["Sample Name", "Target Name"]):
            ct_vals = grp["CT"].dropna()
            if len(ct_vals) < 2:
                continue
            mean = ct_vals.mean()
            outliers = grp[abs(grp["CT"] - mean) > 1.5]
            flagged.extend(outliers["Well Position"].tolist())
        if flagged:
            findings.append(Finding(
                check="ct_outlier_replicate_detection",
                severity=Severity.SOFT,
                section="ct_export",
                message=f"{len(flagged)} well(s) deviate >1.5 Ct from their replicate mean.",
                wells=flagged,
            ))

    # plate_barcode_in_header
    if "plate_barcode_in_header" in check_set:
        if not ct.plate_barcode:
            findings.append(Finding(
                check="plate_barcode_in_header",
                severity=Severity.SOFT,
                section="ct_export",
                message="'# Plate Barcode' not found or empty in Ct export header.",
            ))

    # run_date_in_header
    if "run_date_in_header" in check_set:
        date_ok = bool(ct.run_date and re.match(r"^\d{4}-\d{2}-\d{2}$", ct.run_date))
        if not date_ok:
            findings.append(Finding(
                check="run_date_in_header",
                severity=Severity.SOFT,
                section="ct_export",
                message=(f"Run Date missing or not YYYY-MM-DD format "
                         f"(got: '{ct.run_date or 'empty'}')."),
            ))

    return findings
