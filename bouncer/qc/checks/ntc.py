"""Section 2: NTC contamination checks."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..finding import Finding, Severity

if TYPE_CHECKING:
    from ...parsers import CtExport


def run(ct: "CtExport", checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = ct.data
    check_set = {c["check"] for c in checks}

    ntc_rows = df[df["Task"].str.upper() == "NTC"] if "Task" in df.columns else pd.DataFrame()

    # ntc_ct_undetermined
    if "ntc_ct_undetermined" in check_set and not ntc_rows.empty:
        amplified = ntc_rows[ntc_rows["CT"].notna()]
        if not amplified.empty:
            findings.append(Finding(
                check="ntc_ct_undetermined",
                severity=Severity.HARD,
                section="ntc",
                message=(f"{len(amplified)} NTC well(s) show amplification "
                         f"(numeric CT) — contamination detected."),
                wells=amplified["Well Position"].tolist(),
            ))

    # ntc_ct_high_threshold
    nt_cfg = next((c for c in checks if c["check"] == "ntc_ct_high_threshold"), None)
    if nt_cfg and not ntc_rows.empty:
        threshold = nt_cfg.get("hard_max_ct_for_ntc", 34.9)
        high_contam = ntc_rows[ntc_rows["CT"].notna() & (ntc_rows["CT"] < threshold)]
        if not high_contam.empty:
            findings.append(Finding(
                check="ntc_ct_high_threshold",
                severity=Severity.HARD,
                section="ntc",
                message=(f"{len(high_contam)} NTC well(s) have CT < {threshold} — "
                         f"high-level contamination."),
                wells=high_contam["Well Position"].tolist(),
            ))

    # ntc_present_per_target
    if "ntc_present_per_target" in check_set and "Target Name" in df.columns:
        all_targets = set(df["Target Name"].dropna().unique())
        ntc_targets = set(ntc_rows["Target Name"].dropna().unique()) if not ntc_rows.empty else set()
        missing_ntc = all_targets - ntc_targets
        if missing_ntc:
            findings.append(Finding(
                check="ntc_present_per_target",
                severity=Severity.WARNING,
                section="ntc",
                message=(f"No NTC well found for {len(missing_ntc)} target(s): "
                         f"{', '.join(sorted(missing_ntc))}."),
            ))

    return findings
