"""Section 6: Experimental design advisory checks (never block)."""

from __future__ import annotations

import pandas as pd

from ..finding import Finding, Severity


def run(samplesheet: pd.DataFrame, checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = samplesheet
    check_set = {c["check"] for c in checks}

    is_ntc = df["target_type"].astype(str) == "ntc" if "target_type" in df.columns \
        else pd.Series([False] * len(df))
    non_ntc = df[~is_ntc]

    # min_biological_replicates_per_condition
    bio_rep_cfg = next((c for c in checks if c["check"] == "min_biological_replicates_per_condition"), None)
    if bio_rep_cfg and "condition" in non_ntc.columns and "biological_replicate" in non_ntc.columns:
        min_val = bio_rep_cfg.get("value", 2)
        recommended = bio_rep_cfg.get("recommended", 3)
        counts = non_ntc.groupby("condition", observed=True)["biological_replicate"].nunique()
        below = counts[counts < min_val]
        warn = counts[(counts >= min_val) & (counts < recommended)]
        if not below.empty:
            findings.append(Finding(
                check="min_biological_replicates_per_condition",
                severity=Severity.WARNING,
                section="design",
                message=(f"Condition(s) with < {min_val} biological replicates: "
                         f"{', '.join(f'{c}({n})' for c, n in below.items())}."),
            ))
        if not warn.empty:
            findings.append(Finding(
                check="min_biological_replicates_per_condition",
                severity=Severity.WARNING,
                section="design",
                message=(f"Condition(s) with < {recommended} (recommended) biological replicates: "
                         f"{', '.join(f'{c}({n})' for c, n in warn.items())}."),
            ))

    # min_technical_replicates_per_well
    tech_rep_cfg = next((c for c in checks if c["check"] == "min_technical_replicates_per_well"), None)
    if tech_rep_cfg and "technical_replicate" in non_ntc.columns and "target_name" in non_ntc.columns:
        min_val = tech_rep_cfg.get("value", 2)
        recommended = tech_rep_cfg.get("recommended", 3)
        grp_col = "biological_replicate" if "biological_replicate" in non_ntc.columns else None
        if grp_col:
            counts = non_ntc.groupby([grp_col, "target_name"])["technical_replicate"].nunique()
            below = counts[counts < min_val]
            if not below.empty:
                findings.append(Finding(
                    check="min_technical_replicates_per_well",
                    severity=Severity.WARNING,
                    section="design",
                    message=f"{len(below)} sample×target group(s) have < {min_val} technical replicates.",
                ))
            warn = counts[(counts >= min_val) & (counts < recommended)]
            if not warn.empty:
                findings.append(Finding(
                    check="min_technical_replicates_per_well",
                    severity=Severity.WARNING,
                    section="design",
                    message=(f"{len(warn)} sample×target group(s) have < {recommended} "
                             f"(recommended) technical replicates."),
                ))

    # control_condition_present
    ctrl_cfg = next((c for c in checks if c["check"] == "control_condition_present"), None)
    if ctrl_cfg and "condition" in non_ntc.columns:
        accepted = {x.lower() for x in ctrl_cfg.get("accepted_labels", [])}
        conditions = {str(x).lower() for x in non_ntc["condition"].dropna().unique()}
        if not (accepted & conditions):
            findings.append(Finding(
                check="control_condition_present",
                severity=Severity.WARNING,
                section="design",
                message=(f"No recognisable control condition found. "
                         f"Present conditions: {', '.join(sorted(conditions))}. "
                         f"Expected one of: {', '.join(sorted(accepted)[:10])}."),
            ))

    # batch_column_present
    if "batch_column_present" in check_set:
        if "batch" not in df.columns or df["batch"].isna().all():
            findings.append(Finding(
                check="batch_column_present",
                severity=Severity.WARNING,
                section="design",
                message="'batch' column absent or all-null — batch effects cannot be tracked.",
            ))

    # sex_column_present
    if "sex_column_present" in check_set:
        if "sex" not in df.columns or df["sex"].isna().all():
            findings.append(Finding(
                check="sex_column_present",
                severity=Severity.WARNING,
                section="design",
                message="'sex' column absent or all-null — sex covariate cannot be tracked.",
            ))

    # batch_not_confounded_with_condition
    if "batch_not_confounded_with_condition" in check_set and \
            "batch" in non_ntc.columns and "condition" in non_ntc.columns:
        pivot = non_ntc.dropna(subset=["batch", "condition"]) \
                       .groupby(["batch", "condition"], observed=True).size().unstack(fill_value=0)
        non_zero = (pivot > 0).sum(axis=1)  # conditions per batch
        if (non_zero == 1).all() and len(pivot) > 1:
            findings.append(Finding(
                check="batch_not_confounded_with_condition",
                severity=Severity.WARNING,
                section="design",
                message="Batch is perfectly confounded with condition — cannot separate effects.",
            ))

    # high_passage_number
    hp_cfg = next((c for c in checks if c["check"] == "high_passage_number"), None)
    if hp_cfg and "passage_number" in df.columns:
        threshold = hp_cfg.get("threshold", 50)
        high = df[df["passage_number"].notna() & (df["passage_number"] > threshold)]
        if not high.empty:
            findings.append(Finding(
                check="high_passage_number",
                severity=Severity.WARNING,
                section="design",
                message=(f"{len(high)} row(s) have passage_number > {threshold} — "
                         f"risk of phenotypic drift."),
                samples=high["sample_id"].astype(str).tolist()[:10],
            ))

    # all_conditions_have_replicates
    if "all_conditions_have_replicates" in check_set and \
            "condition" in non_ntc.columns and "biological_replicate" in non_ntc.columns:
        counts = non_ntc.groupby("condition", observed=True)["biological_replicate"].nunique()
        single = counts[counts < 2]
        if not single.empty:
            findings.append(Finding(
                check="all_conditions_have_replicates",
                severity=Severity.WARNING,
                section="design",
                message=(f"Single-replicate condition(s) — no power for significance testing: "
                         f"{', '.join(str(c) for c in single.index)}."),
            ))

    # ct_sd_outlier_samples
    return findings
