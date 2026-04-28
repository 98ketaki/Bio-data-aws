"""Section 3: Reference gene checks."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..finding import Finding, Severity

if TYPE_CHECKING:
    from ...parsers import CtExport


def run(ct: "CtExport", samplesheet: pd.DataFrame, checks: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    df = ct.data
    check_set = {c["check"] for c in checks}

    ref_wells = samplesheet[samplesheet["target_type"].astype(str) == "reference"] \
        if "target_type" in samplesheet.columns else pd.DataFrame()

    ref_targets = set(ref_wells["target_name"].dropna().unique()) \
        if not ref_wells.empty and "target_name" in ref_wells.columns else set()

    # min_reference_genes
    min_cfg = next((c for c in checks if c["check"] == "min_reference_genes"), None)
    if min_cfg:
        hard_min = min_cfg.get("hard_min", 1)
        soft_min = min_cfg.get("soft_min", 2)
        n = len(ref_targets)
        if n < hard_min:
            findings.append(Finding(
                check="min_reference_genes",
                severity=Severity.HARD,
                section="reference_genes",
                message=f"No reference genes found. At least {hard_min} required.",
            ))
        elif n < soft_min:
            findings.append(Finding(
                check="min_reference_genes",
                severity=Severity.SOFT,
                section="reference_genes",
                message=(f"Only {n} reference gene(s) found; "
                         f"{soft_min} recommended for stability cross-validation."),
            ))

    if ref_targets and "Target Name" in df.columns:
        ref_ct_rows = df[df["Target Name"].isin(ref_targets)]

        # reference_gene_ct_range
        rng_cfg = next((c for c in checks if c["check"] == "reference_gene_ct_range"), None)
        if rng_cfg and not ref_ct_rows.empty:
            hard_min = rng_cfg.get("hard_min", 10.0)
            hard_max = rng_cfg.get("hard_max", 30.0)
            soft_min = rng_cfg.get("soft_min", 14.0)
            soft_max = rng_cfg.get("soft_max", 25.0)

            hard_fail = ref_ct_rows[
                (ref_ct_rows["CT"] < hard_min) | (ref_ct_rows["CT"] > hard_max)
            ]
            if not hard_fail.empty:
                findings.append(Finding(
                    check="reference_gene_ct_range",
                    severity=Severity.HARD,
                    section="reference_genes",
                    message=(f"{len(hard_fail)} reference gene well(s) outside hard CT range "
                             f"[{hard_min}, {hard_max}]."),
                    wells=hard_fail["Well Position"].tolist(),
                ))

            soft_fail = ref_ct_rows[
                ((ref_ct_rows["CT"] < soft_min) | (ref_ct_rows["CT"] > soft_max)) &
                ~((ref_ct_rows["CT"] < hard_min) | (ref_ct_rows["CT"] > hard_max))
            ]
            if not soft_fail.empty:
                findings.append(Finding(
                    check="reference_gene_ct_range",
                    severity=Severity.SOFT,
                    section="reference_genes",
                    message=(f"{len(soft_fail)} reference gene well(s) outside soft CT range "
                             f"[{soft_min}, {soft_max}]."),
                    wells=soft_fail["Well Position"].tolist(),
                ))

        # reference_gene_stability_cv
        cv_cfg = next((c for c in checks if c["check"] == "reference_gene_stability_cv"), None)
        if cv_cfg and not ref_ct_rows.empty:
            hard_max = cv_cfg.get("hard_max", 0.25)
            soft_max = cv_cfg.get("soft_max", 0.15)
            for tgt, grp in ref_ct_rows.groupby("Target Name"):
                vals = grp["CT"].dropna()
                if len(vals) < 2:
                    continue
                cv = vals.std() / vals.mean()
                if cv > hard_max:
                    findings.append(Finding(
                        check="reference_gene_stability_cv",
                        severity=Severity.HARD,
                        section="reference_genes",
                        message=(f"Reference gene '{tgt}' CV={cv:.2%} > hard limit {hard_max:.0%}. "
                                 f"Not a stable normaliser."),
                    ))
                elif cv > soft_max:
                    findings.append(Finding(
                        check="reference_gene_stability_cv",
                        severity=Severity.SOFT,
                        section="reference_genes",
                        message=(f"Reference gene '{tgt}' CV={cv:.2%} > soft limit {soft_max:.0%} "
                                 f"— investigate stability."),
                    ))

    # reference_gene_not_in_targets
    if "reference_gene_not_in_targets" in check_set and "target_type" in samplesheet.columns:
        tgt_wells = samplesheet[samplesheet["target_type"].astype(str) == "target"]
        tgt_names = set(tgt_wells["target_name"].dropna().unique()) \
            if "target_name" in tgt_wells.columns else set()
        dual = ref_targets & tgt_names
        if dual:
            findings.append(Finding(
                check="reference_gene_not_in_targets",
                severity=Severity.HARD,
                section="reference_genes",
                message=(f"Gene(s) listed as both reference and target: "
                         f"{', '.join(sorted(dual))}."),
            ))

    return findings
