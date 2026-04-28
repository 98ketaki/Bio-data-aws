"""Converts a validated QCResult into an AnnData object ready for the feature store."""

from __future__ import annotations

import hashlib
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ..qc.engine import QCResult


def build_anndata(result: QCResult) -> "anndata.AnnData":
    """Build AnnData from a QCResult that has passed hard QC.

    .X        — mean Ct per (sample × target), shape (n_samples, n_targets)
    .layers   — {'ct_sd': SD of technical replicates}
    .obs      — sample metadata (one row per biological_sample_id)
    .var      — assay/target metadata (one row per target_name)
    .uns      — plate metadata + provenance
    """
    import anndata as ad

    ct = result.ct
    ss = result.samplesheet

    df = ct.data.copy()
    non_ntc_mask = df["Task"].str.upper() != "NTC" if "Task" in df.columns \
        else pd.Series([True] * len(df))
    df = df[non_ntc_mask].copy()

    # merge samplesheet metadata into ct rows by well_position
    if "well_position" in ss.columns and "Well Position" in df.columns:
        df = df.merge(
            ss.drop(columns=[c for c in ss.columns if c in ("task",)], errors="ignore"),
            left_on="Well Position", right_on="well_position", how="left",
        )

    # pivot: rows = Sample Name (biological sample), cols = Target Name
    id_col = "Sample Name"
    target_col = "Target Name"

    ct_mean = df.groupby([id_col, target_col])["CT"].mean().unstack(target_col)
    ct_sd = df.groupby([id_col, target_col])["CT"].std().unstack(target_col)

    sample_ids = ct_mean.index.tolist()
    target_names = ct_mean.columns.tolist()

    X = ct_mean.values.astype(np.float32)

    # .obs — one row per biological sample, metadata from samplesheet
    obs_cols = [
        "organism", "tissue", "cell_type", "condition",
        "biological_replicate", "batch", "sex", "passage_number",
        "experiment_id", "plate_id",
        "cdna_input_ng", "cdna_concentration_ng_uL",
        "reverse_transcription_kit", "master_mix", "protocol_version",
    ]
    obs = (
        df.drop_duplicates(subset=[id_col])
        .set_index(id_col)
        [[c for c in obs_cols if c in df.columns]]
        .reindex(sample_ids)
    )

    # .var — one row per target (from merged df which has samplesheet cols)
    var_cols = ["target_type", "primer_assay_id"]
    var = (
        df[["Target Name"] + [c for c in var_cols if c in df.columns]]
        .drop_duplicates("Target Name")
        .set_index("Target Name")
        .reindex(target_names)
    ) if "Target Name" in df.columns else pd.DataFrame(index=target_names)

    # .uns — plate provenance
    uns = {
        "assay_type": result.assay_type,
        "schema_version": result.schema.version if result.schema else "",
        "registration_timestamp": datetime.utcnow().isoformat(),
        "plate_barcode": ct.plate_barcode,
        "run_date": ct.run_date,
        "operator": ct.operator,
        "experiment_name": ct.experiment_name,
        "qc_warnings": [
            {"check": f.check, "severity": f.severity.value, "message": f.message}
            for f in result.soft_failures + result.warnings
        ],
    }

    adata = ad.AnnData(X=X, obs=obs, var=var, uns=uns)
    adata.layers["ct_sd"] = ct_sd.reindex(index=sample_ids, columns=target_names).values.astype(np.float32)
    adata.obs_names = sample_ids
    adata.var_names = target_names

    return adata
