"""DuckDB-backed feature store for registered experiments.

Schema (flat table for querying):
  experiments   — one row per registered AnnData (provenance metadata)
  observations  — one row per sample_id × target_name (Ct value + all metadata)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


_DEFAULT_DB = Path(__file__).parent.parent.parent / "bouncer_store.duckdb"


def _connect(db_path: str | Path | None = None):
    import duckdb
    path = str(db_path or _DEFAULT_DB)
    conn = duckdb.connect(path)
    _init_schema(conn)
    return conn


def _init_schema(conn) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS experiments (
            experiment_id   VARCHAR PRIMARY KEY,
            assay_type      VARCHAR,
            schema_version  VARCHAR,
            plate_barcode   VARCHAR,
            run_date        VARCHAR,
            operator        VARCHAR,
            registered_at   TIMESTAMP DEFAULT current_timestamp
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS observations (
            experiment_id              VARCHAR,
            sample_id                  VARCHAR,
            target_name                VARCHAR,
            target_type                VARCHAR,
            ct_mean                    FLOAT,
            ct_sd                      FLOAT,
            -- sample metadata
            assay_type                 VARCHAR,
            organism                   VARCHAR,
            tissue                     VARCHAR,
            cell_type                  VARCHAR,
            condition                  VARCHAR,
            biological_replicate       INTEGER,
            batch                      VARCHAR,
            sex                        VARCHAR,
            passage_number             INTEGER,
            plate_id                   VARCHAR,
            cdna_input_ng              FLOAT,
            cdna_concentration_ng_uL   FLOAT,
            reverse_transcription_kit  VARCHAR,
            master_mix                 VARCHAR,
            primer_assay_id            VARCHAR,
            protocol_version           VARCHAR,
            PRIMARY KEY (experiment_id, sample_id, target_name)
        )
    """)


def register(adata, db_path: str | Path | None = None) -> None:
    """Push an AnnData object to the DuckDB feature store."""
    import numpy as np
    conn = _connect(db_path)

    uns = adata.uns
    exp_id = str(uns.get("experiment_name", "unknown"))

    # upsert experiment row
    conn.execute("""
        INSERT OR REPLACE INTO experiments
          (experiment_id, assay_type, schema_version, plate_barcode, run_date, operator)
        VALUES (?, ?, ?, ?, ?, ?)
    """, [
        exp_id,
        uns.get("assay_type", ""),
        uns.get("schema_version", ""),
        uns.get("plate_barcode", ""),
        uns.get("run_date", ""),
        uns.get("operator", ""),
    ])

    # build long-form observations dataframe
    rows = []
    obs = adata.obs
    var = adata.var
    ct_sd_layer = adata.layers.get("ct_sd")

    for i, sample_id in enumerate(adata.obs_names):
        obs_row = obs.iloc[i]
        for j, target_name in enumerate(adata.var_names):
            ct_mean_val = float(adata.X[i, j]) if not np.isnan(adata.X[i, j]) else None
            ct_sd_val = float(ct_sd_layer[i, j]) \
                if ct_sd_layer is not None and not np.isnan(ct_sd_layer[i, j]) else None

            def g(col, default=None):
                return obs_row[col] if col in obs.columns else default

            rows.append({
                "experiment_id": exp_id,
                "sample_id": sample_id,
                "target_name": target_name,
                "target_type": var.loc[target_name, "target_type"]
                    if "target_type" in var.columns else None,
                "ct_mean": ct_mean_val,
                "ct_sd": ct_sd_val,
                "assay_type": uns.get("assay_type"),
                "organism": g("organism"),
                "tissue": g("tissue"),
                "cell_type": g("cell_type"),
                "condition": g("condition"),
                "biological_replicate": g("biological_replicate"),
                "batch": g("batch"),
                "sex": g("sex"),
                "passage_number": g("passage_number"),
                "plate_id": g("plate_id"),
                "cdna_input_ng": g("cdna_input_ng"),
                "cdna_concentration_ng_uL": g("cdna_concentration_ng_uL"),
                "reverse_transcription_kit": g("reverse_transcription_kit"),
                "master_mix": g("master_mix"),
                "primer_assay_id": var.loc[target_name, "primer_assay_id"]
                    if "primer_assay_id" in var.columns else None,
                "protocol_version": g("protocol_version"),
            })

    if rows:
        df = pd.DataFrame(rows)
        conn.execute("DELETE FROM observations WHERE experiment_id = ?", [exp_id])
        conn.execute("INSERT INTO observations SELECT * FROM df")

    conn.close()
    print(f"Registered experiment '{exp_id}' — "
          f"{len(adata.obs_names)} samples × {len(adata.var_names)} targets.")


def pull_data(
    db_path: str | Path | None = None,
    assay: list[str] | None = None,
    condition: list[str] | None = None,
    organism: list[str] | None = None,
    tissue: list[str] | None = None,
    sex: list[str] | None = None,
    batch: list[str] | None = None,
    experiment_id: list[str] | None = None,
    target_name: list[str] | None = None,
    target_type: list[str] | None = None,
    as_anndata: bool = False,
) -> pd.DataFrame:
    """Query the feature store. Returns a long-form DataFrame or AnnData.

    Examples
    --------
    bouncer.pull_data(assay=['qpcr'])
    bouncer.pull_data(assay=['qpcr'], condition=['treated', 'control'])
    bouncer.pull_data(target_name=['GAPDH', 'MYC'], as_anndata=True)
    """
    conn = _connect(db_path)

    clauses: list[str] = []
    params: list[Any] = []

    def _add(col: str, values: list[str] | None) -> None:
        if values:
            placeholders = ", ".join(["?"] * len(values))
            clauses.append(f"LOWER({col}) IN ({placeholders})")
            params.extend(v.lower() for v in values)

    _add("assay_type", assay)
    _add("condition", condition)
    _add("organism", organism)
    _add("tissue", tissue)
    _add("sex", sex)
    _add("batch", batch)
    _add("experiment_id", experiment_id)
    _add("target_name", target_name)
    _add("target_type", target_type)

    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
    query = f"SELECT * FROM observations {where} ORDER BY experiment_id, sample_id, target_name"

    df = conn.execute(query, params).df()
    conn.close()

    if df.empty:
        print("No data matched your query.")
        return df

    if as_anndata:
        return _df_to_anndata(df)

    return df


def _df_to_anndata(df: pd.DataFrame):
    import anndata as ad
    import numpy as np

    pivot_ct = df.pivot_table(index="sample_id", columns="target_name", values="ct_mean", aggfunc="first")
    pivot_sd = df.pivot_table(index="sample_id", columns="target_name", values="ct_sd", aggfunc="first")

    sample_ids = pivot_ct.index.tolist()
    target_names = pivot_ct.columns.tolist()

    obs_cols = [
        "experiment_id", "assay_type", "organism", "tissue", "cell_type",
        "condition", "biological_replicate", "batch", "sex", "passage_number",
        "plate_id", "cdna_input_ng", "cdna_concentration_ng_uL",
        "reverse_transcription_kit", "master_mix", "protocol_version",
    ]
    obs = (
        df.drop_duplicates("sample_id")
        .set_index("sample_id")
        [[c for c in obs_cols if c in df.columns]]
        .reindex(sample_ids)
    )
    var_cols = ["target_type", "primer_assay_id"]
    var = (
        df.drop_duplicates("target_name")
        .set_index("target_name")
        [[c for c in var_cols if c in df.columns]]
        .reindex(target_names)
    )

    adata = ad.AnnData(
        X=pivot_ct.values.astype(np.float32),
        obs=obs,
        var=var,
    )
    adata.layers["ct_sd"] = pivot_sd.reindex(
        index=sample_ids, columns=target_names
    ).values.astype(np.float32)
    adata.obs_names = sample_ids
    adata.var_names = target_names
    return adata


def list_experiments(db_path: str | Path | None = None) -> pd.DataFrame:
    conn = _connect(db_path)
    df = conn.execute("SELECT * FROM experiments ORDER BY registered_at DESC").df()
    conn.close()
    return df
