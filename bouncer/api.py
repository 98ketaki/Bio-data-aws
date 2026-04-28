"""Bouncer public API — check, register, pull_data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .qc.engine import QCResult, run_qc
from .qc.report import print_report, generate_report
from .features.builder import build_anndata
from .features.store import (
    register as _store_register,
    pull_data as _store_pull,
    list_experiments as _store_list,
)


def check(
    ct_export: str | Path,
    samplesheet: str | Path,
    schema: str | Path,
    qc_contract: str | Path,
    protocol: str | Path | None = None,
    verbose: bool = True,
) -> QCResult:
    """Run all QC checks and print a report.

    Parameters
    ----------
    ct_export:    Path to QuantStudio Ct export CSV.
    samplesheet:  Path to plate samplesheet CSV.
    schema:       Path to bouncer schema YAML.
    qc_contract:  Path to bouncer QC contract YAML.
    protocol:     Path to protocol PDF (optional, used for context).
    verbose:      Print the report immediately.

    Returns
    -------
    QCResult with .passed, .hard_failures, .soft_failures, .warnings.
    """
    result = run_qc(
        ct_export_path=ct_export,
        samplesheet_path=samplesheet,
        schema_path=schema,
        qc_contract_path=qc_contract,
        protocol_path=protocol,
    )
    if verbose:
        print_report(result)
    return result


def register(
    result: QCResult,
    db_path: str | Path | None = None,
    force: bool = False,
) -> "anndata.AnnData":
    """Convert a QC-passed result to AnnData and push to the feature store.

    Parameters
    ----------
    result:   A QCResult returned by bouncer.check().
    db_path:  Path to the DuckDB file. Defaults to bouncer_store.duckdb next to the package.
    force:    If True, register even with soft failures (hard failures always block).

    Returns
    -------
    The registered AnnData object.
    """
    if not result.passed:
        failures = "\n".join(f"  - {f.check}: {f.message}" for f in result.hard_failures)
        raise RuntimeError(
            f"Registration blocked — {len(result.hard_failures)} hard failure(s):\n{failures}\n"
            "Fix the issues above and re-run bouncer.check()."
        )

    adata = build_anndata(result)
    _store_register(adata, db_path=db_path)
    return adata


def pull_data(
    assay: list[str] | None = None,
    condition: list[str] | None = None,
    organism: list[str] | None = None,
    tissue: list[str] | None = None,
    sex: list[str] | None = None,
    batch: list[str] | None = None,
    experiment_id: list[str] | None = None,
    target_name: list[str] | None = None,
    target_type: list[str] | None = None,
    db_path: str | Path | None = None,
    as_anndata: bool = False,
) -> pd.DataFrame:
    """Pull data from the feature store with optional filters.

    Examples
    --------
    >>> import bouncer
    >>> df = bouncer.pull_data(assay=['qpcr'])
    >>> df = bouncer.pull_data(assay=['qpcr'], condition=['treated', 'control'])
    >>> adata = bouncer.pull_data(target_name=['MYC'], as_anndata=True)
    """
    return _store_pull(
        db_path=db_path,
        assay=assay,
        condition=condition,
        organism=organism,
        tissue=tissue,
        sex=sex,
        batch=batch,
        experiment_id=experiment_id,
        target_name=target_name,
        target_type=target_type,
        as_anndata=as_anndata,
    )


def list_experiments(db_path: str | Path | None = None) -> pd.DataFrame:
    """List all registered experiments in the feature store."""
    return _store_list(db_path=db_path)
