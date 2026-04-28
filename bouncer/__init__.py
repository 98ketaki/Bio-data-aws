"""Bouncer — context-aware QC and feature registration for bioinformatics data."""

from .api import check, register, pull_data, list_experiments
from .qc.report import generate_report, print_report
from .qc.engine import QCResult

__version__ = "0.1.0"
__all__ = [
    "check",
    "register",
    "pull_data",
    "list_experiments",
    "generate_report",
    "print_report",
    "QCResult",
]
