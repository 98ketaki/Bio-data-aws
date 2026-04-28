"""Single QC finding (one issue found by one check)."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(str, Enum):
    HARD = "hard"
    SOFT = "soft"
    WARNING = "warning"


@dataclass
class Finding:
    check: str
    severity: Severity
    section: str          # ct_export | ntc | reference_genes | samplesheet | cross_file | design
    message: str
    wells: list[str] = None   # affected well positions, if applicable
    samples: list[str] = None # affected sample_ids

    def __post_init__(self):
        if self.wells is None:
            self.wells = []
        if self.samples is None:
            self.samples = []

    @property
    def blocks_registration(self) -> bool:
        return self.severity == Severity.HARD
