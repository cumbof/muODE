"""Bracken ingestion -- the quantitative abundance contract.

`Bracken <https://github.com/jenniferlu717/Bracken>`_ re-estimates per-taxon
read abundances from a Kraken2 classification.  In muODE the Kraken2 (and
Bracken) databases are built *from the MetaSBT database* (see
:mod:`muode.metasbt` and :mod:`muode.quantify`), so a Bracken report assigns
reads to the very species clusters MetaSBT defines -- making it the community's
quantitative profiler.

A Bracken report is a TSV with the header::

    name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads  new_est_reads  fraction_total_reads

This module reads one into ``{name: abundance}`` (``fraction_total_reads`` by
default, or the re-estimated read count), optionally filtered to a single
``taxonomy_lvl`` (``S`` = species) and renormalised over the kept rows.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

#: single-letter Bracken/Kraken rank codes -> muODE level names.
LEVEL_CODES = {
    "D": "domain", "K": "kingdom", "P": "phylum", "C": "class",
    "O": "order", "F": "family", "G": "genus", "S": "species",
}

_VALUE_COLUMNS = {
    "fraction": "fraction_total_reads",
    "reads": "new_est_reads",
    "est_reads": "new_est_reads",
    "kraken_reads": "kraken_assigned_reads",
}


@dataclass
class BrackenProfile:
    """Parsed Bracken report: per-taxon abundance keyed by ``name``."""

    abundances: Dict[str, float]
    taxonomy_ids: Dict[str, str] = field(default_factory=dict)
    level: Optional[str] = None       # taxonomy_lvl code kept (e.g. "S"), or None
    value: str = "fraction"           # which column the abundances came from

    def taxa(self):
        return list(self.abundances)

    def __len__(self) -> int:
        return len(self.abundances)


def read_bracken(
    path: str | Path,
    level: Optional[str] = "S",
    value: str = "fraction",
    normalize: bool = True,
) -> BrackenProfile:
    """Read a Bracken report into a :class:`BrackenProfile`.

    Parameters
    ----------
    level:
        Keep only rows whose ``taxonomy_lvl`` equals this code (``"S"`` for
        species, the default).  Pass ``None`` to keep every row.
    value:
        Which column to use as the abundance: ``"fraction"``
        (``fraction_total_reads``, the default), ``"reads"``/``"est_reads"``
        (``new_est_reads``) or ``"kraken_reads"`` (``kraken_assigned_reads``).
    normalize:
        Rescale the kept abundances to sum to 1.  Useful when ``level`` filters
        out some rows (the original fractions then no longer sum to 1) or when a
        read-count column is requested.
    """
    import pandas as pd

    column = _VALUE_COLUMNS.get(value)
    if column is None:
        raise ValueError(
            f"unknown value '{value}'; use one of {sorted(_VALUE_COLUMNS)}"
        )

    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    required = {"name", "taxonomy_lvl", column}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: not a Bracken report (missing columns: {sorted(missing)})"
        )

    if level is not None:
        df = df[df["taxonomy_lvl"].astype(str).str.upper() == level.upper()]

    names = df["name"].astype(str).str.strip()
    vals = pd.to_numeric(df[column], errors="coerce").fillna(0.0).clip(lower=0.0)
    abundances: Dict[str, float] = {}
    for n, v in zip(names, vals):
        if n and n.lower() != "nan":
            abundances[n] = abundances.get(n, 0.0) + float(v)

    total = sum(abundances.values())
    if normalize and total > 0:
        abundances = {n: v / total for n, v in abundances.items()}

    taxonomy_ids: Dict[str, str] = {}
    if "taxonomy_id" in df.columns:
        for n, t in zip(names, df["taxonomy_id"].astype(str)):
            if n and n.lower() != "nan":
                taxonomy_ids[n] = t.strip()

    return BrackenProfile(abundances, taxonomy_ids, level=level, value=value)
