"""Quantitative-profile ingestion -- the upstream abundance contract.

.. warning::

   **MetaSBT itself does not produce abundances.**  MetaSBT is a genome
   *clustering / taxonomic-characterization* framework: its ``profile`` module
   reports a genome's closest genome/species/genus/family in the database, not
   how much of it is present.  Relative abundance for a community must therefore
   come from a *separate* quantitative step (read mapping / coverage).  This
   module reads such a quantitative ``{mag_id: abundance}`` table; proper ingestion
   of MetaSBT's *taxonomy* characterization (for labelling MAGs and known/unknown
   cluster assignment) is a separate contract still to be built during M1–M3
   hardening (see the dev branch at https://github.com/cumbof/MetaSBT).

The reader is deliberately *permissive*: quantitative exports vary, so it
auto-detects the identifier, abundance and (optional) taxonomy columns by common
header names and lets the caller override any of them.  Raw counts are normalised
to relative abundance.  Pass columns explicitly (``id_column=...`` /
``abundance_column=...``) when auto-detection is wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

# candidate header names (lower-cased), in priority order
_ID_HINTS = ("mag", "genome", "bin", "cluster", "mag_id", "genome_id", "input", "id", "name", "#")
_ABUNDANCE_HINTS = (
    "relative_abundance", "rel_abundance", "relative abundance", "abundance",
    "rel_ab", "fraction", "proportion", "count", "counts", "reads", "coverage", "depth",
)
_TAXONOMY_HINTS = ("taxonomy", "lineage", "classification", "taxon", "taxa", "tax")


@dataclass
class MetaSBTProfile:
    """Parsed MetaSBT profile: normalised abundances + optional taxonomy."""

    abundances: Dict[str, float]
    taxonomy: Dict[str, str] = field(default_factory=dict)
    id_column: Optional[str] = None
    abundance_column: Optional[str] = None

    def mags(self):
        return list(self.abundances)

    def to_abundance_tsv(self, path: str | Path) -> None:
        from muode.io_utils import write_abundance

        write_abundance(self.abundances, path)


def _pick(columns, hints) -> Optional[str]:
    lower = {c.lower().strip(): c for c in columns}
    for h in hints:
        if h in lower:
            return lower[h]
    # substring fallback (e.g. "estimated_relative_abundance")
    for h in hints:
        for lc, original in lower.items():
            if h in lc:
                return original
    return None


def read_metasbt_profile(
    path: str | Path,
    id_column: Optional[str] = None,
    abundance_column: Optional[str] = None,
    taxonomy_column: Optional[str] = None,
    normalize: bool = True,
    sep: Optional[str] = None,
) -> MetaSBTProfile:
    """Read a MetaSBT-style profile into a :class:`MetaSBTProfile`.

    Columns are auto-detected from common header names unless given explicitly.
    Abundances are coerced to floats; non-positive/blank values become 0 and,
    when ``normalize`` is set, the column is rescaled to sum to 1.
    """
    import pandas as pd

    path = Path(path)
    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    if df.empty or df.shape[1] < 2:
        raise ValueError(f"{path}: expected a table with an id and an abundance column")

    id_col = id_column or _pick(df.columns, _ID_HINTS) or df.columns[0]
    ab_col = abundance_column or _pick(df.columns, _ABUNDANCE_HINTS)
    if ab_col is None:
        # fall back to the first numeric column that is not the id column
        for c in df.columns:
            if c != id_col and pd.api.types.is_numeric_dtype(pd.to_numeric(df[c], errors="coerce")):
                ab_col = c
                break
    if ab_col is None:
        raise ValueError(
            f"{path}: could not find an abundance column; pass abundance_column=... "
            f"(columns: {list(df.columns)})"
        )
    tax_col = taxonomy_column or _pick(df.columns, _TAXONOMY_HINTS)

    ids = df[id_col].astype(str).str.strip()
    vals = pd.to_numeric(df[ab_col], errors="coerce").fillna(0.0).clip(lower=0.0)
    abundances: Dict[str, float] = {}
    for i, v in zip(ids, vals):
        if i and i.lower() != "nan":
            abundances[i] = abundances.get(i, 0.0) + float(v)

    total = sum(abundances.values())
    if normalize and total > 0:
        abundances = {i: v / total for i, v in abundances.items()}

    taxonomy: Dict[str, str] = {}
    if tax_col is not None:
        for i, t in zip(ids, df[tax_col].astype(str)):
            if i and t and t.lower() != "nan":
                taxonomy[i] = t.strip()

    return MetaSBTProfile(abundances, taxonomy, id_column=id_col, abundance_column=ab_col)
