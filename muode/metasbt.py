"""MetaSBT ingestion -- taxonomic *characterization* of MAGs (not abundance).

.. note::

   **MetaSBT does not produce abundances.**  MetaSBT (https://github.com/cumbof/MetaSBT)
   is a framework for *clustering genomes into known and unknown taxonomic
   clusters* (species clusters, and also higher ranks).  Its ``profile`` module
   *characterizes* a genome: for each taxonomic level it reports the closest
   cluster in the database (its label, the ANI distance and a confidence).  That
   is taxonomic **identity**, not **quantity**.

   Relative abundance comes from a *separate* quantitative step.  In muODE the
   plan is: MetaSBT database -> custom Kraken2 DB (``metasbt kraken``) -> custom
   Bracken DB -> Bracken read assignment, parsed by :mod:`muode.bracken`.  The
   two are joined in :mod:`muode.quantify`: a MAG's MetaSBT species cluster is
   matched to its Bracken abundance.

This module parses MetaSBT ``profile`` outputs.  One file per genome
(``{genome}.txt`` or ``{genome}.split.txt``) with a header
``# level<TAB>closest<TAB>ani<TAB>confidence`` and one row per taxonomic level
(plus a ``genome`` row for the single closest genome).
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

#: the seven standard taxonomic levels MetaSBT characterizes, coarse -> fine.
LEVELS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")

#: column order of a MetaSBT profile table (used when the header is absent).
_COLUMNS = ("level", "closest", "ani", "confidence")


@dataclass(frozen=True)
class ProfileMatch:
    """A single ``profile`` row: the closest cluster at one taxonomic level."""

    level: str
    closest: str
    ani: Optional[float] = None
    confidence: Optional[float] = None


@dataclass
class MetaSBTProfile:
    """A genome's MetaSBT characterization: closest cluster per taxonomic level.

    ``matches`` is keyed by level name; it may also contain a ``genome`` entry
    for the single closest genome (MetaSBT emits one).  Use the level
    properties (:attr:`species`, :attr:`lineage`, ...) for the common joins.
    """

    genome_id: str
    matches: Dict[str, ProfileMatch] = field(default_factory=dict)

    # -- per-level accessors ------------------------------------------------
    def label(self, level: str = "species") -> Optional[str]:
        """Closest cluster label at ``level`` (``None`` if not characterized)."""
        m = self.matches.get(level)
        return m.closest if m else None

    @property
    def species(self) -> Optional[str]:
        return self.label("species")

    @property
    def genus(self) -> Optional[str]:
        return self.label("genus")

    @property
    def closest_genome(self) -> Optional[str]:
        return self.label("genome")

    @property
    def lineage(self) -> "OrderedDict[str, str]":
        """Ordered ``{level: closest}`` for the standard ranks present."""
        return OrderedDict(
            (lvl, self.matches[lvl].closest) for lvl in LEVELS if lvl in self.matches
        )

    @property
    def lineage_string(self) -> str:
        """Pipe-separated lineage (e.g. ``k__Bacteria|...|s__E_coli``)."""
        return "|".join(self.lineage.values())

    def confidence(self, level: str = "species") -> Optional[float]:
        m = self.matches.get(level)
        return m.confidence if m else None


def _coerce_float(value: str) -> Optional[float]:
    value = (value or "").strip()
    if not value or value.lower() in ("na", "nan", "none", "-"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _genome_id_from_path(path: Path) -> str:
    name = path.name
    for suffix in (".split.txt", ".profile.txt", ".txt", ".tsv", ".profile"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_metasbt_profile(
    path: str | Path,
    genome_id: Optional[str] = None,
) -> MetaSBTProfile:
    """Parse a single MetaSBT ``profile`` file into a :class:`MetaSBTProfile`.

    The header (``# level<TAB>closest<TAB>ani<TAB>confidence``) is honoured when
    present to locate columns; otherwise the default order is assumed.  Blank
    and comment-only lines are ignored.  ``genome_id`` defaults to the file stem
    (with MetaSBT's ``.split``/``.txt`` suffixes stripped).
    """
    path = Path(path)
    if genome_id is None:
        genome_id = _genome_id_from_path(path)

    header: Optional[List[str]] = None
    matches: "OrderedDict[str, ProfileMatch]" = OrderedDict()

    for raw in path.read_text().splitlines():
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        is_comment = line.lstrip().startswith("#")
        cells = [c.strip() for c in line.lstrip("#").strip().split("\t")]
        if is_comment:
            lowered = [c.lower() for c in cells]
            if "level" in lowered and "closest" in lowered:
                header = lowered
            continue  # any other comment line is metadata, skip it

        cols = header or list(_COLUMNS)
        row = dict(zip(cols, cells))
        level = (row.get("level") or (cells[0] if cells else "")).strip().lower()
        closest = (row.get("closest") or (cells[1] if len(cells) > 1 else "")).strip()
        if not level or not closest:
            continue
        matches[level] = ProfileMatch(
            level=level,
            closest=closest,
            ani=_coerce_float(row.get("ani", "")),
            confidence=_coerce_float(row.get("confidence", "")),
        )

    return MetaSBTProfile(genome_id=genome_id, matches=matches)


def read_metasbt_profiles(
    source: str | Path,
    pattern: str = "*.txt",
) -> Dict[str, MetaSBTProfile]:
    """Read every MetaSBT profile in a directory (or a single file).

    Returns ``{genome_id: MetaSBTProfile}``.  ``*.split.txt`` files (MetaSBT's
    split variant) take precedence over a plain ``*.txt`` for the same genome.
    """
    source = Path(source)
    if source.is_file():
        prof = read_metasbt_profile(source)
        return {prof.genome_id: prof}

    profiles: Dict[str, MetaSBTProfile] = {}
    # plain files first, then let .split.txt override (it is the refined call)
    files = sorted(source.glob(pattern), key=lambda p: (".split." in p.name, p.name))
    for f in files:
        prof = read_metasbt_profile(f)
        profiles[prof.genome_id] = prof
    return profiles


def mag_taxa(
    profiles: Dict[str, MetaSBTProfile],
    level: str = "species",
) -> Dict[str, str]:
    """Map ``{genome_id: closest cluster label}`` at ``level`` (skipping misses)."""
    out: Dict[str, str] = {}
    for gid, prof in profiles.items():
        lab = prof.label(level)
        if lab:
            out[gid] = lab
    return out
