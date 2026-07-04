"""Organism traits: domain (kingdom) and oxygen relationship.

The dynamic-FBA engine (:mod:`muode.dfba`) is **paradigm-defined, not taxon-
defined**: it only ever drives the :class:`~muode.organism.OrganismModel`
protocol, so it has no idea whether a model is a bacterium, an archaeon or a
fungus.  A genome-scale model of *Candida albicans* satisfies the same interface
as one of *Bacteroides* and simulates with the *same* code.

What does differ between kingdoms is (a) **how the GEM is reconstructed** -- the
default CarveMe path is prokaryote-only -- and (b) a handful of **physiological
traits** that the metabolic LP does not encode on its own, chiefly the oxygen
relationship.  This module captures those as lightweight, optional metadata so
the rest of muODE (reconstruction routing, the :mod:`muode.oxygen` layer) can act
on them without the simulation core having to know anything about taxonomy.

Viruses (phages) are included in :class:`Domain` for completeness, but a phage is
**not** an :class:`OrganismModel` -- it has no metabolism.  It is modelled as a
coupled infection process in :mod:`muode.phage`, not as a member of the FBA
community.  See ``docs/kingdoms.md`` and ``docs/LIMITATIONS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class Domain(str, Enum):
    """Top-level lineage of a community member."""

    BACTERIA = "bacteria"
    ARCHAEA = "archaea"
    EUKARYOTE = "eukaryote"   # fungi, protists
    VIRUS = "virus"           # phages -- not an FBA organism (see muode.phage)


class OxygenTolerance(str, Enum):
    """Relationship to molecular oxygen.

    This is the trait that most often makes "kingdom" matter dynamically: the gut
    lumen is anaerobic, and facultative members (many yeasts, Enterobacteriaceae)
    are its principal oxygen scavengers -- their consumption keeps the niche
    anaerobic for the obligate anaerobes.  Consumed by :mod:`muode.oxygen`.
    """

    OBLIGATE_ANAEROBE = "obligate_anaerobe"   # poisoned by O2
    AEROTOLERANT = "aerotolerant"             # indifferent to O2, ferments anyway
    FACULTATIVE = "facultative"               # uses O2 when present, ferments when not
    OBLIGATE_AEROBE = "obligate_aerobe"       # requires O2


@dataclass
class MicrobeTraits:
    """Optional per-organism metadata layered beside its metabolic model.

    Attributes
    ----------
    domain:
        which kingdom the organism belongs to (drives reconstruction routing).
    oxygen:
        its oxygen relationship (drives the :mod:`muode.oxygen` layer).
    cell_mass_pg:
        approximate dry mass of one cell (pg) -- eukaryotic cells are ~10-100x
        heavier than bacteria; recorded for documentation / unit reasoning only.
    notes:
        free-form provenance (e.g. the GEM source or reference).
    """

    domain: Domain = Domain.BACTERIA
    oxygen: OxygenTolerance = OxygenTolerance.OBLIGATE_ANAEROBE
    cell_mass_pg: Optional[float] = None
    notes: str = ""

    @property
    def is_metabolic(self) -> bool:
        """Whether this member is an FBA organism (everything except a virus)."""
        return self.domain is not Domain.VIRUS


# ---------------------------------------------------------------------------
# reconstruction routing
# ---------------------------------------------------------------------------

#: recommended GEM-reconstruction route per domain.  ``engine`` is the muODE
#: reconstruction engine name where one applies; ``note`` is shown to the user.
_RECONSTRUCTION = {
    Domain.BACTERIA: ("carveme", "CarveMe (BiGG namespace); universe gramneg/grampos."),
    Domain.ARCHAEA: ("carveme", "CarveMe with `universe=archaea`."),
    Domain.EUKARYOTE: (
        None,
        "CarveMe does NOT support eukaryotes. muODE routes them separately: "
        "reconstruct_mag(engine='carvefungi') for fungi, or engine='eukaryote_generic' "
        "for other eukaryotes (MetaEuk gene calling + CarveFungi / eggNOG-mapper+"
        "ModelSEEDpy, gated by euk_ref_db). Or supply a curated template (e.g. Yeast8 "
        "for ascomycetes) and load the SBML into the community as a CobraOrganism.",
    ),
    Domain.VIRUS: (
        None,
        "A phage has no metabolism and no GEM. Do not reconstruct it; model it as "
        "a muode.phage.PhageInfection layer coupled to its host species.",
    ),
}


def reconstruction_route(domain: Domain) -> tuple[Optional[str], str]:
    """Return ``(engine_or_None, human_note)`` describing how to get a model.

    Used by the assembly/reconstruction guidance so that a non-bacterial member
    is routed correctly (or flagged as not-an-FBA-organism) instead of being
    silently handed to CarveMe, which would fail or produce nonsense.
    """
    return _RECONSTRUCTION[Domain(domain)]
