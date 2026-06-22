"""The microbial community: a set of organism models sharing one environment.

muODE follows the *compartmentalised / shared-pool* community formulation used
by dynamic community-FBA tools (e.g. COMETS): every species keeps its own
genome-scale model and they interact only through a common pool of extracellular
metabolites.  This avoids materialising one giant merged LP, parallelises
trivially across species, and maps one-to-one onto the README's community ODEs

    dX_i/dt = mu_i * X_i
    dM_j/dt = sum_i v_(j,i) * X_i

A :class:`Community` therefore holds:

* the per-species :class:`~muode.organism.OrganismModel` objects,
* their relative abundances (e.g. from a MetaSBT profile), and
* the initial biomass per species (total biomass split by abundance).

It is solver-agnostic: organisms may be ``CobraOrganism`` (real GEMs) or
``LinprogOrganism`` (toy/test models).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional

from muode.organism import OrganismModel


@dataclass
class Community:
    organisms: List[OrganismModel]
    #: relative abundance per organism id (need not be normalised; it will be)
    abundances: Dict[str, float] = field(default_factory=dict)
    #: total community biomass at t=0 (gDW/L)
    total_biomass: float = 0.01

    def __post_init__(self) -> None:
        ids = [o.id for o in self.organisms]
        if len(set(ids)) != len(ids):
            raise ValueError(f"duplicate organism ids in community: {ids}")
        if not self.abundances:
            self.abundances = {i: 1.0 for i in ids}
        # keep only known organisms and normalise
        self.abundances = {i: float(self.abundances.get(i, 0.0)) for i in ids}
        total = sum(self.abundances.values())
        if total <= 0:
            self.abundances = {i: 1.0 / len(ids) for i in ids}
        else:
            self.abundances = {i: v / total for i, v in self.abundances.items()}

    # -- views --------------------------------------------------------------
    @property
    def organism_ids(self) -> List[str]:
        return [o.id for o in self.organisms]

    def get(self, organism_id: str) -> OrganismModel:
        for o in self.organisms:
            if o.id == organism_id:
                return o
        raise KeyError(organism_id)

    def initial_biomass(self) -> Dict[str, float]:
        """Initial biomass per species = total biomass * relative abundance."""
        return {i: self.total_biomass * a for i, a in self.abundances.items()}

    def environment_metabolites(self) -> List[str]:
        """All extracellular metabolites any species can exchange."""
        mets: set[str] = set()
        for o in self.organisms:
            mets.update(o.exchange_metabolites())
        return sorted(mets)

    # -- construction -------------------------------------------------------
    @classmethod
    def from_models(
        cls,
        model_dir: str | Path,
        abundance: Optional[Mapping[str, float]] = None,
        total_biomass: float = 0.01,
        pattern: str = "*.xml",
    ) -> "Community":
        """Build a community from a directory of SBML GEMs.

        Requires the ``cobra`` extra.  The organism id is the file stem, which
        must match the ids used in the abundance profile.
        """
        from muode.organism import CobraOrganism

        model_dir = Path(model_dir)
        paths = sorted(model_dir.glob(pattern))
        if not paths:
            raise FileNotFoundError(f"no models matching {pattern} in {model_dir}")
        organisms: List[OrganismModel] = [
            CobraOrganism.from_file(p, id=p.stem) for p in paths
        ]
        return cls(organisms, dict(abundance or {}), total_biomass)
