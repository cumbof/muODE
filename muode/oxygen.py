"""Oxygen as a shared environmental variable -- the cross-kingdom mechanism.

The gut lumen is anaerobic, yet it is *kept* anaerobic biologically: facultative
members (many yeasts and Enterobacteriaceae) consume the oxygen that diffuses in
from the mucosa, and that consumption is what lets obligate anaerobes -- the bulk
of the community -- survive.  This is the single most important reason a
multi-kingdom sample behaves differently from a bacteria-only one, because the
gut's main facultative oxygen scavengers include **fungi**.

The dynamic-FBA core does not capture this on its own: a toy model has no O2
exchange at all, and even a real GEM's O2 uptake says nothing about whether the
*species* is poisoned by oxygen or requires it.  This layer adds, on top of the
FBA core:

* **tolerance gating** of growth by the current O2 level --
    - obligate aerobe   : grows only with O2     (``o2/(km+o2)``),
    - obligate anaerobe : poisoned by O2          (``1/(1+o2/ki)``),
    - facultative / aerotolerant : indifferent    (factor 1);
* **oxygen consumption** by aerobes and facultatives (first order in their
  biomass, Monod in O2), which draws the pool down and *creates* the anaerobic
  niche -- the mechanism by which a facultative fungus protects the anaerobes.

For real GEMs that already exchange ``o2_e`` through their own FBA, set
``consume=False`` so uptake is not double-counted; the layer then contributes
only the tolerance gating.  For toy / coarse models with no O2 reaction, keep
``consume=True`` and the layer is the sole O2 bookkeeper.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping, Set

from muode.ecology import EcologyLayer
from muode.traits import Domain, OxygenTolerance

OXYGEN = "o2_e"


@dataclass
class OxygenSensitivity(EcologyLayer):
    """Oxygen tracking, tolerance gating and scavenging.

    Parameters
    ----------
    obligate_aerobes, obligate_anaerobes, facultative:
        organism-id sets by oxygen relationship.  Aerotolerant species need not
        be listed (they are unaffected).  ``facultative`` species are both
        O2-indifferent for growth *and* O2 consumers.
    oxygen, o2_initial, o2_influx:
        the O2 pool metabolite id, its starting concentration (mmol/L) and a
        constant influx (mmol/L/h) standing in for diffusion from the mucosa.
    km_o2:
        half-saturation for both aerobic growth and O2 uptake (mmol/L).
    ki_o2:
        O2 concentration that halves an obligate anaerobe's growth (mmol/L).
    consumption:
        max specific O2 uptake of aerobes/facultatives (mmol O2 / gDW / h).
    consume:
        whether the layer itself draws O2 down.  Set ``False`` when the GEMs
        already exchange O2 via FBA (avoids double counting).
    """

    name: str = "oxygen"
    obligate_aerobes: Set[str] = field(default_factory=set)
    obligate_anaerobes: Set[str] = field(default_factory=set)
    facultative: Set[str] = field(default_factory=set)
    oxygen: str = OXYGEN
    o2_initial: float = 0.25
    o2_influx: float = 0.0
    km_o2: float = 0.01
    ki_o2: float = 0.02
    consumption: float = 8.0
    consume: bool = True

    def extra_metabolites(self) -> Iterable[str]:
        return (self.oxygen,)

    def initial_concentrations(self) -> Mapping[str, float]:
        return {self.oxygen: float(self.o2_initial)}

    def _o2(self, M: Mapping[str, float]) -> float:
        return max(0.0, M.get(self.oxygen, 0.0))

    def growth_factor(self, organism_id, t, M, X) -> float:
        o2 = self._o2(M)
        if organism_id in self.obligate_aerobes:
            return o2 / (self.km_o2 + o2) if o2 > 0 else 0.0
        if organism_id in self.obligate_anaerobes:
            return 1.0 / (1.0 + o2 / self.ki_o2) if self.ki_o2 > 0 else 1.0
        return 1.0  # facultative / aerotolerant: indifferent

    def metabolite_rates(self, t, M, X) -> Mapping[str, float]:
        if not self.consume:
            return {}
        o2 = self._o2(M)
        if o2 <= 0:
            return {}
        monod = o2 / (self.km_o2 + o2)
        scavengers = self.obligate_aerobes | self.facultative
        draw = sum(self.consumption * monod * max(0.0, X.get(s, 0.0)) for s in scavengers)
        return {self.oxygen: -draw} if draw else {}

    def observables(self, t, M, X) -> Dict[str, float]:
        return {"oxygen": self._o2(M)}

    def describe(self) -> str:
        return (f"{self.name}: aerobes={sorted(self.obligate_aerobes)} "
                f"anaerobes={sorted(self.obligate_anaerobes)} "
                f"facultative={sorted(self.facultative)}")

    # -- convenience --------------------------------------------------------
    @classmethod
    def from_traits(cls, traits: Mapping[str, "object"], **kwargs) -> "OxygenSensitivity":
        """Build the species sets from a ``{organism_id: MicrobeTraits}`` map."""
        aer: Set[str] = set()
        ana: Set[str] = set()
        fac: Set[str] = set()
        for oid, tr in traits.items():
            ox = getattr(tr, "oxygen", None)
            if getattr(tr, "domain", None) is Domain.VIRUS:
                continue
            if ox is OxygenTolerance.OBLIGATE_AEROBE:
                aer.add(oid)
            elif ox is OxygenTolerance.OBLIGATE_ANAEROBE:
                ana.add(oid)
            elif ox is OxygenTolerance.FACULTATIVE:
                fac.add(oid)
        return cls(obligate_aerobes=aer, obligate_anaerobes=ana, facultative=fac, **kwargs)
