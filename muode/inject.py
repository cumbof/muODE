"""Timed biomass injection -- transplants, probiotic doses, inoculation events.

A :class:`Perturbation` (:mod:`muode.perturb`) changes the *bounds* of organisms
already in a community.  An :class:`Injection` is the complementary primitive: it
changes the *state* of the community at a chosen time by introducing biomass of
one or more organisms mid-run.  This is what turns a static "who is present"
community into one that can model a discrete ecological event:

* a **fecal microbiota transplant (FMT)** -- a whole donor community delivered as
  a bolus into a dysbiotic recipient at time ``t``;
* a **probiotic dose** -- one or a few strains added at a chosen time;
* a staged **inoculation** -- seeding a second species once the first has
  conditioned the medium.

Mechanistically nothing special happens at the injection: the introduced
organisms are ordinary members of the :class:`~muode.community.Community` that
simply start at *zero* biomass and have it set to the inoculum amount when the
event fires.  The dynamic-FBA loop already idles any species below its biomass
floor, so an injected organism contributes no flux until its bolus arrives, then
competes for the shared metabolite pool like everyone else.  Whether it
*engrafts* -- and whether it displaces an incumbent (e.g. donor commensals
restoring colonisation resistance against a pathogen) -- is an emergent property
of the simulation, not something the injection hard-codes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Mapping, Tuple

if TYPE_CHECKING:  # avoid an import cycle at runtime
    from muode.community import Community


@dataclass
class Injection:
    """A timed introduction of biomass into a running simulation.

    Parameters
    ----------
    time:
        Simulated time (h) at which the bolus is added.  It fires at the first
        integration step whose time is ``>= time``; the injected biomass is
        visible in the recorded time course from that step on (a vertical jump).
    biomass:
        Mapping of organism id -> biomass added (gDW/L).  Each id must belong to
        the community being simulated (injecting into an unknown id is ignored).
        If the id is already present the amount is *added* to its current biomass
        (a boost); otherwise the dormant zero-biomass member is seeded.
    name:
        Human-readable label (e.g. ``"FMT"``).
    """

    time: float
    biomass: Dict[str, float] = field(default_factory=dict)
    name: str = "injection"

    @classmethod
    def from_abundances(
        cls,
        time: float,
        abundances: Mapping[str, float],
        total_biomass: float,
        name: str = "injection",
    ) -> "Injection":
        """Build an injection from *relative* abundances and a total inoculum.

        The abundances need not be normalised; they are scaled so the injected
        biomass sums to ``total_biomass`` (gDW/L).
        """
        total = float(sum(abundances.values())) or 1.0
        return cls(
            time=float(time),
            biomass={k: total_biomass * float(v) / total for k, v in abundances.items()},
            name=name,
        )

    def total(self) -> float:
        """Total biomass introduced by this injection (gDW/L)."""
        return float(sum(self.biomass.values()))

    def describe(self) -> str:
        return f"{self.name}@{self.time:g}h (+{len(self.biomass)} taxa, {self.total():.4g} gDW/L)"

    def to_dict(self) -> dict:
        return {"name": self.name, "time": self.time, "biomass": dict(self.biomass)}


def merge_for_injection(
    base: "Community",
    donor: "Community",
    time: float,
    name: str = "injection",
) -> Tuple["Community", Injection]:
    """Combine a ``donor`` community into a ``base`` community as a timed event.

    The returned community contains every member of ``base`` at its original
    initial biomass plus every donor-only member at *zero* initial biomass.  The
    returned :class:`Injection` introduces each donor member's biomass
    (``donor.initial_biomass()``) at ``time``.

    A donor member whose id already exists in ``base`` is *not* duplicated: the
    base organism's model is kept and the injection simply boosts that species'
    biomass at ``time`` (a realistic transplant that reinforces a strain the
    recipient already carries at low abundance).

    This is the composable building block behind an FMT example: assemble the
    dysbiotic recipient and the healthy donor as two independent communities,
    then transplant one into the other.
    """
    from muode.community import Community

    base_ids = set(base.organism_ids)
    organisms = list(base.organisms)
    for o in donor.organisms:
        if o.id not in base_ids:
            organisms.append(o)

    # base.abundances is already normalised (sums to 1); donor-only ids default
    # to 0 in Community.__post_init__, so the residents keep their exact initial
    # biomass and the donors start dormant.
    merged = Community(organisms, dict(base.abundances), total_biomass=base.total_biomass)
    event = Injection(time=float(time), biomass=dict(donor.initial_biomass()), name=name)
    return merged, event
