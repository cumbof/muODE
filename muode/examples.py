"""A tiny, dependency-free cross-feeding community used for tests and ``muode demo``.

It encodes the canonical microbial cross-feeding motif with two species and a
shared extracellular pool:

* **A_glucose** -- a glucose specialist. It takes up glucose and secretes
  acetate as a fermentation by-product.
* **B_acetate** -- an acetate specialist. It cannot use glucose; it grows only
  on the acetate that A releases.

So B depends entirely on A through the shared metabolite pool.  This is exactly
the situation muODE is built to capture: run the community and B blooms *after*
A has produced enough acetate; remove A (a perturbation) and B suffers a
*secondary extinction* even though nothing touched B directly.

Everything here is built from :class:`~muode.organism.LinprogOrganism`, so it
runs with only numpy/scipy -- no cobra, no solver licence, no GEM files.
"""

from __future__ import annotations

from muode.community import Community
from muode.diet import Diet
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism

# biomass yields (objective weight = gDW produced per unit substrate flux)
_YIELD_A = 0.10
_YIELD_B = 0.15


def build_glucose_specialist(id: str = "A_glucose") -> LinprogOrganism:
    """Glucose -> acetate fermenter (objective ~ glucose throughput)."""
    return LinprogOrganism(
        id=id,
        reactions=[
            ("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0),  # uptake set by engine
            ("EX_ac_e", {"ac_e": -1.0}, 0.0, 1000.0),         # acetate secretion only
            ("GROW", {"glc_e": -1.0, "ac_e": 1.0}, 0.0, 1000.0),
        ],
        objective={"GROW": _YIELD_A},
        exchanges={"glc_e": "EX_glc_e", "ac_e": "EX_ac_e"},
    )


def build_acetate_specialist(id: str = "B_acetate") -> LinprogOrganism:
    """Acetate consumer; cannot use glucose, so depends on A's secretion."""
    return LinprogOrganism(
        id=id,
        reactions=[
            ("EX_ac_e", {"ac_e": -1.0}, -1000.0, 1000.0),  # acetate uptake
            ("GROW", {"ac_e": -1.0}, 0.0, 1000.0),
        ],
        objective={"GROW": _YIELD_B},
        exchanges={"ac_e": "EX_ac_e"},
    )


def toy_kinetics() -> KineticParameters:
    """Per-metabolite Michaelis-Menten constants for the toy system."""
    return KineticParameters(
        metabolite_defaults={
            "glc_e": (10.0, 0.5),  # (Vmax, Km)
            "ac_e": (10.0, 0.5),
        }
    )


def toy_diet() -> Diet:
    """A bolus of glucose and no acetate; acetate must be cross-fed from A."""
    return Diet(concentrations={"glc_e": 20.0, "ac_e": 0.0}, name="toy_glucose")


def build_toy_community(total_biomass: float = 0.02) -> Community:
    """Two-species cross-feeding community, A more abundant than B at t=0."""
    return Community(
        organisms=[build_glucose_specialist(), build_acetate_specialist()],
        abundances={"A_glucose": 0.7, "B_acetate": 0.3},
        total_biomass=total_biomass,
    )
