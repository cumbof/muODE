"""Phage therapy / predation -- a worked scenario, built from toy models.

This is an EXAMPLE, not part of the muODE package.  The reusable mechanism it
exercises -- :class:`muode.phage.PhageInfection`, a lytic predator-prey layer -- lives
in the package; the *specific story* (two toy growers competing for one carbon source,
one of them phage-susceptible) lives here, because a specific community on a specific
diet is a study, not a tool.

A fast pathogen out-blooms a commensal for a shared carbon source; a phage specific to
the pathogen crashes the bloom and releases the commensal.  The predation, the phage
amplification and the competitive rebound are all emergent -- the only difference
between the two arms is whether the phage is present.

The organism yields here are illustrative toy values, not measurements.  See
examples/fmt_cdiff for the direction this is heading: scenarios rebuilt on genome-scale
models so the yields come from stoichiometry rather than from a dial.
"""

from __future__ import annotations

from typing import Optional

from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism
from muode.phage import PhageInfection

PATHOGEN = "pathogen_bloom"       # a fast carbon grower that would otherwise dominate
COMMENSAL = "commensal"           # a slower grower it out-competes for the same carbon


def _grower(id: str, product: Optional[str] = None, yld: float = 0.12) -> LinprogOrganism:
    """A glucose grower; optionally secreting a fermentation product 1:1."""
    grow_stoich = {"glc_e": -1.0}
    reactions = [("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0)]
    exchanges = {"glc_e": "EX_glc_e"}
    if product is not None:
        grow_stoich[product] = 1.0
        reactions.append((f"EX_{product}", {product: -1.0}, 0.0, 1000.0))
        exchanges[product] = f"EX_{product}"
    reactions.append(("GROW", grow_stoich, 0.0, 1000.0))
    return LinprogOrganism(id=id, reactions=reactions, objective={"GROW": yld}, exchanges=exchanges)


def phage_predation_scenario(therapy: bool, t_end: float = 48.0, dt: float = 0.02,
                             dose_titer: float = 0.5) -> SimulationResult:
    """A fast pathogen out-blooms a commensal; a phage reins it back in.

    Two species compete for the same carbon: ``PATHOGEN`` grows faster and, left
    alone, takes over while ``COMMENSAL`` is suppressed.  With ``therapy=True`` a
    phage specific to the pathogen is present; its lytic cycle crashes the bloom
    and *releases* the commensal -- the predation, the amplification of the phage
    and the competitive rebound are all emergent, not scripted.

    The only difference between the two arms is the phage titer; inspect
    ``result.biomass`` and ``result.environment['phage[T4-like]']``.
    """
    organisms = [_grower(PATHOGEN, yld=0.30), _grower(COMMENSAL, yld=0.18)]
    community = Community(organisms, {PATHOGEN: 0.5, COMMENSAL: 0.5}, total_biomass=0.02)
    diet = Diet(concentrations={"glc_e": 8.0}, influx={"glc_e": 1.0}, name="carbon_chemostat")
    kinetics = KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})

    layers = []
    if therapy:
        layers.append(PhageInfection(
            host=PATHOGEN, name="T4-like", adsorption_rate=12.0, burst_size=60.0,
            latent_period=0.4, decay_rate=0.1, initial_titer=dose_titer))
    ecology = EcologyModel(layers)
    return DynamicFBA(t_end=t_end, dt=dt).run(community, diet, kinetics, ecology=ecology)
