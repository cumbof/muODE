"""Worked mechanistic scenarios built from dependency-light toy models.

These tie the :mod:`muode.ecology` layers together into a complete, runnable
story that needs only numpy/scipy -- no GEMs, no CarveMe -- so the *mechanisms*
can be exercised and tested on any machine.  The genome-scale version of the same
story is the ``examples/fmt_cdiff`` workflow.

:func:`cdi_scenario` reproduces the central clinical fact about recurrent
*C. difficile* infection and why FMT cures it, as an **emergent** outcome of the
ecology layers rather than a scripted result:

* a spore reservoir survives an antibiotic course (spores are not killed);
* afterwards, in a wiped gut, the bile pool is germinant-rich and
  secondary-bile-acid-poor, so the spores **germinate and the infection recurs**;
* an FMT introduces a community that (a) competes for carbon, (b) acidifies via
  SCFA, and (c) restores 7α-dehydroxylation (cholate -> deoxycholate); the
  deoxycholate both **blocks germination** and **inhibits** any vegetative
  *C. difficile*, so the spores stay dormant and the infection clears.

The only difference between the recurrence and cure runs is the FMT injection.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from muode.antibiotic import Antibiotic
from muode.bile import BileAcidInhibition, BileAcidTransform
from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.inject import Injection
from muode.kinetics import KineticParameters
from muode.lifecycle import SporeForming
from muode.organism import LinprogOrganism
from muode.ph import WeakAcidInhibition
from muode.phage import PhageInfection

CDIFF = "C_difficile"
COMPETITOR = "donor_competitor"   # SCFA-producing carbon competitor
BAI = "donor_baiguild"            # 7α-dehydroxylating (bile-transforming) guild


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


def cdi_kinetics() -> KineticParameters:
    return KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})


def cdi_diet() -> Diet:
    """Host-supplied carbon, germinant (taurocholate) and primary bile (cholate)."""
    return Diet(
        concentrations={"glc_e": 5.0, "tca_e": 2.0, "ca_e": 2.0, "but_e": 0.0, "dca_e": 0.0},
        influx={"glc_e": 1.5, "tca_e": 0.3, "ca_e": 0.3},
        name="gut_with_bile",
    )


def cdi_community(total_biomass: float = 0.05) -> Community:
    """Recipient: vegetative *C. difficile* plus the (dormant) donor taxa.

    The donor species are members from t=0 so the bile-transform layer can refer
    to them, but they carry zero biomass until an FMT injection seeds them.
    """
    organisms = [_grower(CDIFF, yld=0.10),
                 _grower(COMPETITOR, product="but_e", yld=0.16),
                 _grower(BAI, yld=0.14)]
    return Community(organisms, {CDIFF: 1.0, COMPETITOR: 0.0, BAI: 0.0}, total_biomass=total_biomass)


def cdi_ecology(initial_spores: float = 0.05) -> EcologyModel:
    """The full mechanistic stack shared by both arms."""
    return EcologyModel([
        # SCFA acidification; C. difficile is made extra acid-sensitive.
        WeakAcidInhibition(buffer_capacity=25.0, default_ki=12.0, ki={CDIFF: 3.0}),
        # the donor bai-guild turns cholate into deoxycholate (only once present)
        BileAcidTransform(bai_producers={BAI}, vmax_bai=3.0),
        # deoxycholate inhibits vegetative C. difficile growth
        BileAcidInhibition(targets={CDIFF}, ki=0.05),
        # spore reservoir: survives antibiotics; germination gated by bile acids
        SporeForming(species={CDIFF}, initial_spores={CDIFF: initial_spores},
                     k_germination=0.6, k_sporulation=0.8, mu_stress=0.15),
        # a vancomycin-like course: strong early killing, then cleared
        Antibiotic(susceptible={CDIFF}, dose_times=tuple(float(t) for t in range(0, 10, 2)),
                   dose=3.0, half_life=2.0, emax=5.0, ec50=0.4),
    ])


def cdi_scenario(fmt: bool, t_end: float = 96.0, dt: float = 0.05,
                 fmt_time: float = 12.0, fmt_biomass: float = 0.05) -> SimulationResult:
    """Run the recurrent-CDI scenario with (``fmt=True``) or without an FMT.

    Returns the full :class:`~muode.dfba.SimulationResult`; inspect
    ``result.biomass[CDIFF]`` and ``result.spores[CDIFF]`` to see recurrence vs
    clearance, and ``result.environment`` for pH / germination-signal traces.
    """
    community = cdi_community()
    ecology = cdi_ecology()
    injections = None
    if fmt:
        injections = [Injection.from_abundances(
            fmt_time, {COMPETITOR: 0.6, BAI: 0.4}, total_biomass=fmt_biomass, name="FMT")]
    engine = DynamicFBA(t_end=t_end, dt=dt)
    return engine.run(community, cdi_diet(), cdi_kinetics(),
                      injections=injections, ecology=ecology)


# ---------------------------------------------------------------------------
# phage therapy / predation
# ---------------------------------------------------------------------------

PATHOGEN = "pathogen_bloom"       # a fast carbon grower that would otherwise dominate
COMMENSAL = "commensal"           # a slower grower it out-competes for the same carbon


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
