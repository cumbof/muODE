"""Rich, guild-structured recurrent-CDI / FMT scenario builder.

The companion to :mod:`muode.scenarios.cdi` (the minimal 3-member mechanistic
sketch). This module defines a **defined 12-guild gut community** together with
the full composable ecology stack used to study recurrent *Clostridioides
difficile* infection (rCDI) and its resolution: antibiotic-driven collapse and
spore-driven relapse, FMT rescue, rational design of a defined therapeutic
consortium, and the mechanistic bottlenecks of failed transplants.

Like the rest of :mod:`muode.scenarios`, every member is a dependency-light
:class:`~muode.organism.LinprogOrganism`, so the whole study runs on any machine
with only numpy/scipy -- no GEMs, no solver licence. Swap the members for
``CobraOrganism`` GEMs reconstructed from real MAGs (see ``examples/fmt_cdiff``)
to run the identical analysis on genome-scale models; the engine and ecology
stack are unchanged.

The worked, plotted analysis that drives these builders lives in
``examples/fmt_cdiff/dynamics/``.

Everything is expressed in **days** (rates, dosing times, the time step).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.inject import Injection
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism
from muode.ph import WeakAcidInhibition
from muode.bile import BileAcidTransform, BileAcidInhibition
from muode.lifecycle import SporeForming
from muode.antibiotic import Antibiotic
from muode.traits import Domain, MicrobeTraits, OxygenTolerance


# ===========================================================================
# 1. THE SHARED EXTRACELLULAR POOL (metabolite identifiers)
# ===========================================================================
# muODE couples members solely through a shared pool of extracellular
# metabolites (Methods, "Community members and dynamics"). We use BiGG-style
# `*_e` identifiers so a single namespace is enforced across every member --
# the prerequisite for emergent cross-feeding.

MUCIN = "mucin_e"     # complex mucosal glycans / dietary polysaccharide (host-supplied)
GLC = "glc_e"         # monomeric sugars liberated by primary degraders (shared carbon currency)
PRO = "pro_e"         # proline  -- a Stickland-fermentation substrate
GLY = "gly_e"         # glycine  -- a Stickland-fermentation substrate
AC = "ac_e"           # acetate  -- a short-chain fatty acid (SCFA), cross-fed and acidifying
BUT = "but_e"         # butyrate -- a SCFA (colonocyte fuel), acidifying
SUCC = "succ_e"       # succinate
LAC = "lac__L_e"      # L-lactate -- a SCFA, cross-fed

# Bile acids. The primary -> secondary transformation is performed by the ecology
# layer (BileAcidTransform), keyed on guild membership, NOT by metabolic
# reactions -- exactly as in muode.bile. We only need to name the pool species.
TCA = "tca_e"         # taurocholate  -- conjugated primary bile acid (spore germinant)
CA = "ca_e"           # cholate       -- free primary bile acid
DCA = "dca_e"         # deoxycholate  -- secondary bile acid (germination + growth inhibitor)
LCA = "lca_e"         # lithocholate  -- secondary bile acid (inhibitor)

PRIMARY_BILE = (TCA, CA)
SECONDARY_BILE = (DCA, LCA)

# The pathogen id is referenced by nearly every ecology layer.
CDIFF = "Clostridioides_difficile"


# ===========================================================================
# 2. MEMBER FACTORY -- a small mass-balanced metabolic model per organism
# ===========================================================================
# Every community member is a tiny stoichiometric model solved with
# scipy/HiGHS via muode.organism.LinprogOrganism. The genome-scale
# GEMs would slot in here unchanged (they satisfy the same OrganismModel API).
#
# A member is defined by:
#   * a set of *uptake* metabolites (exchange lower bound < 0: import allowed),
#   * a set of *product* metabolites (exchange lower bound = 0: secretion only),
#   * one or more *growth modes*, each a reaction that consumes substrates and
#     (optionally) secretes products, carrying a biomass yield (objective weight
#     = grams dry weight produced per unit reaction flux).
#
# With several growth modes the LP picks the most productive feasible one given
# the current medium, so substrate preference (e.g. C. difficile favouring
# Stickland fermentation) emerges rather than being hard-coded.

# Flux convention (COBRA): an exchange reaction `met_e -->` (stoichiometry
# {met: -1}) carries positive flux for SECRETION and negative flux for UPTAKE.


@dataclass
class GrowthMode:
    """One biomass-producing reaction: {metabolite: signed coefficient} + yield.

    Negative coefficients are consumed substrates; positive coefficients are
    secreted products. ``yield_`` is the objective weight (gDW per unit flux).
    """

    stoich: Dict[str, float]
    yield_: float


def make_member(
    member_id: str,
    modes: Sequence[GrowthMode],
    uptakes: Iterable[str],
    products: Iterable[str] = (),
) -> LinprogOrganism:
    """Assemble a :class:`LinprogOrganism` from growth modes and exchange lists.

    Parameters
    ----------
    member_id:
        Organism identifier (used everywhere: abundances, ecology layers, plots).
    modes:
        One or more :class:`GrowthMode`. Each becomes a ``GROW_<k>`` reaction and
        contributes its yield to the biomass objective.
    uptakes:
        Metabolites the member may import (exchange lower bound = -1000).
    products:
        Metabolites the member may only secrete (exchange lower bound = 0). A
        metabolite that is both consumed and secreted should be listed in
        ``uptakes`` (import is the more permissive bound).
    """
    uptakes = list(dict.fromkeys(uptakes))      # de-duplicate, keep order
    products = [p for p in dict.fromkeys(products) if p not in uptakes]

    reactions: List[tuple] = []
    exchanges: Dict[str, str] = {}

    # Boundary/exchange reactions. Convention: EX_<met> has stoichiometry {met:-1}.
    for met in uptakes:
        reactions.append((f"EX_{met}", {met: -1.0}, -1000.0, 1000.0))   # import allowed
        exchanges[met] = f"EX_{met}"
    for met in products:
        reactions.append((f"EX_{met}", {met: -1.0}, 0.0, 1000.0))       # secretion only
        exchanges[met] = f"EX_{met}"

    # Growth (biomass) reactions -- one per mode -- plus the weighted objective.
    objective: Dict[str, float] = {}
    for k, mode in enumerate(modes):
        rxn_id = f"GROW_{k}"
        reactions.append((rxn_id, dict(mode.stoich), 0.0, 1000.0))
        objective[rxn_id] = float(mode.yield_)

    return LinprogOrganism(
        id=member_id, reactions=reactions, objective=objective, exchanges=exchanges
    )


# ===========================================================================
# 3. THE GUILD REGISTRY  (the 12-member consortium + the pathogen)
# ===========================================================================
# Each entry records the functional guild, oxygen relationship, domain, and the
# metabolic wiring that realises the guild's stated role. The "role in pool"
# comments tie each member back to its functional-guild rationale.

@dataclass
class GuildSpec:
    """Static description of one community member (its biology + wiring)."""

    member_id: str
    guild: str                              # human-readable functional guild
    modes: List[GrowthMode]
    uptakes: List[str]
    products: List[str] = field(default_factory=list)
    domain: Domain = Domain.BACTERIA
    oxygen: OxygenTolerance = OxygenTolerance.OBLIGATE_ANAEROBE
    bsh: bool = False                       # carries bile-salt hydrolase?
    bai: bool = False                       # carries the bai operon (7a-dehydroxylation)?
    spore_former: bool = False
    vanco_susceptible: bool = False
    note: str = ""

    def build(self) -> LinprogOrganism:
        return make_member(self.member_id, self.modes, self.uptakes, self.products)


# --- yields (gDW per unit growth-reaction flux). Tuned so that fast guilds reach
# --- steady state within a few simulated days and the pathogen can bloom in an
# --- empty lumen but is out-competed in a restored one. ---------------------
Y_DEGRADER = 0.20
Y_BAI = 0.22
Y_SCAVENGER = 0.30       # aggressive Stickland competitor (depletes pro/gly fast)
Y_BUTYRATE = 0.20
Y_SINK = 0.18
Y_FUNGUS = 0.12
Y_CDIFF_STICKLAND = 0.34  # pathogen strongly prefers Stickland fermentation
Y_CDIFF_SUGAR = 0.18


# ---- The pathogen ---------------------------------------------------------
CDIFF_SPEC = GuildSpec(
    member_id=CDIFF,
    guild="Pathogen (spore-forming; Stickland fermenter)",
    # Two growth modes: preferred Stickland fermentation (proline + glycine as a
    # coupled electron donor/acceptor pair) and a fallback on monomeric sugar.
    modes=[
        GrowthMode({PRO: -1.0, GLY: -1.0}, Y_CDIFF_STICKLAND),
        GrowthMode({GLC: -1.0}, Y_CDIFF_SUGAR),
    ],
    uptakes=[PRO, GLY, GLC],
    spore_former=True,
    vanco_susceptible=True,
    note="Germinates on primary bile (taurocholate); blocked by secondary bile; "
         "acid-sensitive; killed (vegetative only) by vancomycin.",
)


# ---- The 12-member designed consortium -------------------------------------
# Each member declares the extracellular pool it targets; the guild set
# is realised by the growth modes / products below.
CONSORTIUM_SPECS: List[GuildSpec] = [
    # -- Secondary bile-acid effectors (the therapeutic core) ---------------
    GuildSpec(
        "Clostridium_scindens", "Secondary effector (7a-dehydroxylation)",
        modes=[GrowthMode({GLC: -1.0, AC: 1.0}, Y_BAI)],
        uptakes=[GLC], products=[AC], bai=True,
        note="Taurocholate->deoxycholate. Needs monomeric sugar from degraders.",
    ),
    GuildSpec(
        "Extibacter_muris", "Secondary effector (bai + bile-salt hydrolase)",
        modes=[GrowthMode({GLC: -1.0, AC: 1.0}, Y_BAI)],
        uptakes=[GLC], products=[AC], bai=True, bsh=True,
        note="Supplementary 7a-dehydroxylation + BSH; cholate->lithocholate.",
    ),
    # -- Primary polysaccharide degraders (the carbon pipeline) -------------
    GuildSpec(
        "Bacteroides_thetaiotaomicron", "Primary degrader (mucosal glycans)",
        modes=[GrowthMode({MUCIN: -1.0, GLC: 0.8}, Y_DEGRADER)],
        uptakes=[MUCIN], products=[GLC],
        note="Broad glycan degradation; secretes monomeric sugars.",
    ),
    GuildSpec(
        "Bacteroides_ovatus", "Primary degrader (dietary polysaccharide)",
        modes=[GrowthMode({MUCIN: -1.0, GLC: 0.8}, Y_DEGRADER)],
        uptakes=[MUCIN], products=[GLC],
        note="Complex polysaccharide (e.g. xylan) degradation -> monomeric sugars.",
    ),
    GuildSpec(
        "Akkermansia_muciniphila", "Primary degrader (mucin)",
        modes=[GrowthMode({MUCIN: -1.0, GLC: 0.5, AC: 0.4}, Y_DEGRADER)],
        uptakes=[MUCIN], products=[GLC, AC],
        note="Mucin degradation; supports the basal trophic chain during fasting.",
    ),
    # -- Metabolic sink / trophic bridges -----------------------------------
    GuildSpec(
        "Blautia_hydrogenotrophica", "Metabolic sink (acetogen)",
        modes=[GrowthMode({GLC: -1.0, AC: 1.0}, Y_SINK)],
        uptakes=[GLC], products=[AC],
        note="Acetogenesis proxy; prevents thermodynamic inhibition, feeds acetate.",
    ),
    # -- Amino-acid scavenger (direct pathogen competitor) ------------------
    GuildSpec(
        "Peptostreptococcus_anaerobius", "Amino-acid scavenger (Stickland)",
        modes=[GrowthMode({PRO: -1.0, GLY: -1.0}, Y_SCAVENGER)],
        uptakes=[PRO, GLY],
        note="Aggressive Stickland fermenter; depletes free proline/glycine, "
             "directly starving C. difficile of its preferred substrates.",
    ),
    # -- Trophic bridges / butyrate producers -------------------------------
    GuildSpec(
        "Eubacterium_hallii", "Trophic bridge (acetate/lactate -> butyrate)",
        modes=[GrowthMode({AC: -1.0, BUT: 1.0}, Y_BUTYRATE),
               GrowthMode({LAC: -1.0, BUT: 1.0}, Y_BUTYRATE)],
        uptakes=[AC, LAC], products=[BUT],
        note="Consumes acetate/lactate, secretes butyrate (colonocyte fuel).",
    ),
    GuildSpec(
        "Faecalibacterium_prausnitzii", "Trophic bridge (primary butyrate producer)",
        modes=[GrowthMode({GLC: -1.0, BUT: 1.0}, Y_BUTYRATE)],
        uptakes=[GLC], products=[BUT], oxygen=OxygenTolerance.AEROTOLERANT,
        note="Primary butyrate producer; high-affinity oxygen scavenger.",
    ),
    GuildSpec(
        "Roseburia_intestinalis", "Trophic bridge (xylan -> butyrate)",
        modes=[GrowthMode({GLC: -1.0, BUT: 1.0}, Y_BUTYRATE)],
        uptakes=[GLC], products=[BUT],
        note="Xylan degradation and robust butyrate production.",
    ),
    GuildSpec(
        "Bifidobacterium_longum", "Trophic bridge (oligosaccharide fermenter)",
        modes=[GrowthMode({GLC: -1.0, AC: 0.7, LAC: 0.5}, Y_BUTYRATE)],
        uptakes=[GLC], products=[AC, LAC],
        note="Specialised oligosaccharide fermentation; secretes acetate/lactate.",
    ),
    # -- Modulator ----------------------------------------------------------
    GuildSpec(
        "Parabacteroides_distasonis", "Modulator (succinate producer)",
        modes=[GrowthMode({GLC: -1.0, SUCC: 1.0}, Y_SINK)],
        uptakes=[GLC], products=[SUCC],
        note="Succinate production; broad secondary-bile-acid tolerance.",
    ),
]

# Antibiotic effect. Although oral vancomycin's *direct* spectrum is Gram-
# positive, clinical microbiome studies of vancomycin-treated CDI report a broad,
# catastrophic collapse of the gut community -- Bacteroidetes and Verrucomicrobia
# included -- not merely a loss of Clostridiales. The narrative
# depends on exactly this: a depauperate post-antibiotic lumen "saturated with
# unconsumed mucosal carbohydrates and free amino acids" that the pathogen then
# exploits. We therefore model the 10-day course as broadly suppressive of the
# resident bacterial community (the resistant fungal opportunist persists), which
# is what makes the post-antibiotic FMT-arm lumen equivalent to the depauperate
# consortium-arm lumen -- so that the transplant is the *only* difference from
# the untreated relapse arm.
for _s in CONSORTIUM_SPECS:
    _s.vanco_susceptible = True
_VANCO_SUSCEPTIBLE_GUILDS = {s.member_id for s in CONSORTIUM_SPECS}

#: The 12 designed-consortium member ids, in guild order.
CONSORTIUM_12: List[str] = [s.member_id for s in CONSORTIUM_SPECS]

#: Vancomycin-susceptible commensal guilds (the "susceptible obligate anaerobes",
#: e.g. Ruminococcaceae / Lachnospiraceae, whose crash the relapse arm tracks).
SUSCEPTIBLE_COMMENSALS: List[str] = [
    s.member_id for s in CONSORTIUM_SPECS if s.vanco_susceptible
]

#: The "minimal" naive consortium: the bai-positive effectors only.
MINIMAL_CONSORTIUM: List[str] = [
    s.member_id for s in CONSORTIUM_SPECS if s.bai
]

#: The primary polysaccharide degraders (the carbon pipeline whose absence is the
#: "donor insufficiency" bottleneck of the failed-FMT autopsy).
PRIMARY_DEGRADERS: List[str] = [
    "Bacteroides_thetaiotaomicron", "Bacteroides_ovatus", "Akkermansia_muciniphila",
]


# ---- Optional multi-kingdom incumbent: an antibiotic-tolerant fungus ------
# The pre-intervention rCDI landscape is enriched for
# "Candida amino-acid fermenters capable of surviving high-antibiotic regimes".
# We include one as a benign incumbent (present in the recipient, not vancomycin-
# susceptible) so the multi-kingdom expansion is represented. It ferments sugar
# and amino acids weakly and never dominates.
FUNGUS_SPEC = GuildSpec(
    member_id="Candida_albicans",
    guild="Fungal opportunist (sugar fermenter)",
    # A weak glucose fermenter. We deliberately do NOT let it consume the
    # pathogen's Stickland substrates (proline/glycine): as a vancomycin-tolerant
    # incumbent it would otherwise scavenge that proline and spuriously starve
    # C. difficile, masking the loss of colonisation resistance. It therefore
    # only nibbles the scarce free glucose and never dominates.
    modes=[GrowthMode({GLC: -1.0}, Y_FUNGUS)],
    uptakes=[GLC],
    domain=Domain.EUKARYOTE, oxygen=OxygenTolerance.FACULTATIVE,
    vanco_susceptible=False,
    note="Eukaryotic opportunist that persists through the antibiotic course.",
)


# All specs, indexed by id, for convenient lookup by the scenario scripts.
ALL_SPECS: Dict[str, GuildSpec] = {
    s.member_id: s for s in [CDIFF_SPEC, FUNGUS_SPEC, *CONSORTIUM_SPECS]
}


# ===========================================================================
# 4. COMMUNITY CONSTRUCTION
# ===========================================================================
# One master community is used throughout. The pathogen and fungal incumbent are
# always present at t=0. The twelve commensal guild members are ALWAYS part of
# the roster (so the ecology layers can reference them and a timed bolus can seed
# them), but their starting biomass is a parameter:
#
#   * ``commensal_level > 0``  -> a "pre-antibiotic patient microbiome": the
#     commensal guilds are present and functional. The vancomycin course then
#     crashes the susceptible (Firmicutes) fraction -- this is the setup for the
#     relapse and FMT-rescue scenarios.
#   * ``commensal_level == 0`` -> a "depauperate post-antibiotic lumen": only the
#     pathogen (and fungus) are present; every commensal is dormant and can only
#     appear via a timed injection. This is the setup for the designed-
#     consortium design (a naive bai-only bolus starves; the full 12-member bolus
#     builds a self-sustaining trophic chain).
#
# A dormant member (0 biomass) contributes no flux -- muODE idles anything below
# its biomass floor -- until a bolus raises it, after which it competes normally.

CDIFF_INITIAL = 0.04           # pathogen biomass at t=0 (gDW/L)
FUNGUS_INITIAL = 0.02          # fungal incumbent biomass at t=0 (gDW/L)
COMMENSAL_LEVEL = 0.60         # total resident commensal biomass, relapse/FMT arms (gDW/L)


def build_master_community(
    commensal_level: float = 0.0, include_fungus: bool = True
) -> Community:
    """Assemble the recipient community used by every scenario.

    Parameters
    ----------
    commensal_level:
        Total biomass (gDW/L) of the resident commensal guilds at t=0, split
        evenly across the twelve members. ``0`` gives a depauperate lumen
        (commensals dormant, seeded only by injection); a positive value gives a
        functional pre-antibiotic microbiome that the vancomycin course crashes.
    include_fungus:
        Include the antibiotic-tolerant fungal incumbent (multi-kingdom flavour).
    """
    incumbents = [CDIFF_SPEC] + ([FUNGUS_SPEC] if include_fungus else [])
    donors = CONSORTIUM_SPECS
    organisms = [s.build() for s in (incumbents + donors)]

    # Absolute initial biomass per member (gDW/L). We pass these as "abundances"
    # together with total_biomass = their sum, so Community's internal
    # normalise-then-rescale returns exactly these absolute values.
    per_commensal = (commensal_level / len(donors)) if commensal_level > 0 else 0.0
    biomass: Dict[str, float] = {CDIFF: CDIFF_INITIAL}
    if include_fungus:
        biomass[FUNGUS_SPEC.member_id] = FUNGUS_INITIAL
    for s in donors:
        biomass[s.member_id] = per_commensal

    total = sum(biomass.values())
    traits = {
        s.member_id: MicrobeTraits(domain=s.domain, oxygen=s.oxygen)
        for s in (incumbents + donors)
    }
    return Community(organisms, biomass, total_biomass=total, traits=traits)


# ===========================================================================
# 5. THE DEFINED GUT MEDIUM  (Methods: "defined in-silico gut medium")
# ===========================================================================
# A Western-diet-like lumen supplemented with host-derived mucin glycans and
# primary bile acids (cholate and taurocholate), plus the Stickland substrates
# proline and glycine. Concentrations are mmol/L; influx (open-system renewal,
# e.g. host secretion) is mmol/L/day.

def gut_medium() -> Diet:
    """The defined gut medium shared by every scenario."""
    return Diet(
        concentrations={
            MUCIN: 8.0,     # complex glycans -- only degraders can use them
            GLC: 1.0,       # little free sugar; most must be liberated by degraders
            PRO: 4.0, GLY: 4.0,     # Stickland substrates (host/dietary protein)
            TCA: 3.0, CA: 3.0,      # primary bile acids (host-derived germinants)
            DCA: 0.0, LCA: 0.0,     # secondary bile acids -- made by the bai guild
            AC: 0.0, BUT: 0.0, LAC: 0.0, SUCC: 0.0,
        },
        influx={
            MUCIN: 1.0,             # continuous host mucin turnover
            PRO: 0.4, GLY: 0.4,     # continuous protein supply
            TCA: 0.4, CA: 0.4,      # continuous hepatic primary-bile secretion
        },
        name="defined_gut_medium",
    )


# ===========================================================================
# 6. UPTAKE KINETICS  (Michaelis-Menten, per day)
# ===========================================================================
# Saturable uptake couples each member's LP to the shared pool. Capacities are
# mmol/gDW/day; affinities (Km) are mmol/L. High capacity + small Km => uptake
# saturates quickly once a substrate is present, so nutrient depletion and
# competitive exclusion emerge from the shared pool.

def cdi_kinetics() -> KineticParameters:
    """Per-metabolite Michaelis-Menten constants for the study (day units)."""
    vmax, km = 14.0, 0.5
    mets = [MUCIN, GLC, PRO, GLY, AC, BUT, LAC, SUCC]
    return KineticParameters(
        default_vmax=vmax, default_km=km,
        metabolite_defaults={m: (vmax, km) for m in mets},
    )


# ===========================================================================
# 7. THE COMPOSABLE ECOLOGY LAYER
# ===========================================================================
# The non-metabolic processes that drive rCDI: SCFA acidification, the two-step
# bile-acid transformation, secondary-bile inhibition, sporulation/germination,
# and (for the antibiotic scenarios) time-varying vancomycin pharmacokinetics.
# Each layer is a real muODE EcologyLayer; they compose without order-dependence.

# --- Antibiotic (vancomycin) pharmacokinetics/dynamics ---------------------
# 10-day course: one dose per day on days 0..9. Single-compartment PK with
# first-order elimination (half-life in days), driving a saturating (Emax) kill
# rate on *vegetative* susceptible biomass. Spores are insulated (SporeForming).
ANTIBIOTIC_DAYS = 10
VANCOMYCIN_DOSE_TIMES: Tuple[float, ...] = tuple(float(d) for d in range(ANTIBIOTIC_DAYS))


def vancomycin_layer(susceptible: Iterable[str]) -> Antibiotic:
    """A 10-day vancomycin course targeting the susceptible vegetative biomass."""
    return Antibiotic(
        name="vancomycin",
        susceptible=set(susceptible),
        dose_times=VANCOMYCIN_DOSE_TIMES,
        dose=4.0,               # per-dose luminal concentration (arbitrary consistent unit)
        half_life=0.4,          # days (luminal residence)
        emax=8.0,               # maximum kill rate (1/day)
        ec50=0.5,               # concentration at half-maximal killing
    )


def build_ecology(community: Community, with_antibiotic: bool = True) -> EcologyModel:
    """Assemble the full ecology stack for a given community.

    The bai / bsh / spore-former / susceptibility sets are read straight off the
    :class:`GuildSpec` registry, so the ecology always matches the members that
    are present.
    """
    present = set(community.organism_ids)

    def ids(pred) -> set:
        return {mid for mid in present if mid in ALL_SPECS and pred(ALL_SPECS[mid])}

    bai_producers = ids(lambda s: s.bai)
    bsh_producers = ids(lambda s: s.bsh)
    spore_formers = ids(lambda s: s.spore_former)
    susceptible = ids(lambda s: s.vanco_susceptible)

    layers: List = [
        # (1) SCFA acidification + weak-acid inhibition. This is a *supporting*
        #     arm of colonisation resistance here: the central rCDI mechanism in
        #     the model is the secondary-bile shield plus Stickland-substrate
        #     starvation (layers 3-4 and the amino-acid scavenger), so we keep the
        #     pathogen only mildly acid-sensitive (large Ki) and use a well-
        #     buffered medium. Acidification still makes dense fermentation self-
        #     limiting and modestly disfavours the pathogen.
        WeakAcidInhibition(
            buffer_capacity=60.0, default_ki=20.0, ki={CDIFF: 20.0},
            acids=(AC, BUT, LAC),
        ),
        # (2) Two-step microbial bile-acid transformation in the shared pool:
        #     taurocholate --BSH--> cholate --7a-dehydroxylase--> deoxycholate.
        #     Driven by the biomass of the guilds carrying each enzyme.
        #     The Vmax values are high because the colon is a flowing system
        #     (see DILUTION): the transformed pool is continuously washed out, so
        #     a thriving guild must dehydroxylate briskly to hold a steady-state
        #     secondary-bile shield. A washed-out / starving guild produces
        #     essentially none -- exactly the fail/cure separation we want.
        BileAcidTransform(
            bsh_producers=bsh_producers, bai_producers=bai_producers,
            vmax_bsh=25.0, vmax_bai=20.0, tca=TCA, ca=CA, dca=DCA,
        ),
        # (3) Secondary bile acids inhibit vegetative C. difficile growth. The
        #     inhibition is deliberately *graded* (a mmol-scale Ki, not a razor-
        #     edge one): a mere trickle of secondary bile from a starving bai
        #     guild is not enough to clear the pathogen. Durable clearance needs a
        #     *thriving*, degrader-fed bai guild (a large secondary-bile pool)
        #     AND depletion of the pathogen's Stickland substrates by the amino-
        #     acid scavenger -- the "dual constraint". This is why
        #     the naive bai-only consortium fails while the full
        #     12-member consortium succeeds.
        BileAcidInhibition(targets={CDIFF}, ki=1.0, secondary=SECONDARY_BILE),
        # (4) Sporulation / germination. Spores survive antibiotics; germination
        #     is triggered by primary bile (taurocholate) and blocked by
        #     secondary bile acids -- coupling the bile pool to the pathogen life
        #     cycle. Rates are per day.
        SporeForming(
            species=spore_formers,
            initial_spores={CDIFF: 0.06},
            germinant=TCA, km_germinant=0.15,
            inhibitor=SECONDARY_BILE, ki_inhibitor=0.03,
            # Germination is deliberately slow (per day): the spore reservoir must
            # survive the 10-day vancomycin course (spores are drug-insulated but
            # any that germinate into vegetative cells are killed) and only bloom
            # once the drug clears (~Day 12) in a primary-bile-rich, secondary-
            # bile-poor lumen. Sporulation under nutrient stress is fast.
            k_germination=0.12, k_sporulation=1.5, mu_stress=0.3,
            # Spores are shed from the flowing lumen at a slow rate. This makes
            # the cure a genuine extinction: once the vegetative pathogen is gone
            # and no fresh spores are made, the residual spore pool washes out
            # rather than persisting forever. In the relapse arm, continuous
            # re-sporulation from the active bloom keeps the reservoir topped up.
            k_spore_decay=0.06,
        ),
    ]
    if with_antibiotic and susceptible:
        layers.append(vancomycin_layer(susceptible))

    return EcologyModel(layers)


# ===========================================================================
# 8. SIMULATION TIMELINE  (days)
# ===========================================================================
T_END = 42.0          # total simulated time (days): 10-day course + 8-week-ish follow-up
DT = 0.02             # integration step (days) -- well under the depletion/CFL scale
FMT_DAY = 12.0        # timed intervention: drug has decayed below EC50 by ~Day 12

# Chemostat-style washout (1/day). The colon is an open, flowing system: biomass
# and luminal metabolites are continuously diluted by transit. Without it, the
# constant nutrient influx would drive biomass and bile acids up without bound.
# A modest dilution (~5-day residence) bounds every field to a steady state,
# which is both more physiological and far more legible. It applies to the
# vegetative pool and the extracellular pool; the (adherent) spore reservoir is
# not washed out.
DILUTION = 0.2


# ===========================================================================
# 9. DATA HELPERS
# ===========================================================================


def biomass_sum(result, member_ids: Iterable[str]):
    """Total biomass over time across a set of members (a pandas Series)."""
    cols = [m for m in member_ids if m in result.biomass.columns]
    return result.biomass[cols].sum(axis=1) if cols else result.biomass.iloc[:, 0] * 0.0


# ===========================================================================
# 10. SCENARIO RUNNER + TIME-SERIES HELPERS
# ===========================================================================
# A single entry point that assembles a community, its ecology stack, and an
# optional timed bolus, then runs the muODE dynamic-FBA engine. Every figure
# script calls this so the exact same engine/ecology configuration is used
# throughout -- the only thing that changes between scenarios is (a) whether the
# resident commensals start present, and (b) which members (if any) are injected
# at FMT_DAY.

def run_scenario(
    commensal_level: float,
    inject_members: Optional[Sequence[str]] = None,
    with_antibiotic: bool = True,
    inoculum_biomass: float = 0.05,
    extra_layers: Optional[Sequence] = None,
    t_end: float = T_END,
    dt: float = DT,
) -> Tuple[SimulationResult, Optional[Antibiotic]]:
    """Build and run one rCDI scenario with the real muODE engine.

    Parameters
    ----------
    commensal_level:
        Resident commensal biomass at t=0 (see :func:`build_master_community`).
    inject_members:
        Member ids delivered as a timed bolus (FMT / consortium) at ``FMT_DAY``,
        seeded in proportion to an even relative abundance. ``None`` = no bolus.
    with_antibiotic:
        Apply the 10-day vancomycin course.
    inoculum_biomass:
        Total biomass (gDW/L) of the injected bolus.
    extra_layers:
        Additional :class:`~muode.ecology.EcologyLayer` objects to compose on top
        (e.g. a :class:`~muode.phage.PhageInfection` for the failed-FMT autopsy).

    Returns
    -------
    (result, antibiotic_layer) -- the :class:`SimulationResult` and the
    :class:`~muode.antibiotic.Antibiotic` layer (or ``None``), the latter so the
    drug-concentration curve can be reconstructed for plotting.
    """
    community = build_master_community(commensal_level=commensal_level)
    ecology = build_ecology(community, with_antibiotic=with_antibiotic)
    for layer in (extra_layers or []):
        ecology.add(layer)

    injections = None
    if inject_members:
        injections = [Injection.from_abundances(
            FMT_DAY, {m: 1.0 for m in inject_members},
            total_biomass=inoculum_biomass, name="bolus",
        )]

    result = DynamicFBA(t_end=t_end, dt=dt, dilution_rate=DILUTION).run(
        community, gut_medium(), cdi_kinetics(),
        injections=injections, ecology=ecology,
    )
    antibiotic = next(
        (L for L in ecology.layers if isinstance(L, Antibiotic)), None
    )
    return result, antibiotic


def primary_bile(result):
    """Summed primary bile-acid concentration over time (pandas Series)."""
    m = result.metabolites
    return m[[x for x in PRIMARY_BILE if x in m.columns]].sum(axis=1)


def secondary_bile(result):
    """Summed secondary bile-acid concentration over time (pandas Series)."""
    m = result.metabolites
    return m[[x for x in SECONDARY_BILE if x in m.columns]].sum(axis=1)


def drug_curve(antibiotic: Optional[Antibiotic], times) -> "list":
    """Reconstruct the luminal drug concentration C(t) at each recorded time."""
    if antibiotic is None:
        return [0.0 for _ in times]
    return [antibiotic.concentration(float(t)) for t in times]

