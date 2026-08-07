"""The rCDI / FMT scenario on REAL genome-scale models.

The point of this scenario is that nothing in the community's *outcome* is a dial.
Each of the five members is a gapseq genome-scale reconstruction, so its growth yield
is the biomass stoichiometry that came out of its genome -- there is nothing left to
tune.  When ``C_difficile`` is out-competed here, it is because the reconstructions say
so, or it is not.  (This is the fix for the circularity ``tests/test_provenance.py``
records: an earlier sketch of this scenario ran on toy organisms whose yields were
hand-set to produce clearance -- ``_grower(CDIFF, yld=0.10)`` -- and no amount of
parameter provenance rescues a yield chosen to give the answer.  That sketch is gone.)

Why gapseq / ModelSEED and not CarveMe / BiGG
---------------------------------------------
CarveMe carves from the BiGG universe, which has no *connected* butyrate pathway for gut
anaerobes: the CarveMe Roseburia and F. prausnitzii carried no butyrate exchange at all,
so the SCFA / acidification arm of colonisation resistance could not emerge from
stoichiometry -- it would have had to be a dial too.  gapseq builds bottom-up from
MetaCyc and, on the *same* genomes, reconstructs models that secrete butyrate and acetate
natively (``EX_cpd00211_e0``, ``EX_cpd00029_e0``).  So the whole community is gapseq, in
one namespace (ModelSEED), and the members cross-feed with no id translation.  muODE
keeps CarveMe as its default engine; this scenario chooses gapseq because its result
rests on a pathway CarveMe cannot supply (see ``gems/PROVENANCE.md``).

Readiness is by namespace, not filenames
-----------------------------------------
A member is *ready* only if its GEM loads AND is in the ModelSEED namespace.  This is
what lets :func:`build_community` refuse to assemble a "CDI scenario" that carries no
*C. difficile* -- reporting an infection cleared from a community that never contained
the pathogen would be the worst kind of fake, and this module will not produce one.
All five members are reconstructed and in ``gems/`` (see ``gems/PROVENANCE.md``); the
readiness gate remains so the refusal is enforced, not assumed.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
GEMS = HERE / "gems"

# Run against the checkout this example ships in, not whatever ``muode`` happens to be
# installed.  A script's sys.path[0] is its own directory, so without this the repo root
# is absent and an older installed muode wins -- silently, and only where it differs
# (an outdated Diet, a missing ecology layer).  A checkout is authoritative for its own
# examples; if muode is not vendored here, fall through to the installed package.
_ROOT = HERE.parents[1]
if (_ROOT / "muode" / "__init__.py").exists() and str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from muode.community import Community  # noqa: E402
from muode.diet import Diet  # noqa: E402
from muode.organism import CobraOrganism  # noqa: E402

#: ModelSEED ids for the metabolites the scenario reasons about.  (The bile pool is an
#: ecology-layer concern in muODE ids, namespace-independent -- see :func:`cdi_diet` --
#: so it is deliberately NOT here: gapseq GEMs need no bile exchanges for the bile layer.)
BUTYRATE = "cpd00211_e0"
ACETATE = "cpd00029_e0"

#: How we tell a gapseq/ModelSEED GEM from a leftover CarveMe/BiGG one: ModelSEED
#: extracellular metabolite ids look like ``cpd00001_e0``.  Water is in every model.
MODELSEED_WATER = "cpd00001_e0"


@dataclass(frozen=True)
class Member:
    """One community member and the claim its reconstruction has to back."""

    slug: str
    species: str
    accession: str
    role: str
    #: ModelSEED metabolite this member is claimed to SECRETE from stoichiometry, or
    #: None if its role is not defined by a secreted product (the degrader, the
    #: pathogen).  verify() checks the claim rather than trusting it.
    secretes: Optional[str] = None
    is_pathogen: bool = False

    def gem_path(self) -> Optional[Path]:
        """This member's GEM, or None if it has not been reconstructed yet.

        The filename does not certify the reconstructor -- every model here is gapseq
        (see ``gems/PROVENANCE.md``), so naming one of them ``.gapseq.`` would only
        imply the others are something else.  What a file *is* is settled by
        :func:`load`, which rejects anything outside the ModelSEED namespace; this
        only has to find it.
        """
        p = GEMS / f"{self.slug}.xml.gz"
        return p if p.exists() else None


#: The designed 5-member FMT study.  Roles map to the mechanisms in :func:`cdi_ecology`:
#: carbon competition (all, from the FBA core), butyrate/acidification (the two
#: producers), and 7a-dehydroxylation (C. scindens, via markers.py -- the bai step is an
#: ecology layer, not an FBA reaction, so it does not need to be in the GEM).
MEMBERS: List[Member] = [
    Member("C_difficile_630", "Clostridioides difficile 630", "GCF_000009205.2",
           "THE PATHOGEN: Stickland fermenter; spore-former", is_pathogen=True),
    Member("R_intestinalis_L182", "Roseburia intestinalis L1-82", "GCA_900537995.1",
           "butyrate producer -- the SCFA/acidification arm", secretes=BUTYRATE),
    Member("F_prausnitzii_A2165", "Faecalibacterium prausnitzii A2-165", "GCA_002734145.1",
           "butyrate producer + published-GEM cross-check", secretes=BUTYRATE),
    Member("B_thetaiotaomicron_VPI5482", "Bacteroides thetaiotaomicron VPI-5482",
           "GCF_000011065.1", "generalist carbohydrate degrader -- broad carbon competition"),
    Member("C_scindens_ATCC35704", "Clostridium scindens ATCC 35704", "GCA_004295125.1",
           "bai+ 7a-dehydroxylase (cholate -> deoxycholate), via markers.py"),
]

MODELSEED_DIET = GEMS / "western_gut_modelseed.csv"


def load_modelseed_diet() -> Diet:
    """The western_gut medium translated into the gapseq (ModelSEED) namespace."""
    return Diet.from_csv(MODELSEED_DIET)


def _is_modelseed(org: CobraOrganism) -> bool:
    return MODELSEED_WATER in set(org.exchange_metabolites())


def load(member: Member) -> Optional[CobraOrganism]:
    """Load a member as a :class:`CobraOrganism`, or None if it is not ready.

    "Ready" means: a GEM file exists AND it is in the ModelSEED namespace.  A leftover
    CarveMe/BiGG file loads fine but is rejected here -- it is the wrong namespace for
    this community, and treating it as ready would mix BiGG and ModelSEED ids in one
    shared pool, so nobody would cross-feed.
    """
    path = member.gem_path()
    if path is None:
        return None
    org = CobraOrganism.from_file(str(path), id=member.slug)
    if not _is_modelseed(org):
        return None
    return org


def readiness() -> Dict[str, str]:
    """Per-member status, for a human deciding whether the scenario can run yet.

    ``ready`` | ``awaiting gapseq (have only old CarveMe/BiGG)`` | ``not reconstructed yet``.
    """
    status: Dict[str, str] = {}
    for m in MEMBERS:
        path = m.gem_path()
        if path is None:
            status[m.slug] = "not reconstructed yet"
        elif load(m) is None:
            status[m.slug] = "awaiting gapseq (have only old CarveMe/BiGG)"
        else:
            status[m.slug] = "ready"
    return status


@dataclass
class MemberReport:
    """What a reconstruction actually does -- the yields that used to be dials."""

    slug: str
    growth: float                 # /h on the ModelSEED western_gut, WITH flux bounds
    secretes: Optional[str]       # the metabolite it was claimed to secrete
    #: Max export flux of that product when the pathway is exercised (medium open).
    #: This is a property of the RECONSTRUCTION -- does it encode a connected pathway
    #: at all -- and is the claim that separates gapseq from CarveMe.  It is NOT the
    #: on-diet secretion rate: secretion is growth-coupled, so a fastidious member that
    #: barely grows on this (incomplete) medium secretes ~0 while still carrying the
    #: pathway.  We report the capability, and the on-diet growth, as separate facts.
    can_secrete: Optional[float]
    pathway_present: Optional[bool]   # None if no claim; else does the pathway carry flux?


def verify(member: Member, diet: Optional[Diet] = None) -> MemberReport:
    """Does this reconstruction grow, and does it encode the pathway its role claims?

    TWO SEPARATE FACTS, deliberately not conflated:

    * ``growth`` -- the rate through the diet's flux bounds (``growth_on_diet``), a
      physiological number, not the open-medium fiction.  F. prausnitzii ~1e-3/h
      (fastidious, consistent with Clark); R. intestinalis ~0.03/h.

    * ``can_secrete`` -- for a member whose role is a secreted product, the MAXIMUM
      export flux the reconstruction can carry when the pathway is exercised.  This is
      the claim that justified moving to gapseq: CarveMe's Roseburia/F. prausnitzii had
      NO butyrate exchange, so this would be None; gapseq's carry a connected pathway,
      so it is large.  Crucially this is *capability*, measured with the objective set
      to the export reaction -- it does not depend on the member growing well on the
      current, still-incomplete ModelSEED medium.
    """
    from muode.media import growth_on_diet

    diet = diet or load_modelseed_diet()
    org = load(member)
    if org is None:
        raise FileNotFoundError(
            f"{member.slug} is not ready ({readiness()[member.slug]}); cannot verify it"
        )

    model = org.model
    mu = growth_on_diet(model, diet)

    cap: Optional[float] = None
    present: Optional[bool] = None
    if member.secretes is not None:
        ex = f"EX_{member.secretes}"
        exch_ids = {r.id for r in model.exchanges}
        if ex not in exch_ids:
            present = False        # claims a product it has no exchange for at all
        else:
            # Exercise the pathway: open what the diet supplies, then MAXIMISE export.
            # Capability, not the growth-coupled rate.
            supplied = {f"EX_{met}": diet.uptake_limit(met) or 10.0
                        for met in diet.metabolites()
                        if diet.initial_concentration(met) > 0 and f"EX_{met}" in exch_ids}
            with model:
                model.medium = supplied
                model.objective = ex
                cap = float(model.slim_optimize() or 0.0)
            present = cap > 1e-6

    return MemberReport(member.slug, mu, member.secretes, cap, present)


def build_community(
    diet: Optional[Diet] = None,
    abundances: Optional[Dict[str, float]] = None,
    total_biomass: float = 0.05,
    require_pathogen: bool = True,
) -> Community:
    """Assemble the community from the members whose gapseq GEMs are ready.

    **Refuses to build a CDI community without the pathogen.**  A scenario that reports
    a *C. difficile* infection cleared while carrying no *C. difficile* would be the
    worst kind of fake, so this raises rather than quietly dropping it.  Pass
    ``require_pathogen=False`` only to inspect the donor community in isolation.
    """
    # Check the pathogen FIRST, and cheaply -- load() short-circuits to None the moment
    # the file is absent, so this parses at most one GEM, not the whole community.  If we
    # cannot run the scenario there is no point loading the donors, and a CDI study that
    # silently drops C. difficile is exactly the failure this gate exists to prevent.
    if require_pathogen:
        pathogen = next(m for m in MEMBERS if m.is_pathogen)
        if load(pathogen) is None:
            raise RuntimeError(
                f"cannot assemble the CDI community: the pathogen {pathogen.slug} is "
                f"{readiness()[pathogen.slug]}. A CDI scenario without C. difficile is not "
                "a CDI scenario. Run gems/reconstruct_gapseq.sh on the workstation, or "
                "pass require_pathogen=False to inspect the donor community alone."
            )

    organisms = [org for org in (load(m) for m in MEMBERS) if org is not None]
    return Community(organisms, abundances or {}, total_biomass=total_biomass)


#: The pathogen and the four donors, by the roles the ecology layers key on.  This
#: membership is genome-derived, not tuned: C. scindens ATCC 35704 is the type strain
#: for 7a-dehydroxylation (the bai arm), C. difficile is the spore-former and the
#: vancomycin target.  (Ideally bai/bsh membership comes from markers.py run on the
#: proteomes; C. scindens is declared here because it is definitional for this strain.)
PATHOGEN = "C_difficile_630"

#: Both guilds are read off the MODELS' OWN ANNOTATIONS, not asserted here -- the same
#: discipline the ModelSEED diet follows.  Searching the five GEMs for the relevant EC
#: numbers gives exactly one answer each:
#:
#:   bai 7a-dehydratase (EC 4.2.1.106) -> C. scindens only          (rxn05066_c0, ...)
#:   BSH               (EC 3.5.1.24)  -> R. intestinalis, F. prausnitzii  (rxn02795_c0, ...)
#:
#: B. thetaiotaomicron carries NO bile-acid reaction at all, which is worth knowing:
#: the obvious guess (Bacteroides deconjugate bile) is wrong for THIS genome, and the
#: guess is what a hand-written guild would have encoded.
#:
#: BSH was previously never passed to BileAcidTransform at all, so bsh_producers was an
#: empty set: taurocholate was never deconjugated, the tca -> ca -> dca cascade was
#: broken at step 1, and the germinant accumulated unopposed (influx, no consumer, no
#: washout) to ~31 mM by t=96 -- 15x physiological.  The bai arm was acting only on the
#: 2 mM of cholate the diet happened to supply directly.
BAI_GUILD = {"C_scindens_ATCC35704"}
BSH_GUILD = {"R_intestinalis_L182", "F_prausnitzii_A2165"}
DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
          "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]


#: Colonic washout, 1/h.  DERIVED, and the one number here with a solid basis: colonic
#: transit is 24-48 h, so D = 1/transit = 0.021-0.042/h; 0.025 is a 40 h transit.
#:
#: This was 0.0 -- i.e. the colon was modelled as a sealed batch culture.  That is not a
#: conservative simplification, it is an actively wrong one that compounds with time:
#: anything the community does not consume integrates without bound.  Taurocholate
#: (influx 0.3, no exchange, and -- until BSH_GUILD was wired up -- no consumer) reached
#: ~31 mM by t=96, 15x physiological, dragging the germination signal from 0.11 to 0.66
#: on an artifact.  With washout it instead sits at a steady state of influx/D.
#:
#: It also makes colonization resistance EXPRESSIBLE.  At D=0, any organism with mu>0
#: persists forever and no community can exclude anything; with washout, a member has to
#: outgrow transit to stay, which is what resistance to colonization means.  Every member
#: here clears D=0.025 on its own in monoculture except C. scindens and F. prausnitzii
#: (0.0008/h), which must be cross-fed to hold their place at all -- so this rate is not
#: a passive parameter, it is what decides who the community can carry.
DILUTION_RATE = 0.025


#: Bile-acid concentration in the caecum, mM.  Ramirez & Abel-Santos 2011.
CAECAL_BILE_MM = 2.0


def cdi_diet(diet: Optional[Diet] = None,
             dilution_rate: float = DILUTION_RATE) -> Diet:
    """ModelSEED western_gut PLUS the bile-acid pool the ecology layers act on.

    The bile ids (tca_e/ca_e/dca_e) are muODE ecology-layer ids, NOT ModelSEED: the
    bile transformation is a phenomenological pool process keyed on guild membership,
    so the GEMs need no bile exchanges.

    The influx is DERIVED from the concentration rather than picked: at steady state an
    unconsumed pool sits at ``influx / dilution_rate``, so supplying it at ``C * D``
    holds it at the physiological ~2 mM instead of letting it drift.  It used to be a
    flat 0.3 mmol/L/h, which -- with no washout and (through the bsh_producers bug) no
    consumer -- meant taurocholate ramped to ~31 mM by t=96, 15x physiological, dragging
    the germination signal from 0.11 to 0.66 on nothing but an accumulation artifact.
    The germinant's own affinity is measured (km 15.9 mM); letting its concentration
    float was quietly overriding that measurement.
    """
    from muode.bile import CHOLATE, TAUROCHOLATE

    base = diet or load_modelseed_diet()
    conc = dict(base.concentrations)
    conc[TAUROCHOLATE] = CAECAL_BILE_MM   # the germinant
    conc[CHOLATE] = CAECAL_BILE_MM        # primary bile: the bai substrate
    influx = dict(base.influx)
    influx[TAUROCHOLATE] = CAECAL_BILE_MM * dilution_rate
    influx[CHOLATE] = CAECAL_BILE_MM * dilution_rate
    return Diet(concentrations=conc, influx=influx, max_uptake=dict(base.max_uptake),
                source=dict(base.source), name="modelseed_western_gut+bile")


def cdi_ecology(ablate: str = "", vmax_bai=None):
    """The mechanistic stack.  ``ablate`` selects which arms are present.

    There are two kinds of number here, and the distinction is the point.

    **Library parameters** -- germination km (15.9 mM, MEASURED), inhibition ki (0.5 mM,
    DERIVED), mu_stress (0.01/h, ASSUMED), k_sporulation / k_germination (INVENTED) --
    are NOT passed.  They come from the layer defaults, which are the registered values
    in :mod:`muode.provenance`, so there is exactly one place to audit them.

    This docstring used to claim that ALL parameters worked that way, while quietly
    passing eleven overrides.  One of them decided the study: ``mu_stress=0.15`` is over
    3x faster than any member of this community can grow on a diet bounded by measured
    intake (the fastest manages 0.042/h), so ``growth < mu_stress`` was true at every
    step for every organism and the pathogen sporulated unconditionally -- every arm
    reported CLEARED, the untreated control included.  See
    :data:`muode.lifecycle.MU_STRESS`.

    **Scenario parameters** stay explicit below, because they describe THIS clinical
    story rather than the mechanism, and no library default could be right for them:

    * the vancomycin course (which drug, what dose, how often);
    * the colonic titration (``buffer_capacity``) and how acid-sensitive THIS pathogen
      is (``ki``).

    All of them are INVENTED -- nobody measured them for this system, and the honest
    treatment of an invented mechanism is to make it removable and show what it does.
    Each is covered by an ``ablate`` arm (``abx``, ``ph``), so its contribution is
    measured rather than asserted.  If an arm turns out to carry the result, that number
    needs evidence before the result can be published.
    """
    from muode.antibiotic import Antibiotic
    from muode.bile import BileAcidInhibition, BileAcidTransform
    from muode.ecology import EcologyModel
    from muode.lifecycle import SporeForming
    from muode.ph import WeakAcidInhibition

    layers = []
    if "ph" not in ablate:
        # INVENTED, both of them: buffer_capacity 25 (vs the library's generic 60) says
        # colonic content titrates more sharply than a lab medium, and ki 3.0 (vs the
        # default 12.0) says C. difficile is ~4x more acid-sensitive than the community
        # around it.  Directionally supported -- SCFAs do inhibit C. difficile -- but
        # neither number is measured for this system.  The `ph` arm is what tests them.
        layers.append(WeakAcidInhibition(buffer_capacity=25.0, default_ki=12.0,
                                         ki={PATHOGEN: 3.0}))
    if "bile" not in ablate:
        # vmax_bai default (1.0) is a library INVENTED param; an override lets us anchor
        # it to the measured in-vivo deoxycholate concentration rather than assert it.
        bat_kw = {} if vmax_bai is None else {"vmax_bai": float(vmax_bai)}
        layers.append(BileAcidTransform(bai_producers=set(BAI_GUILD),
                                        bsh_producers=set(BSH_GUILD), **bat_kw))
        layers.append(BileAcidInhibition(targets={PATHOGEN}))
    layers.append(SporeForming(species={PATHOGEN}, initial_spores={PATHOGEN: 0.05}))
    if "abx" not in ablate:
        # Vancomycin 125 mg PO qid is the standard rCDI course; dosed here at t=0..8 h
        # so the drug is cleared (half-life 2 h) long before t_end and any rebound is
        # the community's doing, not the drug's.
        layers.append(Antibiotic(susceptible={PATHOGEN},
                                 dose_times=tuple(float(t) for t in range(0, 10, 2)),
                                 dose=3.0, half_life=2.0, emax=5.0, ec50=0.4))
    return EcologyModel(layers)




def build_scenario(fmt: bool, ablate: str = "", t_end: float = 120.0, dt: float = 0.05,
                   fmt_time: float = 12.0, fmt_biomass: float = 0.05,
                   dilution_rate: float = DILUTION_RATE, vmax_bai=None):
    """Run the rCDI/FMT scenario on the real gapseq community.

    The recipient carries the (bloomed) pathogen; the four donors are members from t=0
    at zero biomass so the bile layer can refer to them, and an FMT injection seeds them
    at ``fmt_time``.  The ONLY difference between the treatment and control arms is the
    injection.

    ``ablate`` ("", "bile", "ph", "abx", or any space-separated combination) drops
    ecology layers so the mechanism can be decomposed.  ``abx`` matters most: with the
    drug in every arm there is no untreated baseline, so "cleared" measures what
    vancomycin did, not what the community did.

    ``t_end`` is 120 h, about 3 colonic transits -- long enough to reach the washout
    steady state and to see any post-drug rebound (the course ends at t=8 and the drug
    is gone by t~24), and NOT long enough to pretend we are simulating the weeks over
    which real rCDI recurs.  That timescale is governed by immune recovery and mucosal
    refuge, none of which a well-mixed dFBA represents; see docs/LIMITATIONS.md.

    Returns a :class:`~muode.dfba.SimulationResult`.  ``result.meta["spore_latched"]``
    reports any species whose growth never once reached ``mu_stress`` -- see
    :meth:`muode.lifecycle.SporeForming.latched`.  If it is non-empty the sporulation
    trigger never disengaged and no arm of that run can be read.
    """
    from muode.dfba import DynamicFBA
    from muode.inject import Injection
    from muode.kinetics import KineticParameters
    from muode.lifecycle import SporeForming

    community = build_community(
        abundances={PATHOGEN: 1.0, **{d: 0.0 for d in DONORS}})
    injections = None
    if fmt:
        injections = [Injection.from_abundances(
            fmt_time, {d: 1.0 / len(DONORS) for d in DONORS},
            total_biomass=fmt_biomass, name="FMT")]

    ecology = cdi_ecology(ablate, vmax_bai=vmax_bai)
    # The SAME dilution rate must reach the diet: the bile influx is derived from it so
    # the pool holds at the physiological concentration. Passing one and defaulting the
    # other would put the germinant at the wrong steady state, silently.
    result = DynamicFBA(t_end=t_end, dt=dt, dilution_rate=dilution_rate).run(
        community, cdi_diet(dilution_rate=dilution_rate), KineticParameters(),
        injections=injections, ecology=ecology)

    # Surface the latch check on the result: a silently latched trigger looks exactly
    # like success (the pathogen vanishes), so it must travel with the numbers.
    for layer in ecology.layers:
        if isinstance(layer, SporeForming):
            result.meta["spore_latched"] = layer.latched()
    return result


if __name__ == "__main__":
    print("FMT genome-scale community -- reconstruction readiness:\n")
    for slug, status in readiness().items():
        mark = "OK " if status == "ready" else "-- "
        print(f"  {mark} {slug:32} {status}")

    diet = load_modelseed_diet()
    print(f"\nModelSEED western_gut: {len(diet.metabolites())} metabolites "
          f"({sum(1 for m in diet.metabolites() if diet.initial_concentration(m) > 0)} supplied)\n")
    print("Per-member verification (growth is WITH flux bounds -- the yields that used "
          "to be dials):\n")
    for m in MEMBERS:
        if readiness()[m.slug] != "ready":
            continue
        r = verify(m, diet)
        sec = ("" if r.secretes is None
               else f"  {r.secretes} pathway: {'CONNECTED' if r.pathway_present else 'ABSENT'}"
                    f" (max export {r.can_secrete:.3g})")
        print(f"  {r.slug:32} growth={r.growth:.4f}/h{sec}")
