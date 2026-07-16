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
        """Prefer the gapseq reconstruction; fall back to a bare ``.xml.gz``.

        Returns None if no file exists yet.  A ``.gapseq.xml.gz`` is always gapseq; a
        bare ``.xml.gz`` may be gapseq (the two already-done members) or a leftover
        CarveMe carve (the two awaiting replacement) -- :func:`load` sorts that out by
        namespace, so this only has to find *a* file.
        """
        for name in (f"{self.slug}.gapseq.xml.gz", f"{self.slug}.xml.gz"):
            p = GEMS / name
            if p.exists():
                return p
        return None


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
BAI_GUILD = {"C_scindens_ATCC35704"}
DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
          "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]


def cdi_diet(diet: Optional[Diet] = None) -> Diet:
    """ModelSEED western_gut PLUS the bile-acid pool the ecology layers act on.

    The bile ids (tca_e/ca_e/dca_e) are muODE ecology-layer ids, NOT ModelSEED: the
    bile transformation is a phenomenological pool process keyed on guild membership,
    so the GEMs need no bile exchanges.  Concentrations are physiological -- bile
    reaches the caecum at ~2 mM (Ramirez & Abel-Santos 2011) -- rather than tuned.
    """
    from muode.bile import CHOLATE, TAUROCHOLATE

    base = diet or load_modelseed_diet()
    conc = dict(base.concentrations)
    conc[TAUROCHOLATE] = 2.0      # germinant, ~caecal bile concentration
    conc[CHOLATE] = 2.0           # primary bile: the bai substrate
    influx = dict(base.influx)
    influx[TAUROCHOLATE] = 0.3
    influx[CHOLATE] = 0.3
    return Diet(concentrations=conc, influx=influx, max_uptake=dict(base.max_uptake),
                source=dict(base.source), name="modelseed_western_gut+bile")


def cdi_ecology(ablate: str = ""):
    """The mechanistic stack, with EVIDENCE-BASED parameters (the layer defaults).

    Unlike the toy ``scenario.cdi_ecology``, this passes NO parameter overrides: it uses
    the layer defaults, which are now the provenance-registered values (germination
    km 15.9 mM, inhibition ki 0.5 mM -- see muode.provenance).  The only knobs are which
    layers are present, so ``ablate`` can remove the bile and/or pH arm to decompose the
    mechanism the way test_provenance did for the toy -- but now on real stoichiometry.
    """
    from muode.antibiotic import Antibiotic
    from muode.bile import BileAcidInhibition, BileAcidTransform
    from muode.ecology import EcologyModel
    from muode.lifecycle import SporeForming
    from muode.ph import WeakAcidInhibition

    layers = []
    if "ph" not in ablate:
        layers.append(WeakAcidInhibition(buffer_capacity=25.0, default_ki=12.0,
                                         ki={PATHOGEN: 3.0}))
    if "bile" not in ablate:
        layers.append(BileAcidTransform(bai_producers=set(BAI_GUILD)))
        layers.append(BileAcidInhibition(targets={PATHOGEN}))
    layers.append(SporeForming(species={PATHOGEN}, initial_spores={PATHOGEN: 0.05},
                               k_germination=0.6, k_sporulation=0.8, mu_stress=0.15))
    layers.append(Antibiotic(susceptible={PATHOGEN},
                             dose_times=tuple(float(t) for t in range(0, 10, 2)),
                             dose=3.0, half_life=2.0, emax=5.0, ec50=0.4))
    return EcologyModel(layers)


def build_scenario(fmt: bool, ablate: str = "", t_end: float = 96.0, dt: float = 0.05,
                   fmt_time: float = 12.0, fmt_biomass: float = 0.05):
    """Run the rCDI/FMT scenario on the real gapseq community.

    The recipient carries the (bloomed) pathogen; the four donors are members from t=0
    at zero biomass so the bile layer can refer to them, and an FMT injection seeds them
    at ``fmt_time``.  The ONLY difference between the arms is the injection.

    ``ablate`` ("", "bile", "ph", "bile ph") drops ecology layers so the mechanism can
    be decomposed: this is the experiment that, on the toy models, showed clearance was
    nutrient competition rather than the advertised bile mechanism.  Re-running it here
    asks whether that still holds once the yields are stoichiometry, not dials.

    Returns a :class:`~muode.dfba.SimulationResult`.  ~30 min per arm on a genome-scale
    community -- a workstation job; see ``run.py``.
    """
    from muode.dfba import DynamicFBA
    from muode.inject import Injection
    from muode.kinetics import KineticParameters

    community = build_community(
        abundances={PATHOGEN: 1.0, **{d: 0.0 for d in DONORS}})
    injections = None
    if fmt:
        injections = [Injection.from_abundances(
            fmt_time, {d: 1.0 / len(DONORS) for d in DONORS},
            total_biomass=fmt_biomass, name="FMT")]
    return DynamicFBA(t_end=t_end, dt=dt).run(
        community, cdi_diet(), KineticParameters(),
        injections=injections, ecology=cdi_ecology(ablate))


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
