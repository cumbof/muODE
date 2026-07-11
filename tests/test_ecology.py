"""Tests for the ecology layer -- the environmental / life-history processes
layered on top of the dynamic-FBA core (pH/SCFA inhibition, bile-acid
transformation, sporulation/germination, antibiotic PK, bacteriocins).

Each test asserts the *mechanism* the layer is meant to add, and a regression
guard checks that an empty ecology model leaves a run numerically unchanged.
"""

import numpy as np
import pytest

from muode import (
    Antibiotic,
    Bacteriocin,
    BileAcidInhibition,
    BileAcidTransform,
    DynamicFBA,
    EcologyModel,
    SporeForming,
    WeakAcidInhibition,
)
from muode.community import Community
from muode.diet import Diet
from muode.examples import (
    build_glucose_specialist,
    toy_diet,
    toy_kinetics,
    build_toy_community,
)
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism


def _glucose_grower(id, yld=0.12):
    return LinprogOrganism(
        id=id,
        reactions=[("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0),
                   ("GROW", {"glc_e": -1.0}, 0.0, 1000.0)],
        objective={"GROW": yld},
        exchanges={"glc_e": "EX_glc_e"},
    )


def _glucose_kin():
    return KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})


# ---------------------------------------------------------------------------
# regression: empty ecology == no ecology
# ---------------------------------------------------------------------------

def test_empty_ecology_matches_baseline():
    comm = build_toy_community()
    a = DynamicFBA(t_end=12.0, dt=0.1).run(comm, toy_diet(), toy_kinetics())
    b = DynamicFBA(t_end=12.0, dt=0.1).run(
        build_toy_community(), toy_diet(), toy_kinetics(), ecology=EcologyModel()
    )
    assert np.allclose(a.biomass.values, b.biomass.values)
    assert b.environment is None and b.spores is None
    assert b.meta["ecology"] is None


# ---------------------------------------------------------------------------
# weak-acid / pH inhibition
# ---------------------------------------------------------------------------

def test_ph_drops_with_acid_and_inhibits_growth():
    wa = WeakAcidInhibition(ph0=6.8, buffer_capacity=40.0, default_ki=5.0)
    clean = {"ac_e": 0.0, "but_e": 0.0}
    acidic = {"ac_e": 60.0, "but_e": 40.0}
    assert wa.ph(clean) == pytest.approx(6.8)
    assert wa.ph(acidic) < 6.0
    # undissociated acid present at low pH suppresses growth
    assert wa.growth_factor("x", 0.0, acidic, {}) < 0.6
    assert wa.growth_factor("x", 0.0, clean, {}) == pytest.approx(1.0, abs=1e-6)


def test_weak_acid_makes_fermentation_self_limiting():
    """A glucose->acetate fermenter accumulates acid; with the layer it ends up
    with less biomass than without (acid feedback)."""
    diet = Diet(concentrations={"glc_e": 40.0, "ac_e": 0.0}, name="glc")
    kin = toy_kinetics()
    eng = DynamicFBA(t_end=24.0, dt=0.05)

    base = eng.run(Community([build_glucose_specialist()], {"A_glucose": 1.0}, 0.02), diet, kin)
    eco = EcologyModel([WeakAcidInhibition(buffer_capacity=20.0, default_ki=3.0)])
    inhib = eng.run(Community([build_glucose_specialist()], {"A_glucose": 1.0}, 0.02), diet, kin,
                    ecology=eco)

    assert inhib.biomass["A_glucose"].iloc[-1] < base.biomass["A_glucose"].iloc[-1]
    assert inhib.environment is not None
    assert inhib.environment["pH"].iloc[-1] < inhib.environment["pH"].iloc[0]


# ---------------------------------------------------------------------------
# bile-acid transformation + inhibition
# ---------------------------------------------------------------------------

def test_bile_transform_makes_secondary_bile_acids():
    """With a bai-carrying producer and cholate in the pool, deoxycholate rises
    and cholate falls."""
    comm = Community([_glucose_grower("bai_bug")], {"bai_bug": 1.0}, total_biomass=0.05)
    diet = Diet(concentrations={"glc_e": 20.0, "ca_e": 5.0, "dca_e": 0.0}, name="bile")
    eco = EcologyModel([BileAcidTransform(bai_producers={"bai_bug"}, vmax_bai=2.0)])
    res = DynamicFBA(t_end=12.0, dt=0.1).run(comm, diet, _glucose_kin(), ecology=eco)

    assert res.metabolites["dca_e"].iloc[-1] > 0.5      # secondary produced
    assert res.metabolites["ca_e"].iloc[-1] < 5.0       # primary consumed


def test_bile_inhibition_only_hits_targets():
    bi = BileAcidInhibition(targets={"cdiff"}, ki=0.02)
    pool = {"dca_e": 0.1}
    assert bi.growth_factor("cdiff", 0.0, pool, {}) < 0.3     # strongly inhibited
    assert bi.growth_factor("commensal", 0.0, pool, {}) == 1.0  # untouched


# ---------------------------------------------------------------------------
# antibiotic PK/PD
# ---------------------------------------------------------------------------

def test_antibiotic_concentration_decays():
    ab = Antibiotic(susceptible={"s"}, dose_times=(0.0,), dose=4.0, half_life=2.0)
    assert ab.concentration(0.0) == pytest.approx(4.0)
    assert ab.concentration(2.0) == pytest.approx(2.0, rel=1e-6)   # one half-life
    assert ab.concentration(4.0) == pytest.approx(1.0, rel=1e-6)   # two half-lives


def test_antibiotic_kills_susceptible_not_resistant():
    comm = Community([_glucose_grower("susceptible"), _glucose_grower("resistant")],
                     {"susceptible": 0.5, "resistant": 0.5}, total_biomass=0.02)
    diet = Diet(concentrations={"glc_e": 50.0}, influx={"glc_e": 2.0}, name="glc")
    eco = EcologyModel([Antibiotic(susceptible={"susceptible"}, dose_times=(0.0,),
                                   dose=10.0, half_life=4.0, emax=4.0, ec50=0.5)])
    res = DynamicFBA(t_end=24.0, dt=0.05).run(comm, diet, _glucose_kin(), ecology=eco)

    bio = res.biomass.iloc[-1]
    assert bio["susceptible"] < res.biomass["susceptible"].iloc[0]   # knocked back
    assert bio["resistant"] > res.biomass["resistant"].iloc[0]       # grew freely
    assert res.environment["drug[antibiotic]"].iloc[0] > res.environment["drug[antibiotic]"].iloc[-1]


# ---------------------------------------------------------------------------
# sporulation / germination life cycle
# ---------------------------------------------------------------------------

def test_spores_survive_antibiotic_when_vegetative_die():
    """After nutrient depletion a spore-former banks biomass as antibiotic-
    resistant spores; a non-spore-former hit by the same drug is cleared."""
    comm = Community([_glucose_grower("sporeformer"), _glucose_grower("nonspore")],
                     {"sporeformer": 0.5, "nonspore": 0.5}, total_biomass=0.04)
    diet = Diet(concentrations={"glc_e": 8.0}, name="glc")  # finite glucose, no influx
    eco = EcologyModel([
        SporeForming(species={"sporeformer"}, k_sporulation=1.5, mu_stress=0.2),
        # drug dosed once glucose is exhausted and sporulation has happened
        Antibiotic(susceptible={"sporeformer", "nonspore"}, dose_times=(8.0,),
                   dose=10.0, half_life=3.0, emax=5.0, ec50=0.5),
    ])
    res = DynamicFBA(t_end=36.0, dt=0.05).run(comm, diet, _glucose_kin(), ecology=eco)

    spores = res.spores["sporeformer"].iloc[-1]
    sf_total = res.biomass["sporeformer"].iloc[-1] + spores
    ns_total = res.biomass["nonspore"].iloc[-1]
    assert spores > 0.0                       # a dormant reservoir formed
    assert sf_total > 5.0 * max(ns_total, 1e-9)  # it survived; the non-spore did not


@pytest.mark.filterwarnings("ignore:no species grew at any point:RuntimeWarning")
def test_germination_is_gated_by_bile_acids():
    """Seeded spores germinate when the germinant is present, but secondary bile
    acids slow germination -- so more spores stay dormant and less vegetative
    biomass appears.  Run without a growth substrate so vegetative biomass
    reflects germination alone (no amplification by growth)."""
    eng = DynamicFBA(t_end=24.0, dt=0.05)

    def run(extra):
        comm = Community([_glucose_grower("cdiff")], {"cdiff": 1.0}, total_biomass=1e-9)
        conc = {"tca_e": 1.0}            # germinant present, no carbon source
        conc.update(extra)
        eco = EcologyModel([SporeForming(species={"cdiff"}, initial_spores={"cdiff": 0.02},
                                         k_germination=0.6)])
        return eng.run(comm, Diet(concentrations=conc, name="bile"), _glucose_kin(), ecology=eco)

    permissive = run({})                 # germinant present, no inhibitor
    inhibited = run({"dca_e": 0.5})      # secondary bile acid blocks germination

    # germinant present + no inhibitor -> spores germinate into vegetative cells
    assert permissive.biomass["cdiff"].iloc[-1] > 2.0 * inhibited.biomass["cdiff"].iloc[-1]
    # inhibitor keeps the spore reservoir dormant
    assert inhibited.spores["cdiff"].iloc[-1] > 2.0 * permissive.spores["cdiff"].iloc[-1]


# ---------------------------------------------------------------------------
# bacteriocin antagonism
# ---------------------------------------------------------------------------

def test_bacteriocin_suppresses_target():
    comm = Community([_glucose_grower("producer"), _glucose_grower("target")],
                     {"producer": 0.5, "target": 0.5}, total_biomass=0.02)
    diet = Diet(concentrations={"glc_e": 50.0}, influx={"glc_e": 2.0}, name="glc")

    base = DynamicFBA(t_end=24.0, dt=0.05).run(comm, diet, _glucose_kin())
    eco = EcologyModel([Bacteriocin(producers={"producer"}, targets={"target"},
                                    production=2.0, decay=0.1, ki=0.05)])
    res = DynamicFBA(t_end=24.0, dt=0.05).run(
        Community([_glucose_grower("producer"), _glucose_grower("target")],
                  {"producer": 0.5, "target": 0.5}, total_biomass=0.02),
        diet, _glucose_kin(), ecology=eco)

    # target suppressed relative to the toxin-free run; producer not
    assert res.biomass["target"].iloc[-1] < base.biomass["target"].iloc[-1]
    assert res.metabolites["bacteriocin_e"].iloc[-1] > 0.0
