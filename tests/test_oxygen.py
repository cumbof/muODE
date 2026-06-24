"""Oxygen tolerance gating and scavenging (the cross-kingdom mechanism)."""

import pytest

from muode.community import Community
from muode.dfba import DynamicFBA
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism
from muode.oxygen import OXYGEN, OxygenSensitivity
from muode.traits import Domain, MicrobeTraits, OxygenTolerance


def _grower(id, yld=0.2):
    return LinprogOrganism(
        id=id,
        reactions=[("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0),
                   ("GROW", {"glc_e": -1.0}, 0.0, 1000.0)],
        objective={"GROW": yld},
        exchanges={"glc_e": "EX_glc_e"},
    )


def test_growth_factor_by_tolerance():
    layer = OxygenSensitivity(obligate_aerobes={"aer"}, obligate_anaerobes={"ana"},
                              facultative={"fac"}, km_o2=0.01, ki_o2=0.02)
    hi = {OXYGEN: 0.25}
    lo = {OXYGEN: 0.0}
    # aerobe: needs O2
    assert layer.growth_factor("aer", 0, hi, {}) > 0.9
    assert layer.growth_factor("aer", 0, lo, {}) == 0.0
    # anaerobe: poisoned by O2, fine without it
    assert layer.growth_factor("ana", 0, hi, {}) < 0.1
    assert layer.growth_factor("ana", 0, lo, {}) == pytest.approx(1.0)
    # facultative and unlisted: indifferent
    assert layer.growth_factor("fac", 0, hi, {}) == 1.0
    assert layer.growth_factor("other", 0, hi, {}) == 1.0


def test_scavengers_consume_oxygen():
    layer = OxygenSensitivity(facultative={"fac"}, obligate_anaerobes={"ana"},
                              consumption=8.0)
    rates = layer.metabolite_rates(0, {OXYGEN: 0.25}, {"fac": 1.0, "ana": 5.0})
    assert rates[OXYGEN] < 0          # facultative draws O2 down
    # anaerobes do not scavenge
    only_ana = layer.metabolite_rates(0, {OXYGEN: 0.25}, {"ana": 5.0})
    assert only_ana == {}
    # with consume=False the layer is pure gating, no O2 sink
    pure = OxygenSensitivity(facultative={"fac"}, consume=False)
    assert pure.metabolite_rates(0, {OXYGEN: 0.25}, {"fac": 1.0}) == {}


def test_from_traits_builds_sets():
    traits = {
        "yeast": MicrobeTraits(domain=Domain.EUKARYOTE, oxygen=OxygenTolerance.FACULTATIVE),
        "bact": MicrobeTraits(oxygen=OxygenTolerance.OBLIGATE_ANAEROBE),
        "myco": MicrobeTraits(oxygen=OxygenTolerance.OBLIGATE_AEROBE),
        "phage": MicrobeTraits(domain=Domain.VIRUS),
    }
    layer = OxygenSensitivity.from_traits(traits)
    assert layer.facultative == {"yeast"}
    assert layer.obligate_anaerobes == {"bact"}
    assert layer.obligate_aerobes == {"myco"}
    # a virus is not a metabolic member and is skipped entirely
    assert "phage" not in (layer.facultative | layer.obligate_anaerobes | layer.obligate_aerobes)


def test_facultative_scavenger_rescues_anaerobe():
    """A facultative O2 consumer should let an obligate anaerobe grow."""
    diet = Diet(concentrations={"glc_e": 10.0}, influx={"glc_e": 1.0}, name="ox")
    kin = KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})

    anaerobe = _grower("ana")
    facultative = _grower("fac")

    # with a facultative scavenger that knocks O2 down, the anaerobe grows
    comm_pair = Community([anaerobe, _grower("fac")], {"ana": 0.5, "fac": 0.5},
                          total_biomass=0.02)
    eco_pair = EcologyModel([OxygenSensitivity(
        obligate_anaerobes={"ana"}, facultative={"fac"},
        o2_initial=0.25, o2_influx=0.05, consumption=20.0)])
    res_pair = DynamicFBA(t_end=24, dt=0.05).run(comm_pair, diet, kin, ecology=eco_pair)

    # alone, with constant O2 influx and no scavenger, the anaerobe is suppressed
    comm_solo = Community([_grower("ana")], {"ana": 1.0}, total_biomass=0.02)
    eco_solo = EcologyModel([OxygenSensitivity(
        obligate_anaerobes={"ana"}, o2_initial=0.25, o2_influx=0.05, consumption=20.0)])
    res_solo = DynamicFBA(t_end=24, dt=0.05).run(comm_solo, diet, kin, ecology=eco_solo)

    assert res_pair.final_biomass()["ana"] > 2.0 * res_solo.final_biomass()["ana"]
    # the scavenger drove oxygen down
    assert res_pair.environment["oxygen"].iloc[-1] < res_solo.environment["oxygen"].iloc[-1]


def test_obligate_aerobe_needs_oxygen():
    diet = Diet(concentrations={"glc_e": 10.0}, influx={"glc_e": 1.0}, name="ox")
    kin = KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})
    comm = Community([_grower("aer")], {"aer": 1.0}, total_biomass=0.02)

    with_o2 = EcologyModel([OxygenSensitivity(obligate_aerobes={"aer"},
                                              o2_initial=0.5, o2_influx=0.5, consume=False)])
    no_o2 = EcologyModel([OxygenSensitivity(obligate_aerobes={"aer"},
                                            o2_initial=0.0, o2_influx=0.0)])
    grew = DynamicFBA(t_end=12, dt=0.05).run(comm, diet, kin, ecology=with_o2)
    starved = DynamicFBA(t_end=12, dt=0.05).run(comm, diet, kin, ecology=no_o2)
    assert grew.final_biomass()["aer"] > 5.0 * starved.final_biomass()["aer"]
