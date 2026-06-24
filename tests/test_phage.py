"""Bacteriophage predation as a coupled infection ODE."""

import pytest

from muode.community import Community
from muode.dfba import DynamicFBA
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism
from muode.phage import PhageInfection
from muode.scenarios import COMMENSAL, PATHOGEN, phage_predation_scenario


def _grower(id, yld=0.2):
    return LinprogOrganism(
        id=id,
        reactions=[("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0),
                   ("GROW", {"glc_e": -1.0}, 0.0, 1000.0)],
        objective={"GROW": yld},
        exchanges={"glc_e": "EX_glc_e"},
    )


def test_phage_amplifies_and_removes_host():
    p = PhageInfection(host="h", adsorption_rate=10.0, burst_size=50.0,
                       latent_period=0.4, decay_rate=0.1, initial_titer=0.5)
    p.reset(None)
    X = {"h": 1.0}
    for _ in range(50):
        p.integrate(0.0, 0.05, X, {}, {"h": 0.0})
    # host has been driven down and phage has amplified above the inoculum
    assert X["h"] < 1.0
    assert p.titer > 0.5


def test_no_phage_means_no_effect():
    p = PhageInfection(host="h", initial_titer=0.0)
    p.reset(None)
    X = {"h": 1.0}
    for _ in range(50):
        p.integrate(0.0, 0.05, X, {}, {"h": 0.0})
    assert X["h"] == pytest.approx(1.0)
    assert p.titer == 0.0


def test_lysogeny_blocks_predation():
    """A fully temperate phage (lysogeny=1) lyses nothing and does not amplify."""
    p = PhageInfection(host="h", adsorption_rate=10.0, burst_size=50.0,
                       initial_titer=1.0, lysogeny_fraction=1.0)
    p.reset(None)
    X = {"h": 1.0}
    for _ in range(50):
        p.integrate(0.0, 0.05, X, {}, {"h": 0.0})
    # host is not lytically removed and no burst occurs
    assert X["h"] == pytest.approx(1.0)
    assert p.infected_biomass == pytest.approx(0.0)
    assert p.titer <= 1.0  # only decays, never grows


def test_cfl_cap_never_oversubscribes_host():
    p = PhageInfection(host="h", adsorption_rate=1e6, burst_size=10.0,
                       initial_titer=1e6, latent_period=0.4)
    p.reset(None)
    X = {"h": 1.0}
    p.integrate(0.0, 0.05, X, {}, {"h": 0.0})
    assert X["h"] >= 0.0  # never goes negative despite an enormous adsorption term


def test_missing_host_is_a_noop():
    p = PhageInfection(host="absent", initial_titer=1.0)
    p.reset(None)
    X = {"present": 1.0}
    p.integrate(0.0, 0.05, X, {}, {})  # must not raise or touch other species
    assert X["present"] == 1.0


def test_phage_layer_runs_in_engine_and_records():
    comm = Community([_grower("h")], {"h": 1.0}, total_biomass=0.05)
    diet = Diet(concentrations={"glc_e": 8.0}, influx={"glc_e": 1.0}, name="c")
    kin = KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})
    eco = EcologyModel([PhageInfection(host="h", name="vir", adsorption_rate=12.0,
                                       burst_size=60.0, initial_titer=0.5)])
    res = DynamicFBA(t_end=24, dt=0.02).run(comm, diet, kin, ecology=eco)
    assert res.environment is not None
    assert "phage[vir]" in res.environment.columns
    assert res.environment["phage[vir]"].max() > 0.5


def test_phage_therapy_scenario_controls_bloom():
    no = phage_predation_scenario(therapy=False)
    yes = phage_predation_scenario(therapy=True)
    # phage knocks the pathogen well below its phage-free level ...
    assert yes.final_biomass()[PATHOGEN] < 0.25 * no.final_biomass()[PATHOGEN]
    # ... and releases the commensal it was out-competing
    assert yes.final_biomass()[COMMENSAL] > no.final_biomass()[COMMENSAL]
