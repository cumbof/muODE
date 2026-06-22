"""Integration tests for the dynamic-FBA engine on the toy cross-feeding system.

These assert the *biology* the engine is supposed to reproduce, not just that
code runs: a glucose specialist must feed an acetate specialist, and removing
the feeder must drive the dependent species to a secondary extinction.
"""

import numpy as np
import pytest

from muode import DynamicFBA, Perturbation
from muode.examples import build_toy_community, toy_diet, toy_kinetics


@pytest.fixture
def sim():
    return DynamicFBA(t_end=24.0, dt=0.05)


def test_single_organism_fba_respects_uptake_bound():
    from muode.examples import build_glucose_specialist

    org = build_glucose_specialist()
    org.reset_bounds()
    org.set_uptake_bound("glc_e", 4.0)
    sol = org.optimize()
    assert sol.feasible
    # growth = yield (0.10) * glucose throughput, throughput capped at 4.0
    assert sol.growth_rate == pytest.approx(0.10 * 4.0, rel=1e-6)
    # 1:1 acetate secretion
    assert sol.exchange_fluxes["ac_e"] == pytest.approx(4.0, rel=1e-6)
    assert sol.exchange_fluxes["glc_e"] == pytest.approx(-4.0, rel=1e-6)


def test_cross_feeding_supports_growth(sim):
    community = build_toy_community()
    result = sim.run(community, toy_diet(), toy_kinetics())

    init = result.biomass.iloc[0]
    final = result.biomass.iloc[-1]

    # A grows on glucose
    assert final["A_glucose"] > init["A_glucose"]
    # B grows even though it cannot touch glucose -- only via cross-fed acetate
    assert final["B_acetate"] > init["B_acetate"] * 1.5

    # glucose is consumed
    assert result.metabolites["glc_e"].iloc[-1] < result.metabolites["glc_e"].iloc[0]
    # acetate is produced at some point (cross-feeding currency)
    assert result.metabolites["ac_e"].max() > 1e-3


def test_b_blooms_after_a(sim):
    """B's growth should lag A's: it needs acetate to exist first."""
    result = sim.run(build_toy_community(), toy_diet(), toy_kinetics())
    mu = result.growth_rates
    # first time each species reaches a small positive growth rate
    a_on = np.argmax(mu["A_glucose"].values > 1e-4)
    b_on = np.argmax(mu["B_acetate"].values > 1e-4)
    assert b_on >= a_on


def test_cross_feeding_inference(sim):
    result = sim.run(build_toy_community(), toy_diet(), toy_kinetics())
    cf = result.cross_feeding()
    pairs = set(zip(cf["producer"], cf["metabolite"], cf["consumer"]))
    assert ("A_glucose", "ac_e", "B_acetate") in pairs


def test_metabolites_never_negative(sim):
    result = sim.run(build_toy_community(), toy_diet(), toy_kinetics())
    assert (result.metabolites.values >= -1e-9).all()


def test_secondary_extinction_on_feeder_removal(sim):
    """Removing A must starve B even though the perturbation never touches B."""
    community = build_toy_community()
    pert = Perturbation.remove_species(["A_glucose"])
    result = sim.run(community, toy_diet(), toy_kinetics(), perturbation=pert)

    init = result.biomass.iloc[0]
    final = result.biomass.iloc[-1]

    # A produces no biomass
    assert final["A_glucose"] <= init["A_glucose"] + 1e-9
    # B fails to grow without cross-fed acetate
    assert final["B_acetate"] <= init["B_acetate"] * 1.01
    # and essentially no acetate is ever produced
    assert result.metabolites["ac_e"].max() < 1e-6
    assert "B_acetate" in result.extinct(abs_threshold=0.0, rel_threshold=1.5)


def test_reproducible(sim):
    a = sim.run(build_toy_community(), toy_diet(), toy_kinetics())
    b = sim.run(build_toy_community(), toy_diet(), toy_kinetics())
    assert np.allclose(a.biomass.values, b.biomass.values)
