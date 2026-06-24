"""Tests for timed biomass injection (transplants / probiotic doses).

These assert the *mechanism* an FMT-style simulation relies on: an organism that
is absent (zero biomass) until an injection time, then introduced as a bolus and
allowed to grow on the shared metabolite pool -- and that injecting a donor into
an existing community is correctly assembled by ``merge_for_injection``.
"""

import numpy as np
import pytest

from muode import DynamicFBA, Injection, merge_for_injection
from muode.community import Community
from muode.examples import (
    build_acetate_specialist,
    build_glucose_specialist,
    toy_diet,
    toy_kinetics,
)


def test_from_abundances_normalises_to_total():
    inj = Injection.from_abundances(12.0, {"a": 3.0, "b": 1.0}, total_biomass=0.08, name="FMT")
    assert inj.time == 12.0
    assert inj.name == "FMT"
    assert inj.biomass["a"] == pytest.approx(0.06)
    assert inj.biomass["b"] == pytest.approx(0.02)
    assert inj.total() == pytest.approx(0.08)


def test_merge_for_injection_keeps_residents_and_seeds_donors():
    base = Community([build_glucose_specialist()], {"A_glucose": 1.0}, total_biomass=0.02)
    donor = Community([build_acetate_specialist()], {"B_acetate": 1.0}, total_biomass=0.01)

    merged, event = merge_for_injection(base, donor, time=10.0, name="FMT")

    assert set(merged.organism_ids) == {"A_glucose", "B_acetate"}
    # resident keeps its initial biomass; donor is dormant at t=0
    assert merged.initial_biomass()["A_glucose"] == pytest.approx(0.02)
    assert merged.initial_biomass()["B_acetate"] == 0.0
    # the event introduces the donor's biomass at the chosen time
    assert event.time == 10.0
    assert event.biomass["B_acetate"] == pytest.approx(0.01)


def test_injected_species_dormant_then_blooms():
    """B is absent until t=10 h, then injected and grows on accumulated acetate."""
    base = Community([build_glucose_specialist()], {"A_glucose": 1.0}, total_biomass=0.02)
    donor = Community([build_acetate_specialist()], {"B_acetate": 1.0}, total_biomass=0.01)
    merged, event = merge_for_injection(base, donor, time=10.0)

    result = DynamicFBA(t_end=24.0, dt=0.05).run(
        merged, toy_diet(), toy_kinetics(), injections=[event]
    )
    bio = result.biomass

    # B contributes nothing before the injection
    pre = bio["B_acetate"][bio.index < 10.0 - 1e-9]
    assert (pre.values == 0.0).all()

    # the bolus appears exactly at the injection time
    at_inject = bio["B_acetate"].loc[bio.index[np.argmin(np.abs(bio.index - 10.0))]]
    assert at_inject == pytest.approx(0.01, abs=1e-6)

    # and B grows afterwards on the acetate A accumulated while B was absent
    assert bio["B_acetate"].iloc[-1] > 0.01

    # the standing acetate pool that built up while B was absent is drawn down
    # after the injection -- B feeds on A's leftover fermentation product
    ac = result.metabolites["ac_e"]
    ac_at_inject = ac.loc[ac.index[np.argmin(np.abs(ac.index - 10.0))]]
    assert ac_at_inject > 1e-3
    assert ac.iloc[-1] < ac_at_inject


def test_injection_boosts_existing_species_without_duplicating():
    base = Community([build_glucose_specialist()], {"A_glucose": 1.0}, total_biomass=0.02)
    donor = Community([build_glucose_specialist()], {"A_glucose": 1.0}, total_biomass=0.01)

    merged, event = merge_for_injection(base, donor, time=5.0)

    # no duplicate organism id
    assert merged.organism_ids == ["A_glucose"]
    # the injection boosts the resident strain
    assert event.biomass["A_glucose"] == pytest.approx(0.01)

    # at the injection step the biomass jumps by the bolus on top of growth
    result = DynamicFBA(t_end=10.0, dt=0.05).run(
        merged, toy_diet(), toy_kinetics(), injections=[event]
    )
    bio = result.biomass["A_glucose"]
    idx = np.argmin(np.abs(bio.index - 5.0))
    jump = bio.iloc[idx] - bio.iloc[idx - 1]
    # the step-to-step change includes the +0.01 bolus, so it exceeds plain growth
    assert jump > 0.01


def test_injection_recorded_in_meta():
    base = Community([build_glucose_specialist()], {"A_glucose": 1.0}, total_biomass=0.02)
    donor = Community([build_acetate_specialist()], {"B_acetate": 1.0}, total_biomass=0.01)
    merged, event = merge_for_injection(base, donor, time=8.0, name="FMT")
    result = DynamicFBA(t_end=12.0, dt=0.1).run(
        merged, toy_diet(), toy_kinetics(), injections=[event]
    )
    assert result.meta["injections"]
    assert "FMT" in result.meta["injections"][0]


def test_no_injection_matches_baseline():
    """injections=None must not change a run (regression guard)."""
    base = Community(
        [build_glucose_specialist(), build_acetate_specialist()],
        {"A_glucose": 0.7, "B_acetate": 0.3},
        total_biomass=0.02,
    )
    a = DynamicFBA(t_end=12.0, dt=0.1).run(base, toy_diet(), toy_kinetics())
    b = DynamicFBA(t_end=12.0, dt=0.1).run(base, toy_diet(), toy_kinetics(), injections=None)
    assert np.allclose(a.biomass.values, b.biomass.values)
    assert a.meta["injections"] is None
