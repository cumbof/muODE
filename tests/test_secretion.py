"""Tests for GrowthCoupledSecretion -- growth-proportional metabolite release.

Asserts the mechanism (a species adds yield * biomass-increment to the pool each
step), that only positive increments secrete, that declared metabolites are tracked,
and that in a real dFBA run the product accumulates in proportion to the producer's
biomass gain.
"""

import numpy as np

from muode import DynamicFBA, EcologyModel, GrowthCoupledSecretion
from muode.examples import build_toy_community, toy_diet, toy_kinetics


def test_extra_metabolites_declares_products():
    layer = GrowthCoupledSecretion(yields={"A_glucose": {"sig_e": 5.0}, "B_acetate": {"sig_e": 1.0}})
    assert set(layer.extra_metabolites()) == {"sig_e"}


def test_metabolite_rate_is_yield_times_increment_over_dt():
    comm = build_toy_community()
    layer = GrowthCoupledSecretion(yields={"A_glucose": {"sig_e": 5.0}}, dt=0.1)
    layer.reset(comm)
    # first call: dX = 1.0 - 0.0 (reset baseline) -> rate = 5.0 * 1.0 / 0.1
    rates = layer.metabolite_rates(0.0, {}, {"A_glucose": 1.0, "B_acetate": 0.0})
    assert np.isclose(rates["sig_e"], 50.0)
    # second call, no biomass change -> no secretion
    rates2 = layer.metabolite_rates(0.1, {}, {"A_glucose": 1.0, "B_acetate": 0.0})
    assert rates2.get("sig_e", 0.0) == 0.0


def test_shrinking_population_does_not_secrete():
    comm = build_toy_community()
    layer = GrowthCoupledSecretion(yields={"A_glucose": {"sig_e": 5.0}}, dt=0.1)
    layer.reset(comm)
    layer.metabolite_rates(0.0, {}, {"A_glucose": 2.0})            # grow to 2.0
    rates = layer.metabolite_rates(0.1, {}, {"A_glucose": 1.0})    # shrink to 1.0
    assert rates.get("sig_e", 0.0) == 0.0


def test_product_accumulates_in_a_dfba_run():
    layer = GrowthCoupledSecretion(yields={"A_glucose": {"sig_e": 5.0}}, dt=0.1)
    res = DynamicFBA(t_end=12.0, dt=0.1).run(
        build_toy_community(), toy_diet(), toy_kinetics(), ecology=EcologyModel([layer])
    )
    sig_final = res.metabolites["sig_e"].iloc[-1]
    a_final = res.biomass["A_glucose"].iloc[-1]
    assert sig_final > 0.0
    # telescoping sum of positive increments from the reset baseline (0) == final biomass,
    # so the accumulated product is yield * final producer biomass.
    assert np.isclose(sig_final, 5.0 * a_final, rtol=0.05)
