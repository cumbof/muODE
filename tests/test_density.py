"""Tests for the logistic carrying-capacity ecology layer.

The layer adds Verhulst density-dependence: growth slows to zero as an organism's
biomass nears its carrying capacity K.  These check the factor math, that a single
grower saturates near K, and -- the reason the layer exists -- that per-organism
capacities lift the excluded minority of a two-species community without a shared
uniform slowdown (which would cancel in relative abundance).
"""

import numpy as np
import pytest

from muode import DynamicFBA, EcologyModel, LogisticCarryingCapacity
from muode.community import Community
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


def _diet(glc=50.0):
    from muode.diet import Diet
    return Diet(name="glc", concentrations={"glc_e": glc})


# ---------------------------------------------------------------------------
# factor math
# ---------------------------------------------------------------------------
def test_growth_factor_is_logistic():
    lay = LogisticCarryingCapacity(capacity={"a": 2.0})
    assert lay.growth_factor("a", 0.0, {}, {"a": 0.0}) == pytest.approx(1.0)
    assert lay.growth_factor("a", 0.0, {}, {"a": 1.0}) == pytest.approx(0.5)
    assert lay.growth_factor("a", 0.0, {}, {"a": 2.0}) == pytest.approx(0.0)
    # floor clamps: above K stays at floor (default 0.0), not negative
    assert lay.growth_factor("a", 0.0, {}, {"a": 3.0}) == pytest.approx(0.0)


def test_unknown_organism_is_unconstrained():
    lay = LogisticCarryingCapacity(capacity={"a": 2.0})
    assert lay.growth_factor("b", 0.0, {}, {"b": 99.0}) == pytest.approx(1.0)
    # ...unless a default capacity is given
    lay2 = LogisticCarryingCapacity(capacity={}, default_capacity=1.0)
    assert lay2.growth_factor("b", 0.0, {}, {"b": 1.0}) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# a single grower saturates near its carrying capacity
# ---------------------------------------------------------------------------
def test_single_grower_saturates_at_capacity():
    K = 0.3
    comm = Community([_glucose_grower("a")], {"a": 1.0}, 0.01)
    eco = EcologyModel([LogisticCarryingCapacity(capacity={"a": K})])
    res = DynamicFBA(t_end=48.0, dt=0.1).run(comm, _diet(), _glucose_kin(), ecology=eco)
    final = float(res.biomass.iloc[-1]["a"])
    # logistic arrests growth at K: endpoint should sit just below K, never above
    assert final <= K + 1e-6
    assert final > 0.5 * K   # and it did climb toward it


# ---------------------------------------------------------------------------
# the point of the layer: lift the excluded minority
# ---------------------------------------------------------------------------
def test_capacity_lifts_the_minority():
    fast, slow = _glucose_grower("fast", yld=0.20), _glucose_grower("slow", yld=0.10)
    kin, diet = _glucose_kin(), _diet()

    def minority(eco):
        comm = Community([_glucose_grower("fast", 0.20), _glucose_grower("slow", 0.10)],
                         {"fast": 0.5, "slow": 0.5}, 0.01)
        res = DynamicFBA(t_end=48.0, dt=0.1).run(comm, diet, kin, ecology=eco)
        b = res.biomass.iloc[-1].to_dict()
        tot = b["fast"] + b["slow"]
        return b["slow"] / tot if tot > 0 else 0.0

    base = minority(None)
    # equal per-capita capacities let the slower grower hold a larger share
    capped = minority(EcologyModel([LogisticCarryingCapacity(
        capacity={"fast": 0.3, "slow": 0.3})]))
    assert capped > base + 0.05
