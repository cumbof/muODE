"""Tests for MuMaxCapped -- the post-solve max-growth-rate ceiling wrapper.

Asserts the mechanism (cap the growth rate and scale the exchange fluxes with it),
the pass-through when below the cap / infeasible, delegation of the other organism
methods, and that in a real dFBA run a wrapped fast grower never exceeds the cap.
"""

import numpy as np

from muode import DynamicFBA, MuMaxCapped
from muode.community import Community
from muode.examples import build_glucose_specialist, build_acetate_specialist, toy_diet, toy_kinetics
from muode.organism import OrganismSolution


class _FakeOrg:
    """Minimal OrganismModel stand-in with a fixed solution, to unit-test the wrapper."""

    def __init__(self, growth, fluxes, status="optimal"):
        self.id = "X"
        self._sol = OrganismSolution(growth, dict(fluxes), status)
        self.reset_calls = 0

    def optimize(self):
        return OrganismSolution(self._sol.growth_rate, dict(self._sol.exchange_fluxes), self._sol.status)

    def copy(self):
        return _FakeOrg(self._sol.growth_rate, self._sol.exchange_fluxes, self._sol.status)

    def reset_bounds(self):
        self.reset_calls += 1

    def exchange_metabolites(self):
        return ("glc_e", "ac_e")


def test_caps_growth_and_scales_fluxes():
    inner = _FakeOrg(1.0, {"glc_e": -10.0, "ac_e": 8.0})
    capped = MuMaxCapped(inner, mumax=0.4)
    sol = capped.optimize()
    assert sol.growth_rate == 0.4                     # capped
    # fluxes scaled by mumax/growth = 0.4
    assert np.isclose(sol.exchange_fluxes["glc_e"], -4.0)
    assert np.isclose(sol.exchange_fluxes["ac_e"], 3.2)


def test_no_change_below_cap():
    inner = _FakeOrg(0.2, {"glc_e": -2.0})
    sol = MuMaxCapped(inner, mumax=0.4).optimize()
    assert sol.growth_rate == 0.2
    assert np.isclose(sol.exchange_fluxes["glc_e"], -2.0)


def test_nonpositive_mumax_disables_cap():
    inner = _FakeOrg(1.5, {"glc_e": -5.0})
    sol = MuMaxCapped(inner, mumax=0.0).optimize()
    assert sol.growth_rate == 1.5                     # uncapped


def test_infeasible_passes_through():
    inner = _FakeOrg(0.0, {}, status="infeasible")
    sol = MuMaxCapped(inner, mumax=0.4).optimize()
    assert not sol.feasible and sol.growth_rate == 0.0


def test_delegates_other_methods_and_id_and_copy():
    inner = _FakeOrg(1.0, {"glc_e": -1.0})
    capped = MuMaxCapped(inner, mumax=0.4)
    assert capped.id == "X"
    assert capped.exchange_metabolites() == ("glc_e", "ac_e")   # delegated
    capped.reset_bounds()                                        # delegated
    assert inner.reset_calls == 1
    dup = capped.copy()
    assert isinstance(dup, MuMaxCapped) and dup._mumax == 0.4


def test_cap_holds_in_a_dfba_run():
    mumax = 0.05
    comm = Community(
        organisms=[MuMaxCapped(build_glucose_specialist(), mumax), build_acetate_specialist()],
        abundances={"A_glucose": 0.7, "B_acetate": 0.3},
        total_biomass=0.02,
    )
    res = DynamicFBA(t_end=12.0, dt=0.1).run(comm, toy_diet(), toy_kinetics())
    # the wrapped grower's realised growth rate never exceeds the ceiling
    assert res.growth_rates["A_glucose"].max() <= mumax + 1e-9
