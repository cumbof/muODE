"""Solver backends: the default path is unchanged, batched backends agree with
per-LP solves under the same flux rule, and ensembles equal independent runs."""

import numpy as np
import pandas as pd
import pytest

from muode import DynamicFBA
from muode.backends import LegacyBackend, SolveRequest, make_backend
from muode.diet import Diet
from muode.examples import build_toy_community, toy_diet, toy_kinetics
from muode.organism import MuMaxCapped

manylp = pytest.importorskip("manylp")


def _run(backend=None, t_end=12.0, dt=0.1, community=None, diet=None):
    return DynamicFBA(t_end=t_end, dt=dt).run(
        community or build_toy_community(), diet or toy_diet(), toy_kinetics(), backend=backend)


def test_default_is_the_legacy_path_bit_for_bit():
    a = _run()
    b = _run(LegacyBackend(n_jobs=1))
    c = _run(LegacyBackend(n_jobs=4))
    for other in (b, c):
        pd.testing.assert_frame_equal(a.biomass, other.biomass, check_exact=True)
        pd.testing.assert_frame_equal(a.metabolites, other.metabolites, check_exact=True)
        pd.testing.assert_frame_equal(a.growth_rates, other.growth_rates, check_exact=True)


@pytest.mark.parametrize("rule", ["vertex", "pfba", "pfba-unique"])
def test_batched_matches_per_lp_highs_under_the_same_rule(rule):
    batched = _run(make_backend("manylp", flux_rule=rule))
    per_lp = _run(make_backend("highs", flux_rule=rule))
    np.testing.assert_allclose(batched.biomass.values, per_lp.biomass.values, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(batched.metabolites.values, per_lp.metabolites.values, rtol=1e-9, atol=1e-12)


def test_batched_growth_matches_legacy_solver():
    # growth rates are unique, whatever vertex a solver returns
    legacy = _run()
    batched = _run(make_backend("manylp", flux_rule="vertex"))
    np.testing.assert_allclose(batched.growth_rates.values, legacy.growth_rates.values,
                               rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(batched.biomass.values, legacy.biomass.values, rtol=1e-7, atol=1e-12)


def test_batched_backend_certifies_almost_everything():
    be = make_backend("manylp")
    _run(be, t_end=24.0, dt=0.05)
    st = be.stats()
    assert st["lps"] > 800
    assert st["highs_solves"] < 0.05 * st["lps"]


def test_ensemble_members_equal_independent_runs():
    diets = [Diet(concentrations={"glc_e": g, "ac_e": 0.0}, name=f"glc{g}") for g in (5.0, 12.0, 20.0)]
    eng = DynamicFBA(t_end=10.0, dt=0.1)
    be = make_backend("manylp")
    ens = eng.run_ensemble(build_toy_community(), diets, toy_kinetics(), backend=be)
    assert len(ens) == 3
    for d, res in zip(diets, ens):
        solo = eng.run(build_toy_community(), d, toy_kinetics(), backend=make_backend("manylp"))
        np.testing.assert_allclose(res.biomass.values, solo.biomass.values, rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(res.metabolites.values, solo.metabolites.values, rtol=1e-10, atol=1e-14)


def test_ensemble_with_legacy_backend_matches_too():
    diets = [Diet(concentrations={"glc_e": g, "ac_e": 0.0}, name=f"glc{g}") for g in (5.0, 20.0)]
    eng = DynamicFBA(t_end=6.0, dt=0.1)
    ens = eng.run_ensemble(build_toy_community(), diets, toy_kinetics())
    for d, res in zip(diets, ens):
        solo = eng.run(build_toy_community(), d, toy_kinetics())
        pd.testing.assert_frame_equal(res.biomass, solo.biomass, check_exact=True)


def test_mumax_wrapper_is_honoured():
    def capped():
        c = build_toy_community()
        c.organisms = [MuMaxCapped(o, mumax=0.05) for o in c.organisms]
        return c

    legacy = _run(community=capped())
    batched = _run(make_backend("manylp", flux_rule="vertex"), community=capped())
    assert legacy.growth_rates.values.max() <= 0.05 + 1e-12
    np.testing.assert_allclose(batched.growth_rates.values, legacy.growth_rates.values, rtol=1e-7, atol=1e-9)


def test_unbatchable_organism_falls_back_to_legacy():
    class Opaque:
        """An organism the backend cannot extract an LP from."""

        def __init__(self, inner):
            self._o = inner
            self.id = inner.id

        def __getattr__(self, k):
            return getattr(self._o, k)

    comm = build_toy_community()
    comm.organisms = [Opaque(o) for o in comm.organisms]
    be = make_backend("manylp")
    res = _run(be, community=comm)
    ref = _run()
    np.testing.assert_allclose(res.biomass.values, ref.biomass.values, rtol=1e-12)
    assert be.stats()["fallback_lps"] > 0


def test_solve_step_returns_one_solution_per_request():
    comm = build_toy_community()
    reqs = [SolveRequest(o, {"glc_e": 5.0, "ac_e": 1.0}, member=m) for m in range(3) for o in comm.organisms]
    sols = make_backend("manylp").solve_step(reqs)
    assert len(sols) == len(reqs)
    legacy = LegacyBackend().solve_step(reqs)
    for a, b in zip(sols, legacy):
        assert a.status == b.status
        assert a.growth_rate == pytest.approx(b.growth_rate, rel=1e-8, abs=1e-12)


def test_atlas_makes_the_second_run_simplex_free(tmp_path):
    be1 = make_backend("manylp", atlas_dir=str(tmp_path))
    a = _run(be1)
    be1.close()
    assert any(tmp_path.iterdir())
    be2 = make_backend("manylp", atlas_dir=str(tmp_path))
    b = _run(be2)
    assert be2.atlas_loaded > 0
    assert be2.stats()["highs_solves"] == 0
    np.testing.assert_allclose(a.biomass.values, b.biomass.values, rtol=1e-12, atol=1e-15)


def test_policies_give_an_alternative_optima_envelope():
    pols = ["canonical", "max:ac_e", "min:ac_e", "random:0", "random:1"]
    be = make_backend("manylp", policies=pols)
    eng = DynamicFBA(t_end=12.0, dt=0.1)
    runs = eng.run_ensemble(build_toy_community(), toy_diet(), toy_kinetics(), n_members=len(pols),
                            backend=be)
    assert len(runs) == len(pols)
    # growth is unique under every policy, so each trajectory is a valid dFBA solution;
    # the canonical member equals a plain single run
    solo = _run(make_backend("manylp"))
    np.testing.assert_allclose(runs[0].biomass.values, solo.biomass.values, rtol=1e-10, atol=1e-14)


def test_gpu_backend_matches_cpu_on_a_large_ensemble():
    from manylp.backend import gpu_available

    if not gpu_available():
        pytest.skip("needs a GPU")
    diets = [Diet(concentrations={"glc_e": g, "ac_e": 0.0}, name=f"g{g:.2f}")
             for g in np.linspace(2.0, 30.0, 64)]
    eng = DynamicFBA(t_end=8.0, dt=0.1)
    cpu = eng.run_ensemble(build_toy_community(), diets, toy_kinetics(), backend=make_backend("manylp-cpu"))
    gpu = eng.run_ensemble(build_toy_community(), diets, toy_kinetics(), backend=make_backend("manylp-gpu"))
    for a, b in zip(cpu, gpu):
        np.testing.assert_allclose(a.biomass.values, b.biomass.values, rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(a.metabolites.values, b.metabolites.values, rtol=1e-10, atol=1e-14)


def test_spatial_batched_backend_matches_per_lp_highs_and_legacy_growth():
    from muode.spatial import SpatialDynamicFBA

    def run(backend):
        return SpatialDynamicFBA(nx=6, ny=2, t_end=2.0, dt=0.05, backend=backend).run(
            build_toy_community(), toy_diet(), toy_kinetics())

    legacy = SpatialDynamicFBA(nx=6, ny=2, t_end=2.0, dt=0.05).run(
        build_toy_community(), toy_diet(), toy_kinetics())
    batched = run(make_backend("manylp", flux_rule="pfba-unique"))
    per_lp = run(make_backend("highs", flux_rule="pfba-unique"))
    for k in batched.biomass:
        np.testing.assert_allclose(batched.biomass[k], per_lp.biomass[k], rtol=1e-9, atol=1e-14)
    # the toy LPs have unique growth, so biomass agrees with the historical path too
    vertex = run(make_backend("manylp", flux_rule="vertex"))
    for k in legacy.biomass:
        np.testing.assert_allclose(vertex.biomass[k], legacy.biomass[k], rtol=1e-6, atol=1e-12)
