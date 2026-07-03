"""Parallel per-step solves must be numerically identical to sequential ones.

The dynamic-FBA engine can fan the independent per-species LP solves of each step
across worker threads (``n_jobs``).  Because the solves share only the read-only
medium and biomass snapshots within a step and mutate only their own organism's
bounds, the result must not depend on ``n_jobs`` at all -- parallelism buys speed
on genome-scale models, never a different answer.  These tests pin that contract.
"""

import pandas as pd

from muode.dfba import DynamicFBA
from muode.ecology import EcologyModel
from muode.examples import build_toy_community, toy_diet, toy_kinetics
from muode.ph import WeakAcidInhibition


def _run(n_jobs, ecology=None):
    comm = build_toy_community()
    engine = DynamicFBA(t_end=10.0, dt=0.1, n_jobs=n_jobs)
    return engine.run(comm, toy_diet(), toy_kinetics(), ecology=ecology)


def test_threaded_matches_sequential():
    seq = _run(1)
    par = _run(4)
    pd.testing.assert_frame_equal(seq.biomass, par.biomass)
    pd.testing.assert_frame_equal(seq.metabolites, par.metabolites)
    pd.testing.assert_frame_equal(seq.growth_rates, par.growth_rates)


def test_all_cores_matches_sequential():
    seq = _run(1)
    par = _run(-1)          # -1 == all available cores
    pd.testing.assert_frame_equal(seq.biomass, par.biomass)
    pd.testing.assert_frame_equal(seq.metabolites, par.metabolites)


def test_cross_feeding_is_identical():
    key = ["producer", "metabolite", "consumer"]
    seq = _run(1).cross_feeding().sort_values(key).reset_index(drop=True)
    par = _run(4).cross_feeding().sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(seq, par)


def test_threaded_matches_sequential_with_ecology():
    # An active ecology layer adds read-only per-step hooks (here, SCFA-driven
    # acidification from the community's secreted acetate); threading must leave
    # both the trajectories and the environmental observables unchanged.
    eco = lambda: EcologyModel([WeakAcidInhibition(buffer_capacity=25.0, default_ki=12.0)])
    seq = _run(1, ecology=eco())
    par = _run(4, ecology=eco())
    pd.testing.assert_frame_equal(seq.biomass, par.biomass)
    pd.testing.assert_frame_equal(seq.metabolites, par.metabolites)
    assert seq.environment is not None
    pd.testing.assert_frame_equal(seq.environment, par.environment)
