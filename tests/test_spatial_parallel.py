"""Parallel per-cell solves in the spatial engine must match sequential ones.

The 2D reaction-diffusion engine can fan the independent per-cell FBA solves of
each step across worker threads (``n_jobs``).  Because every populated cell reads
only its own slice of the biomass and concentration fields and each worker holds
its own private model copies, the fields must not depend on ``n_jobs`` -- only the
wall-clock time does.  These tests pin that invariant on the toy cross-feeding
gradient (a multi-cell, multi-species run that exercises the parallel path).

Kept deliberately short.  This is an *equivalence* test: it needs enough steps for
the cells to grow, cross-feed and diffuse into each other, and not one step more.
It used to run 300 steps per simulation and compute the sequential run twice (once
per test), which cost 109s of a 513s suite to prove something 50 steps prove just as
well.  ``test_the_run_is_not_trivially_short`` is the guard that stops anyone --
including a future me, chasing the clock -- from shortening it into vacuity.
"""

import numpy as np
import pytest

from muode.examples import build_toy_community, toy_diet, toy_kinetics
from muode.spatial import SpatialDynamicFBA, halves_inoculum

#: 200 steps, down from 300.  Not chosen to be round: this is the point at which the
#: SECOND species starts growing, i.e. the point at which cross-feeding has actually
#: happened.  Below it (t_end=2.0, 3.0) only the glucose eater grows, the acetate eater
#: never gets any acetate, and the "cross-feeding gradient" this module claims to test
#: is not yet a cross-feeding gradient at all.  test_the_run_is_not_trivially_short is
#: what holds this honest -- it fails at t_end=1.0.
T_END = 4.0
DT = 0.02
AMOUNT = 0.04


def _run(n_jobs):
    comm = build_toy_community()
    shape = (1, 8)
    inoc = halves_inoculum(shape, "A_glucose", "B_acetate", amount=AMOUNT, axis=1)
    engine = SpatialDynamicFBA(nx=8, ny=1, dx=1.0, t_end=T_END, dt=DT,
                               default_diffusivity=2.0, n_jobs=n_jobs)
    return engine.run(comm, toy_diet(), toy_kinetics(), inoculum=inoc)


@pytest.fixture(scope="module")
def sequential():
    """The n_jobs=1 reference, computed ONCE and shared by both comparisons."""
    return _run(1)


def _assert_fields_equal(a, b):
    for sp in a.biomass:
        np.testing.assert_array_equal(a.biomass[sp], b.biomass[sp])
    for m in a.metabolites:
        np.testing.assert_array_equal(a.metabolites[m], b.metabolites[m])


@pytest.mark.slow
def test_the_run_is_not_trivially_short(sequential):
    """The equivalence tests below are only worth anything if the run DID something.

    Two fields of zeros are equal to two other fields of zeros.  So before comparing
    parallel against sequential, establish that the simulation actually exercised the
    thing it claims to: **both** species grew, which on this inoculum can only happen
    if the acetate eater received acetate the glucose eater excreted and diffused to
    it.  That is the cross-feeding gradient, and it is the workload the parallel path
    has to reproduce.

    This guard is not decoration.  Shortening T_END to 1.0 to make the suite faster
    makes it fail -- only the glucose eater grows, and the parallel tests would then be
    comparing two runs in which nothing interesting had happened yet.
    """
    grew = [sp for sp, field in sequential.biomass.items() if field.max() > AMOUNT * 1.01]
    assert len(grew) >= 2, (
        f"only {grew} grew -- with one species growing there is no cross-feeding, and "
        "the parallel comparison is testing an empty simulation"
    )

    spread = [
        m for m, field in sequential.metabolites.items()
        if np.count_nonzero(field) > 1 and field.max() > 0
    ]
    assert spread, "no metabolite reached more than one cell -- nothing diffused"


@pytest.mark.slow
def test_spatial_threaded_matches_sequential(sequential):
    _assert_fields_equal(sequential, _run(4))


@pytest.mark.slow
def test_spatial_all_cores_matches_sequential(sequential):
    _assert_fields_equal(sequential, _run(-1))
