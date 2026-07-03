"""Parallel per-cell solves in the spatial engine must match sequential ones.

The 2D reaction-diffusion engine can fan the independent per-cell FBA solves of
each step across worker threads (``n_jobs``).  Because every populated cell reads
only its own slice of the biomass and concentration fields and each worker holds
its own private model copies, the fields must not depend on ``n_jobs`` -- only the
wall-clock time does.  These tests pin that invariant on the toy cross-feeding
gradient (a multi-cell, multi-species run that exercises the parallel path).
"""

import numpy as np

from muode.examples import build_toy_community, toy_diet, toy_kinetics
from muode.spatial import SpatialDynamicFBA, halves_inoculum


def _run(n_jobs):
    comm = build_toy_community()
    shape = (1, 8)
    inoc = halves_inoculum(shape, "A_glucose", "B_acetate", amount=0.04, axis=1)
    engine = SpatialDynamicFBA(nx=8, ny=1, dx=1.0, t_end=6.0, dt=0.02,
                               default_diffusivity=2.0, n_jobs=n_jobs)
    return engine.run(comm, toy_diet(), toy_kinetics(), inoculum=inoc)


def _assert_fields_equal(a, b):
    for sp in a.biomass:
        np.testing.assert_array_equal(a.biomass[sp], b.biomass[sp])
    for m in a.metabolites:
        np.testing.assert_array_equal(a.metabolites[m], b.metabolites[m])


def test_spatial_threaded_matches_sequential():
    _assert_fields_equal(_run(1), _run(4))


def test_spatial_all_cores_matches_sequential():
    _assert_fields_equal(_run(1), _run(-1))
