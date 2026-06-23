"""M4 -- spatiotemporal dynamic FBA (2D reaction-diffusion colony/biofilm).

Dependency-light: the toy linprog community is solved per grid cell, so these run
without cobra or a solver licence.
"""

import numpy as np
import pytest

from muode.community import Community
from muode.diet import Diet
from muode.examples import (
    build_acetate_specialist,
    build_glucose_specialist,
    toy_kinetics,
)
from muode.spatial import (
    SpatialDynamicFBA,
    diffuse_step,
    halves_inoculum,
    point_inoculum,
)


def test_diffusion_conserves_mass_and_smooths():
    field = np.zeros((1, 11))
    field[0, 5] = 10.0
    total0 = field.sum()
    for _ in range(80):
        field = diffuse_step(field, diffusivity=1.0, dx=1.0, dt=0.1)  # D*dt/dx^2 = 0.1
    assert field.sum() == pytest.approx(total0)          # no-flux walls conserve mass
    assert field[0, 5] < 10.0                            # peak spread out
    assert field[0, 0] > 0.0                             # reached the boundary


def test_colony_grows_and_depletes_local_nutrient():
    comm = Community([build_glucose_specialist()])
    diet = Diet({"glc_e": 20.0, "ac_e": 0.0}, name="glc")
    shape = (1, 9)
    inoc = point_inoculum(shape, {"A_glucose": (0, 4, 0.05)})
    eng = SpatialDynamicFBA(nx=9, ny=1, dx=1.0, t_end=4.0, dt=0.02,
                            diffusivity={"glc_e": 2.0, "ac_e": 2.0})
    res = eng.run(comm, diet, toy_kinetics(), inoculum=inoc)

    tb = res.total_biomass()["A_glucose"]
    assert tb.iloc[-1] > tb.iloc[0]                       # the colony grew
    glc = res.final_metabolite("glc_e")
    assert glc[0, 4] < glc[0, 0]                          # glucose drawn down at the colony


def test_spatial_cross_feeding_gradient():
    comm = Community([build_glucose_specialist(), build_acetate_specialist()])
    diet = Diet({"glc_e": 20.0, "ac_e": 0.0}, name="glc")
    shape = (1, 8)
    inoc = halves_inoculum(shape, "A_glucose", "B_acetate", amount=0.04, axis=1)
    eng = SpatialDynamicFBA(nx=8, ny=1, dx=1.0, t_end=6.0, dt=0.02,
                            diffusivity={"glc_e": 2.0, "ac_e": 3.0})
    res = eng.run(comm, diet, toy_kinetics(), inoculum=inoc)

    B = res.biomass["B_acetate"]
    assert B[-1].sum() > B[0].sum()                       # B grew on cross-fed acetate
    final_B = B[-1][0]
    assert final_B[4] > final_B[7]                        # more growth near the A interface


def test_unstable_diffusion_warns():
    comm = Community([build_glucose_specialist()])
    diet = Diet({"glc_e": 10.0}, name="glc")
    eng = SpatialDynamicFBA(nx=5, ny=1, dx=1.0, t_end=0.1, dt=0.1,
                            diffusivity={"glc_e": 10.0})   # D*dt/dx^2 = 1.0 >> 0.25
    with pytest.warns(RuntimeWarning, match="unstable"):
        eng.run(comm, diet, toy_kinetics(),
                inoculum=point_inoculum((1, 5), {"A_glucose": (0, 2, 0.01)}))
