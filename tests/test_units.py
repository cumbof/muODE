"""Unit tests for the supporting components (kinetics, diet, community)."""

import math

import pytest

from muode import Community, Diet, KineticParameters
from muode.diet import available_presets, load_preset
from muode.examples import build_glucose_specialist, build_acetate_specialist
from muode.kinetics import vmax_from_kcat


# -- kinetics ---------------------------------------------------------------


def test_michaelis_menten_monotonic_and_saturating():
    kp = KineticParameters(default_vmax=10.0, default_km=1.0)
    v_low = kp.michaelis_menten("o", "m", 0.1)
    v_mid = kp.michaelis_menten("o", "m", 1.0)
    v_hi = kp.michaelis_menten("o", "m", 1e6)
    assert 0 < v_low < v_mid < v_hi
    assert v_mid == pytest.approx(5.0)         # at M == Km, v == Vmax/2
    assert v_hi == pytest.approx(10.0, rel=1e-3)  # saturates at Vmax


def test_kinetics_override_precedence():
    kp = KineticParameters(metabolite_defaults={"glc_e": (8.0, 0.2)})
    assert kp.get("any", "glc_e") == (8.0, 0.2)
    kp.set("orgX", "glc_e", 99.0, 0.9)
    assert kp.get("orgX", "glc_e") == (99.0, 0.9)   # override beats default
    assert kp.get("orgY", "glc_e") == (8.0, 0.2)


def test_zero_concentration_gives_zero_uptake():
    kp = KineticParameters()
    assert kp.michaelis_menten("o", "m", 0.0) == 0.0


def test_vmax_from_kcat_units():
    # 10 /s * 3600 s/h * 0.001 mmol/gDW = 36 mmol/gDW/h
    assert vmax_from_kcat(10.0, 0.001) == pytest.approx(36.0)


# -- diet -------------------------------------------------------------------


def test_diet_presets_exist():
    assert "western_gut" in available_presets()
    d = load_preset("western_gut")
    assert d.initial_concentration("glc_e") > 0


def test_diet_csv_roundtrip(tmp_path):
    d = Diet({"glc_e": 10.0, "ac_e": 0.0}, influx={"glc_e": 1.0}, name="t")
    path = tmp_path / "diet.csv"
    d.to_csv(path)
    d2 = Diet.from_csv(path)
    assert d2.initial_concentration("glc_e") == 10.0
    assert d2.influx_rate("glc_e") == 1.0


# -- community --------------------------------------------------------------


def test_abundance_normalisation_and_initial_biomass():
    comm = Community(
        organisms=[build_glucose_specialist(), build_acetate_specialist()],
        abundances={"A_glucose": 3.0, "B_acetate": 1.0},
        total_biomass=0.04,
    )
    assert math.isclose(sum(comm.abundances.values()), 1.0)
    init = comm.initial_biomass()
    assert init["A_glucose"] == pytest.approx(0.03)
    assert init["B_acetate"] == pytest.approx(0.01)


def test_duplicate_ids_rejected():
    with pytest.raises(ValueError):
        Community(organisms=[build_glucose_specialist(), build_glucose_specialist()])


def test_environment_metabolites_union():
    comm = Community([build_glucose_specialist(), build_acetate_specialist()])
    assert set(comm.environment_metabolites()) == {"glc_e", "ac_e"}
