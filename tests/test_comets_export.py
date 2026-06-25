"""Tests for the muODE -> COMETS export bridge (`muode.comets`).

Pure stdlib: the exporter operates on a community manifest + a Diet and writes
text/JSON, so it needs neither cobra nor COMETS and runs on any machine.
"""

import json

from muode.comets import export_comets, initial_biomass
from muode.diet import Diet


def _manifest():
    return {
        "models": ["/models/A_glucose.xml", "/models/B_acetate.xml"],
        "abundances": {"A_glucose": 3.0, "B_acetate": 1.0},
        "total_biomass": 0.02,
        "diet": "custom",
    }


def test_initial_biomass_splits_total_by_abundance():
    bm = initial_biomass(_manifest())
    # 3:1 abundance over total 0.02 gDW
    assert bm["A_glucose"] == 0.015
    assert bm["B_acetate"] == 0.005
    assert sum(bm.values()) == 0.02


def test_export_writes_layout_params_driver(tmp_path):
    diet = Diet({"glc_e": 10.0, "ac_e": 0.0}, name="toy")
    paths = export_comets(_manifest(), tmp_path, diet=diet, t_end=24.0, dt=0.1)

    for key in ("layout", "params", "driver", "readme"):
        assert paths[key].exists()

    layout = json.loads(paths["layout"].read_text())
    names = [name for _path, name in layout["models"]]
    assert names == ["A_glucose", "B_acetate"]
    assert layout["initial_biomass"]["A_glucose"] == 0.015
    assert layout["media"]["glc_e"] == 10.0
    assert layout["max_cycles"] == 240          # 24 h / 0.1 h
    assert layout["time_step"] == 0.1

    params = paths["params"].read_text()
    assert "maxCycles = 240" in params
    assert "timeStep = 0.1" in params

    driver = paths["driver"].read_text()
    assert "import cometspy" in driver
    assert "comets_layout.json" in driver


def test_export_resolves_manifest_diet_when_none(tmp_path):
    # diet=None -> resolve the manifest's preset name ("western_gut")
    manifest = dict(_manifest(), diet="western_gut")
    paths = export_comets(manifest, tmp_path, diet=None)
    layout = json.loads(paths["layout"].read_text())
    assert layout["diet"] == "western_gut_demo"
    assert "glc_e" in layout["media"]
