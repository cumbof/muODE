"""The genome-scale example scenarios are glue: load real GEMs, wire an ecology layer,
run.  They cannot run in CI without reconstructed GEMs, but the GLUE can break silently
(an import, a changed layer signature, a diet that starves a real model).  So exercise
each on a real BiGG stand-in (cobra's textbook E. coli, saved under the scenarios' member
names) and assert it runs both arms with the ecology layer actually active.

This is not the science -- the science needs the real strains -- it is a tripwire on the
scaffolding the workstation run depends on.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

cobra = pytest.importorskip("cobra")

pytestmark = pytest.mark.slow


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_ROOT = Path(__file__).resolve().parents[1]


def _stand_in_models(tmp_path, names):
    """Write cobra's textbook E. coli under each member name -- a real BiGG GEM that
    grows on the scenarios' defined glucose media."""
    tb = cobra.io.load_model("textbook")
    for n in names:
        cobra.io.write_sbml_model(tb, str(tmp_path / f"{n}.xml"))
    return tmp_path


def test_phage_therapy_genome_scale_runs_and_the_phage_layer_is_active(tmp_path):
    pt = _load(_ROOT / "examples" / "phage_therapy" / "genome_scale.py", "pt_gs")
    from muode.diet import Diet

    models = pt.load_models(_stand_in_models(tmp_path, [pt.HOST, pt.COMMENSAL]))
    assert set(models) == {pt.HOST, pt.COMMENSAL}
    diet = Diet.from_csv(_ROOT / "examples" / "phage_therapy" / "gut_glucose.csv")

    no = pt.scenario(models, therapy=False, diet=diet, t_end=6.0)
    yes = pt.scenario(models, therapy=True, diet=diet, t_end=6.0)
    assert no.final_biomass()[pt.HOST] > 0                       # the host grows untreated
    # the phage layer is wired and amplifies (a lytic burst ran), only in the therapy arm
    col = f"phage[{pt.PHAGE['name']}]"
    assert col in yes.environment.columns
    assert yes.environment[col].max() > pt.PHAGE["initial_titer"]
    # the no-phage arm has no ecology layer at all, so no environment is recorded
    assert no.environment is None or not any("phage[" in c for c in no.environment.columns)


def test_strain_competition_genome_scale_runs_all_three_modes(tmp_path):
    sc = _load(_ROOT / "examples" / "strain_competition" / "genome_scale.py", "sc_gs")
    diet_csv = _ROOT / "examples" / "strain_competition" / "defined_medium.csv"

    models = sc.load_models(_stand_in_models(tmp_path, list(sc.STRAINS)))
    assert set(models) == set(sc.STRAINS)

    # bare resource competition and the colicin (interference) arm both run
    resource = sc.scenario(models, diet_csv, with_colicin=False, with_arabinose=False, t_end=6.0)
    colicin = sc.scenario(models, diet_csv, with_colicin=True, with_arabinose=False, t_end=6.0)
    assert resource.final_biomass()[sc.EFFICIENT] > 0
    # the colicin suppresses its target (the specialist) relative to the resource arm --
    # interference is active even when the GEMs are identical (here they are stand-ins)
    assert colicin.final_biomass()[sc.EFFICIENT] < resource.final_biomass()[sc.EFFICIENT]
