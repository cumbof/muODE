"""Which of a medium's numbers actually matter?

`western_gut` has 126 bounded rows and 71 of them are the source table's *fill
value* -- not measurements.  The instinct is to go and find 71 literature values.
Most of that work would be wasted, because **nothing grows on zinc**: a cell needs
a trace of it and is saturated long before the bound bites, so whether the row says
0.1 or 0.001 cannot change a prediction.

The rows that DO matter are the ones growth is sensitive to, and which ones those
are depends on the genomes -- so it is computed, not remembered.  These tests pin
the three verdicts on models small enough that the right answer is known by
construction, then check the real medium behaves the way the design assumes.
"""

from __future__ import annotations

import cobra
import pytest

from muode.diet import Diet, load_diet
from muode.qc import nutrient_sensitivity


def _toy():
    """A cell that eats carbon, needs a trace of a metal, and ignores a third nutrient.

    biomass = 1 carbon + 0.001 metal.  `junk` is importable and goes nowhere.
    """
    m = cobra.Model("toy")
    bm = cobra.Metabolite("bm_c", compartment="c")
    bio = cobra.Reaction("BIOMASS")
    bio.bounds = (0, 1000)
    m.add_reactions([bio])

    for nutrient, coeff in (("c_e", -1.0), ("metal_e", -0.001), ("junk_e", 0.0)):
        met = cobra.Metabolite(nutrient, compartment="e")
        ex = cobra.Reaction(f"EX_{nutrient}")
        ex.add_metabolites({met: -1})
        ex.bounds = (-1000, 1000)
        m.add_reactions([ex])
        if coeff:
            bio.add_metabolites({met: coeff})

    bio.add_metabolites({bm: 1})
    sink = cobra.Reaction("DM_bm")
    sink.add_metabolites({bm: -1})
    sink.bounds = (0, 1000)
    m.add_reactions([sink])
    m.objective = "BIOMASS"
    return m


@pytest.fixture(scope="module")
def toy_report():
    diet = Diet(
        concentrations={"c_e": 100.0, "metal_e": 100.0, "junk_e": 100.0},
        max_uptake={"c_e": 1.0, "metal_e": 0.1, "junk_e": 0.1},
        name="toy",
    )
    return nutrient_sensitivity(_toy(), diet)


def test_the_carbon_source_is_load_bearing(toy_report):
    """Growth is 1:1 with the carbon bound, so the number IS the prediction."""
    c = toy_report["nutrients"]["EX_c_e"]
    assert c["verdict"] == "load_bearing"
    assert c["essential"] and c["limiting"]
    # and the mechanism: the cell eats every bit of it, right up to the ceiling
    assert c["uptake_at_optimum"] == pytest.approx(c["bound"])


def test_the_metal_is_needed_but_its_number_is_irrelevant(toy_report):
    """The zinc case, in miniature -- and the reason this tool exists.

    Biomass needs the metal at a 1:1000 stoichiometry, so the cell is saturated far
    below the 0.1 bound: removing it kills growth, but its exact value cannot move a
    prediction.  A fill value here is HARMLESS, and this is what says so.
    """
    metal = toy_report["nutrients"]["EX_metal_e"]
    assert metal["verdict"] == "essential_trace"
    assert metal["essential"], "biomass demands it -- removing it must stop growth"
    assert not metal["limiting"], "but 10x more of it must not buy any growth"
    assert metal["uptake_at_optimum"] < metal["bound"], "the bound is slack, not binding"


def test_a_nutrient_the_model_cannot_use_is_called_unused(toy_report):
    junk = toy_report["nutrients"]["EX_junk_e"]
    assert junk["verdict"] == "unused"
    assert not junk["essential"] and not junk["limiting"]


def test_a_dead_model_is_refused_rather_than_mischaracterised():
    """Every nutrient in a dead model looks irrelevant, because nothing depends on
    anything.  Reporting 71 `unused` rows would be worse than useless -- it would be
    a confident wrong answer.  diagnose_no_growth is the tool for a zero; say so.
    """
    m = cobra.Model("dead")
    bm = cobra.Metabolite("bm_c", compartment="c")
    bio = cobra.Reaction("BIOMASS")
    bio.add_metabolites({bm: -1})       # needs bm; nothing makes bm
    bio.bounds = (0, 1000)
    m.add_reactions([bio])
    m.objective = "BIOMASS"

    r = nutrient_sensitivity(m, Diet(concentrations={"c_e": 10.0}, name="whatever"))
    assert r["verdict"] == "no_growth"
    assert r["nutrients"] == {}


def test_the_minerals_in_western_gut_are_not_load_bearing():
    """The design assumption behind leaving the mineral rows at the fill value.

    If a metal ever came back `load_bearing`, the medium would be deciding growth
    with someone else's default -- and that row would need a real number urgently.
    This is the guard that would catch it.
    """
    model = cobra.io.load_model("iJO1366")
    r = nutrient_sensitivity(model, load_diet("western_gut"))
    assert r["verdict"] == "ok", "iJO1366 must grow on western_gut for this to mean anything"

    minerals = ["EX_zn2_e", "EX_cu2_e", "EX_mn2_e", "EX_ca2_e", "EX_mg2_e",
                "EX_k_e", "EX_cl_e", "EX_cobalt2_e", "EX_mobd_e"]
    offenders = {
        ex: r["nutrients"][ex]["verdict"]
        for ex in minerals
        if ex in r["nutrients"] and r["nutrients"][ex]["limiting"]
    }
    assert not offenders, (
        f"a mineral bound is limiting growth, so its fill value is NOT harmless: {offenders}"
    )
