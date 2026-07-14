"""A zero growth rate owes us a reason.

Two causes look identical in a results table:

* the genome genuinely lacks the pathways -- a true negative, and the answer;
* the medium lacks a nutrient the model's biomass demands -- an artefact, where the
  model was dead on arrival and its zero says nothing about the organism.

Fe(III) on DM38 is the instance we tripped over, but the point is that the framework
takes ARBITRARY genomes: which nutrient goes missing depends entirely on what you
feed it, so it cannot be a fact anyone remembers to check.  It has to be computed,
every run, on whatever models are in hand.
"""

from __future__ import annotations

import cobra
import pytest

from muode.diet import Diet, load_diet
from muode.qc import (
    FORBIDDEN_RESCUES,
    diagnose_no_growth,
    sanity_check_model,
)


@pytest.mark.slow
def test_a_model_that_grows_is_not_diagnosed(ijo):
    rich = Diet(concentrations={m[len("EX_"):]: 10.0 for m in ijo.medium},
                name="its own medium")
    assert diagnose_no_growth(ijo, rich)["verdict"] == "grows"


@pytest.mark.slow
def test_dm38_is_diagnosed_as_a_medium_gap_and_names_the_iron(dm38_diagnosis):
    """The real case.  iJO1366 scores 0.0 on DM38 -- and that zero is not biology."""
    d = dm38_diagnosis

    assert d["verdict"] == "medium_gap", (
        "iJO1366 on DM38 is a medium artefact, not a finding about E. coli"
    )
    assert "EX_fe3_e" in d["rescued_by_any_of"], (
        f"ferric iron should be on the menu; got {d['rescued_by_any_of']}"
    )
    assert d["growth"] <= 1e-6 < d["growth_if_rescued"]


@pytest.mark.slow
def test_the_menu_exposes_the_artefact_route_instead_of_silently_picking_it(dm38_diagnosis):
    """The reason ``rescued_by_any_of`` is a LIST and not a single minimal answer.

    iJO1366's only anaerobic route to cytoplasmic Fe(III) other than importing it is
    FESD2s -- ``4fe4s_c + no_c -> 3fe4s_c + fe3_c`` -- in which nitric oxide destroys
    the cell's own iron-sulfur clusters to liberate iron.  FBA is delighted to use it.

    A shrink-to-minimal diagnosis returns exactly one rescue and picks between fe3_e
    and no_e on iteration order; ours returned NO, which is a true minimal rescue, a
    solver artefact, and useless to a human.  Enumerating the whole equivalence class
    puts both on the table so a person can tell chemistry from damage.  The tool must
    not make that call, but it must not hide it either.
    """
    d = dm38_diagnosis
    assert {"EX_fe3_e", "EX_no_e"} <= set(d["rescued_by_any_of"]), (
        "both the honest fix and the artefact fix must be visible, not one of them"
    )


@pytest.mark.slow
def test_the_menu_is_short_enough_to_read(dm38_diagnosis):
    """DM38 leaves dozens of iJO1366's exchanges closed; listing them all says nothing."""
    d = dm38_diagnosis
    assert 0 < len(d["rescued_by_any_of"]) < d["n_candidates"] / 4


def test_redundant_routes_do_not_make_the_gap_vanish():
    """The bug a naive implementation ships.

    Give a model two ways to satisfy one requirement -- x OR y -- and ask "is any
    single missing nutrient *essential*?"  The answer is NO for both, because each is
    dispensable while the other is available.  A diagnosis built on essentiality
    reports an empty set and declares the medium innocent while the model sits at
    zero.  Asking "does this one ALONE rescue it?" gets both.
    """
    m = cobra.Model("redundant")
    a, b, biomass = (cobra.Metabolite(i, compartment="c") for i in ("a_c", "b_c", "bm_c"))

    # biomass needs `a`, and `a` can be made from EITHER of two nutrients
    for nutrient in ("x_e", "y_e"):
        met = cobra.Metabolite(nutrient, compartment="e")
        ex = cobra.Reaction(f"EX_{nutrient}")
        ex.add_metabolites({met: -1}); ex.bounds = (0, 1000)
        conv = cobra.Reaction(f"to_a_from_{nutrient}")
        conv.add_metabolites({met: -1, a: 1}); conv.bounds = (0, 1000)
        m.add_reactions([ex, conv])

    # ...and it also needs `b`, which has only one source: the one the diet DOES supply
    z = cobra.Metabolite("z_e", compartment="e")
    exz = cobra.Reaction("EX_z_e"); exz.add_metabolites({z: -1}); exz.bounds = (0, 1000)
    toz = cobra.Reaction("to_b"); toz.add_metabolites({z: -1, b: 1}); toz.bounds = (0, 1000)
    bio = cobra.Reaction("BIOMASS")
    bio.add_metabolites({a: -1, b: -1, biomass: 1}); bio.bounds = (0, 1000)
    sink = cobra.Reaction("EX_bm_c"); sink.add_metabolites({biomass: -1}); sink.bounds = (0, 1000)
    m.add_reactions([exz, toz, bio, sink])
    m.objective = "BIOMASS"

    diet = Diet(concentrations={"z_e": 10.0}, max_uptake={"z_e": 10.0}, name="missing x and y")
    d = diagnose_no_growth(m, diet)

    assert d["verdict"] == "medium_gap", "the gap must not vanish just because it has two fixes"
    assert set(d["rescued_by_any_of"]) == {"EX_x_e", "EX_y_e"}, (
        f"BOTH routes fix it and both must be listed; got {d['rescued_by_any_of']}"
    )


@pytest.mark.slow
def test_oxygen_is_never_proposed_as_a_rescue(dm38_diagnosis):
    """A diagnosis may complete an environment.  It may not overturn one.

    Oxygen rescues almost any GEM, so if it were a candidate it would be the answer
    every time and drown every real finding -- and it would silently turn an
    anaerobic gut simulation aerobic.  A diet that omits o2_e is making a claim about
    the environment, and the diagnosis has to respect it.
    """
    assert "o2_e" in FORBIDDEN_RESCUES
    d = dm38_diagnosis
    assert "EX_o2_e" not in d["rescued_by_any_of"]


def test_a_model_that_cannot_grow_on_anything_exonerates_the_medium():
    """Blame the reconstruction, not the diet, when nothing would help."""
    m = cobra.Model("broken")
    bm = cobra.Metabolite("bm_c", compartment="c")
    bio = cobra.Reaction("BIOMASS")
    bio.add_metabolites({bm: -1})          # needs bm, and nothing on earth makes it
    bio.bounds = (0, 1000)
    m.add_reactions([bio])
    m.objective = "BIOMASS"

    d = diagnose_no_growth(m, Diet(concentrations={"z_e": 10.0}, name="anything"))
    assert d["verdict"] == "model_cannot_grow"
    assert d["rescued_by_any_of"] == [] and d["rescued_by_all_of"] == []


@pytest.mark.slow
def test_the_pipeline_surfaces_the_verdict_without_being_asked(ijo):
    """sanity_check_model diagnoses a zero automatically when it is given the diet.

    If this has to be remembered, it will not be remembered.
    """
    report = sanity_check_model(ijo, growth_rate=0.0, diet=load_diet("dm38"))
    assert report["no_growth"]["verdict"] == "medium_gap"
    assert "EX_fe3_e" in report["no_growth"]["rescued_by_any_of"]


@pytest.mark.slow
def test_a_growing_model_is_not_diagnosed_by_the_pipeline(ijo):
    report = sanity_check_model(ijo, growth_rate=0.5, diet=load_diet("dm38"))
    assert "no_growth" not in report
