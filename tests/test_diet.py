"""The `western_gut` diet must be able to feed a real genome-scale model.

Regression tests for the bug that made the first real CarveMe community (89 MAGs,
all of them healthy and growing on their own complete medium) produce a perfectly
flat 24 h run: the `western_gut` preset was a six-metabolite toy with no nitrogen,
phosphate, sulfur or ion source, and its glucose id (`glc_e`) was not even valid
BiGG. Every organism got a zero uptake bound, so FBA correctly returned zero
growth -- and nothing in the engine said a word about it.
"""

import warnings

import pytest

from muode.dfba import _warn_if_medium_cannot_feed
from muode.diet import Diet, load_preset


class _FakeOrganism:
    """Minimal stand-in exposing just the exchange surface the check reads."""

    def __init__(self, oid, mets):
        self.id = oid
        self._mets = list(mets)

    def exchange_metabolites(self):
        return self._mets


# -- the diet itself --------------------------------------------------------


def test_western_gut_is_a_complete_medium():
    """A genome-scale biomass reaction needs N, P, S and trace ions, not just sugar."""
    d = load_preset("western_gut")
    for essential in ("nh4_e", "pi_e", "so4_e", "k_e", "mg2_e", "ca2_e", "fe2_e", "h2o_e"):
        assert d.initial_concentration(essential) > 0, f"{essential} missing from the diet"
    # carbon, amino acids and vitamins too
    assert d.initial_concentration("glc__D_e") > 0
    assert d.initial_concentration("ala__L_e") > 0
    assert d.initial_concentration("btn_e") > 0


def test_western_gut_is_anaerobic():
    d = load_preset("western_gut")
    assert d.initial_concentration("o2_e") == 0.0
    assert d.influx_rate("o2_e") == 0.0


def test_fermentation_products_are_cross_fed_not_supplied():
    """SCFAs must be produced by the community, not handed to it."""
    d = load_preset("western_gut")
    for product in ("ac_e", "but_e", "ppa_e", "lac__D_e", "succ_e"):
        assert d.initial_concentration(product) == 0.0
        assert d.influx_rate(product) == 0.0


def test_diet_csv_tolerates_comments_and_blank_lines(tmp_path):
    path = tmp_path / "d.csv"
    path.write_text(
        "metabolite,concentration,influx\n"
        "# a comment\n"
        "\n"
        "glc__D_e,10.0,1.0\n"
        "  # indented comment\n"
        "nh4_e,5.0,0.0\n"
    )
    d = Diet.from_csv(path)
    assert d.initial_concentration("glc__D_e") == 10.0
    assert d.influx_rate("glc__D_e") == 1.0
    assert d.initial_concentration("nh4_e") == 5.0
    assert len(d.concentrations) == 2


# -- the guard that should have caught it -----------------------------------


def test_warns_when_diet_namespace_does_not_match_models():
    """The exact failure: a bare `glc_e` diet against BiGG (`glc__D_e`) models."""
    orgs = [_FakeOrganism("bin.1", ["glc__D_e", "nh4_e", "pi_e"])]
    broken = Diet({"glc_e": 10.0}, name="broken")
    with pytest.warns(RuntimeWarning, match="different namespaces"):
        _warn_if_medium_cannot_feed(orgs, broken)


def test_warns_when_diet_has_carbon_but_no_nitrogen_phosphorus_or_sulfur():
    """The actual bug: sugars present, so nothing looks obviously wrong, but a
    genome-scale biomass reaction cannot fire without N/P/S -> flat run."""
    orgs = [_FakeOrganism("bin.1", ["fru_e", "lcts_e", "nh4_e", "pi_e", "so4_e"])]
    sugars_only = Diet({"fru_e": 5.0, "lcts_e": 2.0}, name="sugars_only")
    with pytest.warns(RuntimeWarning) as rec:
        _warn_if_medium_cannot_feed(orgs, sugars_only)
    msgs = " ".join(str(w.message) for w in rec)
    assert "no nitrogen source" in msgs
    assert "no phosphorus source" in msgs
    assert "no sulfur source" in msgs


def test_obligate_cross_feeder_is_not_flagged():
    """A species fed only by the community (acetate is produced, not supplied) is
    legitimate -- it must not be warned about."""
    orgs = [
        _FakeOrganism("producer", ["glc__D_e", "nh4_e", "pi_e", "so4_e", "ac_e"]),
        _FakeOrganism("cross_feeder", ["ac_e"]),
    ]
    diet = Diet(
        {"glc__D_e": 10.0, "nh4_e": 10.0, "pi_e": 5.0, "so4_e": 2.0, "ac_e": 0.0},
        name="complete",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _warn_if_medium_cannot_feed(orgs, diet)  # must not raise


def test_toy_models_without_npS_exchanges_are_not_flagged():
    """A toy model with a lumped biomass has no nh4/pi/so4 uptake route at all;
    the absence of those nutrients from its diet is not an error."""
    orgs = [_FakeOrganism("toy", ["glc_e", "ac_e"])]
    toy = Diet({"glc_e": 10.0, "ac_e": 0.0}, name="toy_glucose")
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _warn_if_medium_cannot_feed(orgs, toy)  # must not raise


def test_warns_when_diet_supplies_nothing_at_all():
    orgs = [_FakeOrganism("bin.1", ["glc__D_e"])]
    empty = Diet({"ac_e": 0.0, "but_e": 0.0}, name="all_zero")
    with pytest.warns(RuntimeWarning, match="no nutrient"):
        _warn_if_medium_cannot_feed(orgs, empty)


def test_real_diet_against_bigg_models_is_silent():
    """The fixed diet must feed a BiGG model without complaint."""
    d = load_preset("western_gut")
    orgs = [_FakeOrganism("bin.1", list(d.concentrations))]
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _warn_if_medium_cannot_feed(orgs, d)  # must not raise
