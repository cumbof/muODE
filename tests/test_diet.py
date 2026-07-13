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
    for essential in ("pi_e", "so4_e", "k_e", "mg2_e", "ca2_e", "fe2_e", "h2o_e"):
        assert d.initial_concentration(essential) > 0, f"{essential} missing from the diet"
    assert d.initial_concentration("btn_e") > 0                # vitamins
    # nitrogen -- from amino acids, NOT ammonium (see below)
    assert d.initial_concentration("ala__L_e") > 0
    # carbon -- from what actually reaches the colon
    carbon = ("starch1200_e", "amylose300_e", "lcts_e", "malt_e", "sucr_e", "fru_e")
    assert any(d.initial_concentration(c) > 0 for c in carbon)


def test_the_colonic_medium_has_no_free_glucose_and_no_ammonium():
    """Both absences are correct biology, and both surprise people.

    Glucose is absorbed in the small intestine: what reaches the colon is starch,
    amylose, pullulan, lactose and fibre.  A model that cannot degrade a
    polysaccharide therefore cannot grow on this medium *alone* -- in the gut it
    lives on sugars cross-fed by primary degraders, and in muODE it must do the
    same.  Nitrogen likewise arrives as amino acids, not as NH4+.

    If a future edit reintroduces either, it is smuggling the small intestine into
    the colon and every cross-feeding result becomes suspect.
    """
    d = load_preset("western_gut")
    assert d.initial_concentration("glc__D_e") == 0.0
    assert d.initial_concentration("nh4_e") == 0.0
    assert sum(1 for m in d.metabolites() if m.endswith("__L_e")) >= 15   # the N source


def test_western_gut_is_anaerobic():
    d = load_preset("western_gut")
    assert d.initial_concentration("o2_e") == 0.0
    assert d.influx_rate("o2_e") == 0.0
    assert d.uptake_limit("o2_e") is None


def test_fermentation_products_are_cross_fed_not_supplied():
    """SCFAs must be produced by the community, not handed to it.

    The published source medium *does* supply acetate, formate and H2 (a dietary
    flux of 0.1 mmol/gDW/h).  Carrying that through would start the vessel at
    ~48 mM acetate -- about the colonic steady-state the community is supposed to
    produce -- and make SCFA prediction unfalsifiable.  Same trap as lactate in
    DM38.  The derive script forces them to zero.
    """
    d = load_preset("western_gut")
    for product in ("ac_e", "but_e", "ppa_e", "lac__D_e", "succ_e", "for_e", "h2_e"):
        assert d.initial_concentration(product) == 0.0
        assert d.influx_rate(product) == 0.0
        # ...and uncapped: a butyrate producer eats acetate made by ANOTHER SPECIES,
        # so throttling it at the *dietary* flux would throttle cross-feeding itself.
        assert d.uptake_limit(product) is None, (
            f"{product} is capped at the dietary flux, which would cap cross-feeding"
        )


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
