"""Integration test for the mechanistic recurrent-CDI / FMT scenario.

This exercises the whole ecology stack at once (antibiotic PK, spore survival,
bile-acid transformation + germination gating, SCFA/pH inhibition, nutrient
competition) and asserts the emergent clinical contrast: without an FMT the
spore reservoir germinates and the infection recurs; with an FMT the restored
community keeps the spores dormant and clears the pathogen.
"""

import pytest

from muode.scenarios import BAI, CDIFF, COMPETITOR, cdi_scenario


@pytest.fixture(scope="module")
def runs():
    return cdi_scenario(fmt=False), cdi_scenario(fmt=True)


def test_antibiotic_clears_vegetative_in_both_arms(runs):
    rec, fmt = runs
    # the drug knocks vegetative C. difficile right down early in both arms
    assert rec.biomass[CDIFF].min() < 0.01
    assert fmt.biomass[CDIFF].min() < 0.01


def test_recurrence_without_fmt(runs):
    rec, _ = runs
    # spores survive the antibiotic and regrow: final >> the knocked-down minimum
    assert rec.biomass[CDIFF].iloc[-1] > 0.5
    assert rec.biomass[CDIFF].iloc[-1] > 50 * rec.biomass[CDIFF].min()
    # no donor guild -> no secondary bile acids -> germination stays permissive
    assert rec.metabolites["dca_e"].iloc[-1] < 1e-6
    assert rec.environment["germination_signal"].iloc[-1] > 0.5


def test_fmt_clears_infection(runs):
    _, fmt = runs
    # C. difficile does not recover after the transplant
    assert fmt.biomass[CDIFF].iloc[-1] < 0.05
    # donor community engrafted
    assert fmt.biomass[COMPETITOR].iloc[-1] > 0.5
    assert fmt.biomass[BAI].iloc[-1] > 0.5
    # mechanisms fired: secondary bile acids made, germination shut off, pH dropped
    assert fmt.metabolites["dca_e"].iloc[-1] > 1.0
    assert fmt.environment["germination_signal"].iloc[-1] < 0.05
    assert fmt.environment["pH"].iloc[-1] < fmt.environment["pH"].iloc[0]


def test_fmt_beats_recurrence(runs):
    rec, fmt = runs
    assert rec.biomass[CDIFF].iloc[-1] > 10 * (fmt.biomass[CDIFF].iloc[-1] + 1e-9)
