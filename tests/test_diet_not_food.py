"""Some rows in a published medium are not food, and one of them was doubling growth.

Three failure modes, all present in `western_gut` at the source table's fill value,
all invisible in a CSV of a hundred plausible-looking metabolite ids:

* **A microbial product supplied as a nutrient.**  Indole and H2S are what a gut
  community MAKES -- from tryptophan and cysteine, which this medium already
  supplies.  Handing them to the community is the acetate/lactate trap: you cannot
  predict a metabolite you were given, and H2S in particular is a headline output
  (colonocyte damage, IBD), so being *given* it destroys the finding.

* **An intracellular intermediate supplied as a nutrient.**  No food delivers free
  chorismate to a colon, and cells detoxify formaldehyde rather than eat it.

* **A real nutrient at a thousand times its real amount.**  Nitrite's bound of 0.1
  mmol/gDW/h is 2.2 GRAMS a day; humans eat 1.2 mg.  Nothing about the row looked
  wrong, and it was the most influential number in the file -- relaxing it 10x
  doubled growth.

These tests pin each one, because the next person to regenerate the medium from the
source will get all three back unless something objects.
"""

from __future__ import annotations

import sys
from pathlib import Path

import cobra
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "diets"))

from derive_western_gut import MISLABELLED_AS_FOOD, PRODUCTS, SKIPPED  # noqa: E402

from muode.diet import load_diet  # noqa: E402
from muode.qc import nutrient_sensitivity  # noqa: E402


@pytest.fixture(scope="module")
def diet():
    return load_diet("western_gut")


@pytest.mark.parametrize("met", MISLABELLED_AS_FOOD)
def test_a_microbial_catabolite_must_be_predicted_not_supplied(diet, met):
    """Starts at zero, no uptake ceiling: somebody has to MAKE it before anybody eats it.

    Uncapped matters as much as zero.  A cross-feeder consuming indole is eating what
    another species excreted, not what arrived in food, so its uptake must be governed
    by its own kinetics rather than by a dietary bound that has no meaning here.
    """
    assert met in PRODUCTS
    assert met in MISLABELLED_AS_FOOD, "these three arrived as *nutrients* -- that is the bug"
    assert diet.initial_concentration(met) == 0.0, f"{met} must not seed the vessel"
    assert diet.uptake_limit(met) is None, f"{met} must be uncapped -- cross-feeding is not dietary"
    assert diet.provenance(met) == "product"


def test_the_substrates_of_those_products_are_still_supplied(diet):
    """The community can only make them if the medium feeds the pathway.

    Moving a product out of the medium is only honest if the precursor is still there:
    tryptophanase needs tryptophan, cysteine desulfhydrase needs cysteine, and
    dissimilatory sulfate reduction needs sulfate.  Otherwise we have not fixed the
    trap, we have just deleted the metabolite.
    """
    for precursor in ("trp__L_e", "cys__L_e", "so4_e"):
        assert diet.initial_concentration(precursor) > 0, (
            f"{precursor} is the substrate for a product we just stopped supplying"
        )


@pytest.mark.parametrize("met", ["chor_e", "fald_e"])
def test_a_non_food_is_not_in_the_medium_at_all(diet, met):
    """Not zeroed -- absent.  And the reason is written down, not just the deletion."""
    assert met in SKIPPED and SKIPPED[met].strip()
    assert met not in diet.metabolites()


def test_nitrite_is_bounded_by_what_people_eat(diet):
    """The quantitative offender: a plausible row with an implausible number.

    0.1 mmol/gDW/h * 20 gDW/L * 24 h = 48 mmol/day = 2.2 g of nitrite a day.  Measured
    intake is 1.2 mg/day (2.6 on a cured-meat-heavy diet).  The medium was ~1,000x
    over, and it mattered more than any other row in the file.
    """
    bound = diet.uptake_limit("no2_e")
    assert diet.provenance("no2_e") == "dri"
    assert bound == pytest.approx(6.8e-5, rel=0.05)

    grams_per_day = bound * 20.0 * 24.0 * 46.01 / 1000.0
    assert grams_per_day < 0.005, f"the medium supplies {grams_per_day:.3g} g/day of nitrite"


@pytest.mark.slow
def test_nitrite_no_longer_dominates_the_prediction(probe_sensitivity):
    """It was `load_bearing` and relaxing it 10x DOUBLED growth (0.052 -> 0.101/h).

    At a real dietary bound iJO1366 does not bother with it at all.  That is the whole
    argument for `nutrient_sensitivity` in one row: the number was carrying the
    prediction, and nobody had chosen it.

    On the probe diet, because E. coli cannot live on the real western_gut once the
    phantom carbon is gone (conftest explains why that is the right outcome).
    """
    r = probe_sensitivity
    assert r["verdict"] == "ok"
    assert not r["nutrients"]["EX_no2_e"]["limiting"], (
        "nitrite is limiting growth again -- check its bound before trusting any result"
    )
