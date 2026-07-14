"""A vitamin must not be a food.

`western_gut` shipped with every vitamin bounded at 0.1 mmol/gDW/h -- the fill value
of the table it came from, and, for biotin and cobalamin, a number an earlier version
of this script put there under a comment claiming they were "catalytic, carrying
neither carbon nor nitrogen".

Cobalamin is C72H100CoN18O17P.  At 0.1 mmol/gDW/h it was contributing 7.2 mmol
C/gDW/h -- **the single largest carbon source in the medium**, ahead of the starch,
the fibre and the sugars.  Twelve vitamin rows between them supplied 18 of the
medium's 51.8 mmol C/gDW/h.  FBA does not know a vitamin from a sugar; it knows a
carbon skeleton and an uptake bound, and it will eat whatever the bound allows.

The defence is dimensional, not vigilant: bound a vitamin by its *dietary reference
intake* (micrograms to milligrams a day) and its carbon vanishes on its own.  These
tests pin the arithmetic that does it, and the property that makes it right.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import cobra
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "diets"))

from derive_western_gut import DRI, SUPPLEMENT, flux_from_intake  # noqa: E402

from muode.diet import load_diet  # noqa: E402


@pytest.fixture(scope="module")
def diet():
    return load_diet("western_gut")


@pytest.fixture(scope="module")
def carbon():
    """mmol C per mmol, for every extracellular metabolite iJO1366 knows a formula for."""
    model = cobra.io.load_model("iJO1366")
    out = {}
    for m in model.metabolites:
        if m.compartment != "e" or not m.formula:
            continue
        hits = re.findall(r"C(\d*)(?![a-z])", m.formula)
        out[m.id] = sum(int(h) if h else 1 for h in hits) if hits else 0
    return out


def test_the_intake_conversion_reproduces_the_rows_already_published():
    """The 1000 that was missing from the docstring.

    ``g / (g/mol)`` is MOLES.  The formula written in this script's comments divided
    by residence time and biomass and called the result millimoles, which is short by
    three orders of magnitude.  The published bounds were right -- they were computed
    with the factor -- but the formula printed beside them was not, and the first
    person to derive a new row *from that formula* got a bound 1000x too small.

    So the conversion is a function now, and this is the test that it agrees with the
    numbers already in the file.
    """
    assert flux_from_intake(3.5, 60.06) == pytest.approx(0.12, rel=0.02)    # urea
    assert flux_from_intake(2.0, 528.5) == pytest.approx(0.0079, rel=0.02)  # xylan4
    assert flux_from_intake(2.0, 1057.0) == pytest.approx(0.0039, rel=0.02)  # xylan8


def test_a_vitamin_is_never_a_major_carbon_source(diet, carbon):
    """The property the DRI bound buys, stated as a property rather than a number.

    No vitamin may contribute as much carbon as the medium's real foods.  The yardstick
    is starch, which is what a colon actually runs on: 1e-4 mmol/gDW/h of a 1200-mer of
    glucose is 1e-4 * 1200 * 6 = 0.72 mmol C/gDW/h.  (Note the 6.  iJO1366 has no
    formula for starch1200, so the carbon has to be counted by hand, and counting the
    glucose units and calling them carbon undercounts it six-fold.)

    Cobalamin at the old fill value was 7.2 mmol C/gDW/h -- **ten times the starch**.

    Choline is deliberately the closest of the survivors, and that is not a defect: its
    AI is 550 mg/day, a thousand-fold more than the microgram B vitamins, and it is a
    genuine colonic substrate (it is the precursor gut bacteria turn into trimethylamine).
    It should be a real if minor carbon source, and it is -- ~8% of the starch.
    """
    starch_c = diet.uptake_limit("starch1200_e") * 1200 * 6
    assert starch_c == pytest.approx(0.72, rel=0.05), "the yardstick moved; re-read this test"

    for met in DRI:
        c = diet.uptake_limit(met) * carbon.get(met, 0)
        assert c < starch_c / 10, (
            f"{met} supplies {c:.3g} mmol C/gDW/h against starch's {starch_c:.3g} -- "
            "a vitamin bounded like a food is not a conservative choice, it is a gift"
        )


def test_cobalamin_specifically_is_no_longer_the_biggest_carbon_source(diet, carbon):
    """The headline case, pinned by name because it is the one that got through."""
    b12 = diet.uptake_limit("adocbl_e") * carbon["adocbl_e"]
    biggest = max(
        diet.uptake_limit(m) * carbon.get(m, 0)
        for m in diet.metabolites()
        if m != "h2o_e" and diet.uptake_limit(m) is not None
    )
    assert b12 < biggest / 1000, (
        f"B12 supplies {b12:.3g} mmol C/gDW/h; the medium's largest source is {biggest:.3g}. "
        "It was once 7.2 and it WAS the largest."
    )


def test_vitamins_are_bounded_by_intake_and_minerals_by_the_median(diet):
    """The two dicts mean different things, and conflating them is what caused the bug.

    SUPPLEMENT is for MINERALS: catalytic, no carbon, no nitrogen, and parked at the
    median because their exact bound provably cannot move a prediction (see
    test_nutrient_sensitivity.py).  DRI is for organic cofactors, where the bound is
    the whole point.  An entry in both would be a vitamin at a mineral's bound, which
    is precisely how cobalamin ended up outweighing the starch.
    """
    assert not (set(SUPPLEMENT) & set(DRI))
    for met in DRI:
        assert diet.provenance(met) == "dri", f"{met} should carry its intake provenance"
    for met in SUPPLEMENT:
        assert diet.provenance(met) in ("supplement", "source_default")


def test_fixing_the_vitamins_did_not_change_what_a_prototroph_predicts(diet):
    """A 300x cut to the vitamin carbon must not move E. coli, and must not break it.

    iJO1366 synthesises its own B vitamins, so `nutrient_sensitivity` calls every one
    of these rows `unused` for it: the carbon they were donating was pure artefact,
    and taking it away changes nothing it predicts.  (An auxotroph is a different
    story, and *should* be -- it must now get its cobalamin by cross-feeding from a
    producer, which is what actually happens in a colon.)
    """
    from muode.qc import nutrient_sensitivity

    r = nutrient_sensitivity(cobra.io.load_model("iJO1366"), diet)
    assert r["verdict"] == "ok"
    assert r["growth"] == pytest.approx(0.0516, abs=5e-3), (
        "E. coli's growth on western_gut should be untouched by the vitamin fix"
    )
    for met in ("adocbl_e", "btn_e", "thm_e", "fol_e"):
        ex = f"EX_{met}"
        if ex in r["nutrients"]:
            assert not r["nutrients"][ex]["limiting"], (
                f"{met} is limiting growth -- a vitamin bound should never do that"
            )
