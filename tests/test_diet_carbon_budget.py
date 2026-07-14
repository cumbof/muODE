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
def carbon(ijo):
    """mmol C per mmol, for every extracellular metabolite iJO1366 knows a formula for."""
    out = {}
    for m in ijo.metabolites:
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


def test_the_nucleosides_are_bounded_by_dietary_nucleic_acid(diet, carbon):
    """The rows that were feeding E. coli an RNA diet.

    Nine nucleosides plus AMP, each at 0.1 mmol/gDW/h -- around 9 mmol C/gDW/h once you
    count the ~10 carbons in a nucleoside, which after the vitamin fix made them the
    largest carbon source left in the medium.  Humans eat 0.1-1 g/day of nucleic acid.
    """
    nucleosides = ["adn_e", "gsn_e", "cytd_e", "uri_e", "amp_e",
                   "dad_2_e", "dgsn_e", "dcyt_e", "thymd_e"]
    total_c = sum(diet.uptake_limit(m) * carbon.get(m, 0) for m in nucleosides)

    for met in nucleosides:
        assert diet.provenance(met) == "derived"
    assert total_c < 0.05, (
        f"the nucleosides supply {total_c:.3g} mmol C/gDW/h; at the old fill value they "
        "supplied ~9, and that carbon was the only thing keeping iJO1366 alive"
    )

    # and the mass balance: the whole medium's nucleoside flux must not exceed what a
    # person eats.  0.5 g/day was the assumption; allow a factor of 2 either way.
    g_per_day = sum(
        diet.uptake_limit(m) * 20.0 * 24.0 * mw / 1000.0
        for m, mw in [("adn_e", 267.24), ("gsn_e", 283.24), ("cytd_e", 243.22),
                      ("uri_e", 244.20), ("amp_e", 347.22), ("dad_2_e", 251.24),
                      ("dgsn_e", 267.24), ("dcyt_e", 227.22), ("thymd_e", 242.23)]
    )
    assert 0.2 < g_per_day < 1.5, (
        f"the medium supplies {g_per_day:.2f} g/day of nucleoside; people eat 0.1-1 g/day "
        "of nucleic acid"
    )


def test_mucin_feeds_the_mucolytic_specialists_at_a_host_derived_rate(diet):
    """2-3 g/day of mucin enters the large bowel; MUC2 is ~80% O-glycan by weight.

    This is the substrate Akkermansia and many Bacteroides actually live on, and its
    amino sugars were sitting at the fill value like everything else.  The derived
    bounds have to be small -- but not zero, because mucin is a genuine and important
    colonic carbon source, and a gut medium that omits it is missing a whole guild's
    niche.
    """
    assert diet.provenance("acgam_e") == "derived"
    assert diet.provenance("acnam_e") == "derived"
    assert 1e-3 < diet.uptake_limit("acgam_e") < 2e-2, "GlcNAc from ~2 g/day mucin glycan"
    assert 1e-4 < diet.uptake_limit("acnam_e") < 1e-2, "sialic acid, ~8 mol% of the glycan"

    # ManNAc and glucosamine are what bacteria MAKE from those two; they are not food.
    for downstream in ("acmana_e", "gam_e"):
        assert diet.uptake_limit(downstream) is None
        assert diet.initial_concentration(downstream) == 0.0


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


@pytest.mark.slow
def test_ecoli_was_living_on_the_fill_value_and_now_starves(diet, ijo):
    """The single most consequential thing the medium audit turned up.

    iJO1366 grew at 0.0516/h on western_gut, and **all of it came from nine nucleoside
    rows sitting at the source table's fill value**.  Restore those rows to 0.1 and it
    lives; bound them by real dietary nucleic acid (0.1-1 g/day, so ~200x less) and it
    cannot make its maintenance ATP at all:

        max ATP available to it:     1.97 mmol/gDW/h
        its ATPM demand:             3.15 mmol/gDW/h   -> the LP is INFEASIBLE

    E. coli was eating an RNA diet that nobody had chosen, and it was the only thing
    keeping it alive on a "western gut" medium.

    The starvation is correct biology, which is why this test asserts it rather than
    working around it.  E. coli ferments sugar; a colonic medium's sugar is locked in
    fibre polymers it cannot open; in a gut it survives by cross-feeding on what the
    fibre degraders release.  A monoculture on colonic contents should struggle.

    If this test ever starts failing because iJO1366 grows again, something has put
    phantom carbon back into the medium.  Find out what before believing the result.
    """
    from muode.qc import diagnose_no_growth, nutrient_sensitivity

    r = nutrient_sensitivity(ijo, diet)
    assert r["verdict"] == "no_growth", (
        f"iJO1366 grows at {r.get('growth', 0):.4f}/h on western_gut again -- it should "
        "not be able to.  Check whether a fill value has crept back in."
    )

    # and the framework must SAY it is a medium gap rather than a fact about E. coli
    d = diagnose_no_growth(ijo, diet)
    assert d["verdict"] == "medium_gap"
    assert "EX_glc__D_e" in d["rescued_by_any_of"], (
        "free sugar should be on the rescue menu -- that is the whole diagnosis: it is "
        "short of fermentable carbon, not broken"
    )


@pytest.mark.slow
def test_the_vitamin_carbon_was_pure_artefact(ijo, western_gut_probe, probe_sensitivity):
    """The causal test: delete the vitamins entirely and a prototroph does not notice.

    iJO1366 synthesises its own B vitamins.  So if the 18 mmol C/gDW/h those rows used
    to donate had been doing any real work, removing them would show up in growth --
    and it does not.  The carbon was artefact, and the medium is not quietly running on
    it.  (An auxotroph is a different story, and *should* be: it must now get its
    cobalamin by cross-feeding from a producer, which is what happens in a colon.)

    Deliberately NOT a pin on an absolute growth rate.  Other rows in this medium move
    that number for their own good reasons, and a test that conflates them would fail
    for the wrong cause and get "fixed" by editing the constant.

    Run on the probe diet (see conftest): E. coli cannot live on the real western_gut
    at all, so there would otherwise be no growth to compare.
    """
    from muode.diet import Diet
    from muode.qc import nutrient_sensitivity

    with_vitamins = probe_sensitivity
    assert with_vitamins["verdict"] == "ok"

    without = Diet(
        concentrations={m: c for m, c in western_gut_probe.concentrations.items()
                        if m not in DRI},
        influx={m: v for m, v in western_gut_probe.influx.items() if m not in DRI},
        max_uptake={m: v for m, v in western_gut_probe.max_uptake.items() if m not in DRI},
        name="probe minus every vitamin",
    )
    stripped = nutrient_sensitivity(ijo, without)
    assert stripped["verdict"] == "ok", "a prototroph must not need the vitamin rows at all"
    assert stripped["growth"] == pytest.approx(with_vitamins["growth"], rel=1e-3), (
        "removing every vitamin changed growth -- so the medium IS running on their "
        "carbon, and the bound on those rows is doing work it should not be doing"
    )

    for met in ("adocbl_e", "btn_e", "thm_e", "fol_e"):
        ex = f"EX_{met}"
        if ex in with_vitamins["nutrients"]:
            assert not with_vitamins["nutrients"][ex]["limiting"], (
                f"{met} is limiting growth -- a vitamin bound should never do that"
            )
