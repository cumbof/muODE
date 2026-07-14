"""Dietary flux bounds: what the diet supplies, not just what it contains.

The bug these pin: `diet_medium` turned a diet into uptake bounds using only
`Vmax * C / (Km + C)` with a *uniform* default Vmax of 10 mmol/gDW/h.  Every
nutrient in a large medium therefore saturated at the same bound, so a
genome-scale model imported ~100 nutrients at 10 mmol/gDW/h each.  Measured on
iJO1366 that is 1129 mmol/gDW/h of total uptake and a growth rate of 12.8/h --
a three-minute doubling.

The published western gut diet instead specifies a flux bound per metabolite
(median 0.1 mmol/gDW/h).  Same model, same nutrients, growth 0.055/h.
"""

import cobra
import pytest

from muode.diet import Diet, load_diet
from muode.kinetics import KineticParameters
from muode.media import diet_medium, growth_on_diet
from muode.qc import MAX_PLAUSIBLE_GROWTH


def test_the_published_diet_carries_flux_bounds():
    """The median used to be 0.1 -- the source table's *fill value* -- and is now 0.018.

    That drop is the medium audit, not a regression.  0.1 was the number the source
    used for every row its intake data did not itemise, and 56% of the file sat on it.
    Bounding those rows by what people actually eat (vitamins by DRI, nucleosides by
    dietary nucleic acid, nitrite by measured intake) pulled the median down to the
    scale of the medium's *real* foods -- sucrose is 0.0074, xylose 0.018.

    The assertion that matters has not changed: the bounds are far below the uniform
    Vmax of 10 that this file exists to have replaced.
    """
    d = load_diet("western_gut")
    assert d.max_uptake, "western_gut must carry dietary flux bounds"
    bounds = sorted(d.max_uptake.values())
    median = bounds[len(bounds) // 2]
    assert median == pytest.approx(0.018, rel=0.2), (
        f"median bound is {median:g}; it should sit with the medium's real foods "
        "(sucrose 0.0074, xylose 0.018), not on the source's 0.1 fill value"
    )
    assert median * 50 < 10.0, "...which is far below muODE's old uniform Vmax of 10"


def test_uptake_is_the_min_of_kinetics_and_dietary_availability():
    """A cell may be slower than the diet allows; it may never be faster."""
    kin = KineticParameters()
    model = cobra.io.load_model("textbook")

    generous = Diet(concentrations={"glc__D_e": 100.0}, max_uptake={"glc__D_e": 0.1},
                    name="capped")
    stingy = Diet(concentrations={"glc__D_e": 100.0}, max_uptake={"glc__D_e": 1e6},
                  name="uncapped-by-diet")

    # diet is the binding constraint -> the diet's number wins
    assert diet_medium(model, generous, kin)["EX_glc__D_e"] == pytest.approx(0.1)
    # kinetics are the binding constraint -> a huge dietary ceiling changes nothing
    kinetic_rate = kin.michaelis_menten(model.id, "glc__D_e", 100.0)
    assert diet_medium(model, stingy, kin)["EX_glc__D_e"] == pytest.approx(kinetic_rate)
    assert kinetic_rate < 1e6


def test_a_diet_without_bounds_is_unchanged():
    """The cap is additive: no max_uptake, no behaviour change."""
    kin = KineticParameters()
    model = cobra.io.load_model("textbook")
    d = Diet(concentrations={"glc__D_e": 20.0}, name="no-bounds")
    assert d.uptake_limit("glc__D_e") is None
    assert diet_medium(model, d, kin)["EX_glc__D_e"] == pytest.approx(
        kin.michaelis_menten(model.id, "glc__D_e", 20.0)
    )


def test_zero_supply_and_unbounded_uptake_are_not_confused():
    """An empty max_uptake cell means "uncapped", not "not supplied"."""
    d = Diet(concentrations={"ac_e": 0.0, "glc__D_e": 5.0}, max_uptake={"glc__D_e": 0.1},
             name="mixed")
    assert d.uptake_limit("ac_e") is None        # uncapped -- a cross-fed product
    assert d.uptake_limit("glc__D_e") == 0.1     # capped by the diet
    assert d.initial_concentration("ac_e") == 0.0


def test_flux_bounds_survive_a_csv_round_trip(tmp_path):
    d = Diet(concentrations={"glc__D_e": 5.0, "ac_e": 0.0},
             max_uptake={"glc__D_e": 0.1}, name="rt")
    path = tmp_path / "d.csv"
    d.to_csv(path)
    back = Diet.from_csv(path)
    assert back.max_uptake == {"glc__D_e": 0.1}
    assert back.uptake_limit("ac_e") is None     # blank stayed blank, did not become 0.0


def test_with_metabolites_does_not_drop_the_bounds():
    """A regression guard: `max_uptake` was added as a new positional field."""
    d = Diet(concentrations={"glc__D_e": 5.0}, max_uptake={"glc__D_e": 0.1}, name="x")
    extended = d.with_metabolites({"but_e": 0.0})
    assert extended.max_uptake == {"glc__D_e": 0.1}
    assert extended.name == "x"
    assert "but_e" in extended.concentrations


@pytest.mark.slow
def test_the_bounds_keep_a_genome_scale_model_inside_biology():
    """The whole point, measured end to end on a real genome-scale GEM.

    E. coli core has only 20 exchanges, so it could never expose this bug -- it
    physically cannot import 100 nutrients.  It takes a genome-scale model to show
    it, which is why the 89-MAG run found it and the test suite did not.
    """
    model = cobra.io.load_model("iJO1366")
    kin = KineticParameters()
    base = load_diet("western_gut")

    # Supplement what iJO1366 needs and a colonic medium lacks: a sugar it can actually
    # ferment, an N source, its trace ions.  The sugar is load-bearing here -- E. coli
    # cannot live on western_gut at all now that the phantom nucleoside carbon is gone
    # (see tests/conftest.py), because a colon's sugar is locked in fibre polymers it
    # has no enzymes for.  1.0 mmol/gDW/h of glucose clears its maintenance demand; the
    # rest ride at the medium's own median.
    missing = ["nh4_e", "na1_e", "ni2_e", "sel_e", "slnt_e", "tungs_e", "cbl1_e", "co2_e"]
    conc = dict(base.concentrations)
    lim = dict(base.max_uptake)
    for met in missing:
        conc[met], lim[met] = 48.0, 0.1
    conc["glc__D_e"], lim["glc__D_e"] = 480.0, 1.0

    capped = Diet(concentrations=conc, max_uptake=lim, name="published-bounds")
    uncapped = Diet(concentrations=conc, name="uniform-Vmax-10")   # the old behaviour

    fast = growth_on_diet(model, uncapped, kin)
    sane = growth_on_diet(model, capped, kin)

    assert fast > 5.0, "the old behaviour really was this broken"
    assert fast > MAX_PLAUSIBLE_GROWTH, "...and the QC ceiling catches it"
    assert 0.0 < sane < MAX_PLAUSIBLE_GROWTH, "the published bounds restore biology"
    assert sane < fast / 50, "the bounds are not a rounding correction"
