"""Why DM38 carries no ``max_uptake`` -- and what would break if it did.

`western_gut` needs a dietary flux ceiling: it is a flow system, the dietary flux
IS the measured quantity, and without the ceiling a GEM imports ~100 nutrients at a
uniform Vmax and "grows" at 6.5/h.  DM38 is the opposite case, and copying the fix
across would be a mistake:

* DM38 is a **sealed batch culture with published concentrations**.  25 mM glucose
  is genuinely available.  There is no dietary flux to bound, so a ``max_uptake``
  would be a number we invented, not one we measured.
* DM38 is the medium the Clark benchmark **calibrates Vmax against**.  Uptake is
  ``min(Michaelis-Menten, max_uptake)``, so a dietary ceiling clamps uptake *below*
  Vmax, flattens growth-vs-Vmax, and makes Vmax **unidentifiable** -- destroying the
  very quantity the benchmark exists to fit.

So the exposure is real (see the first test: the default Vmax really does produce an
implausible growth rate) but the instrument that fixes it is Vmax, not the diet.
These tests pin that reasoning so nobody -- including a future me -- "fixes" DM38 by
adding bounds to it.
"""

from __future__ import annotations

import cobra
import pytest

from muode.diet import Diet, load_diet
from muode.kinetics import DEFAULT_VMAX, KineticParameters
from muode.media import growth_on_diet
from muode.qc import MAX_PLAUSIBLE_GROWTH

#: DM38 lacks these three, and iJO1366 (E. coli -- NOT one of Clark's 25 strains)
#: needs them.  Supplied here only so a real GEM grows at all and the Vmax effect is
#: visible.  Iron is the load-bearing one AMONG THESE FOUR -- not the only possible
#: rescue; see tests/test_no_growth_diagnosis.py.
IJO_NEEDS = ("co2_e", "fe3_e", "sel_e", "tungs_e")


@pytest.fixture(scope="module")
def model(ijo):
    return ijo


@pytest.fixture(scope="module")
def dm38_for_ijo():
    """DM38, plus what E. coli needs to be alive on it.  Still anaerobic: no o2_e."""
    conc = dict(load_diet("dm38").concentrations)
    for m in IJO_NEEDS:
        conc[m] = 1.0
    assert "o2_e" not in conc, "DM38 is anaerobic; adding oxygen would fake the result"
    return conc


def _growth(model, conc, vmax):
    return growth_on_diet(model, Diet(concentrations=conc, name="dm38"),
                          KineticParameters(default_vmax=vmax))


def test_dm38_carries_no_dietary_ceiling():
    """Deliberate, not an oversight.  The rest of this module is the justification."""
    assert not load_diet("dm38").max_uptake


@pytest.mark.slow
def test_the_default_vmax_really_does_produce_an_impossible_organism(model, dm38_for_ijo):
    """The exposure is real: DM38 + the uniform default Vmax is above the QC ceiling.

    This is the test that must NOT be silenced by raising MAX_PLAUSIBLE_GROWTH.  It
    is here to keep the problem visible until the Clark calibration supplies a
    measured Vmax to replace the default with.
    """
    g = _growth(model, dm38_for_ijo, DEFAULT_VMAX)
    assert g > MAX_PLAUSIBLE_GROWTH, (
        "DM38 at the default Vmax used to grow iJO1366 at ~2.9/h.  If this now "
        "passes, either the default Vmax was fixed (good -- update this test) or the "
        "QC ceiling was raised to hide it (bad)."
    )


@pytest.mark.slow
def test_vmax_is_identifiable_on_dm38_which_is_why_the_diet_must_not_cap_uptake(
    model, dm38_for_ijo
):
    """Growth is strictly monotone in Vmax over the plausible range -- no plateau.

    That monotonicity is what makes the Clark calibration possible at all: a plateau
    would mean the data cannot distinguish one Vmax from another, which is exactly
    the `unidentifiable` verdict muode.calibrate.assess() is built to return.
    """
    grid = [0.5, 1.0, 2.0, 5.0, 10.0]
    growth = [_growth(model, dm38_for_ijo, v) for v in grid]

    assert all(b > a for a, b in zip(growth, growth[1:])), (
        f"growth must increase with Vmax for it to be fittable; got {growth}"
    )
    # and it must MOVE, not creep -- a 20x Vmax range has to span a real growth range
    assert growth[-1] > 10 * growth[0]


@pytest.mark.slow
def test_a_dietary_ceiling_would_make_vmax_unidentifiable(model, dm38_for_ijo):
    """The concrete reason DM38 must not get max_uptake.

    Cap uptake the way western_gut does and growth stops responding to Vmax: every
    value in the grid gives the same answer, so no amount of Clark data could ever
    tell us what Vmax is.  We would have broken the calibration to fix the medium.
    """
    capped = {m: 0.1 for m in dm38_for_ijo}          # a western_gut-style ceiling
    growth = [
        growth_on_diet(model,
                       Diet(concentrations=dm38_for_ijo, max_uptake=capped, name="capped"),
                       KineticParameters(default_vmax=v))
        for v in (1.0, 2.0, 5.0, 10.0)
    ]
    assert max(growth) - min(growth) < 1e-6, (
        f"a dietary cap should flatten growth-vs-Vmax entirely; got {growth}"
    )


@pytest.mark.slow
def test_dm38_supplies_ferrous_iron_only_and_iJO1366_starves_on_it(model, dm38_for_ijo):
    """DM38 has fe2_e and no fe3_e, and iJO1366 needs *cytoplasmic* Fe(III).

    Chemically the medium is right -- Fe(III) is not stable in a reducing anaerobic
    broth -- so this is not a derivation error.  It is a trap: the model is dead on
    arrival and its zero says nothing about E. coli.

    Note the careful claim.  Ferric iron is load-bearing *among the four nutrients
    supplied here*, which is what this test checks.  It is NOT the only thing that
    would rescue the model: nitric oxide does too, via FESD2s, by wrecking the cell's
    own Fe-S clusters to liberate iron.  Naming one nutrient as "the" cause is exactly
    the mistake muode.qc.diagnose_no_growth exists to prevent -- see
    tests/test_no_growth_diagnosis.py, which enumerates the whole equivalence class
    instead of picking a favourite.
    """
    dm38 = load_diet("dm38")
    assert dm38.initial_concentration("fe2_e") > 0
    assert "fe3_e" not in dm38.metabolites()

    without_fe3 = {m: c for m, c in dm38_for_ijo.items() if m != "fe3_e"}
    assert _growth(model, without_fe3, DEFAULT_VMAX) < 1e-6, (
        "among co2/sel/tungs, iron is the one that matters; if this now grows, the "
        "trap has changed shape and the diagnosis in test_no_growth_diagnosis.py "
        "should be re-read rather than this test deleted"
    )
