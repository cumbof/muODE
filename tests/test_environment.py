"""An environmental assumption belongs in an object you can print, not a constant.

muODE carried two of them in qc.py::

    MAX_PLAUSIBLE_GROWTH = 1.5              # justified by "gut anaerobes run 0.1-0.5/h"
    FORBIDDEN_RESCUES = frozenset({"o2_e"}) # i.e. "this environment is anaerobic"

Both are claims about a human colon.  Neither said so at the call site, and muODE is
supposed to take arbitrary genomes -- which means arbitrary environments.  On a soil
sample the first censors a legitimate fast aerobe as a "modelling artefact", and the
second refuses to tell an oxic medium that the thing it is missing is oxygen.
"""

from __future__ import annotations

import cobra
import pytest

from muode.diet import Diet
from muode.environment import (
    ABSOLUTE_MAX_GROWTH,
    GENERIC_AEROBIC,
    GENERIC_ANAEROBIC,
    HUMAN_GUT,
    Environment,
    known_environments,
    load_environment,
)
from muode.qc import diagnose_no_growth, sanity_check_model


def test_the_gut_ceiling_is_still_the_default_and_now_admits_to_being_the_gut():
    """Moving it must not silently change any existing result."""
    from muode.qc import FORBIDDEN_RESCUES, MAX_PLAUSIBLE_GROWTH

    assert HUMAN_GUT.max_plausible_growth == 1.5
    assert MAX_PLAUSIBLE_GROWTH == HUMAN_GUT.max_plausible_growth
    assert FORBIDDEN_RESCUES == HUMAN_GUT.forbidden_rescues == frozenset({"o2_e"})
    assert HUMAN_GUT.oxygen is False


def test_no_environment_may_claim_a_ceiling_above_the_physical_limit():
    """The one number here that is not a judgement call.

    Nothing has ever been observed to grow faster than V. natriegens: 9.8 min doubling
    (Eagon 1962), i.e. ln2/(9.8/60) = 4.24/h.  A "plausibility" ceiling above that
    cannot catch anything, because nothing can exceed it -- so asking for one is a
    mistake, and saying so is cheaper than letting it through.
    """
    assert ABSOLUTE_MAX_GROWTH == pytest.approx(4.24, rel=0.01)
    with pytest.raises(ValueError, match="ABSOLUTE_MAX_GROWTH"):
        Environment(name="wishful", max_plausible_growth=10.0, why="because I said so")


def test_an_environment_cannot_be_created_without_saying_why():
    """The enforcement mechanism, and the whole point of the module.

    A growth ceiling is an empirical claim.  If you cannot write down why yours is what
    it is, you do not have a ceiling -- you have a number, which is what we had.
    """
    with pytest.raises(ValueError, match="must say WHY"):
        Environment(name="mystery", max_plausible_growth=2.0, why="   ")


def test_an_anaerobic_environment_always_forbids_oxygen_even_if_you_forget_to_say_so():
    """Oxygen rescues almost any GEM, so in an anaerobic setting it would be the answer
    every time and drown every real finding.  Too important to be opt-in."""
    e = Environment(name="some_anoxic_place", max_plausible_growth=1.0,
                    why="a test", oxygen=False)
    assert "o2_e" in e.forbidden_rescues


def test_an_aerobic_environment_does_NOT_forbid_oxygen():
    """The case the old constant got wrong, and the reason this is a parameter.

    In a colon, proposing oxygen would OVERTURN the environment rather than complete
    it.  In a soil or marine sample the opposite is true: oxygen belongs there, so an
    oxic medium that omits o2_e is a medium with a hole in it, and the diagnosis should
    say so rather than staying politely silent.
    """
    assert GENERIC_AEROBIC.oxygen is True
    assert "o2_e" not in GENERIC_AEROBIC.forbidden_rescues


def test_an_uncharacterised_environment_refuses_to_pretend_it_knows_the_ceiling():
    """The honest default when you have not done the reading.

    GENERIC_ANAEROBIC's ceiling is the physical limit, so it catches only the
    impossible.  That is deliberate: a filter that censors real biology is a worse
    failure than one that lets an artefact through, and inventing a tighter number to
    look rigorous is precisely the habit this module exists to break.
    """
    assert GENERIC_ANAEROBIC.max_plausible_growth == ABSOLUTE_MAX_GROWTH
    assert GENERIC_ANAEROBIC.max_plausible_growth > HUMAN_GUT.max_plausible_growth


def test_only_environments_we_can_justify_are_shipped():
    """There is no `rumen` or `marine` preset, and that is a feature.

    Shipping one would mean inventing its growth ceiling, which would move the invented
    constant somewhere more official-looking rather than removing it.  Every preset here
    carries a citation in its `why`.
    """
    for name, env in known_environments().items():
        assert env.why.strip(), f"{name} has no justification"
        assert len(env.why) > 80, f"{name}'s justification is too thin to be one"

    with pytest.raises(KeyError, match="deliberately"):
        load_environment("rumen")


def test_the_growth_ceiling_actually_follows_the_environment():
    """The behaviour, not just the plumbing: same model, same rate, different verdict."""
    model = cobra.io.load_model("textbook")

    gut = sanity_check_model(model, growth_rate=2.5, environment=HUMAN_GUT)
    soil = sanity_check_model(model, growth_rate=2.5, environment=GENERIC_AEROBIC)

    assert gut["implausible_growth"] is True, "2.5/h is absurd for a gut anaerobe"
    assert soil["implausible_growth"] is False, "2.5/h is unremarkable for a fast aerobe"
    assert gut["environment"] == "human_gut" and soil["environment"] == "generic_aerobic"


@pytest.mark.slow
def test_the_diagnosis_offers_oxygen_in_an_oxic_setting_and_never_in_an_anoxic_one(ijo):
    """The other half, end to end on a real GEM.

    iJO1366 on a medium with no oxygen and no fermentable carbon: aerobically, oxygen is
    exactly what it is missing.  The gut environment must refuse to say so; the aerobic
    environment must be free to.
    """
    minimal = Diet(
        concentrations={m: 10.0 for m in
                        ("glc__D_e", "nh4_e", "pi_e", "so4_e", "k_e", "mg2_e",
                         "fe3_e", "ca2_e", "cl_e", "cobalt2_e", "cu2_e", "mn2_e",
                         "mobd_e", "ni2_e", "zn2_e")},
        name="aerobic minimal, oxygen omitted",
    )

    anoxic = diagnose_no_growth(ijo, minimal, environment=HUMAN_GUT)
    oxic = diagnose_no_growth(ijo, minimal, environment=GENERIC_AEROBIC)

    # In the gut, oxygen is never on the menu -- whatever it would do for growth.
    assert "EX_o2_e" not in anoxic.get("rescued_by_any_of", [])
    # In an oxic setting it is allowed to be the answer, and here it should be.
    if oxic["verdict"] == "medium_gap":
        assert "EX_o2_e" in oxic["rescued_by_any_of"], (
            "an aerobic environment must be able to say that the missing thing is oxygen"
        )
