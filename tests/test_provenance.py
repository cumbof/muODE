"""Where did that number come from -- and does the answer depend on it?

The medium already answers the first question: every diet row carries a `source`
column, and that column is what let an audit discover that vitamin B12 was the largest
carbon source in `western_gut` and that nine nucleoside rows were feeding E. coli an
RNA diet nobody eats.

The kinetic and ecology parameters had no such column, and they needed one more, because
the rCDI/FMT prediction -- the headline result -- rests on them entirely.  What the audit
found there is worse than "unjustified", and this module pins all of it.
"""

from __future__ import annotations

import pytest

from muode.bile import KI_BILE, VMAX_BAI, VMAX_BSH
from muode.lifecycle import KI_INHIBITOR, KM_GERMINANT
from muode.provenance import (
    Evidence,
    Parameter,
    invented,
    registry,
    report,
    sensitivity,
    unpublishable,
)


def test_a_parameter_must_say_where_it_came_from():
    with pytest.raises(ValueError, match="where its value came from"):
        Parameter(name="x", value=1.0, units="mM", evidence=Evidence.ASSUMED, why="  ")


def test_a_measurement_without_a_citation_is_just_an_assertion():
    """The one rule that makes MEASURED mean anything."""
    with pytest.raises(ValueError, match="cites nothing"):
        Parameter(name="x", value=1.0, units="mM", evidence=Evidence.MEASURED,
                  why="somebody measured it once")


def test_every_registered_parameter_declares_its_evidence():
    for key, p in registry().items():
        assert p.why.strip(), f"{key} has no provenance"
        assert p.units, f"{key} has no units"
        if p.evidence is Evidence.MEASURED:
            assert p.citation.strip(), f"{key} claims MEASURED but cites nothing"


def test_the_bile_parameters_that_could_be_looked_up_HAVE_been(monkeypatch):
    """Two of them were not merely unjustified -- they were WRONG, by 1-2 orders of
    magnitude, against numbers that have been in the literature for over a decade.

        km_germinant  0.1 mM  ->  15.9 mM   (159x)  taurocholate germination EC50
        ki_inhibitor  0.02 mM ->  0.5 mM    (25x)   deoxycholate growth inhibition

    And both errors pushed the SAME WAY: spores germinating far too readily, secondary
    bile acids suppressing far too strongly.  That is the entire rCDI mechanism, biased
    toward the answer we were hoping to see.  This is what an evidence class is for.
    """
    assert KM_GERMINANT.value == pytest.approx(15.9)
    assert KM_GERMINANT.evidence is Evidence.MEASURED
    assert "Ramirez" in KM_GERMINANT.citation

    assert KI_INHIBITOR.value == pytest.approx(0.5)
    assert KI_INHIBITOR.evidence is Evidence.DERIVED
    assert KI_BILE.value == pytest.approx(0.5)
    assert "Usui" in KI_BILE.citation


def test_the_invented_parameters_are_named_and_the_list_only_shrinks():
    """A pinned inventory of what nobody chose.

    This list is allowed to get SHORTER (someone did the reading) and must never get
    longer without a deliberate decision.  A new INVENTED parameter appearing here is a
    new way for a result to be indefensible.

    lifecycle.k_sporulation and lifecycle.k_germination joined the list in the rebuild,
    and that is the registry doing its job rather than a regression.  They were ALWAYS
    invented; they were merely unregistered, so they sat in the layer's defaults looking
    like facts, and the FMT scenario then overrode them (0.8 and 0.6) behind a docstring
    claiming it passed no overrides.  Registering an existing invention does not create
    a new one -- it stops it hiding.  The list got longer because the code got honest.
    """
    assert set(invented()) == {
        "bile.km_bai",     # community-level Km for 7a-dehydroxylation: not measured
        "bile.km_bsh",     # community-level Km for BSH: not measured
        "bile.vmax_bai",   # per-gram-of-community turnover: nobody has measured it
        "bile.vmax_bsh",   # ditto
        # published work reports sporulation FREQUENCIES per generation in vitro, which
        # is not a per-hour rate in a colon; the germinant's affinity is measured
        # (km_germinant) but the rate it saturates at is not
        "lifecycle.k_sporulation",
        "lifecycle.k_germination",
    }
    assert unpublishable() == invented()


def test_the_report_refuses_to_let_an_invented_parameter_pass_quietly():
    text = report()
    assert "INVENTED" in text and "nobody chose these numbers" in text
    assert "bile.vmax_bsh" in text
    # and it must say what to DO about it, not merely tut
    assert "sensitivity" in text


def test_sensitivity_separates_a_conclusion_that_survives_from_one_that_does_not():
    """The move that makes an invented parameter survivable -- or fatal.

    You cannot fix twenty numbers by reading twenty papers.  You CAN find out which of
    them your result actually turns on.  A conclusion that does not move when a
    parameter sweeps across the range you cannot rule out is *robust to not knowing it*,
    and the sweep is the evidence you cite in place of the measurement you lack.
    """
    robust = sensitivity(lambda v: 10.0 + 0.001 * v, "x", [0.1, 1.0, 10.0, 100.0])
    assert robust.verdict == "robust" and robust.robust

    fragile = sensitivity(lambda v: v * 5.0, "x", [0.1, 1.0, 10.0])
    assert fragile.verdict == "sensitive" and not fragile.robust
    assert fragile.spread > 0.9


def test_a_conclusion_that_moves_with_an_invented_number_is_not_a_conclusion():
    """The rule, asserted as a rule.  If the sweep says `sensitive` and the parameter
    says `INVENTED`, no amount of prose rescues the claim -- only a measurement does.
    """
    p = Parameter(name="made_up", value=1.0, units="1/h", evidence=Evidence.INVENTED,
                  why="nobody chose this")
    s = sensitivity(lambda v: v * 3.0, "made_up", [0.1, 1.0, 10.0])
    assert not p.publishable
    assert s.verdict == "sensitive"
    # both true at once == the result must not be reported.  There is no third option.


def test_the_bile_decomposition_is_now_a_genome_scale_STUDY_not_a_toy_unit_test():
    """THE FINDING has moved to the real scenario -- deliberately, and here is the trail.

    On the *toy* rCDI models this module used to prove that the FMT arm clears the
    pathogen even with the bile mechanism removed entirely -- i.e. that clearance was
    nutrient competition and the advertised bile mechanism (Buffie 2015; Theriot 2014)
    was decorative.  But that finding rested on hand-tuned yields whose own comment
    admitted the competitive outcome was tuned in, so it was a toy result about a toy.

    The toy scenario is gone (examples/fmt_cdiff is now genome-scale only).  The finding
    is re-asked properly by ``examples/fmt_cdiff/run.py``, which ablates the bile and pH
    layers on the REAL gapseq community and reports whether clearance survives -- now on
    stoichiometry, not dials.  That is a research result (a ~30 min/arm workstation run),
    not a code invariant, so it is not pinned as a unit test.

    What IS still a code invariant -- that the decomposition machinery works, i.e. that
    ablating a layer actually removes it -- is pinned in
    tests/test_fmt_scenario.py::test_ablation_removes_exactly_the_named_layers.  This
    test just guards the trail so the finding cannot be quietly lost in the refactor.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "fmt_cdiff"))
    import run as fmt_run  # noqa: E402

    # the study still runs the ablation arms that decompose the mechanism
    assert "fmt_full" in fmt_run.ARMS
    assert "fmt_competition" in fmt_run.ARMS      # bile AND pH removed -> pure competition
    assert fmt_run.ARMS["fmt_competition"] == (True, "bile ph")

    # ...AND it contains an arm where nothing is done at all.  The first genome-scale
    # run dosed vancomycin in every arm, so its "control" was drug-without-FMT; the drug
    # kills at ~4.5/h against a pathogen growing at 0.036/h, so every arm reported
    # CLEARED and the study measured the drug rather than the community.  A CDI study
    # with no untreated arm cannot detect colonization resistance, only pharmacology.
    assert "untreated" in fmt_run.ARMS
    assert fmt_run.ARMS["untreated"] == (False, "abx")
    assert "abx_only" in fmt_run.ARMS
    assert "untreated" in fmt_run.CORE_ARMS
