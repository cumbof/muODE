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
    """
    assert set(invented()) == {
        "bile.km_bai",     # community-level Km for 7a-dehydroxylation: not measured
        "bile.km_bsh",     # community-level Km for BSH: not measured
        "bile.vmax_bai",   # per-gram-of-community turnover: nobody has measured it
        "bile.vmax_bsh",   # ditto
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


@pytest.mark.slow
def test_the_fmt_result_does_NOT_come_from_the_bile_mechanism_it_advertises():
    """THE FINDING, and the reason this module exists.

    muode/bile.py opens by explaining that bile acids are "the second major arm of
    colonization resistance against C. difficile", that FMT works by restoring the rare
    bai guild, and that this is "the mechanistic basis of recurrence, and of why FMT
    cures it" (Buffie 2015 Nature; Theriot 2014 Nat Commun).  All of that is real
    biology and it is properly cited.

    **The simulation does not do it.**

    Strip the bile layers out of the rCDI scenario entirely and the FMT arm still clears
    the pathogen.  Strip the pH layer out too and it still clears it.  The clearance is
    nutrient competition: the injected competitor out-eats C. difficile for proline and
    glycine.  And the yields that decide THAT race carry this comment in rcdi.py:

        "Tuned so that fast guilds reach steady state within a few simulated days and
         the pathogen can bloom in an empty lumen but is out-competed in a restored one."

    i.e. the competitive outcome was tuned in.  So the scenario is circular -- tuned to
    produce its conclusion -- and it reaches that conclusion through a mechanism it does
    not advertise, while the mechanism it DOES advertise is decorative.

    That is not a bug in the code; every layer does what it says.  It is a bug in the
    CLAIM, and it is exactly the kind of thing that survives peer review and should not.

    This test pins the defect so it cannot be forgotten.  When the bile mechanism is
    made load-bearing -- with a real bile pool in the diet, at physiological
    concentrations (bile reaches the caecum at ~2 mM), and rates that are measured
    rather than invented -- this test SHOULD fail, and that failure is the fix.
    """
    from muode.antibiotic import Antibiotic
    from muode.bile import BileAcidInhibition, BileAcidTransform
    from muode.dfba import DynamicFBA
    from muode.ecology import EcologyModel
    from muode.inject import Injection
    import sys
    from pathlib import Path

    from muode.lifecycle import SporeForming
    from muode.ph import WeakAcidInhibition

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "fmt_cdiff"))
    from scenario import (  # noqa: E402
        BAI,
        CDIFF,
        COMPETITOR,
        cdi_community,
        cdi_diet,
        cdi_kinetics,
    )

    def run(with_bile: bool) -> float:
        layers = [WeakAcidInhibition(buffer_capacity=25.0, default_ki=12.0,
                                     ki={CDIFF: 3.0})]
        if with_bile:
            layers += [BileAcidTransform(bai_producers={BAI}, vmax_bai=3.0),
                       BileAcidInhibition(targets={CDIFF}, ki=0.05)]
        layers += [
            SporeForming(species={CDIFF}, initial_spores={CDIFF: 0.05},
                         k_germination=0.6, k_sporulation=0.8, mu_stress=0.15),
            Antibiotic(susceptible={CDIFF},
                       dose_times=tuple(float(t) for t in range(0, 10, 2)),
                       dose=3.0, half_life=2.0, emax=5.0, ec50=0.4),
        ]
        fmt = [Injection.from_abundances(12.0, {COMPETITOR: 0.6, BAI: 0.4},
                                         total_biomass=0.05, name="FMT")]
        res = DynamicFBA(t_end=96.0, dt=0.05).run(
            cdi_community(), cdi_diet(), cdi_kinetics(),
            injections=fmt, ecology=EcologyModel(layers),
        )
        return float(res.biomass[CDIFF].iloc[-1])

    with_bile = run(True)
    without_bile = run(False)

    assert with_bile < 0.1, "the FMT arm should clear the pathogen (it does)"
    assert without_bile < 0.1, (
        "REMOVING the entire bile mechanism should ALSO clear it -- and it does. "
        "If this now fails, the bile layer has become load-bearing, which is the "
        "outcome we want: delete this assertion and celebrate."
    )
