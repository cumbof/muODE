"""The growth ceiling: a GEM must not be allowed to grow faster than any cell.

This exists because of a real run.  89 CarveMe models were reconstructed against
the `western_gut` diet; all 89 passed QC; 25 of them grew faster than 2.5/h and
the fastest hit 6.53/h -- a 6-minute doubling time, quicker than any organism ever
measured.  Nothing in the pipeline objected.  A flat community looks broken and
gets investigated; an absurdly fast one looks like a *result*.
"""

import json

import cobra
import pytest

from muode.qc import MAX_PLAUSIBLE_GROWTH, sanity_check_model, summarize_reconstruction


def _qc_json(path, mag, growth, implausible):
    path.write_text(json.dumps({
        "model": mag,
        "grew_initially": True,
        "reactions_added": [],
        "grows_now": True,
        "growth_on_diet": growth,
        "qc": {
            "n_mass_unbalanced": 0,
            "energy_generating_cycle": False,
            "growth_rate": growth,
            "implausible_growth": implausible,
            "passed": not implausible,
        },
    }))
    return path


def test_the_ceiling_sits_above_every_real_gut_organism():
    """A slow-growing commensal must never be flagged; only the impossible is."""
    assert MAX_PLAUSIBLE_GROWTH > 1.0        # E. coli in rich aerobic medium
    assert MAX_PLAUSIBLE_GROWTH < 2.5        # ...but below a 17-minute doubling


@pytest.mark.parametrize("rate, flagged", [
    (0.21, False),   # E. coli, anaerobic, glucose-limited
    (0.45, False),   # B. thetaiotaomicron -- a normal gut anaerobe
    (1.20, False),   # fast, but real organisms do this
    (2.27, True),    # from the real run
    (6.53, True),    # the fastest model in the real run: 6-minute doubling
])
def test_sanity_check_flags_only_impossible_rates(rate, flagged):
    model = cobra.io.load_model("textbook")
    report = sanity_check_model(model, growth_rate=rate)
    assert report["implausible_growth"] is flagged
    assert report["passed"] is not flagged


def test_a_healthy_model_on_a_sane_medium_still_passes():
    """The ceiling must not become a blanket failure: no growth_rate, no verdict."""
    model = cobra.io.load_model("textbook")
    report = sanity_check_model(model)
    assert "implausible_growth" not in report
    assert report["passed"] is True


def test_implausible_models_are_flagged_but_never_silently_dropped(tmp_path):
    """The medium is at fault, not the model.

    Dropping the 25 offenders would delete a quarter of the community and hide
    the cause -- the simulation would just be wrong more quietly.  So
    `implausible_growth` is surfaced and `simulatable` is left alone.
    """
    paths = [
        _qc_json(tmp_path / "a.qc.json", "a", growth=0.42, implausible=False),
        _qc_json(tmp_path / "b.qc.json", "b", growth=6.53, implausible=True),
    ]
    df = summarize_reconstruction(paths).set_index("mag")

    assert bool(df.loc["b", "implausible_growth"]) is True
    assert bool(df.loc["b", "qc_passed"]) is False      # it is reported as broken
    assert bool(df.loc["b", "simulatable"]) is True     # ...but it is NOT deleted
    assert bool(df.loc["a", "qc_passed"]) is True


def test_the_real_run_would_have_been_caught():
    """The 89-model run that motivated all this: 25 models over the ceiling."""
    observed = [1.06, 2.27, 0.04, 0.08, 0.47, 1.45, 0.08, 3.98, 1.27, 0.37, 2.13,
                6.00, 5.86, 3.91, 6.53, 6.35, 4.47, 6.30, 5.37, 5.08, 3.27, 5.84]
    flagged = [g for g in observed if g > MAX_PLAUSIBLE_GROWTH]
    assert len(flagged) >= 14, "the ceiling must fire on the run that motivated it"
    assert max(observed) == 6.53
