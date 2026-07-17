"""The multi-kingdom demo teaches two cross-kingdom mechanisms -- pin both directions.

`examples/multikingdom/mechanistic_demo.py` runs on any machine (toy LinprogOrganism
models, no eukaryote reconstruction needed) and is what a reader executes to see why a
sample that is *not* bacteria-only behaves differently.  It makes two falsifiable,
ablatable claims:

  1. a facultative fungus scavenges mucosal O2 and PROTECTS the obligate-anaerobe
     keystone that O2 would otherwise poison;
  2. a lytic phage crashes a facultative pathobiont's bloom.

Nothing guarded them.  This pins each mechanism by its ablation -- and the oxygen one is
pinned through its CAUSE (the O2 pool), not just its effect, so a change that helped
the anaerobe for the wrong reason still fails.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

# See test_strain_competition_demo: both examples name their module `mechanistic_demo`,
# so load by path under a unique name to avoid a sys.modules collision between the two.
_DEMO = Path(__file__).resolve().parents[1] / "examples" / "multikingdom" / "mechanistic_demo.py"
_spec = importlib.util.spec_from_file_location("multikingdom_demo", _DEMO)
mk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mk)

# Three 36 h dynamic-FBA integrations; slow for the reason the test needs (it runs the
# demo's real dynamics), so it joins the -m "not slow" exclusion.
pytestmark = pytest.mark.slow


def _snapshot(with_fungus: bool, with_phage: bool) -> dict:
    res = mk.run(with_fungus=with_fungus, with_phage=with_phage)
    return {"bio": res.final_biomass(), "o2": float(res.environment["oxygen"].iloc[-1])}


@pytest.fixture(scope="module")
def outcomes():
    return {
        "full": _snapshot(True, True),      # fungus + phage
        "no_fungus": _snapshot(False, True),
        "no_phage": _snapshot(True, False),
    }


def test_the_fungus_scavenges_oxygen(outcomes):
    """The CAUSE: pulling the fungus leaves the mucosal O2 pool high instead of drawn down."""
    assert outcomes["no_fungus"]["o2"] > outcomes["full"]["o2"]
    # scavenged to near zero with the fungus; left high without it -- not a marginal shift
    assert outcomes["full"]["o2"] < 0.5
    assert outcomes["no_fungus"]["o2"] > 1.0


def test_that_oxygen_scavenging_protects_the_obligate_anaerobe(outcomes):
    """The EFFECT: the O2 the fungus removes is what the anaerobe keystone would be poisoned by."""
    bact_full = outcomes["full"]["bio"][mk.BACTEROIDES]
    bact_bare = outcomes["no_fungus"]["bio"][mk.BACTEROIDES]
    assert bact_full > bact_bare
    # cross-kingdom protection is large, not a nudge: the anaerobe is suppressed without it
    assert bact_full > 10 * bact_bare


def test_a_lytic_phage_checks_the_pathobiont_bloom(outcomes):
    """Ablating the phage lets the facultative pathobiont bloom unchecked."""
    kleb_checked = outcomes["full"]["bio"][mk.KLEBSIELLA]
    kleb_bloom = outcomes["no_phage"]["bio"][mk.KLEBSIELLA]
    assert kleb_bloom > kleb_checked
    # a real bloom vs a crashed host, not a small difference
    assert kleb_bloom > 5.0
    assert kleb_checked < 0.5
