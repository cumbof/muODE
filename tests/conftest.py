"""Shared fixtures -- and the reason the western_gut tests need a scaffold at all.

**iJO1366 cannot grow on western_gut, and that is a result, not an obstacle.**

Until the medium was audited it grew at 0.0516/h, and every bit of that came from
nine nucleoside rows sitting at the source table's fill value of 0.1 mmol/gDW/h.  A
nucleoside is a ~10-carbon molecule, so those rows were feeding E. coli roughly 9 mmol
C/gDW/h -- an RNA diet nobody chose.  Bound them by real dietary nucleic acid (0.1-1
g/day) and they fall ~200x, and E. coli stops growing entirely:

    max ATP it can make on the honest medium:  1.97 mmol/gDW/h
    its maintenance demand (ATPM):             3.15 mmol/gDW/h

so the LP is *infeasible* rather than merely zero.  That is correct biology.  E. coli
ferments sugar, and a colonic medium's sugar is locked inside fibre polymers -- starch,
xylan, cellobiose, pullulan -- that it has no enzymes to open.  In a real gut E. coli
is a **cross-feeder**, living on mucus sugars and on what the fibre degraders release;
a monoculture on colonic contents is not a situation it is adapted to.  Its ATPM of
3.15 is fitted on aerobic rich medium, and is high for an anaerobe.

That the medium is now too lean for iJO1366 is a sign it is in the right physiological
range, not the wrong one: it supplies ~2-3 mmol ATP/gDW/h per gram of community
biomass, which is exactly why real gut communities grow at 0.03-0.05/h rather than at
E. coli's textbook 0.5/h.

But the medium's *other* rows still need exercising -- are the mineral bounds slack? is
a vitamin bound limiting? -- and that needs a model that is alive.  Hence a scaffold.
"""

from __future__ import annotations

import cobra
import pytest

from muode.diet import Diet, load_diet

# ---------------------------------------------------------------------------
# Shared heavy objects.
#
# These fixtures exist for SPEED as much as for tidiness, and the difference is
# large enough to be worth stating.  `nutrient_sensitivity` runs two FBAs per
# nutrient -- about 240 on a genome-scale model against a 113-row diet -- and
# `diagnose_no_growth` runs one per candidate rescue, about 226.  Recomputing either
# once per test costs 10-30s a time, and five modules were each loading their own
# copy of iJO1366 on top of that.
#
# Session scope makes each of them happen exactly once.  Nothing here is mutated:
# cobra's `with model:` context manager rolls back every bound change, and the diet
# objects are read-only in practice.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def bigg():
    """Load a BiGG model by name, at most once per test session.

    ``cobra.io.load_model`` does NOT cache.  Every call re-parses (iJO1366: ~4.5s) or,
    for anything not bundled with cobra, **re-downloads from the BiGG API** (iCN900:
    ~3.3s and a network round-trip).  Five test modules were each loading their own
    iJO1366, which is ~20s of the suite spent parsing the same file five times -- and on
    a slow or throttled CI network the download path turns that into minutes.

    Sharing one instance is safe because nothing here mutates a model persistently:
    every consumer changes bounds inside cobra's ``with model:`` block, which rolls them
    back on exit (pinned by test_growth_on_diet_does_not_mutate_the_model).
    """
    cache: dict = {}

    def _load(name: str):
        if name not in cache:
            cache[name] = cobra.io.load_model(name)
        return cache[name]

    return _load


@pytest.fixture(scope="session")
def ijo(bigg):
    """iJO1366, loaded once for the whole suite."""
    return bigg("iJO1366")


@pytest.fixture(scope="session")
def western_gut_probe():
    """western_gut plus free glucose: a stand-in for the community E. coli cannot be.

    **This is a test scaffold and not a claim about the colon.**  Its only job is to
    put a model above its maintenance threshold so the medium's other rows can be
    interrogated; nothing derived from it should be reported as a property of the diet.

    The glucose is not arbitrary.  The medium supplies starch1200, pullulan1200,
    amylose300 and cellobiose, and a real community's amylases and glucanases hydrolyse
    them:

        1e-4 * 1200  +  1e-4 * 1200  +  1e-4 * 300  +  0.0745 * 2  =  0.42 mmol glc/gDW/h

    That is sugar which genuinely exists in this medium and which E. coli genuinely
    cannot reach on its own.  Even handed all of it, iJO1366 reaches 3.02 mmol ATP/gDW/h
    against the 3.15 it needs -- still short.  So the scaffold rounds up to 1.0, which
    clears maintenance with room to grow, and the gap between 0.42 and 1.0 is the size
    of the fudge.  It is written down here rather than hidden in a constant.
    """
    diet = load_diet("western_gut")
    conc = dict(diet.concentrations)
    upt = dict(diet.max_uptake)
    conc["glc__D_e"] = 1.0 * 20.0 * 24.0
    upt["glc__D_e"] = 1.0
    return Diet(
        concentrations=conc,
        influx=dict(diet.influx),
        max_uptake=upt,
        source=dict(diet.source),
        name="western_gut + cross-fed glucose (TEST SCAFFOLD)",
    )


@pytest.fixture(scope="session")
def probe_sensitivity(ijo, western_gut_probe):
    """``nutrient_sensitivity`` on the probe diet -- ~240 FBAs, so computed once."""
    from muode.qc import nutrient_sensitivity

    return nutrient_sensitivity(ijo, western_gut_probe)


@pytest.fixture(scope="session")
def dm38_nofe3():
    """DM38 with ferric iron removed -- reconstructs the historical Fe(III) trap.

    DM38 as shipped now supplies ``fe3_e`` (commit 2233c82: the Clark GEMs' biomass
    demands cytoplasmic Fe(III) and nothing converts fe2->fe3, so the medium offers
    the same iron pool in both oxidation states).  That deliberately DISARMS the
    trap -- iJO1366 grows on DM38 now.  But the trap is the canonical case this
    diagnosis machinery exists for, so we reconstruct it by dropping ``fe3_e``,
    which is exactly the medium these tests were first written against.
    """
    dm38 = load_diet("dm38")
    return Diet(
        concentrations={m: c for m, c in dm38.concentrations.items() if m != "fe3_e"},
        influx=dict(dm38.influx),
        max_uptake=dict(dm38.max_uptake),
        source=dict(dm38.source),
        name="DM38 (fe3 removed -- reconstructs the Fe(III) trap)",
    )


@pytest.fixture(scope="session")
def dm38_nofe3_diagnosis(ijo, dm38_nofe3):
    """``diagnose_no_growth`` on iJO1366 / fe3-less DM38 -- ~226 FBAs, computed once.

    Five tests in test_no_growth_diagnosis.py interrogate different facets of this
    one report (the verdict, the iron, the nitric-oxide artefact, the menu length,
    the forbidden oxygen).  They were each recomputing it.
    """
    from muode.qc import diagnose_no_growth

    return diagnose_no_growth(ijo, dm38_nofe3)
