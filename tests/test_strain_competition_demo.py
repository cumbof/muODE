"""The strain-competition demo teaches three ecological laws -- pin the ORDERING.

`examples/strain_competition/mechanistic_demo.py` is one of the example parts that
actually runs on any machine (toy LinprogOrganism models, no CarveMe), so it is what a
reader executes to see the mechanism.  It makes three falsifiable claims, one per
competition mode, and nothing guarded them -- a change to the engine, the Bacteriocin
layer or the kinetics could silently invert a winner and the demo would still "run".

This test pins the *ordering* (who wins in each mode), not the exact biomass: the
ordering IS the ecology (Gause / Tilman R*, interference beating exploitation, Freter's
nutrient niche), and it must survive numeric drift while a flipped winner must fail.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_DEMO_DIR = Path(__file__).resolve().parents[1] / "examples" / "strain_competition"
sys.path.insert(0, str(_DEMO_DIR))
import mechanistic_demo as scd  # noqa: E402

# Runs three 48 h dynamic-FBA integrations (~2 min); slow for the reason the test needs
# -- it exercises the demo's real dynamics -- so it joins the -m "not slow" exclusion.
pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def outcomes():
    """Final biomass for each of the three scenarios (deterministic; run once)."""
    return {
        "resource": scd.run(with_colicin=False, with_arabinose=False).final_biomass(),
        "interference": scd.run(with_colicin=True, with_arabinose=False).final_biomass(),
        "niche": scd.run(with_colicin=True, with_arabinose=True).final_biomass(),
    }


def test_resource_competition_excludes_all_but_the_best_glucose_strain(outcomes):
    """Mode 1: one limiting sugar -> the lowest-R* strain competitively excludes the rest."""
    fb = outcomes["resource"]
    assert fb[scd.EFFICIENT] > fb[scd.COLICIN]
    assert fb[scd.EFFICIENT] > fb[scd.NICHE]
    # exclusion, not a photo-finish: the specialist dominates by a wide margin
    assert fb[scd.EFFICIENT] > 3 * fb[scd.COLICIN]


def test_a_colicin_overturns_the_metabolic_winner(outcomes):
    """Mode 2: interference beats exploitation -- the weaker PRODUCER wins, specialist collapses."""
    resource, interference = outcomes["resource"], outcomes["interference"]
    # the flip: the producer now outbiomasses the strain that dominated on nutrients alone
    assert interference[scd.COLICIN] > interference[scd.EFFICIENT]
    # ...and it is a collapse of the specialist, not merely the producer catching up
    assert interference[scd.EFFICIENT] < resource[scd.EFFICIENT]


def test_a_private_niche_buys_coexistence(outcomes):
    """Mode 3: a private substrate lets the niche strain escape exclusion and dominate."""
    resource, niche = outcomes["resource"], outcomes["niche"]
    assert niche[scd.NICHE] > niche[scd.EFFICIENT]
    assert niche[scd.NICHE] > niche[scd.COLICIN]
    # the private sugar is what does it: the niche strain ends far higher than when
    # glucose was the only carbon source
    assert niche[scd.NICHE] > resource[scd.NICHE]
