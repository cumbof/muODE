"""The fungal-GEM exchange harmonizer: rename exchange metabolites to their BiGG ids.

Without this the multi-kingdom fungus (Yeast8, `s_0565` = glucose) and the CarveMe
bacteria (`glc__D_e`) never touch the same molecule, so the fungus cannot scavenge the
shared O2 pool -- its whole role.  The logic is small (read a bigg.metabolite annotation,
rename the metabolite) but load-bearing, so pin its three cases: rename, collision, and
no-annotation.  Fast -- toy cobra models, no Yeast8 download needed.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

cobra = pytest.importorskip("cobra")

_MOD = Path(__file__).resolve().parents[1] / "examples" / "multikingdom" / "harmonize_fungal_gem.py"
_spec = importlib.util.spec_from_file_location("mk_harmonize", _MOD)
hz = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hz)


def _exchange(model, met_id, bigg, compartment="e"):
    met = cobra.Metabolite(met_id, compartment=compartment)
    if bigg is not None:
        met.annotation = {"bigg.metabolite": bigg}
    rxn = cobra.Reaction(f"EX_{met_id}")
    rxn.add_metabolites({met: -1.0})
    rxn.lower_bound = -10.0
    model.add_reactions([rxn])
    return met


def test_an_annotated_exchange_metabolite_is_renamed_to_its_bigg_id():
    m = cobra.Model("t")
    _exchange(m, "s_0565", "glc__D")     # Yeast8's glucose
    _exchange(m, "s_1277", "o2")         # Yeast8's oxygen
    renamed, skipped, collisions = hz.harmonize(m)

    ids = {met.id for met in m.metabolites}
    assert "glc__D_e" in ids and "o2_e" in ids
    assert renamed == 2 and skipped == 0 and collisions == 0


def test_an_unannotated_exchange_is_left_native_and_counted():
    m = cobra.Model("t")
    _exchange(m, "s_9999", None)         # a yeast aroma ester: no bigg annotation
    renamed, skipped, collisions = hz.harmonize(m)

    assert {met.id for met in m.metabolites} == {"s_9999"}   # unchanged
    assert skipped == 1 and renamed == 0


def test_a_collision_keeps_the_first_and_leaves_the_second_native():
    """Two exchanges annotated to the same bigg id must not both become glc__D_e (that is
    not a valid model); the second is reported and left alone."""
    m = cobra.Model("t")
    _exchange(m, "s_0565", "glc__D")
    _exchange(m, "s_0566", "glc__D")     # a second thing claiming glucose
    renamed, skipped, collisions = hz.harmonize(m)

    ids = {met.id for met in m.metabolites}
    assert "glc__D_e" in ids            # the first won
    assert "s_0566" in ids              # the second kept its native id
    assert collisions == 1
