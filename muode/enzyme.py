"""Enzyme-constrained bounds -- the *intracellular* role of kcat (Phase 3b).

A turnover number does not enter a substrate uptake bound; it caps how fast the
enzyme that catalyses a reaction can run: ``Vmax = kcat * [E]``.  Enzyme-
constrained models (GECKO, sMOMENT, ecModels) impose exactly this.  muODE
implements the tractable, widely-used *per-reaction cap* form here:

    |v_r| <= kcat_r * [E]_r        (mmol / gDW / h)

applied to every internal reaction for which a kcat is known, leaving exchange
and biomass pseudo-reactions untouched.  This is a deliberate simplification of
full GECKO -- it does **not** yet impose a shared protein-*pool* budget
(``sum_r MW_r * [E]_r <= P``), which needs per-enzyme molecular weights and a
calibrated proteome fraction.  The pool constraint is the natural next step once
real kcat (DLKcat) and proteomics are wired in; until then, kcat values from the
:mod:`muode.predict` heuristic are placeholders and this layer is opt-in.
"""

from __future__ import annotations

from typing import Optional

from muode.kinetics import KineticParameters


def _is_internal(rxn) -> bool:
    """A reaction we may enzyme-constrain (not exchange/boundary/biomass/objective)."""
    return (
        not rxn.boundary
        and "biomass" not in rxn.id.lower()
        and getattr(rxn, "objective_coefficient", 0.0) == 0.0
    )


def apply_enzyme_constraints(
    model,
    kinetics: KineticParameters,
    organism_id: str,
    enzyme_concentration: Optional[float] = None,
) -> dict:
    """Cap intracellular reaction velocities of a cobra ``model`` from kcat.

    For every internal reaction with a known kcat, the forward capacity is capped
    at ``Vmax = kcat * [E]`` and (if the reaction is reversible) the reverse
    capacity at ``-Vmax``.  Reactions without a kcat are left as-is.  Mutates the
    model in place and returns a small report.
    """
    constrained = []
    for rxn in model.reactions:
        if not _is_internal(rxn):
            continue
        vmax = kinetics.enzyme_vmax(organism_id, rxn.id, enzyme_concentration)
        if vmax is None:
            continue
        vmax = abs(float(vmax))
        new_ub = min(rxn.upper_bound, vmax)
        new_lb = max(rxn.lower_bound, -vmax)
        if new_ub != rxn.upper_bound or new_lb != rxn.lower_bound:
            rxn.upper_bound, rxn.lower_bound = new_ub, new_lb
            constrained.append(rxn.id)
    e = kinetics.default_enzyme_concentration if enzyme_concentration is None else enzyme_concentration
    return {
        "organism": organism_id,
        "n_constrained": len(constrained),
        "enzyme_concentration": e,
        "reactions": constrained[:20],
    }
