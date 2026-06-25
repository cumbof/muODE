"""Enzyme-constrained bounds -- the *intracellular* role of kcat (Phase 3b).

A turnover number does not enter a substrate uptake bound; it caps how fast the
enzyme that catalyses a reaction can run: ``Vmax = kcat * [E]``.  Enzyme-
constrained models (GECKO, sMOMENT, ecModels) impose exactly this.  muODE
implements the tractable, widely-used *per-reaction cap* form here:

    |v_r| <= kcat_r * [E]_r        (mmol / gDW / h)

applied to every internal reaction for which a kcat is known, leaving exchange
and biomass pseudo-reactions untouched.

On top of the per-reaction caps muODE also offers the shared protein-*pool*
budget of GECKO/sMOMENT (:func:`apply_protein_pool_constraint`)::

    sum_r ( MW_r / (kcat_r * 3600) ) * |v_r|  <=  P            (g enzyme / gDW)

i.e. every reaction draws on one finite enzyme mass budget ``P`` (the fraction of
the proteome available to the constrained reactions), so the cell must *allocate*
its limited enzyme between competing pathways rather than running them all at the
per-reaction maximum.  Molecular weights (kDa = g/mmol) default to a typical
enzyme mass and can be overridden per reaction.  This is still GECKO-*lite* — a
single global budget with default masses, not a measured per-enzyme proteome
allocation — but it captures the central proteome-constraint phenomenon
(overflow metabolism, substrate hierarchies) that per-reaction caps alone cannot.
Both layers are opt-in; with the heuristic predictor the kcat/MW values are
placeholders, so calibrate before any quantitative claim.
"""

from __future__ import annotations

from typing import Mapping, Optional

from muode.kinetics import KineticParameters

#: Default enzyme molecular weight (kDa = g/mmol) when none is supplied; ~the
#: median soluble-enzyme subunit mass.  Only the *ratio* MW/kcat matters for the
#: relative enzyme cost of one reaction versus another.
DEFAULT_ENZYME_MW_KDA = 40.0


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


def pool_coefficient(mw_kda: float, kcat_per_s: float) -> float:
    """Enzyme-mass cost (g enzyme / gDW) per unit flux (mmol/gDW/h) of a reaction.

    ``coeff = MW / (kcat * 3600)`` with ``MW`` in kDa (= g/mmol) and ``kcat`` in
    1/s.  Multiplying ``coeff`` by a reaction's flux gives the enzyme mass that
    flux requires; summing the products over all enzyme-constrained reactions is
    the quantity the protein-pool budget caps.
    """
    denom = float(kcat_per_s) * 3600.0
    if denom <= 0.0:
        return 0.0
    return float(mw_kda) / denom


def apply_protein_pool_constraint(
    model,
    kinetics: KineticParameters,
    organism_id: str,
    pool_budget: float,
    default_mw: float = DEFAULT_ENZYME_MW_KDA,
    mw: Optional[Mapping[str, float]] = None,
) -> dict:
    """Impose one shared enzyme-mass budget across all kcat-constrained reactions.

    Adds the GECKO/sMOMENT proteome-pool constraint to a cobra ``model``::

        sum_r ( MW_r / (kcat_r * 3600) ) * |v_r|  <=  pool_budget

    Reversible reactions are handled exactly via an auxiliary non-negative usage
    variable ``a_r >= |v_r|`` (two linear constraints), so no reaction splitting
    or model rewrite is needed.  The constraint is added to the solver and
    therefore *survives* the per-step ``reset_bounds`` of the dynamic loop (which
    only restores reaction bounds), so it stays in force for the whole run.

    Parameters
    ----------
    pool_budget:
        Total enzyme mass available, in g enzyme / gDW (a fraction of the
        proteome, e.g. ~0.2-0.5).  Smaller values force tighter allocation.
    default_mw, mw:
        Per-reaction enzyme molecular weight (kDa = g/mmol); ``mw`` overrides
        ``default_mw`` for specific reaction ids.

    Returns a report dict (number of reactions drawing on the pool, the budget).
    """
    mw = dict(mw or {})
    prob = model.problem
    from optlang.symbolics import Zero

    pool_expr = Zero
    aux_cons_vars = []
    n = 0
    for rxn in model.reactions:
        if not _is_internal(rxn):
            continue
        kcat = kinetics.get_kcat(organism_id, rxn.id)
        if kcat is None:
            continue
        coeff = pool_coefficient(mw.get(rxn.id, default_mw), kcat)
        if coeff <= 0.0:
            continue
        usage = prob.Variable(f"mu_prot_{rxn.id}", lb=0.0)
        # a_r >= v_r  and  a_r >= -v_r   =>   a_r >= |v_r|
        c_pos = prob.Constraint(usage - rxn.flux_expression, lb=0.0, name=f"mu_ec_pos_{rxn.id}")
        c_neg = prob.Constraint(usage + rxn.flux_expression, lb=0.0, name=f"mu_ec_neg_{rxn.id}")
        aux_cons_vars += [usage, c_pos, c_neg]
        pool_expr = pool_expr + coeff * usage
        n += 1

    if n == 0:
        return {"organism": organism_id, "n_pooled": 0, "pool_budget": float(pool_budget)}

    budget = prob.Constraint(pool_expr, ub=float(pool_budget), name="mu_prot_pool")
    model.add_cons_vars(aux_cons_vars + [budget])
    model.solver.update()
    return {"organism": organism_id, "n_pooled": n, "pool_budget": float(pool_budget)}
