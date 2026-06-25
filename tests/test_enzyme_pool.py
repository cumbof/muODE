"""Tests for the GECKO/sMOMENT protein-pool budget (`muode.enzyme`).

The pool coefficient is pure arithmetic and is tested directly.  The constraint
itself needs a real LP, so that test uses cobra's bundled `textbook` E. coli
model and is skipped where cobra is unavailable (core-only environments).
"""

import pytest

from muode.enzyme import DEFAULT_ENZYME_MW_KDA, _is_internal, pool_coefficient
from muode.kinetics import KineticParameters


def test_pool_coefficient_units():
    # coeff = MW / (kcat * 3600); 40 kDa at 10/s
    assert pool_coefficient(40.0, 10.0) == pytest.approx(40.0 / (10.0 * 3600.0))
    # a faster enzyme costs proportionally less mass per unit flux
    assert pool_coefficient(40.0, 20.0) == pytest.approx(0.5 * pool_coefficient(40.0, 10.0))
    # guard against division by zero / non-positive kcat
    assert pool_coefficient(40.0, 0.0) == 0.0


def test_pool_budget_binds_and_limits_growth():
    cobra = pytest.importorskip("cobra")

    def fresh():
        return cobra.io.load_model("textbook")

    kcat_per_s = 10.0
    coeff = pool_coefficient(DEFAULT_ENZYME_MW_KDA, kcat_per_s)

    # assign the same kcat to every internal reaction
    def kinetics_for(model):
        kin = KineticParameters()
        for rxn in model.reactions:
            if _is_internal(rxn):
                kin.set_kcat("ecoli", rxn.id, kcat_per_s)
        return kin

    # unconstrained optimum and the enzyme mass it would demand
    base = fresh()
    kin = kinetics_for(base)
    sol0 = base.optimize()
    mu0 = float(sol0.objective_value)
    demand0 = sum(
        coeff * abs(float(sol0.fluxes[r.id]))
        for r in base.reactions
        if _is_internal(r) and kin.get_kcat("ecoli", r.id) is not None
    )
    assert mu0 > 0.0 and demand0 > 0.0

    from muode.enzyme import apply_protein_pool_constraint

    # a generous budget (10x the demand) must not bite: growth is preserved
    loose = fresh()
    rep = apply_protein_pool_constraint(loose, kinetics_for(loose), "ecoli", pool_budget=10.0 * demand0)
    assert rep["n_pooled"] > 0
    mu_loose = float(loose.slim_optimize())
    assert mu_loose == pytest.approx(mu0, rel=1e-3)

    # a tight budget (half the demand) must force a lower maximum growth rate
    tight = fresh()
    apply_protein_pool_constraint(tight, kinetics_for(tight), "ecoli", pool_budget=0.5 * demand0)
    mu_tight = float(tight.slim_optimize())
    assert mu_tight < mu0 * 0.999
