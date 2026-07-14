"""Fitting Vmax to Clark's measurements -- and refusing to when the data cannot.

The trap this guards: hand an optimiser an endpoint-only objective and it will
always return a confident number.  In a 48 h batch culture that runs to substrate
exhaustion the endpoint is set by *yield*, not by *rate*, so that number is fitted
noise on a plateau.  A calibration that cannot say "the data does not determine
this" is worse than no calibration, because it launders an invented parameter into
a measured-looking one.
"""

import warnings

import cobra
import pytest

from muode.calibrate import Identifiability, assess, fit_vmax, profile
from muode.community import Community
from muode.dfba import DynamicFBA
from muode.diet import load_diet
from muode.kinetics import KineticParameters
from muode.organism import CobraOrganism

INOCULUM = 0.01


@pytest.fixture(scope="module")
def predict_dm38():
    """Endpoint of E. coli core in DM38 after 48 h, as a function of Vmax."""
    model = cobra.io.load_model("textbook")
    diet = load_diet("dm38")

    def predict(vmax: float):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            org = CobraOrganism(model.copy(), id="e")
            comm = Community([org], abundances={"e": 1.0}, total_biomass=INOCULUM)
            res = DynamicFBA(t_end=48.0, dt=0.2).run(comm, diet, KineticParameters(
                default_vmax=float(vmax)))
        return {
            "biomass": float(res.biomass["e"].iloc[-1]),
            "ac_e": float(res.metabolites["ac_e"].iloc[-1]),
        }

    return predict


# --- the real thing: a batch endpoint does not identify a rate -----------------

@pytest.mark.slow
def test_endpoint_data_bounds_vmax_but_does_not_identify_it(predict_dm38):
    """The finding, on the real DM38 medium and a real BiGG model.

    Below ~4 mmol/gDW/h the cell cannot pay its ATP maintenance and dies.  Above
    ~6 the glucose is gone well before 48 h and the endpoint stops moving.  So the
    measurement says "Vmax is at least this big" and nothing more.
    """
    ident, prof = fit_vmax(
        predict_dm38,
        observed={"ac_e": 21.0},                 # any endpoint on the plateau
        vmax_grid=[2, 4, 6, 10, 20, 50, 100],
        inoculum=INOCULUM,
    )
    assert ident.verdict == "bounded_below"
    assert ident.best is None, "a point estimate here would be fitting noise"
    assert ident.lower_bound is not None and ident.lower_bound >= 4
    assert "NOT identified" in ident.summary()

    # and the plateau is real: a 10x change in Vmax barely moves the observable
    on_plateau = [p for v, p in zip(prof.vmax, prof.predicted["ac_e"]) if v >= 10]
    assert max(on_plateau) - min(on_plateau) < 0.10 * max(on_plateau)


@pytest.mark.slow
def test_a_starved_model_is_reported_as_dead_not_as_a_bad_fit(predict_dm38):
    """No growth at any Vmax indicts the model or the medium -- never the kinetics."""
    ident, _ = fit_vmax(
        predict_dm38, observed={"ac_e": 21.0},
        vmax_grid=[0.01, 0.1, 0.5],              # far too low to pay maintenance
        inoculum=INOCULUM,
    )
    assert ident.verdict == "no_growth"
    assert ident.best is None
    assert "fix the model or the medium" in ident.summary()


# --- the classifier itself, on synthetic profiles ------------------------------

def _prof(grid, biomass, target, observed):
    return profile(
        lambda v: {"biomass": biomass[grid.index(v)], "y": target[grid.index(v)]},
        grid, observed=observed,
    )


def test_an_interior_optimum_is_reported_as_identified():
    """When the data DOES pin Vmax down, say so and give the number."""
    grid = [1.0, 2.0, 4.0, 8.0, 16.0]
    biomass = [0.5] * 5                       # grows everywhere
    target = [1.0, 3.0, 5.0, 3.0, 1.0]        # peaks at vmax=4
    ident = assess(_prof(grid, biomass, target, {"y": 5.0}), inoculum=0.01)
    assert ident.verdict == "identified"
    assert ident.best == 4.0
    assert ident.fitted


def test_a_monotone_objective_is_a_bound_not_an_optimum():
    """The objective falling to the edge of the sweep is not a minimum."""
    grid = [1.0, 2.0, 4.0, 8.0, 16.0]
    biomass = [0.5] * 5
    target = [1.0, 2.0, 3.0, 4.0, 5.0]        # still climbing at the last point
    ident = assess(_prof(grid, biomass, target, {"y": 5.0}), inoculum=0.01)
    assert ident.verdict == "bounded_below"
    assert ident.best is None


def test_a_flat_prediction_is_unidentifiable_not_a_perfect_fit():
    """A parameter that moves nothing must never be reported as fitted."""
    grid = [1.0, 2.0, 4.0, 8.0]
    biomass = [0.5] * 4
    target = [3.00, 3.01, 2.99, 3.00]         # inert
    ident = assess(_prof(grid, biomass, target, {"y": 3.0}), inoculum=0.01)
    assert ident.verdict == "unidentifiable"
    assert ident.best is None
    assert ident.sensitivity < 0.10


def test_dead_cultures_do_not_count_as_evidence_that_vmax_matters():
    """A zero prediction from a dead cell is not sensitivity to the parameter.

    Without this, the huge jump from "dead" to "alive" reads as a strong dependence
    on Vmax and the whole sweep is misclassified as `identified`.
    """
    grid = [1.0, 2.0, 4.0, 8.0]
    biomass = [0.01, 0.01, 0.5, 0.5]          # dead, dead, alive, alive
    target = [0.0, 0.0, 3.0, 3.0]             # the 0s are death, not a rate effect
    ident = assess(_prof(grid, biomass, target, {"y": 3.0}), inoculum=0.01)
    assert ident.verdict == "bounded_below"
    assert ident.lower_bound == 4.0           # the Vmax at which it first lived
    assert ident.best is None


def test_a_culture_still_running_is_not_evidence_that_vmax_matters():
    """The bug this caught, and the reason `assess` looks for a plateau SUFFIX.

    At a low-but-viable Vmax the culture is alive and still climbing when the 48 h
    ends: its endpoint is low because it did not *finish*, not because the rate
    changed the outcome.  Measuring sensitivity across the whole live range reads
    that climb as "the parameter matters", finds an interior optimum, and hands back
    a fitted Vmax -- which is precisely the laundering this module exists to stop.
    """
    grid = [4.0, 6.0, 10.0, 20.0, 50.0, 100.0]
    biomass = [0.04, 0.39, 0.53, 0.61, 0.60, 0.48]     # all alive
    target = [3.9, 22.4, 21.5, 20.7, 20.9, 20.4]       # climbing, THEN flat
    ident = assess(_prof(grid, biomass, target, {"y": 21.0}), inoculum=0.01)

    assert ident.verdict == "bounded_below"
    assert ident.best is None
    assert ident.plateau_from == 6.0        # flat from here up
    assert ident.sensitivity < 0.10


def test_summary_never_launders_a_bound_into_a_value():
    for verdict in ("bounded_below", "unidentifiable", "no_growth"):
        ident = Identifiability(verdict, None, 4.0, 4.0, 0.02)
        assert ident.best is None
        assert not ident.fitted
        assert "identified: " not in ident.summary()


def test_an_optimum_on_the_edge_of_the_grid_is_not_an_optimum():
    """If the best fit is the lowest Vmax we tried, we never looked low enough.

    Reporting it as `identified` would turn an artefact of the grid into a measured
    parameter -- the same laundering, one step removed.
    """
    grid = [2.0, 5.0, 10.0, 20.0, 50.0]
    biomass = [0.4] * 5                       # alive throughout
    target = [3.0, 8.0, 9.0, 9.05, 9.02]      # rises, then plateaus from 10
    # observation matches the LOWEST grid point -> optimum sits on the edge
    ident = assess(_prof(grid, biomass, target, {"y": 3.0}), inoculum=0.01)
    assert ident.verdict != "identified"
    assert ident.best is None
