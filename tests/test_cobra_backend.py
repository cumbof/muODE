"""Production-backend tests: dynamic FBA on a real genome-scale model.

Skipped automatically when cobra is not installed (e.g. core-only environments).
The assertions reproduce the canonical batch-growth dFBA result for the textbook
E. coli core model: growth on glucose with acetate overflow followed by
re-uptake once glucose is exhausted (Mahadevan et al., Biophys. J. 2002).
"""

import warnings

import pytest

cobra = pytest.importorskip("cobra")

from muode.community import Community
from muode.dfba import DynamicFBA
from muode.diet import Diet
from muode.kinetics import KineticParameters
from muode.organism import CobraOrganism


@pytest.fixture(scope="module")
def ecoli_community():
    warnings.filterwarnings("ignore")
    model = cobra.io.load_model("textbook")
    org = CobraOrganism(model, id="Ecoli")
    return Community([org], abundances={"Ecoli": 1.0}, total_biomass=0.01)


def test_batch_growth_and_acetate_switch(ecoli_community):
    diet = Diet(
        concentrations={
            "glc__D_e": 15.0, "o2_e": 1000.0, "nh4_e": 1000.0,
            "pi_e": 1000.0, "h2o_e": 1000.0, "h_e": 1000.0,
        },
        name="glc_aerobic",
    )
    kin = KineticParameters(
        default_vmax=10.0, default_km=0.1,
        metabolite_defaults={"glc__D_e": (10.0, 0.5), "o2_e": (15.0, 0.1)},
    )
    res = DynamicFBA(t_end=10.0, dt=0.1).run(ecoli_community, diet, kin)

    biomass = res.biomass["Ecoli"]
    glucose = res.metabolites["glc__D_e"]
    acetate = res.metabolites["ac_e"]

    assert biomass.iloc[-1] > 2 * biomass.iloc[0]      # substantial growth
    assert glucose.iloc[-1] < 1.0                      # glucose consumed
    assert acetate.max() > 1.0                         # overflow metabolism
    assert acetate.iloc[-1] < acetate.max()            # re-uptake after switch


def test_pathway_perturbation_resolves_subsystem():
    # The bundled textbook model ships without subsystem annotations, so set a
    # couple explicitly and confirm case-insensitive substring resolution.
    model = cobra.io.load_model("textbook")
    model.reactions.PFK.subsystem = "Glycolysis/Gluconeogenesis"
    model.reactions.PGI.subsystem = "Glycolysis/Gluconeogenesis"
    org = CobraOrganism(model, id="Ecoli")
    rxns = org.reactions_in_subsystem("glycolysis")
    assert set(rxns) == {"PFK", "PGI"}


def test_pathway_attenuation_reduces_capacity():
    from muode.perturb import Perturbation

    model = cobra.io.load_model("textbook")
    model.reactions.PFK.subsystem = "Glycolysis"
    org = CobraOrganism(model, id="Ecoli")
    comm = Community([org], total_biomass=0.01)
    Perturbation.antibiotic("Glycolysis", efficacy=0.5).apply(comm)
    org.reset_bounds()
    # PFK upper bound (default 1000) should be halved by the 50%-efficacy hit
    assert org.model.reactions.PFK.upper_bound == pytest.approx(500.0)
