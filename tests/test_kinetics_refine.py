"""Phase-3 kinetics refinement: the heuristic predictor, the kcat store and its
persistence, and the GECKO-lite enzyme-constraint layer.

The store / predictor-scalar tests are dependency-light; the model-level tests
(refinement, enzyme caps, namespace context) require cobra and use the bundled
textbook *E. coli* core model.
"""

import pytest

from muode.kinetics import KineticParameters, vmax_from_kcat
from muode.predict import HeuristicPredictor, get_predictor


# ---------------------------------------------------------------------------
# Heuristic predictor (no cobra)
# ---------------------------------------------------------------------------


def test_heuristic_km_table_and_jitter():
    p = HeuristicPredictor()
    assert p.predict_km_for("glc__D_e") == 0.5         # curated literature value
    assert p.predict_km_for("o2_e") == 0.024
    # unknown metabolite -> deterministic jitter in a sane range
    km = p.predict_km_for("weird_metabolite_e")
    assert 0.01 <= km < 1.0
    assert p.predict_km_for("weird_metabolite_e") == km   # reproducible


def test_heuristic_kcat_deterministic_and_varied():
    p = HeuristicPredictor()
    k1 = p.predict_kcat_for("orgA", "PFK")
    k2 = p.predict_kcat_for("orgA", "PGI")
    assert k1 == p.predict_kcat_for("orgA", "PFK")        # deterministic
    assert k1 != k2                                       # varies by reaction
    assert p.predict_kcat_for("orgB", "PFK") != k1        # varies by organism
    # within ~1 decade of the genome-wide median by construction
    assert 1.0 < k1 < 140.0


def test_get_predictor_registry():
    assert isinstance(get_predictor("heuristic"), HeuristicPredictor)
    with pytest.raises(ValueError):
        get_predictor("nope")


# ---------------------------------------------------------------------------
# kcat store, enzyme Vmax, persistence, merge (no cobra)
# ---------------------------------------------------------------------------


def test_enzyme_vmax_none_without_kcat():
    kin = KineticParameters(default_enzyme_concentration=0.01)
    assert kin.enzyme_vmax("o", "R1") is None
    kin.set_kcat("o", "R1", 10.0)
    assert kin.enzyme_vmax("o", "R1") == pytest.approx(vmax_from_kcat(10.0, 0.01))


def test_kinetics_json_roundtrip(tmp_path):
    kin = KineticParameters(metabolite_defaults={"glc_e": (8.0, 0.2)})
    kin.set("orgA", "ac_e", 9.0, 0.7)
    kin.set_kcat("orgA", "PFK", 22.0)

    kin2 = KineticParameters.from_dict(kin.to_dict())
    assert kin2.get("orgA", "ac_e") == (9.0, 0.7)
    assert kin2.get_kcat("orgA", "PFK") == 22.0
    assert kin2.get("anyone", "glc_e") == (8.0, 0.2)     # metabolite default survives

    path = tmp_path / "k.json"
    kin.to_json(path)
    kin3 = KineticParameters.from_json(path)
    assert kin3.get_kcat("orgA", "PFK") == 22.0


def test_kinetics_merge():
    a = KineticParameters()
    a.set("o", "m1", 1.0, 2.0)
    a.set_kcat("o", "r1", 5.0)
    b = KineticParameters()
    b.set("o", "m2", 3.0, 4.0)
    b.set_kcat("o", "r2", 6.0)
    a.merge(b)
    assert a.get("o", "m2") == (3.0, 4.0)
    assert a.get_kcat("o", "r2") == 6.0
    assert a.get_kcat("o", "r1") == 5.0


# ---------------------------------------------------------------------------
# Model-level refinement + enzyme constraints (cobra)
# ---------------------------------------------------------------------------


def test_refine_kinetics_fills_km_and_kcat():
    cobra = pytest.importorskip("cobra")
    from muode.predict import refine_kinetics

    model = cobra.io.load_model("textbook")
    kin = refine_kinetics(model, "Ecoli")
    # every internal reaction got a kcat; exchange metabolites got a Km
    assert kin.get_kcat("Ecoli", "PFK") is not None
    assert kin.get_kcat("Ecoli", "EX_glc__D_e") is None          # exchanges excluded
    vmax, km = kin.get("Ecoli", "glc__D_e")
    assert km == 0.5                                             # from the curated table


def test_enzyme_constraints_cap_flux_and_reduce_growth():
    cobra = pytest.importorskip("cobra")
    from muode.enzyme import apply_enzyme_constraints

    model = cobra.io.load_model("textbook")
    base_growth = model.slim_optimize()

    kin = KineticParameters(default_enzyme_concentration=0.01)
    kin.set_kcat("Ecoli", "PFK", 0.001)                          # a tiny turnover number
    report = apply_enzyme_constraints(model, kin, "Ecoli")

    assert report["n_constrained"] == 1 and "PFK" in report["reactions"]
    assert model.reactions.PFK.upper_bound == pytest.approx(vmax_from_kcat(0.001, 0.01))
    assert (model.slim_optimize() or 0.0) < base_growth          # PFK bottleneck bites


def test_cobra_organism_enzyme_caps_survive_reset():
    cobra = pytest.importorskip("cobra")
    from muode.organism import CobraOrganism

    org = CobraOrganism(cobra.io.load_model("textbook"), id="Ecoli")
    kin = KineticParameters(default_enzyme_concentration=0.01)
    kin.set_kcat("Ecoli", "PFK", 0.001)
    org.apply_enzyme_constraints(kin)

    capped = vmax_from_kcat(0.001, 0.01)
    assert org.model.reactions.PFK.upper_bound == pytest.approx(capped)
    org.set_uptake_bound("glc__D_e", 5.0)   # a per-step tweak
    org.reset_bounds()                       # dynamic loop resets each step
    assert org.model.reactions.PFK.upper_bound == pytest.approx(capped)   # cap persists


def test_build_enzyme_context_resolves_sequences(tmp_path):
    cobra = pytest.importorskip("cobra")
    from cobra import Metabolite, Model, Reaction

    from muode.predict import build_enzyme_context

    model = Model("m")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    a.annotation = {"smiles": "C"}
    r = Reaction("R1")
    r.add_metabolites({a: -1.0, b: 1.0})
    model.add_reactions([r])
    model.reactions.R1.gene_reaction_rule = "g1"

    faa = tmp_path / "p.faa"
    faa.write_text(">g1\nMSEQVENCE\n")
    ctx = build_enzyme_context(model, str(faa))
    assert ctx["R1"]["genes"] == ["g1"]
    assert ctx["R1"]["sequences"]["g1"] == "MSEQVENCE"
    assert ctx["R1"]["substrates"]["a_c"] == "C"


@pytest.mark.filterwarnings("ignore:Solver status")  # over-tight caps -> handled infeasible steps
def test_enzyme_constraints_change_dynamics():
    """A restrictive kcat propagates through the dFBA loop to slower growth."""
    cobra = pytest.importorskip("cobra")
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.diet import Diet
    from muode.organism import CobraOrganism

    diet = Diet(
        {"glc__D_e": 15.0, "o2_e": 1000.0, "nh4_e": 1000.0,
         "pi_e": 1000.0, "h2o_e": 1000.0, "h_e": 1000.0}, name="glc",
    )
    uptake = KineticParameters(
        default_vmax=10.0, default_km=0.1,
        metabolite_defaults={"glc__D_e": (10.0, 0.5), "o2_e": (15.0, 0.1)},
    )

    free = Community([CobraOrganism(cobra.io.load_model("textbook"), id="Ecoli")], total_biomass=0.01)
    res_free = DynamicFBA(t_end=8.0, dt=0.1).run(free, diet, uptake)

    org = CobraOrganism(cobra.io.load_model("textbook"), id="Ecoli")
    caps = KineticParameters(default_enzyme_concentration=0.01)
    caps.set_kcat("Ecoli", "PFK", 0.1)             # throttle the committed glycolytic step
    org.apply_enzyme_constraints(caps)
    res_con = DynamicFBA(t_end=8.0, dt=0.1).run(Community([org], total_biomass=0.01), diet, uptake)

    assert res_con.biomass["Ecoli"].iloc[-1] < res_free.biomass["Ecoli"].iloc[-1]


def test_ml_predictors_are_optin():
    from muode.external import MuodeToolError
    from muode.predict import DLKcatPredictor

    try:
        p = DLKcatPredictor()
    except MuodeToolError:
        return                                   # torch absent -> opt-in gate fired
    with pytest.raises(MuodeToolError):          # torch present -> scoring not bundled
        p.predict_kcat("MSEQ", "C(C(=O)O)")
