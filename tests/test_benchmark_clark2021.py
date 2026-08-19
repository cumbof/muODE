"""The Clark et al. 2021 benchmark: real measurements, real medium, real strains.

These tests run against a committed excerpt of the authors' measurement table
(`tests/data/clark2021_excerpt.csv`, CC-BY 4.0) -- every monoculture, a few
communities, plus contaminated wells and a COMM* pool so each filtering branch is
exercised.  They guard the things that would quietly corrupt the benchmark:

* dropping a strain because its code was not in our table (HB vs BH),
* scoring lactate as an endpoint when 28.3 mM of it is *in the medium*,
* scoring contaminated wells,
* penalising a model for not growing a strain that did not grow in vitro either.
"""

import sys
from pathlib import Path

import pytest

pytest.importorskip("pandas")

# the clark2021 harness lives with its benchmark scripts (not shipped in the package)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "benchmarks" / "clark2021"))
import clark2021 as ck  # noqa: E402

DATA = Path(__file__).parent / "data" / "clark2021_excerpt.csv"


@pytest.fixture(scope="module")
def df():
    return ck.load(DATA)


# --- the medium ------------------------------------------------------------

def test_dm38_is_the_published_medium():
    dm38 = ck.dm38()
    # the four bulk carbon sources, at the concentrations in Supplementary Data 4
    assert dm38.initial_concentration("glc__D_e") == pytest.approx(24.98, abs=0.01)
    assert dm38.initial_concentration("arab__L_e") == pytest.approx(21.31, abs=0.01)
    assert dm38.initial_concentration("malt_e") == pytest.approx(4.38, abs=0.01)
    # anaerobic: the colon is, and so was the incubator
    assert dm38.initial_concentration("o2_e") == 0.0
    # batch culture, not a chemostat -- nothing is replenished
    assert all(dm38.influx_rate(m) == 0.0 for m in dm38.metabolites())


# --- namespace switching (CarveMe/BiGG vs gapseq/ModelSEED) -----------------

def test_metabolite_map_switches_namespace():
    assert ck.metabolite_map("bigg")["Butyrate"] == "but_e"
    assert ck.metabolite_map("modelseed")["Butyrate"] == "cpd00211_e0"
    assert ck.metabolite_map("modelseed")["Lactate"] == "cpd00159_e0"
    with pytest.raises(ValueError, match="bigg.*modelseed"):
        ck.metabolite_map("nonsense")


def test_dm38_modelseed_accessor_is_a_workstation_artifact_and_says_so_when_absent():
    """It is derived from the gapseq GEMs, so it is not committed; the accessor must not
    fail obscurely -- it points at the script that builds it."""
    from muode.diet import _DIET_DIR

    if (_DIET_DIR / "dm38_modelseed.csv").exists():
        pytest.skip("dm38_modelseed.csv is present (built on a workstation)")
    with pytest.raises(FileNotFoundError, match="derive_dm38_modelseed"):
        ck.medium("modelseed")


def test_the_measured_side_keys_by_the_chosen_namespace_and_subtracts_that_baseline(df):
    """A gapseq prediction is keyed cpd00211_e0; the measured truth it is scored against
    must be keyed the same way and have the SAME namespace's medium baseline removed.
    Pass the ModelSEED map + a ModelSEED medium and the observation flips namespace with
    the lactate-baseline correction intact."""
    from muode.diet import Diet

    ms = ck.metabolite_map("modelseed")
    # a minimal ModelSEED DM38: only the lactate baseline matters for net (the SCFAs are
    # products, absent from the medium -> baseline 0), and it is the 28.3 mM trap.
    ms_diet = Diet(concentrations={"cpd00159_e0": 28.3082}, name="dm38_ms_min")

    obs = ck.observations(df, diet=ms_diet, metabolites=ms)
    o = next(x for x in obs if x.net_metabolites)
    assert set(o.net_metabolites) <= set(ms.values())        # ModelSEED ids, not BiGG
    if "cpd00159_e0" in o.metabolites:                       # lactate baseline removed
        assert o.net_metabolites["cpd00159_e0"] < o.metabolites["cpd00159_e0"]


def test_dm38_feeds_a_real_bigg_model_anaerobically():
    """If DM38 cannot grow E. coli core, the mapping to BiGG dropped a nutrient."""
    cobra = pytest.importorskip("cobra")
    from muode.media import growth_on_diet

    core = cobra.io.load_model("textbook")
    assert growth_on_diet(core, ck.dm38()) > 0.05


def test_lactate_is_a_medium_component_not_only_a_product():
    """The trap: 28.3 mM lactate is supplied, so the endpoint is mostly medium."""
    assert ck.dm38().initial_concentration("lac__L_e") == pytest.approx(28.31, abs=0.01)


# --- loading ---------------------------------------------------------------

def test_contaminated_wells_are_excluded(df):
    raw = ck.load(DATA, drop_contaminated=False)
    assert (raw["Contamination?"] != "No").any()      # the fixture has some
    assert (df["Contamination?"] == "No").all()       # the loader drops them


def test_holdemanella_is_not_silently_dropped(df):
    """HB (Holdemanella biformis) is a 26th strain, distinct from BH (Blautia).

    Keying only on the authors' 25-species design space would discard every
    community containing HB -- silently, because the row simply would not match.
    """
    assert "HB" in ck.STRAINS and "BH" in ck.STRAINS
    assert ck.STRAINS["HB"].species == "Holdemanella biformis"
    assert ck.STRAINS["BH"].species == "Blautia hydrogenotrophica"
    assert "HB" not in ck.DESIGN_SPECIES          # ...but not in the design space
    assert "BH-HB" in set(df["Treatment"])        # and it survives loading


def test_named_pools_are_excluded_because_membership_is_unknown(df):
    assert not any(str(t).startswith("COMM") for t in df["Treatment"])


# --- observations ----------------------------------------------------------

def test_observations_report_net_production_against_the_dm38_baseline(df):
    obs = {o.community: o for o in ck.observations(df)}
    ac = obs["AC"]                                    # Anaerostipes caccae, monoculture

    assert ac.richness == 1 and ac.species == ("AC",)
    assert ac.n_replicates > 1
    # butyrate is absent from DM38, so endpoint == net production
    assert ac.metabolites["but_e"] == pytest.approx(ac.net_metabolites["but_e"])
    assert ac.net_metabolites["but_e"] > 10.0
    # lactate is supplied, so the two differ by exactly the medium concentration
    assert (ac.metabolites["lac__L_e"] - ac.net_metabolites["lac__L_e"]) == pytest.approx(28.31, abs=0.01)


# --- tier 1: the phenotype ground truth ------------------------------------

def test_monoculture_phenotypes_recover_the_known_butyrate_producers(df):
    """Derived from the measurements, cross-checked against the authors' own list.

    The one disagreement is FP, and it is the medium's doing, not an error --
    see the next test.
    """
    pheno = ck.monoculture_phenotypes(df)
    producers = {c for c, p in pheno.items() if p["but_e"]}

    assert producers == {"AC", "CC", "ER", "RI"}
    # every Bacteroides makes no butyrate; that is the sharp, falsifiable half
    assert not any(pheno[c]["but_e"] for c in ("BT", "BO", "BU", "BV", "BF", "BC"))
    # ...and they do make succinate
    assert all(pheno[c]["succ_e"] for c in ("BT", "BO", "BV", "BF", "BC", "PC"))

    # the authors' designation agrees except on FP
    assert set(ck.AUTHORS_BUTYRATE_PRODUCERS) - producers == {"FP"}


def test_f_prausnitzii_did_not_grow_in_dm38_so_it_must_be_excludable(df):
    """FP is a canonical butyrate producer that reaches OD ~0.02 in DM38 alone.

    A model that correctly grows FP would be scored *wrong* against its measured
    (zero) butyrate.  The benchmark has to be able to say so out loud.
    """
    assert ck.monoculture_growth(df)["FP"] < 0.1
    assert "FP" in ck.non_growers(df)

    truth = ck.monoculture_phenotypes(df)
    perfect = {c: dict(p) for c, p in truth.items()}
    perfect["FP"]["but_e"] = True                     # a model that gets FP right

    penalised = ck.score_phenotypes(perfect, truth, "but_e")
    assert penalised["false_positive"] == 1           # punished for being correct

    fair = ck.score_phenotypes(perfect, truth, "but_e", exclude=ck.non_growers(df))
    assert fair["false_positive"] == 0
    assert fair["accuracy"] == 1.0
    assert "FP" in fair["excluded"]                   # the exclusion travels with the score


def test_score_phenotypes_catches_a_bacteroides_that_secretes_butyrate(df):
    """The failure mode tier 1 exists to catch: a chemically impossible GEM."""
    truth = ck.monoculture_phenotypes(df)
    bad = {c: dict(p) for c, p in truth.items()}
    bad["BT"]["but_e"] = True                         # B. theta does not do this

    score = ck.score_phenotypes(bad, truth, "but_e")
    assert score["false_positive"] == 1
    assert score["accuracy"] < 1.0
    assert any("BT" in m for m in score["misclassified"])


# --- scoring metabolites across communities --------------------------------

def test_score_metabolites_is_perfect_on_a_perfect_prediction(df):
    obs = ck.observations(df)
    perfect = {o.community: dict(o.net_metabolites) for o in obs}

    score = ck.score_metabolites(perfect, obs, "but_e", net=True)
    assert score["n"] == len(obs)
    assert score["mae"] == pytest.approx(0.0, abs=1e-9)
    assert score["pearson_r"] == pytest.approx(1.0, abs=1e-6)
    assert score["bias"] == pytest.approx(0.0, abs=1e-9)


def test_scoring_lactate_on_the_endpoint_flatters_a_useless_model(df):
    """Why `net=True` is the default, in one assertion.

    A model that predicts *nothing happens to lactate* -- it just reports the
    medium back -- scores a tiny error on the raw endpoint, because ~95% of the
    measured lactate was there at t=0.  On net production, the same model is
    exposed.
    """
    obs = ck.observations(df)
    baseline = ck.dm38().initial_concentration("lac__L_e")
    useless = {o.community: {"lac__L_e": baseline} for o in obs}       # "nothing happened"

    raw = ck.score_metabolites(useless, obs, "lac__L_e", net=False)
    net = ck.score_metabolites(
        {c: {"lac__L_e": 0.0} for c in useless}, obs, "lac__L_e", net=True
    )
    # same model, same predictions -- the raw-endpoint framing hides the error
    assert raw["mae"] == pytest.approx(net["mae"])
    assert raw["observed_mean"] > 25.0        # dominated by the medium
    assert net["observed_mean"] < 10.0        # the real signal is much smaller
    # relative error is what changes: tiny against the endpoint, large against net
    assert raw["mae"] / raw["observed_mean"] < 0.5
    assert net["mae"] / abs(net["observed_mean"]) > 0.5


def test_score_metabolites_reports_directional_bias(df):
    obs = ck.observations(df)
    over = {o.community: {"but_e": o.net_metabolites["but_e"] + 5.0} for o in obs}

    score = ck.score_metabolites(over, obs, "but_e")
    assert score["bias"] == pytest.approx(5.0)        # +ve: over-predicts
    assert score["mae"] == pytest.approx(5.0)
    assert score["pearson_r"] == pytest.approx(1.0, abs=1e-6)   # correlated but biased
