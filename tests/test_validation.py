"""M3 -- validation & scale: MetaSBT ingestion, abundance subsampling, and the
benchmark/validation framework.

These are dependency-light (numpy/pandas/scipy only); the end-to-end validation
test reuses the bundled toy cross-feeding community.
"""

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# MetaSBT profile ingestion
# ---------------------------------------------------------------------------


def test_metasbt_autodetect_and_normalise(tmp_path):
    from muode.metasbt import read_metasbt_profile

    prof = tmp_path / "profile.tsv"
    prof.write_text(
        "genome\ttaxonomy\trelative_abundance\n"
        "MAG_1\tk__Bacteria|p__Firmicutes\t0.6\n"
        "MAG_2\tk__Bacteria|p__Bacteroidetes\t0.4\n"
    )
    p = read_metasbt_profile(prof)
    assert p.abundances == {"MAG_1": pytest.approx(0.6), "MAG_2": pytest.approx(0.4)}
    assert p.taxonomy["MAG_1"].startswith("k__Bacteria")
    assert p.abundance_column == "relative_abundance"


def test_metasbt_counts_normalised_and_override(tmp_path):
    from muode.metasbt import read_metasbt_profile

    prof = tmp_path / "counts.tsv"
    prof.write_text("bin\treads\nA\t30\nB\t10\n")     # raw counts -> relative
    p = read_metasbt_profile(prof)
    assert p.abundances == {"A": pytest.approx(0.75), "B": pytest.approx(0.25)}

    # explicit columns + no normalisation
    p2 = read_metasbt_profile(prof, id_column="bin", abundance_column="reads", normalize=False)
    assert p2.abundances == {"A": 30.0, "B": 10.0}


# ---------------------------------------------------------------------------
# Abundance-aware subsampling
# ---------------------------------------------------------------------------


def test_select_by_abundance_policies():
    from muode.subsample import select_by_abundance

    ab = {"a": 50, "b": 30, "c": 15, "d": 5}     # rel: .50 .30 .15 .05
    assert select_by_abundance(ab, top_n=2) == ["a", "b"]
    assert set(select_by_abundance(ab, min_abundance=0.1)) == {"a", "b", "c"}
    assert select_by_abundance(ab, coverage=0.8) == ["a", "b"]      # .50+.30 = .80


def test_community_subsample_renormalises():
    from muode.community import Community
    from muode.examples import build_acetate_specialist, build_glucose_specialist

    comm = Community(
        [build_glucose_specialist(), build_acetate_specialist()],
        abundances={"A_glucose": 0.9, "B_acetate": 0.1},
    )
    small = comm.subsample(top_n=1)
    assert small.organism_ids == ["A_glucose"]
    assert small.abundances["A_glucose"] == pytest.approx(1.0)   # renormalised

    with pytest.raises(ValueError):
        comm.subsample(min_abundance=2.0)                        # removes everything


# ---------------------------------------------------------------------------
# Validation framework
# ---------------------------------------------------------------------------


def test_compare_passes_and_fails():
    from muode.validate import BenchmarkExpectation, compare

    exp = BenchmarkExpectation(
        relative_abundances={"A": 0.4, "B": 0.6},
        metabolites={"glc_e": 0.0},
        cross_feeding=[("A", "ac_e", "B")],
    )
    good = compare({"A": 4.0, "B": 6.0}, {"glc_e": 0.0}, [("A", "ac_e", "B")], exp)
    assert good.passed
    assert good.metrics["cross_feeding"]["f1"] == 1.0

    bad = compare({"A": 9.0, "B": 1.0}, {"glc_e": 0.0}, [], exp)   # wrong composition + missed edge
    assert not bad.passed
    assert bad.metrics["relative_abundance"]["passed"] is False
    assert bad.metrics["cross_feeding"]["recall"] == 0.0


def test_validate_toy_community_against_benchmark():
    from muode.dfba import DynamicFBA
    from muode.examples import build_toy_community, toy_diet, toy_kinetics
    from muode.validate import BenchmarkExpectation, validate

    res = DynamicFBA(t_end=24.0, dt=0.05).run(build_toy_community(), toy_diet(), toy_kinetics())
    exp = BenchmarkExpectation.from_file(REPO / "examples/benchmarks/toy_cross_feeding.yaml")
    rep = validate(res, exp)
    assert rep.passed, rep.metrics
    assert rep.metrics["cross_feeding"]["f1"] == 1.0


def test_validate_outputs_roundtrip(tmp_path):
    from muode.dfba import DynamicFBA
    from muode.examples import build_toy_community, toy_diet, toy_kinetics
    from muode.validate import BenchmarkExpectation, validate_outputs

    res = DynamicFBA(t_end=24.0, dt=0.05).run(build_toy_community(), toy_diet(), toy_kinetics())
    res.to_csv(tmp_path)
    exp = BenchmarkExpectation.from_file(REPO / "examples/benchmarks/toy_cross_feeding.yaml")
    rep = validate_outputs(tmp_path, exp)
    assert rep.passed
