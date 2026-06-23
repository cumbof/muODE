"""Tests for Phase-1 reconstruction: the stub engine, command construction,
the MAG quality gate, and the reconstruction-at-scale QC summary.

The stub-model tests need cobra (they emit and re-read SBML); the command-
construction and summary tests are dependency-light.
"""

import json

import pytest


# ---------------------------------------------------------------------------
# Stub engine: produces simulatable placeholder GEMs and a cross-feeding pair
# ---------------------------------------------------------------------------

cobra = pytest.importorskip("cobra")

from muode.reconstruct import _stub_role, stub_reconstruct  # noqa: E402


def _write_mag(path, role):
    path.write_text(f">contig_1 [muode-stub:{role}] toy\nATGCATGCATGCATGCATGC\n")
    return path


def test_stub_role_parsing(tmp_path):
    assert _stub_role(_write_mag(tmp_path / "f.fna", "fermenter")) == "fermenter"
    assert _stub_role(_write_mag(tmp_path / "c.fna", "consumer")) == "consumer"
    # untagged genome defaults to generic
    plain = tmp_path / "p.fna"
    plain.write_text(">just_a_contig\nATGC\n")
    assert _stub_role(plain) == "generic"


def test_stub_fermenter_grows_and_secretes_acetate(tmp_path):
    out = stub_reconstruct(_write_mag(tmp_path / "ferm.fna", "fermenter"), tmp_path / "ferm.xml")
    model = cobra.io.read_sbml_model(str(out))
    assert {"EX_glc_e", "EX_ac_e"} <= {r.id for r in model.reactions}
    sol = model.optimize()
    assert sol.objective_value > 1e-6                 # grows
    assert sol.fluxes["EX_ac_e"] > 1e-6               # secretes acetate


def test_stub_consumer_cannot_use_glucose(tmp_path):
    out = stub_reconstruct(_write_mag(tmp_path / "cons.fna", "consumer"), tmp_path / "cons.xml")
    model = cobra.io.read_sbml_model(str(out))
    assert "EX_glc_e" not in {r.id for r in model.reactions}
    assert "EX_ac_e" in {r.id for r in model.reactions}
    # no acetate available on a closed medium -> no growth
    model.reactions.EX_ac_e.lower_bound = 0.0
    assert (model.slim_optimize() or 0.0) < 1e-6


def test_stub_pair_cross_feeds(tmp_path):
    """Two stub MAGs assemble into the canonical cross-feeding community."""
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.diet import Diet
    from muode.kinetics import KineticParameters
    from muode.organism import CobraOrganism

    ferm = stub_reconstruct(_write_mag(tmp_path / "ferm.fna", "fermenter"), tmp_path / "ferm.xml")
    cons = stub_reconstruct(_write_mag(tmp_path / "cons.fna", "consumer"), tmp_path / "cons.xml")
    comm = Community(
        [CobraOrganism.from_file(ferm, id="ferm"), CobraOrganism.from_file(cons, id="cons")],
        abundances={"ferm": 0.7, "cons": 0.3},
        total_biomass=0.02,
    )
    diet = Diet({"glc_e": 20.0, "ac_e": 0.0}, name="glc")
    kin = KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5), "ac_e": (10.0, 0.5)})
    res = DynamicFBA(t_end=12.0, dt=0.1).run(comm, diet, kin)

    # both grow, and the acetate consumer depends entirely on the fermenter
    assert res.biomass["ferm"].iloc[-1] > res.biomass["ferm"].iloc[0]
    assert res.biomass["cons"].iloc[-1] > res.biomass["cons"].iloc[0]
    assert not res.cross_feeding().empty


# ---------------------------------------------------------------------------
# CarveMe command construction (no binary needed -- the shell-out is mocked)
# ---------------------------------------------------------------------------


def test_carveme_command_construction(tmp_path, monkeypatch):
    import muode.reconstruct as recon

    captured = {}
    monkeypatch.setattr(recon, "require", lambda binary, env_hint="": binary)
    monkeypatch.setattr(recon, "run", lambda cmd, log=None: captured.setdefault("cmd", list(map(str, cmd))))

    recon.carveme(tmp_path / "p.faa", tmp_path / "m.xml", universe="grampos", gapfill_media="M9")
    cmd = captured["cmd"]
    assert cmd[0] == "carve"
    assert "-u" in cmd and "grampos" in cmd
    assert "-g" in cmd and "M9" in cmd


# ---------------------------------------------------------------------------
# MAG quality gate + reconstruction-at-scale summary
# ---------------------------------------------------------------------------


def test_filter_mags(tmp_path):
    from muode.qc import filter_mags

    report = tmp_path / "quality_report.tsv"
    report.write_text(
        "Name\tCompleteness\tContamination\n"
        "good_bin\t95.0\t1.0\n"
        "fragmented\t40.0\t1.0\n"
        "contaminated\t90.0\t25.0\n"
    )
    keep = filter_mags(report, min_completeness=50.0, max_contamination=10.0)
    assert keep == ["good_bin"]


def _qc_json(path, mag, grows, egc=False, added=0):
    path.write_text(json.dumps({
        "model": mag,
        "grew_initially": grows and added == 0,
        "reactions_added": list(range(added)),
        "grows_now": grows,
        "qc": {"n_mass_unbalanced": 0, "energy_generating_cycle": egc, "passed": not egc},
    }))
    return path


def test_summarize_reconstruction_and_filtering(tmp_path):
    from muode.qc import simulatable_mags, summarize_reconstruction

    paths = [
        _qc_json(tmp_path / "a.qc.json", "a", grows=True, added=3),
        _qc_json(tmp_path / "b.qc.json", "b", grows=False),            # dead -> not simulatable
        _qc_json(tmp_path / "c.qc.json", "c", grows=True, egc=True),   # energy cycle -> excluded
    ]
    df = summarize_reconstruction(paths)
    assert list(df["mag"]) == ["a", "b", "c"]
    assert df.set_index("mag").loc["a", "n_reactions_added"] == 3
    assert set(df.loc[df["simulatable"], "mag"]) == {"a"}

    summary = tmp_path / "reconstruction_summary.tsv"
    df.to_csv(summary, sep="\t", index=False)
    assert simulatable_mags(summary) == ["a"]
