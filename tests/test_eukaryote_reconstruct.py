"""Tests for the eukaryote reconstruction routes.

The external tools (MetaEuk, CarveFungi, eggNOG-mapper, ModelSEEDpy) are not
available here, so the shell-outs are mocked.  We test what muODE owns: command
construction, the gene-call -> reconstruct dispatch by engine, ref-db guarding,
and the eggNOG annotation parser.
"""

import pytest

import muode.reconstruct as recon
from muode.eukaryote_draft import parse_eggnog_kos
from muode.external import MuodeToolError


# ---------------------------------------------------------------------------
# MetaEuk eukaryote gene calling
# ---------------------------------------------------------------------------


def test_metaeuk_command_construction(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, log=None):
        captured["cmd"] = list(map(str, cmd))
        # MetaEuk writes <prefix>.fas; emulate that so the wrapper can rename it.
        prefix = cmd[cmd.index("easy-predict") + 3]
        from pathlib import Path
        Path(f"{prefix}.fas").write_text(">g1\nMAAA\n")

    monkeypatch.setattr(recon, "require", lambda binary, env_hint="": binary)
    monkeypatch.setattr(recon, "run", fake_run)

    out = recon.call_genes_eukaryote(
        tmp_path / "euk.fna", tmp_path / "euk.faa", tmp_path / "refdb", threads=4
    )
    cmd = captured["cmd"]
    assert cmd[0] == "metaeuk" and cmd[1] == "easy-predict"
    assert "--threads" in cmd and "4" in cmd
    assert out.read_text().startswith(">g1")


# ---------------------------------------------------------------------------
# CarveFungi command construction
# ---------------------------------------------------------------------------


def test_carvefungi_default_command(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, log=None):
        captured["cmd"] = list(map(str, cmd))
        out = cmd[cmd.index("-o") + 1]
        from pathlib import Path
        Path(out).write_text("<sbml/>")

    monkeypatch.setattr(recon, "require", lambda binary, env_hint="": binary)
    monkeypatch.setattr(recon, "run", fake_run)

    recon.carvefungi(tmp_path / "p.faa", tmp_path / "m.xml", threads=2)
    cmd = captured["cmd"]
    assert cmd[0] == "carvefungi"
    assert "-o" in cmd and "--threads" in cmd and "2" in cmd


def test_carvefungi_cmd_template_override(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, log=None):
        captured["cmd"] = list(map(str, cmd))
        from pathlib import Path
        Path(tmp_path / "m.xml").write_text("<sbml/>")

    # require() must NOT be called when a template is supplied.
    monkeypatch.setattr(recon, "require", lambda *a, **k: pytest.fail("require() should not run"))
    monkeypatch.setattr(recon, "run", fake_run)

    recon.carvefungi(
        tmp_path / "p.faa", tmp_path / "m.xml",
        cmd_template="myfungi --in {proteins} --out {output} -t {threads}",
    )
    cmd = captured["cmd"]
    assert cmd[0] == "myfungi"
    assert str(tmp_path / "p.faa") in cmd and str(tmp_path / "m.xml") in cmd


# ---------------------------------------------------------------------------
# Engine dispatch in reconstruct_mag
# ---------------------------------------------------------------------------


def test_reconstruct_mag_routes_fungus_to_carvefungi(tmp_path, monkeypatch):
    calls = {}
    monkeypatch.setattr(recon, "call_genes_eukaryote",
                        lambda g, p, db, threads=1: calls.setdefault("genes", (str(g), str(db))) or p)
    monkeypatch.setattr(recon, "carvefungi",
                        lambda p, o, threads=1, cmd_template=None: calls.setdefault("carvefungi", True) or o)
    monkeypatch.setattr(recon, "eukaryote_generic",
                        lambda *a, **k: pytest.fail("generic should not run for carvefungi"))

    recon.reconstruct_mag(tmp_path / "f.fna", tmp_path / "f.xml",
                          engine="carvefungi", ref_db=tmp_path / "db")
    assert calls.get("carvefungi") is True
    assert calls["genes"][1].endswith("db")


def test_reconstruct_mag_routes_other_euk_to_generic(tmp_path, monkeypatch):
    calls = {}
    monkeypatch.setattr(recon, "call_genes_eukaryote", lambda g, p, db, threads=1: p)
    monkeypatch.setattr(recon, "eukaryote_generic",
                        lambda p, o, eggnog_data_dir=None, threads=1: calls.setdefault("generic", True) or o)
    monkeypatch.setattr(recon, "carvefungi",
                        lambda *a, **k: pytest.fail("carvefungi should not run for generic"))

    recon.reconstruct_mag(tmp_path / "e.fna", tmp_path / "e.xml",
                          engine="eukaryote_generic", ref_db=tmp_path / "db")
    assert calls.get("generic") is True


def test_reconstruct_mag_euk_requires_ref_db(tmp_path):
    with pytest.raises(MuodeToolError, match="reference database"):
        recon.reconstruct_mag(tmp_path / "f.fna", tmp_path / "f.xml", engine="carvefungi")


def test_reconstruct_mag_unknown_engine(tmp_path):
    with pytest.raises(ValueError, match="unknown reconstruction engine"):
        recon.reconstruct_mag(tmp_path / "f.fna", tmp_path / "f.xml", engine="bogus")


# ---------------------------------------------------------------------------
# eggNOG annotation parsing
# ---------------------------------------------------------------------------


def test_parse_eggnog_kos(tmp_path):
    ann = tmp_path / "x.emapper.annotations"
    ann.write_text(
        "## emapper version blah\n"
        "#query\tseed_ortholog\tKEGG_ko\tEC\n"
        "gene1\tortho1\tko:K00001,ko:K00002\t1.1.1.1\n"
        "gene2\tortho2\t-\t-\n"
        "gene3\tortho3\tko:K12345\t-\n"
        "## done\n"
    )
    ko_map = parse_eggnog_kos(ann)
    assert ko_map == {"gene1": ["K00001", "K00002"], "gene3": ["K12345"]}


def test_parse_eggnog_kos_empty(tmp_path):
    ann = tmp_path / "empty.emapper.annotations"
    ann.write_text(
        "#query\tseed_ortholog\tKEGG_ko\tEC\n"
        "geneA\to\t-\t-\n"
    )
    assert parse_eggnog_kos(ann) == {}
