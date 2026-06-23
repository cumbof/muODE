"""Upstream ingestion: MetaSBT taxonomy + Bracken abundance + the join.

MetaSBT *characterizes* taxonomy (closest species cluster per MAG); Bracken
gives quantitative abundance (from a Kraken2/Bracken DB built off the MetaSBT
database). These tests cover the parsers, the join, and the heavy-tool command
*builders* (constructed, never executed -- no kraken2/bracken/MetaSBT here).
"""

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# MetaSBT profile parser
# ---------------------------------------------------------------------------

_PROFILE = (
    "# level\tclosest\tani\tconfidence\n"
    "kingdom\tk__Bacteria\t0.0\t1.0\n"
    "phylum\tp__Firmicutes\t2.1\t0.95\n"
    "class\tc__Clostridia\t3.0\t0.9\n"
    "order\to__Eubacteriales\t4.0\t0.88\n"
    "family\tf__Lachnospiraceae\t5.0\t0.8\n"
    "genus\tg__Roseburia\t6.0\t0.75\n"
    "species\ts__Roseburia_intestinalis\t7.0\t0.7\n"
    "genome\tGCF_000209855\t1.2\t0.99\n"
)


def test_metasbt_profile_parses_levels(tmp_path):
    from muode.metasbt import read_metasbt_profile

    f = tmp_path / "MAG_1.txt"
    f.write_text(_PROFILE)
    prof = read_metasbt_profile(f)

    assert prof.genome_id == "MAG_1"
    assert prof.species == "s__Roseburia_intestinalis"
    assert prof.genus == "g__Roseburia"
    assert prof.closest_genome == "GCF_000209855"
    assert prof.confidence("species") == pytest.approx(0.7)
    assert prof.matches["phylum"].ani == pytest.approx(2.1)
    # lineage covers the 7 standard ranks (genome row excluded), coarse->fine
    assert list(prof.lineage) == [
        "kingdom", "phylum", "class", "order", "family", "genus", "species",
    ]
    assert prof.lineage_string.startswith("k__Bacteria|p__Firmicutes")


def test_metasbt_profile_split_suffix_and_no_header(tmp_path):
    from muode.metasbt import read_metasbt_profile

    # the .split.txt suffix is stripped from the genome id; header is optional
    f = tmp_path / "bin.42.split.txt"
    f.write_text("species\ts__Bacteroides_fragilis\t1.0\t0.9\n")
    prof = read_metasbt_profile(f)
    assert prof.genome_id == "bin.42"
    assert prof.species == "s__Bacteroides_fragilis"


def test_read_metasbt_profiles_directory(tmp_path):
    from muode.metasbt import mag_taxa, read_metasbt_profiles

    (tmp_path / "A.txt").write_text("species\ts__Species_one\t1\t0.9\n")
    (tmp_path / "B.txt").write_text("species\ts__Species_two\t1\t0.8\n")
    profiles = read_metasbt_profiles(tmp_path)
    assert set(profiles) == {"A", "B"}
    assert mag_taxa(profiles, level="species") == {
        "A": "s__Species_one", "B": "s__Species_two",
    }


# ---------------------------------------------------------------------------
# Bracken parser
# ---------------------------------------------------------------------------

_BRACKEN = (
    "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
    "Roseburia intestinalis\t166486\tS\t600\t100\t700\t0.70\n"
    "Bacteroides fragilis\t817\tS\t200\t100\t300\t0.30\n"
    "Clostridia\t186801\tC\t50\t0\t50\t0.05\n"
)


def test_read_bracken_species_fraction(tmp_path):
    from muode.bracken import read_bracken

    f = tmp_path / "sample.bracken.tsv"
    f.write_text(_BRACKEN)
    brk = read_bracken(f, level="S", value="fraction")
    # class-level row dropped; the two species renormalised to sum to 1
    assert set(brk.abundances) == {"Roseburia intestinalis", "Bacteroides fragilis"}
    assert brk.abundances["Roseburia intestinalis"] == pytest.approx(0.7)
    assert sum(brk.abundances.values()) == pytest.approx(1.0)
    assert brk.taxonomy_ids["Bacteroides fragilis"] == "817"


def test_read_bracken_reads_value_and_all_levels(tmp_path):
    from muode.bracken import read_bracken

    f = tmp_path / "sample.bracken.tsv"
    f.write_text(_BRACKEN)
    brk = read_bracken(f, level=None, value="reads", normalize=False)
    assert brk.abundances["Roseburia intestinalis"] == 700.0
    assert "Clostridia" in brk.abundances          # level=None keeps every row

    with pytest.raises(ValueError):
        read_bracken(f, value="bogus")


# ---------------------------------------------------------------------------
# The join: MetaSBT identity x Bracken abundance
# ---------------------------------------------------------------------------


def _profiles(species_by_mag):
    from muode.metasbt import MetaSBTProfile, ProfileMatch

    return {
        mag: MetaSBTProfile(mag, {"species": ProfileMatch("species", sp)})
        for mag, sp in species_by_mag.items()
    }


def test_join_matches_exact_and_normalised_labels():
    from muode.bracken import BrackenProfile
    from muode.quantify import mag_abundance_from_bracken, match_profiles_to_bracken

    brk = BrackenProfile({"Roseburia intestinalis": 0.7, "Bacteroides fragilis": 0.3})
    profiles = _profiles({
        "MAG_exact": "Roseburia intestinalis",          # exact name
        "MAG_prefixed": "s__Bacteroides_fragilis",       # rank prefix + underscores
        "MAG_missing": "s__Unknown_species",             # no Bracken taxon
    })
    matched = match_profiles_to_bracken(profiles, brk, level="species")
    assert matched["MAG_exact"] == "Roseburia intestinalis"
    assert matched["MAG_prefixed"] == "Bacteroides fragilis"
    assert matched["MAG_missing"] is None

    abund = mag_abundance_from_bracken(brk, profiles, level="species", normalize=True)
    assert abund["MAG_missing"] == 0.0
    assert abund["MAG_exact"] == pytest.approx(0.7)
    assert abund["MAG_prefixed"] == pytest.approx(0.3)
    assert sum(abund.values()) == pytest.approx(1.0)


def test_join_splits_shared_taxon_across_mags():
    from muode.bracken import BrackenProfile
    from muode.quantify import mag_abundance_from_bracken

    brk = BrackenProfile({"Escherichia coli": 0.8, "Bacteroides fragilis": 0.2})
    profiles = _profiles({
        "bin1": "Escherichia coli",     # two bins of the same species cluster
        "bin2": "Escherichia coli",
        "bin3": "Bacteroides fragilis",
    })
    abund = mag_abundance_from_bracken(brk, profiles, split_shared=True, normalize=False)
    assert abund["bin1"] == pytest.approx(0.4)   # 0.8 split between bin1 and bin2
    assert abund["bin2"] == pytest.approx(0.4)
    assert abund["bin3"] == pytest.approx(0.2)

    whole = mag_abundance_from_bracken(brk, profiles, split_shared=False, normalize=False)
    assert whole["bin1"] == pytest.approx(0.8)   # each bin inherits the full taxon
    assert whole["bin2"] == pytest.approx(0.8)


# ---------------------------------------------------------------------------
# Heavy-tool command builders (constructed only -- never executed here)
# ---------------------------------------------------------------------------


def test_command_builders_are_pure():
    from muode import quantify as q

    prof = q.metasbt_profile_cmd("genomes.txt", "/db/index", "/out/profiles")
    assert prof[:2] == ["MetaSBT", "profile"]
    assert "--output-dir" in prof and "/out/profiles" in prof
    assert "--expand" in prof

    krk = q.metasbt_kraken_cmd("/db", "/out/kraken_db", genomes="g.txt",
                               names_dmp="names.dmp", nodes_dmp="nodes.dmp")
    assert krk[:2] == ["MetaSBT", "kraken"]
    assert "--genomes" in krk and "g.txt" in krk

    bb = q.bracken_build_cmd("/out/kraken_db", read_len=100)
    assert bb[0] == "bracken-build" and "100" in bb

    k2 = q.kraken2_cmd("/out/kraken_db", ["r1.fq", "r2.fq"], "rep.txt", paired=True)
    assert k2[0] == "kraken2" and "--paired" in k2
    assert k2[-2:] == ["r1.fq", "r2.fq"]            # reads come last

    bk = q.bracken_cmd("/out/kraken_db", "rep.txt", "out.bracken", level="S")
    assert bk[0] == "bracken" and "out.bracken" in bk


def test_runner_requires_missing_tool(monkeypatch):
    """The runners locate the binary first, so a missing tool fails loudly."""
    from muode import quantify as q
    from muode.external import MuodeToolError

    monkeypatch.setattr(q, "require", lambda *a, **k: (_ for _ in ()).throw(
        MuodeToolError("required tool 'MetaSBT' not found on PATH")))
    with pytest.raises(MuodeToolError):
        q.metasbt_kraken("/db", "/out")
