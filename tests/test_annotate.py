"""The trait-calling logic, pinned.

The searches shell out to HMMER and DIAMOND and are not tested here.  Everything
that carries a *judgement* -- how many markers make a call, what a hit means, what
an empty guild means -- is pure, and is tested.
"""

from __future__ import annotations

from muode import annotate
from muode.markers import BAI, BSH, TRAITS, MarkerSource, accessions, call_all


# -- parsing ----------------------------------------------------------------

def test_pfam_version_suffixes_are_stripped():
    """A registry pins a family (PF02275), not a release (PF02275.20)."""
    tbl = (
        "#  target  accession  query  accession\n"
        "prot_1 - cbah PF02275.20 1.2e-40 138.0\n"
    )
    assert annotate.parse_hmmsearch_tblout(tbl) == {"PF02275"}


def test_hmm_libraries_without_accessions_do_not_produce_a_phantom_hit():
    tbl = "prot_1 - somehmm - 1e-40 138.0\n"
    assert annotate.parse_hmmsearch_tblout(tbl) == set()


def test_diamond_subject_ids_resolve_whether_or_not_they_are_uniprot_formatted():
    tab = "q1\tsp|P19412|BAIE_CLOSD\t74.0\t250.0\nq2\tP19409\t68.0\t210.0\n"
    assert annotate.parse_diamond_tab(tab) == {"P19412", "P19409"}


def test_weak_alignments_are_not_evidence():
    """A short, low-identity hit to a promiscuous family must not call a trait."""
    tab = "q1\tsp|P19337|BAIA2_CLOSD\t22.0\t35.0\n"   # below both floors
    assert annotate.parse_diamond_tab(tab) == set()


# -- calling ----------------------------------------------------------------

def test_one_curated_pfam_hit_is_enough_for_bsh():
    """PF02275 is specific to the chemistry, so a single hit calls it."""
    assert BSH.call(frozenset({"PF02275"}))
    assert not BSH.call(frozenset())


def test_a_lone_dehydrogenase_hit_does_not_make_a_7alpha_dehydroxylator():
    """BaiA2 is an SDR.  Any SDR-rich genome will hit it; that is not the bai operon.

    This is the whole reason min_markers exists, and getting it wrong would invent
    a secondary-bile-acid producer out of an ordinary gut anaerobe.
    """
    assert not BAI.call(frozenset({"P19337"}))
    assert not BAI.call(frozenset({"P19337", "P19409"}))          # still only 2


def test_three_bile_specific_proteins_together_do_call_bai():
    assert BAI.call(frozenset({"P19409", "P19412", "P19337"}))


def test_bai_needs_more_evidence_than_bsh():
    """The evidence bar tracks how promiscuous the marker family is, not the trait."""
    assert BSH.min_markers == 1
    assert BAI.min_markers >= 3


# -- guilds -----------------------------------------------------------------

def test_guilds_invert_into_the_shape_the_ecology_layers_want():
    called = call_all({
        "mag_a": frozenset({"PF02275"}),                          # bsh only
        "mag_b": frozenset({"P19409", "P19412", "P19337"}),       # bai only
        "mag_c": frozenset(),                                     # neither
    })
    g = annotate.guilds(called)
    assert g["bsh"] == {"mag_a"}
    assert g["bai"] == {"mag_b"}


def test_a_community_with_no_bai_carrier_yields_an_empty_guild_not_an_error():
    """No bai operon means no deoxycholate.  That is a finding to report.

    It must not raise, must not fall back to a default guild, and must not be
    silently filled in -- all three would manufacture a producer that the genomes
    say is not there.
    """
    g = annotate.guilds(call_all({"mag_a": frozenset({"PF02275"})}))
    assert g["bai"] == set()
    assert g["bsh"] == {"mag_a"}


# -- the registry's own contract --------------------------------------------

def test_every_trait_says_why_the_gem_cannot_encode_it():
    """The rule for belonging in this module, enforced.

    If a trait cannot answer this, it is ordinary metabolism and belongs in the
    medium or the reconstruction -- as urease did.
    """
    for t in TRAITS:
        assert t.why_not_fba.strip(), f"{t.name} does not justify its existence"
        assert t.layer.strip(), f"{t.name} is read by no ecology layer"
        assert t.citation.strip(), f"{t.name} cites nothing"


def test_pfam_and_refprot_accessions_are_dispatched_to_different_tools():
    pfam, refprot = accessions(MarkerSource.PFAM), accessions(MarkerSource.REFPROT)
    assert "PF02275" in pfam and "P19412" in refprot
    assert not set(pfam) & set(refprot)


# -- round trip -------------------------------------------------------------

def test_traits_tsv_round_trips_and_keeps_the_columns_the_pipeline_reads(tmp_path):
    called = call_all({
        "mag_a": frozenset({"PF02275"}),
        "mag_b": frozenset({"P19409", "P19412", "P19337"}),
    })
    p = annotate.write_traits_tsv(tmp_path / "traits.tsv", called,
                                  domains={"mag_a": "bacteria"})
    header = [ln for ln in p.read_text().splitlines() if not ln.startswith("#")][0]
    # the Snakefile routes on these two; they must stay, and stay first
    assert header.split("\t")[:2] == ["mag_id", "domain"]

    back = annotate.read_traits_tsv(p)
    assert back["mag_a"]["bsh"] is True
    assert back["mag_a"]["bai"] is False
    assert back["mag_b"]["bai"] is True
