"""Strain matching for the Clark et al. benchmark genomes.

Reconstructing the wrong strain is the quietest possible way to invalidate the
whole benchmark: metabolic phenotype is strain-specific, every downstream number
still computes, and nothing looks broken.  These tests pin the two ways it nearly
happened.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "benchmarks" / "clark2021"))

from fetch_genomes import _norm, strains_match  # noqa: E402

from muode.benchmarks.clark2021 import STRAINS  # noqa: E402


def test_strain_designations_normalise_across_punctuation():
    assert _norm("A2-165") == _norm("A2_165") == _norm("A2 165") == "a2165"
    assert _norm("VPI-5482") == _norm("VPI 5482")


def test_matching_is_exact_or_confidently_fuzzy():
    assert strains_match(_norm("VPI-5482"), _norm("VPI 5482"))          # punctuation
    assert strains_match(_norm("DSM 14662"), _norm("DSM 14662 4"))      # suffixed
    assert not strains_match(_norm("ATCC 8482"), _norm("ATCC 8483"))    # near-miss
    assert not strains_match("", _norm("A2-165"))                       # no evidence


def test_a_short_strain_token_never_fuzzy_matches():
    """The bug this caught: NCBI reports F. prausnitzii A2-165's strain as "1".

    A naive substring test finds "1" inside "a2165", declares a confident match,
    and hands back a genome from a different organism -- with no warning.
    """
    assert not strains_match(_norm("A2-165"), _norm("1"))
    assert not strains_match(_norm("L1-82"), _norm("82"))
    # ...but a long, genuinely discriminating token still matches
    assert strains_match(_norm("A2-165"), _norm("A2165"))


def test_reclassified_organisms_carry_an_ncbi_taxon_override():
    """F. prausnitzii A2-165 is now *Faecalibacterium duncaniae* A2-165.

    Searching NCBI under the name in the 2021 paper silently returns the species
    reference (M21/2) instead -- a different strain, with different metabolism.
    """
    fp = STRAINS["FP"]
    assert fp.strain == "A2-165"
    assert fp.taxon == "Faecalibacterium duncaniae"     # what NCBI calls it now
    assert fp.species == "Faecalibacterium prausnitzii"  # what the paper called it


def test_taxon_defaults_to_species_when_not_reclassified():
    bt = STRAINS["BT"]
    assert bt.ncbi_taxon == "" and bt.taxon == bt.species


def test_every_strain_has_a_designation_to_match_on():
    """A blank strain would silently fall back to the species reference genome."""
    assert len(STRAINS) == 26
    assert all(s.strain.strip() for s in STRAINS.values())
