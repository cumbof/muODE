"""Traits / domain metadata and reconstruction routing."""

from muode.community import Community
from muode.organism import LinprogOrganism
from muode.traits import (
    Domain,
    MicrobeTraits,
    OxygenTolerance,
    reconstruction_route,
)


def _grower(id):
    return LinprogOrganism(
        id=id,
        reactions=[("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0),
                   ("GROW", {"glc_e": -1.0}, 0.0, 1000.0)],
        objective={"GROW": 0.1},
        exchanges={"glc_e": "EX_glc_e"},
    )


def test_reconstruction_routing_per_domain():
    assert reconstruction_route(Domain.BACTERIA)[0] == "carveme"
    assert reconstruction_route(Domain.ARCHAEA)[0] == "carveme"
    # eukaryotes and viruses have no CarveMe route -- flagged, not silently carved
    assert reconstruction_route(Domain.EUKARYOTE)[0] is None
    assert "CarveMe does NOT" in reconstruction_route(Domain.EUKARYOTE)[1]
    assert reconstruction_route(Domain.VIRUS)[0] is None
    assert "PhageInfection" in reconstruction_route(Domain.VIRUS)[1]


def test_route_accepts_plain_string_domain():
    # Domain is a str-enum, so a bare string round-trips
    assert reconstruction_route("eukaryote")[0] is None


def test_virus_is_not_metabolic():
    assert MicrobeTraits(domain=Domain.VIRUS).is_metabolic is False
    assert MicrobeTraits(domain=Domain.EUKARYOTE).is_metabolic is True
    assert MicrobeTraits().is_metabolic is True  # default bacterium


def test_community_traits_default_and_lookup():
    comm = Community(
        [_grower("a"), _grower("fungus")],
        {"a": 0.5, "fungus": 0.5},
        traits={"fungus": MicrobeTraits(domain=Domain.EUKARYOTE,
                                        oxygen=OxygenTolerance.FACULTATIVE)},
    )
    assert comm.domain_of("fungus") is Domain.EUKARYOTE
    assert comm.traits_of("fungus").oxygen is OxygenTolerance.FACULTATIVE
    # a member without an explicit entry defaults to an anaerobic bacterium
    assert comm.domain_of("a") is Domain.BACTERIA
    assert comm.traits_of("a").oxygen is OxygenTolerance.OBLIGATE_ANAEROBE


def test_subsample_preserves_traits():
    comm = Community(
        [_grower("a"), _grower("b")],
        {"a": 0.9, "b": 0.1},
        traits={"a": MicrobeTraits(domain=Domain.EUKARYOTE)},
    )
    sub = comm.subsample(top_n=1)
    assert sub.organism_ids == ["a"]
    assert sub.domain_of("a") is Domain.EUKARYOTE
