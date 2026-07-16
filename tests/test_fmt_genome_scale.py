"""Step 1 of making the FMT scenario publishable: real GEMs instead of tuned dials.

The toy scenario's growth yields are hand-set (``_grower(CDIFF, yld=0.10)``), and those
dials are the root of the circularity pinned in test_provenance.py -- the competitive
outcome that clears the pathogen was tuned in.  ``examples/fmt_cdiff/genome_scale.py``
removes the dial: each member is a gapseq genome-scale reconstruction, so its yield is
biomass stoichiometry that came out of a genome.

These tests pin what is checkable NOW, while three of the five GEMs (including the
pathogen) are still reconstructing on the workstation.  They run on the two members that
are ready and skip cleanly for the rest, and they pin the two things that must not
regress: the community refuses to assemble without the pathogen, and the leftover
CarveMe/BiGG GEMs are never silently mixed into a ModelSEED community.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "fmt_cdiff"))

import genome_scale as gs  # noqa: E402


def _ready():
    return [m for m in gs.MEMBERS if gs.readiness()[m.slug] == "ready"]


def test_readiness_is_honest_about_what_is_and_is_not_reconstructed():
    """The report a human reads to decide whether the scenario can run yet.

    Every member is in exactly one of three states, and the states are not
    interchangeable: a leftover CarveMe/BiGG file is NOT 'ready' just because a file
    exists, because its ids would not cross-feed with the ModelSEED community.
    """
    status = gs.readiness()
    assert set(status) == {m.slug for m in gs.MEMBERS}
    assert set(status.values()) <= {
        "ready",
        "awaiting gapseq (have only old CarveMe/BiGG)",
        "not reconstructed yet",
    }


def test_the_community_refuses_to_assemble_without_the_pathogen(monkeypatch):
    """THE honesty gate.  A CDI scenario with no C. difficile is not a CDI scenario.

    build_community must raise rather than quietly hand back a donors-only community
    that would 'clear' an infection that was never there.  We simulate the pathogen
    being unavailable (rather than depending on whether its GEM happens to be present)
    so the guard stays covered now that the reconstruction exists.
    """
    pathogen = next(m for m in gs.MEMBERS if m.is_pathogen)
    real_load = gs.load
    monkeypatch.setattr(gs, "load",
                        lambda m: None if m.is_pathogen else real_load(m))

    with pytest.raises(RuntimeError, match="without C. difficile is not"):
        gs.build_community()

    # ...but you may inspect the donor community explicitly, with eyes open.
    donors = gs.build_community(require_pathogen=False)
    assert pathogen.slug not in donors.organism_ids


def test_the_scenario_uses_evidence_based_parameters_not_the_toy_overrides():
    """The genome-scale scenario must not re-introduce the invented bile numbers.

    The toy ``scenario.cdi_ecology`` overrode the germination/inhibition constants with
    the values that made bile look load-bearing (ki=0.05, etc.).  This scenario passes
    NO overrides, so its layers carry the provenance-registered defaults: germination
    km 15.9 mM and inhibition ki 0.5 mM (muode.provenance / test_provenance.py).  If
    someone re-adds an override to force a bile effect, this fails.
    """
    from muode.bile import KI_BILE
    from muode.lifecycle import KI_INHIBITOR, KM_GERMINANT

    layers = {type(layer).__name__: layer for layer in gs.cdi_ecology().layers}
    assert layers["SporeForming"].km_germinant == pytest.approx(KM_GERMINANT.value)
    assert layers["SporeForming"].ki_inhibitor == pytest.approx(KI_INHIBITOR.value)
    assert layers["BileAcidInhibition"].ki == pytest.approx(KI_BILE.value)


def test_guild_membership_is_the_real_strains_not_a_label():
    """bai / spore / antibiotic membership keys on the actual reconstructed members.

    C. scindens ATCC 35704 is the 7a-dehydroxylator; C. difficile is the spore-former
    and the vancomycin target.  These are genome-derived roles, not tunable knobs.
    """
    eco = gs.cdi_ecology()
    layers = {type(layer).__name__: layer for layer in eco.layers}
    assert layers["BileAcidTransform"].bai_producers == {"C_scindens_ATCC35704"}
    assert layers["SporeForming"].species == {gs.PATHOGEN}
    assert layers["Antibiotic"].susceptible == {gs.PATHOGEN}
    assert gs.PATHOGEN == "C_difficile_630"


def test_ablation_removes_exactly_the_named_layers():
    """The decomposition knob: ablate must drop the bile and/or pH arm and nothing else."""
    names = lambda eco: {type(l).__name__ for l in eco.layers}

    full = names(gs.cdi_ecology(""))
    assert {"WeakAcidInhibition", "BileAcidTransform", "BileAcidInhibition"} <= full

    assert "BileAcidTransform" not in names(gs.cdi_ecology("bile"))
    assert "BileAcidInhibition" not in names(gs.cdi_ecology("bile"))
    assert "WeakAcidInhibition" in names(gs.cdi_ecology("bile"))       # pH kept

    assert "WeakAcidInhibition" not in names(gs.cdi_ecology("ph"))
    assert "BileAcidTransform" in names(gs.cdi_ecology("ph"))          # bile kept

    competition = names(gs.cdi_ecology("bile ph"))
    assert "WeakAcidInhibition" not in competition
    assert "BileAcidTransform" not in competition
    # spore reservoir + antibiotic always remain -- they define the clinical setup
    assert {"SporeForming", "Antibiotic"} <= competition


def test_the_diet_carries_a_physiological_bile_pool():
    """The bile pool is added at ~caecal concentration, in muODE ecology ids.

    The GEMs are ModelSEED and have no bile exchanges; the bile pool is a layer concern.
    2 mM matches what reaches the caecum (Ramirez & Abel-Santos 2011) rather than a dial.
    """
    from muode.bile import CHOLATE, TAUROCHOLATE

    diet = gs.cdi_diet()
    assert diet.initial_concentration(TAUROCHOLATE) == pytest.approx(2.0)
    assert diet.initial_concentration(CHOLATE) == pytest.approx(2.0)
    # and it still carries the ModelSEED carbon the pathogen ferments (Stickland)
    assert diet.initial_concentration("cpd00129_e0") > 0    # L-proline
    assert diet.initial_concentration("cpd00033_e0") > 0    # glycine


@pytest.mark.slow
def test_the_old_carveme_gems_are_rejected_not_silently_used():
    """A file existing is not enough; it must be the right namespace.

    B. theta and C. scindens still have their CarveMe/BiGG carves in the tree.  Those
    load fine as cobra models, so the only thing standing between them and a corrupted
    community (BiGG ids that cross-feed with nothing) is the namespace check.  If a
    member reports 'ready', every exchange id it contributes must be ModelSEED.
    """
    any_bigg_pending = any(
        gs.readiness()[m.slug].startswith("awaiting") for m in gs.MEMBERS
    )
    if not any_bigg_pending:
        pytest.skip("no leftover CarveMe/BiGG members remain to reject")

    for m in gs.MEMBERS:
        org = gs.load(m)
        if org is None:
            continue                      # not ready -- correctly not loaded
        assert gs._is_modelseed(org), f"{m.slug} loaded but is not ModelSEED"


@pytest.mark.slow
def test_the_ready_producers_encode_a_connected_butyrate_pathway():
    """The claim that justified moving off CarveMe.

    CarveMe's Roseburia and F. prausnitzii carried NO butyrate exchange -- the SCFA arm
    could not come from stoichiometry.  Every ready member whose role is butyrate
    production must have a butyrate exchange that can carry export flux.  This is a
    property of the RECONSTRUCTION and is checked as capability (objective = export),
    so it holds even for a fastidious member that barely grows on the current medium.
    """
    producers = [m for m in _ready() if m.secretes == gs.BUTYRATE]
    if not producers:
        pytest.skip("no butyrate producer is reconstructed yet")

    diet = gs.load_modelseed_diet()
    for m in producers:
        r = gs.verify(m, diet)
        assert r.pathway_present, (
            f"{m.slug} claims to produce butyrate but the pathway carries no flux -- "
            "this is exactly the gap CarveMe left and gapseq was meant to fix"
        )
        assert r.can_secrete > 1e-6


@pytest.mark.slow
def test_the_genome_scale_scenario_integrates_and_the_antibiotic_bites():
    """A SHORT end-to-end run: the 5-GEM community integrates and behaves.

    Not the full study (that is ~30 min/arm on the workstation -- run_genome_scale.py).
    This is a 4 h smoke run that pins the things that must not silently break: all five
    members are tracked, the bile pool and spore reservoir exist, and the early
    vancomycin course actually knocks the vegetative pathogen down (the clinical setup).
    The FMT/clearance verdict is NOT asserted here -- that needs the full run.
    """
    if gs.readiness()[gs.PATHOGEN] != "ready":
        pytest.skip("pathogen not reconstructed yet")

    res = gs.build_scenario(fmt=True, t_end=4.0, dt=0.1)

    assert set(res.biomass.columns) == {m.slug for m in gs.MEMBERS}
    assert res.spores is not None and gs.PATHOGEN in res.spores.columns
    for bile_id in ("tca_e", "ca_e", "dca_e"):
        assert bile_id in res.metabolites.columns, f"bile pool {bile_id} not tracked"

    path = res.biomass[gs.PATHOGEN]
    assert path.iloc[0] == pytest.approx(0.05, abs=1e-6)
    assert path.min() < 0.01, "the vancomycin course must knock the vegetative pathogen down"


@pytest.mark.slow
def test_build_community_assembles_every_ready_member():
    if any(s != "ready" for s in gs.readiness().values()):
        pytest.skip("not all members reconstructed yet")
    comm = gs.build_community(abundances={m.slug: 1.0 for m in gs.MEMBERS})
    assert set(comm.organism_ids) == {m.slug for m in gs.MEMBERS}


@pytest.mark.slow
def test_ready_members_grow_at_a_physiological_rate_not_the_open_medium_fiction():
    """The yields are stoichiometry now, and they land in the physiological range.

    Growth is measured through the diet's flux bounds, so it must be a rate an organism
    could manage -- below the human-gut ceiling -- not the 8-10/h an open medium
    produces.  A member that barely grows (F. prausnitzii is fastidious) is fine and
    expected; a member ABOVE the ceiling would mean the flux bounds are not being
    applied, which is the uniform-Vmax bug this whole line of work exists to kill.
    """
    from muode.environment import HUMAN_GUT

    ready = _ready()
    if not ready:
        pytest.skip("no members reconstructed yet")

    diet = gs.load_modelseed_diet()
    for m in ready:
        r = gs.verify(m, diet)
        assert 0.0 <= r.growth < HUMAN_GUT.max_plausible_growth, (
            f"{m.slug} grows at {r.growth:.3f}/h -- above the gut ceiling means the diet's "
            "flux bounds are not biting (the uniform-Vmax bug)"
        )
