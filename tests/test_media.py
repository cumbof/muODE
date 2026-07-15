"""The medium bridge: a Diet -> a CarveMe media db, and a Diet -> cobra bounds.

The failure these guard against is the quiet one.  A CarveMe model carved without
a gap-fill medium grows only with *every* exchange open; the engine then closes
everything the diet does not name, the model cannot grow, and the run returns a
flat line that looks like a successful simulation.  So: the medium we carve
against and the medium we simulate on must be the same medium, and the QC growth
check must be run under the diet's bounds rather than the model's own.
"""

import cobra
import pytest

from muode.diet import Diet, load_diet, load_preset
from muode.gapfill import ensure_biomass, grows
from muode.media import (
    diet_medium,
    growth_on_diet,
    supplied_metabolites,
    to_compound,
    write_carveme_mediadb,
)


@pytest.fixture(scope="module")
def core():
    """E. coli core -- a real BiGG model, small enough to solve in milliseconds."""
    return cobra.io.load_model("textbook")


@pytest.fixture(scope="module")
def genome_scale(ijo):
    """iJO1366, from the session fixture -- loaded once for the whole suite.

    Needed because E. coli core's 20 exchanges physically *cannot* import enough
    nutrients to expose a diet's uptake bounds.  The 6.5/h growth rates in the real
    89-MAG run were invisible to every core-based test for exactly this reason.
    """
    return ijo


# --- Diet -> CarveMe media db ---------------------------------------------

def test_compound_ids_drop_the_compartment_suffix():
    assert to_compound("glc__D_e") == "glc__D"
    assert to_compound("h2o_e") == "h2o"
    assert to_compound("nh4") == "nh4"          # already bare: unchanged


def test_mediadb_lists_only_what_the_diet_actually_supplies(tmp_path):
    diet = Diet(
        concentrations={"glc__D_e": 10.0, "nh4_e": 5.0, "ac_e": 0.0, "but_e": 0.0},
        influx={"pi_e": 1.0},          # supplied by influx alone, starts at 0
        name="t",
    )
    path = write_carveme_mediadb(diet, tmp_path / "mediadb.tsv")
    lines = path.read_text().strip().split("\n")

    assert lines[0].split("\t") == ["medium", "description", "compound"]
    compounds = {ln.split("\t")[2] for ln in lines[1:]}
    assert compounds == {"glc__D", "nh4", "pi"}
    # acetate/butyrate start at 0 with no influx: they are cross-fed by the
    # community, not provided -- carving against them would be a lie.
    assert "ac" not in compounds and "but" not in compounds
    assert all(ln.startswith("t\t") for ln in lines[1:])


def test_mediadb_of_the_real_diet_is_carveme_shaped(tmp_path):
    diet = load_preset("western_gut")
    path = write_carveme_mediadb(diet, tmp_path / "mediadb.tsv", medium="western_gut")
    rows = [ln.split("\t") for ln in path.read_text().strip().split("\n")[1:]]

    assert len(rows) == len(supplied_metabolites(diet))
    assert all(len(r) == 3 and r[0] == "western_gut" for r in rows)
    assert all(not r[2].endswith("_e") for r in rows)      # no compartment suffix
    assert "fru" in {r[2] for r in rows}
    assert "o2" not in {r[2] for r in rows}                # the gut medium is anaerobic


# --- Diet -> cobra medium --------------------------------------------------

def test_diet_medium_closes_every_exchange_the_diet_does_not_name(core):
    diet = Diet(concentrations={"glc__D_e": 10.0}, name="glc_only")
    medium = diet_medium(core, diet)

    assert medium["EX_glc__D_e"] > 0.0
    # everything else is shut: this is exactly what the engine does, and why an
    # un-gap-filled model dies on a defined diet
    assert medium["EX_nh4_e"] == 0.0
    assert medium["EX_o2_e"] == 0.0
    assert set(medium) == {ex.id for ex in core.exchanges}


@pytest.mark.slow
def test_the_colonic_medium_feeds_a_gut_anaerobe_and_starves_the_organisms_that_are_not_one(bigg):
    """The medium must be lean, but it must not be too lean for a GUT organism.

    This test used to run on iJO1366 and assert that a *genome-scale* E. coli could
    gather enough from the medium's long tail of minor nutrients to grow.  It could --
    and the tail was phantom.  Nine nucleoside rows were sitting at the source table's
    fill value of 0.1, feeding it roughly 9 mmol C/gDW/h of RNA that nobody eats.  Bound
    them by real dietary nucleic acid and iJO1366 cannot pay its ATP maintenance at all.

    Which raises the only question that matters: is the honest medium too poor for
    *everything*, or just for E. coli?  So ask a real gut anaerobe.

        iCN900  (Clostridioides difficile 630, a gut anaerobe)   0.0036/h
        iJO1366 (Escherichia coli K-12)                          0
        iYO844  (Bacillus subtilis -- a soil organism)           0

    The colonic medium feeds the colonic organism.  E. coli ferments sugar, and a
    colon's sugar is locked inside fibre polymers it has no enzymes to open; in a real
    gut it is a cross-feeder, living on what the primary degraders release.  That is
    exactly what test_e_coli_core_cannot_grow_on_the_colonic_medium (below) has always
    said -- and the genome-scale model turns out to be in the same position, which we
    could not see while the fill value was feeding it.

    HONEST CAVEAT, because it weakens the comparison: iCN900 carries no ATP maintenance
    demand (ATPM = 0) while iJO1366 demands 3.15, so this is not a clean
    metabolism-only contrast.  What it does establish -- which is what the test is for
    -- is that the medium is not too lean to support growth at all.
    """
    from muode.qc import MAX_PLAUSIBLE_GROWTH

    anaerobe = bigg("iCN900")           # NB: downloaded from BiGG -- hence `slow`
    mu = growth_on_diet(anaerobe, load_preset("western_gut"))
    assert mu > 1e-3, "a gut anaerobe must be able to grow on a colonic medium"
    assert mu < MAX_PLAUSIBLE_GROWTH, "...but at a rate an organism could actually manage"


@pytest.mark.slow
def test_a_genome_scale_e_coli_also_starves_on_the_colonic_medium(genome_scale):
    """The other half of the pair, pinned so nobody quietly re-feeds it.

    If this starts failing -- if iJO1366 grows on western_gut again -- then carbon has
    been put back into the medium that is not in a colon.  Check what, before believing
    any downstream result.
    """
    assert growth_on_diet(genome_scale, load_preset("western_gut")) == pytest.approx(
        0.0, abs=1e-6
    )


def test_e_coli_core_cannot_grow_on_the_colonic_medium(core):
    """Not a bug -- the finding that the flux bounds are doing their job.

    E. coli core has 20 exchanges and no way to eat a polysaccharide, so on a
    colonic medium it can reach only fructose plus two amino acids: nowhere near
    enough to pay its 8.39 mmol/gDW/h ATP maintenance.  Which is exactly the
    position E. coli is in in the real colon, where it survives on sugars released
    by primary degraders.  If this ever starts passing, someone has put free
    glucose back into the large intestine.
    """
    assert growth_on_diet(core, load_preset("western_gut")) == pytest.approx(0.0, abs=1e-6)
    assert core.reactions.ATPM.lower_bound > 8.0


def test_growth_on_diet_does_not_mutate_the_model(core):
    before = {ex.id: ex.bounds for ex in core.exchanges}
    growth_on_diet(core, load_preset("western_gut"))
    assert {ex.id: ex.bounds for ex in core.exchanges} == before


def test_model_starves_when_the_diet_lacks_an_element(core):
    """Drop phosphate and a genome-scale biomass reaction cannot fire."""
    diet = load_preset("western_gut")
    starved = Diet(
        {m: c for m, c in diet.concentrations.items() if m != "pi_e"},
        {m: f for m, f in diet.influx.items() if m != "pi_e"},
        name="no_phosphate",
    )
    assert growth_on_diet(core, starved) == pytest.approx(0.0, abs=1e-9)


# --- the QC hole this closes ----------------------------------------------

def test_grows_on_open_medium_but_not_on_the_diet(core):
    """The exact silent failure: QC says fine, the engine says flat.

    `grows()` without a diet asks the model's *stored* medium, which for a
    CarveMe model is wide open -- a check that cannot fail.  With the diet it
    asks the question the simulation will ask.
    """
    no_p = Diet({"glc__D_e": 10.0}, name="carbon_only")    # C but no N/P/S
    assert grows(core) is True                             # "healthy" on its own medium
    assert grows(core, diet=no_p) is False                 # dead on the medium we simulate


def test_ensure_biomass_reports_both_media(core):
    report = ensure_biomass(core, diet=Diet({"glc__D_e": 10.0}, name="carbon_only"))

    assert report["growth_on_complete_medium"] > 0.1
    assert report["growth_on_diet"] == pytest.approx(0.0, abs=1e-9)
    # structurally fine, so it still reaches the community: it may yet be fed by
    # cross-feeding.  The flat-diet growth is the diagnostic, not the gate.
    assert report["grows_now"] is True
    assert report["diet"] == "carbon_only"


def test_ensure_biomass_without_a_diet_is_unchanged(core):
    report = ensure_biomass(core)
    assert report["grows_now"] is True
    assert "growth_on_diet" not in report


# --- diet resolution -------------------------------------------------------

def test_load_diet_accepts_a_preset_or_a_csv(tmp_path):
    assert load_diet("western_gut").name == "western_gut"

    csv = tmp_path / "mine.csv"
    csv.write_text("metabolite,concentration,influx\nglc__D_e,5.0,0.5\n")
    d = load_diet(str(csv))
    assert d.name == "mine" and d.initial_concentration("glc__D_e") == 5.0

    with pytest.raises(FileNotFoundError):
        load_diet(str(tmp_path / "nope.csv"))
    with pytest.raises(KeyError):
        load_diet("not_a_preset")
