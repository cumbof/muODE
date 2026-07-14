"""The `source` column: what each number in a medium is actually worth.

A bound is not evidence just because it is a number.  Over half of `western_gut`
sits at the source table's *fill* value, and before this column existed nothing in
the file distinguished those rows from the ones derived from dietary intake data.
These tests pin the distinction so it cannot quietly rot back out.
"""

from __future__ import annotations

import pytest

from muode.diet import Diet, load_preset


@pytest.fixture(scope="module")
def western():
    return load_preset("western_gut")


def test_the_medium_says_where_every_bounded_number_came_from(western):
    for met in western.metabolites():
        assert western.provenance(met), f"{met} does not say where its number came from"


def test_the_fill_value_is_not_dressed_up_as_dietary_intake(western):
    """0.1 is the source's default.  Rows carrying it must NOT claim `intake`.

    This is the failure this column exists to prevent: 71 rows that look exactly
    as authoritative as the 44 that are.
    """
    for met in western.metabolites():
        limit = western.uptake_limit(met)
        if limit == 0.1 and western.provenance(met) == "intake":
            pytest.fail(f"{met}: bound is the source's fill value but claims to be intake")


def test_most_of_this_medium_is_not_a_measurement_and_says_so(western):
    """Guards the honesty of the headline, in both directions.

    If someone later relabels the defaults as `intake`, this fails.  If someone fixes
    them with real intake data, this also fails -- and that is a good failure, which
    should be celebrated by updating the number here.

    **It has now been celebrated once: 71/126 (56%) -> 35/113 (31%).**  The 36 rows that
    moved did so with a source and an arithmetic:

      * 13 vitamins and cofactors, bounded by Dietary Reference Intake (IOM 1998).
        Cobalamin at the fill value was the largest CARBON source in the medium.
      * 11 nucleosides and mucin amino sugars, bounded by dietary nucleic acid
        (0.1-1 g/day) and by the 2-3 g/day of mucin entering the large bowel.
      * 11 rows that are not food at all -- indole, H2S, the free nucleobases -- moved
        to PRODUCTS, plus chorismate and formaldehyde dropped outright.
      * nitrite, which the fill value supplied at 2.2 g/day against a real intake of
        1.2 mg/day, and which was the single most influential row in the file.

    What is LEFT is the honest remainder: amino acids, polyamines, glutathione, some
    organic acids, and the mineral rows (whose fill value is provably harmless -- see
    tests/test_nutrient_sensitivity.py).  Lowering this ceiling further needs a real
    ileal amino-acid table, which is a literature job and not a code one.
    """
    bounded = [m for m in western.metabolites() if western.uptake_limit(m) is not None]
    defaults = [m for m in bounded if western.provenance(m) == "source_default"]
    assert 0.2 < len(defaults) / len(bounded) < 0.4, (
        f"{len(defaults)}/{len(bounded)} rows are the source's default; if that "
        "changed, say so here rather than moving the goalposts"
    )


def test_the_amino_acid_inconsistency_is_visible_rather_than_hidden(western):
    """Phe sits 6.7x above Leu for no dietary reason.  The column has to show why."""
    assert western.provenance("phe__L_e") == "source_default"
    assert western.provenance("leu__L_e") == "intake"
    assert western.uptake_limit("phe__L_e") > 5 * western.uptake_limit("leu__L_e")


def test_products_are_never_supplied_and_never_capped(western):
    """The rule that keeps SCFA prediction falsifiable, now machine-checkable."""
    for met in western.metabolites():
        if western.provenance(met) != "product":
            continue
        assert western.initial_concentration(met) == 0.0, f"{met} is seeded"
        assert western.uptake_limit(met) is None, f"{met} has a dietary ceiling"


def test_methanol_is_no_longer_a_free_carbon_source(western):
    """4.8 M methanol would sterilise the vessel; it is an override, and marked one."""
    assert western.provenance("meoh_e") == "override"
    assert western.uptake_limit("meoh_e") < 0.1
    assert western.initial_concentration("meoh_e") < 100.0


def test_provenance_survives_a_round_trip(tmp_path):
    d = Diet(concentrations={"glc__D_e": 10.0}, max_uptake={"glc__D_e": 0.5},
             source={"glc__D_e": "intake"}, name="t")
    d.to_csv(tmp_path / "d.csv")
    assert Diet.from_csv(tmp_path / "d.csv").provenance("glc__D_e") == "intake"


def test_a_diet_without_the_column_still_loads(tmp_path):
    """Backwards compatibility: DM38 has no `source` column and must keep working."""
    p = tmp_path / "old.csv"
    p.write_text("metabolite,concentration,influx\nglc__D_e,10,0\n")
    d = Diet.from_csv(p)
    assert d.initial_concentration("glc__D_e") == 10.0
    assert d.provenance("glc__D_e") is None
