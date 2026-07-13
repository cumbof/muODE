#!/usr/bin/env python3
"""Derive the western gut diet as a muODE diet CSV, with *published* flux bounds.

Why this script exists
----------------------
muODE's original `western_gut.csv` was hand-curated: it listed metabolites and
concentrations, and the engine turned those into uptake bounds via a *uniform*
default Vmax of 10 mmol/gDW/h.  A genome-scale model on an 89-metabolite diet
therefore imported 89 nutrients at 10 mmol/gDW/h each -- no cell has ever been in
such a medium.  In a real 89-MAG run this produced growth rates up to **6.5/h**: a
6.4-minute doubling time, faster than any organism ever measured.

Published diets avoid this by specifying a **flux bound per metabolite**, derived
from actual dietary intake.  In this one the median bound is **0.1 mmol/gDW/h**,
not 10.  That is the entire fix.

    python examples/diets/derive_western_gut.py

Source
------
`western_diet_gut_carveme.qza` from the MICOM media collection
(https://github.com/micom-dev/media, Apache-2.0), which maps the VMH/AGORA western
gut diet into the BiGG/CarveMe namespace -- the namespace muODE's CarveMe models
are actually in.  The recipe is `recipes/carveme.ipynb` in that repo.

    Diener C, Gibbons SM, Resendis-Antonio O.  "MICOM: Metagenome-Scale Modeling
    To Infer Metabolic Interactions in the Gut Microbiota."  mSystems 5, e00606-19
    (2020).  https://doi.org/10.1128/mSystems.00606-19

Modelling decisions (deliberate -- change them here, not downstream)
-------------------------------------------------------------------
* **`flux` becomes `max_uptake`, not `concentration`.**  They are different
  quantities in different units (mmol/gDW/h vs mmol/L) and conflating them is the
  bug this script exists to fix.  `max_uptake` caps how fast a cell may import the
  nutrient; `concentration` says how much of it is in the vessel.
* **Concentrations are still needed**, because muODE is a *dynamic* simulation: it
  tracks depletion, and a metabolite with no pool cannot be consumed.  The colon is
  a flow system, so we set each pool from the dietary flux over a residence time
  (`RESIDENCE_H`) at a nominal community density (`BIOMASS_G_PER_L`), and replenish
  it with an influx matched to the same flux.  Both constants are stated here and
  are the two assumptions a reviewer should argue with.
* **Fermentation products start at zero** but must exist as metabolites, or the
  engine has nowhere to put what the community secretes.
* **The published table is reproduced faithfully, including its oddities.**  We do
  not quietly "improve" someone else's medium; the script instead *reports* every
  bound far off the median (see OUTLIER_FACTOR) so they cannot pass unnoticed:
    - `h2o_e` and `meoh_e` are 10 mmol/gDW/h (the source marks both `dilution: 1.0`,
      i.e. free).  Water is uncontroversial.  **Methanol is a free carbon source at
      100x the median bound** -- it is real in the colon (pectin demethylation), but
      any model that can oxidise it gets a large gift.  Watch it.
    - `o2_e` carries a 0.001 mmol/gDW/h trace in the source, and we DROP it: see
      SKIPPED below for why.  The derived medium is anaerobic, with no o2_e row.
"""

from __future__ import annotations

import io
import sys
import urllib.request
import zipfile
from pathlib import Path

MEDIUM_URL = (
    "https://github.com/micom-dev/media/raw/HEAD/media/western_diet_gut_carveme.qza"
)

OUT = Path(__file__).resolve().parents[2] / "muode" / "data" / "diets" / "western_gut.csv"

#: Colonic residence time (h).  Transit through the large intestine is 12-48 h;
#: 24 h is the round number in the middle and the value MICOM's own gut model uses.
RESIDENCE_H = 24.0

#: Microbial density in the colon (gDW/L).  ~10^11 cells/mL at ~10^-12 gDW/cell.
BIOMASS_G_PER_L = 20.0

#: Fermentation products.  These ALWAYS start at zero concentration with **no**
#: uptake cap, overriding the source table where it supplies them (it gives
#: acetate, formate and H2 a dietary flux of 0.1 mmol/gDW/h).  Two reasons, both
#: load-bearing:
#:
#: 1. **They must be predicted, not assumed.**  Sizing an acetate pool from that
#:    dietary flux would start the vessel at 48 mM acetate -- roughly the colonic
#:    steady-state concentration the community is supposed to *produce*.  Scoring
#:    SCFA output against a medium that already contains the SCFAs is the same trap
#:    as lactate in DM38 (see muode/benchmarks/clark2021.py): a model that predicts
#:    nothing happens would look accurate.
#: 2. **Cross-feeding must not be throttled by dietary availability.**  A butyrate
#:    producer eats acetate made by *another species*, not acetate from food.
#:    Capping its uptake at the dietary flux would cap the cross-feeding that is the
#:    whole point of a community simulation.  So products are left uncapped and
#:    their uptake is governed by the cell's kinetics.
#:
#: The dietary contribution of SCFAs is small compared with microbial production,
#: so dropping it costs little; assuming it would cost the experiment.
PRODUCTS = ("ac_e", "but_e", "ppa_e", "lac__L_e", "lac__D_e", "succ_e", "for_e",
            "etoh_e", "h2_e", "ch4_e", "co2_e")

#: Dropped from the source table, with the reason.  Nothing is dropped silently.
SKIPPED = {
    "o2_e": (
        "The colon is anaerobic and the rest of muODE is built on that premise.  "
        "The source supplies a 0.001 mmol/gDW/h trace (the mucosal oxygen gradient) "
        "-- 100x BELOW the median bound, far too little to support aerobic growth, "
        "so dropping it changes no energetics.  Keeping it would, however, force us "
        "to assert a dissolved-O2 pool of 0.48 mM, which is about twice air "
        "saturation (~0.25 mM): an unphysical initial condition."
    ),
}

#: Trace minerals and vitamins ADDED to the source table (bound = MEDIAN_BOUND).
#:
#: These are not in the published medium and are not a measurement -- they are an
#: explicit modelling addition, listed here so a reviewer can delete them and see
#: what changes.  The justification: the source is a *dietary intake* table, and it
#: itemises food, not the trace minerals and vitamins that are ubiquitous in the
#: colon (from bile, host secretion, sloughed epithelium, and microbial synthesis).
#: Their absence is an artefact of what the table was for, not a claim that the gut
#: lacks nickel.
#:
#: Empirically, iJO1366 needs exactly ONE of these to grow on the medium at all --
#: nickel, a urease/hydrogenase cofactor -- and then reaches 0.052/h.  The rest are
#: included because real gut anaerobes are frequently vitamin auxotrophs, and a
#: model failing for want of biotin would be a medium artefact misread as biology.
#:
#: Deliberately NOT added: glucose (absorbed in the small intestine -- adding it
#: would smuggle the small intestine into the colon) and ammonium (in a community
#: it arises from amino-acid deamination; supplying it free would let every model
#: bypass proteolysis).
SUPPLEMENT = {
    # trace minerals -- catalytic, not stoichiometric; no carbon, no nitrogen
    "ni2_e": "nickel: urease / [NiFe]-hydrogenase cofactor (iJO1366 CANNOT grow without it)",
    "sel_e": "selenium: selenocysteine, formate dehydrogenase",
    "slnt_e": "selenite",
    "tungs_e": "tungstate: tungsten-dependent formate dehydrogenases in anaerobes",
    "na1_e": "sodium: Na+-coupled transport, ubiquitous",
    # vitamins / cofactor precursors -- gut anaerobes are commonly auxotrophic
    "btn_e": "biotin (B7)",
    "thm_e": "thiamine (B1)",
    "ribflv_e": "riboflavin (B2)",
    "pnto__R_e": "pantothenate (B5)",
    "nac_e": "nicotinate (B3)",
    "pydxn_e": "pyridoxine (B6)",
    "fol_e": "folate (B9)",
    "adocbl_e": "adenosylcobalamin (B12) -- many gut anaerobes cannot make it",
}

#: The source table's own default bound, used for anything without a measured
#: dietary flux.  We reuse it for SUPPLEMENT rather than inventing a new number.
MEDIAN_BOUND = 0.1

#: Report any bound this many times off the median.  The published table is
#: reproduced as-is; this makes its outliers loud instead of silent.
OUTLIER_FACTOR = 10.0


def fetch_medium() -> "list[tuple[str, float]]":
    """Download the published artifact and return (metabolite_id, flux) pairs."""
    import csv

    print(f"downloading {MEDIUM_URL}")
    with urllib.request.urlopen(MEDIUM_URL, timeout=120) as fh:  # noqa: S310
        blob = fh.read()
    z = zipfile.ZipFile(io.BytesIO(blob))
    root = z.namelist()[0].split("/")[0]
    text = z.read(f"{root}/data/medium.csv").decode()

    # Several source rows can collapse onto one `global_id` (the table is keyed by
    # MICOM's per-compartment `reaction`).  Deduplicate, and refuse to guess if the
    # duplicates disagree: last-wins would silently pick one flux over another.
    fluxes: "dict[str, float]" = {}
    for row in csv.DictReader(io.StringIO(text)):
        # `global_id` is the standard exchange (EX_ac_e); `reaction` is MICOM's
        # gut-compartment variant (EX_ac_m), which no CarveMe model has.
        ex = (row.get("global_id") or "").strip()
        if not ex.startswith("EX_"):
            continue
        met, flux = ex[len("EX_"):], float(row["flux"])
        if met in fluxes and fluxes[met] != flux:
            sys.exit(
                f"{met}: source gives conflicting fluxes {fluxes[met]:g} and {flux:g} "
                "-- resolve it here rather than letting one win silently"
            )
        fluxes[met] = flux
    if not fluxes:
        sys.exit("no exchanges parsed from the artifact -- has its schema changed?")
    return sorted(fluxes.items())


def main() -> int:
    raw = fetch_medium()
    dropped = [(m, f) for m, f in raw if m in SKIPPED]
    supplied = [(m, f) for m, f in raw if m not in SKIPPED and m not in PRODUCTS]
    limits = dict(supplied)

    # A pool sized so that, at the nominal density, the dietary flux would take
    # RESIDENCE_H to consume it.  mmol/gDW/h * gDW/L * h = mmol/L.
    def pool(flux: float) -> float:
        return flux * BIOMASS_G_PER_L * RESIDENCE_H

    lines = [
        "# western_gut -- anaerobic western-diet colonic medium, BiGG/CarveMe namespace",
        "#",
        "# DERIVED FILE -- do not hand-edit.  Regenerate with:",
        "#   python examples/diets/derive_western_gut.py",
        "#",
        "# max_uptake (mmol/gDW/h) is the PUBLISHED dietary flux bound, taken from",
        "# western_diet_gut_carveme.qza in the MICOM media collection",
        "# (github.com/micom-dev/media), which maps the VMH/AGORA western gut diet",
        "# into the BiGG namespace.  Cite:",
        "#   Diener C, Gibbons SM, Resendis-Antonio O.  MICOM: Metagenome-Scale",
        "#   Modeling To Infer Metabolic Interactions in the Gut Microbiota.",
        "#   mSystems 5, e00606-19 (2020).  doi:10.1128/mSystems.00606-19",
        "#",
        "# These bounds are the point of this file.  The median is 0.1 mmol/gDW/h;",
        "# muODE's uniform default Vmax was 10, i.e. 100x too permissive on every one",
        "# of ~90 nutrients at once, which is how GEMs reached 6.5/h (a 6-minute",
        "# doubling).  See muode/qc.py::MAX_PLAUSIBLE_GROWTH.",
        "#",
        f"# concentration (mmol/L) = max_uptake * {BIOMASS_G_PER_L:g} gDW/L * "
        f"{RESIDENCE_H:g} h  -- the pool a flow system of this",
        "# density and residence time sustains.  influx replenishes it at the same",
        "# rate.  Those two constants are assumptions, not measurements: they set the",
        "# size of the pool, NOT how fast a cell may eat from it (max_uptake does).",
        "#",
        "# Anaerobic: o2_e is dropped (see SKIPPED in the derive script).",
        "#",
        "# Fermentation products start at 0 with a BLANK max_uptake (= uncapped), even",
        "# where the source supplies them (it gives acetate, formate and H2 a dietary",
        "# flux).  They must be PREDICTED, not assumed -- seeding 48 mM acetate would",
        "# make SCFA production unfalsifiable -- and cross-feeding on them must be",
        "# governed by the cell's kinetics, not by a dietary ceiling.",
        "#",
        "# Trace minerals and vitamins are ADDED (see SUPPLEMENT in the derive script):",
        "# the source is a dietary-intake table and itemises food, not the ubiquitous",
        "# trace nutrients.  These are a modelling addition, not a measurement.",
        "metabolite,concentration,influx,max_uptake",
    ]

    for met, flux in supplied:
        lines.append(f"{met},{pool(flux):.6g},{flux:.6g},{flux:.6g}")
    for met in sorted(SUPPLEMENT):
        if met in limits:                      # already in the source: leave it alone
            continue
        limits[met] = MEDIAN_BOUND
        lines.append(f"{met},{pool(MEDIAN_BOUND):.6g},{MEDIAN_BOUND:.6g},{MEDIAN_BOUND:.6g}")
    for met in PRODUCTS:
        lines.append(f"{met},0.0,0.0,")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")

    fluxes = sorted(limits.values())
    mid = fluxes[len(fluxes) // 2]
    added = [m for m in sorted(SUPPLEMENT) if m not in dict(supplied)]
    print(f"wrote {OUT}")
    print(f"  {len(supplied)} from the published table + {len(added)} supplemented "
          f"+ {len(PRODUCTS)} product sinks")
    if added:
        print(f"  ADDED (not in the source; a modelling decision): {added}")
    seeded = [m for m, _ in raw if m in PRODUCTS]
    if seeded:
        print(f"  products the source SUPPLIES, forced to 0 and uncapped: {seeded}")
        print("    (else the vessel starts full of the SCFAs we are trying to predict)")
    for met, flux in dropped:
        print(f"  DROPPED {met} (source flux {flux:g}): {SKIPPED[met].split('.')[0]}.")
    print(f"  max_uptake: min={fluxes[0]:g}  median={mid:g}  max={fluxes[-1]:g} mmol/gDW/h")
    print(f"  (muODE's old uniform default Vmax was 10 -- {10 / mid:.0f}x the median "
          "bound, on every nutrient at once)")

    high = {m: f for m, f in limits.items() if f > mid * OUTLIER_FACTOR}
    low = {m: f for m, f in limits.items() if f < mid / OUTLIER_FACTOR}
    if high:
        print(f"\n  bounds >{OUTLIER_FACTOR:g}x the median -- these nutrients are near-free:")
        for m, f in sorted(high.items(), key=lambda kv: -kv[1]):
            note = " <-- a CARBON SOURCE, effectively unlimited" if m != "h2o_e" else ""
            print(f"    {m:10} {f:g}{note}")
    if low:
        print(f"\n  bounds <1/{OUTLIER_FACTOR:g} of the median -- trace only:")
        for m, f in sorted(low.items(), key=lambda kv: kv[1]):
            print(f"    {m:10} {f:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
