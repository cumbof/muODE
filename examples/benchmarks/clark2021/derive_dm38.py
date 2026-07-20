#!/usr/bin/env python3
"""Derive the DM38 defined medium (Clark et al. 2021) as a muODE diet CSV.

Why this is a script and not a hand-typed table
-----------------------------------------------
DM38 is the *published, chemically defined* medium the Clark et al. synthetic gut
communities were actually grown in, so it is the medium muODE must gap-fill and
simulate on to predict those measurements.  Hand-transcribing 67 components would
be unauditable; this script downloads Supplementary Data 4 from the publisher and
applies an explicit, reviewable name -> BiGG mapping, so anyone can re-derive
`muode/data/diets/dm38.csv` and diff it.

    python examples/benchmarks/clark2021/derive_dm38.py

Source
------
Clark RL, Connors BM, Stevenson DM, Hromada SE, Hamilton JJ, Amador-Noguez D,
Venturelli OS. "Design of synthetic human gut microbiome assembly and butyrate
production." Nat Commun 12, 3254 (2021).  https://doi.org/10.1038/s41467-021-22938-y
Supplementary Data 4 (DM38 composition), CC-BY 4.0.

Modelling decisions (all deliberate -- change them here, not downstream)
-----------------------------------------------------------------------
* **Salts are dissociated into ions.**  A GEM does not import "MgSO4"; it imports
  mg2_e and so4_e.  Where several salts supply the same ion the concentrations are
  summed (e.g. sulfate arrives from Na2SO4, MgSO4, MnSO4, ZnSO4, CuSO4, FeSO4).
* **Batch culture, so no influx.**  The experiment is a sealed 48 h incubation,
  not a chemostat: every influx rate is 0.  Nutrients deplete, which is precisely
  what makes the endpoint metabolite concentrations informative.
* **Anaerobic.**  No o2_e at any concentration.
* **Sodium lactate (28.3 mM) is a medium component**, not only a product.  Lactate
  is therefore both consumed and secreted, and any validation against measured
  lactate must score the *change* from this baseline, not the endpoint (see
  `clark2021.py`, which carries this baseline through).
* **Non-metabolites are dropped** with a reason recorded in `SKIPPED`: MOPS is a
  pH buffer, EDTA a chelator; neither is taken up.  Trace ions with no BiGG
  exchange in the standard universe (Al, B, W) are dropped -- they are catalytic,
  not stoichiometric, so no biomass reaction depends on them.
"""

from __future__ import annotations

import io
import sys
import urllib.request
from pathlib import Path

# Supplementary Data 4 (DM38 composition) -- Nature's static asset host.
SUPP_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-021-22938-y/MediaObjects/"
    "41467_2021_22938_MOESM5_ESM.xlsx"
)

OUT = Path(__file__).resolve().parents[3] / "muode" / "data" / "diets" / "dm38.csv"

# --- DM38 component -> BiGG extracellular metabolite(s) ---------------------
# A salt maps to every ion it releases.  Values are stoichiometric multipliers:
# CaCl2 at 1.29 mM gives 1.29 mM ca2_e and 2 x 1.29 mM cl_e.
MAPPING: dict[str, dict[str, float]] = {
    # -- bulk carbon sources
    "D-Glucose": {"glc__D_e": 1},
    "L-Arabinose": {"arab__L_e": 1},
    "Sodium Lactate": {"lac__L_e": 1, "na1_e": 1},
    "Maltose": {"malt_e": 1},
    "Inositol": {"inost_e": 1},
    "Sodium Bicarbonate": {"hco3_e": 1, "na1_e": 1},
    # -- amino acids (also major carbon/nitrogen sources for gut anaerobes)
    "L-alanine": {"ala__L_e": 1},
    "L-arginine": {"arg__L_e": 1},
    "L-asparagine": {"asn__L_e": 1},
    "L-Aspartic acid": {"asp__L_e": 1},
    "L-Cysteine": {"cys__L_e": 1},
    "L-Glutamic acid Sodium": {"glu__L_e": 1, "na1_e": 1},
    "L-glutamine": {"gln__L_e": 1},
    "L-Glycine": {"gly_e": 1},
    "L-histidine": {"his__L_e": 1},
    "L-isoleucine": {"ile__L_e": 1},
    "L-leucine": {"leu__L_e": 1},
    "L-lysine": {"lys__L_e": 1},
    "L-methionine": {"met__L_e": 1},
    "L-phenylalanine": {"phe__L_e": 1},
    "L-proline": {"pro__L_e": 1},
    "L-serine": {"ser__L_e": 1},
    "L-threonine": {"thr__L_e": 1},
    "L-tryptophan": {"trp__L_e": 1},
    "L-tyrosine": {"tyr__L_e": 1},
    "L-valine": {"val__L_e": 1},
    # -- nucleobases / nucleosides
    "L-Adenine": {"ade_e": 1},
    "Guanine": {"gua_e": 1},
    "Cytosine": {"csn_e": 1},
    "Uracil": {"ura_e": 1},
    "Xanthine": {"xan_e": 1},
    "Thymidine": {"thymd_e": 1},
    "Orotic Acid": {"orot_e": 1},
    # -- vitamins / cofactors
    "Biotin": {"btn_e": 1},
    "Ca-D-pantothenate": {"pnto__R_e": 1, "ca2_e": 0.5},
    "Cobalamin (B12)": {"cbl1_e": 1},
    "Folic Acid": {"fol_e": 1},
    "Tetrahydrofolic Acid": {"thf_e": 1},
    "Nicotinamide": {"ncam_e": 1},
    "Riboflavin": {"ribflv_e": 1},
    "Thiamine HCl": {"thm_e": 1, "cl_e": 1},
    "Pyrixodal": {"pydx_e": 1},           # sic: "Pyridoxal", as spelled in the source
    "Pyridoxine HCl": {"pydxn_e": 1, "cl_e": 1},
    "Pyridoxamine": {"pydam_e": 1},
    "p-Aminobenzoic Acid": {"4abz_e": 1},
    "Haemin": {"pheme_e": 1},
    # -- bulk salts: N, P, S, and the major ions
    "Ammonium Chloride": {"nh4_e": 1, "cl_e": 1},
    "Potassium Phosphate": {"pi_e": 1, "k_e": 1},
    "Sodium Sulfate": {"so4_e": 1, "na1_e": 2},
    "NaCl": {"na1_e": 1, "cl_e": 1},
    "Potassium Hydroxide": {"k_e": 1},
    "CaCl2": {"ca2_e": 1, "cl_e": 2},
    "MgSO4": {"mg2_e": 1, "so4_e": 1},
    "Magnesium Chloride": {"mg2_e": 1, "cl_e": 2},
    # -- trace metals
    "FeSO4": {"fe2_e": 1, "so4_e": 1},
    "MnSO4": {"mn2_e": 1, "so4_e": 1},
    "ZnSO4": {"zn2_e": 1, "so4_e": 1},
    "CuSO4": {"cu2_e": 1, "so4_e": 1},
    "Co(NO3)2": {"cobalt2_e": 1, "no3_e": 2},
    "NiCl2": {"ni2_e": 1, "cl_e": 2},
    "Na2MoO4": {"mobd_e": 1, "na1_e": 2},
    "Na2SeO3": {"slnt_e": 1, "na1_e": 2},
}

# Present in DM38 but deliberately not given to the models, with the reason.
SKIPPED: dict[str, str] = {
    "MOPS": "pH buffer -- not metabolised, and pH is handled by muode.ph",
    "EDTA": "chelator -- keeps trace metals soluble; not a nutrient",
    "AlK(SO4)2": "aluminium has no exchange in the BiGG universe (catalytic trace)",
    "H3BO3": "boron has no exchange in the BiGG universe (catalytic trace)",
    "Na2WO4": "tungstate has no exchange in most BiGG GEMs (catalytic trace)",
}

# Always available in an aqueous medium; GEMs expect them open.
IMPLICIT: dict[str, float] = {"h2o_e": 1000.0, "h_e": 1000.0}

# A DELIBERATE DEVIATION from the published medium.  DM38's only iron source is
# FeSO4, so the medium above carries Fe(II) alone -- chemically right for a reducing
# anaerobic culture.  But reconstructions do not agree on which oxidation state the
# cell imports: CarveMe and gapseq both build biomass that demands fe3_e/cpd10516,
# and neither ships a Fe(II)->Fe(III) conversion.  The result is that ALL 26 clark
# strains score exactly 0.0/h on DM38 in both namespaces -- a medium artefact that
# looks like 26 biological predictions of no growth.
#
# So we expose the SAME iron pool in both oxidation states.  This does not double the
# iron: 0.0658 mM is the measured total, and it is offered as either form because the
# GEMs pick a state arbitrarily.  Iron is ~1000x from limiting at this concentration
# (biomass needs ~1e-5 mmol/gDW), so which form is drawn on cannot shift a flux.
#
# Recorded here rather than patched into the CSV so the deviation stays visible: the
# published medium says Fe(II), muODE feeds Fe(II)+Fe(III), and this comment is why.
IRON_STATE_DEVIATION: dict[str, str] = {"fe3_e": "fe2_e"}


def load_dm38(source: str | Path | None = None) -> list[tuple[str, float]]:
    """Return [(component, mM)] from Supplementary Data 4 (URL or local xlsx)."""
    import openpyxl

    if source is None:
        with urllib.request.urlopen(SUPP_URL) as fh:  # noqa: S310  (fixed, https)
            blob = io.BytesIO(fh.read())
    else:
        blob = io.BytesIO(Path(source).read_bytes())

    ws = openpyxl.load_workbook(blob, read_only=True).worksheets[0]
    rows: list[tuple[str, float]] = []
    for name, conc, *_ in ws.iter_rows(min_row=3, max_col=3, values_only=True):
        if name is None or conc is None:
            continue
        rows.append((str(name).strip(), float(conc)))
    return rows


def to_diet(components: list[tuple[str, float]]) -> tuple[dict[str, float], list[str]]:
    """Map DM38 components onto BiGG exchange metabolites, summing shared ions."""
    conc: dict[str, float] = dict(IMPLICIT)
    unmapped: list[str] = []
    for name, mM in components:
        if name in SKIPPED:
            continue
        if name not in MAPPING:
            unmapped.append(name)
            continue
        for met, mult in MAPPING[name].items():
            conc[met] = conc.get(met, 0.0) + mult * mM
    # mirror the iron pool into the other oxidation state (see IRON_STATE_DEVIATION)
    for alias, src in IRON_STATE_DEVIATION.items():
        if src in conc:
            conc.setdefault(alias, conc[src])
    return conc, unmapped


def main() -> int:
    source = sys.argv[1] if len(sys.argv) > 1 else None
    components = load_dm38(source)
    conc, unmapped = to_diet(components)

    if unmapped:
        print(f"ERROR: {len(unmapped)} DM38 components are neither mapped nor "
              f"explicitly skipped: {unmapped}", file=sys.stderr)
        print("Add them to MAPPING or to SKIPPED (with a reason) -- silently "
              "dropping a nutrient would starve the models.", file=sys.stderr)
        return 1

    header = [
        "# DM38 -- the chemically defined medium of Clark et al. 2021, in the BiGG",
        "# namespace.  DERIVED FILE: regenerate with",
        "#   python examples/benchmarks/clark2021/derive_dm38.py",
        "#",
        "# Source: Clark et al., 'Design of synthetic human gut microbiome assembly",
        "# and butyrate production', Nat Commun 12, 3254 (2021), Supplementary Data 4.",
        "# https://doi.org/10.1038/s41467-021-22938-y  (CC-BY 4.0)",
        "#",
        "# Concentrations are mmol/L, exactly as published.  Influx is 0 throughout:",
        "# the experiment is a sealed 48 h batch culture, not a chemostat.",
        "# Anaerobic -- no o2_e.  Salts are dissociated into the ions a GEM imports,",
        "# summing shared ions across salts.",
        "#",
        "# NOTE: L-lactate is a MEDIUM COMPONENT here (28.3 mM sodium lactate), not",
        "# only a fermentation product.  Score measured lactate as a change from this",
        "# baseline, never as an endpoint concentration.",
        "#",
        "# NO max_uptake COLUMN, AND THAT IS DELIBERATE.  western_gut carries dietary",
        "# flux bounds because it is a flow system where the dietary flux is the datum.",
        "# DM38 is the opposite: a sealed batch culture whose CONCENTRATIONS are the",
        "# datum.  25 mM glucose is genuinely available, so a dietary ceiling here would",
        "# be a number we invented rather than one anybody measured.",
        "#",
        "# It would also break the benchmark.  Uptake is min(Michaelis-Menten, max_uptake),",
        "# so a ceiling clamps uptake BELOW Vmax, flattens growth-vs-Vmax, and makes Vmax",
        "# unidentifiable -- and DM38 is the medium the Clark benchmark exists to CALIBRATE",
        "# Vmax against.  See tests/test_dm38_bounds.py, which demonstrates the flattening",
        "# rather than asserting it.",
        "#",
        "# What DOES constrain uptake here is Vmax, and muODE's uniform default of 10 is",
        "# too permissive: on DM38 it grows iJO1366 at ~2.9/h, above muode.qc's ceiling of",
        "# 1.5.  The fix is a CALIBRATED Vmax (examples/benchmarks/clark2021/fit_kinetics.py),",
        "# not an invented dietary bound.",
        "#",
        "# IRON -- A DELIBERATE DEVIATION FROM THE PUBLISHED MEDIUM.  DM38's only iron",
        "# source is FeSO4, i.e. Fe(II), which is chemically right for a reducing",
        "# anaerobic culture.  But CarveMe and gapseq both build biomass demanding",
        "# fe3_e / cpd10516_e0 and ship no Fe(II)->Fe(III) conversion, so on the medium",
        "# AS PUBLISHED all 26 clark strains score exactly 0.0/h in both namespaces --",
        "# a medium artefact that reads like 26 biological predictions of no growth.",
        "# muODE therefore offers the SAME 0.0658 mM iron pool as BOTH fe2_e and fe3_e.",
        "# The iron is not doubled; it is one pool exposed in two oxidation states,",
        "# because the GEMs pick a state arbitrarily.  At ~1000x above the biomass",
        "# requirement iron cannot be limiting, so this cannot shift a predicted flux.",
        f"# {len(components)} components -> {len(conc)} exchange metabolites.",
        "metabolite,concentration,influx",
    ]
    lines = header + [f"{m},{c:.6g},0.0" for m, c in sorted(conc.items())]
    OUT.write_text("\n".join(lines) + "\n")

    print(f"wrote {OUT}")
    print(f"  {len(components)} DM38 components "
          f"({len(SKIPPED)} skipped) -> {len(conc)} BiGG metabolites")
    print(f"  carbon: glucose {conc['glc__D_e']:.1f}, arabinose {conc['arab__L_e']:.1f}, "
          f"lactate {conc['lac__L_e']:.1f}, maltose {conc['malt_e']:.1f} mM")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
