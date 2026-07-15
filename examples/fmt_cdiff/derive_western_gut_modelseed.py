#!/usr/bin/env python3
"""Translate the BiGG `western_gut` diet into the ModelSEED namespace the gapseq
FMT community actually eats -- ONCE, from the models' own annotations.

Why this script exists
----------------------
muODE is BiGG-native: `western_gut.csv`, the kinetics and the marker genes all key
on BiGG ids like `glc__D_e`.  The FMT community, by contrast, is reconstructed with
gapseq and is ModelSEED-native (`cpd00027_e0`).  A ModelSEED community is
internally consistent -- its members cross-feed with NO id translation -- so the
*only* BiGG<->ModelSEED boundary in the whole simulation is the diet.  Rather than
harmonise every model's ~160 exchanges onto BiGG, we translate the one diet across.

How the map is built (and why it is trustworthy)
------------------------------------------------
`muode.media.build_bigg_to_modelseed` reads every gapseq model in `gems/` and, for
each exchange metabolite that carries a `bigg.metabolite` annotation, records
`bigg -> cpd`.  That auto-resolves the bulk of the diet from the models' OWN
annotations -- not from memory.  A small hand-verified layer
(`muode.media._CURATED_BIGG_TO_MODELSEED`) sits on top for the rows the annotations
cannot get right: the secreted SCFAs (no annotation), ammonia (annotated `nh3`, the
diet says `nh4`), lactose (two compounds share the `lcts` annotation), and a couple
of un-annotated trace nutrients read out of the models BY NAME.  Anything neither
layer resolves is DROPPED and reported below -- never guessed.

    python examples/fmt_cdiff/derive_western_gut_modelseed.py

Output
------
`gems/western_gut_modelseed.csv` -- the same diet, same concentrations / influx /
max_uptake / provenance, with ModelSEED ids.  Plus a coverage report to stderr:
what resolved, what did not, and (per model) whether any UNRESOLVED row is one the
model needs to grow -- so the curation of the tail is driven by that, not by taste.

Note on completeness
--------------------
The unresolved tail shrinks as more members are reconstructed: a compound only
B. thetaiotaomicron carries (a polysaccharide, say) cannot resolve until B. theta's
model is present.  Re-run this once all five gapseq models are in `gems/`.
"""
from __future__ import annotations

import glob
import gzip
import logging
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEMS = HERE / "gems"
OUT = GEMS / "western_gut_modelseed.csv"


def _load_models():
    import cobra

    logging.disable(logging.CRITICAL)  # cobra's SBML reader is very chatty
    paths = sorted(glob.glob(str(GEMS / "*.gapseq.xml.gz")))
    # the first two members predate the `.gapseq.` naming convention
    for legacy in ("F_prausnitzii_A2165.xml.gz", "R_intestinalis_L182.xml.gz"):
        p = GEMS / legacy
        if p.exists() and str(p) not in paths:
            paths.append(str(p))
    models = {}
    for p in sorted(paths):
        with gzip.open(p, "rt") as fh:
            models[Path(p).name] = cobra.io.read_sbml_model(fh)
    return models


def main() -> int:
    from muode.diet import _western_gut
    from muode.media import (
        build_bigg_to_modelseed,
        diet_medium,
        translate_diet_to_modelseed,
    )

    models = _load_models()
    if not models:
        print(f"no gapseq models in {GEMS}", file=sys.stderr)
        return 1
    print(f"building map from {len(models)} model(s): {', '.join(models)}",
          file=sys.stderr)

    mapping, report = build_bigg_to_modelseed(models.values())
    diet = _western_gut()
    seed, unresolved = translate_diet_to_modelseed(diet, mapping)

    # --- write the CSV, provenance carried over, with a documented banner ------
    banner = [
        "# western_gut (ModelSEED / gapseq namespace) -- TRANSLATION of western_gut.csv",
        "#",
        "# DERIVED FILE -- do not hand-edit.  Regenerate with:",
        "#   python examples/fmt_cdiff/derive_western_gut_modelseed.py",
        "#",
        "# This is the BiGG western_gut diet with its ids mapped to ModelSEED, so the",
        "# gapseq FMT community can be fed by it.  Concentrations, influx, max_uptake and",
        "# the `source` provenance column are carried over UNCHANGED from western_gut.csv;",
        "# only the metabolite id changes (glc__D_e -> cpd00027_e0).  The mapping is built",
        "# from the gapseq models' own bigg.metabolite annotations plus a small verified",
        "# pin layer (muode.media._CURATED_BIGG_TO_MODELSEED); see the derive script.",
        "#",
        f"# Built from {len(models)} gapseq model(s): {', '.join(sorted(models))}.",
        f"# {len(seed.metabolites())} of {len(diet.metabolites())} diet rows resolved; "
        f"{len(unresolved)} unresolved (dropped, listed below).",
        "# UNRESOLVED (no ModelSEED id in any present model; dropped from this medium):",
    ]
    # wrap the unresolved list at ~4 per comment line so it stays readable
    for i in range(0, len(unresolved), 6):
        banner.append("#   " + ", ".join(unresolved[i:i + 6]))
    banner.append("metabolite,concentration,influx,max_uptake,source")

    lines = list(banner)
    for met in sorted(seed.metabolites()):
        limit = seed.uptake_limit(met)
        lines.append(
            f"{met},{seed.initial_concentration(met):.6g},"
            f"{seed.influx_rate(met):.6g},"
            f"{'' if limit is None else format(limit, '.6g')},"
            f"{seed.provenance(met) or ''}"
        )
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT}  ({len(seed.metabolites())} rows)", file=sys.stderr)

    # --- report: coverage, ambiguity, and essentiality of the unresolved tail --
    print("\n===== COVERAGE =====", file=sys.stderr)
    print(f"resolved   : {len(seed.metabolites())}/{len(diet.metabolites())}",
          file=sys.stderr)
    print(f"ambiguous  : {report['ambiguous'] or 'none'}", file=sys.stderr)
    print(f"unresolved : {len(unresolved)}", file=sys.stderr)
    print("  " + ", ".join(unresolved), file=sys.stderr)

    print("\n===== does any UNRESOLVED row block growth? =====", file=sys.stderr)
    print("(monoculture growth on this lean colonic medium is expected to be low --",
          file=sys.stderr)
    print(" these are fibre cross-feeders; the community, not the diet, feeds them)",
          file=sys.stderr)
    for name, m in models.items():
        default_mu = m.slim_optimize() or 0.0
        with m:
            m.medium = diet_medium(m, seed)
            tr_mu = m.slim_optimize() or 0.0
        print(f"  {name:34s} own medium {default_mu:6.4f}/h   "
              f"translated western_gut {tr_mu:6.4f}/h", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
