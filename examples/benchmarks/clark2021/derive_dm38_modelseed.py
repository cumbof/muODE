#!/usr/bin/env python3
"""Translate the BiGG DM38 medium into the ModelSEED namespace a gapseq clark2021
community eats -- ONCE, from the models' own annotations.

Why this script exists
----------------------
`dm38.csv` is BiGG (`glc__D_e`), because the Clark benchmark reconstructs with CarveMe
by default.  But CarveMe's BiGG universe has NO connected butyrate pathway for gut
anaerobes, so its Tier-1 butyrate producers (AC, CC, ER, RI) reconstruct WITHOUT a
butyrate exchange and score as non-producers -- a reconstruction failure the benchmark
correctly detects but that understates what muODE can do.  gapseq builds bottom-up from
MetaCyc and secretes butyrate natively; to feed a gapseq community you need DM38 in the
ModelSEED namespace, because a gapseq community is ModelSEED-native and internally
consistent (its members cross-feed with no id translation).  The only BiGG<->ModelSEED
boundary is the medium, translated here once -- exactly as
`examples/gut_western/derive_medium_modelseed.py` does for western_gut.

How the map is built (and why it is trustworthy)
------------------------------------------------
`muode.media.build_bigg_to_modelseed` reads every gapseq model in ``--gems`` and, for
each exchange metabolite carrying a `bigg.metabolite` annotation, records `bigg -> cpd`.
That resolves the bulk of DM38 from the models' OWN annotations -- not from memory.  A
small hand-verified pin layer (`muode.media._CURATED_BIGG_TO_MODELSEED`) sits on top for
the rows annotations cannot get right (secreted SCFAs, ammonia nh3-vs-nh4, lactose).
Anything neither layer resolves is DROPPED and reported -- never guessed.

DM38 is simpler than western_gut in two ways that matter here:
  * it is a SEALED BATCH medium: influx is 0 and there is no max_uptake (the cells run
    to substrate exhaustion), so only concentrations translate -- exactly as in dm38.csv.
  * it is chemically defined and carries no fibre polymers (its largest sugars are
    glucose, arabinose and maltose), so the starch degree-of-polymerisation rescale that
    western_gut needs is a no-op here.  The shared translator handles both regardless.

    python examples/benchmarks/clark2021/derive_dm38_modelseed.py \
        --gems results/clark2021/models_gapseq

Output
------
`muode/data/diets/dm38_modelseed.csv` by default (parallel to dm38.csv, so it can be
registered as a preset) -- the same medium, same concentrations / provenance, with
ModelSEED ids.  Plus a coverage report to stderr: which of the 26 clark strains have a
GEM, what resolved, what did not, and per-model growth on the translated medium.

Note on completeness
--------------------
The unresolved tail shrinks as more strains are reconstructed: a compound only one
species carries cannot resolve until that model is present.  RE-RUN once all 26 gapseq
GEMs are in ``--gems`` (the reconstruction is a workstation job), or the medium will be
missing substrates the not-yet-built members would have resolved.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import logging
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parents[2] / "muode" / "data" / "diets" / "dm38_modelseed.csv"


def _load_models(gems_dir: Path):
    import cobra

    logging.disable(logging.CRITICAL)  # cobra's SBML reader is very chatty
    paths = sorted(glob.glob(str(gems_dir / "*.xml.gz"))) + \
        sorted(glob.glob(str(gems_dir / "*.xml")))
    models = {}
    for p in paths:
        opener = gzip.open if p.endswith(".gz") else open
        with opener(p, "rt") as fh:
            models[Path(p).name] = cobra.io.read_sbml_model(fh)
    return models


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gems", type=Path, required=True,
                    help="directory of gapseq GEMs for the clark strains (<CODE>.xml[.gz])")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help=f"output CSV (default: {DEFAULT_OUT})")
    args = ap.parse_args()

    import clark2021 as ck
    from muode.diet import _dm38
    from muode.media import (
        build_bigg_to_modelseed,
        diet_medium,
        translate_diet_to_modelseed,
    )

    models = _load_models(args.gems)
    if not models:
        print(f"no gapseq models in {args.gems}", file=sys.stderr)
        return 1

    # Which of the 26 clark strains does this cover?  GEMs are named by strain code.
    covered = {Path(n).stem.split(".")[0] for n in models}
    known = covered & set(ck.STRAINS)
    missing = sorted(set(ck.STRAINS) - known)
    if missing:
        print(f"WARNING: {len(known)}/{len(ck.STRAINS)} clark strains have a GEM; missing "
              f"{missing}.  The medium will be missing substrates those strains would "
              f"resolve -- re-run once all 26 are reconstructed.", file=sys.stderr)

    mapping, report = build_bigg_to_modelseed(models.values())
    diet = _dm38()
    seed, unresolved = translate_diet_to_modelseed(diet, mapping)

    # --- write the CSV, provenance carried over, with a documented banner ------
    banner = [
        "# DM38 (ModelSEED / gapseq namespace) -- TRANSLATION of dm38.csv",
        "#",
        "# DERIVED FILE -- do not hand-edit.  Regenerate with:",
        "#   python examples/benchmarks/clark2021/derive_dm38_modelseed.py --gems <dir>",
        "#",
        "# The DM38 medium of Clark et al. 2021 (Nat Commun 12, 3254), with its BiGG ids",
        "# mapped to ModelSEED so the gapseq-reconstructed clark community can be fed by it.",
        "# Concentrations and the `source` provenance column are carried over UNCHANGED from",
        "# dm38.csv; only the metabolite id changes (glc__D_e -> cpd00027_e0).  Like dm38.csv",
        "# this is a SEALED BATCH medium: influx is 0 and there is no max_uptake.  The map is",
        "# built from the gapseq models' own bigg.metabolite annotations plus the verified pin",
        "# layer (muode.media._CURATED_BIGG_TO_MODELSEED); see the derive script.",
        "#",
        f"# Built from {len(models)} gapseq model(s) covering {len(known)}/{len(ck.STRAINS)} "
        f"clark strains: {', '.join(sorted(known))}.",
        f"# {len(seed.metabolites())} of {len(diet.metabolites())} medium rows resolved; "
        f"{len(unresolved)} unresolved (dropped, listed below).",
        "# UNRESOLVED (no ModelSEED id in any present model; dropped from this medium):",
    ]
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
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}  ({len(seed.metabolites())} rows)", file=sys.stderr)

    # --- report: coverage, ambiguity, and essentiality of the unresolved tail --
    print("\n===== COVERAGE =====", file=sys.stderr)
    print(f"resolved   : {len(seed.metabolites())}/{len(diet.metabolites())}",
          file=sys.stderr)
    print(f"ambiguous  : {report['ambiguous'] or 'none'}", file=sys.stderr)
    print(f"unresolved : {len(unresolved)}", file=sys.stderr)
    print("  " + ", ".join(unresolved), file=sys.stderr)

    print("\n===== growth on the translated DM38 (sanity) =====", file=sys.stderr)
    print("(sealed batch: endpoint is set by yield, not rate; low mu is expected)",
          file=sys.stderr)
    for name, m in models.items():
        with m:
            m.medium = diet_medium(m, seed)
            tr_mu = m.slim_optimize() or 0.0
        print(f"  {name:34s} DM38 (ModelSEED) {tr_mu:6.4f}/h", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
