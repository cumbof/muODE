#!/usr/bin/env python
"""Gap-fill, QC and (optionally) kinetically refine a single draft GEM (Phase 2/3).

Invoked per-MAG by the workflow so reconstruction + refinement parallelise across
the cluster.  Writes the refined model, a small JSON QC report and -- when
``--predict-kinetics`` is given -- a ``{stem}.kinetics.json`` with predicted Km
and kcat values.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cobra

from muode.diet import load_diet
from muode.gapfill import ensure_biomass
from muode.qc import sanity_check_model


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("draft", help="input draft SBML model")
    p.add_argument("output", help="output refined SBML model")
    p.add_argument("--universal", help="universal model (SBML) for LP gap-filling")
    p.add_argument("--diet", help="diet preset name or CSV: judge growth on the "
                                  "simulation medium instead of the model's open one")
    p.add_argument("--qc-json", help="path to write the QC report")
    p.add_argument("--kinetics-json", help="path to write predicted kinetics")
    p.add_argument("--predict-kinetics", action="store_true",
                   help="predict Km + kcat with the heuristic predictor")
    p.add_argument("--predictor", default="heuristic", help="heuristic | dlkcat | km-ml")
    args = p.parse_args()

    stem = Path(args.draft).stem
    model = cobra.io.read_sbml_model(args.draft)
    universal = cobra.io.read_sbml_model(args.universal) if args.universal else None
    diet = load_diet(args.diet) if args.diet else None

    gap = ensure_biomass(model, universal=universal, diet=diet)
    qc = sanity_check_model(model)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, args.output)

    if args.kinetics_json:
        from muode.kinetics import KineticParameters

        kin = KineticParameters()
        if args.predict_kinetics:
            from muode.predict import get_predictor, refine_kinetics

            kin = refine_kinetics(model, stem, predictor=get_predictor(args.predictor))
        kin.to_json(args.kinetics_json)

    report = {"model": stem, **gap, "qc": qc}
    if args.qc_json:
        Path(args.qc_json).write_text(json.dumps(report, indent=2, default=str))
    keys = [k for k in ("model", "grows_now", "growth_on_diet") if k in report]
    print(json.dumps({k: report[k] for k in keys}, default=str))


if __name__ == "__main__":
    main()
