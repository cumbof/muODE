#!/usr/bin/env python
"""Gap-fill and QC a single draft GEM (Phase 2/3, one MAG).

Invoked per-MAG by the workflow so reconstruction + refinement parallelise across
the cluster.  Writes the refined model and a small JSON QC report.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cobra

from muode.gapfill import ensure_biomass
from muode.qc import sanity_check_model


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("draft", help="input draft SBML model")
    p.add_argument("output", help="output refined SBML model")
    p.add_argument("--universal", help="universal model (SBML) for LP gap-filling")
    p.add_argument("--qc-json", help="path to write the QC report")
    args = p.parse_args()

    model = cobra.io.read_sbml_model(args.draft)
    universal = cobra.io.read_sbml_model(args.universal) if args.universal else None

    gap = ensure_biomass(model, universal=universal)
    qc = sanity_check_model(model)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, args.output)

    report = {"model": Path(args.draft).stem, **gap, "qc": qc}
    if args.qc_json:
        Path(args.qc_json).write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps({k: report[k] for k in ("model", "grows_now")}, default=str))


if __name__ == "__main__":
    main()
