#!/usr/bin/env python
"""Assemble a community manifest from the per-MAG refined models (Phase 4a).

Takes the explicit list of models that passed QC plus an abundance profile and a
diet, and writes the ``community.json`` consumed by ``muode simulate``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from muode.io_utils import read_abundance


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models", nargs="+", required=True, help="refined SBML model paths")
    p.add_argument("--abundance", help="MAG abundance TSV (MetaSBT profile)")
    p.add_argument("--diet", default="western_gut", help="diet preset name or CSV path")
    p.add_argument("--total-biomass", type=float, default=0.01)
    p.add_argument("--out", required=True, help="output community.json")
    args = p.parse_args()

    abund = read_abundance(args.abundance) if args.abundance else {
        Path(m).stem: 1.0 for m in args.models
    }
    manifest = {
        "models": [str(Path(m).resolve()) for m in args.models],
        "abundances": {Path(m).stem: abund.get(Path(m).stem, 0.0) for m in args.models},
        "diet": args.diet,
        "total_biomass": args.total_biomass,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(manifest, indent=2))
    print(f"assembled community with {len(args.models)} model(s) -> {args.out}")


if __name__ == "__main__":
    main()
