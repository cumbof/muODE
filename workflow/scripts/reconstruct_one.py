#!/usr/bin/env python
"""Reconstruct a GEM for a single MAG (Phase 1, one genome).

Used by the workflow's ``reconstruct`` rule for engines driven from Python (the
``stub`` engine for local end-to-end testing, and ``gapseq``).  CarveMe is run
directly from the rule's shell instead, inside its own conda environment.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genome", help="input MAG nucleotide FASTA")
    p.add_argument("output", help="output SBML model")
    p.add_argument("--engine", default="stub", choices=["stub", "gapseq", "carveme"])
    p.add_argument("--universe", default="bacteria")
    p.add_argument("--gapfill-media", default=None)
    args = p.parse_args()

    from muode.reconstruct import reconstruct_mag

    out = reconstruct_mag(
        args.genome, args.output, engine=args.engine,
        universe=args.universe, gapfill_media=args.gapfill_media,
    )
    print(f"reconstructed [{args.engine}] {Path(args.genome).stem} -> {out}")


if __name__ == "__main__":
    main()
