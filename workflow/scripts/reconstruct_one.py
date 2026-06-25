#!/usr/bin/env python
"""Reconstruct a GEM for a single MAG (Phase 1, one genome).

Used by the workflow's ``reconstruct`` rule for engines driven from Python: the
``stub`` engine (local end-to-end testing), ``gapseq``, and the eukaryote engines
``carvefungi`` (fungi) and ``eukaryote_generic`` (other eukaryotes).  CarveMe is
run directly from the rule's shell instead, inside its own conda environment.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genome", help="input MAG nucleotide FASTA")
    p.add_argument("output", help="output SBML model")
    p.add_argument(
        "--engine", default="stub",
        choices=["stub", "gapseq", "carveme", "carvefungi", "eukaryote_generic"],
    )
    p.add_argument("--universe", default="bacteria")
    p.add_argument("--gapfill-media", default=None)
    p.add_argument("--ref-db", default=None,
                   help="MetaEuk protein reference DB (eukaryote engines)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--carvefungi-cmd", default=None,
                   help="override CarveFungi invocation; {proteins}/{output} placeholders")
    p.add_argument("--eggnog-data-dir", default=None,
                   help="eggNOG-mapper data directory (eukaryote_generic engine)")
    args = p.parse_args()

    from muode.reconstruct import reconstruct_mag

    out = reconstruct_mag(
        args.genome, args.output, engine=args.engine,
        universe=args.universe, gapfill_media=args.gapfill_media,
        ref_db=args.ref_db, threads=args.threads,
        carvefungi_cmd=args.carvefungi_cmd, eggnog_data_dir=args.eggnog_data_dir,
    )
    print(f"reconstructed [{args.engine}] {Path(args.genome).stem} -> {out}")


if __name__ == "__main__":
    main()
