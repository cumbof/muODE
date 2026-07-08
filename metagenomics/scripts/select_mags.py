#!/usr/bin/env python3
"""Filter a CheckM2 (or EukCC) quality report and copy passing bins under stable ids.

Reads a quality report TSV, keeps bins with completeness >= --min-completeness and
contamination <= --max-contamination, copies each passing bin FASTA to
    <outdir>/<sample>__<tag>.<k>.fa
(a stable, sample-provenanced, collision-free mag_id), and writes a tidy quality
table keyed by the new mag_id to --quality-out.

Handles both CheckM2 (`Name`, `Completeness`, `Contamination`) and EukCC
(`bin`, `completeness`, `contamination`) column spellings.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

# accepted column spellings -> canonical
_NAME = ("Name", "bin", "Bin", "genome")
_COMP = ("Completeness", "completeness", "Completeness_General", "comp")
_CONT = ("Contamination", "contamination", "cont")


def _pick(row: dict, names) -> str | None:
    for n in names:
        if n in row and row[n] != "":
            return row[n]
    return None


def _resolve_bin(bins: Path, name: str, ext: str) -> Path | None:
    for cand in (bins / name, bins / f"{name}.{ext}", bins / f"{name}.fa",
                 bins / f"{name}.fasta", bins / f"{name}.fna"):
        if cand.is_file():
            return cand
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--bins", required=True, type=Path, help="dir of bin FASTAs")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--quality-out", required=True, type=Path)
    ap.add_argument("--tag", default="bin", help="id infix (bin|euk)")
    ap.add_argument("--ext", default="fa")
    ap.add_argument("--min-completeness", type=float, default=50.0)
    ap.add_argument("--max-contamination", type=float, default=10.0)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.quality_out.parent.mkdir(parents=True, exist_ok=True)

    passing = []
    if args.report.is_file():
        with args.report.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            k = 0
            for row in reader:
                name = _pick(row, _NAME)
                comp = _pick(row, _COMP)
                cont = _pick(row, _CONT)
                if name is None or comp is None or cont is None:
                    continue
                if float(comp) < args.min_completeness:
                    continue
                if float(cont) > args.max_contamination:
                    continue
                src = _resolve_bin(args.bins, name, args.ext)
                if src is None:
                    print(f"warn: {args.sample}: passing bin {name!r} has no FASTA "
                          f"in {args.bins}", file=sys.stderr)
                    continue
                k += 1
                mag_id = f"{args.sample}__{args.tag}.{k}"
                shutil.copy(src, args.outdir / f"{mag_id}.fa")
                passing.append((mag_id, float(comp), float(cont), name))

    with args.quality_out.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["mag_id", "completeness", "contamination", "source_bin", "sample"])
        for mag_id, comp, cont, name in passing:
            w.writerow([mag_id, f"{comp:.2f}", f"{cont:.2f}", name, args.sample])

    print(f"{args.sample}: {len(passing)} MAG(s) passed "
          f"(comp>={args.min_completeness}, cont<={args.max_contamination})",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
