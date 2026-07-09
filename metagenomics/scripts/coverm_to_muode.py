#!/usr/bin/env python3
"""Convert a CoverM genome relative-abundance table into muODE's abundance TSV.

CoverM `-m relative_abundance` emits `Genome<TAB><sample> Relative Abundance (%)`
with an `unmapped` row. muODE's abundance contract is a 2-column TSV
`(mag_id, rel_abundance)` over the modelled community, so we drop `unmapped` and
renormalise the remaining genomes to sum to 1 (community composition among the
catalogue MAGs). Pass --keep-percent to keep raw percentages instead.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coverm", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--keep-percent", action="store_true",
                    help="write raw CoverM percentages instead of renormalised fractions")
    args = ap.parse_args(argv)

    entries = []
    with args.coverm.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) < 2 or not row[0]:
                continue
            genome = row[0].strip()
            if genome.lower() == "unmapped":
                continue
            try:
                val = float(row[1])
            except ValueError:
                continue
            entries.append((genome, val))

    total = sum(v for _, v in entries)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["mag_id", "rel_abundance"])
        for genome, val in entries:
            if args.keep_percent:
                w.writerow([genome, f"{val:.6f}"])
            else:
                frac = (val / total) if total > 0 else 0.0
                w.writerow([genome, f"{frac:.6f}"])

    print(f"{args.out.name}: {len(entries)} MAG(s), total input mass {total:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
