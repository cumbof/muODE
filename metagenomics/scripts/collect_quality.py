#!/usr/bin/env python3
"""Merge per-sample MAG quality tables and emit a dRep genomeInfo CSV.

Concatenates the per-sample `*.quality.tsv` files written by select_mags.py into
one table, and (optionally) writes a dRep-compatible genomeInfo CSV
(`genome,completeness,contamination`) so dRep can skip its internal CheckM run
(essential for eukaryotic MAGs, which CheckM cannot score).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quality", nargs="*", type=Path, default=[],
                    help="per-sample quality TSVs (mag_id, completeness, ...)")
    ap.add_argument("--merged-out", required=True, type=Path)
    ap.add_argument("--genomeinfo-out", type=Path,
                    help="optional dRep genomeInfo CSV")
    ap.add_argument("--ext", default="fa")
    args = ap.parse_args(argv)

    rows = []
    header = ["mag_id", "completeness", "contamination", "source_bin", "sample"]
    for path in args.quality:
        if not path.is_file():
            continue
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if not row.get("mag_id"):
                    continue
                rows.append(row)

    args.merged_out.parent.mkdir(parents=True, exist_ok=True)
    with args.merged_out.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, "") for c in header])

    if args.genomeinfo_out:
        with args.genomeinfo_out.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["genome", "completeness", "contamination"])
            for r in rows:
                w.writerow([f"{r['mag_id']}.{args.ext}",
                            r.get("completeness", "50"),
                            r.get("contamination", "0")])

    print(f"collected {len(rows)} MAG(s) from {len(args.quality)} sample table(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
