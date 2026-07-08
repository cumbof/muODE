#!/usr/bin/env python3
"""Emit muODE's `domains.tsv` (mag_id, domain) + a catalogue quality table.

Domain is decided per catalogue MAG:
  * id contains `__euk.`           -> eukaryote
  * GTDB-Tk lineage is d__Archaea  -> archaea
  * otherwise (prokaryotic MAG)    -> bacteria

`domains.tsv` is drop-in for muODE's `traits:` option (../workflow, config key
`traits`), which routes each MAG to the right reconstruction engine
(bacteria/archaea -> CarveMe/gapseq; eukaryote -> MetaEuk route).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_gtdbtk(path: Path) -> dict:
    lineage = {}
    if path and path.is_file():
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                g = row.get("user_genome")
                c = row.get("classification", "")
                if g:
                    lineage[g] = c
    return lineage


def load_quality(path: Path) -> dict:
    q = {}
    if path and path.is_file():
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                mid = row.get("mag_id")
                if mid:
                    q[mid] = row
    return q


def domain_of(mag_id: str, lineage: dict) -> str:
    if "__euk." in mag_id:
        return "eukaryote"
    tax = lineage.get(mag_id, "")
    if tax.startswith("d__Archaea"):
        return "archaea"
    return "bacteria"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalogue", required=True, type=Path,
                    help="dir of dereplicated catalogue MAG FASTAs")
    ap.add_argument("--gtdbtk", type=Path, help="GTDB-Tk merged summary TSV")
    ap.add_argument("--quality", type=Path, help="merged all_mags.quality.tsv")
    ap.add_argument("--out-domains", required=True, type=Path)
    ap.add_argument("--out-quality", required=True, type=Path)
    ap.add_argument("--ext", default="fa")
    args = ap.parse_args(argv)

    lineage = load_gtdbtk(args.gtdbtk)
    quality = load_quality(args.quality)

    mags = sorted(p.stem for p in args.catalogue.glob(f"*.{args.ext}"))

    args.out_domains.parent.mkdir(parents=True, exist_ok=True)
    counts = {}
    with args.out_domains.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mag_id", "domain"])
        for mid in mags:
            dom = domain_of(mid, lineage)
            counts[dom] = counts.get(dom, 0) + 1
            w.writerow([mid, dom])

    with args.out_quality.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mag_id", "domain", "completeness", "contamination",
                    "classification", "source_sample"])
        for mid in mags:
            dom = domain_of(mid, lineage)
            q = quality.get(mid, {})
            w.writerow([mid, dom, q.get("completeness", ""),
                        q.get("contamination", ""),
                        lineage.get(mid, ""), q.get("sample", "")])

    summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"catalogue: {len(mags)} MAG(s) [{summary or 'empty'}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
