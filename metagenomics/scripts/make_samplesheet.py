#!/usr/bin/env python3
"""Build a metagenomics/ sample sheet from a directory of paired FASTQ files.

Scans a run directory (e.g. an SRA project download like PRJNA426573) for
paired reads named `<id>_1.fastq[.gz]` / `<id>_2.fastq[.gz]` (also `_R1`/`_R2`)
and prints a TSV with header `sample<TAB>r1<TAB>r2` to stdout.

Usage:
    python scripts/make_samplesheet.py /path/to/PRJNA426573 > config/samples.tsv
    python scripts/make_samplesheet.py /path/to/dir --glob '*_R1_001.fastq.gz'
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# forward-read stem -> (sample_id, mate2_stem) via these paired tokens
_PAIR_TOKENS = [("_1", "_2"), ("_R1", "_R2"), ("_R1_001", "_R2_001")]
_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _strip_ext(name: str) -> tuple[str, str]:
    for ext in _EXTS:
        if name.endswith(ext):
            return name[: -len(ext)], ext
    return name, ""


def find_pairs(directory: Path):
    files = {p.name for p in directory.iterdir() if p.is_file()}
    pairs, used = [], set()
    for name in sorted(files):
        if name in used:
            continue
        stem, ext = _strip_ext(name)
        if not ext:
            continue
        for t1, t2 in _PAIR_TOKENS:
            if stem.endswith(t1):
                sample = stem[: -len(t1)]
                mate = f"{sample}{t2}{ext}"
                if mate in files:
                    pairs.append((sample, name, mate))
                    used.update({name, mate})
                break
    return pairs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory", type=Path, help="directory containing paired FASTQ files")
    ap.add_argument("--abs", action="store_true", default=True,
                    help="emit absolute paths (default)")
    ap.add_argument("--relative", dest="abs", action="store_false",
                    help="emit paths as given (relative to the directory)")
    args = ap.parse_args(argv)

    d = args.directory
    if not d.is_dir():
        ap.error(f"not a directory: {d}")

    pairs = find_pairs(d)
    if not pairs:
        print(f"error: no *_1/_2 (or _R1/_R2) FASTQ pairs found in {d}", file=sys.stderr)
        return 1

    root = d.resolve() if args.abs else d
    print("sample\tr1\tr2")
    for sample, r1, r2 in pairs:
        print(f"{sample}\t{root / r1}\t{root / r2}")
    print(f"# {len(pairs)} paired samples found in {d}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
