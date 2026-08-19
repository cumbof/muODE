#!/usr/bin/env python3
"""Resolve and download the 26 Clark et al. strain genomes from NCBI.

    python examples/benchmarks/clark2021/fetch_genomes.py --outdir data/clark2021/genomes

Why this resolves rather than hardcodes accessions
--------------------------------------------------
Accessions typed from memory are exactly the kind of thing that is wrong and
looks right.  This queries the NCBI Datasets API, matches each assembly's
*infraspecific strain* against the strain Clark et al. actually grew, and writes a
manifest recording what it chose and how confident the match was.  Nothing is
silently substituted.

Why the strain matters (and this is the whole point)
----------------------------------------------------
Metabolic phenotype is a strain property.  NCBI's species reference for
*Faecalibacterium prausnitzii* is M21/2; Clark et al. measured **A2-165**.  Carve a
GEM from the reference and you have validated a model of an organism nobody grew.
So an exact strain match is required; a species-level fallback is *reported as
such* (`match=species-only`) and must be justified in the paper or replaced by
hand, not accepted because the script exited zero.

Output
------
* ``<outdir>/<CODE>.fna`` -- one genome per strain code (AC.fna, BT.fna, ...),
  named so the muODE workflow's `mags_dir` can consume them directly.
* ``<outdir>/manifest.tsv`` -- code, species, wanted strain, chosen accession,
  the strain NCBI reports, match quality, assembly level.  **Read this file.**
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from clark2021 import STRAINS

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


#: A substring match on a short token is a coin flip, not evidence.  NCBI reports
#: F. prausnitzii A2-165's strain as "1"; a naive `"1" in "a2165"` test "matches"
#: and hands you the wrong genome with a clean bill of health.  Below this length,
#: only exact equality counts.
MIN_FUZZY_LEN = 4


def _norm(text: str) -> str:
    """Strain designations are written a dozen ways: 'A2-165', 'A2_165', 'A2 165'."""
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def strains_match(want: str, got: str) -> bool:
    """Do two strain designations refer to the same strain?

    Exact after normalisation, or one contains the other *and* the contained token
    is long enough to be discriminating (:data:`MIN_FUZZY_LEN`).
    """
    if not want or not got:
        return False
    if want == got:
        return True
    shorter = min(want, got, key=len)
    if len(shorter) < MIN_FUZZY_LEN:
        return False
    return want in got or got in want


def candidates(species: str, limit: int = 50) -> list[dict]:
    """Every assembly NCBI holds for this species (annotated ones first)."""
    url = (f"{API}/genome/taxon/{urllib.parse.quote(species)}/dataset_report"
           f"?page_size={limit}&filters.assembly_version=current")
    with urllib.request.urlopen(url, timeout=60) as fh:  # noqa: S310
        return json.load(fh).get("reports", []) or []


def pick(species: str, strain: str) -> tuple[dict | None, str]:
    """Choose the assembly matching ``strain``; fall back to the species reference."""
    reports = candidates(species)
    if not reports:
        return None, "not-found"

    want = _norm(strain)
    for report in reports:
        info = report.get("assembly_info", {})
        got = _norm(info.get("biosample", {}).get("strain") or info.get("infraspecific_names", {}).get("strain") or "")
        if strains_match(want, got):
            return report, "strain"

    # The strain field is often blank or wrong; the organism name usually still
    # carries the designation ("Faecalibacterium prausnitzii A2-165").
    for report in reports:
        if strains_match(want, _norm(report.get("organism", {}).get("organism_name", ""))):
            return report, "strain-via-name"

    # No exact strain: prefer a reference/representative assembly, then anything.
    for report in reports:
        if report.get("assembly_info", {}).get("refseq_category") in ("reference genome", "representative genome"):
            return report, "species-only"
    return reports[0], "species-only"


def download(accession: str, dest: Path) -> None:
    """Fetch the genomic FASTA for one accession and write it uncompressed."""
    url = (f"{API}/genome/accession/{accession}/download"
           f"?include_annotation_type=GENOME_FASTA")
    with urllib.request.urlopen(url, timeout=300) as fh:  # noqa: S310
        blob = io.BytesIO(fh.read())

    with zipfile.ZipFile(blob) as zf:
        names = [n for n in zf.namelist() if n.endswith((".fna", ".fna.gz"))]
        if not names:
            raise RuntimeError(f"{accession}: no FASTA in the NCBI archive")
        data = zf.read(names[0])
    if names[0].endswith(".gz"):
        data = gzip.decompress(data)
    dest.write_bytes(data)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=Path("data/clark2021/genomes"))
    ap.add_argument("--dry-run", action="store_true", help="resolve and write the manifest, download nothing")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows, inexact, missing = [], [], []
    for code, s in sorted(STRAINS.items()):
        report, match = pick(s.taxon, s.strain)
        if report is None:
            missing.append(code)
            print(f"  {code}  {s.species:38} NOT FOUND", file=sys.stderr)
            continue

        info = report.get("assembly_info", {})
        accession = report["accession"]
        got_strain = (info.get("biosample", {}).get("strain")
                      or info.get("infraspecific_names", {}).get("strain") or "-")
        level = info.get("assembly_level", "-")

        flag = " " if match.startswith("strain") else "!"
        print(f"{flag} {code}  {s.species:38} want={s.strain:12} got={got_strain:16} "
              f"{accession}  [{match}, {level}]")
        if not match.startswith("strain"):
            inexact.append(f"{code} ({s.species}: wanted {s.strain!r}, got {got_strain!r})")

        if not args.dry_run:
            download(accession, args.outdir / f"{code}.fna")

        rows.append([code, s.species, s.strain, accession, got_strain, match, level])

    manifest = args.outdir / "manifest.tsv"
    manifest.write_text(
        "code\tspecies\twanted_strain\taccession\tncbi_strain\tmatch\tassembly_level\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in rows)
    )
    print(f"\nwrote {manifest} ({len(rows)}/{len(STRAINS)} strains)")

    if inexact:
        print(f"\nWARNING: {len(inexact)} strain(s) resolved to the species reference, "
              f"NOT the strain Clark et al. grew:", file=sys.stderr)
        for line in inexact:
            print(f"  ! {line}", file=sys.stderr)
        print("\nMetabolic phenotype is strain-specific.  Either source the correct\n"
              "assembly by hand and overwrite the FASTA, or state this substitution\n"
              "explicitly when reporting the benchmark.  Do not let it pass silently.",
              file=sys.stderr)
    if missing:
        print(f"\nERROR: no assembly found for {missing}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
