#!/usr/bin/env python3
"""Assemble ONE sample's final genome set: its MAGs plus any fetched references.

This is the set muODE simulates, and the set CoverM quantifies. Building it before
quantification is what makes the abundance table coherent: one tool measures one
genome set against one denominator, so there is never a second table to merge, and
never a CoverM percentage standing next to a MetaPhlAn percentage (a read fraction
next to a marker-normalised cell fraction).

Adding the reference genomes also *improves* the MAG abundances. Without them,
reads from a species that failed to assemble have nowhere correct to go and map
onto its assembled relatives, inflating them.

Emits, alongside the genome directory:
  * genome_sources.tsv   -- (genome_id, source, accession, gtdb_species) provenance;
                            `source` is `mag` or `reference`, so a proxy genome is
                            never mistaken for one assembled from the sample;
  * classification.tsv   -- GTDB-Tk-format lineages for MAGs *and* references, so
                            downstream domain routing treats them identically.
"""
from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def load_classification(path: Path) -> dict:
    out = {}
    if path and path.is_file():
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                genome = (row.get("user_genome") or "").strip()
                if genome:
                    out[genome] = (row.get("classification") or "").strip()
    return out


def species_of(lineage: str) -> str:
    species = ""
    for token in lineage.split(";"):
        token = token.strip()
        if token.startswith("s__"):
            species = token[3:].strip()
    return species


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--catalogue", required=True, type=Path,
                    help="this sample's dereplicated MAG FASTAs")
    ap.add_argument("--gtdbtk", type=Path, help="GTDB-Tk summary for the MAGs")
    ap.add_argument("--references", type=Path, help="dir of fetched reference FASTAs")
    ap.add_argument("--fetched", type=Path, help="fetched.tsv from fetch_reference_genomes.py")
    ap.add_argument("--ext", default="fa")
    ap.add_argument("--out-mags", required=True, type=Path)
    ap.add_argument("--out-sources", required=True, type=Path)
    ap.add_argument("--out-classification", required=True, type=Path)
    args = ap.parse_args(argv)

    out_mags = args.out_mags
    if out_mags.exists():
        shutil.rmtree(out_mags)
    out_mags.mkdir(parents=True)

    mag_lineage = load_classification(args.gtdbtk) if args.gtdbtk else {}
    sources, classification = [], {}

    for fasta in sorted(args.catalogue.glob(f"*.{args.ext}")):
        genome_id = fasta.stem
        shutil.copyfile(fasta, out_mags / fasta.name)
        lineage = mag_lineage.get(genome_id, "")
        classification[genome_id] = lineage
        sources.append({"genome_id": genome_id, "source": "mag", "accession": "",
                        "gtdb_species": species_of(lineage)})

    if args.fetched and args.fetched.is_file():
        ok = {"fetched", "cached"}
        with args.fetched.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("status") not in ok:
                    continue
                genome_id = (row.get("genome_id") or "").strip()
                if not genome_id:
                    continue
                src = (args.references / f"{genome_id}.{args.ext}") if args.references else None
                if not (src and src.is_file()):
                    raise SystemExit(
                        f"build_genome_set: {args.sample}: {row['status']} genome "
                        f"'{genome_id}' has no FASTA at {src}")
                shutil.copyfile(src, out_mags / src.name)
                lineage = (row.get("gtdb_lineage") or "").strip()
                classification[genome_id] = lineage
                sources.append({"genome_id": genome_id, "source": "reference",
                                "accession": (row.get("accession") or "").strip(),
                                "gtdb_species": species_of(lineage)})

    args.out_sources.parent.mkdir(parents=True, exist_ok=True)
    with args.out_sources.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "source", "accession", "gtdb_species"])
        for r in sources:
            w.writerow([r["genome_id"], r["source"], r["accession"], r["gtdb_species"]])

    args.out_classification.parent.mkdir(parents=True, exist_ok=True)
    with args.out_classification.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["user_genome", "classification"])
        for genome_id in sorted(classification):
            w.writerow([genome_id, classification[genome_id]])

    n_mag = sum(1 for r in sources if r["source"] == "mag")
    n_ref = sum(1 for r in sources if r["source"] == "reference")
    print(f"[{args.sample}] genome set: {n_mag} MAG(s) + {n_ref} reference(s) "
          f"= {len(sources)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
