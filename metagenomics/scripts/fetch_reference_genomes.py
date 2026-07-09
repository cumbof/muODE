#!/usr/bin/env python3
"""Fetch NCBI reference genomes for taxa that were PROFILED but never ASSEMBLED.

A quantitative profiler sees species that assembly cannot recover -- too
low-abundance, or too strain-diverse for binning to resolve. Those species are
real members of the community, but muODE cannot simulate them: a GEM is built
from a genome, and there is no genome. This script supplies one, standing in a
public reference genome of the same GTDB species.

Species -> accession is resolved from `gtdb_taxonomy.tsv`, which ships with the
GTDB-Tk reference data the pipeline already needs -- so no extra database, and the
accession is guaranteed to sit in the same taxonomy GTDB-Tk classified the MAGs
with. RefSeq (`RS_GCF_*`) is preferred over GenBank (`GB_GCA_*`); ties break
lexicographically, so the choice is deterministic across reruns.

Genomes are downloaded with the NCBI `datasets` CLI and written as
`<sample>__ref.<k>.fa`, alongside a GTDB-Tk-format classification table so the
downstream domain routing treats them exactly like MAGs.

HONEST SCOPE. The fetched genome is a *proxy*: a reference isolate of the species,
not the strain actually present in this sample. Its accessory-gene content, and
therefore its GEM, may differ from the real organism's. It restores the species'
metabolic presence in the simulation; it does not recover the sample's strain.
Every substitution is recorded in `reconciliation.tsv` with `source=reference` so
this never becomes invisible.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


def load_missing(path: Path, min_abundance: float, max_genomes: int):
    """Profiled-only species, richest first, filtered and capped."""
    rows = []
    if not path.is_file():
        return rows
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            species = (row.get("gtdb_species") or "").strip()
            lineage = (row.get("gtdb_lineage") or "").strip()
            try:
                pct = float(row.get("metaphlan_species_pct") or 0.0)
            except ValueError:
                continue
            if not species or not lineage:
                continue  # unmapped SGB: no GTDB species to resolve
            if pct < min_abundance:
                continue
            rows.append({"gtdb_species": species, "gtdb_lineage": lineage,
                         "metaphlan_species_pct": pct})
    rows.sort(key=lambda d: -d["metaphlan_species_pct"])
    return rows[:max_genomes] if max_genomes > 0 else rows


def resolve_accessions(gtdb_taxonomy: Path, wanted: set) -> dict:
    """{species: accession} from gtdb_taxonomy.tsv, preferring RefSeq.

    The file is ~600k rows (`RS_GCF_000005845.2<TAB>d__Bacteria;...;s__X`), so it
    is streamed and only the species we asked for are retained.
    """
    best: dict = {}
    with gtdb_taxonomy.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc, lineage = parts[0].strip(), parts[1].strip()
            species = ""
            for token in lineage.split(";"):
                token = token.strip()
                if token.startswith("s__"):
                    species = token[3:].strip()
            if species not in wanted:
                continue
            # strip the GTDB source prefix: RS_GCF_... / GB_GCA_... -> GCF_/GCA_
            clean = acc.split("_", 1)[1] if acc[:3] in ("RS_", "GB_") else acc
            rank = 0 if clean.startswith("GCF_") else 1  # RefSeq first
            prev = best.get(species)
            if prev is None or (rank, clean) < (prev[0], prev[1]):
                best[species] = (rank, clean)
    return {sp: acc for sp, (_rank, acc) in best.items()}


def download_genome(accession: str, dest: Path) -> bool:
    """Fetch one genome FASTA with the NCBI `datasets` CLI. True on success."""
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "genome.zip"
        cmd = ["datasets", "download", "genome", "accession", accession,
               "--include", "genome", "--filename", str(zip_path),
               "--no-progressbar"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=900)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
                FileNotFoundError) as exc:
            print(f"  download failed for {accession}: {exc}", file=sys.stderr)
            return False
        try:
            with zipfile.ZipFile(zip_path) as zf:
                members = [m for m in zf.namelist() if m.endswith((".fna", ".fasta"))]
                if not members:
                    print(f"  no FASTA inside {accession}.zip", file=sys.stderr)
                    return False
                zf.extract(members[0], Path(tmp) / "x")
                shutil.copyfile(Path(tmp) / "x" / members[0], dest)
        except (zipfile.BadZipFile, OSError) as exc:
            print(f"  could not unpack {accession}: {exc}", file=sys.stderr)
            return False
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--missing", required=True, type=Path,
                    help="missing_species.tsv from reconcile_profile.py")
    ap.add_argument("--gtdb-taxonomy", required=True, type=Path,
                    help="gtdb_taxonomy.tsv from the GTDB-Tk reference data")
    ap.add_argument("--outdir", required=True, type=Path, help="where FASTAs land")
    ap.add_argument("--out-table", required=True, type=Path)
    ap.add_argument("--out-classification", required=True, type=Path,
                    help="GTDB-Tk-format (user_genome, classification) for the fetched genomes")
    ap.add_argument("--min-abundance", type=float, default=0.1,
                    help="skip profiled species below this MetaPhlAn %% (default 0.1)")
    ap.add_argument("--max-genomes", type=int, default=50,
                    help="cap on genomes fetched per sample; 0 = no cap (default 50)")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve accessions but download nothing")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    wanted = load_missing(args.missing, args.min_abundance, args.max_genomes)

    accessions = {}
    if wanted:
        if not args.gtdb_taxonomy.is_file():
            print(f"fetch_reference_genomes: no GTDB taxonomy at {args.gtdb_taxonomy}",
                  file=sys.stderr)
            return 1
        accessions = resolve_accessions(args.gtdb_taxonomy,
                                        {d["gtdb_species"] for d in wanted})

    results, k = [], 0
    for entry in wanted:
        species = entry["gtdb_species"]
        accession = accessions.get(species, "")
        if not accession:
            results.append({**entry, "genome_id": "", "accession": "",
                            "status": "no_accession"})
            continue
        k += 1
        genome_id = f"{args.sample}__ref.{k}"
        dest = args.outdir / f"{genome_id}.fa"
        if args.dry_run:
            status = "planned"
        elif dest.is_file():
            status = "cached"
        elif download_genome(accession, dest):
            status = "fetched"
        else:
            k -= 1
            results.append({**entry, "genome_id": "", "accession": accession,
                            "status": "download_failed"})
            continue
        results.append({**entry, "genome_id": genome_id, "accession": accession,
                        "status": status})

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    with args.out_table.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "gtdb_species", "accession", "gtdb_lineage",
                    "metaphlan_species_pct", "status"])
        for r in results:
            w.writerow([r["genome_id"], r["gtdb_species"], r["accession"],
                        r["gtdb_lineage"], f"{r['metaphlan_species_pct']:.6f}", r["status"]])

    ok = {"fetched", "cached"}
    args.out_classification.parent.mkdir(parents=True, exist_ok=True)
    with args.out_classification.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["user_genome", "classification"])
        for r in results:
            if r["status"] in ok:
                w.writerow([r["genome_id"], r["gtdb_lineage"]])

    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) or "none"
    print(f"[{args.sample}] reference genomes: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
