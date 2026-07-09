#!/usr/bin/env python3
"""Final per-sample accounting: what was simulated, what was not, and why.

Joins three things for one sample:
  * status.tsv        -- matched / reconstructed_only / profiled_only (reconcile_profile.py)
  * fetched.tsv       -- which profiled_only species got a reference genome
  * abundance.tsv     -- CoverM over the final genome set (the single source of truth)

and answers the only question that matters when reading a simulation: *which
organisms is it actually about?*

Every genome and every profiled taxon appears exactly once, with `simulated` set
from the ground truth -- presence in the abundance table -- rather than inferred.
A `profiled_only` species whose reference genome could not be resolved is reported,
not quietly forgotten.

`pct_profiled_abundance_captured` is the headline number: of everything MetaPhlAn
saw, how much of it is represented by a genome in the simulation. It is routinely
well below 100%, and it belongs in any write-up of the results.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_tsv(path: Path) -> list:
    if not (path and path.is_file()):
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def coverm_unmapped_pct(path: Path) -> float:
    """The `unmapped` row: reads matching no genome in the final set."""
    if not (path and path.is_file()):
        return 0.0
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2 and row[0].strip().lower() == "unmapped":
                try:
                    return float(row[1])
                except ValueError:
                    return 0.0
    return 0.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--status", required=True, type=Path)
    ap.add_argument("--stats", required=True, type=Path)
    ap.add_argument("--abundance", required=True, type=Path)
    ap.add_argument("--coverm", type=Path)
    ap.add_argument("--sources", required=True, type=Path)
    ap.add_argument("--fetched", type=Path)
    ap.add_argument("--out-reconciliation", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    args = ap.parse_args(argv)

    status_rows = load_tsv(args.status)
    sources = {r["genome_id"]: r for r in load_tsv(args.sources)}
    fetched = load_tsv(args.fetched) if args.fetched else []
    abundance = {r["mag_id"]: float(r["rel_abundance"])
                 for r in load_tsv(args.abundance)}
    stats = (load_tsv(args.stats) or [{}])[0]

    # species -> the reference genome that now stands in for it
    ref_by_species = {r["gtdb_species"]: r for r in fetched}
    ok_status = {"fetched", "cached"}

    out_rows = []
    n_recovered = n_missing = 0
    captured_pct = 0.0
    # `metaphlan_species_pct` is a SPECIES-level percentage: when two MAGs share a
    # species they both carry it, so accumulating per genome would double-count.
    # Credit each species' profiled mass at most once.
    counted_species = set()

    def credit(species: str, pct: float) -> None:
        nonlocal captured_pct
        key = species.strip().lower()
        if key and key not in counted_species:
            counted_species.add(key)
            captured_pct += pct

    for row in status_rows:
        status = row["status"]
        species = row["gtdb_species"]
        pct = float(row["metaphlan_species_pct"] or 0.0)

        if status in ("matched", "reconstructed_only"):
            genome_id = row["genome_id"]
            src = sources.get(genome_id, {})
            simulated = genome_id in abundance
            if status == "matched" and simulated:
                credit(species, pct)
            out_rows.append({
                "status": status, "genome_id": genome_id,
                "source": src.get("source", "mag"), "accession": "",
                "gtdb_species": species,
                "metaphlan_species_pct": row["metaphlan_species_pct"],
                "rel_abundance": f"{abundance[genome_id]:.6f}" if simulated else "",
                "simulated": "yes" if simulated else "no",
                "note": row["note"] if not simulated else
                        ("" if status == "matched" else row["note"]),
            })
            continue

        # profiled_only: did a reference genome stand in for it?
        ref = ref_by_species.get(species)
        if ref and ref["status"] in ok_status and ref["genome_id"] in abundance:
            n_recovered += 1
            credit(species, pct)
            out_rows.append({
                "status": "profiled_only", "genome_id": ref["genome_id"],
                "source": "reference", "accession": ref["accession"],
                "gtdb_species": species,
                "metaphlan_species_pct": row["metaphlan_species_pct"],
                "rel_abundance": f"{abundance[ref['genome_id']]:.6f}",
                "simulated": "yes",
                "note": "profiled but not assembled; a reference genome of this "
                        "species stands in for it (proxy, not the sample's strain)",
            })
        else:
            n_missing += 1
            reason = ref["status"] if ref else "not_fetched"
            out_rows.append({
                "status": "profiled_only", "genome_id": "", "source": "",
                "accession": ref["accession"] if ref else "",
                "gtdb_species": species,
                "metaphlan_species_pct": row["metaphlan_species_pct"],
                "rel_abundance": "", "simulated": "no",
                "note": f"profiled but not assembled and no reference genome "
                        f"({reason}); NOT part of the simulated community",
            })

    cols = ["status", "genome_id", "source", "accession", "gtdb_species",
            "metaphlan_species_pct", "rel_abundance", "simulated", "note"]
    args.out_reconciliation.parent.mkdir(parents=True, exist_ok=True)
    with args.out_reconciliation.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + cols)
        order = {"matched": 0, "reconstructed_only": 1, "profiled_only": 2}
        for r in sorted(out_rows, key=lambda r: (order[r["status"]],
                                                 -float(r["rel_abundance"] or 0))):
            w.writerow([args.sample] + [r[c] for c in cols])

    total_profiled = float(stats.get("total_profiled_pct") or 0.0)
    pct_captured = (100.0 * captured_pct / total_profiled) if total_profiled > 0 else 0.0
    n_mags = sum(1 for s in sources.values() if s["source"] == "mag")
    n_refs = sum(1 for s in sources.values() if s["source"] == "reference")

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.out_summary.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "n_genomes_simulated", "n_mags", "n_references_added",
                    "n_profiled_species", "n_matched", "n_reconstructed_only",
                    "n_profiled_only", "n_profiled_only_recovered",
                    "n_profiled_only_missing", "pct_profiled_abundance_captured",
                    "metaphlan_unknown_pct", "coverm_unmapped_pct"])
        w.writerow([args.sample, len(abundance), n_mags, n_refs,
                    stats.get("n_profiled_species", 0), stats.get("n_matched", 0),
                    stats.get("n_reconstructed_only", 0),
                    stats.get("n_profiled_only", 0), n_recovered, n_missing,
                    f"{pct_captured:.2f}",
                    f"{float(stats.get('metaphlan_unknown_pct') or 0.0):.2f}",
                    f"{coverm_unmapped_pct(args.coverm):.2f}"])

    print(f"[{args.sample}] simulated={len(abundance)} "
          f"({n_mags} MAG + {n_refs} ref) captured={pct_captured:.1f}% "
          f"profiled_only_missing={n_missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
