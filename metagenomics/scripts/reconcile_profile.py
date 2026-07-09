#!/usr/bin/env python3
"""Classify ONE sample's genomes as reconstructed, profiled, or both.

Two independent views of the same sample:

  * reconstructed -- the MAGs assembled and binned from this sample's own reads,
    dereplicated within the sample and classified by GTDB-Tk;
  * profiled       -- MetaPhlAn4's SGB-level relative abundances, mapped into GTDB
    taxonomy (via the DB's SGB->GTDB table) so both views share one namespace.

They disagree, and the disagreement is the point. A species can be profiled from a
handful of marker reads yet never assemble (binning needs roughly 5-10x even
coverage; low-abundance and high-microdiversity taxa fail). A novel MAG can
assemble yet carry no markers.

This script only DECIDES WHICH GENOMES EXIST, joining the two sets on GTDB species:

    matched             assembled AND profiled
    reconstructed_only  assembled, not profiled   (novel / below markers)
    profiled_only       profiled, never assembled -> a reference genome can stand
                                                     in for it (fetch_reference_genomes.py)

It deliberately computes NO abundance. Quantification happens once, downstream,
when CoverM maps the sample's reads against the final genome set (MAGs + fetched
references). Mixing CoverM percentages with MetaPhlAn percentages would mix a
read-fraction with a marker-normalised cell fraction -- two different quantities.

Degrades to "everything is reconstructed_only" when there is no profile.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

_WS = re.compile(r"\s+")


# --- taxonomy helpers -------------------------------------------------------
def norm_species(name: str) -> str:
    name = name.strip()
    if name.lower().startswith("s__"):
        name = name[3:]
    return _WS.sub(" ", name).strip().lower()


def species_from_lineage(lineage: str) -> str:
    """Text after the last `s__`; "" when the species rank is unpopulated.

    Handles `;`-separated (GTDB-Tk) and `|`-separated (MetaPhlAn) lineages.
    """
    species = ""
    for token in re.split(r"[;|]", lineage):
        token = token.strip()
        if token.startswith("s__"):
            species = token[3:].strip()
    return species


def norm_sgb(sgb: str) -> str:
    """`t__SGB10068_group` / `SGB10068` -> `SGB10068`."""
    sgb = sgb.strip()
    if sgb.startswith("t__"):
        sgb = sgb[3:]
    if sgb.endswith("_group"):
        sgb = sgb[: -len("_group")]
    return sgb


# --- loaders ----------------------------------------------------------------
def load_gtdbtk(path: Path) -> dict:
    """mag_id -> GTDB lineage."""
    out = {}
    if path and path.is_file():
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                mag = (row.get("user_genome") or "").strip()
                if mag:
                    out[mag] = (row.get("classification") or "").strip()
    return out


def load_sgb2gtdb(path: Path) -> dict:
    """SGB id -> GTDB lineage (the table shipped with the MetaPhlAn DB)."""
    out = {}
    if not (path and path.is_file()):
        return out
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            key, lineage = norm_sgb(fields[0]), fields[-1].strip()
            if not key.upper().startswith("SGB") or "__" not in lineage:
                continue  # header or unusable row
            out[key] = lineage
    return out


def load_metaphlan(path: Path):
    """Parse a MetaPhlAn4 `-t rel_ab` profile -> ({SGB: pct}, unknown_pct)."""
    sgb_ab: dict = {}
    unknown_pct = 0.0
    if not (path and path.is_file()):
        return sgb_ab, unknown_pct

    ab_idx = None
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = [c.strip().lstrip("#") for c in line.split("\t")]
                if "relative_abundance" in header:
                    ab_idx = header.index("relative_abundance")
                continue
            fields = line.split("\t")
            clade = fields[0].strip()

            idx = ab_idx
            if idx is None or idx >= len(fields):
                idx = 2 if len(fields) >= 3 else 1
            if idx >= len(fields):
                continue
            try:
                value = float(fields[idx])
            except ValueError:
                continue

            if clade.upper() in {"UNCLASSIFIED", "UNKNOWN"}:
                unknown_pct += value
                continue
            leaf = clade.split("|")[-1].strip()
            if not leaf.startswith("t__"):
                continue  # higher-rank aggregate; the SGB rows carry the mass
            sgb = norm_sgb(leaf)
            sgb_ab[sgb] = sgb_ab.get(sgb, 0.0) + value
    return sgb_ab, unknown_pct


def profile_to_species(sgb_ab: dict, sgb2gtdb: dict):
    """SGB abundances -> {norm_species: (pct, label, lineage)}, plus unmapped SGBs.

    Several SGBs can share one GTDB species; their mass is summed. SGBs with no
    GTDB mapping (notably eukaryotic SGBs -- GTDB is prokaryote-only) cannot be
    compared with GTDB-Tk and are reported rather than silently dropped.
    """
    species_ab: dict = {}
    unmapped: dict = {}
    for sgb, ab in sgb_ab.items():
        lineage = sgb2gtdb.get(sgb)
        species = species_from_lineage(lineage) if lineage else ""
        if not species:
            unmapped[sgb] = ab
            continue
        key = norm_species(species)
        if key in species_ab:
            pct, label, lin = species_ab[key]
            species_ab[key] = (pct + ab, label, lin)
        else:
            species_ab[key] = (ab, species, lineage)
    return species_ab, unmapped


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--catalogue", required=True, type=Path,
                    help="dir of this sample's dereplicated MAG FASTAs")
    ap.add_argument("--ext", default="fa")
    ap.add_argument("--gtdbtk", type=Path, help="GTDB-Tk summary TSV for this sample")
    ap.add_argument("--profile", type=Path, help="MetaPhlAn4 rel_ab profile")
    ap.add_argument("--sgb2gtdb", type=Path, help="SGB->GTDB table from the MetaPhlAn DB")
    ap.add_argument("--out-status", required=True, type=Path)
    ap.add_argument("--out-missing", required=True, type=Path,
                    help="profiled_only species -> candidates for a reference genome")
    ap.add_argument("--out-stats", required=True, type=Path)
    args = ap.parse_args(argv)

    mags = sorted(p.stem for p in args.catalogue.glob(f"*.{args.ext}"))
    mag_lineage = load_gtdbtk(args.gtdbtk) if args.gtdbtk else {}
    sgb_ab, unknown_pct = load_metaphlan(args.profile) if args.profile else ({}, 0.0)
    sgb2gtdb = load_sgb2gtdb(args.sgb2gtdb) if args.sgb2gtdb else {}
    species_ab, unmapped_sgbs = profile_to_species(sgb_ab, sgb2gtdb)

    # index MAGs by GTDB species
    species_to_mags: dict = {}
    mag_species: dict = {}
    unmatchable = []
    for mag in mags:
        raw = species_from_lineage(mag_lineage.get(mag, ""))
        if not raw:
            unmatchable.append(mag)
            continue
        mag_species[mag] = raw
        species_to_mags.setdefault(norm_species(raw), []).append(mag)

    rows = []
    matched_keys = set()
    for key, group in sorted(species_to_mags.items()):
        hit = key in species_ab
        if hit:
            matched_keys.add(key)
        pct = f"{species_ab[key][0]:.6f}" if hit else ""
        for mag in group:
            rows.append({
                "status": "matched" if hit else "reconstructed_only",
                "genome_id": mag, "gtdb_species": mag_species[mag],
                "gtdb_lineage": mag_lineage.get(mag, ""), "metaphlan_species_pct": pct,
                "note": "" if hit else (
                    "assembled but not profiled (novel / below MetaPhlAn markers / "
                    "DB-version skew)" if species_ab else
                    "no quantitative profile for this sample"),
            })
    for mag in unmatchable:
        rows.append({
            "status": "reconstructed_only", "genome_id": mag, "gtdb_species": "",
            "gtdb_lineage": mag_lineage.get(mag, ""), "metaphlan_species_pct": "",
            "note": ("no species-level GTDB assignment; cannot be matched to a profile"
                     if mag_lineage else "no GTDB-Tk taxonomy for this sample"),
        })

    missing = []
    for key, (pct, label, lineage) in species_ab.items():
        if key in matched_keys:
            continue
        missing.append({"gtdb_species": label, "gtdb_lineage": lineage,
                        "metaphlan_species_pct": pct})
        rows.append({
            "status": "profiled_only", "genome_id": "", "gtdb_species": label,
            "gtdb_lineage": lineage, "metaphlan_species_pct": f"{pct:.6f}",
            "note": "profiled but not assembled; a reference genome can stand in for it",
        })
    for sgb, pct in unmapped_sgbs.items():
        rows.append({
            "status": "profiled_only", "genome_id": "", "gtdb_species": sgb,
            "gtdb_lineage": "", "metaphlan_species_pct": f"{pct:.6f}",
            "note": "profiled SGB with no GTDB mapping (e.g. a eukaryotic SGB); "
                    "cannot be matched to a MAG or resolved to a reference genome",
        })

    # NB `metaphlan_species_pct` is the abundance of the genome's SPECIES, not of
    # the genome. Two MAGs of one species both carry the species' percentage, so
    # summing this column over genomes double-counts. Aggregate by species.
    args.out_status.parent.mkdir(parents=True, exist_ok=True)
    cols = ["status", "genome_id", "gtdb_species", "gtdb_lineage",
            "metaphlan_species_pct", "note"]
    with args.out_status.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + cols)
        order = {"matched": 0, "reconstructed_only": 1, "profiled_only": 2}
        for r in sorted(rows, key=lambda r: (order[r["status"]], r["genome_id"],
                                             -float(r["metaphlan_species_pct"] or 0))):
            w.writerow([args.sample] + [r[c] for c in cols])

    missing.sort(key=lambda d: -d["metaphlan_species_pct"])
    args.out_missing.parent.mkdir(parents=True, exist_ok=True)
    with args.out_missing.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gtdb_species", "gtdb_lineage", "metaphlan_species_pct"])
        for d in missing:
            w.writerow([d["gtdb_species"], d["gtdb_lineage"],
                        f"{d['metaphlan_species_pct']:.6f}"])

    total_profiled = (sum(p for p, _, _ in species_ab.values())
                      + sum(unmapped_sgbs.values()))
    n_matched = sum(1 for r in rows if r["status"] == "matched")
    n_recon = sum(1 for r in rows if r["status"] == "reconstructed_only")
    args.out_stats.parent.mkdir(parents=True, exist_ok=True)
    with args.out_stats.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "n_mags", "n_profiled_species", "n_matched",
                    "n_reconstructed_only", "n_profiled_only",
                    "total_profiled_pct", "metaphlan_unknown_pct"])
        w.writerow([args.sample, len(mags), len(species_ab), n_matched, n_recon,
                    len(missing) + len(unmapped_sgbs), f"{total_profiled:.6f}",
                    f"{unknown_pct:.6f}"])

    print(f"[{args.sample}] mags={len(mags)} matched={n_matched} "
          f"reconstructed_only={n_recon} "
          f"profiled_only={len(missing) + len(unmapped_sgbs)} "
          f"(fetchable={len(missing)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
