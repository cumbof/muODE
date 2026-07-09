#!/usr/bin/env python3
"""Reconcile ONE sample's RECONSTRUCTED genomes with its QUANTITATIVELY PROFILED ones.

Two independent views of the same sample:

  * reconstructed -- the MAGs assembled and binned from this sample's own reads,
    dereplicated within the sample and classified by GTDB-Tk. These are the only
    genomes muODE can simulate, because a GEM is built from a genome.
  * profiled       -- MetaPhlAn4's SGB-level relative abundances, mapped into GTDB
    taxonomy (via the DB's SGB->GTDB table) so both views share one namespace.

They disagree, and the disagreement is informative. A species can be profiled yet
never assemble (too low-abundance, too much microdiversity); a MAG can assemble yet
carry no MetaPhlAn markers (a novel/uncharacterised genome -- MetaPhlAn4 would call
it a uSGB only if it is in the DB). So:

  status=matched             assembled AND profiled -> muODE simulates it
  status=reconstructed_only  assembled, not profiled -> novel / below markers
  status=profiled_only       profiled, not assembled -> swap in a reference genome
                                                        if you need to simulate it

muODE's abundance is written over the INTERSECTION, taking the value from the
quantitative profiler (relative abundance, comparable across samples and studies)
rather than from within-sample assembly coverage. Nothing is dropped silently:
every excluded genome, from either side, lands in the reconciliation table.

Degrades gracefully -- with no profile (or no GTDB-Tk taxonomy) it falls back to
CoverM over all of the sample's MAGs and records that in the summary.

The muODE abundance contract is unchanged: a 2-column TSV `(mag_id, rel_abundance)`
renormalised to sum to 1 over the simulated community.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

_WS = re.compile(r"\s+")


# --- taxonomy helpers -------------------------------------------------------
def norm_species(name: str) -> str:
    """Normalise a GTDB species label for matching."""
    name = name.strip()
    if name.lower().startswith("s__"):
        name = name[3:]
    return _WS.sub(" ", name).strip().lower()


def species_from_lineage(lineage: str) -> str:
    """Text after the last `s__` in a GTDB lineage; "" when species is unpopulated.

    Handles both `;`-separated (GTDB-Tk) and `|`-separated (MetaPhlAn) lineages.
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
    """mag_id -> GTDB lineage (GTDB-Tk classify_wf summary)."""
    out = {}
    if path and path.is_file():
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                mag = (row.get("user_genome") or "").strip()
                if mag:
                    out[mag] = (row.get("classification") or "").strip()
    return out


def load_coverm(path: Path) -> dict:
    """mag_id -> CoverM relative abundance (%); the `unmapped` row is dropped."""
    out = {}
    if path and path.is_file():
        with path.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader, None)  # header
            for row in reader:
                if len(row) < 2 or not row[0]:
                    continue
                mag = row[0].strip()
                if mag.lower() == "unmapped":
                    continue
                try:
                    out[mag] = float(row[1])
                except ValueError:
                    continue
    return out


def load_sgb2gtdb(path: Path) -> dict:
    """SGB id -> GTDB lineage, from the SGB->GTDB table shipped with the MetaPhlAn DB.

    Tolerates a header line and extra columns: the key is the first field and the
    lineage is the last field (a GTDB lineage contains no tabs).
    """
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
                continue  # header row or unusable line
            out[key] = lineage
    return out


def load_metaphlan(path: Path):
    """Parse a MetaPhlAn4 `-t rel_ab` profile.

    Returns (sgb_ab, unknown_pct):
      sgb_ab      -- {SGB id: relative abundance (%)} taken from SGB-level (`t__`)
                     clades only, which is the finest rank MetaPhlAn4 reports;
      unknown_pct -- the UNCLASSIFIED percentage when the profile was produced with
                     `--unclassified_estimation`, else 0.0.
    """
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
                # the column header is a comment: `#clade_name<TAB>NCBI_tax_id<TAB>...`
                header = [c.strip().lstrip("#") for c in line.split("\t")]
                if "relative_abundance" in header:
                    ab_idx = header.index("relative_abundance")
                continue
            fields = line.split("\t")
            clade = fields[0].strip()

            idx = ab_idx
            if idx is None or idx >= len(fields):
                # headerless / unexpected width: MetaPhlAn writes rel_ab in col 2
                # when the NCBI tax id column is present, else col 1.
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
                continue  # a higher-rank aggregate row; SGB rows carry the mass
            sgb = norm_sgb(leaf)
            sgb_ab[sgb] = sgb_ab.get(sgb, 0.0) + value
    return sgb_ab, unknown_pct


def profile_to_species(sgb_ab: dict, sgb2gtdb: dict):
    """Map SGB abundances onto GTDB species. Returns (species_ab, unmapped_sgbs).

    Several SGBs can share one GTDB species (strain-level bins) -- their mass is
    summed. SGBs absent from the mapping table cannot be compared with GTDB-Tk and
    are reported separately rather than silently discarded.
    """
    species_ab: dict = {}
    unmapped = {}
    for sgb, ab in sgb_ab.items():
        lineage = sgb2gtdb.get(sgb)
        species = species_from_lineage(lineage) if lineage else ""
        if not species:
            unmapped[sgb] = ab
            continue
        key = norm_species(species)
        prev = species_ab.get(key)
        species_ab[key] = (prev[0] + ab, prev[1]) if prev else (ab, species)
    return species_ab, unmapped


# --- the join ---------------------------------------------------------------
def reconcile(all_mags: set, mag_lineage: dict, coverm: dict, species_ab: dict):
    """Join MAGs to profiled species on GTDB species.

    Both sides are always fully accounted for: every MAG lands in `matched` or
    `reconstructed_only`, and every profiled species in `matched` or
    `profiled_only` -- including the degenerate cases (no MAGs at all, no
    taxonomy, no profile), where one side is simply empty.

    Returns (matched, reconstructed_only, profiled_only).
    """
    species_to_mags: dict = {}
    mag_species: dict = {}
    unmatchable = []  # MAGs with no species-level GTDB assignment
    for mag in sorted(all_mags):
        raw = species_from_lineage(mag_lineage.get(mag, ""))
        if not raw:
            unmatchable.append(mag)
            continue
        key = norm_species(raw)
        mag_species[mag] = raw
        species_to_mags.setdefault(key, []).append(mag)

    matched, reconstructed_only, profiled_only = [], [], []
    for key, mags in species_to_mags.items():
        if key not in species_ab:
            for mag in mags:
                reconstructed_only.append({
                    "mag_id": mag, "species": mag_species[mag],
                    "coverm": coverm.get(mag, 0.0),
                    "note": "assembled but not profiled "
                            "(novel / below MetaPhlAn markers / DB-version mismatch)"
                    if species_ab else
                    "no quantitative profile for this sample; using CoverM abundance",
                })
            continue
        total, _label = species_ab[key]
        # One species, several MAGs (residual strain-level bins): split the
        # profiled mass by within-sample CoverM coverage; equal split if absent.
        weights = {m: coverm.get(m, 0.0) for m in mags}
        wsum = sum(weights.values())
        for mag in mags:
            share = (weights[mag] / wsum) if wsum > 0 else 1.0 / len(mags)
            matched.append({
                "mag_id": mag, "species": mag_species[mag],
                "metaphlan": total * share, "coverm": coverm.get(mag, 0.0),
            })

    for mag in unmatchable:
        reconstructed_only.append({
            "mag_id": mag, "species": "", "coverm": coverm.get(mag, 0.0),
            "note": "no species-level GTDB assignment; cannot be matched to a profile entry"
            if mag_lineage else
            "no GTDB-Tk taxonomy for this sample; cannot be matched to a profile entry",
        })

    matched_keys = {norm_species(d["species"]) for d in matched}
    for key, (ab, label) in species_ab.items():
        if key not in matched_keys:
            profiled_only.append({"species": label, "metaphlan": ab})

    matched.sort(key=lambda d: d["mag_id"])
    reconstructed_only.sort(key=lambda d: d["mag_id"])
    profiled_only.sort(key=lambda d: -d["metaphlan"])
    return matched, reconstructed_only, profiled_only


# --- writers ----------------------------------------------------------------
def write_abundance(path: Path, pairs):
    """pairs: [(mag_id, mass)]. Renormalise to sum to 1 -> muODE's contract."""
    total = sum(m for _, m in pairs)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mag_id", "rel_abundance"])
        for mag_id, mass in pairs:
            w.writerow([mag_id, f"{(mass / total) if total > 0 else 0.0:.6f}"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--coverm", type=Path, help="CoverM relative_abundance TSV")
    ap.add_argument("--gtdbtk", type=Path, help="GTDB-Tk summary TSV for this sample")
    ap.add_argument("--profile", type=Path, help="MetaPhlAn4 rel_ab profile")
    ap.add_argument("--sgb2gtdb", type=Path,
                    help="SGB->GTDB table shipped with the MetaPhlAn DB")
    ap.add_argument("--source", choices=["reconciled", "coverm"], default="reconciled",
                    help="which abundance muODE consumes (default: reconciled)")
    ap.add_argument("--out-abundance", required=True, type=Path)
    ap.add_argument("--out-reconciliation", required=True, type=Path)
    ap.add_argument("--out-summary", type=Path)
    args = ap.parse_args(argv)

    mag_lineage = load_gtdbtk(args.gtdbtk) if args.gtdbtk else {}
    coverm = load_coverm(args.coverm) if args.coverm else {}
    sgb_ab, unknown_pct = load_metaphlan(args.profile) if args.profile else ({}, 0.0)
    sgb2gtdb = load_sgb2gtdb(args.sgb2gtdb) if args.sgb2gtdb else {}
    species_ab, unmapped_sgbs = profile_to_species(sgb_ab, sgb2gtdb)

    all_mags = set(mag_lineage) | set(coverm)
    matched, reconstructed_only, profiled_only = reconcile(
        all_mags, mag_lineage, coverm, species_ab)

    # --- pick the abundance muODE consumes --------------------------------
    # Every fallback is named in the summary rather than silently taken: a bundle
    # built on CoverM is a different scientific object from a reconciled one.
    if not all_mags:
        used = "coverm(empty-catalogue)"  # no MAGs recovered; nothing to simulate
    elif not species_ab:
        used = ("coverm(fallback:no-sgb-gtdb-mapping)" if sgb_ab
                else "coverm")  # sgb_ab non-empty means --sgb2gtdb was wrong/missing
    elif not mag_lineage:
        used = "coverm(fallback:no-gtdbtk-taxonomy)"
    elif args.source == "coverm":
        used = "coverm"
    elif not matched:
        # Profiling ran but nothing intersected. Handing muODE an empty community
        # would be worse than falling back -- do so loudly.
        used = "coverm(fallback:empty-intersection)"
    else:
        used = "reconciled"

    if used.startswith("coverm"):
        pairs = [(m, coverm.get(m, 0.0)) for m in sorted(all_mags)]
    else:
        pairs = [(d["mag_id"], d["metaphlan"]) for d in matched]
    write_abundance(args.out_abundance, pairs)

    # --- the reconciliation table -----------------------------------------
    args.out_reconciliation.parent.mkdir(parents=True, exist_ok=True)
    with args.out_reconciliation.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "status", "mag_id", "gtdb_species",
                    "metaphlan_pct", "coverm_pct", "note"])
        for d in matched:
            w.writerow([args.sample, "matched", d["mag_id"], d["species"],
                        f"{d['metaphlan']:.6f}", f"{d['coverm']:.6f}", ""])
        for d in reconstructed_only:
            w.writerow([args.sample, "reconstructed_only", d["mag_id"], d["species"],
                        "", f"{d['coverm']:.6f}", d["note"]])
        for d in profiled_only:
            w.writerow([args.sample, "profiled_only", "", d["species"],
                        f"{d['metaphlan']:.6f}", "",
                        "profiled but not assembled; substitute a reference genome "
                        "to simulate it"])
        for sgb, ab in sorted(unmapped_sgbs.items(), key=lambda kv: -kv[1]):
            w.writerow([args.sample, "profiled_only", "", sgb, f"{ab:.6f}", "",
                        "profiled SGB with no GTDB mapping; cannot be matched to a MAG"])

    # --- one-row summary (concatenated across samples downstream) ----------
    total_profiled = sum(ab for ab, _ in species_ab.values()) + sum(unmapped_sgbs.values())
    captured = sum(d["metaphlan"] for d in matched)
    pct_captured = (100.0 * captured / total_profiled) if total_profiled > 0 else 0.0
    if args.out_summary:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        with args.out_summary.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "n_mags", "n_profiled_species", "n_matched",
                        "n_reconstructed_only", "n_profiled_only",
                        "pct_profiled_abundance_captured", "metaphlan_unknown_pct",
                        "muode_abundance_source"])
            w.writerow([args.sample, len(all_mags), len(species_ab), len(matched),
                        len(reconstructed_only),
                        len(profiled_only) + len(unmapped_sgbs),
                        f"{pct_captured:.2f}", f"{unknown_pct:.2f}", used])

    print(f"[{args.sample}] matched={len(matched)} "
          f"reconstructed_only={len(reconstructed_only)} "
          f"profiled_only={len(profiled_only) + len(unmapped_sgbs)} "
          f"captured={pct_captured:.1f}% source={used}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
