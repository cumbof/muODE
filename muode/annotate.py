"""Call functional traits from MAG proteomes, and hand the guilds to the ecology layers.

This closes the gap described in :mod:`muode.markers`: the ecology layers all take
a *set of species ids*, and until now nothing produced those sets from genomes, so
they were hand-written in the designed scenarios and empty on real data.

The pipeline is deliberately split in two, because only one half can be tested on
a machine without HMMER:

* **search** (:func:`hmmsearch`, :func:`diamond_blastp`) -- thin shell-outs.  Not
  unit-tested; they run where the tools are installed.
* **calling** (:func:`parse_hmmsearch_tblout`, :func:`parse_diamond_tab`,
  :func:`muode.markers.call_all`) -- pure functions over a hit table.  These carry
  all the judgement, so these are the ones under test.

Usage::

    python -m muode.annotate --proteomes proteins/ --refs bai_refs.faa \\
        --pfam Pfam-A.hmm --out traits.tsv

``traits.tsv`` is the same side file the Snakemake pipeline already reads for
domain routing (``traits:`` in config); this adds columns to it, and the existing
``csv.DictReader`` ignores the ones it does not know.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Mapping, Set, Tuple

from muode.markers import TRAITS, MarkerSource, accessions, call_all

#: DIAMOND defaults.  The bitscore is a judgement call and we say so: it is
#: deliberately permissive, because the evidence bar for a REFPROT trait is
#: `min_markers` (three separate bai proteins), not the score of any single hit.
#: Tightening this makes bai calls *more* conservative, never less.
MIN_BITSCORE = 60.0
MIN_IDENTITY = 30.0


# -- search (shell-outs; run where the tools live) ---------------------------

def hmmsearch(proteome: Path, pfam_hmm: Path, out: Path) -> Path:
    """``hmmsearch --cut_ga`` one proteome against Pfam-A.

    ``--cut_ga`` uses Pfam's own curated gathering threshold, so the cutoff is
    theirs and not a number we invented.
    """
    subprocess.run(
        ["hmmsearch", "--cut_ga", "--noali", "--tblout", str(out),
         str(pfam_hmm), str(proteome)],
        check=True, stdout=subprocess.DEVNULL,
    )
    return out


def diamond_blastp(proteome: Path, refs_dmnd: Path, out: Path) -> Path:
    """DIAMOND one proteome against the reference proteins (bai operon)."""
    subprocess.run(
        ["diamond", "blastp", "--quiet", "--query", str(proteome),
         "--db", str(refs_dmnd), "--outfmt", "6", "qseqid", "sseqid", "pident",
         "bitscore", "--out", str(out), "--max-target-seqs", "5"],
        check=True,
    )
    return out


# -- parsing + calling (pure; these are the tested ones) ---------------------

def parse_hmmsearch_tblout(text: str) -> Set[str]:
    """Pfam accessions that hit, from ``hmmsearch --tblout``.

    Column 4 is the query (HMM) accession, e.g. ``PF02275.20``.  The version
    suffix is dropped: the registry pins families, not releases.
    """
    hits: Set[str] = set()
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 4:
            continue
        acc = cols[3]
        if acc == "-":                      # some HMM libraries carry no accession
            continue
        hits.add(acc.split(".")[0])
    return hits


def parse_diamond_tab(text: str, min_bitscore: float = MIN_BITSCORE,
                      min_identity: float = MIN_IDENTITY) -> Set[str]:
    """Reference accessions that hit, from ``diamond --outfmt 6 qseqid sseqid pident bitscore``.

    The subject id is the reference protein.  DIAMOND writes it as whatever the
    FASTA header held, so ``sp|P19412|BAIE_CLOSD`` and a bare ``P19412`` both have
    to resolve to ``P19412``.
    """
    hits: Set[str] = set()
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        _q, subject, pident, bitscore = cols[0], cols[1], cols[2], cols[3]
        try:
            if float(bitscore) < min_bitscore or float(pident) < min_identity:
                continue
        except ValueError:
            continue
        parts = [p for p in subject.split("|") if p]
        hits.add(parts[1] if len(parts) >= 3 else parts[-1])
    return hits


def call(hits_by_mag: Mapping[str, FrozenSet]) -> Dict[str, Dict[str, bool]]:
    """``{mag: {accession}}`` -> ``{mag: {trait: present}}``.  See markers.call_all."""
    return call_all(hits_by_mag)


# -- guilds: the whole point ------------------------------------------------

def guilds(traits: Mapping[str, Mapping[str, bool]]) -> Dict[str, Set[str]]:
    """``{mag: {trait: present}}`` -> ``{trait: {mag, ...}}``.

    This is the shape every ecology layer wants::

        g = guilds(call(hits))
        BileAcidTransform(bsh_producers=g["bsh"], bai_producers=g["bai"])

    A trait nobody carries yields an empty set, and an empty set is a *finding*:
    no bai carrier means the community cannot 7a-dehydroxylate, and no modelling
    choice can conjure deoxycholate out of a community that lacks the operon.
    """
    out: Dict[str, Set[str]] = {t.name: set() for t in TRAITS}
    for mag, called in traits.items():
        for name, present in called.items():
            if present:
                out.setdefault(name, set()).add(mag)
    return out


# -- traits.tsv -------------------------------------------------------------

def write_traits_tsv(path: str | Path, traits: Mapping[str, Mapping[str, bool]],
                     domains: Mapping[str, str] | None = None) -> Path:
    """Write the side file the Snakemake pipeline already reads.

    Keeps ``mag_id`` and ``domain`` first so an existing domain-routing run keeps
    working; the trait columns are additive and older readers ignore them.
    """
    path = Path(path)
    domains = domains or {}
    names = [t.name for t in TRAITS]
    with path.open("w", newline="") as fh:
        fh.write("# muODE functional traits -- generated by muode.annotate\n")
        fh.write("# Traits the GEM cannot encode; they drive the ecology layers.\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mag_id", "domain", *names])
        for mag in sorted(traits):
            called = traits[mag]
            w.writerow([mag, domains.get(mag, "bacteria"),
                        *[int(bool(called.get(n, False))) for n in names]])
    return path


def read_traits_tsv(path: str | Path) -> Dict[str, Dict[str, bool]]:
    """Read trait columns back.  Unknown columns are ignored; missing ones are False."""
    names = {t.name for t in TRAITS}
    out: Dict[str, Dict[str, bool]] = {}
    with Path(path).open() as fh:
        reader = csv.DictReader(
            (ln for ln in fh if not ln.strip().startswith("#")), delimiter="\t"
        )
        for row in reader:
            mag = row["mag_id"]
            out[mag] = {
                n: str(row.get(n, "0")).strip() in {"1", "true", "True", "yes"}
                for n in names
            }
    return out


# -- CLI --------------------------------------------------------------------

def _mag_id(proteome: Path) -> str:
    return proteome.name.split(".")[0]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Call functional traits from MAG proteomes")
    ap.add_argument("--proteomes", required=True,
                    help="directory of per-MAG protein FASTAs (*.faa)")
    ap.add_argument("--pfam", help="Pfam-A.hmm (press'd); enables the PFAM traits")
    ap.add_argument("--refs", help="DIAMOND db of reference proteins; enables the REFPROT traits")
    ap.add_argument("--out", default="traits.tsv")
    ap.add_argument("--workdir", default=".muode_annotate")
    args = ap.parse_args(argv)

    if not args.pfam and not args.refs:
        ap.error("give --pfam and/or --refs, or nothing can be called")

    proteomes = sorted(Path(args.proteomes).glob("*.faa"))
    if not proteomes:
        sys.exit(f"no *.faa under {args.proteomes}")
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)

    wanted_pfam = set(accessions(MarkerSource.PFAM))
    wanted_ref = set(accessions(MarkerSource.REFPROT))

    hits: Dict[str, FrozenSet] = {}
    for p in proteomes:
        mag = _mag_id(p)
        found: Set[str] = set()
        if args.pfam:
            tbl = hmmsearch(p, Path(args.pfam), work / f"{mag}.pfam.tbl")
            found |= parse_hmmsearch_tblout(tbl.read_text()) & wanted_pfam
        if args.refs:
            tab = diamond_blastp(p, Path(args.refs), work / f"{mag}.refs.tsv")
            found |= parse_diamond_tab(tab.read_text()) & wanted_ref
        hits[mag] = frozenset(found)

    traits = call(hits)
    write_traits_tsv(args.out, traits)
    g = guilds(traits)

    n = len(proteomes)
    print(f"{n} MAGs -> {args.out}\n")
    for t in TRAITS:
        members = sorted(g.get(t.name, ()))
        print(f"  {t.name:16} {len(members):>3}/{n}  -> {t.layer}")
        if not members:
            print("      NONE.  This community cannot do it; that is a result, not a gap.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
