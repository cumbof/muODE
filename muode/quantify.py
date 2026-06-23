"""Quantitative profiling: join MetaSBT taxonomy with Bracken abundance.

This is the bridge between the two upstream contracts:

* :mod:`muode.metasbt` gives each MAG its **taxonomic identity** -- the closest
  MetaSBT species cluster (and full lineage).
* :mod:`muode.bracken` gives each species cluster its **relative abundance** in
  a sample.

:func:`mag_abundance_from_bracken` joins them into the ``{mag_id: rel_abundance}``
contract the rest of muODE consumes (``assemble``/``simulate``).  The join key is
the MAG's MetaSBT cluster label at a chosen rank (species by default), matched to
Bracken's ``name``; matching is exact first, then on a normalised leaf token so
``s__Escherichia_coli`` and ``Escherichia coli`` reconcile.

It also provides thin wrappers that build (and optionally run) the heavy upstream
commands -- ``MetaSBT profile`` / ``MetaSBT kraken``, ``bracken-build``,
``kraken2`` and ``bracken``.  The command *builders* (``*_cmd``) are pure and
testable; the runners shell out via :func:`muode.external.run` and are meant for
a capable workstation, **not** the aarch64 dev box.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from muode.bracken import BrackenProfile
from muode.external import require, run
from muode.metasbt import MetaSBTProfile

# ---------------------------------------------------------------------------
# label normalisation + the join
# ---------------------------------------------------------------------------

_RANK_PREFIX = re.compile(r"^[a-z]__")


def normalize_label(label: str) -> str:
    """Reduce a taxon label to a comparable leaf token.

    Takes the last component of a pipe-separated lineage, strips a ``s__``-style
    rank prefix, lower-cases and collapses spaces/underscores -- so MetaSBT's
    ``...|s__Escherichia_coli`` matches Bracken's ``Escherichia coli``.
    """
    leaf = str(label).split("|")[-1].strip()
    leaf = _RANK_PREFIX.sub("", leaf)
    return re.sub(r"[\s_]+", " ", leaf).strip().lower()


def match_profiles_to_bracken(
    profiles: Dict[str, MetaSBTProfile],
    bracken: BrackenProfile,
    level: str = "species",
) -> Dict[str, Optional[str]]:
    """Map ``{mag_id: matched Bracken name}`` (``None`` where no taxon matched)."""
    norm_index: Dict[str, str] = {}
    for name in bracken.abundances:
        norm_index.setdefault(normalize_label(name), name)

    result: Dict[str, Optional[str]] = {}
    for mag_id, prof in profiles.items():
        label = prof.label(level)
        if not label:
            result[mag_id] = None
            continue
        if label in bracken.abundances:           # exact name match
            result[mag_id] = label
        else:                                       # normalised leaf match
            result[mag_id] = norm_index.get(normalize_label(label))
    return result


def mag_abundance_from_bracken(
    bracken: BrackenProfile,
    profiles: Dict[str, MetaSBTProfile],
    level: str = "species",
    split_shared: bool = True,
    normalize: bool = True,
) -> Dict[str, float]:
    """Join MetaSBT identity with Bracken abundance into ``{mag_id: abundance}``.

    Each MAG inherits the abundance of the Bracken taxon its MetaSBT ``level``
    cluster matches.  When several MAGs map to the *same* taxon (e.g. multiple
    bins of one species cluster) and ``split_shared`` is set, that taxon's
    abundance is divided equally among them so the community total is conserved.
    MAGs with no matching taxon get ``0.0``.  With ``normalize`` the result is
    rescaled to sum to 1 over the matched MAGs.
    """
    matched = match_profiles_to_bracken(profiles, bracken, level=level)

    shares: Dict[str, int] = {}
    if split_shared:
        for name in matched.values():
            if name is not None:
                shares[name] = shares.get(name, 0) + 1

    abundances: Dict[str, float] = {}
    for mag_id, name in matched.items():
        if name is None:
            abundances[mag_id] = 0.0
            continue
        value = bracken.abundances.get(name, 0.0)
        if split_shared and shares.get(name, 1) > 1:
            value /= shares[name]
        abundances[mag_id] = value

    total = sum(abundances.values())
    if normalize and total > 0:
        abundances = {m: v / total for m, v in abundances.items()}
    return abundances


# ---------------------------------------------------------------------------
# Upstream tool command scaffolds (build then optionally run).
#
# These construct standard MetaSBT / Kraken2 / Bracken invocations.  They are
# resource-intensive and are intended for a capable workstation; only the
# command *builders* are exercised in CI.  `extra_args` is passed through so you
# can supply version-specific flags without editing muODE.
# ---------------------------------------------------------------------------


def _flatten(*parts) -> List[str]:
    out: List[str] = []
    for p in parts:
        if p is None:
            continue
        if isinstance(p, (list, tuple)):
            out.extend(str(x) for x in p)
        else:
            out.append(str(p))
    return out


def metasbt_profile_cmd(
    input_file: str | Path,
    database: str | Path,
    output_dir: str | Path,
    *,
    expand: bool = True,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    """``MetaSBT profile`` -- characterize genomes against a MetaSBT database.

    Writes one ``{genome}.txt`` per input genome under ``output_dir`` (parsed by
    :func:`muode.metasbt.read_metasbt_profiles`).
    """
    cmd = _flatten(
        "MetaSBT", "profile",
        "--input-file", input_file,
        "--tree", database,
        "--output-dir", output_dir,
    )
    if expand:
        cmd.append("--expand")
    return cmd + list(extra_args or [])


def metasbt_kraken_cmd(
    database: str | Path,
    output_db: str | Path,
    *,
    genomes: Optional[str | Path] = None,
    names_dmp: Optional[str | Path] = None,
    nodes_dmp: Optional[str | Path] = None,
    threads: int = 4,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    """``MetaSBT kraken`` -- build a custom Kraken2 DB from a MetaSBT database.

    The Kraken2 DB inherits MetaSBT's known/unknown species clusters, so reads
    classify directly to the clusters MetaSBT characterizes.
    """
    cmd = _flatten(
        "MetaSBT", "kraken",
        "--db-dir", database,
        "--output-dir", output_db,
        "--threads", threads,
    )
    cmd += _flatten("--genomes", genomes) if genomes else []
    cmd += _flatten("--names", names_dmp) if names_dmp else []
    cmd += _flatten("--nodes", nodes_dmp) if nodes_dmp else []
    return cmd + list(extra_args or [])


def bracken_build_cmd(
    kraken_db: str | Path,
    *,
    kmer_len: int = 35,
    read_len: int = 150,
    threads: int = 4,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    """``bracken-build`` -- build the Bracken DB on top of a Kraken2 DB."""
    cmd = _flatten(
        "bracken-build",
        "-d", kraken_db,
        "-k", kmer_len,
        "-l", read_len,
        "-t", threads,
    )
    return cmd + list(extra_args or [])


def kraken2_cmd(
    kraken_db: str | Path,
    reads: Sequence[str | Path],
    report: str | Path,
    *,
    output: Optional[str | Path] = None,
    paired: bool = False,
    threads: int = 4,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    """``kraken2`` -- classify a sample's reads against the custom DB."""
    cmd = _flatten(
        "kraken2",
        "--db", kraken_db,
        "--threads", threads,
        "--report", report,
        "--output", output if output is not None else "-",
    )
    if paired:
        cmd.append("--paired")
    return cmd + list(extra_args or []) + _flatten(reads)


def bracken_cmd(
    kraken_db: str | Path,
    kraken_report: str | Path,
    output: str | Path,
    *,
    read_len: int = 150,
    level: str = "S",
    threshold: int = 10,
    extra_args: Optional[Sequence[str]] = None,
) -> List[str]:
    """``bracken`` -- re-estimate abundances from a Kraken2 report.

    The output TSV is read by :func:`muode.bracken.read_bracken`.
    """
    cmd = _flatten(
        "bracken",
        "-d", kraken_db,
        "-i", kraken_report,
        "-o", output,
        "-r", read_len,
        "-l", level,
        "-t", threshold,
    )
    return cmd + list(extra_args or [])


# -- runners (shell out; for a capable workstation, not aarch64) -------------

def _run_built(binary: str, cmd: List[str], log: Optional[str | Path]):
    require(binary, env_hint="see the muODE upstream docs (MetaSBT/Kraken2/Bracken)")
    return run(cmd, log=log)


def metasbt_profile(input_file, database, output_dir, *, log=None, **kw):
    cmd = metasbt_profile_cmd(input_file, database, output_dir, **kw)
    return _run_built("MetaSBT", cmd, log)


def metasbt_kraken(database, output_db, *, log=None, **kw):
    cmd = metasbt_kraken_cmd(database, output_db, **kw)
    return _run_built("MetaSBT", cmd, log)


def bracken_build(kraken_db, *, log=None, **kw):
    cmd = bracken_build_cmd(kraken_db, **kw)
    return _run_built("bracken-build", cmd, log)


def run_kraken2(kraken_db, reads, report, *, log=None, **kw):
    cmd = kraken2_cmd(kraken_db, reads, report, **kw)
    return _run_built("kraken2", cmd, log)


def run_bracken(kraken_db, kraken_report, output, *, log=None, **kw):
    cmd = bracken_cmd(kraken_db, kraken_report, output, **kw)
    return _run_built("bracken", cmd, log)
