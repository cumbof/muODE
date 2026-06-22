"""Phase 1 -- automated GEM reconstruction from MAGs.

Thin, well-defined wrappers around the established reconstruction tools.  muODE
defaults to **CarveMe** because it is fast (seconds-to-minutes per genome) and
emits models in a single consistent namespace (BiGG), which is a hard
requirement for assembling many MAGs into one community with a shared metabolite
pool.  **gapseq** is offered as an alternative engine (more comprehensive
pathway evidence, ModelSEED namespace, but far slower) -- it should not be mixed
with CarveMe models in the same community.

The actual binaries run inside the workflow's per-rule conda environments; here
we only construct and execute the commands.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from muode.external import MuodeToolError, require, run


def call_genes(
    genome_fna: str | Path,
    proteins_faa: str | Path,
    translation_table: int = 11,
) -> Path:
    """Predict protein-coding genes from a (MAG) nucleotide FASTA.

    Prefers the ``pyrodigal`` Python binding (fast, no subprocess); falls back to
    a ``prodigal`` binary.  Metagenomic mode is used since MAGs are drafts.
    """
    genome_fna, proteins_faa = Path(genome_fna), Path(proteins_faa)
    proteins_faa.parent.mkdir(parents=True, exist_ok=True)
    try:
        import pyrodigal
        from Bio import SeqIO  # type: ignore

        orf_finder = pyrodigal.GeneFinder(meta=True)
        with proteins_faa.open("w") as out:
            for record in SeqIO.parse(str(genome_fna), "fasta"):
                genes = orf_finder.find_genes(bytes(record.seq))
                genes.write_translations(out, sequence_id=record.id)
        return proteins_faa
    except ImportError:
        prodigal = require("prodigal", env_hint="mamba install -c bioconda prodigal")
        run([prodigal, "-i", genome_fna, "-a", proteins_faa, "-p", "meta",
             "-g", str(translation_table), "-q"])
        return proteins_faa


def carveme(
    proteins_faa: str | Path,
    output_xml: str | Path,
    universe: str = "bacteria",
    gapfill_media: Optional[str] = None,
    solver: Optional[str] = None,
    from_dna: bool = False,
) -> Path:
    """Reconstruct a GEM with CarveMe (top-down, BiGG namespace)."""
    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    carve = require("carve", env_hint="mamba install -c bioconda carveme")
    cmd = [carve, str(proteins_faa), "-o", str(output_xml), "-u", universe]
    if from_dna:
        cmd.append("--dna")
    if gapfill_media:
        cmd += ["-g", gapfill_media]
    if solver:
        cmd += ["--solver", solver]
    run(cmd, log=output_xml.with_suffix(".carveme.log"))
    return output_xml


def gapseq(
    genome_fna: str | Path,
    outdir: str | Path,
    media: Optional[str] = None,
) -> Path:
    """Reconstruct a GEM with gapseq (bottom-up, ModelSEED namespace).

    Returns the path to the gap-filled SBML model produced by ``gapseq doall``.
    """
    genome_fna, outdir = Path(genome_fna), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gs = require("gapseq", env_hint="see https://github.com/jotech/gapseq")
    cmd = [gs, "doall", str(genome_fna)]
    if media:
        cmd += [str(media)]
    run(cmd, log=outdir / f"{genome_fna.stem}.gapseq.log")
    produced = genome_fna.with_suffix("").name
    candidates = list(Path.cwd().glob(f"{produced}*.xml")) + list(outdir.glob(f"{produced}*.xml"))
    if not candidates:
        raise MuodeToolError(f"gapseq produced no SBML model for {genome_fna}")
    return candidates[0]


def reconstruct_mag(
    genome_fna: str | Path,
    output_xml: str | Path,
    engine: str = "carveme",
    universe: str = "bacteria",
    gapfill_media: Optional[str] = None,
    solver: Optional[str] = None,
) -> Path:
    """End-to-end single-MAG reconstruction used by the workflow's per-MAG rule."""
    genome_fna, output_xml = Path(genome_fna), Path(output_xml)
    if engine == "carveme":
        proteins = output_xml.with_suffix(".faa")
        call_genes(genome_fna, proteins)
        return carveme(proteins, output_xml, universe=universe,
                       gapfill_media=gapfill_media, solver=solver)
    if engine == "gapseq":
        produced = gapseq(genome_fna, output_xml.parent, media=gapfill_media)
        if produced != output_xml:
            produced.replace(output_xml)
        return output_xml
    raise ValueError(f"unknown reconstruction engine '{engine}' (use 'carveme' or 'gapseq')")
