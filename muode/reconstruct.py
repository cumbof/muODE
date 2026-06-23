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
    if engine == "stub":
        return stub_reconstruct(genome_fna, output_xml)
    raise ValueError(
        f"unknown reconstruction engine '{engine}' (use 'carveme', 'gapseq' or 'stub')"
    )


# ---------------------------------------------------------------------------
# Stub engine -- a dependency-free placeholder so the *workflow* can be run and
# validated end-to-end on machines where CarveMe/gapseq cannot (e.g. aarch64).
# ---------------------------------------------------------------------------

#: roles a stub MAG can declare via a ``muode-stub:<role>`` tag in its first
#: FASTA header, each mapping to a toy metabolism in the ``*_e`` namespace of the
#: built-in diet presets so a set of stubs assembles into a cross-feeding community.
STUB_ROLES = ("fermenter", "consumer", "generic")


def _stub_role(genome_fna: str | Path) -> str:
    """Read the role tag from a stub MAG's first FASTA header (default 'generic')."""
    with Path(genome_fna).open() as fh:
        for line in fh:
            if line.startswith(">"):
                low = line.lower()
                for role in STUB_ROLES:
                    if f"muode-stub:{role}" in low:
                        return role
                break
    return "generic"


def stub_reconstruct(genome_fna: str | Path, output_xml: str | Path) -> Path:
    """Write a tiny, *simulatable* placeholder GEM without any external tool.

    This is **not** a scientific reconstruction -- it exists so the Snakemake
    pipeline (reconstruct -> refine -> QC -> assemble -> simulate) can be exercised
    locally before moving to a machine where CarveMe runs.  The model's metabolism
    is selected by a ``muode-stub:<role>`` tag in the genome's first FASTA header:

    * ``fermenter`` -- takes up glucose, secretes acetate, grows;
    * ``consumer``  -- grows only on acetate (cross-feeds off a fermenter);
    * ``generic``   -- grows on glucose alone (the default).

    Exchange-metabolite ids (``glc_e``, ``ac_e``) match the built-in diet presets,
    so a directory of stub MAGs assembles into a real cross-feeding community.
    Swap ``engine: carveme`` on a capable machine for actual models.
    """
    from cobra import Metabolite, Model, Reaction

    genome_fna, output_xml = Path(genome_fna), Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    role = _stub_role(genome_fna)

    model = Model(f"stub_{genome_fna.stem}")
    glc_e = Metabolite("glc_e", name="D-glucose [e]", compartment="e")
    ac_e = Metabolite("ac_e", name="acetate [e]", compartment="e")

    def _exchange(met: Metabolite, lower: float, upper: float) -> Reaction:
        rxn = Reaction(f"EX_{met.id}", name=f"{met.id} exchange")
        rxn.add_metabolites({met: -1.0})
        rxn.bounds = (lower, upper)
        return rxn

    biomass = Reaction("BIOMASS_stub", name="stub biomass")
    rxns = [biomass]
    if role == "fermenter":
        rxns += [_exchange(glc_e, -10.0, 1000.0), _exchange(ac_e, 0.0, 1000.0)]
        biomass.add_metabolites({glc_e: -1.0, ac_e: 1.0})  # glucose -> acetate + growth
    elif role == "consumer":
        rxns += [_exchange(ac_e, -10.0, 1000.0)]
        biomass.add_metabolites({ac_e: -1.0})              # grows only on acetate
    else:  # generic
        rxns += [_exchange(glc_e, -10.0, 1000.0)]
        biomass.add_metabolites({glc_e: -1.0})
    biomass.bounds = (0.0, 1000.0)

    model.add_reactions(rxns)
    model.objective = "BIOMASS_stub"

    import cobra.io

    cobra.io.write_sbml_model(model, str(output_xml))
    return output_xml
