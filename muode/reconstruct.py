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


def call_genes_eukaryote(
    genome_fna: str | Path,
    proteins_faa: str | Path,
    ref_db: str | Path,
    threads: int = 1,
) -> Path:
    """Predict proteins from a *eukaryotic* (MAG) nucleotide FASTA with MetaEuk.

    prodigal (used for prokaryotes) cannot call eukaryotic genes — introns and
    splicing break its ORF model.  MetaEuk performs reference-based eukaryotic
    gene prediction directly on contigs, which is the right fit for eukaryotic
    genome bins.  ``ref_db`` is a MetaEuk/MMseqs2 protein reference database
    (e.g. UniRef90 or a curated eukaryotic proteome set).

    Returns the predicted-protein FASTA (``proteins_faa``).
    """
    genome_fna, proteins_faa = Path(genome_fna), Path(proteins_faa)
    proteins_faa.parent.mkdir(parents=True, exist_ok=True)
    metaeuk = require(
        "metaeuk", env_hint="mamba install -c bioconda metaeuk"
    )
    prefix = proteins_faa.with_suffix("")          # MetaEuk writes <prefix>.fas
    tmpdir = proteins_faa.parent / f"{proteins_faa.stem}.metaeuk_tmp"
    run(
        [metaeuk, "easy-predict", str(genome_fna), str(ref_db),
         str(prefix), str(tmpdir), "--threads", str(threads)],
        log=proteins_faa.with_suffix(".metaeuk.log"),
    )
    produced = Path(f"{prefix}.fas")
    if not produced.exists():
        raise MuodeToolError(
            f"MetaEuk produced no protein FASTA for {genome_fna} (expected {produced})"
        )
    if produced != proteins_faa:
        produced.replace(proteins_faa)
    return proteins_faa


def carvefungi(
    proteins_faa: str | Path,
    output_xml: str | Path,
    threads: int = 1,
    cmd_template: Optional[str] = None,
) -> Path:
    """Reconstruct a *fungal* GEM from a proteome with CarveFungi.

    CarveFungi (Castillo-Priego et al.) is the eukaryote-fungal analog of
    CarveMe: a deep-learning model predicts reaction presence from protein
    sequences and assembles a gap-filled fungal GEM.  It is the automated route
    for fungi, where CarveMe (prokaryote) and gapseq (crashes on eukaryotic
    contigs) cannot be used.

    The exact CarveFungi entry point varies by install, so the invocation is
    overridable via ``cmd_template`` (a format string with ``{proteins}`` and
    ``{output}`` placeholders, e.g. from the ``carvefungi_cmd`` config key).
    The default assumes a ``carvefungi`` executable on PATH.
    """
    proteins_faa, output_xml = Path(proteins_faa), Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    if cmd_template:
        cmd = cmd_template.format(
            proteins=str(proteins_faa), output=str(output_xml), threads=threads
        ).split()
    else:
        cf = require(
            "carvefungi",
            env_hint="see https://github.com/SandraCastilloPriego/CarveFungi",
        )
        cmd = [cf, str(proteins_faa), "-o", str(output_xml), "--threads", str(threads)]
    run(cmd, log=output_xml.with_suffix(".carvefungi.log"))
    if not output_xml.exists():
        # some builds name the output after the proteome; recover it if so.
        alt = list(output_xml.parent.glob(f"{proteins_faa.stem}*.xml"))
        if not alt:
            raise MuodeToolError(f"CarveFungi produced no SBML model for {proteins_faa}")
        alt[0].replace(output_xml)
    return output_xml


def eukaryote_generic(
    proteins_faa: str | Path,
    output_xml: str | Path,
    eggnog_data_dir: Optional[str | Path] = None,
    template: str = "fungi",
    threads: int = 1,
) -> Path:
    """Draft a GEM for a *non-fungal* eukaryote from a proteome (experimental).

    This is the lower-quality fallback for eukaryotes that CarveFungi does not
    cover (protists, etc.): annotate the proteome with orthology (eggNOG-mapper)
    and build a draft from those annotations with ModelSEEDpy.  The result needs
    substantially more gap-filling and curation than a CarveMe/CarveFungi model
    — prefer a curated model (``eukaryote_models``) when one is available.
    """
    proteins_faa, output_xml = Path(proteins_faa), Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)

    # 1) orthology annotation with eggNOG-mapper
    emapper = require(
        "emapper.py", env_hint="mamba install -c bioconda eggnog-mapper"
    )
    ann_prefix = output_xml.with_suffix("")
    cmd = [emapper, "-i", str(proteins_faa), "-o", ann_prefix.name,
           "--output_dir", str(output_xml.parent), "--itype", "proteins",
           "--cpu", str(threads)]
    if eggnog_data_dir:
        cmd += ["--data_dir", str(eggnog_data_dir)]
    run(cmd, log=output_xml.with_suffix(".emapper.log"))
    annotations = Path(f"{ann_prefix}.emapper.annotations")
    if not annotations.exists():
        raise MuodeToolError(
            f"eggNOG-mapper produced no annotations for {proteins_faa} "
            f"(expected {annotations})"
        )

    # 2) build a draft from the annotations with ModelSEEDpy
    try:
        from muode.eukaryote_draft import build_from_eggnog
    except ImportError as exc:  # pragma: no cover - exercised on the workstation
        raise MuodeToolError(
            "the generic eukaryote engine needs ModelSEEDpy "
            "(`pip install modelseedpy`); or supply a curated model via "
            "`eukaryote_models`"
        ) from exc
    build_from_eggnog(annotations, proteins_faa, output_xml, template=template)
    if not output_xml.exists():
        raise MuodeToolError(f"generic eukaryote build produced no model for {proteins_faa}")
    return output_xml


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
    ref_db: Optional[str | Path] = None,
    threads: int = 1,
    carvefungi_cmd: Optional[str] = None,
    eggnog_data_dir: Optional[str | Path] = None,
) -> Path:
    """End-to-end single-MAG reconstruction used by the workflow's per-MAG rule.

    Engines:
      * ``carveme``           — prokaryotes (prodigal genes + CarveMe, BiGG).
      * ``gapseq``            — prokaryotes (ModelSEED namespace).
      * ``carvefungi``        — fungi (MetaEuk genes + CarveFungi).
      * ``eukaryote_generic`` — other eukaryotes (MetaEuk genes + eggNOG/ModelSEEDpy).
      * ``stub``              — dependency-free placeholder (local DAG testing).

    The eukaryote engines call genes with MetaEuk, which needs a protein
    reference database (``ref_db``).
    """
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
    if engine in ("carvefungi", "eukaryote_generic"):
        if not ref_db:
            raise MuodeToolError(
                f"engine '{engine}' needs a MetaEuk protein reference database; "
                "set `euk_ref_db` in the config (or supply a curated model via "
                "`eukaryote_models`)"
            )
        proteins = output_xml.with_suffix(".faa")
        call_genes_eukaryote(genome_fna, proteins, ref_db, threads=threads)
        if engine == "carvefungi":
            return carvefungi(proteins, output_xml, threads=threads,
                              cmd_template=carvefungi_cmd)
        return eukaryote_generic(proteins, output_xml,
                                 eggnog_data_dir=eggnog_data_dir, threads=threads)
    if engine == "stub":
        return stub_reconstruct(genome_fna, output_xml)
    raise ValueError(
        f"unknown reconstruction engine '{engine}' (use 'carveme', 'gapseq', "
        "'carvefungi', 'eukaryote_generic' or 'stub')"
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
