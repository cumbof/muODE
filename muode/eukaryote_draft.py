"""Generic (non-fungal) eukaryote draft reconstruction from eggNOG annotations.

This is the **experimental, lower-quality** fallback used by
:func:`muode.reconstruct.eukaryote_generic` for eukaryotes that CarveFungi does
not cover (protists, etc.).  CarveMe is prokaryote-only and gapseq crashes on
eukaryotic contigs, so there is no fast push-button tool here; the best automated
option is orthology annotation (eggNOG-mapper) followed by a template-based draft
with ModelSEEDpy.  The resulting model needs substantially more gap-filling and
curation than a CarveMe/CarveFungi model — prefer a curated model
(``eukaryote_models``) whenever one is available.

The eggNOG parsing is self-contained and tested; the ModelSEEDpy build is kept in
one small function so it is easy to adapt to a given ModelSEEDpy version.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from muode.external import MuodeToolError


def parse_eggnog_kos(annotations: str | Path) -> Dict[str, List[str]]:
    """Parse an ``*.emapper.annotations`` file into ``{query: [KO, ...]}``.

    eggNOG-mapper writes a tab-separated table whose comment lines start with
    ``#``; the first non-comment-prefixed line beginning with ``#query`` is the
    header.  The ``KEGG_ko`` column holds values like ``ko:K00001,ko:K00002``
    (or ``-`` when none).  KO ids are returned without the ``ko:`` prefix.
    """
    annotations = Path(annotations)
    ko_map: Dict[str, List[str]] = {}
    header: List[str] = []
    with annotations.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                # the header row is the comment line starting with "#query"
                if line.lstrip("#").startswith("query"):
                    header = line.lstrip("#").split("\t")
                continue
            if not header:
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            query = row.get("query") or fields[0]
            raw = row.get("KEGG_ko", "-")
            if raw in ("", "-"):
                continue
            kos = [k.split(":", 1)[-1] for k in raw.split(",") if k and k != "-"]
            if kos:
                ko_map.setdefault(query, []).extend(kos)
    return ko_map


def build_from_eggnog(
    annotations: str | Path,
    proteins_faa: str | Path,
    output_xml: str | Path,
    template: str = "fungi",
) -> Path:
    """Build a draft GEM from eggNOG annotations with ModelSEEDpy.

    Writes a ``<output>.gene_ko.tsv`` companion (gene -> KO list) alongside the
    model so the orthology evidence behind the draft is inspectable.  Raises
    :class:`MuodeToolError` with actionable guidance if ModelSEEDpy is missing or
    its build API differs on the installed version.
    """
    annotations, proteins_faa = Path(annotations), Path(proteins_faa)
    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)

    ko_map = parse_eggnog_kos(annotations)
    if not ko_map:
        raise MuodeToolError(
            f"no KEGG KO annotations found in {annotations}; cannot draft a model"
        )
    # Persist the orthology evidence regardless of the build outcome.
    ko_tsv = output_xml.with_suffix(".gene_ko.tsv")
    with ko_tsv.open("w") as out:
        out.write("gene\tKO\n")
        for gene, kos in sorted(ko_map.items()):
            out.write(f"{gene}\t{','.join(kos)}\n")

    try:
        from modelseedpy import MSBuilder, MSGenome  # type: ignore
    except ImportError as exc:
        raise MuodeToolError(
            "the generic eukaryote engine needs ModelSEEDpy "
            "(`pip install modelseedpy`); or supply a curated model via "
            "`eukaryote_models`"
        ) from exc

    try:
        genome = MSGenome.from_fasta(str(proteins_faa), split=" ")
        # Attach the eggNOG KO evidence as functional annotations.
        for feature in genome.features:
            for ko in ko_map.get(feature.id, []):
                feature.add_ontology_term("KEGG", ko)
        model = MSBuilder(genome, template).build(
            str(output_xml.stem), index="0", allow_all_non_grp_reactions=True,
        )
        import cobra.io

        cobra.io.write_sbml_model(model, str(output_xml))
    except MuodeToolError:
        raise
    except Exception as exc:  # pragma: no cover - depends on ModelSEEDpy version
        raise MuodeToolError(
            f"ModelSEEDpy draft build failed for {proteins_faa}: {exc}. "
            "The generic eukaryote engine is experimental; consider supplying a "
            "curated model via `eukaryote_models`. The gene->KO evidence was "
            f"written to {ko_tsv}."
        ) from exc
    return output_xml
