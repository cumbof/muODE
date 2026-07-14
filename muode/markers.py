"""Marker genes for the traits a genome-scale model **cannot** encode.

Why this module exists
----------------------
muODE has two halves that do not currently touch.  Reconstruction turns genomes
into GEMs; the ecology layers (:mod:`muode.bile`, :mod:`muode.lifecycle`,
:mod:`muode.antagonism`, :mod:`muode.antibiotic`, :mod:`muode.phage`) turn
*guild membership* into dynamics.  Nothing turns genomes into guild membership.
So every ecology layer is hand-populated in the designed scenarios
(``muode.scenarios.cdi`` sets ``bai=True`` on a species by hand) and is silently
**empty** on a real MAG run.  This module is the missing bridge.

The rule for what belongs here
------------------------------
**Can the GEM already represent it?**  If yes, this is the wrong module.

Urease is the worked example.  It looked like a candidate -- host urea, a real
gut function, present in only some organisms.  But CarveMe already gives 23 of the
89 SRR13844389 MAGs an ``EX_urea_e``, because urea hydrolysis is ordinary
metabolism: it consumes a substrate, releases ammonia, and the FBA stoichiometry
handles it.  Urease needed a *medium* fix (supply urea), not a marker gene.

What lands here is what FBA structurally cannot carry:

* **gratuitous transformations** -- bile-acid deconjugation and 7a-dehydroxylation
  change the shared pool but are not growth-coupled, so a biomass objective will
  never route flux through them, no matter how good the reconstruction;
* **dormancy** -- sporulation is not a metabolic state;
* **warfare** -- a bacteriocin costs its producer and kills its target;
* **host interaction** and **inhibition** -- effects on *other* organisms' growth.

Every :class:`FunctionalTrait` must therefore say, in ``why_not_fba``, why the
reconstruction cannot supply it.  If that field is hard to write, the trait almost
certainly belongs in the medium or in the reconstruction instead.

Evidence, and how much to trust it
----------------------------------
Two kinds of marker, and they are **not** equally trustworthy:

* :attr:`MarkerSource.PFAM` -- an HMM with a curated gathering threshold.  Search
  with ``hmmsearch --cut_ga`` and the cutoff is Pfam's, not ours.  No invented
  numbers.
* :attr:`MarkerSource.REFPROT` -- homology to a named reference protein.  The
  bitscore cutoff *is* a judgement call, so we lean on **multi-marker evidence**
  instead: a trait needs ``min_markers`` distinct proteins before it is called.

That distinction is load-bearing for the *bai* operon.  **Pfam cannot call bai.**
BaiCD and BaiH are Old Yellow Enzyme family; BaiA is a short-chain dehydrogenase;
BaiB and BaiF are CoA transferases.  Those families have thousands of members
doing unrelated chemistry, so a Pfam hit is worthless as evidence of
7a-dehydroxylation.  bai is called from Swiss-Prot reference proteins of
*Clostridium scindens*, requiring several of them together.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Mapping, Tuple


class MarkerSource(str, Enum):
    """Where a marker's evidence comes from -- and so how far to trust one hit."""

    PFAM = "pfam"        # HMM + curated gathering threshold (hmmsearch --cut_ga)
    REFPROT = "refprot"  # homology to a named reference protein (diamond/blastp)


@dataclass(frozen=True)
class Marker:
    accession: str
    gene: str
    source: MarkerSource


@dataclass(frozen=True)
class FunctionalTrait:
    """A capability the GEM cannot encode, called from marker genes.

    Parameters
    ----------
    name:
        Trait id, and the column name in ``traits.tsv``.
    markers:
        Marker genes.  A trait is called when ``min_markers`` distinct ones hit.
    min_markers:
        Evidence bar.  For a specific, curated Pfam family, 1 is enough.  For
        reference proteins drawn from promiscuous enzyme families, it must be
        several -- that requirement, not the bitscore, is what makes the call.
    layer:
        The ecology layer that consumes this trait.  A trait no layer reads is
        dead weight; say which one, or do not add it.
    why_not_fba:
        Why reconstruction cannot supply this.  If you cannot fill this in, the
        trait does not belong in this module.
    citation:
        Where the marker set comes from.
    """

    name: str
    markers: Tuple[Marker, ...]
    min_markers: int
    layer: str
    why_not_fba: str
    citation: str

    def call(self, hits: frozenset) -> bool:
        """Is this trait present, given the accessions that hit in one genome?"""
        return sum(1 for m in self.markers if m.accession in hits) >= self.min_markers


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------

#: Bile-salt hydrolase.  One curated Pfam family, specific to the chemistry:
#: choloylglycine hydrolase (CBAH).  Deconjugation is the near-universal first
#: step -- Heinken et al. find BSH in 204 of 693 gut genomes (29%) -- so expect a
#: substantial fraction of any gut community to be BSH+.
BSH = FunctionalTrait(
    name="bsh",
    markers=(Marker("PF02275", "cbah", MarkerSource.PFAM),),
    min_markers=1,
    layer="muode.bile.BileAcidTransform (bsh_producers)",
    why_not_fba=(
        "Deconjugation transforms the shared bile pool without conserving energy "
        "or yielding biomass, so an FBA biomass objective routes zero flux through "
        "it.  It also has no BiGG reaction: CarveMe's universe carries no bile-acid "
        "exchange at all (0/89 MAGs), so reconstruction cannot represent it however "
        "well the genome is annotated."
    ),
    citation=(
        "Pfam PF02275 (CBAH, linear amide C-N hydrolase / choloylglycine hydrolase). "
        "Prevalence: Heinken et al., Microbiome 7:75 (2019)."
    ),
)

#: The *bai* operon -- 7a-dehydroxylation, primary -> secondary bile acids.
#:
#: RARE.  Heinken et al. find the operon in 7 of 693 gut genomes (~1%), confined to
#: *Clostridioides*, *Lachnoclostridium* and *Eggerthella*.  A community with no
#: bai carrier CANNOT make deoxycholate or lithocholate, and that is a result to
#: report, not a gap to fill.
#:
#: Swiss-Prot reference proteins of *Clostridium scindens* (the type organism for
#: this pathway).  Pfam is useless here -- see the module docstring.  We require 3
#: of the 5, because BaiA2 is a short-chain dehydrogenase and would produce
#: spurious lone hits against any SDR-rich genome; three bile-acid-specific
#: proteins together will not.
BAI = FunctionalTrait(
    name="bai",
    markers=(
        Marker("P19409", "baiB", MarkerSource.REFPROT),   # bile acid--CoA ligase
        Marker("P19410", "baiCD", MarkerSource.REFPROT),  # 3-oxocholoyl-CoA 4-desaturase
        Marker("P19412", "baiE", MarkerSource.REFPROT),   # bile acid 7alpha-dehydratase
        Marker("P19337", "baiA2", MarkerSource.REFPROT),  # 3a-hydroxy BA-CoA 3-dehydrogenase
        Marker("P32369", "baiG", MarkerSource.REFPROT),   # bile acid transporter
    ),
    min_markers=3,
    layer="muode.bile.BileAcidTransform (bai_producers)",
    why_not_fba=(
        "Same as BSH -- no BiGG reaction exists, and the pathway is an electron sink "
        "rather than a growth-coupled one, so the biomass objective would not use it "
        "even if it did."
    ),
    citation=(
        "UniProt/Swiss-Prot, Clostridium scindens bai operon: P19409 (BaiB), "
        "P19410 (BaiCD), P19412 (BaiE), P19337 (BaiA2), P32369 (BaiG). "
        "Prevalence: Heinken et al., Microbiome 7:75 (2019); Ridlon et al."
    ),
)

#: Defence against reactive oxygen species.  Catalase, Fe/Mn superoxide dismutase
#: (two Pfam domains, one enzyme), cytochrome bd oxidase.
#:
#: HONEST LIMIT: this separates "has ROS defence" from "has none".  It does NOT
#: distinguish aerotolerant from facultative from microaerophilic -- gene presence
#: cannot, and pretending otherwise would put a number on a guess.  Use it to flag
#: organisms whose OBLIGATE_ANAEROBE default looks wrong, not to assign a class.
OXYGEN_DEFENCE = FunctionalTrait(
    name="oxygen_defence",
    markers=(
        Marker("PF00199", "katA", MarkerSource.PFAM),   # catalase
        Marker("PF00081", "sodA_N", MarkerSource.PFAM), # Fe/Mn SOD, alpha-hairpin
        Marker("PF02777", "sodA_C", MarkerSource.PFAM), # Fe/Mn SOD, C-terminal
        Marker("PF01654", "cydA", MarkerSource.PFAM),   # cytochrome bd oxidase I
    ),
    min_markers=1,
    layer="muode.traits.OxygenTolerance / muode.oxygen",
    why_not_fba=(
        "Oxygen tolerance is a survival property, not a flux.  A GEM carrying no "
        "oxygen-consuming reaction is indistinguishable from one that is killed by "
        "oxygen; only the ROS-defence genes tell them apart."
    ),
    citation="Pfam PF00199 (catalase), PF00081/PF02777 (Fe/Mn SOD), PF01654 (cyt bd).",
)

#: Every trait muODE can call.  Add here; the caller and the TSV follow automatically.
TRAITS: Tuple[FunctionalTrait, ...] = (BSH, BAI, OXYGEN_DEFENCE)


def by_name() -> Dict[str, FunctionalTrait]:
    return {t.name: t for t in TRAITS}


def accessions(source: MarkerSource) -> Tuple[str, ...]:
    """Every accession to search for with one tool, deduplicated.

    The search step needs this: all PFAM accessions go to one ``hmmsearch``, all
    REFPROT accessions to one ``diamond blastp``, rather than one search per trait.
    """
    seen: Dict[str, None] = {}
    for trait in TRAITS:
        for m in trait.markers:
            if m.source is source:
                seen[m.accession] = None
    return tuple(seen)


def call_all(hits: Mapping[str, frozenset]) -> Dict[str, Dict[str, bool]]:
    """``{mag: {accession, ...}}`` -> ``{mag: {trait: present}}``.

    Pure, so the trait logic is testable without HMMER or DIAMOND on the box.
    """
    return {
        mag: {t.name: t.call(acc) for t in TRAITS}
        for mag, acc in hits.items()
    }
