"""Phase 2 -- gap filling.

Two complementary mechanisms, mirroring the README:

* **Algorithmic LP gap filling** (implemented here with cobra): find the minimal
  set of reactions from a universal model that lets a draft GEM produce biomass
  on a defined medium.  This guarantees every model is simulatable.

* **AI-driven gap filling** (:func:`metapathpredict`, optional wrapper): use
  MetaPathPredict to flag KEGG metabolic modules that are *probably present* but
  unobserved because the MAG is incomplete.  Those predictions bias / seed the LP
  gap filling toward biologically plausible reactions instead of arbitrary ones.

A practical caveat documented in :doc:`/docs/EVALUATION`: MetaPathPredict works
in the KEGG-module namespace, whereas CarveMe models are in BiGG.  The predicted
modules must be mapped (via KEGG<->BiGG/ModelSEED cross-references) before they
can constrain the LP -- this mapping is the real integration cost, so AI gap
filling is an opt-in refinement rather than a default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

from muode.external import require, run


def grows(model, threshold: float = 1e-6) -> bool:
    """Whether the model carries flux through its biomass objective."""
    return (model.slim_optimize() or 0.0) > threshold


def lp_gapfill(
    model,
    universal,
    threshold: float = 1e-3,
    integer: bool = False,
) -> List[str]:
    """Minimal-reaction LP gap fill so ``model`` can produce biomass.

    Returns the ids of reactions added from ``universal``.  Uses cobra's
    gapfilling solver; raises if no gap-filling solution exists on the current
    medium (a strong signal the medium or draft is wrong).
    """
    from cobra.flux_analysis import gapfill as cobra_gapfill

    if grows(model, threshold):
        return []
    solutions = cobra_gapfill(
        model, universal, demand_reactions=False, integer_threshold=1e-9,
        lower_bound=threshold, iterations=1,
    )
    added = []
    for rxn in solutions[0]:
        if rxn.id not in model.reactions:
            model.add_reactions([rxn.copy()])
        added.append(rxn.id)
    return added


def ensure_biomass(model, universal=None, threshold: float = 1e-3) -> dict:
    """Make sure a model grows; gap-fill against ``universal`` if provided."""
    before = bool(grows(model, threshold))
    added: List[str] = []
    if not before and universal is not None:
        added = lp_gapfill(model, universal, threshold=threshold)
    return {
        "grew_initially": before,
        "reactions_added": added,
        "grows_now": bool(grows(model, threshold)),
    }


def metapathpredict(
    annotations: str | Path,
    output_tsv: str | Path,
) -> Path:
    """Run MetaPathPredict to predict KEGG module presence from a MAG.

    ``annotations`` is a KEGG-ortholog annotation table (e.g. from KofamScan /
    DRAM).  This is an optional wrapper; its predictions feed
    :func:`lp_gapfill` after namespace mapping.
    """
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    mpp = require("MetaPathPredict", env_hint="pip install MetaPathPredict")
    run([mpp, "predict", "-i", str(annotations), "-o", str(output_tsv)],
        log=output_tsv.with_suffix(".log"))
    return output_tsv
