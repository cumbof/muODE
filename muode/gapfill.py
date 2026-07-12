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


def grows(model, threshold: float = 1e-6, diet=None) -> bool:
    """Whether the model carries flux through its biomass objective.

    ``diet`` restricts uptake to the simulation medium first.  Without it the
    check runs on the model's *stored* medium, which for a CarveMe model carved
    without ``-g`` is wide open -- a test that essentially cannot fail, and so
    tells you nothing about whether the model will grow in the community.
    """
    if diet is not None:
        from muode.media import growth_on_diet

        return growth_on_diet(model, diet) > threshold
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


def ensure_biomass(model, universal=None, threshold: float = 1e-3, diet=None) -> dict:
    """Make sure a model grows; gap-fill against ``universal`` if provided.

    When ``diet`` is given, the LP gap fill (if a universal model is supplied)
    closes the gaps *for the medium the engine will impose*, and growth on that
    medium is reported next to growth with every exchange open.  A model that
    grows only on the latter is not a healthy model: it is an un-gap-filled one
    leaning on nutrients the medium will never provide.
    """
    if diet is not None:
        from muode.media import diet_medium, growth_on_diet

        open_before = float(model.slim_optimize() or 0.0)
        added: List[str] = []
        if growth_on_diet(model, diet) <= threshold and universal is not None:
            # Gap-fill *under the diet's uptake bounds*, so the reactions added are
            # the ones that make the model viable on the medium we simulate.  This
            # cannot run inside a `with model:` block -- cobra's context manager
            # would revert the added reactions on exit -- so the exchange bounds are
            # saved and restored by hand, leaving the written SBML unconstrained.
            saved = {ex.id: ex.bounds for ex in model.exchanges}
            try:
                model.medium = diet_medium(model, diet)
                added = lp_gapfill(model, universal, threshold=threshold)
            finally:
                for rid, bounds in saved.items():
                    model.reactions.get_by_id(rid).bounds = bounds
        # `grows_now` keeps its structural meaning -- "this model can make biomass
        # when given nutrients" -- because it is what gates a model into the
        # community.  Growth on the *diet* must NOT gate: an obligate cross-feeder
        # (an acetate-eating methanogen, or the stub `consumer`) cannot grow on the
        # medium alone by definition -- the community feeds it.  Filtering on diet
        # growth would silently delete exactly the species the simulation exists to
        # study.  So it is reported, and `reconstruct_report` shouts when hardly
        # anything grows on the diet, which is the un-gap-filled signature.
        return {
            "grew_initially": bool(open_before > threshold),
            "reactions_added": added,
            "grows_now": bool(grows(model, threshold)),
            "growth_on_diet": float(growth_on_diet(model, diet)),
            "growth_on_complete_medium": open_before,
            "diet": getattr(diet, "name", None),
        }

    before = bool(grows(model, threshold))
    added = []
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
