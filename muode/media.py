"""Bridge a muODE :class:`~muode.diet.Diet` to the *media* other tools consume.

A :class:`Diet` is the engine's boundary condition (concentrations + influx).
Two other consumers need the same information in a different shape, and getting
either wrong is silent rather than loud:

* **CarveMe** gap-fills against a *medium* -- a plain set of compound ids in a
  tab-separated "media db".  Carved **without** one (``carve`` with no ``-g``),
  a model is only guaranteed to grow with *every* exchange open: it may import
  finished CoA, glutathione, NMN or peptidoglycan precursors rather than make
  them, because nothing ever forced it to close those gaps.  Such a model looks
  perfect in QC and then flatlines the moment the engine restricts it to a
  defined diet.  :func:`write_carveme_mediadb` renders the diet as the medium to
  carve against, so the gaps are filled for the medium we actually simulate.

* **cobra** judges growth against ``model.medium``.  :func:`diet_medium` builds
  exactly the uptake bounds the dFBA engine imposes at *t=0* -- Michaelis-Menten
  on the diet's initial concentrations, and **zero for anything the diet does not
  name** -- so a QC growth check can ask the question the simulation will ask.

The two are deliberately built from the same :class:`Diet` object: if the medium
you gap-fill on and the medium you simulate on can drift apart, they will.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from muode.diet import Diet


def supplied_metabolites(diet: Diet) -> Dict[str, float]:
    """The metabolites a diet actually provides, mapped to their t=0 concentration.

    A metabolite listed at concentration 0 with no influx (a fermentation product
    that the community is expected to cross-feed) is *not* supplied: it starts at
    zero, so its uptake bound is zero, and no organism can consume it.
    """
    return {
        m: diet.initial_concentration(m)
        for m in diet.metabolites()
        if diet.initial_concentration(m) > 0.0 or diet.influx_rate(m) > 0.0
    }


def to_compound(metabolite_id: str) -> str:
    """Extracellular metabolite id -> CarveMe media-db compound id (``glc__D_e`` -> ``glc__D``)."""
    return metabolite_id[:-2] if metabolite_id.endswith("_e") else metabolite_id


def write_carveme_mediadb(diet: Diet, path: str | Path, medium: Optional[str] = None) -> Path:
    """Write ``diet`` as a CarveMe media-db TSV and return the path.

    CarveMe reads this with ``--mediadb`` and gap-fills on the named medium with
    ``-g``; it expects one row per compound, keyed by medium, with compound ids in
    the BiGG namespace *without* a compartment suffix.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    name = medium or diet.name
    lines = ["medium\tdescription\tcompound"]
    for met in sorted(supplied_metabolites(diet)):
        lines.append(f"{name}\t{name} (muODE diet)\t{to_compound(met)}")
    path.write_text("\n".join(lines) + "\n")
    return path


def diet_medium(model, diet: Diet, kinetics=None) -> Dict[str, float]:
    """The uptake bounds the dFBA engine imposes on ``model`` at t=0, as a cobra medium.

    Keyed by exchange-reaction id, so it can be assigned to ``model.medium``.
    Exchanges for metabolites the diet does not supply get a bound of 0 -- which is
    the whole point: the engine closes them, so a growth check must too.

    The bound is ``min(Michaelis-Menten rate, dietary availability)``.  The kinetic
    term is what the *cell* can transport; ``diet.max_uptake`` is what the *diet*
    supplies.  Neither alone is sufficient: kinetics without a dietary ceiling let a
    model import every nutrient in a large medium at the same uniform Vmax (which
    is how a GEM ends up "growing" at 6.5/h), and a dietary ceiling without kinetics
    ignores that a cell with no transporter for a sugar cannot eat it however much
    is present.  Diets that specify no ``max_uptake`` are unchanged.
    """
    from muode.kinetics import KineticParameters

    kinetics = kinetics or KineticParameters()
    conc = supplied_metabolites(diet)
    medium: Dict[str, float] = {}
    for ex in model.exchanges:
        if len(ex.metabolites) != 1:
            continue
        met = next(iter(ex.metabolites)).id
        rate = kinetics.michaelis_menten(model.id, met, conc.get(met, 0.0))
        limit = diet.uptake_limit(met)
        medium[ex.id] = rate if limit is None else min(rate, limit)
    return medium


def growth_on_diet(model, diet: Diet, kinetics=None) -> float:
    """Maximum biomass flux of ``model`` on ``diet`` (0.0 if infeasible).

    Non-destructive: the model's own bounds are restored afterwards.
    """
    with model:
        model.medium = diet_medium(model, diet, kinetics)
        value = model.slim_optimize()
    return 0.0 if value is None or value != value else float(value)  # NaN => infeasible
