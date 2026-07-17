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


# ---------------------------------------------------------------------------
# BiGG diet -> ModelSEED (gapseq) diet
# ---------------------------------------------------------------------------
#
# muODE is BiGG-native: diets, kinetics and markers all key on BiGG ids like
# ``glc__D_e``.  A gapseq community is ModelSEED-native (``cpd00027_e0``), and it
# is internally consistent -- its members cross-feed with no id translation at
# all.  The *only* BiGG<->ModelSEED boundary is the diet.  So we translate the
# diet across once, rather than harmonising every model's ~160 exchanges.
#
# The map is built from the models' OWN ``bigg.metabolite`` annotations, not from
# memory -- but that alone is wrong for three kinds of metabolite, so a small
# hand-verified layer sits on top.  Every id in it was read out of a gapseq model
# BY METABOLITE NAME (never guessed):
#
#   * SCFAs the community secretes (butyrate, propionate) carry NO
#     ``bigg.metabolite`` annotation, so the auto-map cannot see them at all --
#     yet they are the acidification arm of colonisation resistance and must be
#     tracked pools.
#   * gapseq names ammonia "NH3" and annotates it ``nh3``; the diet's nitrogen row
#     is ``nh4``, so the auto-map never connects them.
#   * lactose appears as two ModelSEED compounds (``cpd00208`` LACT and
#     ``cpd01354`` beta-Lactose) *both* annotated ``lcts``; the auto-map would
#     flag it ambiguous.  beta-Lactose is the transported form.
_CURATED_BIGG_TO_MODELSEED = {
    "but": "cpd00211",    # Butyrate-e0       -- secreted; no bigg.metabolite annotation
    "ppa": "cpd00141",    # Propionate-e0     -- secreted; no bigg.metabolite annotation
    "nh4": "cpd00013",    # NH3-e0            -- gapseq annotates ammonia as nh3, not nh4
    "lcts": "cpd01354",   # beta-Lactose-e0   -- also cpd00208; this is the transported form
    "ni2": "cpd00244",    # Ni2+-e0           -- trace metal; no bigg.metabolite annotation
    "fuc__L": "cpd00751",  # L-Fucose-e0      -- colonic sugar; no bigg.metabolite annotation
}


# A polymer is not a compound -- it is a compound TIMES a chain length, and the two
# namespaces disagree about the length.  BiGG encodes dietary starch as `starch1200`
# (1200 glucose units); gapseq/ModelSEED encodes it as `cpd90003` ("starch (n=27,
# 3xalpha1-6, 23xalpha1-4)").  Neither id resolves to the other: the diet says
# `starch1200`, the models annotate the bare stem `starch`, and that stem is ambiguous
# anyway (it sits on BOTH cpd90003 and cpd90004), so the auto-map correctly refuses it
# and the ambiguity guard drops the row.  The result was a *fibre-free* colonic diet.
#
# That gap was not cosmetic.  All four donors in examples/fmt_cdiff carry starch
# exchanges and C. difficile carries none, so starch is the one substrate in the medium
# the donor community can eat and the pathogen cannot -- and it was the substrate the
# translation silently threw away.  The competition arm was then read as evidence that
# the donors are poor competitors.
#
# TRANSLATING A POLYMER REQUIRES A UNIT CONVERSION, NOT A RENAME.  Carrying max_uptake
# across verbatim would deliver 27/1200 of the dietary carbon -- a 44x cut -- because
# the flux is molar and the molecules are different sizes.  The invariant that IS
# conserved across the boundary is the **monomer** (glucose-equivalent) flux:
#
#     glucose-equiv/gDW/h = max_uptake(starch1200) * 1200 = 0.0001 * 1200 = 0.12
#     max_uptake(cpd90003) = 0.12 / 27                                   = 0.00444
#
# Both degrees of polymerisation are STATED BY THE SOURCE DATA -- 1200 by the BiGG id,
# 27 by the ModelSEED compound's own name -- so the factor is derived, not chosen.
# Concentration and influx scale by the same ratio, for the same reason.
#
# The one judgement call is cpd90003 (n=27) over cpd90004 (n=19): both are annotated
# `starch` and both are carried by the same four donors.  Because the conversion
# conserves glucose-equivalents, the choice moves the molar flux but not the carbon,
# which is what growth is limited by; supplying both would double the fibre.
#
#: BiGG polymer id -> (ModelSEED compound, monomers per BiGG molecule, monomers per
#: ModelSEED molecule).  Fluxes and concentrations are scaled by ``bigg_n / seed_n``.
_CURATED_POLYMER_BIGG_TO_MODELSEED = {
    "starch1200": ("cpd90003", 1200, 27),
}
# Deliberately NOT in this table, and each for a reason worth keeping:
#   amylose300, pullulan1200, lmn30  -- no corresponding compound in ANY gapseq model,
#       so there is nothing to translate to.  Inventing one would invent a nutrient.
#   xylan4, xylan8  -- cpd90021/cpd90022 exist, but C. difficile carries them TOO, so
#       xylan feeds the pathogen alongside the donors.  It is a real gap in the medium
#       and translating it is defensible future work; it is not a competition asymmetry,
#       and it must not be added quietly under cover of the starch fix.


def build_bigg_to_modelseed(
    models, curated: Optional[Dict[str, str]] = None
) -> tuple[Dict[str, str], Dict[str, list]]:
    """Build a ``{bigg_base: modelseed_cpd}`` map from gapseq models' exchange annotations.

    Reads every single-metabolite exchange in ``models`` and, where its metabolite
    carries a ``bigg.metabolite`` annotation, records ``bigg -> cpd`` (the compound
    id with the compartment stripped, e.g. ``cpd00027_e0`` -> ``cpd00027``).  The
    hand-verified ``curated`` pins win over -- and silence the ambiguity of -- the
    annotation-derived entries.

    Returns ``(mapping, report)``.  ``report["ambiguous"]`` lists every BiGG id that
    the models annotate onto more than one compound and that no curated pin
    resolves: those are deliberately left OUT of ``mapping`` rather than guessed.
    """
    curated = dict(_CURATED_BIGG_TO_MODELSEED if curated is None else curated)
    candidates: Dict[str, set] = {}
    for model in models:
        for ex in model.exchanges:
            if len(ex.metabolites) != 1:
                continue
            met = next(iter(ex.metabolites))
            bigg = (met.annotation or {}).get("bigg.metabolite")
            # gapseq annotates some compounds with SEVERAL BiGG synonyms, e.g.
            # Zn2+ is ["HC02172", "zn2"].  Record the compound under EVERY synonym,
            # not just the first -- else the diet's own id (zn2) silently misses,
            # and an essential trace metal drops out of the medium.
            biggs = bigg if isinstance(bigg, list) else [bigg]
            cpd = met.id.rsplit("_", 1)[0]
            for b in biggs:
                if b:
                    candidates.setdefault(b, set()).add(cpd)

    mapping: Dict[str, str] = {}
    ambiguous: Dict[str, list] = {}
    for bigg, cpds in candidates.items():
        if bigg in curated:
            continue                       # curated wins, added below
        if len(cpds) == 1:
            mapping[bigg] = next(iter(cpds))
        else:
            ambiguous[bigg] = sorted(cpds)  # do not pick; report it
    mapping.update(curated)
    return mapping, {"ambiguous": ambiguous}


def translate_diet_to_modelseed(
    diet: Diet,
    mapping: Dict[str, str],
    compartment: str = "_e0",
    polymers: Optional[Dict[str, tuple]] = None,
) -> tuple[Diet, list]:
    """Translate a BiGG-namespace ``diet`` into ModelSEED ids a gapseq community eats.

    Every metabolite's concentration, influx, ``max_uptake`` and provenance are
    carried over verbatim onto its ModelSEED id (``glc__D_e`` -> ``cpd00027_e0``);
    only the id changes.  A row whose BiGG base has no entry in ``mapping`` is
    dropped and returned in ``unresolved`` -- so the loss is reported, never silent.
    (A diet that does not name a metabolite gives it uptake bound 0, so a dropped
    *fermentation product* still starts at 0 exactly as it did; a dropped *nutrient*
    is a real gap in the medium, which is why the caller is handed the list.)

    ``polymers`` handles the rows a rename cannot: the two namespaces encode a polymer
    at *different chain lengths*, so its row is rescaled by the ratio of the degrees of
    polymerisation to conserve **monomer** flux (see
    ``_CURATED_POLYMER_BIGG_TO_MODELSEED``).  Carrying a polymer's molar flux across
    verbatim would silently change how much carbon the diet delivers.

    Returns ``(translated_diet, unresolved)``.
    """
    polymers = _CURATED_POLYMER_BIGG_TO_MODELSEED if polymers is None else polymers
    conc: Dict[str, float] = {}
    influx: Dict[str, float] = {}
    max_uptake: Dict[str, float] = {}
    source: Dict[str, str] = {}
    unresolved: list = []
    for mid in diet.metabolites():
        base = to_compound(mid)
        # A polymer's own pin wins over the annotation map AND over the ambiguity
        # guard, because it carries the chain-length conversion the map cannot express.
        if base in polymers:
            cpd, bigg_n, seed_n = polymers[base]
            scale = bigg_n / seed_n
        else:
            cpd, scale = mapping.get(base), 1.0
        if cpd is None:
            unresolved.append(base)
            continue
        new = f"{cpd}{compartment}"
        conc[new] = diet.initial_concentration(mid) * scale
        if diet.influx_rate(mid):
            influx[new] = diet.influx_rate(mid) * scale
        limit = diet.uptake_limit(mid)
        if limit is not None:
            max_uptake[new] = limit * scale
        prov = diet.provenance(mid)
        if scale != 1.0:
            # This row's numbers are no longer the source's numbers -- they were divided
            # by a chain-length ratio.  The `source` column is the mechanism that exists
            # for saying so, and a rescaled row that still claims `intake` is a lie.
            prov = f"{prov or 'unknown'}+dp_scaled"
        if prov:
            source[new] = prov
    return (
        Diet(conc, influx, max_uptake, source, name=f"{diet.name}_modelseed"),
        sorted(unresolved),
    )


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
