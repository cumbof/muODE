"""Phase 3 -- kinetic-parameter prediction (refinement layer).

Two roles, made precise in :mod:`muode.kinetics`:

* **Km** sharpens the Michaelis-Menten *uptake* bounds in the dynamic loop;
* **kcat** caps *intracellular* reaction velocity via enzyme constraints
  (:mod:`muode.enzyme`).

Predictors implement a single :meth:`refine` method that fills a
:class:`~muode.kinetics.KineticParameters` store for one organism's model:

* :class:`HeuristicPredictor` -- the **default**.  Dependency-free and
  deterministic: literature-style Km for common substrates and a turnover number
  jittered around the genome-wide median kcat (~13.7 /s; Bar-Even et al., 2011).
  It needs no sequences, so the whole pipeline runs everywhere -- but the numbers
  are placeholders, not measurements.
* :class:`DLKcatPredictor` / :class:`KmPredictor` -- opt-in wrappers around the
  deep-learning predictors (DLKcat for kcat; Kroll et al. for Km).  They require
  the ``ml`` extra **and** an enzyme/substrate context (per-reaction protein
  sequences + substrate SMILES), which is the real integration cost: CarveMe
  models are in the BiGG namespace, so metabolites must be mapped to SMILES and
  reactions to their catalysing gene sequences before a sequence-based model can
  score them (see :func:`build_enzyme_context`).
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional

from muode.external import MuodeToolError
from muode.kinetics import KineticParameters, KineticPredictor

# Genome-wide median turnover number (Bar-Even et al., Biochemistry 2011), 1/s.
MEDIAN_KCAT = 13.7

# Rough literature Km values (mmol/L) for common exchange metabolites, keyed by
# BiGG id (with and without the explicit D-stereochemistry suffix).
_KM_TABLE: Dict[str, float] = {
    "glc__D_e": 0.5, "glc_e": 0.5, "fru_e": 0.5, "lcts_e": 1.0,
    "ac_e": 0.7, "but_e": 0.8, "ppa_e": 0.8, "succ_e": 0.5, "lac__D_e": 0.7,
    "etoh_e": 1.0, "for_e": 1.0, "o2_e": 0.024, "co2_e": 0.5,
    "nh4_e": 0.1, "pi_e": 0.05, "so4_e": 0.05, "h2o_e": 5.0, "h_e": 1.0,
}


def _hash_unit(*parts: str) -> float:
    """A deterministic pseudo-random number in [0, 1) from string parts."""
    digest = hashlib.md5("|".join(parts).encode()).hexdigest()
    return (int(digest[:8], 16) % 1_000_000) / 1_000_000.0


def _read_fasta(path) -> Dict[str, str]:
    """Minimal FASTA reader (id = first whitespace-delimited token of the header).

    Dependency-light on purpose -- the heavy gene callers live in the workflow's
    own conda envs, so context building must not require Biopython locally.
    """
    seqs: Dict[str, str] = {}
    current = None
    chunks: List[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(chunks)
            current = line[1:].split()[0] if len(line) > 1 else ""
            chunks = []
        elif current is not None:
            chunks.append(line.strip())
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def _exchange_metabolite_ids(model) -> List[str]:
    ids = []
    for rxn in model.exchanges:
        mets = list(rxn.metabolites)
        if len(mets) == 1:
            ids.append(mets[0].id)
    return ids


def _is_internal(rxn) -> bool:
    return (
        not rxn.boundary
        and "biomass" not in rxn.id.lower()
        and getattr(rxn, "objective_coefficient", 0.0) == 0.0
    )


class HeuristicPredictor:
    """Default, dependency-free kinetic predictor (placeholders, not measurements)."""

    name = "heuristic"

    def __init__(self, median_kcat: float = MEDIAN_KCAT, default_km: float = 0.1,
                 kcat_spread: float = 1.0):
        #: lognormal-ish spread of kcat around the median, in decades
        self.median_kcat = float(median_kcat)
        self.default_km = float(default_km)
        self.kcat_spread = float(kcat_spread)

    def predict_km_for(self, metabolite_id: str) -> float:
        if metabolite_id in _KM_TABLE:
            return _KM_TABLE[metabolite_id]
        # deterministic jitter in [0.01, 1.0) for unknown metabolites
        return round(0.01 + _hash_unit("km", metabolite_id) * 0.99, 4)

    def predict_kcat_for(self, organism_id: str, reaction_id: str) -> float:
        # multiply the median by 10**U(-spread, +spread): a crude lognormal spread
        exponent = (2.0 * _hash_unit("kcat", organism_id, reaction_id) - 1.0) * self.kcat_spread
        return round(self.median_kcat * (10.0 ** exponent), 4)

    def refine(self, model, organism_id: str,
               into: Optional[KineticParameters] = None, **_: object) -> KineticParameters:
        kin = into or KineticParameters()
        for met in _exchange_metabolite_ids(model):
            vmax, _km = kin.get(organism_id, met)          # keep Vmax (uptake capacity)
            kin.set(organism_id, met, vmax, self.predict_km_for(met))
        for rxn in model.reactions:
            if _is_internal(rxn):
                kin.set_kcat(organism_id, rxn.id, self.predict_kcat_for(organism_id, rxn.id))
        return kin


def refine_kinetics(model, organism_id: str, predictor=None,
                    into: Optional[KineticParameters] = None, **kwargs) -> KineticParameters:
    """Refine an organism's kinetics with ``predictor`` (default heuristic)."""
    predictor = predictor or HeuristicPredictor()
    return predictor.refine(model, organism_id, into=into, **kwargs)


# ---------------------------------------------------------------------------
# Namespace mapping for the sequence-based (ML) predictors
# ---------------------------------------------------------------------------


def build_enzyme_context(model, proteins_faa: Optional[str] = None) -> Dict[str, dict]:
    """Map each internal reaction to its catalysing sequences + substrate SMILES.

    This is the BiGG -> (sequence, SMILES) bridge the deep-learning predictors
    need.  Gene ids are taken from each reaction's GPR and resolved to amino-acid
    sequences in ``proteins_faa`` (the per-MAG CarveMe/Prodigal output); substrate
    SMILES are read from metabolite annotations when present.  Returns
    ``{reaction_id: {"genes": [...], "sequences": {...}, "substrates": {...}}}``.
    """
    sequences: Dict[str, str] = {}
    if proteins_faa and Path(proteins_faa).exists():
        sequences = _read_fasta(proteins_faa)

    context: Dict[str, dict] = {}
    for rxn in model.reactions:
        if not _is_internal(rxn):
            continue
        genes = sorted(g.id for g in getattr(rxn, "genes", []))
        substrates = {}
        for met in rxn.metabolites:
            smiles = (met.annotation or {}).get("smiles") if hasattr(met, "annotation") else None
            if smiles:
                substrates[met.id] = smiles
        context[rxn.id] = {
            "genes": genes,
            "sequences": {g: sequences[g] for g in genes if g in sequences},
            "substrates": substrates,
        }
    return context


class DLKcatPredictor(KineticPredictor):
    """Opt-in wrapper around DLKcat (kcat from enzyme sequence + substrate SMILES).

    Requires the ``ml`` extra and a usable DLKcat checkpoint; without an enzyme
    context it cannot score reactions.  Construction fails loudly with guidance.
    """

    name = "dlkcat"

    def __init__(self, model_path: Optional[str] = None):
        try:
            import torch  # noqa: F401
        except ImportError as exc:  # pragma: no cover - exercised only without torch
            raise MuodeToolError(
                "DLKcat needs the 'ml' extra: `pip install muode[ml]` (and a DLKcat checkpoint)."
            ) from exc
        self.model_path = model_path

    def predict_kcat(self, enzyme_sequence: str, substrate_smiles: str) -> float:  # pragma: no cover
        raise MuodeToolError(
            "DLKcat scoring is not bundled; point DLKcatPredictor at a checkpoint and "
            "implement the call, or use the HeuristicPredictor default."
        )

    def refine(self, model, organism_id: str,
               into: Optional[KineticParameters] = None,
               proteins_faa: Optional[str] = None, **_: object) -> KineticParameters:  # pragma: no cover
        context = build_enzyme_context(model, proteins_faa)
        if not any(c["sequences"] for c in context.values()):
            raise MuodeToolError(
                "no enzyme sequences resolved for any reaction; pass proteins_faa from "
                "the MAG's CarveMe/Prodigal output so reactions map to sequences."
            )
        kin = into or KineticParameters()
        for rxn_id, ctx in context.items():
            seq = next(iter(ctx["sequences"].values()), None)
            smi = next(iter(ctx["substrates"].values()), None)
            if seq and smi:
                kin.set_kcat(organism_id, rxn_id, self.predict_kcat(seq, smi))
        return kin


class KmPredictor(KineticPredictor):
    """Opt-in wrapper around the Kroll et al. Km predictor (sequence + substrate)."""

    name = "km-ml"

    def __init__(self, model_path: Optional[str] = None):
        try:
            import torch  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            raise MuodeToolError(
                "The Km predictor needs the 'ml' extra: `pip install muode[ml]`."
            ) from exc
        self.model_path = model_path

    def predict_km(self, enzyme_sequence: str, substrate_smiles: str) -> float:  # pragma: no cover
        raise MuodeToolError(
            "Km-model scoring is not bundled; supply a checkpoint or use HeuristicPredictor."
        )


#: registry so the CLI/workflow can select a predictor by name
PREDICTORS = {
    "heuristic": HeuristicPredictor,
    "dlkcat": DLKcatPredictor,
    "km-ml": KmPredictor,
}


def get_predictor(name: str, **kwargs):
    """Instantiate a predictor by name ('heuristic' default; others need the ml extra)."""
    if name not in PREDICTORS:
        raise ValueError(f"unknown predictor '{name}'; available: {sorted(PREDICTORS)}")
    return PREDICTORS[name](**kwargs)
