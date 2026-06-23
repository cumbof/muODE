"""Kinetic parameters for the dynamic-FBA loop.

This module is where the README's *Phase 3* is made mathematically precise,
because the original description conflated two distinct uses of "kinetics":

1. **Substrate uptake kinetics (Km, Vmax).**  In dynamic FBA the *upper bound*
   of every substrate uptake reaction is recomputed at each time step from the
   extracellular concentration via Michaelis-Menten / Monod kinetics::

       v_uptake(M) = Vmax * M / (Km + M)

   These are the parameters that actually couple the intracellular LP to the
   extracellular ODEs, so they are *required* by the engine.  When no predicted
   value is available we fall back to literature-style defaults (a single
   permissive Vmax and a small Km), which reduces the loop to the classic
   substrate-unlimited dFBA used throughout the field.

2. **Enzyme turnover (kcat).**  A kcat does **not** enter the uptake bound; it
   bounds an *intracellular* reaction's maximal velocity through the enzyme
   it is catalysed by: ``Vmax = kcat * [E]``.  This is the realm of
   enzyme-constrained models (GECKO / sMOMENT).  muODE treats this as an
   optional refinement layer (see :func:`vmax_from_kcat`) rather than a
   prerequisite for a working simulation.

The deep-learning predictors named in the README (DLKcat, the Kroll et al. Km
model, RealKcat, ...) plug in through :class:`KineticPredictor`.  They are
optional: the pipeline runs end-to-end on defaults and the predictions only
*sharpen* the bounds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

# Sensible, deliberately permissive defaults (mmol / gDW / h and mmol / L).
# A large default Vmax recovers the substrate-unlimited dFBA limit; a small Km
# means uptake saturates quickly once a substrate is present.
DEFAULT_VMAX = 10.0
DEFAULT_KM = 0.01


@dataclass
class KineticParameters:
    """Container resolving ``(organism, metabolite) -> (Vmax, Km)``.

    Lookups fall back gracefully: a per-(organism, metabolite) override wins,
    then a per-metabolite default, then the global default.
    """

    default_vmax: float = DEFAULT_VMAX
    default_km: float = DEFAULT_KM
    #: per-metabolite defaults, applied to every organism
    metabolite_defaults: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    #: per-(organism, metabolite) overrides (e.g. from Km/Vmax predictions)
    overrides: Dict[Tuple[str, str], Tuple[float, float]] = field(default_factory=dict)
    #: per-(organism, reaction) turnover numbers kcat (1/s), for enzyme constraints
    kcat: Dict[Tuple[str, str], float] = field(default_factory=dict)
    #: default enzyme abundance (mmol enzyme / gDW) for the kcat -> Vmax conversion
    default_enzyme_concentration: float = 1e-2

    def get(self, organism_id: str, metabolite_id: str) -> Tuple[float, float]:
        """Return ``(Vmax, Km)`` for an organism/metabolite uptake reaction."""
        if (organism_id, metabolite_id) in self.overrides:
            return self.overrides[(organism_id, metabolite_id)]
        if metabolite_id in self.metabolite_defaults:
            return self.metabolite_defaults[metabolite_id]
        return (self.default_vmax, self.default_km)

    def set(
        self,
        organism_id: str,
        metabolite_id: str,
        vmax: float,
        km: float,
    ) -> None:
        self.overrides[(organism_id, metabolite_id)] = (float(vmax), float(km))

    def michaelis_menten(
        self, organism_id: str, metabolite_id: str, concentration: float
    ) -> float:
        """Maximum uptake rate at a given extracellular concentration."""
        if concentration <= 0.0:
            return 0.0
        vmax, km = self.get(organism_id, metabolite_id)
        return vmax * concentration / (km + concentration)

    # -- enzyme turnover (kcat) --------------------------------------------
    def set_kcat(self, organism_id: str, reaction_id: str, kcat: float) -> None:
        self.kcat[(organism_id, reaction_id)] = float(kcat)

    def get_kcat(self, organism_id: str, reaction_id: str) -> Optional[float]:
        return self.kcat.get((organism_id, reaction_id))

    def enzyme_vmax(
        self,
        organism_id: str,
        reaction_id: str,
        enzyme_concentration: Optional[float] = None,
    ) -> Optional[float]:
        """Enzyme-constrained Vmax (mmol/gDW/h) for an intracellular reaction.

        Returns ``None`` when no kcat is known for the reaction, so callers can
        leave such reactions unconstrained.
        """
        kcat = self.get_kcat(organism_id, reaction_id)
        if kcat is None:
            return None
        e = self.default_enzyme_concentration if enzyme_concentration is None else enzyme_concentration
        return vmax_from_kcat(kcat, e)

    # -- combination & persistence -----------------------------------------
    def merge(self, other: "KineticParameters") -> "KineticParameters":
        """Fold another store's per-organism Km/Vmax and kcat into this one."""
        self.metabolite_defaults.update(other.metabolite_defaults)
        self.overrides.update(other.overrides)
        self.kcat.update(other.kcat)
        return self

    def to_dict(self) -> dict:
        """JSON-serialisable view (tuple keys are nested by organism)."""
        overrides: Dict[str, Dict[str, list]] = {}
        for (org, met), vk in self.overrides.items():
            overrides.setdefault(org, {})[met] = list(vk)
        kcat: Dict[str, Dict[str, float]] = {}
        for (org, rxn), k in self.kcat.items():
            kcat.setdefault(org, {})[rxn] = k
        return {
            "default_vmax": self.default_vmax,
            "default_km": self.default_km,
            "default_enzyme_concentration": self.default_enzyme_concentration,
            "metabolite_defaults": {m: list(vk) for m, vk in self.metabolite_defaults.items()},
            "overrides": overrides,
            "kcat": kcat,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "KineticParameters":
        kp = cls(
            default_vmax=data.get("default_vmax", DEFAULT_VMAX),
            default_km=data.get("default_km", DEFAULT_KM),
            default_enzyme_concentration=data.get("default_enzyme_concentration", 1e-2),
            metabolite_defaults={
                m: tuple(vk) for m, vk in data.get("metabolite_defaults", {}).items()
            },
        )
        for org, mets in data.get("overrides", {}).items():
            for met, vk in mets.items():
                kp.overrides[(org, met)] = tuple(vk)
        for org, rxns in data.get("kcat", {}).items():
            for rxn, k in rxns.items():
                kp.kcat[(org, rxn)] = float(k)
        return kp

    def to_json(self, path) -> None:
        import json
        from pathlib import Path

        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def from_json(cls, path) -> "KineticParameters":
        import json
        from pathlib import Path

        return cls.from_dict(json.loads(Path(path).read_text()))


def vmax_from_kcat(kcat: float, enzyme_concentration: float) -> float:
    """Enzyme-constrained Vmax from a turnover number and enzyme abundance.

    ``kcat`` in 1/s, ``enzyme_concentration`` in mmol/gDW -> Vmax in mmol/gDW/h.
    """
    return kcat * 3600.0 * enzyme_concentration


class KineticPredictor:
    """Optional adapter to deep-learning kinetic predictors.

    Concrete subclasses wrap DLKcat / RealKcat (turnover numbers) or the Kroll
    et al. Km model.  They are intentionally *not* imported at module load: the
    heavy ``torch`` dependency lives behind the ``ml`` extra and is only touched
    when a user explicitly opts into prediction.
    """

    name = "abstract"

    def predict_kcat(self, enzyme_sequence: str, substrate_smiles: str) -> float:
        raise NotImplementedError(
            "Install the 'ml' extra and use a concrete predictor (e.g. DLKcatPredictor)."
        )

    def predict_km(self, enzyme_sequence: str, substrate_smiles: str) -> float:
        raise NotImplementedError(
            "Install the 'ml' extra and use a concrete predictor (e.g. KmPredictor)."
        )
