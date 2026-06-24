"""Ecology layer -- environmental & life-history processes around the dFBA core.

The dynamic-FBA engine (:mod:`muode.dfba`) models metabolism: each species is an
LP, and they couple through a shared pool of extracellular metabolites.  That is
the right level for *nutrient competition* and *cross-feeding*, but a real gut
community is shaped by processes that are **not** metabolic fluxes:

* the **chemistry of the medium** -- short-chain fatty acids acidify it, and the
  undissociated acid inhibits growth (a pH effect, not a flux);
* **chemical / enzymatic transformation of the pool** that no single species'
  biomass objective captures -- e.g. bile-acid 7α-dehydroxylation;
* **direct antagonism** -- bacteriocins, not cross-feeding;
* **pharmacology** -- a time-varying antibiotic concentration;
* **life history** -- sporulation and germination, a non-growth state.

These are layered *on top of* the FBA core as an :class:`EcologyModel`: a list of
:class:`EcologyLayer` objects, each hooking into well-defined points of the
integration loop.  Every hook has a no-op default, so an empty ``EcologyModel``
leaves the simulation numerically identical to the bare engine.  Concrete layers
live in :mod:`muode.ph`, :mod:`muode.bile`, :mod:`muode.lifecycle`,
:mod:`muode.antibiotic` and :mod:`muode.antagonism`.

Design intent: every layer is a *real, literature-grounded mechanism* with
tunable parameters, not a fudge factor.  What remains genuinely outside this
paradigm (gene regulation, evolution, agent-level stochasticity) is documented
in ``docs/LIMITATIONS.md`` rather than faked here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Iterable, List, Mapping, Tuple

if TYPE_CHECKING:  # pragma: no cover
    from muode.community import Community


class EcologyLayer:
    """One environmental / life-history process.

    Subclasses override only the hooks they need.  The hooks, and when the dFBA
    loop calls them, are:

    ``extra_metabolites`` / ``initial_concentrations``
        declare pool species the layer needs tracked that may not be exchanged
        by any GEM (e.g. bile acids, a bacteriocin) and their starting levels.
    ``reset(community)``
        called once before the run (initialise internal state, e.g. spore pools).
    ``uptake_factor(org, met, t, M, X)``  -> multiplies a substrate uptake bound.
    ``growth_factor(org, t, M, X)``       -> multiplies a species' growth rate.
    ``extra_death(org, t, M, X)``         -> adds to a species' death rate (1/h).
    ``metabolite_rates(t, M, X)``         -> extra dM/dt (mmol/L/h) for the pool.
    ``integrate(t, dt, X, M, mu)``        -> may move biomass (e.g. veg<->spore).
    ``observables(t, M, X)``              -> scalars to record (e.g. pH).
    ``spores()``                          -> per-species dormant biomass to record.
    """

    name: str = "layer"

    def extra_metabolites(self) -> Iterable[str]:
        return ()

    def initial_concentrations(self) -> Mapping[str, float]:
        return {}

    def reset(self, community: "Community") -> None:
        return None

    def uptake_factor(self, organism_id: str, metabolite_id: str, t: float,
                      M: Mapping[str, float], X: Mapping[str, float]) -> float:
        return 1.0

    def growth_factor(self, organism_id: str, t: float,
                      M: Mapping[str, float], X: Mapping[str, float]) -> float:
        return 1.0

    def extra_death(self, organism_id: str, t: float,
                    M: Mapping[str, float], X: Mapping[str, float]) -> float:
        return 0.0

    def metabolite_rates(self, t: float, M: Mapping[str, float],
                         X: Mapping[str, float]) -> Mapping[str, float]:
        return {}

    def integrate(self, t: float, dt: float, X: Dict[str, float],
                  M: Mapping[str, float], mu: Mapping[str, float]) -> None:
        return None

    def observables(self, t: float, M: Mapping[str, float],
                    X: Mapping[str, float]) -> Mapping[str, float]:
        return {}

    def spores(self) -> Mapping[str, float]:
        return {}

    def describe(self) -> str:
        return self.name


@dataclass
class EcologyModel:
    """A composable stack of :class:`EcologyLayer` objects.

    Aggregates the layers' hooks the way the dFBA loop expects: growth/uptake
    factors *multiply*, death rates and metabolite rates *sum*, declared
    metabolites *union*, observables/spores *merge*.
    """

    layers: List[EcologyLayer] = field(default_factory=list)

    def add(self, layer: EcologyLayer) -> "EcologyModel":
        self.layers.append(layer)
        return self

    def is_empty(self) -> bool:
        return not self.layers

    # -- declarations -------------------------------------------------------
    def extra_metabolites(self) -> set:
        out: set = set()
        for layer in self.layers:
            out |= set(layer.extra_metabolites())
        return out

    def initial_concentrations(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for layer in self.layers:
            out.update(layer.initial_concentrations())
        return out

    def reset(self, community: "Community") -> None:
        for layer in self.layers:
            layer.reset(community)

    # -- per-step modifiers -------------------------------------------------
    def uptake_factor(self, organism_id, metabolite_id, t, M, X) -> float:
        f = 1.0
        for layer in self.layers:
            f *= layer.uptake_factor(organism_id, metabolite_id, t, M, X)
        return f

    def growth_factor(self, organism_id, t, M, X) -> float:
        f = 1.0
        for layer in self.layers:
            f *= layer.growth_factor(organism_id, t, M, X)
        return f

    def extra_death(self, organism_id, t, M, X) -> float:
        return sum(layer.extra_death(organism_id, t, M, X) for layer in self.layers)

    def metabolite_rates(self, t, M, X) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for layer in self.layers:
            for met, rate in layer.metabolite_rates(t, M, X).items():
                out[met] = out.get(met, 0.0) + rate
        return out

    def integrate(self, t, dt, X, M, mu) -> None:
        for layer in self.layers:
            layer.integrate(t, dt, X, M, mu)

    # -- recording ----------------------------------------------------------
    def observables(self, t, M, X) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for layer in self.layers:
            out.update(layer.observables(t, M, X))
        return out

    def spores(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for layer in self.layers:
            out.update(layer.spores())
        return out

    def describe(self) -> List[str]:
        return [layer.describe() for layer in self.layers]
