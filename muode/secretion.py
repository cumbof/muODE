"""Growth-coupled metabolite secretion.

A phenomenological *source* layer: each species releases one or more metabolites in
proportion to its **growth** (its biomass increment), at a per-species yield.  It is
the controlled stand-in for the fermentation / overflow that the FBA biomass objective
does not route to secretion -- give it yields measured from monoculture and the shared
pool accumulates the products realistically, *without* hand-tuning to the community
outcome.

This generalises a pattern that kept being re-invented per study (SCFA production in
the FMT work; ``GrowthCoupledAcid`` in the Clark benchmark): any metabolite, any guild.

Pair it with the consumers/inhibitors that read the pool -- e.g. :class:`muode.ph.
WeakAcidInhibition` for an acidification feedback, or a GEM that takes the product up
as a cross-feeding substrate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping

from muode.ecology import EcologyLayer


@dataclass
class GrowthCoupledSecretion(EcologyLayer):
    """Secrete metabolites in proportion to each species' growth increment.

    ``yields[species_id][metabolite_id]`` is the amount of the metabolite (mmol/L)
    released per unit of that species' biomass.  Each step the layer tracks the
    per-species biomass change ``dX`` since the previous step and adds ``yield * dX``
    to the pool -- implemented as a rate ``yield * dX / dt`` that the engine multiplies
    back by ``dt`` when it integrates ``metabolite_rates``.  Only positive increments
    secrete (a shrinking population does not release product).

    Parameters
    ----------
    yields:
        ``{species_id: {metabolite_id: amount_per_unit_biomass}}``.
    dt:
        the integration step.  Because the returned value is a rate (divided by ``dt``)
        that the engine multiplies by ``dt``, this **must match** the ``dt`` passed to
        :class:`~muode.dfba.DynamicFBA`, or the amount secreted per unit growth is
        mis-scaled.

    Notes
    -----
    ``reset`` seeds the per-species baseline at 0, so the first step's increment
    includes the starting biomass (a one-off startup release of ``yield * X0``); this
    matches the original benchmark behaviour.
    """

    name: str = "growth_coupled_secretion"
    yields: Dict[str, Dict[str, float]] = field(default_factory=dict)
    dt: float = 0.1
    _last: Dict[str, float] = field(default_factory=dict)

    def extra_metabolites(self) -> Iterable[str]:
        mets: set = set()
        for y in self.yields.values():
            mets |= set(y)
        return tuple(sorted(mets))

    def reset(self, community) -> None:
        self._last = {o.id: 0.0 for o in community.organisms}

    def metabolite_rates(self, t, M, X) -> Dict[str, float]:
        rates: Dict[str, float] = {}
        for sp, y in self.yields.items():
            if sp not in X:
                continue
            dX = X[sp] - self._last.get(sp, 0.0)
            self._last[sp] = X[sp]
            if dX <= 0:
                continue
            for met, k in y.items():
                rates[met] = rates.get(met, 0.0) + k * dX / self.dt   # engine multiplies by dt
        return rates
