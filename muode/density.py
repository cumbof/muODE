"""Density-dependent (logistic) self-limitation -- carrying-capacity growth control.

The FBA core sets an *instantaneous* growth rate from the current medium; nothing in
it stops a fast grower from compounding exponentially and monopolizing a shared
substrate.  In a closed community that produces competitive exclusion far sharper
than reality: on the Clark et al. 2021 DM38 panel the per-step max-biomass dFBA
drives the minority member of a pair to ~0.13 of the community where the measured
value is ~0.26, and a controlled experiment (capping each organism's growth rate)
showed the over-exclusion is a *growth-rate* phenomenon -- moderate the winner's
rate and coexistence is restored.

This layer is the principled form of that moderation: classic Verhulst logistic
self-limitation.  Each organism's growth rate is multiplied by

    f_i = max(floor, 1 - X_i / K_i)

so growth slows to zero as the organism's biomass ``X_i`` approaches its carrying
capacity ``K_i`` and the fast grower stops taking everything once it nears its own
capacity.  The natural, data-grounded choice for ``K_i`` is the strain's calibrated
**monoculture yield** -- in a defined medium the endpoint OD *is* the carrying
capacity of that strain on that medium -- so a community settles near the
yield-proportional composition that the measurements actually show.  Because the
capacities are per-organism (not a single shared cap), this moderates each species
by *its own* ceiling rather than slowing everyone uniformly (which cancels in
relative abundance and does nothing).

``floor`` (default 0.0) bounds the factor from below: 0.0 means growth simply
arrests at ``K`` (no density-induced death), which is stable under the explicit
Euler step; a small negative floor would let an over-shooting population decline
back toward ``K`` at the cost of possible oscillation.  Organisms with no capacity
supplied (and no ``default_capacity``) are left unconstrained (factor 1.0).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional

from muode.ecology import EcologyLayer


@dataclass
class LogisticCarryingCapacity(EcologyLayer):
    """Logistic (density-dependent) growth self-limitation, per organism.

    Parameters
    ----------
    capacity:
        Per-organism carrying capacity ``K_i`` in the same units as biomass
        (gDW/L).  In a defined medium set it to the strain's calibrated
        monoculture yield / measured endpoint OD.
    floor:
        Lower bound on the growth factor (default 0.0 -> growth arrests at ``K``;
        a small negative value allows decline above ``K``).
    default_capacity:
        Capacity for organisms absent from ``capacity``.  ``None`` (default)
        leaves them unconstrained (factor 1.0).
    """

    name: str = "logistic_carrying_capacity"
    capacity: Dict[str, float] = field(default_factory=dict)
    floor: float = 0.0
    default_capacity: Optional[float] = None

    def growth_factor(self, organism_id, t, M, X) -> float:
        K = self.capacity.get(organism_id, self.default_capacity)
        if K is None or K <= 0.0:
            return 1.0
        return max(self.floor, 1.0 - X.get(organism_id, 0.0) / K)

    def observables(self, t, M, X) -> Dict[str, float]:
        """Record each organism's saturation fraction X_i / K_i (0 = empty, 1 = full)."""
        out: Dict[str, float] = {}
        for org, K in self.capacity.items():
            if K and K > 0.0:
                out[f"saturation_{org}"] = X.get(org, 0.0) / K
        return out

    def describe(self) -> str:
        n = len(self.capacity)
        return f"{self.name}(n_capacities={n}, floor={self.floor})"
