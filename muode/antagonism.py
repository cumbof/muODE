"""Direct interference competition -- bacteriocins and contact-independent toxins.

Cross-feeding and nutrient competition are *exploitative* interactions, mediated
by shared metabolites the FBA core already tracks.  Communities also engage in
**interference competition**: a producer secretes a toxin (a bacteriocin, a
microcin, an antimicrobial peptide) that suppresses specific competitors without
any metabolic exchange.  Examples relevant to the gut include *Blautia*/
*Ruminococcaceae* lantibiotics implicated in resisting *C. difficile* and
Enterococcus, and *Lactobacillus* bacteriocins.

This layer models a diffusible inhibitor as an extra pool species:

* **production** -- first-order in the summed biomass of the ``producers``;
* **decay** -- first-order turnover of the inhibitor;
* **action** -- Hill (n=1) suppression of the ``targets`` growth rate,
  ``f = 1 / (1 + C / ki)``.

It is deliberately generic: the same layer covers any "species A chemically
suppresses species B" interaction the user can parameterise.  What it does *not*
attempt is contact-dependent (T6SS) killing, which is spatial and would belong in
the :mod:`muode.spatial` engine.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable

from muode.ecology import EcologyLayer


@dataclass
class Bacteriocin(EcologyLayer):
    """A diffusible inhibitor secreted by ``producers`` and suppressing ``targets``.

    Parameters
    ----------
    producers, targets:
        organism ids that secrete / are inhibited by the toxin.
    metabolite:
        pool id used to track the toxin concentration.
    production:
        secretion rate (concentration units / gDW / h), first-order in producer
        biomass.
    decay:
        first-order turnover of the toxin (1/h).
    ki:
        toxin concentration giving half-maximal growth suppression of targets.
    """

    name: str = "bacteriocin"
    producers: set = field(default_factory=set)
    targets: set = field(default_factory=set)
    metabolite: str = "bacteriocin_e"
    production: float = 0.5
    decay: float = 0.2
    ki: float = 0.1

    def extra_metabolites(self) -> Iterable[str]:
        return (self.metabolite,)

    def metabolite_rates(self, t, M, X) -> Dict[str, float]:
        x_prod = sum(X.get(s, 0.0) for s in self.producers)
        produced = self.production * x_prod
        decayed = self.decay * max(0.0, M.get(self.metabolite, 0.0))
        return {self.metabolite: produced - decayed}

    def growth_factor(self, organism_id, t, M, X) -> float:
        if organism_id not in self.targets:
            return 1.0
        c = max(0.0, M.get(self.metabolite, 0.0))
        return 1.0 / (1.0 + c / self.ki) if self.ki > 0 else 1.0
