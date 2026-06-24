"""Bile-acid metabolism and bile-acid inhibition.

Bile acids are the second major arm of colonization resistance against
*Clostridioides difficile*, and they are a *chemical transformation of the
shared pool* driven by specific enzymes -- exactly the kind of process the FBA
biomass objective does not represent.

The relevant pathway:

* the liver secretes **conjugated primary** bile acids (taurocholate, ``tca_e``);
* broad **bile-salt hydrolase (BSH)** activity -- common across *Bacteroides*,
  *Bifidobacterium*, *Lactobacillus*, *Clostridium* -- deconjugates them to
  **free primary** bile acids (cholate, ``ca_e``);
* a *rare* **7α-dehydroxylating** guild carrying the *bai* operon (e.g.
  *Clostridium scindens*) converts free primary into **secondary** bile acids
  (deoxycholate, ``dca_e``);
* **taurocholate is the physiological germinant** for *C. difficile* spores,
  while **secondary bile acids inhibit** both germination and vegetative growth.

So a healthy community that carries the *bai* guild keeps the pool shifted toward
secondary bile acids -- germination is blocked and the pathogen is suppressed.
Antibiotics that wipe the *bai* guild leave a taurocholate-rich, secondary-poor
pool that *permits* germination: the mechanistic basis of recurrence, and of why
FMT (which restores the guild) cures it (Buffie et al. 2015 Nature; Theriot et
al. 2014 Nat. Commun.).

:class:`BileAcidTransform` evolves the bile-acid pool; :class:`BileAcidInhibition`
applies secondary-bile-acid toxicity to vegetative growth.  Germination is gated
by these same pools in :mod:`muode.lifecycle`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping, Tuple

from muode.ecology import EcologyLayer

#: default extracellular ids for the bile-acid pool
TAUROCHOLATE = "tca_e"   # conjugated primary (germinant)
CHOLATE = "ca_e"         # free primary
DEOXYCHOLATE = "dca_e"   # secondary (inhibitor)
LITHOCHOLATE = "lca_e"   # secondary (inhibitor; from chenodeoxycholate, modelled lumped)


@dataclass
class BileAcidTransform(EcologyLayer):
    """Microbial bile-acid transformation in the shared pool.

    Two enzymatic steps, each Michaelis-Menten in the substrate and first-order
    in the (summed) biomass of the guild that carries the enzyme:

        taurocholate --BSH-->            cholate
        cholate      --7α-dehydroxylase->deoxycholate

    Parameters
    ----------
    bsh_producers:
        species ids with bile-salt hydrolase (broad).
    bai_producers:
        species ids with the *bai* operon (rare 7α-dehydroxylating guild).
    """

    name: str = "bile_acid_transform"
    bsh_producers: set = field(default_factory=set)
    bai_producers: set = field(default_factory=set)
    vmax_bsh: float = 5.0
    km_bsh: float = 0.05
    vmax_bai: float = 1.0
    km_bai: float = 0.05
    tca: str = TAUROCHOLATE
    ca: str = CHOLATE
    dca: str = DEOXYCHOLATE

    def extra_metabolites(self) -> Iterable[str]:
        return (self.tca, self.ca, self.dca)

    def metabolite_rates(self, t, M, X) -> Dict[str, float]:
        x_bsh = sum(X.get(s, 0.0) for s in self.bsh_producers)
        x_bai = sum(X.get(s, 0.0) for s in self.bai_producers)
        tca = max(0.0, M.get(self.tca, 0.0))
        ca = max(0.0, M.get(self.ca, 0.0))
        v_bsh = self.vmax_bsh * x_bsh * tca / (self.km_bsh + tca) if tca > 0 else 0.0
        v_bai = self.vmax_bai * x_bai * ca / (self.km_bai + ca) if ca > 0 else 0.0
        return {
            self.tca: -v_bsh,
            self.ca: v_bsh - v_bai,
            self.dca: v_bai,
        }


@dataclass
class BileAcidInhibition(EcologyLayer):
    """Secondary bile acids inhibit vegetative growth of target species.

    ``f = 1 / (1 + sum(secondary) / ki)`` for organisms in ``targets`` (others
    unaffected).  Deoxycholate is potent against *C. difficile* at low
    micromolar-to-millimolar levels, so ``ki`` is small.
    """

    name: str = "bile_acid_inhibition"
    targets: set = field(default_factory=set)
    ki: float = 0.02
    secondary: Tuple[str, ...] = (DEOXYCHOLATE, LITHOCHOLATE)

    def growth_factor(self, organism_id, t, M, X) -> float:
        if organism_id not in self.targets:
            return 1.0
        s = sum(max(0.0, M.get(b, 0.0)) for b in self.secondary)
        return 1.0 / (1.0 + s / self.ki) if self.ki > 0 else 1.0
