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
from muode.provenance import Evidence, Parameter, register

#: default extracellular ids for the bile-acid pool
TAUROCHOLATE = "tca_e"   # conjugated primary (germinant)
CHOLATE = "ca_e"         # free primary
DEOXYCHOLATE = "dca_e"   # secondary (inhibitor)
LITHOCHOLATE = "lca_e"   # secondary (inhibitor; from chenodeoxycholate, modelled lumped)


# ---------------------------------------------------------------------------
# The parameters.  The MECHANISM above is real and cited.  These NUMBERS are not,
# and the difference is the whole point of muode.provenance.
#
# There is no published Vmax for bile-salt hydrolase or for the bai operon expressed
# per gram of community biomass in a gut.  The enzymes are characterised in vitro on
# purified protein; what we need is a whole-community turnover rate, and nobody has
# measured it.  So these are INVENTED -- not "roughly right", not "conservative":
# nobody chose them.
#
# That is survivable ONLY if the conclusion does not depend on them, and that is a
# question with an answer.  muode.provenance.sensitivity() sweeps a parameter across
# the range you cannot rule out and reports whether the result moves.  Any claim that
# rests on these rates must cite that sweep, or it is not a claim.
# ---------------------------------------------------------------------------

VMAX_BSH = register("bile", Parameter(
    name="vmax_bsh",
    value=5.0,
    units="mmol/gDW/h",
    evidence=Evidence.INVENTED,
    why=(
        "Community-level deconjugation rate by bile-salt hydrolase. NOBODY CHOSE THIS "
        "NUMBER. BSH is well characterised in vitro on purified enzyme, but no one has "
        "measured a per-gram-of-community turnover in a gut, which is what this layer "
        "needs. Known qualitatively: BSH is COMMON (Heinken 2019 finds it in 204/693 "
        "gut genomes, 29%), so deconjugation is fast relative to 7a-dehydroxylation -- "
        "which is why vmax_bsh > vmax_bai. The ORDERING is defensible; the VALUES are not."
    ),
))

VMAX_BAI = register("bile", Parameter(
    name="vmax_bai",
    value=1.0,
    units="mmol/gDW/h",
    evidence=Evidence.INVENTED,
    why=(
        "Community-level 7a-dehydroxylation rate (bai operon). NOBODY CHOSE THIS NUMBER. "
        "Set 5x below vmax_bsh because the guild is RARE -- the bai operon appears in "
        "only 7/693 gut genomes, ~1% (Heinken 2019) -- so the flux to secondary bile "
        "acids is the bottleneck, which is the entire mechanism of colonisation "
        "resistance. Again: the ordering is the claim, not the value."
    ),
))

KM_BSH = register("bile", Parameter(
    name="km_bsh", value=0.05, units="mM", evidence=Evidence.INVENTED,
    why="Half-saturation for BSH on taurocholate at the community level. Not measured.",
))

KM_BAI = register("bile", Parameter(
    name="km_bai", value=0.05, units="mM", evidence=Evidence.INVENTED,
    why="Half-saturation for 7a-dehydroxylation on cholate at the community level. "
        "Not measured.",
))

#: The one bile parameter that IS anchored -- see lifecycle.KI_INHIBITOR for the twin.
KI_BILE = register("bile", Parameter(
    name="ki",
    value=0.5,
    units="mM",
    evidence=Evidence.DERIVED,
    why=(
        "Half-maximal inhibition of C. difficile vegetative growth by secondary bile "
        "acids. Deoxycholate delays growth at 0.01% w/v (0.24 mM) and abolishes it at "
        "0.05% (1.2 mM), so the IC50 lies between; 0.5 mM is the log midpoint. Was "
        "0.02 mM -- 25x too potent, which made colonisation resistance look far more "
        "robust than the chemistry supports."
    ),
    citation="Usui Y et al. Heliyon 6:e03717 (2020), doi:10.1016/j.heliyon.2020.e03717",
))


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
    vmax_bsh: float = VMAX_BSH.value
    km_bsh: float = KM_BSH.value
    vmax_bai: float = VMAX_BAI.value
    km_bai: float = KM_BAI.value
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
    ki: float = KI_BILE.value
    secondary: Tuple[str, ...] = (DEOXYCHOLATE, LITHOCHOLATE)

    def growth_factor(self, organism_id, t, M, X) -> float:
        if organism_id not in self.targets:
            return 1.0
        s = sum(max(0.0, M.get(b, 0.0)) for b in self.secondary)
        return 1.0 / (1.0 + s / self.ki) if self.ki > 0 else 1.0
