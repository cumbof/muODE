"""Sporulation / germination life cycle.

A growing-cell-only model cannot explain the central clinical fact about
*C. difficile*: it **recurs**.  The reason is a dormant **spore** state that

* does not grow and consumes nothing,
* is **resistant to antibiotics** (and to acid/bile), and
* **germinates** back into a vegetative cell when it senses a germinant
  (taurocholate + glycine) *in the absence* of secondary-bile-acid inhibitors.

This layer splits a spore-former's biomass into the vegetative pool (the ordinary
``X`` the FBA core grows and the antibiotic layer kills) and an internal **spore
pool** that the engine never exposes to growth or death.  Per step:

    sporulation:  veg --(k_spo, when growth is stressed)--> spore
    germination:  spore --(k_germ * germination_signal)----> veg

The ``germination_signal`` in [0, 1] rises with the germinant concentration
(Michaelis-Menten) and is shut off by secondary bile acids -- so the same bile
pool that :mod:`muode.bile` controls also gates germination here.  Together with
:mod:`muode.antibiotic` this reproduces the recurrence cycle: antibiotics clear
vegetative cells but not spores; if the post-antibiotic pool is germinant-rich
and secondary-bile-acid-poor (a wiped community), spores germinate and the
infection returns; an FMT that restores the bile-transforming guild keeps the
pool inhibitory and the spores dormant.

References: Sorg & Sonenshein (2008) J. Bacteriol. (taurocholate germinant);
Theriot et al. (2014); Buffie et al. (2015).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Tuple

from muode.bile import DEOXYCHOLATE, LITHOCHOLATE, TAUROCHOLATE
from muode.ecology import EcologyLayer
from muode.provenance import Evidence, Parameter, register


#: WHERE THESE NUMBERS COME FROM -- and what they used to be.
#:
#: km_germinant was **0.1 mM** and ki_inhibitor was **0.02 mM**.  Neither was measured;
#: neither was cited; both were floats that made the scenario run.  And both are wrong
#: against numbers that have been in the literature for over a decade -- wrong by 159x
#: and ~25x respectively, and wrong in the SAME DIRECTION: spores germinating far too
#: readily, secondary bile acids suppressing far too strongly.  That is the whole rCDI
#: mechanism, biased toward the answer we were hoping to see.
#:
#: This is what muode.provenance exists to make impossible to do quietly.
KM_GERMINANT = register("lifecycle", Parameter(
    name="km_germinant",
    value=15.9,
    units="mM",
    evidence=Evidence.MEASURED,
    why=(
        "Half-maximal germination of C. difficile spores by taurocholate. Measured "
        "EC50 = 15.9 mM. NB this is taurocholate ALONE: the same work shows strong "
        "synergy with glycine and Ca2+, which at physiological levels lowers the "
        "effective threshold substantially, so 15.9 is an UPPER bound on the true in "
        "vivo half-saturation and the model is correspondingly conservative about "
        "germination. Bile reaches the caecum at ~2 mM. Was 0.1 mM -- 159x too low, "
        "i.e. spores germinated on almost any trace of germinant."
    ),
    citation="Ramirez N, Abel-Santos E. J Bacteriol 193:531-539 (2011); "
             "Shrestha & Sorg, mSphere 3:e00335-18 (2018) for the synergy",
))

KI_INHIBITOR = register("lifecycle", Parameter(
    name="ki_inhibitor",
    value=0.5,
    units="mM",
    evidence=Evidence.DERIVED,
    why=(
        "Half-maximal block of germination/growth by secondary bile acids. Derived "
        "from a dose-response: deoxycholate delays C. difficile vegetative growth at "
        "0.01% w/v (0.24 mM) and abolishes it completely at 0.05% (1.2 mM), so the "
        "half-inhibitory concentration lies between them and 0.5 mM is the midpoint "
        "on a log scale. Was 0.02 mM -- ~25x too potent, so a trace of DCA shut the "
        "pathogen down and colonisation resistance looked far more robust than the "
        "chemistry supports."
    ),
    citation="Usui Y et al. Heliyon 6:e03717 (2020), doi:10.1016/j.heliyon.2020.e03717",
))


@dataclass
class SporeForming(EcologyLayer):
    """Two-compartment (vegetative / spore) life cycle for spore-forming species.

    Parameters
    ----------
    species:
        spore-forming organism ids.
    initial_spores:
        starting dormant biomass (gDW/L) per species -- e.g. an ingested spore
        inoculum, or the dormant reservoir left by a prior antibiotic course.
    k_sporulation, mu_stress:
        vegetative cells convert to spores at rate ``k_sporulation`` (1/h) when
        their growth rate falls below ``mu_stress`` (a stress trigger).
    k_germination:
        max germination rate (1/h), scaled by the germination signal.
    k_spore_decay:
        first-order spore loss (1/h); spores are durable, so ~0 by default.
    germinant, km_germinant:
        pool metabolite that triggers germination and its half-saturation.
    inhibitor, ki_inhibitor:
        pool metabolites (secondary bile acids) that block germination.
    """

    name: str = "spore_forming"
    species: set = field(default_factory=set)
    initial_spores: Dict[str, float] = field(default_factory=dict)
    k_sporulation: float = 0.5
    mu_stress: float = 0.05
    k_germination: float = 0.4
    k_spore_decay: float = 0.0
    germinant: str = TAUROCHOLATE
    km_germinant: float = KM_GERMINANT.value
    inhibitor: Tuple[str, ...] = (DEOXYCHOLATE, LITHOCHOLATE)
    ki_inhibitor: float = KI_INHIBITOR.value
    _spores: Dict[str, float] = field(default_factory=dict, init=False, repr=False)

    def extra_metabolites(self) -> Tuple[str, ...]:
        """The germinant and its inhibitors: read from the medium, so declare them.

        No organism *exchanges* these -- they gate germination rather than being
        consumed -- but they are not inert, and an undeclared metabolite looks to
        the engine like a diet entry that feeds nothing.
        """
        return (self.germinant, *self.inhibitor)

    def reset(self, community) -> None:
        self._spores = {s: float(self.initial_spores.get(s, 0.0)) for s in self.species}

    def spores(self) -> Dict[str, float]:
        return dict(self._spores)

    def germination_signal(self, M: Mapping[str, float]) -> float:
        g = max(0.0, M.get(self.germinant, 0.0))
        trigger = g / (self.km_germinant + g) if g > 0 else 0.0
        inh = sum(max(0.0, M.get(b, 0.0)) for b in self.inhibitor)
        block = 1.0 / (1.0 + inh / self.ki_inhibitor) if self.ki_inhibitor > 0 else 1.0
        return trigger * block

    def integrate(self, t, dt, X, M, mu) -> None:
        signal = self.germination_signal(M)
        for s in self.species:
            veg = X.get(s, 0.0)
            spore = self._spores.get(s, 0.0)
            growth = mu.get(s, 0.0)

            d_spo = (self.k_sporulation * veg * dt) if growth < self.mu_stress else 0.0
            d_germ = self.k_germination * spore * signal * dt
            d_decay = self.k_spore_decay * spore * dt

            X[s] = max(0.0, veg - d_spo + d_germ)
            self._spores[s] = max(0.0, spore + d_spo - d_germ - d_decay)

    def observables(self, t, M, X) -> Dict[str, float]:
        return {"germination_signal": self.germination_signal(M)}
