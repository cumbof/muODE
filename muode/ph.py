"""pH dynamics and weak-acid (SCFA) growth inhibition.

Short-chain fatty acids are the dominant fermentation products of a gut
community, and they are **weak acids**.  Two coupled effects follow, neither of
which is a metabolic flux (so neither is captured by the FBA core):

1. **Acidification.**  As SCFAs accumulate they titrate the medium and lower the
   pH.  We use a buffered linear titration model: ``pH = pH0 - (total acid) /
   buffer_capacity`` (clamped at a floor).  This is the coarse-grained model used
   in dynamic gut-microbiome models; ``pH0`` and ``buffer_capacity`` are tunable
   to a given compartment (proximal colon pH ~5.5-6.5, distal ~6.5-7).

2. **Weak-acid inhibition.**  Only the *undissociated* acid HA crosses the
   membrane; inside, at near-neutral cytoplasmic pH, it dissociates, dumping
   protons and dissipating the proton-motive force.  The undissociated fraction
   follows Henderson-Hasselbalch, ``HA = A_total / (1 + 10^(pH - pKa))``, so the
   same acid is far more toxic at low pH.  Growth is scaled by

       f_acid = 1 / (1 + sum_i HA_i / Ki)            (per-species Ki)

   and by a cardinal-pH window (Rosso et al. 1995) ``f_pH`` so that growth also
   falls off outside a species' viable pH range.  The species growth rate from
   the FBA solve is multiplied by ``f_acid * f_pH``.

This is what makes butyrate/acetate accumulation *self-limiting* and lets an
acid-sensitive pathogen (e.g. *C. difficile*, given a low ``Ki``) be suppressed
by a dense SCFA-producing community -- the metabolic arm of colonization
resistance that pure nutrient competition under-represents.

References: Russell (1992) J. Appl. Bacteriol.; Rosso et al. (1995) J. Theor.
Biol.; Cherrington et al. (1991) on undissociated SCFA antibacterial action.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Tuple

from muode.ecology import EcologyLayer

#: pKa of common fermentation acids (BiGG extracellular ids).
SCFA_PKA: Dict[str, float] = {
    "ac_e": 4.76,        # acetate
    "ppa_e": 4.87,       # propionate
    "but_e": 4.82,       # butyrate
    "lac__L_e": 3.86,    # L-lactate
    "lac__D_e": 3.86,    # D-lactate
    "for_e": 3.75,       # formate
    "succ_e": 4.21,      # succinate (1st pKa)
}


@dataclass
class WeakAcidInhibition(EcologyLayer):
    """SCFA-driven acidification and weak-acid growth inhibition.

    Parameters
    ----------
    ph0, buffer_capacity, ph_floor:
        Titration model ``pH = max(ph_floor, ph0 - total_acid / buffer_capacity)``
        (mmol/L per pH unit for the buffer).
    default_ki, ki:
        Undissociated-acid sensitivity (mmol/L) giving half-maximal inhibition;
        ``ki`` is a per-species override (lower = more acid-sensitive).
    ph_opt, ph_min, ph_max:
        Cardinal-pH growth window applied to all species.
    acids:
        Which pool metabolites are treated as weak acids.
    """

    name: str = "weak_acid_inhibition"
    ph0: float = 6.8
    buffer_capacity: float = 60.0
    ph_floor: float = 4.3
    default_ki: float = 10.0
    ki: Dict[str, float] = field(default_factory=dict)
    ph_opt: float = 6.8       # default = ph0 so a clean medium is at the optimum
    ph_min: float = 4.3
    ph_max: float = 8.5
    pka: Dict[str, float] = field(default_factory=lambda: dict(SCFA_PKA))
    acids: Tuple[str, ...] = ("ac_e", "ppa_e", "but_e", "lac__L_e", "lac__D_e", "succ_e")

    def total_acid(self, M: Mapping[str, float]) -> float:
        return sum(max(0.0, M.get(a, 0.0)) for a in self.acids)

    def ph(self, M: Mapping[str, float]) -> float:
        return max(self.ph_floor, self.ph0 - self.total_acid(M) / self.buffer_capacity)

    def undissociated(self, M: Mapping[str, float], pH: float) -> float:
        tot = 0.0
        for a in self.acids:
            conc = M.get(a, 0.0)
            if conc <= 0.0:
                continue
            tot += conc / (1.0 + 10.0 ** (pH - self.pka.get(a, 4.8)))
        return tot

    def _cardinal_ph(self, pH: float) -> float:
        """Rosso et al. (1995) cardinal-pH model, in [0, 1].

        In the viable window ``ph_min < pH < ph_max`` both the numerator
        ``(pH-min)(pH-max)`` and the denominator are negative, so their ratio is
        positive and equals 1 at ``ph_opt``.
        """
        if pH <= self.ph_min or pH >= self.ph_max:
            return 0.0
        num = (pH - self.ph_min) * (pH - self.ph_max)
        den = num - (pH - self.ph_opt) ** 2
        if den == 0.0:
            return 0.0
        return max(0.0, min(1.0, num / den))

    def growth_factor(self, organism_id, t, M, X) -> float:
        pH = self.ph(M)
        ha = self.undissociated(M, pH)
        ki = self.ki.get(organism_id, self.default_ki)
        f_acid = 1.0 / (1.0 + ha / ki) if ki > 0 else 1.0
        return f_acid * self._cardinal_ph(pH)

    def observables(self, t, M, X) -> Dict[str, float]:
        pH = self.ph(M)
        return {"pH": pH, "undissociated_acid": self.undissociated(M, pH)}
