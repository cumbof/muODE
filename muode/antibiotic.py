"""Antibiotic pharmacokinetics / pharmacodynamics (time-varying killing).

A :class:`~muode.perturb.Perturbation` applies a *constant* bound change for the
whole run -- fine for a gene knockout, wrong for a drug.  A real antibiotic has a
**concentration that rises with each dose and decays between doses**, and a
**kill rate that saturates** with that concentration.  This layer models both:

* **PK** -- one-compartment, first-order elimination.  After each dose at
  ``dose_times`` the concentration is the superposition of decaying boluses,
  ``C(t) = sum_d dose * exp(-k_el (t - t_d))`` with ``k_el = ln2 / half_life``.
* **PD** -- an Emax kill model on the *susceptible* species' vegetative biomass:
  ``death(t) = Emax * C / (EC50 + C)`` (1/h), added to their death rate.

Crucially the extra death applies only to the vegetative pool ``X``; spores in
:mod:`muode.lifecycle` are held separately and are therefore **not** killed --
which is exactly why an antibiotic course clears the active infection but leaves
a dormant reservoir that can recur.

This turns "antibiotic" from an instantaneous switch into a dosing schedule, so
muODE can contrast *vancomycin alone* (knocks back vegetative *C. difficile*,
spores survive, it returns) with *vancomycin + FMT* (spores stay dormant because
the restored community keeps the bile pool inhibitory).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Mapping, Tuple

from muode.ecology import EcologyLayer


@dataclass
class Antibiotic(EcologyLayer):
    """Time-varying antibiotic killing with one-compartment PK and an Emax PD.

    Parameters
    ----------
    susceptible:
        organism ids the drug kills (its spectrum).
    dose_times, dose, half_life:
        dosing schedule (h), per-dose concentration (arbitrary consistent unit,
        e.g. mg/L) and elimination half-life (h).
    emax, ec50:
        maximum kill rate (1/h) and the concentration at half-maximal killing.
    """

    name: str = "antibiotic"
    susceptible: set = field(default_factory=set)
    dose_times: Tuple[float, ...] = (0.0,)
    dose: float = 1.0
    half_life: float = 2.0
    emax: float = 2.0
    ec50: float = 0.5

    def concentration(self, t: float) -> float:
        k_el = math.log(2.0) / self.half_life if self.half_life > 0 else 0.0
        c = 0.0
        for td in self.dose_times:
            if t >= td:
                c += self.dose * math.exp(-k_el * (t - td))
        return c

    def extra_death(self, organism_id, t, M, X) -> float:
        if organism_id not in self.susceptible:
            return 0.0
        c = self.concentration(t)
        return self.emax * c / (self.ec50 + c) if c > 0 else 0.0

    def observables(self, t, M, X) -> Dict[str, float]:
        return {f"drug[{self.name}]": self.concentration(t)}

    def describe(self) -> str:
        return (f"{self.name}: {len(self.dose_times)} dose(s) of {self.dose} "
                f"(t1/2={self.half_life} h) vs {sorted(self.susceptible)}")
