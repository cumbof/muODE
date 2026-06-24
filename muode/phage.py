"""Bacteriophage predation -- a coupled infection ODE, not an FBA organism.

A phage has no metabolism: no biomass objective, no exchange fluxes, nothing for
Flux Balance Analysis to optimise.  Representing it as an :class:`OrganismModel`
would be a category error.  What a phage *does* is **population dynamics**: it
adsorbs to a susceptible host, the infected cell stops dividing, and after a
latent period it lyses, releasing a burst of new virions.  That is the
Levin-Stewart / Lotka-Volterra paradigm, and it couples cleanly to the metabolic
ODEs because the only thing it changes in the FBA world is **host biomass**.

This layer therefore lives in the ecology stack (:mod:`muode.ecology`) and uses
the same ``integrate`` hook that :class:`~muode.lifecycle.SporeForming` uses to
move biomass between pools.  It maintains two internal state variables:

* ``P`` -- free phage titer (arbitrary but self-consistent units, e.g. PFU/mL),
* ``I`` -- infected host biomass (gDW/L), which neither grows nor is counted as
  live host until it lyses.

Per Euler step of size ``dt`` (one susceptible host species per layer)::

    adsorption   a = k_ads * P * X_host                 (capped at X_host/dt)
    X_host  -=  (1 - lysogeny) * a * dt                  # lytically infected leave X
    I       +=  (1 - lysogeny) * a * dt  -  lysis * dt   # lysis = I / latent_period
    P       +=  burst * lysis * dt  -  a * dt  -  decay * P * dt

The lytic loop amplifies phage and crashes a host bloom; ``lysogeny_fraction``
diverts that fraction of adsorptions into surviving carriers (a coarse temperate
model) so the predation is weaker.  What is deliberately **not** modelled here --
because it belongs to a different paradigm and would be faked otherwise -- is the
*evolution* of phage host-range or of host resistance (see ``docs/LIMITATIONS.md``
section 4); the host range is fixed for the run.

Reference: Levin, Stewart & Chao (1977) *Am. Nat.*; Weitz (2016) *Quantitative
Viral Ecology*.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping

from muode.ecology import EcologyLayer


@dataclass
class PhageInfection(EcologyLayer):
    """Lytic (optionally temperate) phage preying on one host species.

    Parameters
    ----------
    host:
        organism id of the susceptible host (one host per layer; use several
        layers, each with its own phage pool, for a multi-host or cocktail case).
    adsorption_rate:
        infection rate constant ``k_ads`` ((gDW/L)^-1 h^-1): the per-host,
        per-titer probability of a productive adsorption.
    burst_size:
        effective virions released per unit host biomass lysed (lumps the
        per-cell burst with the cell-mass-to-titer conversion).
    latent_period:
        mean time from infection to lysis (h); lysis is first-order at
        ``1/latent_period`` (an exponential-latency approximation of the delay).
    decay_rate:
        first-order free-phage decay (1/h).
    initial_titer:
        free phage present at t=0 -- a resident virome, or set by a dose event.
    lysogeny_fraction:
        fraction of adsorptions that lysogenise (host survives as a carrier)
        instead of lysing; 0 = strictly lytic (default), 1 = no predation.
    name:
        label used for recorded observables (``phage[name]``, ``infected[host]``).
    """

    host: str
    adsorption_rate: float = 5.0
    burst_size: float = 50.0
    latent_period: float = 0.5
    decay_rate: float = 0.1
    initial_titer: float = 0.0
    lysogeny_fraction: float = 0.0
    name: str = "phage"

    _P: float = field(default=0.0, init=False, repr=False)
    _I: float = field(default=0.0, init=False, repr=False)

    # -- EcologyLayer hooks -------------------------------------------------
    def reset(self, community) -> None:
        self._P = float(self.initial_titer)
        self._I = 0.0

    def integrate(self, t, dt, X, M, mu) -> None:
        x = max(0.0, X.get(self.host, 0.0))
        # adsorption flux, capped so a single step cannot infect more host than
        # is present (a CFL-style guard mirroring the engine's uptake cap).
        a = self.adsorption_rate * self._P * x
        if dt > 0:
            a = min(a, x / dt)
        lytic = (1.0 - self.lysogeny_fraction) * a
        lysis = self._I / self.latent_period if self.latent_period > 0 else 0.0

        X[self.host] = max(0.0, x - lytic * dt)
        self._I = max(0.0, self._I + (lytic - lysis) * dt)
        dP = self.burst_size * lysis - a - self.decay_rate * self._P
        self._P = max(0.0, self._P + dP * dt)

    def observables(self, t, M, X) -> Dict[str, float]:
        return {f"phage[{self.name}]": self._P, f"infected[{self.host}]": self._I}

    def describe(self) -> str:
        mode = "temperate" if self.lysogeny_fraction > 0 else "lytic"
        return (f"{self.name}: {mode} phage on {self.host} "
                f"(burst={self.burst_size}, latent={self.latent_period} h)")

    # -- introspection ------------------------------------------------------
    @property
    def titer(self) -> float:
        return self._P

    @property
    def infected_biomass(self) -> float:
        return self._I
