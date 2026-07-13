"""Calibrating uptake kinetics against measurements -- and knowing when you cannot.

muODE's Vmax values are invented (a uniform default, jittered by a hash).  The
obvious remedy is to fit them to Clark et al.'s 1,850 measured communities.  This
module does that, and -- more importantly -- refuses to when the data will not
support it.

Why the refusal matters
-----------------------
Clark's measurements are *endpoints*: one OD and one set of HPLC concentrations
after 48 h of sealed batch culture.  In a batch culture that runs to substrate
exhaustion, the endpoint is determined by **yield** (how much biomass and product
the stoichiometry gets per mole of substrate), not by **rate** (how fast the
substrate was taken up).  Doubling Vmax makes the culture finish sooner; it does
not make it finish differently.

Measured on E. coli core in DM38, endpoint acetate moves from 22.4 to 20.4 mM
while Vmax ranges over 6 -> 100 mmol/gDW/h.  That is a 16-fold change in the
parameter producing a 9% change in the observable, most of it integration noise.
Any optimiser handed that objective will report a confident minimum somewhere on
the plateau, and it will be fitting nothing.

What the data *does* determine is a **lower bound**: Vmax must be large enough to
exhaust the substrate within the incubation.  Below that the culture is still
growing when the experiment ends (and far below it, the cell cannot even pay its
ATP maintenance and dies).  So:

* :func:`profile` sweeps Vmax and records what the model predicts at each value.
* :func:`assess` classifies the result -- ``identified``, ``bounded_below``,
  ``unidentifiable`` or ``no_growth`` -- and returns a point estimate *only* in the
  first case.
* :func:`fit_vmax` is the two together.

The corollary is good news for the benchmark: if the predictions are insensitive to
Vmax above the bound, then muODE's Clark endpoint scores test the *reconstructions*
and their stoichiometry, and the invented kinetics do not contaminate them.  Where
kinetics genuinely bite is competition (who wins a shared substrate, i.e. community
composition) and any time-resolved prediction -- neither of which an endpoint
measures.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

#: Below this relative spread, a prediction is flat: the parameter does not move it.
FLAT_TOL = 0.10

#: Biomass must exceed the inoculum by this factor to count as "grew at all".
GROWTH_FACTOR = 1.5


@dataclass(frozen=True)
class Profile:
    """What the model predicts across a sweep of Vmax."""

    vmax: Tuple[float, ...]
    #: target name -> predicted value at each Vmax
    predicted: Dict[str, Tuple[float, ...]]
    #: objective (weighted SSE against the observations) at each Vmax; None if no
    #: observations were supplied
    objective: Optional[Tuple[float, ...]] = None

    def spread(self, target: str) -> float:
        """Relative spread of a prediction over the *growing* part of the sweep."""
        values = [v for v in self.predicted[target] if v is not None]
        if not values:
            return 0.0
        lo, hi = min(values), max(values)
        scale = max(abs(hi), abs(lo))
        return 0.0 if scale <= 0 else (hi - lo) / scale


@dataclass(frozen=True)
class Identifiability:
    """The honest verdict on whether the data pins Vmax down.

    ``verdict`` is one of:

    ``no_growth``
        The model never grew, at any Vmax.  Nothing can be fitted; the model or the
        medium is at fault, not the kinetics.
    ``identified``
        Predictions genuinely vary with Vmax and the objective has a minimum inside
        the swept range.  ``best`` is meaningful.
    ``bounded_below``
        The model needs a minimum Vmax to grow / finish, and above that the
        predictions plateau.  The data gives a bound, not a value.  ``best`` is None
        *on purpose*.
    ``unidentifiable``
        Predictions are flat across the whole sweep.  Vmax does not matter here.
    """

    verdict: str
    best: Optional[float]
    lower_bound: Optional[float]
    plateau_from: Optional[float]
    #: largest relative spread of any prediction across the growing range
    sensitivity: float

    @property
    def fitted(self) -> bool:
        return self.verdict == "identified"

    def summary(self) -> str:
        if self.verdict == "no_growth":
            return "no growth at any Vmax -- fix the model or the medium, not the kinetics"
        if self.verdict == "identified":
            return f"Vmax identified: {self.best:g} mmol/gDW/h"
        if self.verdict == "bounded_below":
            return (
                f"Vmax NOT identified -- only bounded: >= {self.lower_bound:g} mmol/gDW/h. "
                f"Above {self.plateau_from:g} the predictions plateau "
                f"({self.sensitivity:.0%} spread), so the endpoint carries no rate "
                f"information. Reporting a point estimate here would be fitting noise."
            )
        return (
            f"Vmax unidentifiable: predictions vary by only {self.sensitivity:.0%} across "
            "the whole sweep. The measurement does not constrain the rate."
        )


def profile(
    predict: Callable[[float], Mapping[str, float]],
    vmax_grid: Sequence[float],
    observed: Optional[Mapping[str, float]] = None,
    weights: Optional[Mapping[str, float]] = None,
) -> Profile:
    """Sweep Vmax and record the prediction at each value.

    ``predict(vmax)`` returns a mapping of target -> predicted value (e.g.
    ``{"but_e": 3.1, "biomass": 0.42}``).  ``observed`` enables the objective.
    """
    grid = tuple(float(v) for v in vmax_grid)
    rows = [dict(predict(v)) for v in grid]
    targets = sorted({k for r in rows for k in r})
    predicted = {t: tuple(r.get(t) for r in rows) for t in targets}

    objective: Optional[Tuple[float, ...]] = None
    if observed:
        w = weights or {}
        obj = []
        for r in rows:
            sse = 0.0
            for t, o in observed.items():
                p = r.get(t)
                if p is None:
                    continue
                # scale-free: an error of 1 mM in butyrate (~3 mM) must not be
                # swamped by an error of 1 mM in acetate (~30 mM)
                scale = max(abs(float(o)), 1e-6)
                sse += w.get(t, 1.0) * ((float(p) - float(o)) / scale) ** 2
            obj.append(sse)
        objective = tuple(obj)

    return Profile(vmax=grid, predicted=predicted, objective=objective)


def assess(
    prof: Profile,
    biomass_key: str = "biomass",
    inoculum: Optional[float] = None,
    flat_tol: float = FLAT_TOL,
) -> Identifiability:
    """Decide whether the sweep actually pins Vmax down.  See :class:`Identifiability`."""
    grid = prof.vmax
    growth = prof.predicted.get(biomass_key)

    # which Vmax values produced a culture that grew at all?
    if growth is not None and inoculum is not None:
        alive = [
            i for i, g in enumerate(growth)
            if g is not None and g > inoculum * GROWTH_FACTOR
        ]
    else:
        alive = list(range(len(grid)))
    if not alive:
        return Identifiability("no_growth", None, None, None, 0.0)

    lower_bound = grid[alive[0]]
    targets = [t for t in prof.predicted if t != biomass_key]

    def _spread(idxs) -> float:
        """Largest relative spread of any prediction over these grid points."""
        worst = 0.0
        for t in targets:
            values = [prof.predicted[t][i] for i in idxs if prof.predicted[t][i] is not None]
            if len(values) < 2:
                continue
            lo, hi = min(values), max(values)
            scale = max(abs(hi), abs(lo))
            if scale > 0:
                worst = max(worst, (hi - lo) / scale)
        return worst

    # Find where the plateau starts: the earliest point from which *every* prediction
    # stops moving.  This must be a suffix, not a global spread.  A culture that is
    # alive but has not yet exhausted its substrate is still climbing -- its low
    # endpoint is incompleteness, not sensitivity to the rate -- so measuring spread
    # across the whole live range would mistake "still running" for "the parameter
    # matters" and hand back a fitted number.
    plateau_k: Optional[int] = None
    for k in range(len(alive) - 1):
        if _spread(alive[k:]) < flat_tol:
            plateau_k = k
            break

    if plateau_k is None:                      # never settles: the rate does move it
        sensitivity = _spread(alive)
        if prof.objective is None:
            return Identifiability("bounded_below", None, lower_bound, None, sensitivity)
        best_i = min(alive, key=lambda i: prof.objective[i])
        if best_i in (alive[0], alive[-1]):    # optimum at the edge is not an optimum
            return Identifiability("bounded_below", None, lower_bound, None, sensitivity)
        return Identifiability("identified", grid[best_i], lower_bound, None, sensitivity)

    plateau_from = grid[alive[plateau_k]]
    sensitivity = _spread(alive[plateau_k:])   # how flat the plateau really is

    if plateau_k == 0 and len(alive) == len(grid):
        # flat across the entire sweep, and it grew everywhere: the parameter is inert
        return Identifiability("unidentifiable", None, None, None, sensitivity)

    if prof.objective is not None:
        best_i = min(alive, key=lambda i: prof.objective[i])
        # The optimum must sit BELOW the plateau (where predictions genuinely move)
        # *and* strictly inside the swept range.  An optimum on the lowest grid point
        # is not an optimum -- the real one may lie below it, and we never looked.
        if alive[0] < best_i < alive[plateau_k]:
            return Identifiability("identified", grid[best_i], lower_bound, plateau_from,
                                   sensitivity)

    # the best fit lies on the plateau -- i.e. anywhere up there fits equally well
    return Identifiability("bounded_below", None, lower_bound, plateau_from, sensitivity)


def fit_vmax(
    predict: Callable[[float], Mapping[str, float]],
    observed: Mapping[str, float],
    vmax_grid: Sequence[float],
    inoculum: Optional[float] = None,
    biomass_key: str = "biomass",
    weights: Optional[Mapping[str, float]] = None,
    flat_tol: float = FLAT_TOL,
) -> Tuple[Identifiability, Profile]:
    """Profile Vmax against measurements and report what the data actually supports."""
    prof = profile(predict, vmax_grid, observed=observed, weights=weights)
    return assess(prof, biomass_key=biomass_key, inoculum=inoculum, flat_tol=flat_tol), prof
