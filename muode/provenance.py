"""Where did that number come from?

muODE's medium already answers this.  Every row of a diet carries a ``source`` column
-- ``intake`` (a real dietary figure), ``dri`` (a reference intake), ``derived`` (our
arithmetic, written out), ``source_default`` (someone else's fill value, which is *not*
a measurement and says so).  That column is what let an audit find that vitamin B12 was
the largest carbon source in the medium, and that nine nucleoside rows were feeding
E. coli an RNA diet nobody eats.

**The kinetic and ecology parameters have no such column, and they need one more.**

    bile.py         vmax_bsh=5.0, km_bsh=0.05, vmax_bai=1.0, km_bai=0.05, ki=0.02
    antibiotic.py   dose=1.0, half_life=2.0, emax=2.0, ec50=0.5
    antagonism.py   production=0.5, decay=0.2, ki=0.1
    kinetics.py     DEFAULT_VMAX=10.0, DEFAULT_KM=0.01

Not one of those is measured.  They are bare floats with no units at the call site, no
citation, and nothing that would stop one reaching a figure in a paper.  The rCDI/FMT
prediction -- the headline result -- rests on the first two rows entirely.

The distinction that matters
----------------------------
``ecology.py`` says the layers are "real, literature-grounded mechanisms, not fudge
factors", and **that is true of the mechanisms and false of the numbers**.  The BSH ->
cholate -> deoxycholate pathway is real chemistry with real citations; the *rate* at
which we run it is a number somebody typed.  Those are separable, and only the second is
a problem.  This module keeps them separable, and makes the second one impossible to
ignore.

How this is meant to be used
----------------------------
1. **Declare** every parameter with an :class:`Evidence` class and a reason.
   :data:`Evidence.INVENTED` is not a slur, it is a status: it means *nobody chose this
   number*, and it is the honest label for most of what is listed above.
2. **Report.**  :func:`report` renders the whole registry as a table -- the parameter
   table of a methods section, generated rather than remembered.
3. **Compute what depends on it.**  An invented parameter is only fatal if a conclusion
   *moves* when you change it.  :func:`sensitivity` re-runs a scenario across a range of
   values and reports whether the answer holds.  This is the same trick as
   :func:`muode.qc.nutrient_sensitivity`, for the same reason: you cannot fix 20 numbers
   by reading 20 papers, but you can find out which ones the result actually turns on.

**A conclusion that is INSENSITIVE to an invented parameter is publishable, and the
sensitivity sweep is the evidence.  A conclusion that is SENSITIVE to one is not, and no
amount of prose fixes that -- only a measurement does.**
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from typing import Callable, Dict, Iterable, List, Optional, Sequence


class Evidence(str, Enum):
    """How much do we actually know about this number?

    Ordered from strongest to weakest.  The line that matters for publication sits
    between ASSUMED and INVENTED: an assumption is a choice you can defend, an invention
    is a number nobody chose.
    """

    #: A specific published measurement of this quantity, in a comparable system.
    #: `citation` is mandatory and must be a real one.
    MEASURED = "measured"

    #: Computed from measured quantities, with the arithmetic written down.  The diet's
    #: xylan and urea bounds are the model: a g/day intake and a molecular weight, and a
    #: conversion you can check.
    DERIVED = "derived"

    #: A deliberate modelling choice, defensible and stated, but NOT a measurement.  A
    #: residence time of 24 h for the colon is an assumption; so is splitting one niacin
    #: RDA across two vitamers.  Legitimate -- as long as it says so.
    ASSUMED = "assumed"

    #: **Nobody chose this number.**  It is a placeholder that made the code run.
    #: It is not a measurement, not a derivation, and not a defensible assumption.
    #: A result that depends on one of these is not a result.
    INVENTED = "invented"

    @property
    def is_evidence(self) -> bool:
        return self in (Evidence.MEASURED, Evidence.DERIVED)


@dataclass(frozen=True)
class Parameter:
    """One number, and everything we know about where it came from."""

    name: str
    value: float
    units: str
    evidence: Evidence
    why: str
    citation: str = ""

    def __post_init__(self) -> None:
        if not self.why.strip():
            raise ValueError(f"{self.name}: a parameter must say where its value came from")
        if self.evidence is Evidence.MEASURED and not self.citation.strip():
            raise ValueError(
                f"{self.name}: claims to be MEASURED but cites nothing.  A measurement "
                "without a citation is an assertion."
            )

    @property
    def publishable(self) -> bool:
        """False if a conclusion resting on this number would be indefensible."""
        return self.evidence is not Evidence.INVENTED


#: Every parameter muODE uses, keyed ``module.parameter``.
_REGISTRY: Dict[str, Parameter] = {}


def register(module: str, param: Parameter) -> Parameter:
    key = f"{module}.{param.name}"
    if key in _REGISTRY and _REGISTRY[key] != param:
        raise ValueError(f"{key} is already registered with a different value")
    _REGISTRY[key] = param
    return param


def registry() -> Dict[str, Parameter]:
    return dict(_REGISTRY)


def invented() -> List[str]:
    """Every parameter nobody chose.  This list should only ever get shorter."""
    return sorted(k for k, p in _REGISTRY.items() if p.evidence is Evidence.INVENTED)


def unpublishable(keys: Optional[Iterable[str]] = None) -> List[str]:
    """Of the parameters in play, which would sink a paper?

    Pass the keys a given scenario actually uses; omit for the whole framework.
    """
    ks = sorted(keys) if keys is not None else sorted(_REGISTRY)
    return [k for k in ks if k in _REGISTRY and not _REGISTRY[k].publishable]


def report(keys: Optional[Iterable[str]] = None) -> str:
    """The parameter table of a methods section, generated rather than remembered."""
    ks = sorted(keys) if keys is not None else sorted(_REGISTRY)
    rows = [(k, _REGISTRY[k]) for k in ks if k in _REGISTRY]
    if not rows:
        return "no parameters registered\n"

    w = max(len(k) for k, _ in rows)
    out = [
        f"{'parameter'.ljust(w)}  {'value':>9}  {'units':<12}  {'evidence':<9}  why",
        f"{'-' * w}  {'-' * 9}  {'-' * 12}  {'-' * 9}  {'-' * 40}",
    ]
    for k, p in rows:
        cite = f"  [{p.citation}]" if p.citation else ""
        out.append(
            f"{k.ljust(w)}  {p.value:>9.4g}  {p.units:<12}  "
            f"{p.evidence.value:<9}  {p.why}{cite}"
        )

    bad = [k for k, p in rows if not p.publishable]
    n = len(rows)
    ev = sum(1 for _, p in rows if p.evidence.is_evidence)
    out += [
        "",
        f"{ev}/{n} parameters rest on evidence (measured or derived).",
    ]
    if bad:
        out += [
            "",
            f"!! {len(bad)} INVENTED -- nobody chose these numbers:",
            *(f"     {k}" for k in bad),
            "",
            "   A conclusion that MOVES when one of these changes is not a result.",
            "   Run muode.provenance.sensitivity() to find out whether yours does.",
        ]
    return "\n".join(out) + "\n"


@dataclass(frozen=True)
class Sensitivity:
    """Does the conclusion survive not knowing this number?"""

    parameter: str
    values: Sequence[float]
    outcomes: Sequence[float]
    verdict: str            # "robust" | "sensitive"
    spread: float           # relative spread of the outcome across the sweep

    @property
    def robust(self) -> bool:
        return self.verdict == "robust"


def sensitivity(
    run: Callable[[float], float],
    parameter: str,
    values: Sequence[float],
    tol: float = 0.10,
) -> Sensitivity:
    """Sweep one parameter and ask whether the conclusion holds regardless.

    ``run(value) -> outcome`` re-runs the scenario with the parameter set to ``value``
    and returns the number the claim rests on (a final biomass, a butyrate
    concentration, a cleared-or-not flag).

    An invented parameter is only fatal **if the answer depends on it**.  Sweep it over
    the range you cannot rule out; if the outcome moves by less than ``tol`` in relative
    terms across that whole range, the conclusion is *robust to not knowing the number*,
    and this object is the evidence you cite instead of the measurement you do not have.
    If it moves more, the conclusion is an artefact of a guess, and the only fix is to go
    and measure it.

    This is deliberately the same move as :func:`muode.qc.nutrient_sensitivity`.  You
    cannot fix twenty numbers by reading twenty papers.  You *can* find out which of the
    twenty your result actually turns on -- and usually it is two or three.
    """
    outcomes = [float(run(v)) for v in values]
    lo, hi = min(outcomes), max(outcomes)
    scale = max(abs(lo), abs(hi))
    spread = 0.0 if scale == 0 else (hi - lo) / scale
    return Sensitivity(
        parameter=parameter,
        values=list(values),
        outcomes=outcomes,
        verdict="robust" if spread <= tol else "sensitive",
        spread=spread,
    )


def with_value(layer, field: str, value: float):
    """A copy of a dataclass ecology layer with one field changed -- for sweeping."""
    return replace(layer, **{field: value})
