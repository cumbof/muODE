"""Worked mechanistic scenarios built from dependency-light toy models.

These tie the :mod:`muode.ecology` layers together into complete, runnable
stories that need only numpy/scipy -- no GEMs, no CarveMe -- so the *mechanisms*
can be exercised and tested on any machine. The genome-scale version of the same
story is the ``examples/fmt_cdiff`` workflow.

Two levels of detail are provided:

* :mod:`muode.scenarios.cdi` -- the minimal recurrent-CDI / FMT sketch (a
  3-member community) and a 2-member phage-predation scenario. Its public API
  (:func:`cdi_scenario`, :func:`phage_predation_scenario` and the role
  constants) is re-exported here for convenience and backwards compatibility.
* :mod:`muode.scenarios.rcdi` -- a richer, guild-structured 12-member community
  and the full ecology stack behind the ``examples/fmt_cdiff/dynamics`` study.
"""

from __future__ import annotations

from muode.scenarios import rcdi  # noqa: F401  (rich 12-guild rCDI scenario module)
from muode.scenarios.cdi import (
    BAI,
    CDIFF,
    COMMENSAL,
    COMPETITOR,
    PATHOGEN,
    cdi_community,
    cdi_diet,
    cdi_ecology,
    cdi_kinetics,
    cdi_scenario,
    phage_predation_scenario,
)

__all__ = [
    # minimal CDI / phage sketch (muode.scenarios.cdi)
    "CDIFF",
    "COMPETITOR",
    "BAI",
    "PATHOGEN",
    "COMMENSAL",
    "cdi_kinetics",
    "cdi_diet",
    "cdi_community",
    "cdi_ecology",
    "cdi_scenario",
    "phage_predation_scenario",
    # richer guild-structured rCDI scenario builder (muode.scenarios.rcdi)
    "rcdi",
]
