"""
common.py -- Presentation helpers + a single import surface for the scenario
             scripts in this folder.

The *model* (guild registry, gut medium, ecology stack, timeline and the
``run_scenario`` driver) now lives in the muODE package itself, at
:mod:`muode.scenarios.rcdi`, so this example and the library share one source of
truth. This module simply re-exports that model and adds the plotting style used
by the four scenario scripts (``antibiotic_relapse.py``, ``fmt_resolution.py``,
``designed_consortium.py``, ``failed_fmt_autopsy.py``), which import everything
they need from here as ``import common as C``.

To inspect or tune the model -- guilds, yields, medium, ecology constants,
timeline -- edit ``muode/scenarios/rcdi.py``.
"""

from __future__ import annotations

from pathlib import Path

# Re-export the whole scenario model (metabolite ids, guild registry, medium,
# kinetics, ecology builder, timeline constants, run_scenario, bile/drug series
# helpers, biomass_sum) so the scripts can reach it as C.<name>.
from muode.scenarios.rcdi import *          # noqa: F401,F403  (model + scenario API)
from muode.scenarios.rcdi import ANTIBIOTIC_DAYS


# ===========================================================================
# PLOTTING  (a small, consistent, colour-blind-friendly style)
# ===========================================================================
# All figures share one look. Colours are drawn from the Okabe-Ito palette
# (colour-blind safe). Matplotlib is imported lazily so importing this module on
# a headless node never requires a display.

# Okabe-Ito qualitative palette.
PALETTE = {
    "pathogen": "#D55E00",      # vermillion  -- C. difficile
    "commensal": "#0072B2",     # blue        -- donor / commensal biomass
    "primary": "#E69F00",       # orange      -- primary bile acids
    "secondary": "#009E73",     # bluish green -- secondary bile acids
    "drug": "#555555",          # grey        -- antibiotic concentration
    "accent": "#CC79A7",        # reddish purple
    "muted": "#999999",
}


def _mpl():
    """Return a configured pyplot (Agg backend, shared rcParams)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 150,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
        "axes.titlesize": 11,
    })
    return plt


def output_dir() -> Path:
    """Create and return the results directory (next to these scripts)."""
    out = Path(__file__).resolve().parent / "results"
    out.mkdir(parents=True, exist_ok=True)
    return out


def shade_dosing(ax, days: float = ANTIBIOTIC_DAYS) -> None:
    """Shade the antibiotic dosing window on a time-axis plot."""
    ax.axvspan(0.0, days, color=PALETTE["drug"], alpha=0.08, lw=0,
               label="_dosing")
