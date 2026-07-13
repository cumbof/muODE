"""Diet / growth-medium definitions.

A :class:`Diet` sets the *extracellular* boundary conditions of a simulation:
the initial concentration of every available nutrient, optionally a constant
influx (for an open/chemostat system), and optionally a **maximum uptake flux**
per nutrient.

Units: concentrations in mmol/L; influx in mmol/L/h; ``max_uptake`` in
mmol/gDW/h.

Why ``max_uptake`` exists
-------------------------
Without it, the uptake bound the engine hands a model is purely
``Vmax * C / (Km + C)`` -- and with a *uniform* default Vmax, every nutrient a
diet supplies saturates at the same bound.  A genome-scale model presented with
an 89-metabolite diet then imports 89 nutrients at full tilt simultaneously,
which is not a medium any cell has ever been in.  In a real run this produced
growth rates up to 6.5/h (a 6-minute doubling; the fastest organism ever measured
manages ~10 min, aerobically, in rich medium).

Published diets solve this by specifying a **flux bound per metabolite**, derived
from actual dietary intake: in the VMH/AGORA western gut diet the median bound is
0.1 mmol/gDW/h, not 10.  ``max_uptake`` carries those bounds, and
:func:`muode.media.diet_medium` applies them as a *cap* on the kinetic rate, so
uptake is ``min(Michaelis-Menten, dietary availability)`` -- a cell can be slower
than the diet allows (low affinity, low Vmax) but never faster.

A diet with no ``max_uptake`` behaves exactly as before, so this is additive.

Realistic diets are large curated tables; load them from CSV with
:meth:`Diet.from_csv`.  The small built-in presets here are intended for the toy
example and for smoke-testing the engine, not as biologically complete media.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, Mapping, Optional

_DIET_DIR = Path(__file__).parent / "data" / "diets"


def _strip_comments(lines: Iterable[str]) -> Iterator[str]:
    """Drop blank lines and ``#`` comments so a diet CSV can be documented."""
    for line in lines:
        if line.strip() and not line.lstrip().startswith("#"):
            yield line


@dataclass
class Diet:
    """Extracellular medium for a community simulation.

    Parameters
    ----------
    concentrations:
        Initial concentration (mmol/L) per *extracellular metabolite id*.
    influx:
        Optional constant replenishment rate (mmol/L/h) per metabolite.
    max_uptake:
        Optional ceiling on uptake flux (mmol/gDW/h) per metabolite -- how much of
        this nutrient the diet actually makes available to a cell.  Where given, it
        caps the kinetic uptake rate; where absent, the kinetics are unconstrained.
    name:
        Human-readable label.
    """

    concentrations: Dict[str, float] = field(default_factory=dict)
    influx: Dict[str, float] = field(default_factory=dict)
    max_uptake: Dict[str, float] = field(default_factory=dict)
    name: str = "custom"

    def metabolites(self) -> tuple[str, ...]:
        """Every metabolite the diet mentions, including influx-only ones.

        A nutrient can be supplied purely by influx (initial concentration 0, fed
        in continuously).  Listing only ``concentrations`` would hide it from the
        engine's environment and from the medium exporters, so it would silently
        never be provided.
        """
        return tuple(self.concentrations) + tuple(
            m for m in self.influx if m not in self.concentrations
        )

    def initial_concentration(self, metabolite_id: str) -> float:
        return float(self.concentrations.get(metabolite_id, 0.0))

    def influx_rate(self, metabolite_id: str) -> float:
        return float(self.influx.get(metabolite_id, 0.0))

    def uptake_limit(self, metabolite_id: str) -> Optional[float]:
        """Dietary ceiling on uptake flux (mmol/gDW/h), or None if unconstrained."""
        limit = self.max_uptake.get(metabolite_id)
        return None if limit is None else float(limit)

    def with_metabolites(self, extra: Mapping[str, float]) -> "Diet":
        """Return a copy ensuring ``extra`` metabolites exist (default conc 0)."""
        merged = dict(self.concentrations)
        for m, c in extra.items():
            merged.setdefault(m, c)
        return Diet(merged, dict(self.influx), dict(self.max_uptake), self.name)

    # -- IO -----------------------------------------------------------------
    @classmethod
    def from_csv(cls, path: str | Path, name: Optional[str] = None) -> "Diet":
        """Load a diet from CSV.

        Expected columns (header, case-insensitive): ``metabolite``,
        ``concentration`` and optionally ``influx`` and ``max_uptake``.  This is
        intentionally permissive so that VMH/MICOM-style diet exports can be
        adapted with a light column rename.

        An empty ``max_uptake`` cell means *unconstrained*, which is not the same
        as ``0`` (not supplied at all) -- so blanks are skipped rather than
        coerced.
        """
        import csv

        path = Path(path)
        conc: Dict[str, float] = {}
        influx: Dict[str, float] = {}
        max_uptake: Dict[str, float] = {}
        with path.open(newline="") as fh:
            reader = csv.DictReader(_strip_comments(fh))
            fields = {f.lower(): f for f in (reader.fieldnames or [])}
            met_col = fields.get("metabolite") or fields.get("reaction") or fields.get("id")
            conc_col = fields.get("concentration") or fields.get("flux") or fields.get("amount")
            influx_col = fields.get("influx")
            uptake_col = fields.get("max_uptake") or fields.get("max_flux")
            if met_col is None or conc_col is None:
                raise ValueError(
                    f"{path}: need 'metabolite' and 'concentration' columns, got {reader.fieldnames}"
                )
            for row in reader:
                met = row[met_col].strip()
                if not met:
                    continue
                conc[met] = float(row[conc_col])
                if influx_col and row.get(influx_col):
                    influx[met] = float(row[influx_col])
                if uptake_col and (row.get(uptake_col) or "").strip():
                    max_uptake[met] = float(row[uptake_col])
        return cls(conc, influx, max_uptake, name or path.stem)

    def to_csv(self, path: str | Path) -> None:
        import csv

        with Path(path).open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["metabolite", "concentration", "influx", "max_uptake"])
            for met, c in sorted(self.concentrations.items()):
                limit = self.max_uptake.get(met)
                writer.writerow([met, c, self.influx.get(met, 0.0),
                                 "" if limit is None else limit])


# ---------------------------------------------------------------------------
# Built-in presets (demo / smoke-test scale only)
# ---------------------------------------------------------------------------


def _glucose_minimal() -> Diet:
    """A minimal medium with a single fermentable sugar (toy example)."""
    return Diet(
        concentrations={"glc_e": 20.0, "ac_e": 0.0},
        name="glucose_minimal",
    )


def _western_gut_demo() -> Diet:
    """Tiny stand-in for a Western gut diet -- TOY MEDIUM, sugars only.

    Six metabolites, no nitrogen/phosphate/sulfur/ion source, so a *genome-scale*
    biomass reaction cannot fire on it: real GEMs will report zero growth.  It is
    here for the toy `LinprogOrganism` examples and engine smoke-tests only.  For
    real GEMs use the ``western_gut`` preset (:func:`_western_gut`), which is a
    complete medium.
    """
    return Diet(
        concentrations={
            "glc__D_e": 10.0,  # glucose (BiGG id -- NOT `glc_e`, which matches nothing)
            "fru_e": 5.0,      # fructose
            "lcts_e": 2.0,     # lactose
            "ac_e": 0.0,       # acetate (cross-fed)
            "but_e": 0.0,      # butyrate
            "ppa_e": 0.0,      # propionate
        },
        influx={"glc__D_e": 1.0, "fru_e": 0.5},
        name="western_gut_demo",
    )


def _western_gut() -> Diet:
    """Western-diet colonic medium in the BiGG namespace (CarveMe-compatible).

    **Published, not curated.**  Derived from the VMH/AGORA western gut diet as
    mapped to BiGG by the MICOM media collection (Diener et al., mSystems 2020);
    regenerate with ``examples/diets/derive_western_gut.py``.  It carries
    ``max_uptake`` -- the dietary flux bound per metabolite, median 0.1 mmol/gDW/h
    -- which is what stops a genome-scale model importing every nutrient at once
    and "growing" faster than any organism alive.

    Two properties surprise people, and both are correct colonic biology:

    * **No free glucose.**  Glucose is absorbed in the small intestine; what
      reaches the colon is starch, amylose, pullulan, lactose and fibre.  A model
      that cannot degrade a polysaccharide therefore cannot grow on this medium
      alone -- in the real gut it lives on sugars cross-fed by primary degraders,
      and in a muODE simulation it must do the same.
    * **Microaerobic, not anaerobic.**  o2_e is supplied at 0.001 mmol/gDW/h (100x
      *below* the median bound) -- the mucosal oxygen gradient.  Far too little to
      support aerobic growth, but it is not zero.
    """
    return Diet.from_csv(_DIET_DIR / "western_gut.csv", name="western_gut")


def _dm38() -> Diet:
    """DM38 -- the chemically defined medium of Clark et al. 2021 (Nat Commun).

    Unlike :func:`_western_gut`, this is a *published* medium with an exact
    composition, and it is the medium in which 1,850 synthetic gut communities
    were actually measured.  It is therefore the medium to use whenever the point
    is to compare a prediction against those measurements
    (:mod:`muode.benchmarks.clark2021`).

    Anaerobic batch culture: no oxygen, and no influx of anything.  Note that it
    supplies 28.3 mM L-lactate, so lactate is a substrate here as well as a
    fermentation product.
    """
    return Diet.from_csv(_DIET_DIR / "dm38.csv", name="DM38")


_PRESETS = {
    "glucose_minimal": _glucose_minimal,
    "western_gut": _western_gut,
    "western_gut_demo": _western_gut_demo,
    "dm38": _dm38,
}


def load_preset(name: str) -> Diet:
    """Return a built-in diet preset by name."""
    if name not in _PRESETS:
        raise KeyError(f"unknown diet preset '{name}'; available: {sorted(_PRESETS)}")
    return _PRESETS[name]()


def load_diet(spec: str | Path) -> Diet:
    """Resolve a diet given as a preset name *or* a path to a CSV."""
    text = str(spec)
    if text.endswith(".csv"):
        path = Path(text)
        if not path.exists():
            raise FileNotFoundError(f"diet CSV not found: {path}")
        return Diet.from_csv(path)
    return load_preset(text)


def available_presets() -> tuple[str, ...]:
    return tuple(sorted(_PRESETS))
