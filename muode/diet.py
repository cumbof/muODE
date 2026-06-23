"""Diet / growth-medium definitions.

A :class:`Diet` sets the *extracellular* boundary conditions of a simulation:
the initial concentration of every available nutrient and, optionally, a
constant influx (for an open/chemostat system) so that staple nutrients are
replenished rather than exhausted within a few hours of simulated time.

Concentrations are in mmol/L; influx rates in mmol/L/h.

Realistic diets (e.g. the VMH "Western" gut diet used by MICOM/AGORA) are large
curated tables keyed by exchange-reaction id.  Those are loaded from CSV with
:meth:`Diet.from_csv`.  The small built-in presets here are intended for the
toy example and for smoke-testing the engine, not as biologically complete
media.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Mapping, Optional


@dataclass
class Diet:
    """Extracellular medium for a community simulation.

    Parameters
    ----------
    concentrations:
        Initial concentration (mmol/L) per *extracellular metabolite id*.
    influx:
        Optional constant replenishment rate (mmol/L/h) per metabolite.
    name:
        Human-readable label.
    """

    concentrations: Dict[str, float] = field(default_factory=dict)
    influx: Dict[str, float] = field(default_factory=dict)
    name: str = "custom"

    def metabolites(self) -> tuple[str, ...]:
        return tuple(self.concentrations)

    def initial_concentration(self, metabolite_id: str) -> float:
        return float(self.concentrations.get(metabolite_id, 0.0))

    def influx_rate(self, metabolite_id: str) -> float:
        return float(self.influx.get(metabolite_id, 0.0))

    def with_metabolites(self, extra: Mapping[str, float]) -> "Diet":
        """Return a copy ensuring ``extra`` metabolites exist (default conc 0)."""
        merged = dict(self.concentrations)
        for m, c in extra.items():
            merged.setdefault(m, c)
        return Diet(merged, dict(self.influx), self.name)

    # -- IO -----------------------------------------------------------------
    @classmethod
    def from_csv(cls, path: str | Path, name: Optional[str] = None) -> "Diet":
        """Load a diet from CSV.

        Expected columns (header, case-insensitive): ``metabolite``,
        ``concentration`` and optionally ``influx``.  This is intentionally
        permissive so that VMH/MICOM-style diet exports can be adapted with a
        light column rename.
        """
        import csv

        path = Path(path)
        conc: Dict[str, float] = {}
        influx: Dict[str, float] = {}
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = {f.lower(): f for f in (reader.fieldnames or [])}
            met_col = fields.get("metabolite") or fields.get("reaction") or fields.get("id")
            conc_col = fields.get("concentration") or fields.get("flux") or fields.get("amount")
            influx_col = fields.get("influx")
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
        return cls(conc, influx, name or path.stem)

    def to_csv(self, path: str | Path) -> None:
        import csv

        with Path(path).open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["metabolite", "concentration", "influx"])
            for met, c in sorted(self.concentrations.items()):
                writer.writerow([met, c, self.influx.get(met, 0.0)])


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
    """Tiny stand-in for a Western gut diet.

    The real Western gut diet is a ~200-metabolite VMH table; load it with
    :meth:`Diet.from_csv`.  This demo keeps just a few representative inputs so
    the engine has something to chew on out of the box.
    """
    return Diet(
        concentrations={
            "glc_e": 10.0,   # glucose
            "fru_e": 5.0,    # fructose
            "lcts_e": 2.0,   # lactose
            "ac_e": 0.0,     # acetate (cross-fed)
            "but_e": 0.0,    # butyrate
            "ppa_e": 0.0,    # propionate
        },
        influx={"glc_e": 1.0, "fru_e": 0.5},
        name="western_gut_demo",
    )


_PRESETS = {
    "glucose_minimal": _glucose_minimal,
    "western_gut": _western_gut_demo,
    "western_gut_demo": _western_gut_demo,
}


def load_preset(name: str) -> Diet:
    """Return a built-in diet preset by name."""
    if name not in _PRESETS:
        raise KeyError(f"unknown diet preset '{name}'; available: {sorted(_PRESETS)}")
    return _PRESETS[name]()


def available_presets() -> tuple[str, ...]:
    return tuple(sorted(_PRESETS))
