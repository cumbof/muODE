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
    name:
        Human-readable label.
    """

    concentrations: Dict[str, float] = field(default_factory=dict)
    influx: Dict[str, float] = field(default_factory=dict)
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
            reader = csv.DictReader(_strip_comments(fh))
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
    """Anaerobic Western-style gut medium in the BiGG namespace (CarveMe-compatible).

    A *complete* medium: carbon sources plus the nitrogen, phosphate, sulfur,
    ions, trace metals, amino acids, nucleobases and vitamins a genome-scale
    biomass reaction needs.  Curated and literature-informed -- it is **not** the
    official VMH Western-diet table; export that and use :meth:`Diet.from_csv` if
    you need a published diet.
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
