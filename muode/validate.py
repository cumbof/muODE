"""Validation & benchmarking -- does a simulation reproduce a known community?

A dynamic community model is only trustworthy if it recovers things we already
know about reference communities: the **composition** (relative abundances), the
**short-chain fatty acids** and other key metabolites it should produce, and the
**cross-feeding** structure.  This module compares a :class:`SimulationResult`
against a :class:`BenchmarkExpectation` and returns a :class:`ValidationReport`
with per-component metrics and an overall pass/fail under explicit tolerances.

Metrics
-------
* **Relative abundance** -- mean absolute error over the expected taxa, plus a
  Spearman rank correlation when at least two taxa are compared.
* **Metabolites** -- mean absolute error (mmol/L) and Spearman correlation over
  the expected metabolites (e.g. acetate/butyrate/propionate).
* **Cross-feeding** -- precision / recall / F1 of the predicted producer→
  metabolite→consumer edges against the expected ones.

Only the components an expectation actually specifies are scored and gated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

Edge = Tuple[str, str, str]  # (producer, metabolite, consumer)


@dataclass
class BenchmarkExpectation:
    """Known/expected observables for a reference community simulation."""

    relative_abundances: Dict[str, float] = field(default_factory=dict)
    metabolites: Dict[str, float] = field(default_factory=dict)
    cross_feeding: List[Edge] = field(default_factory=list)
    #: tolerances (component is "passed" at or below/above these)
    abundance_mae: float = 0.1        # max mean abs error of relative abundance
    metabolite_mae: float = 5.0       # max mean abs error (mmol/L)
    min_f1: float = 0.5               # min cross-feeding edge F1
    name: str = "benchmark"

    @classmethod
    def from_dict(cls, data: dict) -> "BenchmarkExpectation":
        edges = [tuple(e) for e in data.get("cross_feeding", [])]
        tol = data.get("tolerances", {})
        return cls(
            relative_abundances={str(k): float(v) for k, v in data.get("relative_abundances", {}).items()},
            metabolites={str(k): float(v) for k, v in data.get("metabolites", {}).items()},
            cross_feeding=edges,  # type: ignore[arg-type]
            abundance_mae=float(tol.get("abundance_mae", 0.1)),
            metabolite_mae=float(tol.get("metabolite_mae", 5.0)),
            min_f1=float(tol.get("min_f1", 0.5)),
            name=str(data.get("name", "benchmark")),
        )

    @classmethod
    def from_file(cls, path: str | Path) -> "BenchmarkExpectation":
        import yaml

        return cls.from_dict(yaml.safe_load(Path(path).read_text()) or {})


@dataclass
class ValidationReport:
    metrics: dict
    passed: bool
    name: str = "benchmark"

    def to_dict(self) -> dict:
        return {"name": self.name, "passed": self.passed, "metrics": self.metrics}

    def to_json(self, path: str | Path) -> None:
        import json

        Path(path).write_text(json.dumps(self.to_dict(), indent=2))


def _spearman(pred: Sequence[float], obs: Sequence[float]) -> Optional[float]:
    if len(pred) < 2:
        return None
    try:
        import warnings

        from scipy.stats import spearmanr

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # constant-input arrays -> undefined corr
            rho = float(spearmanr(pred, obs).correlation)
        return None if np.isnan(rho) else rho
    except Exception:
        return None


def _relative_abundance(biomass_final: Dict[str, float]) -> Dict[str, float]:
    total = sum(max(0.0, v) for v in biomass_final.values())
    if total <= 0:
        return {k: 0.0 for k in biomass_final}
    return {k: max(0.0, v) / total for k, v in biomass_final.items()}


def _edge_scores(predicted: set, expected: set) -> dict:
    tp = len(predicted & expected)
    precision = tp / len(predicted) if predicted else (1.0 if not expected else 0.0)
    recall = tp / len(expected) if expected else 1.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return {
        "precision": precision, "recall": recall, "f1": f1,
        "true_positive": tp, "n_predicted": len(predicted), "n_expected": len(expected),
        "missed": sorted("->".join(e) for e in (expected - predicted)),
        "spurious": sorted("->".join(e) for e in (predicted - expected)),
    }


def compare(
    predicted_abundances: Dict[str, float],
    predicted_metabolites: Dict[str, float],
    predicted_edges: Sequence[Edge],
    expectation: BenchmarkExpectation,
) -> ValidationReport:
    """Score predicted observables against an expectation (pure / testable core)."""
    metrics: dict = {}
    passed = True

    # -- relative abundance --------------------------------------------------
    if expectation.relative_abundances:
        rel = _relative_abundance(predicted_abundances)
        species = list(expectation.relative_abundances)
        pred = [rel.get(s, 0.0) for s in species]
        obs = [expectation.relative_abundances[s] for s in species]
        mae = float(np.mean(np.abs(np.array(pred) - np.array(obs))))
        ok = mae <= expectation.abundance_mae
        metrics["relative_abundance"] = {
            "mae": mae, "tolerance": expectation.abundance_mae, "passed": ok,
            "spearman": _spearman(pred, obs),
            "predicted": {s: round(rel.get(s, 0.0), 4) for s in species},
        }
        passed = passed and ok

    # -- metabolites ---------------------------------------------------------
    if expectation.metabolites:
        mets = list(expectation.metabolites)
        pred = [float(predicted_metabolites.get(m, 0.0)) for m in mets]
        obs = [expectation.metabolites[m] for m in mets]
        mae = float(np.mean(np.abs(np.array(pred) - np.array(obs))))
        ok = mae <= expectation.metabolite_mae
        metrics["metabolites"] = {
            "mae": mae, "tolerance": expectation.metabolite_mae, "passed": ok,
            "spearman": _spearman(pred, obs),
            "predicted": {m: round(float(predicted_metabolites.get(m, 0.0)), 4) for m in mets},
        }
        passed = passed and ok

    # -- cross-feeding edges -------------------------------------------------
    if expectation.cross_feeding:
        scores = _edge_scores({tuple(e) for e in predicted_edges},
                              {tuple(e) for e in expectation.cross_feeding})
        ok = scores["f1"] >= expectation.min_f1
        scores.update({"tolerance": expectation.min_f1, "passed": ok})
        metrics["cross_feeding"] = scores
        passed = passed and ok

    return ValidationReport(metrics=metrics, passed=passed, name=expectation.name)


def edges_from_cross_feeding(df) -> List[Edge]:
    """Extract ``(producer, metabolite, consumer)`` edges from a cross_feeding table."""
    if df is None or len(df) == 0:
        return []
    return [(str(r["producer"]), str(r["metabolite"]), str(r["consumer"])) for _, r in df.iterrows()]


def validate(result, expectation: BenchmarkExpectation) -> ValidationReport:
    """Validate a :class:`~muode.dfba.SimulationResult` against an expectation."""
    biomass_final = result.biomass.iloc[-1].to_dict()
    met_final = result.metabolites.iloc[-1].to_dict()
    edges = edges_from_cross_feeding(result.cross_feeding())
    return compare(biomass_final, met_final, edges, expectation)


def validate_outputs(outdir: str | Path, expectation: BenchmarkExpectation) -> ValidationReport:
    """Validate a results directory written by ``to_csv`` (biomass/metabolites/cross_feeding)."""
    import pandas as pd

    outdir = Path(outdir)
    biomass = pd.read_csv(outdir / "biomass.csv", index_col=0).iloc[-1].to_dict()
    met = pd.read_csv(outdir / "metabolites.csv", index_col=0).iloc[-1].to_dict()
    cf_path = outdir / "cross_feeding.csv"
    edges = edges_from_cross_feeding(pd.read_csv(cf_path)) if cf_path.exists() else []
    return compare(biomass, met, edges, expectation)
