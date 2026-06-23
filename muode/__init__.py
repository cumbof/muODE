"""muODE (µODE) — automated ODE-based simulation engine for microbial community dynamics.

The package is organised around a small, dependency-light core (the dynamic-FBA
integrator, the community container and the perturbation engine) plus a set of
optional wrappers around external bioinformatics / deep-learning tools that make
up the upstream reconstruction phases.

Public, stable entry points are re-exported here for convenience.
"""

from __future__ import annotations

__version__ = "0.1.0"

from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.kinetics import KineticParameters
from muode.organism import LinprogOrganism, OrganismModel, OrganismSolution
from muode.perturb import Perturbation
from muode.predict import HeuristicPredictor, refine_kinetics
from muode.validate import BenchmarkExpectation, ValidationReport, validate

__all__ = [
    "__version__",
    "Community",
    "DynamicFBA",
    "SimulationResult",
    "Diet",
    "KineticParameters",
    "LinprogOrganism",
    "OrganismModel",
    "OrganismSolution",
    "Perturbation",
    "HeuristicPredictor",
    "refine_kinetics",
    "BenchmarkExpectation",
    "ValidationReport",
    "validate",
]
