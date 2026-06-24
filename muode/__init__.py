"""muODE (µODE) — automated ODE-based simulation engine for microbial community dynamics.

The package is organised around a small, dependency-light core (the dynamic-FBA
integrator, the community container and the perturbation engine) plus a set of
optional wrappers around external bioinformatics / deep-learning tools that make
up the upstream reconstruction phases.

Public, stable entry points are re-exported here for convenience.
"""

from __future__ import annotations

__version__ = "0.1.0"

from muode.antagonism import Bacteriocin
from muode.antibiotic import Antibiotic
from muode.bile import BileAcidInhibition, BileAcidTransform
from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.ecology import EcologyLayer, EcologyModel
from muode.inject import Injection, merge_for_injection
from muode.kinetics import KineticParameters
from muode.lifecycle import SporeForming
from muode.oxygen import OxygenSensitivity
from muode.ph import WeakAcidInhibition
from muode.phage import PhageInfection
from muode.scenarios import cdi_scenario, phage_predation_scenario
from muode.organism import LinprogOrganism, OrganismModel, OrganismSolution
from muode.perturb import Perturbation
from muode.traits import Domain, MicrobeTraits, OxygenTolerance, reconstruction_route
from muode.predict import HeuristicPredictor, refine_kinetics
from muode.spatial import SpatialDynamicFBA, SpatialResult
from muode.validate import BenchmarkExpectation, ValidationReport, validate

__all__ = [
    "__version__",
    "Community",
    "DynamicFBA",
    "SimulationResult",
    "Diet",
    "Injection",
    "merge_for_injection",
    "KineticParameters",
    "EcologyLayer",
    "EcologyModel",
    "WeakAcidInhibition",
    "BileAcidTransform",
    "BileAcidInhibition",
    "SporeForming",
    "Antibiotic",
    "Bacteriocin",
    "OxygenSensitivity",
    "PhageInfection",
    "Domain",
    "MicrobeTraits",
    "OxygenTolerance",
    "reconstruction_route",
    "cdi_scenario",
    "phage_predation_scenario",
    "LinprogOrganism",
    "OrganismModel",
    "OrganismSolution",
    "Perturbation",
    "HeuristicPredictor",
    "refine_kinetics",
    "SpatialDynamicFBA",
    "SpatialResult",
    "BenchmarkExpectation",
    "ValidationReport",
    "validate",
]
