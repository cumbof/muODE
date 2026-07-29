"""Organism-level metabolic models and the FBA solve abstraction.

The dynamic-FBA engine (:mod:`muode.dfba`) never talks to a concrete modelling
library directly.  Instead it drives the small :class:`OrganismModel` protocol
defined here.  Two backends implement it:

* :class:`CobraOrganism` -- wraps a genome-scale ``cobra.Model`` (the production
  backend used for real MAG-derived GEMs).  ``cobra`` is imported lazily so the
  rest of the package, and the whole test-suite, work without it.
* :class:`LinprogOrganism` -- a self-contained stoichiometric model solved with
  ``scipy.optimize.linprog`` (HiGHS).  It depends only on numpy/scipy and is used
  for unit tests, the bundled toy cross-feeding example, and as a zero-install
  fallback solver.

Flux sign convention follows COBRA: an *exchange* reaction ``met_e -->`` carries
positive flux for **secretion** and negative flux for **uptake**.  A substrate's
maximum uptake rate is therefore encoded as the (negative) lower bound of its
exchange reaction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping, Optional, Protocol, runtime_checkable

import numpy as np


@dataclass
class OrganismSolution:
    """Result of a single FBA optimisation for one organism."""

    growth_rate: float
    #: exchange flux per *extracellular* metabolite id (secretion > 0, uptake < 0)
    exchange_fluxes: Dict[str, float]
    status: str = "optimal"

    @property
    def feasible(self) -> bool:
        return self.status == "optimal"


@runtime_checkable
class OrganismModel(Protocol):
    """Minimal interface the dynamic-FBA engine relies on.

    A backend must be able to (a) report which extracellular metabolites it can
    exchange, (b) accept a per-step maximum uptake bound for each of them,
    (c) reset to its baseline bounds, and (d) solve the FBA LP.
    """

    id: str

    def exchange_metabolites(self) -> Iterable[str]:
        """Extracellular metabolite ids this organism can take up / secrete."""

    def set_uptake_bound(self, metabolite_id: str, max_uptake: float) -> None:
        """Constrain the maximum uptake rate (>= 0) of ``metabolite_id``."""

    def reset_bounds(self) -> None:
        """Restore baseline reaction bounds (undo per-step uptake settings)."""

    def optimize(self) -> OrganismSolution:
        """Maximise the biomass objective under the current bounds."""


# ---------------------------------------------------------------------------
# scipy/linprog backend (dependency-light, used for tests & the toy example)
# ---------------------------------------------------------------------------


@dataclass
class _Reaction:
    id: str
    stoichiometry: Dict[str, float]
    lower_bound: float = 0.0
    upper_bound: float = 1000.0


class LinprogOrganism:
    """A small mass-balanced metabolic model solved with ``scipy.linprog``.

    Parameters
    ----------
    id:
        Organism identifier.
    reactions:
        Iterable of reaction specs ``(id, {metabolite: coeff}, lb, ub)``.
    objective:
        Reaction id (or mapping of id->weight) used as the biomass objective.
    exchanges:
        Mapping ``extracellular_metabolite_id -> exchange_reaction_id``.  Each
        referenced reaction must have the metabolite as its single participant
        with coefficient ``-1`` (COBRA exchange convention).
    """

    def __init__(
        self,
        id: str,
        reactions: Iterable[_Reaction | tuple],
        objective: str | Mapping[str, float],
        exchanges: Mapping[str, str],
    ) -> None:
        self.id = id
        self._reactions: list[_Reaction] = [
            r if isinstance(r, _Reaction) else _Reaction(r[0], r[1], r[2], r[3])
            for r in reactions
        ]
        self._exchanges = dict(exchanges)

        metabolites = sorted({m for r in self._reactions for m in r.stoichiometry})
        self._met_index = {m: i for i, m in enumerate(metabolites)}
        self._rxn_index = {r.id: j for j, r in enumerate(self._reactions)}

        n_met, n_rxn = len(metabolites), len(self._reactions)
        self._S = np.zeros((n_met, n_rxn))
        for j, r in enumerate(self._reactions):
            for m, coeff in r.stoichiometry.items():
                self._S[self._met_index[m], j] = coeff

        self._lb = np.array([r.lower_bound for r in self._reactions], dtype=float)
        self._ub = np.array([r.upper_bound for r in self._reactions], dtype=float)
        # Baseline copies that perturbations / resets restore to.
        self._base_lb = self._lb.copy()
        self._base_ub = self._ub.copy()

        self._c = np.zeros(n_rxn)
        weights = {objective: 1.0} if isinstance(objective, str) else dict(objective)
        for rxn_id, w in weights.items():
            self._c[self._rxn_index[rxn_id]] = w

    # -- introspection ------------------------------------------------------
    def exchange_metabolites(self) -> Iterable[str]:
        return tuple(self._exchanges)

    def reaction_ids(self) -> Iterable[str]:
        return tuple(self._rxn_index)

    # -- bound management ---------------------------------------------------
    def set_uptake_bound(self, metabolite_id: str, max_uptake: float) -> None:
        rxn_id = self._exchanges.get(metabolite_id)
        if rxn_id is None:
            return
        self._lb[self._rxn_index[rxn_id]] = -abs(float(max_uptake))

    def reset_bounds(self) -> None:
        self._lb = self._base_lb.copy()
        self._ub = self._base_ub.copy()

    def knock_out(self, reaction_ids: Iterable[str]) -> None:
        """Permanently clamp reactions to zero flux (applied to the baseline)."""
        for rid in reaction_ids:
            j = self._rxn_index.get(rid)
            if j is not None:
                self._base_lb[j] = 0.0
                self._base_ub[j] = 0.0
        self.reset_bounds()

    def scale_bounds(self, reaction_ids: Iterable[str], factor: float) -> None:
        """Scale the flux capacity of reactions by ``factor`` (e.g. antibiotic)."""
        for rid in reaction_ids:
            j = self._rxn_index.get(rid)
            if j is not None:
                self._base_lb[j] *= factor
                self._base_ub[j] *= factor
        self.reset_bounds()

    def copy(self) -> "LinprogOrganism":
        """Return an independent deep copy (for per-thread parallel solves)."""
        import copy as _copy

        return _copy.deepcopy(self)

    # -- solve --------------------------------------------------------------
    def optimize(self) -> OrganismSolution:
        from scipy.optimize import linprog

        bounds = list(zip(self._lb, self._ub))
        res = linprog(
            c=-self._c,
            A_eq=self._S,
            b_eq=np.zeros(self._S.shape[0]),
            bounds=bounds,
            method="highs",
        )
        if not res.success or res.x is None:
            return OrganismSolution(0.0, {m: 0.0 for m in self._exchanges}, "infeasible")
        v = res.x
        growth = float(self._c @ v)
        ex = {
            met: float(v[self._rxn_index[rxn]]) for met, rxn in self._exchanges.items()
        }
        return OrganismSolution(growth, ex, "optimal")


# ---------------------------------------------------------------------------
# cobra backend (production GEMs)
# ---------------------------------------------------------------------------


class CobraOrganism:
    """Wraps a genome-scale ``cobra.Model`` behind the :class:`OrganismModel` API.

    ``cobra`` is imported lazily, so importing :mod:`muode` never requires it.
    """

    def __init__(self, model, id: Optional[str] = None) -> None:
        self.model = model
        self.id = id or model.id
        # GLPK's simplex can *cycle indefinitely* on the degenerate LPs that
        # community GEMs routinely produce -- observed as a single solve pinning
        # a core for >21 h with no iteration cap.  Presolve collapses most of that
        # degeneracy, and a per-solve wall-clock ceiling turns any residual
        # pathological LP into a non-optimal status (handled as zero growth in
        # ``optimize``) instead of an unbounded hang.  Both are solver-config, so
        # they ride along through ``model.copy()`` used in the dynamic loop.
        try:
            cfg = model.solver.configuration
            cfg.presolve = True
            cfg.timeout = 60  # seconds per LP; normal solves are sub-second
        except Exception:
            pass
        # Map extracellular metabolite id -> exchange reaction id.
        self._exchanges: Dict[str, str] = {}
        for rxn in model.exchanges:
            mets = list(rxn.metabolites)
            if len(mets) == 1:
                self._exchanges[mets[0].id] = rxn.id
        self._base_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}

    @classmethod
    def from_file(cls, path: str, id: Optional[str] = None) -> "CobraOrganism":
        import cobra

        model = cobra.io.read_sbml_model(str(path))
        return cls(model, id=id)

    def exchange_metabolites(self) -> Iterable[str]:
        return tuple(self._exchanges)

    def reaction_ids(self) -> Iterable[str]:
        return tuple(self._base_bounds)

    def reactions_in_subsystem(self, subsystem: str) -> list[str]:
        """Reaction ids whose subsystem/pathway matches ``subsystem``.

        Matching is case-insensitive and substring-based so that loose pathway
        names from the CLI (e.g. "Folate Biosynthesis") resolve against the
        GEM's annotated subsystems.
        """
        needle = subsystem.lower()
        return [
            r.id
            for r in self.model.reactions
            if r.subsystem and needle in r.subsystem.lower()
        ]

    def set_uptake_bound(self, metabolite_id: str, max_uptake: float) -> None:
        rxn_id = self._exchanges.get(metabolite_id)
        if rxn_id is not None:
            self.model.reactions.get_by_id(rxn_id).lower_bound = -abs(float(max_uptake))

    def reset_bounds(self) -> None:
        for rid, (lb, ub) in self._base_bounds.items():
            rxn = self.model.reactions.get_by_id(rid)
            rxn.lower_bound, rxn.upper_bound = lb, ub

    def knock_out(self, reaction_ids: Iterable[str]) -> None:
        for rid in reaction_ids:
            if rid in self._base_bounds:
                self._base_bounds[rid] = (0.0, 0.0)
        self.reset_bounds()

    def scale_bounds(self, reaction_ids: Iterable[str], factor: float) -> None:
        for rid in reaction_ids:
            if rid in self._base_bounds:
                lb, ub = self._base_bounds[rid]
                self._base_bounds[rid] = (lb * factor, ub * factor)
        self.reset_bounds()

    def copy(self) -> "CobraOrganism":
        """Return an independent copy with its own cobra model.

        Gives each worker thread private models for parallel per-cell solves.
        ``cobra.Model.copy`` duplicates the solver (so an enzyme-pool constraint
        carries over); the exchange map and baseline bounds are copied verbatim
        so the copy resets to the same state as the original.
        """
        new = CobraOrganism.__new__(CobraOrganism)
        new.model = self.model.copy()
        new.id = self.id
        new._exchanges = dict(self._exchanges)
        new._base_bounds = dict(self._base_bounds)
        return new

    def apply_enzyme_constraints(self, kinetics, organism_id: Optional[str] = None,
                                 enzyme_concentration: Optional[float] = None,
                                 pool_budget: Optional[float] = None) -> dict:
        """Cap intracellular reaction velocities from kcat (GECKO-lite).

        The per-reaction caps become part of this organism's *baseline*, so they
        persist across the per-step ``reset_bounds`` of the dynamic loop.  When
        ``pool_budget`` (g enzyme / gDW) is given, a single shared protein-pool
        constraint is added on top via :func:`muode.enzyme.apply_protein_pool_constraint`;
        it lives on the solver, so it likewise survives every per-step reset.
        """
        from muode.enzyme import apply_enzyme_constraints, apply_protein_pool_constraint

        report = apply_enzyme_constraints(
            self.model, kinetics, organism_id or self.id, enzyme_concentration
        )
        # fold the new caps into the baseline the dFBA loop resets to each step
        self._base_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in self.model.reactions}
        if pool_budget is not None:
            report.update(apply_protein_pool_constraint(
                self.model, kinetics, organism_id or self.id, pool_budget
            ))
        return report

    def optimize(self) -> OrganismSolution:
        sol = self.model.optimize()
        if sol.status != "optimal":
            return OrganismSolution(0.0, {m: 0.0 for m in self._exchanges}, sol.status)
        growth = float(sol.objective_value or 0.0)
        ex = {
            met: float(sol.fluxes.get(rxn, 0.0)) for met, rxn in self._exchanges.items()
        }
        return OrganismSolution(growth, ex, "optimal")
