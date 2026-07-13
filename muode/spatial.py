"""Spatiotemporal dynamic FBA -- 2D colony / biofilm simulation (M4, reach).

The well-mixed engine (:mod:`muode.dfba`) treats the whole community as one
beaker.  Real colonies and biofilms are *spatial*: nutrients and secreted
metabolites form gradients, and who-grows-where depends on local conditions.
This module generalises the Static Optimization Approach to a 2D lattice, the
discretisation used by spatial community-FBA tools such as COMETS:

* the world is an ``ny x nx`` grid of cells, each a small well-mixed FBA reactor;
* every cell holds a local biomass per species and a local concentration per
  extracellular metabolite;
* at each step every populated cell solves its species' FBA under local Michaelis
  -Menten uptake bounds, biomass grows locally, and metabolites are produced /
  consumed locally;
* metabolites then **diffuse** between neighbouring cells (a finite-difference
  reaction-diffusion step, ``dM/dt = D nabla^2 M + ...``), with optional biomass
  spreading, under no-flux (Neumann) boundaries so mass is conserved.

The per-cell LP solves make this a proof-of-concept for genome-scale communities
(it is cheap with the linprog backend used by the toy example); scaling to large
GEMs on large grids is an optimisation problem in its own right.
"""

from __future__ import annotations

import os
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy as _deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from muode.community import Community
from muode.diet import Diet
from muode.kinetics import KineticParameters

#: default metabolite diffusivity (mm^2 / h), order-of-magnitude for small
#: molecules in aqueous gel; override per metabolite via ``diffusivity=``.
DEFAULT_DIFFUSIVITY = 1.0


def laplacian(field: np.ndarray, dx: float) -> np.ndarray:
    """5-point finite-difference Laplacian with no-flux (Neumann) boundaries."""
    p = np.pad(field, 1, mode="edge")  # edge padding -> zero gradient at walls
    return (
        p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:] - 4.0 * p[1:-1, 1:-1]
    ) / (dx * dx)


def diffuse_step(field: np.ndarray, diffusivity: float, dx: float, dt: float) -> np.ndarray:
    """One explicit diffusion step (returns a new field)."""
    return field + diffusivity * laplacian(field, dx) * dt


# ---------------------------------------------------------------------------
# inoculum helpers (initial spatial biomass per species)
# ---------------------------------------------------------------------------


def uniform_inoculum(shape: Tuple[int, int], ids: Iterable[str], amount: float) -> Dict[str, np.ndarray]:
    """Spread ``amount`` gDW of each species evenly over the grid."""
    ny, nx = shape
    per_cell = amount / (ny * nx)
    return {i: np.full(shape, per_cell) for i in ids}


def point_inoculum(shape: Tuple[int, int], placements: Dict[str, Tuple[int, int, float]]) -> Dict[str, np.ndarray]:
    """Place each species as a point colony at ``(row, col, amount)``."""
    fields = {sp: np.zeros(shape) for sp in placements}
    for sp, (i, j, amount) in placements.items():
        fields[sp][i, j] = amount
    return fields


def halves_inoculum(shape: Tuple[int, int], left_id: str, right_id: str,
                    amount: float, axis: int = 1) -> Dict[str, np.ndarray]:
    """Inoculate one species in one half of the grid and another in the other."""
    ny, nx = shape
    left, right = np.zeros(shape), np.zeros(shape)
    mid = (nx if axis == 1 else ny) // 2
    if axis == 1:
        left[:, :mid] = amount / (ny * mid)
        right[:, mid:] = amount / (ny * (nx - mid))
    else:
        left[:mid, :] = amount / (mid * nx)
        right[mid:, :] = amount / ((ny - mid) * nx)
    return {left_id: left, right_id: right}


# ---------------------------------------------------------------------------
# result container
# ---------------------------------------------------------------------------


@dataclass
class SpatialResult:
    """Recorded frames of a spatial run: fields of shape ``(n_frames, ny, nx)``."""

    times: np.ndarray
    biomass: Dict[str, np.ndarray]
    metabolites: Dict[str, np.ndarray]
    meta: dict = field(default_factory=dict)

    def total_biomass(self):
        """Total biomass per species over time (summed across the grid)."""
        import pandas as pd

        return pd.DataFrame(
            {sp: arr.reshape(arr.shape[0], -1).sum(axis=1) for sp, arr in self.biomass.items()},
            index=self.times,
        )

    def final_biomass(self, species: str) -> np.ndarray:
        return self.biomass[species][-1]

    def final_metabolite(self, metabolite: str) -> np.ndarray:
        return self.metabolites[metabolite][-1]

    def to_npz(self, outdir: str | Path) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            outdir / "spatial.npz",
            times=self.times,
            **{f"biomass__{k}": v for k, v in self.biomass.items()},
            **{f"metabolite__{k}": v for k, v in self.metabolites.items()},
        )
        self.total_biomass().to_csv(outdir / "spatial_total_biomass.csv", index_label="time_h")


# ---------------------------------------------------------------------------
# engine
# ---------------------------------------------------------------------------


@dataclass
class SpatialDynamicFBA:
    """2D reaction-diffusion dynamic-FBA integrator.

    ``n_jobs`` sets the number of worker threads for the per-cell FBA solves
    (``1`` = sequential, the default; ``-1`` = all cores).  Populated cells solve
    independently, so any value yields identical fields; the speed-up grows with
    the grid size and the model size.
    """

    nx: int = 20
    ny: int = 1
    dx: float = 1.0                 # cell length (mm)
    t_end: float = 12.0
    dt: float = 0.05
    diffusivity: Dict[str, float] = field(default_factory=dict)
    default_diffusivity: float = DEFAULT_DIFFUSIVITY
    biomass_diffusivity: float = 0.0
    death_rate: float = 0.0
    min_biomass: float = 1e-12
    record_every: int = 1
    n_jobs: int = 1

    def _D(self, metabolite: str) -> float:
        return float(self.diffusivity.get(metabolite, self.default_diffusivity))

    def run(
        self,
        community: Community,
        diet: Diet,
        kinetics: Optional[KineticParameters] = None,
        inoculum: Optional[Dict[str, np.ndarray]] = None,
    ) -> SpatialResult:
        kinetics = kinetics or KineticParameters()
        shape = (self.ny, self.nx)
        organisms = community.organisms
        env_mets = sorted(set(community.environment_metabolites()) | set(diet.metabolites()))

        if inoculum is None:
            inoculum = uniform_inoculum(shape, [o.id for o in organisms], community.total_biomass)
        X = {o.id: np.array(inoculum.get(o.id, np.zeros(shape)), dtype=float) for o in organisms}
        M = {m: np.full(shape, diet.initial_concentration(m), dtype=float) for m in env_mets}

        # diffusion stability (explicit scheme): D*dt/dx^2 must stay < 0.25 in 2D
        max_courant = max((self._D(m) for m in env_mets), default=0.0) * self.dt / (self.dx ** 2)
        if max_courant > 0.25:
            warnings.warn(
                f"diffusion may be unstable (D*dt/dx^2 = {max_courant:.2f} > 0.25); "
                "reduce dt or increase dx.",
                RuntimeWarning, stacklevel=2,
            )

        n_steps = int(round(self.t_end / self.dt))
        bio_frames: Dict[str, list] = {o.id: [] for o in organisms}
        met_frames: Dict[str, list] = {m: [] for m in env_mets}
        rec_times: list = []

        org_ids = [o.id for o in organisms]
        shared_orgs = {o.id: o for o in organisms}

        # Populated cells solve independently -- each reads only its own cell's
        # concentrations and biomass -- so they can run on a thread pool.  Unlike
        # the well-mixed engine, the *same* organism object is shared across cells,
        # so each worker thread is given private copies of the models (thread-local)
        # to avoid races; results are applied on the main thread, making a parallel
        # run identical to a sequential one.  n_jobs=1 keeps the in-place path.
        if self.n_jobs < 0:
            n_workers = os.cpu_count() or 1
        else:
            n_workers = self.n_jobs
        n_workers = max(1, n_workers)
        pool = ThreadPoolExecutor(max_workers=n_workers) if n_workers > 1 else None
        _tls = threading.local()

        def _worker_orgs():
            local = getattr(_tls, "orgs", None)
            if local is None:
                local = {o.id: (o.copy() if hasattr(o, "copy") else _deepcopy(o))
                         for o in organisms}
                _tls.orgs = local
            return local

        for step in range(n_steps + 1):
            if step % self.record_every == 0 or step == n_steps:
                rec_times.append(step * self.dt)
                for o in organisms:
                    bio_frames[o.id].append(X[o.id].copy())
                for m in env_mets:
                    met_frames[m].append(M[m].copy())
            if step == n_steps:
                break

            total = np.zeros(shape)
            for o in organisms:
                total += X[o.id]
            dM = {m: np.zeros(shape) for m in env_mets}

            # Solve one populated cell: every resident species under the cell's
            # local uptake bounds.  Reads only this cell's slice of the read-only
            # X and M snapshots and its worker's private models, so distinct cells
            # solve concurrently; the caller applies the returned updates, so all
            # writes are deferred to after the join (no cell reads another's cell).
            def _solve_cell(cell):
                i, j = cell
                orgs = _worker_orgs() if pool is not None else shared_orgs
                new_x: Dict[str, float] = {}
                dm_cell = {m: 0.0 for m in env_mets}
                for oid in org_ids:
                    xij = X[oid][i, j]
                    if xij <= self.min_biomass:
                        continue
                    o = orgs[oid]
                    o.reset_bounds()
                    for m in env_mets:
                        conc = M[m][i, j]
                        mm = kinetics.michaelis_menten(oid, m, conc)
                        limit = diet.uptake_limit(m)      # dietary availability
                        if limit is not None:
                            mm = min(mm, limit)
                        cap = conc / (xij * self.dt) if self.dt > 0 else np.inf
                        o.set_uptake_bound(m, min(mm, cap))
                    sol = o.optimize()
                    mu = sol.growth_rate if sol.feasible else 0.0
                    new_x[oid] = max(0.0, xij + (mu - self.death_rate) * xij * self.dt)
                    for m in env_mets:
                        dm_cell[m] += sol.exchange_fluxes.get(m, 0.0) * xij
                return i, j, new_x, dm_cell

            cells = [tuple(c) for c in np.argwhere(total > self.min_biomass)]
            if pool is None or len(cells) <= 1:
                results = [_solve_cell(c) for c in cells]
            else:
                results = list(pool.map(_solve_cell, cells))

            for i, j, new_x, dm_cell in results:
                for oid, xv in new_x.items():
                    X[oid][i, j] = xv
                for m in env_mets:
                    dM[m][i, j] += dm_cell[m]

            # reaction + influx, then diffusion (Neumann BC) for each metabolite
            for m in env_mets:
                M[m] += dM[m] * self.dt + diet.influx_rate(m) * self.dt
                D = self._D(m)
                if D > 0:
                    M[m] += D * laplacian(M[m], self.dx) * self.dt
                np.clip(M[m], 0.0, None, out=M[m])

            if self.biomass_diffusivity > 0:
                for o in organisms:
                    X[o.id] += self.biomass_diffusivity * laplacian(X[o.id], self.dx) * self.dt
                    np.clip(X[o.id], 0.0, None, out=X[o.id])

        if pool is not None:
            pool.shutdown()

        return SpatialResult(
            times=np.array(rec_times),
            biomass={k: np.stack(v) for k, v in bio_frames.items()},
            metabolites={k: np.stack(v) for k, v in met_frames.items()},
            meta={
                "grid": [self.ny, self.nx], "dx": self.dx, "t_end": self.t_end, "dt": self.dt,
                "diffusivity_default": self.default_diffusivity,
                "biomass_diffusivity": self.biomass_diffusivity,
                "n_species": len(organisms), "diet": diet.name,
            },
        )
