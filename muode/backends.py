"""LP solver backends for the dynamic-FBA engines.

Every time step of :class:`muode.dfba.DynamicFBA` needs one FBA LP per active
species (and, in :meth:`~muode.dfba.DynamicFBA.run_ensemble`, per ensemble
member).  The engine no longer calls ``organism.optimize()`` itself: it collects
the step's LPs as :class:`SolveRequest` objects and hands the *whole batch* to a
:class:`SolverBackend`.  That is what lets a batched solver see thousands of LPs
at once instead of one at a time.

Backends
--------
``LegacyBackend`` (default)
    Exactly the historical behaviour: per request ``reset_bounds()``, one
    ``set_uptake_bound`` per environment metabolite, ``optimize()``, optionally
    on a thread pool.  Outputs are identical to muODE before the refactor.
``ManyLPBackend``
    The organisms' LPs are extracted once and solved by the ``manylp`` package:
    ``engine="batched"`` uses its certify-and-repair batched solver (CPU or
    GPU); ``engine="highs"`` solves every LP individually with a persistent,
    warm-started HiGHS model -- the strong CPU reference, with the same
    flux-uniqueness rule.

Flux uniqueness
---------------
FBA fixes the growth rate, not the flux vector; the exchange fluxes that drive
the ODE are then whatever vertex a solver happens to return.  ``ManyLPBackend``
applies an explicit rule (``flux_rule``):

``"vertex"``        plain FBA (an optimal vertex; may differ between solvers)
``"pfba"``          max growth, then min sum |v|         (parsimonious FBA)
``"pfba-unique"``   ... then a fixed generic tie-break  (default: unique, deterministic)

Only LPs the backend can express *exactly* are extracted; anything else (e.g. a
cobra model carrying extra solver constraints such as an enzyme pool) falls back
to the legacy path for that organism, so results are never silently changed.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence

import numpy as np

from muode.organism import CobraOrganism, LinprogOrganism, MuMaxCapped, OrganismSolution


@dataclass
class SolveRequest:
    """One FBA LP of one time step."""

    organism: object
    #: maximum uptake rate (>= 0) per environment metabolite, in the engine's order;
    #: may be ``None`` when ``mets``/``vec`` are given (built on demand)
    uptake: Optional[Mapping[str, float]]
    #: ensemble member (warm starts are tracked per organism and member)
    member: int = 0
    #: vectorised form: metabolites (the organism's exchange order) and their uptake bounds
    mets: Optional[tuple] = None
    vec: Optional[np.ndarray] = None

    def uptake_items(self):
        if self.uptake is not None:
            return self.uptake.items()
        return zip(self.mets, self.vec.tolist())


class SolverBackend:
    """Solve all LPs of one time step.  Subclasses override :meth:`solve_step`."""

    name = "abstract"

    def solve_step(self, requests: Sequence[SolveRequest]) -> List[OrganismSolution]:
        raise NotImplementedError

    def close(self) -> None:
        pass

    def stats(self) -> dict:
        return {}


# ---------------------------------------------------------------------------
# legacy: one optimize() per request (bit-identical to the historical engine)
# ---------------------------------------------------------------------------


class LegacyBackend(SolverBackend):
    """Historical per-organism path; ``n_jobs`` threads across organisms."""

    name = "legacy"

    def __init__(self, n_jobs: int = 1) -> None:
        self.n_jobs = n_jobs
        self._pool: Optional[ThreadPoolExecutor] = None

    @staticmethod
    def _solve_one(req: SolveRequest) -> OrganismSolution:
        o = req.organism
        o.reset_bounds()
        for m, v in req.uptake_items():
            o.set_uptake_bound(m, v)
        return o.optimize()

    def solve_step(self, requests):
        workers = (os.cpu_count() or 1) if self.n_jobs < 0 else self.n_jobs
        # requests that share an organism object must run sequentially (they
        # mutate its bounds), so parallelise over organisms, not requests
        by_org: Dict[int, List[int]] = {}
        for i, r in enumerate(requests):
            by_org.setdefault(id(r.organism), []).append(i)
        out: List[Optional[OrganismSolution]] = [None] * len(requests)

        def run(idx):
            for i in idx:
                out[i] = self._solve_one(requests[i])

        groups = list(by_org.values())
        if workers > 1 and len(groups) > 1:
            if self._pool is None:
                self._pool = ThreadPoolExecutor(max_workers=workers)
            list(self._pool.map(run, groups))
        else:
            for idx in groups:
                run(idx)
        return out  # type: ignore[return-value]

    def close(self):
        if self._pool is not None:
            self._pool.shutdown()
            self._pool = None


# ---------------------------------------------------------------------------
# LP extraction
# ---------------------------------------------------------------------------


class NotBatchable(Exception):
    """The organism's LP cannot be expressed exactly as ``S v = 0, lb <= v <= ub``."""


@dataclass
class OrganismLP:
    S: object                         # scipy.sparse (m x n)
    lb: np.ndarray
    ub: np.ndarray
    c: np.ndarray
    ex_mets: List[str]                # extracellular metabolite of each exchange
    ex_idx: np.ndarray                # reaction index of each exchange
    post: List[Callable] = field(default_factory=list)


def _cap_growth(mumax: float) -> Callable[[OrganismSolution], OrganismSolution]:
    def post(sol: OrganismSolution) -> OrganismSolution:
        g = sol.growth_rate
        if sol.feasible and g > mumax > 0.0:
            s = mumax / g
            return OrganismSolution(mumax, {m: v * s for m, v in sol.exchange_fluxes.items()},
                                    sol.status)
        return sol
    return post


def extract_lp(o) -> OrganismLP:
    """The organism's baseline LP (after any perturbation applied to its baseline)."""
    import scipy.sparse as sp

    if isinstance(o, MuMaxCapped):
        inner = extract_lp(o._inner)
        inner.post.append(_cap_growth(o._mumax))   # same post-solve rule as MuMaxCapped
        return inner
    if isinstance(o, LinprogOrganism):
        mets = list(o._exchanges)
        return OrganismLP(S=sp.csc_matrix(o._S), lb=o._base_lb.copy(), ub=o._base_ub.copy(),
                          c=o._c.copy(), ex_mets=mets,
                          ex_idx=np.array([o._rxn_index[o._exchanges[m]] for m in mets], dtype=np.int64))
    if isinstance(o, CobraOrganism):
        from cobra.util.array import create_stoichiometric_matrix
        from cobra.util.solver import linear_reaction_coefficients

        model = o.model
        if len(model.constraints) != len(model.metabolites):
            raise NotBatchable(f"{o.id}: solver carries constraints beyond mass balance")
        rxns = model.reactions
        idx = {r.id: i for i, r in enumerate(rxns)}
        lb = np.array([o._base_bounds[r.id][0] for r in rxns], dtype=float)
        ub = np.array([o._base_bounds[r.id][1] for r in rxns], dtype=float)
        c = np.zeros(len(rxns))
        for r, w in linear_reaction_coefficients(model).items():
            c[idx[r.id]] = float(w)
        if model.objective_direction == "min":
            c = -c
        mets = list(o._exchanges)
        return OrganismLP(S=sp.csc_matrix(create_stoichiometric_matrix(model, array_type="dok")),
                          lb=lb, ub=ub, c=c, ex_mets=mets,
                          ex_idx=np.array([idx[o._exchanges[m]] for m in mets], dtype=np.int64))
    raise NotBatchable(f"{getattr(o, 'id', o)!r}: unknown organism type {type(o).__name__}")


# ---------------------------------------------------------------------------
# manylp: batched certify-and-repair (CPU / GPU) or per-LP warm HiGHS
# ---------------------------------------------------------------------------


class ManyLPBackend(SolverBackend):
    """Solve each step's LPs with ``manylp`` (see module docstring)."""

    def __init__(self, device: str = "cpu", engine: str = "batched", flux_rule: str = "pfba-unique",
                 n_workers: int = 16, seed: int = 0, atlas_dir: Optional[str] = None,
                 policies: Optional[Sequence[str]] = None) -> None:
        """
        atlas_dir:
            Directory of persistent basis atlases (one file per distinct LP).  Loaded
            when an organism is first seen and written back on :meth:`close`, so later
            simulations with the same GEMs skip nearly all simplex work.  Atlases hold
            certificates, not answers: every LP is still certified against its own
            bounds, so they can only save time, never change a result.
        policies:
            Alternative-optima selection policies, one per ensemble member
            (``member % len(policies)``): ``"canonical"``, ``"random:<seed>"``,
            ``"max:<metabolite>"`` or ``"min:<metabolite>"``.  Each is applied after
            the growth objective and before the flux rule's own stages, so every member
            still gets a unique, certified flux.  Integrating the same diet with many
            policies (``DynamicFBA.run_ensemble``) yields the *alternative-optima
            envelope*: the part of a dFBA prediction that FBA itself leaves undetermined.
        """
        from manylp import BatchLPSolver

        if engine not in ("batched", "highs"):
            raise ValueError("engine must be 'batched' or 'highs'")
        mode = {"vertex": "fba", "pfba": "pfba", "pfba-unique": "pfba-unique"}.get(flux_rule)
        if mode is None:
            raise ValueError("flux_rule must be 'vertex', 'pfba' or 'pfba-unique'")
        self.mode = mode
        self.engine = engine
        self.seed = seed
        self.name = f"manylp-{engine}-{device}-{flux_rule}"
        self.solver = BatchLPSolver(device=device, n_workers=n_workers) if engine == "batched" else None
        self.n_workers = n_workers
        self._groups: Dict[int, tuple] = {}        # id(organism) -> (spec, problem, handle | None)
        self._warm: Dict[int, Dict[int, object]] = {}
        self._legacy = LegacyBackend()
        self._pool: Optional[ThreadPoolExecutor] = None
        self.atlas_dir = os.path.expanduser(atlas_dir) if atlas_dir else None
        self.policies = list(policies) if policies else ["canonical"]
        self._specs: Dict[int, Optional[OrganismLP]] = {}
        self.n_lps = 0
        self.n_fallback = 0
        self.atlas_loaded = 0

    def _spec(self, o) -> Optional[OrganismLP]:
        key = id(o)
        if key not in self._specs:
            try:
                self._specs[key] = extract_lp(o)
            except NotBatchable:
                self._specs[key] = None
        return self._specs[key]

    def _policy_objectives(self, spec: OrganismLP, q: int) -> list:
        kind, _, arg = self.policies[q].partition(":")
        n = spec.c.size
        if kind == "canonical":
            return []
        if kind == "random":
            rng = np.random.default_rng(int(arg) + 1_000_003)
            return [(rng.normal(size=n), "max")]
        if kind in ("max", "min"):
            if arg not in spec.ex_mets:
                return []                       # this organism cannot exchange it
            e = np.zeros(n)
            e[spec.ex_idx[spec.ex_mets.index(arg)]] = 1.0   # > 0 is secretion
            return [(e, kind)]
        raise ValueError(f"unknown policy {self.policies[q]!r}")

    def _group(self, o, q: int = 0):
        key = (id(o), q)
        g = self._groups.get(key)
        if g is None:
            from manylp.fba import FBAModel, compile_fba

            spec = self._spec(o)
            if spec is None:
                g = self._groups[key] = (None, None, None)
                return g
            model = FBAModel(S=spec.S, lb=spec.lb, ub=spec.ub, c=spec.c, exchanges=spec.ex_idx,
                             name=getattr(o, "id", ""))
            prob = compile_fba(model, mode=self.mode, seed=self.seed,
                               extra_objectives=self._policy_objectives(spec, q))
            handle = None
            if self.solver is not None:
                handle = self.solver.register_group(prob.lp, out_z=prob.out_z(model.exchanges))
                if self.atlas_dir:
                    self.atlas_loaded += self.solver.load_atlas(handle, self._atlas_path(handle))
            g = self._groups[key] = (spec, prob, handle)
            self._warm[key] = {}
        return g

    def _atlas_path(self, handle) -> str:
        from manylp.solver import lp_fingerprint

        os.makedirs(self.atlas_dir, exist_ok=True)
        return os.path.join(self.atlas_dir, lp_fingerprint(handle.lp)[:24] + ".npz")

    def save_atlas(self) -> int:
        """Write every group's certified bases to ``atlas_dir``."""
        if not (self.atlas_dir and self.solver is not None):
            return 0
        return sum(self.solver.save_atlas(h, self._atlas_path(h))
                   for _, _, h in self._groups.values() if h is not None)

    def _ex_lb(self, spec: OrganismLP, reqs: Sequence[SolveRequest]) -> np.ndarray:
        """Exchange lower bounds per request, in ``spec.ex_mets`` order."""
        base = spec.lb[spec.ex_idx]
        out = np.empty((len(reqs), spec.ex_idx.size))
        mets = spec.ex_mets
        tmets = tuple(mets)
        for b, r in enumerate(reqs):
            if r.vec is not None and r.mets == tmets:
                out[b] = -np.abs(r.vec)            # already in this organism's exchange order
                continue
            row = base.copy()
            up = r.uptake if r.uptake is not None else dict(r.uptake_items())
            for i, m in enumerate(mets):
                v = up.get(m)
                if v is not None:
                    row[i] = -abs(float(v))
            out[b] = row
        return out

    def solve_step(self, requests):
        out: List[Optional[OrganismSolution]] = [None] * len(requests)
        by_org: Dict[tuple, List[int]] = {}
        P = len(self.policies)
        for i, r in enumerate(requests):
            by_org.setdefault((id(r.organism), r.member % P), []).append(i)
        for key, idx in by_org.items():
            reqs = [requests[i] for i in idx]
            spec, prob, handle = self._group(reqs[0].organism, key[1])
            if spec is None:
                self.n_fallback += len(reqs)
                for i, s in zip(idx, self._legacy.solve_step(reqs)):
                    out[i] = s
                continue
            ex_lb = self._ex_lb(spec, reqs)
            members = [r.member for r in reqs]
            growth, flux, status = self._solve_group(key, prob, handle, ex_lb, members)
            for j, i in enumerate(idx):
                if status[j]:
                    sol = OrganismSolution(float(growth[j]),
                                           {m: float(v) for m, v in zip(spec.ex_mets, flux[j])}, "optimal")
                else:
                    sol = OrganismSolution(0.0, {m: 0.0 for m in spec.ex_mets}, "infeasible")
                for post in spec.post:
                    sol = post(sol)
                out[i] = sol
            self.n_lps += len(reqs)
        return out  # type: ignore[return-value]

    def _solve_group(self, key, prob, handle, ex_lb, members):
        warm = self._warm[key]
        Lp, Up = prob.param_bounds(ex_lb)
        ex = prob.model.exchanges
        if handle is not None:
            ws = np.array([warm.get(mb, -1) for mb in members], dtype=np.int64)
            be = self.solver.backend
            if be.is_gpu and len(members) >= self.solver.gpu_min_batch:
                # device-resident: bounds built and fluxes recovered on the GPU
                with be.device_ctx():
                    Ld, Ud = prob.param_bounds(be.asarray(ex_lb))
                    sol = self.solver.solve_batch(handle, Ld, Ud, warm_start=ws, device_out=True)
                    growth = be.to_host(sol.objective[:, 0])
                    flux = be.to_host(prob.fluxes(sol.z, ex))
            else:
                sol = self.solver.solve_batch(handle, Lp, Up, warm_start=ws)
                growth, flux = sol.objective[:, 0], prob.fluxes(sol.z, ex)
            for mb, bid in zip(members, sol.basis_id):
                warm[mb] = int(bid)
            return growth, flux, sol.status == 1
        # engine == "highs": one warm-started lexicographic solve per LP
        from manylp.repair import HighsLexSolver

        if "solver" not in warm:
            warm["solver"] = HighsLexSolver(prob.lp)
        hs = warm["solver"]
        oz = prob.out_z(ex)
        B = Lp.shape[0]
        growth = np.zeros(B)
        Z = np.zeros((B, oz.size))
        ok = np.zeros(B, dtype=bool)
        for b in range(B):
            lb, ub = prob.lp.full_bounds(Lp[b], Up[b])
            r = hs.solve(lb, ub, warm_status=warm.get(("basis", members[b])))
            if r.status == 1:
                ok[b] = True
                growth[b] = r.stage_obj[0]
                Z[b] = r.z[oz]
                warm[("basis", members[b])] = r.zstatus
        return growth, prob.fluxes(Z, ex), ok

    def stats(self) -> dict:
        st = {"lps": self.n_lps, "fallback_lps": self.n_fallback}
        if self.solver is not None:
            tot: Dict[str, float] = {}
            for spec, prob, h in self._groups.values():
                if h is not None:
                    for k, v in h.totals.items():
                        tot[k] = tot.get(k, 0) + v
            st.update(tot)
        return st

    def close(self):
        self.save_atlas()
        if self.solver is not None:
            self.solver.close()
        self._legacy.close()


def make_backend(spec: Optional[str] = None, n_jobs: int = 1, **kw) -> SolverBackend:
    """``None``/``"legacy"``, ``"manylp"``/``"manylp-gpu"``/``"manylp-cpu"``, or ``"highs"``."""
    if spec is None or spec == "legacy":
        return LegacyBackend(n_jobs=n_jobs)
    if spec in ("manylp", "manylp-cpu"):
        return ManyLPBackend(device="cpu", **kw)
    if spec in ("manylp-gpu", "manylp-cuda"):
        return ManyLPBackend(device="cuda", **kw)
    if spec == "highs":
        return ManyLPBackend(engine="highs", **kw)
    raise ValueError(f"unknown backend {spec!r}")
