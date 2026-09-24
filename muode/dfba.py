"""Dynamic Flux Balance Analysis -- the muODE simulation engine.

This implements the **Static Optimization Approach** (SOA) to dynamic FBA
(Mahadevan, Edwards & Doyle, *Biophys. J.* 2002), extended to communities of
several genome-scale models sharing one extracellular pool.

At each time step of size ``dt``:

1. For every species *i* and every extracellular metabolite *j*, the maximum
   uptake rate is set from the current concentration via Michaelis-Menten
   kinetics, ``v_max,ij = Vmax * M_j / (Km + M_j)`` (:mod:`muode.kinetics`),
   additionally capped so a species cannot consume more of *j* than is present
   within the step (``M_j / (X_i * dt)``) -- a CFL-style stability bound.
2. Each species' intracellular LP is solved for maximum biomass, yielding a
   growth rate ``mu_i`` and a set of exchange fluxes ``v_(j,i)``.
3. The extracellular ODEs are integrated one Euler step::

       X_i  += (mu_i - death) * X_i * dt
       M_j  += (sum_i v_(j,i) * X_i + influx_j - D * M_j) * dt

   with optional first-order cell death and chemostat dilution ``D``.  Negative
   concentrations are clamped to zero.

This is the formulation the README's "Mathematical Framework" section describes;
the engine is what turns a *steady-state* community model (a snapshot of who
grows on what) into a *time course* of biomass and metabolite concentrations.

Why not MICOM directly?  MICOM computes a single regularised steady-state
community solution -- excellent for one snapshot of cross-feeding, but it is not
a time integrator.  muODE's SOA loop *is* the dynamics; a MICOM/cobra community
model can be plugged in as the per-step solver.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Generator, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from muode.backends import LegacyBackend, SolveRequest, SolverBackend
from muode.community import Community
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.inject import Injection
from muode.kinetics import KineticParameters
from muode.perturb import Perturbation


@dataclass
class SimulationResult:
    """Time-course output of a dynamic-FBA run."""

    times: np.ndarray
    biomass: pd.DataFrame          # index=time, columns=organism ids (gDW/L)
    metabolites: pd.DataFrame      # index=time, columns=metabolite ids (mmol/L)
    growth_rates: pd.DataFrame     # index=time, columns=organism ids (1/h)
    exchange_fluxes: Dict[str, pd.DataFrame] = field(default_factory=dict)
    #: environmental observables from the ecology layer (e.g. pH); may be empty
    environment: Optional[pd.DataFrame] = None
    #: dormant (spore) biomass per species from a life-cycle layer; may be empty
    spores: Optional[pd.DataFrame] = None
    meta: dict = field(default_factory=dict)

    # -- summaries ----------------------------------------------------------
    def final_biomass(self) -> Dict[str, float]:
        return self.biomass.iloc[-1].to_dict()

    def extinct(self, abs_threshold: float = 1e-6, rel_threshold: float = 1e-3) -> List[str]:
        """Species whose final biomass collapsed.

        A species counts as extinct if its final biomass is below
        ``abs_threshold`` (gDW/L) or has fallen below ``rel_threshold`` times its
        initial biomass.
        """
        out = []
        first, last = self.biomass.iloc[0], self.biomass.iloc[-1]
        for sp in self.biomass.columns:
            if last[sp] < abs_threshold or last[sp] < rel_threshold * max(first[sp], 1e-12):
                out.append(sp)
        return out

    def cross_feeding(self, flux_threshold: float = 1e-6) -> pd.DataFrame:
        """Infer cross-feeding interactions from the recorded exchange fluxes.

        Returns one row per ``(producer, metabolite, consumer)`` for which, at
        any recorded time point, the producer secreted a metabolite that the
        consumer simultaneously took up.  The reported strength is the maximum
        co-occurring secretion flux over the trajectory.
        """
        if not self.exchange_fluxes:
            return pd.DataFrame(columns=["producer", "metabolite", "consumer", "strength"])
        species = list(self.exchange_fluxes)
        mets = list(next(iter(self.exchange_fluxes.values())).columns)
        rows = []
        for met in mets:
            for prod in species:
                sec = self.exchange_fluxes[prod][met]
                for cons in species:
                    if cons == prod:
                        continue
                    upt = self.exchange_fluxes[cons][met]
                    coupled = (sec > flux_threshold) & (upt < -flux_threshold)
                    if coupled.any():
                        rows.append(
                            {
                                "producer": prod,
                                "metabolite": met,
                                "consumer": cons,
                                "strength": float(sec[coupled].max()),
                            }
                        )
        return pd.DataFrame(rows, columns=["producer", "metabolite", "consumer", "strength"])

    # -- abundance snapshots ------------------------------------------------
    def abundance_profile(
        self,
        times: Optional[List[float]] = None,
        interval: Optional[float] = None,
    ) -> pd.DataFrame:
        """Relative abundance (sum-to-1) sampled at chosen time points.

        Parameters
        ----------
        times:
            Explicit time points in hours.  Each is matched to the nearest
            recorded simulation step.  Takes precedence over ``interval``.
        interval:
            Sample every ``interval`` hours from 0 to t_end.  Ignored when
            ``times`` is provided.  If neither is given, every recorded step
            is returned.

        Returns
        -------
        DataFrame with time index (h), one column per species, values are
        relative abundances in [0, 1] that sum to 1 at each row (or 0 if
        total biomass is zero at that time point).
        """
        bm = self.biomass
        if times is not None:
            idx = [int(bm.index.get_indexer([t], method="nearest")[0]) for t in times]
            bm = bm.iloc[idx]
        elif interval is not None:
            t_max = float(bm.index[-1])
            sample_times = np.arange(0.0, t_max + interval * 1e-6, interval)
            idx = [int(bm.index.get_indexer([t], method="nearest")[0]) for t in sample_times]
            bm = bm.iloc[idx]
        total = bm.sum(axis=1).replace(0.0, np.nan)
        return bm.div(total, axis=0).fillna(0.0)

    # -- IO -----------------------------------------------------------------
    def to_csv(
        self,
        outdir: str | Path,
        snapshot_interval: Optional[float] = None,
        snapshot_times: Optional[List[float]] = None,
    ) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.biomass.to_csv(outdir / "biomass.csv", index_label="time_h")
        self.metabolites.to_csv(outdir / "metabolites.csv", index_label="time_h")
        self.growth_rates.to_csv(outdir / "growth_rates.csv", index_label="time_h")
        cf = self.cross_feeding()
        if not cf.empty:
            cf.to_csv(outdir / "cross_feeding.csv", index=False)
        if self.environment is not None and not self.environment.empty:
            self.environment.to_csv(outdir / "environment.csv", index_label="time_h")
        if self.spores is not None and not self.spores.empty:
            self.spores.to_csv(outdir / "spores.csv", index_label="time_h")
        if snapshot_interval is not None or snapshot_times is not None:
            snap = self.abundance_profile(times=snapshot_times, interval=snapshot_interval)
            snap.to_csv(outdir / "abundance_snapshots.tsv", sep="\t", index_label="time_h")


# Uptake routes (BiGG) by the element a biomass reaction cannot be built without.
# A diet supplying none of a row, for models that *can* take one up, cannot grow.
_ELEMENTAL_SOURCES: tuple = (
    ("nitrogen", {"nh4_e", "no3_e", "no2_e", "urea_e", "gln__L_e", "glu__L_e",
                  "ala__L_e", "asp__L_e", "asn__L_e", "arg__L_e", "ser__L_e"}),
    ("phosphorus", {"pi_e", "ppi_e", "g3pe_e", "glyc3p_e", "pep_e"}),
    ("sulfur", {"so4_e", "so3_e", "h2s_e", "tsul_e", "cys__L_e", "met__L_e"}),
)


def _warn_if_medium_cannot_feed(organisms, diet: Diet, eco_mets=()) -> None:
    """Warn when the diet cannot actually feed the models.

    A metabolite absent from the diet starts at concentration zero, and a
    zero concentration gives a zero uptake bound -- so an organism can only ever
    consume what the diet names.  A diet whose ids do not match the models'
    exchange namespace (BiGG `glc__D_e` vs a bare `glc_e`), or one missing the
    N/P/S and ion sources a genome-scale biomass reaction needs, therefore yields
    growth = 0 for every species and a perfectly flat, plausible-looking run.
    That silent failure is worse than a loud one, so say it out loud.
    """
    supplied = {
        m for m in diet.metabolites()
        if diet.initial_concentration(m) > 0.0 or diet.influx_rate(m) > 0.0
    }
    if not supplied:
        warnings.warn(
            f"diet '{diet.name}' supplies no nutrient at a nonzero concentration or "
            "influx: every species will have zero uptake and cannot grow.",
            RuntimeWarning,
            stacklevel=3,
        )
        return

    known: set = set()
    for o in organisms:
        known |= set(o.exchange_metabolites())

    # An ecology layer (bile acids, antibiotics, oxygen) may own a metabolite no
    # organism exchanges -- it acts on growth directly, so it is not inert.
    accounted = known | set(eco_mets)

    # (1) Namespace mismatch: diet ids that name no exchange in any model. These
    # are silently inert, so a diet can look full and feed nothing.
    unmatched = supplied - accounted
    if unmatched and len(unmatched) == len(supplied):
        warnings.warn(
            f"diet '{diet.name}': none of its {len(supplied)} nutrients match an exchange "
            "in any model -- the diet and the models are probably in different namespaces "
            f"(CarveMe emits BiGG, e.g. 'glc__D_e'). Unmatched: {sorted(unmatched)[:5]}...",
            RuntimeWarning,
            stacklevel=3,
        )
    elif unmatched:
        warnings.warn(
            f"diet '{diet.name}': {len(unmatched)} of {len(supplied)} nutrients match no "
            f"exchange in any model and are inert: {sorted(unmatched)[:5]}"
            f"{'...' if len(unmatched) > 5 else ''}",
            RuntimeWarning,
            stacklevel=3,
        )

    # (2) Missing elemental source.  Biomass is not made of carbon alone: with no
    # N, P or S source the LP is feasible but its objective is 0, so every species
    # sits at exactly its seeded biomass and the run looks "successful" while
    # being flat.  Only complain when the models *do* have an uptake route for
    # that element -- otherwise this is a toy model with a lumped biomass, not a
    # broken diet.  (Note we deliberately do NOT flag a species that shares no
    # metabolite with the diet: an obligate cross-feeder is fed by the community,
    # not by the medium, and that is legitimate.)
    for element, sources in _ELEMENTAL_SOURCES:
        routes = sources & known
        if routes and not (sources & supplied):
            warnings.warn(
                f"diet '{diet.name}' supplies no {element} source, but the models can take "
                f"one up ({', '.join(sorted(routes)[:4])}). A genome-scale biomass reaction "
                f"cannot fire without {element}: expect zero growth for every species and a "
                "flat run. Add one to the diet.",
                RuntimeWarning,
                stacklevel=3,
            )


@dataclass
class DynamicFBA:
    """Configurable dynamic-FBA community integrator.

    Parameters
    ----------
    t_end:
        Total simulated time (h).
    dt:
        Integration step (h).  Smaller is more accurate and reduces the chance
        of substrate over-depletion warnings.
    death_rate:
        First-order biomass decay (1/h), applied to every species.
    dilution_rate:
        Chemostat dilution ``D`` (1/h); also washes out biomass and metabolites.
    min_biomass:
        Biomass floor (gDW/L); species below it stop contributing flux (a numeric
        proxy for local extinction) but can recover if it rises again.
    record_fluxes:
        Keep per-species exchange-flux trajectories (needed for cross-feeding
        inference).  Disable for very large communities to save memory.
    n_jobs:
        Number of worker threads for the per-step FBA solves.  ``1`` (default)
        runs sequentially; ``-1`` uses all available cores.  The species' solves
        are independent within a step, so any value yields identical results --
        the speed-up is realised on genome-scale models, whose solver calls
        dominate each step and release the GIL.  Only used by the default
        (legacy) backend.
    backend:
        The :class:`~muode.backends.SolverBackend` that solves each step's LPs
        (``None``: the historical per-organism path, bit-identical to earlier
        muODE).  The engine hands the backend *every* LP of a step at once --
        all species, and in :meth:`run_ensemble` all ensemble members -- which is
        what a batched (e.g. GPU) solver needs.
    """

    t_end: float = 24.0
    dt: float = 0.1
    death_rate: float = 0.0
    dilution_rate: float = 0.0
    min_biomass: float = 1e-9
    record_fluxes: bool = True
    n_jobs: int = 1
    backend: Optional[SolverBackend] = None

    # -- drivers ------------------------------------------------------------------
    def run(
        self,
        community: Community,
        diet: Diet,
        kinetics: Optional[KineticParameters] = None,
        perturbation: Optional[Perturbation] = None,
        injections: Optional[Sequence[Injection]] = None,
        ecology: Optional[EcologyModel] = None,
        backend: Optional[SolverBackend] = None,
    ) -> SimulationResult:
        """Integrate one trajectory (see the class docstring)."""
        # Apply any perturbation to the baseline bounds *once*, up front, so it
        # persists across every step (an antibiotic does not wear off mid-run).
        if perturbation is not None:
            perturbation.apply(community)
        gen = self._integrate(community, diet, kinetics, perturbation, injections, ecology, member=0)
        be, owned = self._resolve_backend(backend)
        try:
            return _drive([gen], be)[0]
        finally:
            if owned:
                be.close()

    def run_ensemble(
        self,
        community: Community,
        diets: Union[Diet, Sequence[Diet]],
        kinetics: Union[None, KineticParameters, Sequence[KineticParameters]] = None,
        n_members: Optional[int] = None,
        perturbation: Optional[Perturbation] = None,
        injections: Optional[Sequence[Injection]] = None,
        ecology: Optional[Callable[[], EcologyModel]] = None,
        backend: Optional[SolverBackend] = None,
    ) -> List[SimulationResult]:
        """Integrate an ensemble of trajectories *in lockstep*.

        Members share the community's organism models (and any ``perturbation``,
        applied once) and differ in diet and/or kinetics -- the parameter-sweep /
        Monte-Carlo setting.  Every step, the LPs of all species of all members
        go to the backend as one batch.  Each member's result is identical to
        running :meth:`run` on it alone (the members never interact).

        ``diets`` / ``kinetics``: one object shared by all members, or one per
        member.  ``ecology``: a zero-argument factory, since ecology layers carry
        per-trajectory state.
        """
        diets_l = list(diets) if isinstance(diets, (list, tuple)) else None
        kin_l = list(kinetics) if isinstance(kinetics, (list, tuple)) else None
        E = n_members or len(diets_l or kin_l or [0])
        if E < 1:
            raise ValueError("ensemble needs at least one member")
        if perturbation is not None:
            perturbation.apply(community)
        gens = []
        for e in range(E):
            d = diets_l[e] if diets_l is not None else diets
            k = kin_l[e] if kin_l is not None else kinetics
            eco = ecology() if ecology is not None else None
            gens.append(self._integrate(community, d, k, perturbation, injections, eco, member=e))
        be, owned = self._resolve_backend(backend)
        try:
            return _drive(gens, be)
        finally:
            if owned:
                be.close()

    def _resolve_backend(self, backend: Optional[SolverBackend]):
        be = backend or self.backend
        if be is not None:
            return be, False
        return LegacyBackend(n_jobs=self.n_jobs), True

    # -- the integrator (a coroutine: yields each step's LPs, receives solutions) ----
    def _integrate(
        self,
        community: Community,
        diet: Diet,
        kinetics: Optional[KineticParameters],
        perturbation: Optional[Perturbation],
        injections: Optional[Sequence[Injection]],
        ecology: Optional[EcologyModel],
        member: int,
    ) -> Generator[List[SolveRequest], list, SimulationResult]:
        kinetics = kinetics or KineticParameters()

        # Timed biomass injections (transplants / probiotic doses).  They are
        # state events, not bound changes: each fires at the first step reaching
        # its time, adding biomass to community members (see muode.inject).
        injections = list(injections or [])
        injected = [False] * len(injections)

        # Optional environmental / life-history layer (pH, bile acids, spores,
        # antibiotic PK, bacteriocins).  An empty model is a no-op, leaving the
        # run numerically identical to the bare engine (see muode.ecology).
        ecology = ecology or EcologyModel()
        ecology.reset(community)
        record_eco = not ecology.is_empty()

        organisms = community.organisms
        env_mets = sorted(
            set(community.environment_metabolites())
            | set(diet.metabolites())
            | set(ecology.extra_metabolites())
        )

        # state vectors
        X = dict(community.initial_biomass())
        eco_init = ecology.initial_concentrations()
        M = {}
        for m in env_mets:
            c = diet.initial_concentration(m)
            if c == 0.0 and m in eco_init:
                c = float(eco_init[m])
            M[m] = c

        _warn_if_medium_cannot_feed(organisms, diet, ecology.extra_metabolites())

        n_steps = int(round(self.t_end / self.dt))
        times = np.linspace(0.0, n_steps * self.dt, n_steps + 1)

        # history buffers
        bio_hist = {o.id: np.empty(n_steps + 1) for o in organisms}
        mu_hist = {o.id: np.empty(n_steps + 1) for o in organisms}
        met_hist = {m: np.empty(n_steps + 1) for m in env_mets}
        flux_hist: Dict[str, Dict[str, np.ndarray]] = (
            {o.id: {m: np.zeros(n_steps + 1) for m in env_mets} for o in organisms}
            if self.record_fluxes
            else {}
        )
        obs_hist: List[dict] = []      # ecology observables per step (e.g. pH)
        spore_hist: List[dict] = []    # dormant biomass per species per step

        depletion_warned = False
        anything_grew = False

        # Fast path: without ecology layers every ecology hook is the identity, so the
        # per-metabolite Python loops below can run as NumPy vector operations.  They
        # perform the same IEEE operations in the same order, so results are bit-for-bit
        # identical to the scalar path; metabolites an organism cannot exchange are
        # skipped (``set_uptake_bound`` ignores them anyway).
        fast = ecology.is_empty()
        if fast:
            env_index = {m: i for i, m in enumerate(env_mets)}
            org_mets, org_idx, org_vmax, org_km = [], [], [], []
            for o in organisms:
                mets = tuple(m for m in o.exchange_metabolites() if m in env_index)
                org_mets.append(mets)
                org_idx.append(np.array([env_index[m] for m in mets], dtype=np.int64))
                pars = [kinetics.get(o.id, m) for m in mets]
                org_vmax.append(np.array([p[0] for p in pars], dtype=float))
                org_km.append(np.array([p[1] for p in pars], dtype=float))
            lim_all = np.array([np.inf if diet.uptake_limit(m) is None else diet.uptake_limit(m)
                                for m in env_mets], dtype=float)
            influx_all = np.array([diet.influx_rate(m) for m in env_mets], dtype=float)
            Marr = np.array([M[m] for m in env_mets], dtype=float)
            met_arr = np.empty((n_steps + 1, len(env_mets)))
            flux_arr = ({o.id: np.zeros((n_steps + 1, len(env_mets))) for o in organisms}
                        if self.record_fluxes else {})
            org_pos = {o.id: k for k, o in enumerate(organisms)}

        for step in range(n_steps + 1):
            t_now = times[step]

            # --- fire any pending injections (bolus appears at this time) ----
            # Done before recording/solving so the introduced biomass shows in
            # the time course and is metabolically active from this step on.
            for k, inj in enumerate(injections):
                if not injected[k] and t_now + 1e-9 >= inj.time:
                    for oid, amount in inj.biomass.items():
                        if oid in X:
                            X[oid] = X[oid] + float(amount)
                    injected[k] = True

            # record current state
            for o in organisms:
                bio_hist[o.id][step] = X[o.id]
            if fast:
                met_arr[step] = Marr
            else:
                for m in env_mets:
                    met_hist[m][step] = M[m]
            if record_eco:
                obs_hist.append(ecology.observables(t_now, M, X))
                spore_hist.append(dict(ecology.spores()))

            # --- solve each species' FBA under current medium ---------------
            # Uptake bound of every environment metabolite for every active
            # species; the LPs themselves go to the backend as one batch.
            step_growth: Dict[str, float] = {}
            step_flux: Dict[str, Dict[str, float]] = {}

            # Species below the biomass floor are dormant: no flux, no solve.
            active = [o for o in organisms if X[o.id] > self.min_biomass]
            for o in organisms:
                if X[o.id] <= self.min_biomass:
                    step_growth[o.id] = 0.0
                    if not fast:
                        step_flux[o.id] = {m: 0.0 for m in env_mets}
                    mu_hist[o.id][step] = 0.0

            if fast:
                F = np.zeros((len(organisms), len(env_mets)))
                requests = []
                for o in active:
                    k = org_pos[o.id]
                    idx = org_idx[k]
                    conc = Marr[idx]
                    vmax, km = org_vmax[k], org_km[k]
                    with np.errstate(divide="ignore", invalid="ignore"):
                        mm = np.where(conc <= 0.0, 0.0, vmax * conc / (km + conc))
                        mm = np.minimum(mm, lim_all[idx])
                        cap = conc / (X[o.id] * self.dt) if self.dt > 0 else np.full(conc.shape, np.inf)
                    requests.append(SolveRequest(o, None, member, mets=org_mets[k], vec=np.minimum(mm, cap)))
                sols = yield requests
                for o, sol in zip(active, sols):
                    k = org_pos[o.id]
                    raw = sol.growth_rate if sol.feasible else 0.0
                    step_growth[o.id] = raw
                    mu_hist[o.id][step] = raw
                    ex = sol.exchange_fluxes
                    F[k, org_idx[k]] = [ex.get(m, 0.0) for m in org_mets[k]]
                if any(step_growth[o.id] > 0.0 for o in active):
                    anything_grew = True
                if self.record_fluxes:
                    for o in organisms:
                        flux_arr[o.id][step] = F[org_pos[o.id]]
                if step == n_steps:
                    break
                for o in organisms:
                    mu = step_growth[o.id]
                    death = self.death_rate + self.dilution_rate + ecology.extra_death(o.id, t_now, M, X)
                    X[o.id] = max(0.0, X[o.id] + (mu - death) * X[o.id] * self.dt)
                ecology.integrate(t_now, self.dt, X, M, step_growth)
                dM = influx_all - self.dilution_rate * Marr + 0.0
                for o in organisms:
                    dM += F[org_pos[o.id]] * X[o.id]
                new = Marr + dM * self.dt
                neg = new < 0.0
                if neg.any():
                    new[neg] = 0.0
                    if not depletion_warned:
                        warnings.warn(
                            "A metabolite was fully depleted within a step; consider a "
                            "smaller dt for better accuracy.",
                            RuntimeWarning,
                            stacklevel=2,
                        )
                        depletion_warned = True
                Marr = new
                continue

            requests = []
            for o in active:
                uptake = {}
                for m in env_mets:
                    conc = M[m]
                    mm = kinetics.michaelis_menten(o.id, m, conc)
                    # what the *diet* makes available, if it says.  Without this a
                    # large medium feeds every nutrient at the same uniform Vmax and
                    # the cell grows faster than any organism ever measured.
                    limit = diet.uptake_limit(m)
                    if limit is not None:
                        mm = min(mm, limit)
                    mm *= ecology.uptake_factor(o.id, m, t_now, M, X)
                    # CFL-style cap: do not let one species take more than exists.
                    cap = conc / (X[o.id] * self.dt) if self.dt > 0 else np.inf
                    uptake[m] = min(mm, cap)
                requests.append(SolveRequest(o, uptake, member))
            sols = yield requests

            solved = []
            for o, sol in zip(active, sols):
                raw = sol.growth_rate if sol.feasible else 0.0
                # environmental growth modifiers (pH, bile, bacteriocins, ...)
                mu_eff = raw * ecology.growth_factor(o.id, t_now, M, X)
                flux = {m: sol.exchange_fluxes.get(m, 0.0) for m in env_mets}
                solved.append((o, (mu_eff, flux)))

            for o, (mu_eff, flux) in solved:
                step_growth[o.id] = mu_eff
                step_flux[o.id] = flux
                mu_hist[o.id][step] = mu_eff

            if any(step_growth[o.id] > 0.0 for o in active):
                anything_grew = True

            if self.record_fluxes:
                for o in organisms:
                    for m in env_mets:
                        flux_hist[o.id][m][step] = step_flux[o.id][m]

            if step == n_steps:
                break  # state recorded; no integration past the horizon

            # --- integrate one Euler step -----------------------------------
            for o in organisms:
                mu = step_growth[o.id]
                death = self.death_rate + self.dilution_rate + ecology.extra_death(o.id, t_now, M, X)
                X[o.id] = max(
                    0.0,
                    X[o.id] + (mu - death) * X[o.id] * self.dt,
                )

            # life-history transitions (e.g. vegetative <-> spore) may move
            # biomass between the live pool X and an internal dormant pool.
            ecology.integrate(t_now, self.dt, X, M, step_growth)

            eco_rates = ecology.metabolite_rates(t_now, M, X)
            for m in env_mets:
                dM = diet.influx_rate(m) - self.dilution_rate * M[m] + eco_rates.get(m, 0.0)
                for o in organisms:
                    dM += step_flux[o.id][m] * X[o.id]
                new = M[m] + dM * self.dt
                if new < 0.0:
                    new = 0.0
                    if not depletion_warned:
                        warnings.warn(
                            "A metabolite was fully depleted within a step; consider a "
                            "smaller dt for better accuracy.",
                            RuntimeWarning,
                            stacklevel=2,
                        )
                        depletion_warned = True
                M[m] = new

        # Not one species grew at any point in the whole run: biomass is exactly
        # what was seeded, and the output is a flat line that looks like a
        # successful simulation. Nearly always a setup error (a medium that does
        # not feed these models, a namespace mismatch, models that never
        # gap-filled) rather than biology, so refuse to let it pass silently.
        if organisms and not anything_grew:
            warnings.warn(
                f"no species grew at any point on diet '{diet.name}': every growth rate was "
                "zero for the entire run, so biomass is flat and identical to the inoculum. "
                "Check that the diet feeds these models -- ids in the same namespace (CarveMe "
                "emits BiGG, e.g. 'glc__D_e'), and a nitrogen, phosphate and sulfur source "
                "present -- and that the models grow standalone.",
                RuntimeWarning,
                stacklevel=2,
            )

        if fast:
            for i, m in enumerate(env_mets):
                met_hist[m] = met_arr[:, i].copy()
                M[m] = float(Marr[i])
            if self.record_fluxes:
                for o in organisms:
                    for i, m in enumerate(env_mets):
                        flux_hist[o.id][m] = flux_arr[o.id][:, i].copy()

        # --- assemble result -------------------------------------------------
        biomass_df = pd.DataFrame(bio_hist, index=times)
        mu_df = pd.DataFrame(mu_hist, index=times)
        met_df = pd.DataFrame(met_hist, index=times)
        exch = (
            {oid: pd.DataFrame(flux_hist[oid], index=times) for oid in flux_hist}
            if self.record_fluxes
            else {}
        )
        env_df = (
            pd.DataFrame(obs_hist, index=times) if record_eco and obs_hist and obs_hist[0] else None
        )
        spore_df = (
            pd.DataFrame(spore_hist, index=times) if record_eco and spore_hist and spore_hist[0] else None
        )
        return SimulationResult(
            times=times,
            biomass=biomass_df,
            metabolites=met_df,
            growth_rates=mu_df,
            exchange_fluxes=exch,
            environment=env_df,
            spores=spore_df,
            meta={
                "t_end": self.t_end,
                "dt": self.dt,
                "death_rate": self.death_rate,
                "dilution_rate": self.dilution_rate,
                "diet": diet.name,
                "n_species": len(organisms),
                "perturbation": None if perturbation is None else perturbation.describe(),
                "injections": [inj.describe() for inj in injections] or None,
                "ecology": ecology.describe() or None,
            },
        )


def _drive(gens: List[Generator], backend: SolverBackend) -> List[SimulationResult]:
    """Advance integrator coroutines in lockstep, one batched solve per step."""
    results: List[Optional[SimulationResult]] = [None] * len(gens)
    pending: Dict[int, List[SolveRequest]] = {}
    for i, g in enumerate(gens):
        try:
            pending[i] = next(g)
        except StopIteration as stop:
            results[i] = stop.value
    while pending:
        order = list(pending)
        flat = [r for i in order for r in pending[i]]
        sols = backend.solve_step(flat) if flat else []
        pos = 0
        nxt: Dict[int, List[SolveRequest]] = {}
        for i in order:
            k = len(pending[i])
            try:
                nxt[i] = gens[i].send(sols[pos:pos + k])
            except StopIteration as stop:
                results[i] = stop.value
            pos += k
        pending = nxt
    return results  # type: ignore[return-value]
