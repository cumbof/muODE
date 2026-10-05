# µODE — Dynamic-FBA engine

The engine is the scientific core of µODE. It implements **dynamic Flux Balance
Analysis (dFBA)** through the **Static Optimization Approach** (SOA; Mahadevan,
Edwards & Doyle, *Biophys. J.* 2002), extended to multi-species communities that
share one extracellular metabolite pool.

---

## Conceptual overview

Ordinary FBA solves a linear program at a single snapshot in time — one growth
condition, no dynamics. Dynamic FBA connects a series of such snapshots into a
trajectory: at each discrete time step the engine

1. sets the substrate-uptake bounds from the current extracellular concentrations
   using Michaelis–Menten kinetics,
2. solves each species' FBA independently (one LP per organism),
3. integrates the extracellular ODEs one step forward with the resulting flux
   vectors (Euler method),
4. repeats.

What makes this a *community* simulation is that every species draws from and
secretes into the **same** extracellular pool, so cross-feeding metabolites
produced by species A are immediately available to species B.

---

## Mathematical framework

### Per-species FBA (intracellular)

At each time step, for every organism *i*, the engine solves:

$$\text{Maximize } Z_i = c_i^T v_i$$
$$\text{Subject to } S_i \cdot v_i = 0, \quad v_{i,\min} \le v_i \le v_{i,\max}$$

where $S_i$ is the stoichiometric matrix, $v_i$ the flux vector, and $c_i$ the
biomass objective. The solution yields the growth rate $\mu_i$ (the LP objective
value) and the exchange fluxes $v_{j,i}$ of every extracellular metabolite $j$.

### Michaelis–Menten uptake bounds

The upper bound of each substrate-uptake reaction is set from the extracellular
concentration $M_j$ and the per-species kinetic parameters:

$$v^{up}_{j,i} = \min\!\left(\frac{V_{\max}\, M_j}{K_m + M_j},\; \frac{M_j}{X_i \cdot \Delta t}\right)$$

The second term is a **CFL-style stability cap** that prevents a single step from
over-depleting a metabolite below zero. When no kinetics are supplied µODE falls
back to safe defaults ($V_{\max} = 10$, $K_m = 0.01$ mmol/gDW/h), which recover
substrate-unlimited FBA for abundant metabolites.

### Extracellular ODEs

After solving all species' FBAs, the environment evolves:

$$\frac{dX_i}{dt} = (\mu_i - \delta - D)\, X_i$$

$$\frac{dM_j}{dt} = \sum_{i} v_{j,i}\, X_i + \phi_j - D\, M_j$$

- $X_i$ — biomass of species $i$ (gDW/L)
- $M_j$ — extracellular concentration of metabolite $j$ (mmol/L)
- $\delta$ — first-order death rate (1/h; default 0)
- $D$ — chemostat dilution rate (1/h; default 0, i.e. batch)
- $\phi_j$ — nutrient influx from the diet (mmol/L/h)

Concentrations are clamped to $\ge 0$ after each step.

---

## Solver abstraction

The engine never imports a metabolic-modeling library directly. Instead it drives
the `OrganismModel` protocol — a minimal interface with four methods:

| Method | What it does |
|--------|-------------|
| `exchange_metabolites()` | Returns the list of extracellular metabolite ids |
| `set_uptake_bound(met_id, ub)` | Sets the upper bound of one exchange reaction |
| `reset_bounds()` | Restores all bounds to their pre-step values |
| `optimize()` | Solves the LP; returns `(growth_rate, {met_id: flux})` |

Two backends implement this protocol:

### `CobraOrganism` — real GEMs

Wraps a `cobra.Model` loaded from an SBML file. Used whenever genome-scale
metabolic models are available (CarveMe output, curated models, etc.). Requires
COBRApy and a solver (HiGHS via `scipy`/`optlang` by default; Gurobi or CPLEX are
drop-in alternatives for large communities).

```python
from muode.organism import CobraOrganism
org = CobraOrganism.from_file("path/to/model.xml", id="Bacteroides_fragilis")
```

### `LinprogOrganism` — dependency-light backend

A self-contained stoichiometric model solved with `scipy.optimize.linprog`
(HiGHS backend). It carries its own stoichiometry and bounds as plain NumPy
arrays — no SBML, no COBRApy — so the engine runs on any machine. It backs the
built-in demo and the offline reconstruction engine.

```python
from muode.examples import build_toy_community
comm = build_toy_community()   # a Community of LinprogOrganisms
```

---

## The `DynamicFBA` class

```python
from muode.dfba import DynamicFBA, SimulationResult

engine = DynamicFBA(
    t_end=24.0,       # simulated hours
    dt=0.1,           # Euler step size (h)
    death_rate=0.0,   # first-order death (1/h)
    dilution_rate=0.0,# chemostat D (1/h)
    n_jobs=1,         # per-step solver threads (-1 = all cores)
)
result: SimulationResult = engine.run(community, diet, kinetics,
                                      perturbation=None,  # bound changes (see perturbation.md)
                                      injections=None)    # timed biomass events (see injection.md)
```

`perturbation=` applies bound changes once up front; `injections=` is a list of
timed [`Injection`](injection.md) state events (transplants, probiotic doses) that
introduce biomass mid-run.

### Parallel per-step solves (`n_jobs`)

Within each time step the species' FBA problems are independent: each reads only
the shared medium and biomass snapshots and writes only its own bounds, so they
can be solved concurrently. `n_jobs` sets the number of worker threads — `1`
(default) is sequential, `-1` uses all cores. Results are keyed by species and
order-independent, so **the trajectory is identical for any `n_jobs`**; the
speed-up is realized on genome-scale models, whose solver calls dominate each step
and release the GIL (a moderate LP already runs ~2.4× faster on 4 threads). From
the CLI: `muode simulate --community community.json -j -1`.

### `SimulationResult`

| Attribute | Type | Content |
|-----------|------|---------|
| `result.biomass` | `pd.DataFrame` | rows = time points, columns = species ids |
| `result.metabolites` | `pd.DataFrame` | rows = time points, columns = metabolite ids |
| `result.growth_rates` | `pd.DataFrame` | rows = time points, columns = species ids |
| `result.meta` | `dict` | run parameters (t_end, dt, n_steps, ...) |

```python
result.to_csv("results/")           # writes biomass.csv, metabolites.csv, etc.
result.cross_feeding()              # DataFrame of producer→metabolite→consumer edges
result.extinct()                    # list of species that reached zero biomass
```

---

## Cross-feeding inference

After the run, `result.cross_feeding()` analyzes the time series: a
producer→consumer edge exists when species A has a positive net secretion flux
for some metabolite and species B has a positive uptake flux for the same
metabolite. The `strength` column is the mean product of their fluxes, a proxy for
interaction intensity.

This is an *emergent* property. The engine never hard-codes interactions; they
appear because the shared metabolite pool creates a genuine dependency.

---

## Operational notes

- **Integration step.** The SOA uses first-order (Euler) integration; use
  `dt ≤ 0.1 h` for typical gut communities, where the per-step accuracy is
  `O(dt)`.
- **Pseudo-steady-state fluxes.** The SOA assumes intracellular fluxes reach
  steady state within each step — the standard assumption in community dFBA,
  appropriate when the intracellular timescale is much faster than the
  extracellular one.
- **COMETS interoperability.** A community can be exported to COMETS (the
  independent dynamic community-FBA engine) for cross-validation via
  `muode export-comets` / `muode.comets.export_comets`, which writes a layout,
  parameters and a `cometspy` driver. The `OrganismModel` protocol is the
  integration point for additional per-step solvers.

---

## References

- Mahadevan, Edwards & Doyle (2002) *Dynamic flux balance analysis of diauxic
  growth in Escherichia coli.* Biophys. J. 83(3):1331–1340.
