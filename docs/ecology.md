# µODE — Ecology Layer (pH, bile acids, spores, antibiotics, antagonism)

The dynamic-FBA core ([engine.md](engine.md)) models **metabolism**: who eats
what, and the cross-feeding that emerges from a shared metabolite pool. Real
communities are also shaped by processes that are *not* metabolic fluxes — the
chemistry of the medium, chemical transformation of the pool, pharmacology,
direct antagonism, and dormancy. The **ecology layer** adds these as composable,
optional modifiers on top of the FBA core.

> An empty ecology model is a strict no-op: the run is numerically identical to
> the bare engine. Layers only ever *add* mechanisms.

See [LIMITATIONS.md](LIMITATIONS.md) for how these fit the full picture of what
muODE can and cannot model, and `examples/fmt_cdiff/` for an end-to-end use.

---

## How it plugs in

An `EcologyModel` is a list of `EcologyLayer` objects. Each layer overrides only
the hooks it needs, and the dynamic-FBA loop calls them at defined points:

| Hook | Effect in the loop |
|------|--------------------|
| `extra_metabolites` / `initial_concentrations` | track pool species no GEM exchanges (bile acids, a toxin) and seed them |
| `uptake_factor(org, met, t, M, X)` | multiplies a substrate uptake bound (e.g. O₂ gating) |
| `growth_factor(org, t, M, X)` | multiplies a species' growth rate (pH, bile, bacteriocin) |
| `extra_death(org, t, M, X)` | adds to a species' death rate (antibiotic killing) |
| `metabolite_rates(t, M, X)` | extra dM/dt for the pool (bile transformation, toxin secretion) |
| `integrate(t, dt, X, M, mu)` | moves biomass between pools (vegetative ↔ spore) |
| `observables` / `spores` | recorded into `result.environment` / `result.spores` |

```python
from muode import DynamicFBA, EcologyModel, WeakAcidInhibition

eco = EcologyModel([WeakAcidInhibition()])
result = DynamicFBA(t_end=48).run(community, diet, kinetics, ecology=eco)
result.environment["pH"]      # recorded pH trajectory
```

Growth/uptake factors **multiply**, death and metabolite rates **sum**, so layers
compose without ordering surprises.

---

## Layers

### `WeakAcidInhibition` — pH and SCFA toxicity (`muode.ph`)

Short-chain fatty acids are weak acids. As they accumulate they (1) acidify the
medium — buffered linear titration `pH = pH0 − total_acid / buffer_capacity` —
and (2) inhibit growth, because only the **undissociated** acid
(`HA = A / (1 + 10^(pH − pKa))`, Henderson–Hasselbalch) crosses the membrane and
dissipates the proton-motive force. Growth is scaled by
`1/(1 + ΣHA/Ki)` times a Rosso cardinal-pH window. Per-species `Ki` makes an
acid-sensitive pathogen (low `Ki`) suppressible by a dense SCFA producer.

```python
WeakAcidInhibition(ph0=6.8, buffer_capacity=25.0,
                   default_ki=12.0, ki={"C_difficile": 3.0})
```

This makes fermentation **self-limiting** and is the metabolic arm of
colonization resistance. Records `pH` and `undissociated_acid`.

### `BileAcidTransform` + `BileAcidInhibition` — bile metabolism (`muode.bile`)

Two enzymatic steps transform the bile-acid pool, each Michaelis–Menten in the
substrate and first-order in the carrying guild's biomass:

```
taurocholate --BSH (broad)----------> cholate
cholate      --7α-dehydroxylase (rare, bai operon)--> deoxycholate
```

`BileAcidTransform(bsh_producers=…, bai_producers=…)` evolves the pool;
`BileAcidInhibition(targets=…, ki=…)` applies secondary-bile-acid toxicity to
target species. Taurocholate is the *C. difficile* germinant; deoxycholate
inhibits it — the bile arm of colonization resistance and the reason a *bai*-
carrying donor community resists the pathogen.

### `SporeForming` — dormancy life cycle (`muode.lifecycle`)

Splits a spore-former's biomass into the vegetative pool (`X`, grown and killed
as usual) and an internal **spore** pool that does not grow and is **not killed**
by antibiotics. Per step: vegetative cells **sporulate** under nutrient stress
(growth below `mu_stress`), and spores **germinate** at a rate gated by a
germinant (taurocholate, Michaelis–Menten) and **blocked** by secondary bile
acids. This is what makes infection recur (spores survive a drug course) and what
an FMT shuts off (restored bile metabolism keeps germination blocked). Records
`germination_signal` and per-species spore biomass.

### `Antibiotic` — pharmacokinetics / dynamics (`muode.antibiotic`)

A drug is not a constant. One-compartment PK (`C(t) = Σ dose·e^{−k(t−t_d)}` over
the dosing schedule) feeds an Emax kill model `Emax·C/(EC50+C)` on the
**vegetative** biomass of susceptible species. Because spores live in a separate
pool, they survive — the mechanistic basis of recurrence and of why a drug alone
fails where a drug **+** FMT succeeds.

```python
Antibiotic(susceptible={"C_difficile"}, dose_times=(0,2,4,6,8),
           dose=3.0, half_life=2.0, emax=5.0, ec50=0.4)
```

### `Bacteriocin` — interference competition (`muode.antagonism`)

A diffusible toxin produced first-order in `producers` biomass, decaying
first-order, suppressing `targets` growth by `1/(1 + C/ki)`. Models bacteriocin /
antimicrobial-peptide antagonism that is not mediated by shared nutrients.

---

## Outputs

When an ecology model is active, `SimulationResult` gains two optional frames,
also written by `to_csv`:

- `result.environment` → `environment.csv` — observables over time (pH,
  undissociated acid, germination signal, drug concentration, …).
- `result.spores` → `spores.csv` — dormant biomass per species over time.

`result.meta["ecology"]` lists the active layers.

---

## A complete worked example

`muode.scenarios.cdi_scenario(fmt=…)` composes **all** of the above into the
recurrent-*C. difficile* / FMT story on dependency-light toy models (runs
anywhere). The only difference between the two arms is a timed donor
[injection](injection.md); the recurrence-vs-clearance outcome is emergent:

```bash
python examples/fmt_cdiff/mechanistic_demo.py
```

The genome-scale version is the rest of `examples/fmt_cdiff/`.
