# µODE — Architecture

muODE is two things that fit together:

1. **`muode` — an installable Python library + CLI** holding the science: the
   dynamic-FBA engine, the community container, the perturbation engine, kinetics,
   diets, QC and the external-tool wrappers.
2. **`workflow/` — a Snakemake pipeline** that drives the library at scale,
   fanning the per-MAG phases across a cluster and chaining them into the
   integrative simulation.

This mirrors the `nf-core` / `snakemake-workflows` convention: a reusable library
of logic + a thin, reproducible orchestration layer. Either can be used alone.

```
                 ┌──────────────────────── workflow/ (Snakemake) ───────────────────────┐
   MAGs (*.fna) ─┤ checkm2 ─▶ select_mags(checkpoint) ─▶ reconstruct ─▶ refine ─▶ assemble ─▶ simulate ─▶ perturb
                 └───────────────┬───────────────────────────┬─────────────┬──────────────┬──────────────────────┘
                                 │ calls                      │ calls       │ calls        │ calls
                 ┌───────────────▼────────────────────────────▼─────────────▼──────────────▼──────────┐
                 │                              muode (library + CLI)                                  │
                 │  qc · reconstruct · gapfill · kinetics · community · diet · dfba · perturb · viz     │
                 └──────────────────────────────────────────────────────────────────────────────────┘
```

## The solver abstraction (why both backends exist)

The engine never imports a modelling library directly; it drives the
`OrganismModel` protocol (`exchange_metabolites`, `set_uptake_bound`,
`reset_bounds`, `optimize`). Two backends implement it:

- **`CobraOrganism`** — wraps a genome-scale `cobra.Model` (production, real MAGs).
- **`LinprogOrganism`** — a self-contained stoichiometric model solved with
  `scipy.optimize.linprog` (HiGHS). Zero heavy dependencies, so the engine and the
  whole test-suite run anywhere, and the bundled toy example needs no GEM files.

This is what lets the core be tested and demonstrated without a solver licence,
while the same engine runs real BiGG-namespace community models unchanged.

## The dynamic-FBA loop (`muode.dfba`)

Static Optimization Approach (Mahadevan 2002), per step `dt`:

1. For each species *i* and metabolite *j*: uptake bound = `min(Vmax·M/(Km+M),
   M/(X_i·dt))` — Michaelis–Menten capped by availability (a CFL-style bound).
2. Solve each species' FBA for `µ_i` and exchange fluxes `v_(j,i)`.
3. Integrate `X_i += (µ_i−death−D)·X_i·dt` and
   `M_j += (Σ_i v_(j,i)·X_i + influx_j − D·M_j)·dt`; clamp `M_j ≥ 0`.

Species share one extracellular pool (the COMETS-style compartmentalised
formulation), which is what makes cross-feeding emerge and scales by adding
species rather than rebuilding one giant LP.

## Data contracts (library ⇄ workflow)

- **Abundance:** 2-column TSV `mag_id, abundance` (e.g. a MetaSBT profile).
- **Diet:** preset name or CSV `metabolite, concentration[, influx]`.
- **Per-MAG QC** (`{mag}.qc.json`): gap-fill + sanity-check report emitted by
  `refine`, one per model.
- **Per-MAG kinetics** (`{mag}.kinetics.json`): predicted Km (per uptake
  metabolite) + kcat (per internal reaction), emitted by `refine`; consumed by
  `simulate --kinetics` for Michaelis–Menten bounds and (with
  `--enzyme-constraints`) GECKO-lite kcat caps. Serialised via
  `KineticParameters.to_json`.
- **Reconstruction summary** (`reconstruction_summary.tsv`): one row per MAG
  (`grows_now`, `n_reactions_added`, `energy_generating_cycle`, `simulatable`, …),
  aggregated by `reconstruct_report`. `assemble` reads its `simulatable` column to
  drop pathological models, so a single bad MAG never aborts a large run.
- **Community manifest** (`community.json`): `{models[], abundances{}, diet,
  total_biomass}` — produced by `assemble`, consumed by `simulate`/`perturb`.
- **Results:** `biomass.csv`, `metabolites.csv`, `growth_rates.csv`,
  `cross_feeding.csv`, `meta.json`, and PNG figures.
- **Abundance (MetaSBT):** a MetaSBT profile is ingested by `muode.metasbt`
  (auto-detected id/abundance/taxonomy columns) into the same `{mag_id:
  rel_abundance}` contract; `assemble --metasbt` does this in the workflow.
- **Benchmark expectation** (YAML, `examples/benchmarks/`): expected
  `relative_abundances`, `metabolites`, `cross_feeding` edges + tolerances.
  Consumed by `muode validate` / the `validate` rule, which emits
  `validation/report.json` (per-component metrics + overall pass/fail).
- **Spatial results** (`muode spatial`): `spatial.npz` (per-species biomass and
  per-metabolite fields, shape `n_frames × ny × nx`), `spatial_total_biomass.csv`
  and PNG fields/time-course. The spatial engine reuses the same
  `OrganismModel`/`Community`/`KineticParameters`, solving FBA per grid cell.

See [EVALUATION.md](EVALUATION.md) for the scientific rationale behind these
choices.
