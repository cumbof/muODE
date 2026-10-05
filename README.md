# µODE

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Solver](https://img.shields.io/badge/solver-HiGHS%20(free)-brightgreen.svg)](#installation)

**An ODE-based simulation engine for microbial community dynamics, built directly from metagenome-assembled genomes.**

µODE (pronounced *micro-O-D-E*) takes you from compositional metagenomics — knowing *who* is in a sample — to functional, predictive systems biology — knowing *what* the community does and *how* it changes over time. It reconstructs genome-scale metabolic models from MAGs, couples them through a shared extracellular pool, and integrates the whole community forward in time with **dynamic Flux Balance Analysis (dFBA)**.

On top of the metabolic core, µODE layers the chemistry and ecology that actually shape a gut community — pH and short-chain-fatty-acid inhibition, bile-acid transformation, sporulation, antibiotic pharmacokinetics, bacteriocins, oxygen gradients, and bacteriophage predation — and runs bacteria, archaea, fungi and viruses through one engine.

## Contents

- [Why µODE](#why-µode)
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [The pipeline](#the-pipeline)
- [Command-line interface](#command-line-interface)
- [Ecology and multi-kingdom modeling](#ecology-and-multi-kingdom-modeling)
- [Inputs and outputs](#inputs-and-outputs)
- [Mathematical framework](#mathematical-framework)
- [Worked examples](#worked-examples)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

## Why µODE

Metagenomics tells you the roster of a community. It does not tell you how that roster behaves — who feeds whom, which species a perturbation will knock out, how an intervention propagates through the cross-feeding network. µODE closes that gap with a mechanistic, dynamic model assembled automatically from the genomes themselves:

- **From genomes, not assumptions.** Metabolic capabilities come from the MAGs. Interactions are *emergent* — cross-feeding, competitive exclusion and secondary extinctions arise from a shared metabolite pool, never from hard-coded rules.
- **Dynamic, not static.** A native dynamic-FBA integrator produces full time courses of biomass, metabolites and growth rates, not a single steady-state snapshot.
- **Metabolism *and* ecology.** A composable ecology layer adds the non-metabolic processes that govern real communities, so µODE can model colonization resistance, recurrence after antibiotics, and phage therapy — not just carbon flux.
- **Runs anywhere.** The core engine uses the free HiGHS solver and has no commercial-license requirement; a dependency-light backend runs the entire pipeline offline on any architecture.

## Features

- **Automated GEM reconstruction** from MAGs: CarveMe (default, BiGG namespace) and gapseq (ModelSEED) for prokaryotes, a MetaEuk→CarveFungi route for fungi, and a dependency-free offline engine for portable runs.
- **Kinetic parameterization**: Michaelis–Menten uptake bounds ($K_m$, $V_{max}$) couple the intracellular linear program to the extracellular ODEs; GECKO-lite enzyme constraints cap intracellular velocities from $k_{cat}$, with an optional shared protein-pool budget that drives overflow metabolism. A built-in predictor assigns parameters with no external dependencies; deep-learning predictors (DLKcat, Kroll $K_m$) are available through the `ml` extra.
- **Community dynamic FBA**: a native Static-Optimization-Approach integrator solves community dynamics — cross-feeding, nutrient depletion, biomass growth — over a shared extracellular pool, with parallel per-species solves.
- **Ecology layer**: pH / SCFA inhibition, bile-acid transformation and secondary-bile-acid toxicity, sporulation and germination, time-varying antibiotic PK/PD, and bacteriocin antagonism — composable modifiers on top of the metabolic core.
- **Multi-kingdom**: bacteria, archaea, fungi and protists run through the same engine; oxygen-tolerance gating and facultative O₂ scavenging couple the kingdoms; bacteriophages enter as Levin–Stewart infection ODEs.
- **Perturbations and interventions**: antibiotic pathway attenuation, reaction knockouts, species removal, and timed biomass injection (fecal microbiota transplant, probiotic dosing) — all producing emergent community responses.
- **Spatial dynamic FBA**: a 2D reaction–diffusion colony/biofilm engine that resolves spatial cross-feeding gradients.
- **Validation and interoperability**: a benchmark framework that scores a run against known expectations (relative-abundance error, metabolite error, cross-feeding-edge F1), a self-contained HTML report, and a COMETS export bridge for independent cross-validation.
- **Scale**: a Snakemake workflow that fans reconstruction out per MAG (quality-control gating and a SLURM profile included) and abundance-aware subsampling for large communities.

## Installation

µODE's core engine runs on the free, open-source **HiGHS** solver — no commercial license required. Gurobi or CPLEX are drop-in alternatives for very large communities.

**Prerequisites:** Python 3.10+, Conda or Mamba. *(Optional: a Gurobi/CPLEX academic license for large-scale FBA.)*

```bash
git clone https://github.com/cumbof/muODE.git
cd muODE

# Create and activate the environment (core engine, CLI, workflow driver)
mamba env create -f environment.yml
conda activate muode

# Install the package
pip install -e .

# Verify
muode demo
```

The external per-phase tools (CarveMe, CheckM2, gapseq, CarveFungi) are provisioned automatically by Snakemake from `workflow/envs/` when you run with `--use-conda`.

## Quick start

The built-in two-species cross-feeding simulation runs immediately, with only numpy and scipy:

```bash
muode demo --outdir results/demo
```

`A_glucose` ferments glucose and secretes acetate; `B_acetate` grows only on that acetate. B blooms *after* A has produced its substrate — the canonical cross-feeding motif. Remove the feeder and watch the dependent species go secondarily extinct:

```bash
muode demo --outdir results/demo_ko --remove A_glucose
```

## The pipeline

µODE orchestrates a phased, reproducible pipeline from raw MAGs to a dynamic community simulation.

**Phase 0 — MAG quality control.** CheckM2 gates MAGs on completeness and contamination before any modeling effort is spent. Controlled by `run_checkm2`; when disabled, all MAGs pass.

**Phase 1 — Network reconstruction.** Genes are called with Prodigal/Pyrodigal (prokaryotes) or MetaEuk (eukaryotes). CarveMe is the default reconstructor (BiGG namespace, consistent across MAGs, which a shared community pool requires); gapseq is available as an alternative (ModelSEED namespace), never mixed with CarveMe in one community. Domain-aware routing (`muode.traits.reconstruction_route`) sends bacteria and archaea to CarveMe, fungi to MetaEuk + CarveFungi, other eukaryotes to a MetaEuk + eggNOG/ModelSEEDpy draft, and viruses to the ecology layer as `PhageInfection` models rather than GEMs.

**Phase 2 — Gap filling.** A minimal-cardinality LP gap-fill (cobra) ensures every model can produce biomass on a defined medium. MetaPathPredict can optionally flag probable KEGG modules in incomplete MAGs.

**Phase 3 — Kinetics.** $K_m$ and $V_{max}$ constrain substrate uptake through Michaelis–Menten bounds — the coupling between the intracellular LP and the extracellular ODEs. $k_{cat}$ constrains intracellular velocities through enzyme-constrained (GECKO-lite) bounds, with an optional shared protein-pool budget that forces the cell to allocate a finite enzyme mass between pathways. The default predictor is dependency-free; DLKcat and Kroll $K_m$ predictors are available through the `ml` extra.

**Phase 4 — Community simulation.** A native Static-Optimization-Approach integrator (Mahadevan et al., 2002) updates extracellular concentrations and per-species biomass at each time step, with species interacting through one shared extracellular pool under a defined diet. A MICOM/cobra community model can serve as the per-step solver; COMETS is supported as an export target for independent cross-validation.

**Phase 5 — Perturbation, intervention and ecology.** Antibiotic pathway attenuation, reaction knockouts, species removal and timed biomass injections (FMT, probiotics) drive the run; the ecology layer adds pH, bile, spores, antibiotic PK, bacteriocins, oxygen and phages. Secondary extinctions from disrupted cross-feeding emerge naturally.

## Command-line interface

µODE exposes one subcommand per pipeline phase. Each step notes what it **consumes** and **produces** so the chain is unambiguous.

**Try the engine immediately** — no data needed:

```bash
muode demo --outdir results/demo
```

**1. Build draft GEMs from MAGs.** *Consumes* MAG FASTA files; *produces* one SBML model per MAG.

```bash
muode build --mags ./data/raw_mags/ --outdir ./models/draft_gems/ --engine carveme
```

For dozens to thousands of MAGs, use the Snakemake workflow: it fans reconstruction out per MAG with CheckM2 gating and domain-aware routing.

**2. Gap-fill and predict kinetics.** *Consumes* draft GEMs; *produces* gap-filled GEMs and one `{stem}.kinetics.json` per model.

```bash
muode refine --models ./models/draft_gems/ --outdir ./models/kinetic_gems/ --predict-kinetics
# deep-learning kinetics (needs the `ml` extra):
# muode refine ... --predictor dlkcat
```

**3. Assemble the community.** *Consumes* gap-filled GEMs and an abundance TSV; *produces* `community.json`, consumed by all downstream commands.

```bash
muode assemble \
  --models    ./models/kinetic_gems/ \
  --abundance ./data/abundance.tsv \
  --diet      western_gut \
  --outdir    ./simulation_env/
```

The abundance TSV is a simple 2-column file (`mag_id`, `rel_abundance`) produced by whatever profiler you ran — MetaSBT, Bracken, Kraken2, mOTUs, MetaPhlAn, or a manual estimate. µODE is agnostic to how the sample was profiled; the TSV is the interface. Omit it for equal abundance. For large communities, subsample by abundance:

```bash
muode assemble ... --coverage 0.95       # fewest top MAGs covering 95% of abundance
muode assemble ... --max-species 100     # at most 100 MAGs
muode assemble ... --min-abundance 0.001 # drop anything below 0.1%
```

**4. Run the dynamic simulation.** *Consumes* `community.json`; *produces* `biomass.csv`, `metabolites.csv`, `growth_rates.csv`, `cross_feeding.csv` and figures.

```bash
# basic simulation:
muode simulate --community ./simulation_env/community.json \
               --time 24 --step 0.1 --outdir ./results/

# with predicted kinetics, GECKO-lite enzyme constraints and a protein-pool budget:
muode simulate --community ./simulation_env/community.json \
               --kinetics ./models/kinetic_gems/ --enzyme-constraints --protein-pool 0.2 \
               --time 24 --step 0.1 --outdir ./results/

# transplant a donor community mid-run (FMT) at t = 24 h:
muode simulate --community ./recipient_env/community.json \
               --inject ./donor_env/community.json --inject-time 24 \
               --time 96 --step 0.1 --outdir ./results_fmt/
```

**5. Validate against a known community.** *Consumes* a results directory and a benchmark YAML.

```bash
muode validate --results ./results/ \
               --expected examples/benchmarks/toy_cross_feeding.yaml
```

**6. Spatial (colony / biofilm) simulation.** Standalone; runs the built-in cross-feeding community on a 2D reaction–diffusion grid.

```bash
muode spatial --nx 24 --time 12 --outdir ./results/spatial/
```

**7. Report and COMETS export.**

```bash
muode report --results ./results/                 # one self-contained report.html
muode export-comets --community ./simulation_env/community.json --outdir ./comets_run/
```

**Run the whole pipeline with Snakemake** (recommended for many MAGs):

```bash
conda activate muode
snakemake --use-conda --cores 8 --configfile config/config.yaml
```

A dependency-free configuration runs the entire DAG (reconstruct → refine → QC → assemble → simulate) offline, on any architecture, with no external tools:

```bash
snakemake --cores 4 --configfile config/config.demo.yaml
```

## Ecology and multi-kingdom modeling

The dFBA core models metabolism. Real communities are also shaped by chemistry, pharmacology, antagonism, dormancy and predation, which µODE adds as composable **ecology layers** (`muode.ecology.EcologyModel`). An empty model is a strict no-op — numerically identical to the bare engine.

| Layer | Models | Module |
|-------|--------|--------|
| `WeakAcidInhibition` | pH dynamics and undissociated-SCFA toxicity (self-limiting fermentation; colonization resistance) | `muode.ph` |
| `BileAcidTransform` + `BileAcidInhibition` | BSH and 7α-dehydroxylation of the bile pool; secondary-bile-acid growth inhibition | `muode.bile` |
| `SporeForming` | vegetative↔spore life cycle; spores survive antibiotics and germinate on a bile germinant | `muode.lifecycle` |
| `Antibiotic` | one-compartment PK and Emax kill on vegetative cells | `muode.antibiotic` |
| `Bacteriocin` | diffusible-toxin interference competition | `muode.antagonism` |
| `OxygenSensitivity` | O₂-tolerance growth gating and facultative O₂ scavenging | `muode.oxygen` |
| `PhageInfection` | Levin–Stewart lytic/temperate predation, coupled to host biomass | `muode.phage` |

**Multi-kingdom.** µODE's core is paradigm-defined, not taxon-defined: any GEM satisfying the `OrganismModel` interface — bacterial, archaeal, fungal or protist — simulates with the same code. What differs between kingdoms is the reconstruction tooling (Phase 1 routing) and the oxygen relationship (the `OxygenSensitivity` layer). Phages have no metabolism, so they are infection ODEs rather than metabolic species. See [docs/ecology.md](docs/ecology.md) and [docs/kingdoms.md](docs/kingdoms.md).

## Inputs and outputs

**Inputs**

- **Genomes:** `.fasta`/`.fna` MAG files, one per bin.
- **Abundance profile:** a 2-column TSV (`mag_id`, `rel_abundance`). Any upstream profiler works; omit it for equal abundance.
- **Traits (optional):** a side file (`mag_id`, `domain`, `oxygen`) carrying kingdom and oxygen relationship, keeping the abundance TSV a plain 2-column file.
- **Diet / medium:** a preset name (e.g. `western_gut`) or a CSV of metabolite concentrations (mmol/L) with an optional influx column.

**Outputs**

- **Metabolic models:** standard SBML (`.xml`) / JSON per species.
- **Time-course data:** `biomass.csv`, `metabolites.csv`, `growth_rates.csv`; with the ecology layer, also `environment.csv` (pH, drug, germination signal) and `spores.csv`.
- **Interactions:** `cross_feeding.csv` (producer → metabolite → consumer).
- **Visualization and report:** stacked-area composition charts, metabolite trajectories, cross-feeding network graphs, and a one-file `report.html`.

## Mathematical framework

µODE implements **dynamic Flux Balance Analysis** through the Static Optimization Approach. Inside each cell, a pseudo-steady-state linear program maximizes the biomass objective $Z$:

$$ \text{Maximize } Z = c^T v $$
$$ \text{Subject to } S \cdot v = 0 $$
$$ v_{min} \leq v \leq v_{max} $$

where $S$ is the stoichiometric matrix and $v$ the flux vector. Outside the cell, the environment evolves by ODEs. For biomass $X_i$ and metabolite $M_j$:

$$ \frac{dX_i}{dt} = (\mu_i - \delta - D)\, X_i $$
$$ \frac{dM_j}{dt} = \sum_{i} v_{j,i}\, X_i + \phi_j - D\, M_j $$

Here $\mu_i$ is the growth rate of species $i$ (the LP objective value), $v_{j,i}$ the exchange flux of metabolite $j$ by species $i$, $\delta$ an optional death rate, $D$ a chemostat dilution rate, and $\phi_j$ a nutrient influx. At each step the substrate uptake bound is set by Michaelis–Menten kinetics, $v^{up}_{j,i} = V_{max}\,M_j/(K_m + M_j)$, capped so a step cannot over-deplete a metabolite. Ecology layers compose by multiplying growth/uptake factors and summing death/metabolite rates. See [docs/engine.md](docs/engine.md) and [docs/ecology.md](docs/ecology.md).

## Worked examples

Each example ships a `README.md`, an abundance TSV, a config, a `download_genomes.sh`, and (where relevant) a dependency-light `mechanistic_demo.py` that runs anywhere.

| Example | What it shows |
|---------|---------------|
| [`examples/gut_western/`](examples/gut_western/) | A dense Western-diet gut community: carbon competition and SCFA cross-feeding. |
| [`examples/fmt_cdiff/`](examples/fmt_cdiff/) | Fecal microbiota transplant for recurrent *C. difficile*: control and treatment differ only by a timed donor injection, and the full ecology stack (bile acids, spores, antibiotic PK, pH) reproduces the recurrence-vs-cure contrast. |
| [`examples/multikingdom/`](examples/multikingdom/) | Bacteria + fungus + phage: dropping the fungus collapses the obligate anaerobe (oxygen), dropping the phage unleashes the pathobiont (predation) — both emergent. |
| [`examples/strain_competition/`](examples/strain_competition/) | One species, three strains: on glucose the best grower excludes the rest; a colicin producer overturns it (interference); a private-niche strain coexists. |
| [`examples/phage_therapy/`](examples/phage_therapy/) | A lytic phage crashes a *Klebsiella* bloom and releases the resources a commensal needs (predation and resource-release emerge together). |
| [`examples/benchmarks/`](examples/benchmarks/) | Benchmark-expectation YAMLs for `muode validate`. |

```bash
# run a mechanistic demo on any machine (no GEMs needed):
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/multikingdom/mechanistic_demo.py
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/strain_competition/mechanistic_demo.py
```

## Documentation

Per-feature documentation lives in [docs/](docs/):

| Document | Covers |
|----------|--------|
| [engine.md](docs/engine.md) | SOA integrator, solver abstraction, cross-feeding |
| [reconstruction.md](docs/reconstruction.md) | CarveMe / gapseq / offline / eukaryote engines, gene calling, namespaces, gap-fill, CheckM2 |
| [kinetics.md](docs/kinetics.md) | Km/Vmax uptake, kcat enzyme constraints, protein-pool budget, predictors |
| [assembly.md](docs/assembly.md) | Community manifest, abundance TSV, diet, subsampling |
| [perturbation.md](docs/perturbation.md) | Antibiotic / knockout / species-removal engine |
| [injection.md](docs/injection.md) | Timed biomass injection (FMT, probiotics) |
| [ecology.md](docs/ecology.md) | pH/SCFA, bile acids, spores, antibiotic PK, bacteriocins, oxygen, phage |
| [kingdoms.md](docs/kingdoms.md) | Fungi/protists, phages, archaea |
| [validation.md](docs/validation.md) | Benchmark framework, metrics, pass/fail report |
| [spatial.md](docs/spatial.md) | 2D reaction–diffusion colony/biofilm engine |
| [workflow.md](docs/workflow.md) | Snakemake pipeline, rules, config reference, SLURM |
| [visualization.md](docs/visualization.md) | Static figures and the self-contained HTML report |

## Citation

A manuscript describing µODE is in preparation. If you use µODE in your research in the meantime, please cite this repository:

```bibtex
@software{muode,
  author = {Cumbo, Fabio},
  title  = {{µODE}: an ODE-based simulation engine for microbial community dynamics},
  url    = {https://github.com/cumbof/muODE},
  year   = {2026}
}
```

## License

µODE is released under the MIT License. See [LICENSE](LICENSE) for the full text.

© 2026 Fabio Cumbo.
