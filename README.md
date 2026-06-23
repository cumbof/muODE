# µODE (muODE) 🦠📉

An automated, ODE-based simulation engine for modeling large-scale microbial community dynamics, AI-predicted kinetics, and metabolic cross-feeding directly from metagenome-assembled genomes.

µODE (pronounced micro-O-D-E) bridges the gap between compositional metagenomics (knowing *who* is there) and functional, predictive systems biology (knowing *what* they are doing, and *how* they will interact over time).

Designed to seamlessly ingest species classifications and Metagenome-Assembled Genomes (MAGs) from upstream binning/taxonomy pipelines (e.g., [MetaSBT](https://github.com/cumbof/MetaSBT)), µODE automates the construction of Genome-Scale Metabolic Models (GEMs), uses deep learning to predict missing kinetic parameters, and simulates dynamic community evolution using **dynamic Flux Balance Analysis (dFBA)**.

> **Status — v0.1 (core engine working).** The dynamic-FBA community engine, the
> perturbation engine, the diet/kinetics layer, the CLI and the Snakemake workflow
> scaffold are implemented and tested (toy cross-feeding community **and** a real
> *E. coli* core dFBA that reproduces the textbook acetate-overflow result). The
> reconstruction/AI phases are wired as tool wrappers + workflow rules and are the
> focus of the next milestone. See **[docs/EVALUATION.md](docs/EVALUATION.md)** for
> the scientific assessment and **[docs/DESIGN.md](docs/DESIGN.md)** for the
> architecture.

## 📑 Table of Contents

- [Key Features](#-key-features)
- [What Works Today](#-what-works-today)
- [The µODE Pipeline Architecture](#-the-µode-pipeline-architecture)
- [Installation](#️-installation)
- [Quick Start & CLI Usage](#-quick-start--cli-usage)
- [Inputs and Outputs](#-inputs-and-outputs)
- [Mathematical Framework](#-mathematical-framework)
- [Perturbation & Antibiotic Modeling](#-perturbation--antibiotic-modeling)
- [Roadmap](#️-roadmap)
- [Acknowledgments](#-acknowledgments)

## 🚀 Key Features

- **Automated GEM Reconstruction:** high-throughput extraction of context-specific metabolic models directly from known and unknown MAGs.
- **AI-Driven Gap Filling:** deep learning to flag probable metabolic modules in fragmented metagenomic assemblies, seeding algorithmic LP gap filling.
- **Kinetic Parameter Prediction:** Graph Neural Networks and protein language models to predict species-specific $k_{cat}$ and $K_m$ values without experimental data (an *optional refinement layer*).
- **Massive-Scale dynamic FBA:** solves community dynamics — cross-feeding, nutrient depletion, biomass growth — over time with a native Static-Optimization-Approach integrator.
- **Perturbation Sandbox:** simulate environmental stressors, targeted antibiotic knockouts, and probiotic interventions, and watch **secondary extinctions** cascade through the cross-feeding network.

## ✅ What Works Today

```bash
muode demo                       # two-species cross-feeding sim, no data needed
muode demo --remove A_glucose    # watch the dependent species go secondarily extinct
```

- A **native dynamic-FBA engine** (`muode.dfba`) implementing the Static
  Optimization Approach over a community sharing one extracellular pool.
- Two solver backends behind one interface: **cobra** (real GEMs) and a
  dependency-light **scipy/HiGHS** backend (tests, toy demo, no licence).
- **Perturbation engine**: reaction knockouts, pathway attenuation by *efficacy*,
  and whole-species removal.
- **Michaelis–Menten kinetics + diets**, cross-feeding inference, time-course CSVs
  and figures.
- **Kinetics refinement** (`muode.predict`): a dependency-free heuristic predictor
  for $K_m$/$k_{cat}$ (opt-in DLKcat/Kroll wrappers behind the `ml` extra) and a
  **GECKO-lite enzyme-constraint layer** (`muode.enzyme`) that caps intracellular
  reaction velocities from $k_{cat}$.
- **Validation & scale:** a **benchmark framework** (`muode validate`) scoring a
  run against a known community (relative-abundance error, SCFA/metabolite error,
  cross-feeding-edge F1), **abundance-aware subsampling** for huge communities,
  and a quantitative-profile reader (`muode.metasbt`).
- **Spatial dynamic FBA** (`muode.spatial`, `muode spatial`): a 2D reaction-
  diffusion colony/biofilm engine — per-cell community FBA with metabolite
  diffusion — reproducing **spatial cross-feeding gradients**.
- A **Snakemake workflow** that fans reconstruction/refinement out per MAG (QC
  checkpoint + SLURM profile included), aggregates a per-run
  `reconstruction_summary.tsv`, and **drops non-simulatable models** instead of
  aborting the whole run.
- A dependency-free **`stub` reconstruction engine** so the *entire* pipeline
  (reconstruct → refine → QC → assemble → simulate) can be run and validated on
  any machine — including ones where CarveMe/CheckM2 cannot run — before scaling
  out on a cluster:

  ```bash
  snakemake --cores 4 --configfile config/config.demo.yaml   # full DAG on toy MAGs
  ```

- A **test suite** asserting the biology (cross-feeding-supported growth; secondary
  extinction; the *E. coli* acetate switch; stub-MAG cross-feeding through SBML).

## 🧬 The µODE Pipeline Architecture

µODE orchestrates a phased pipeline, integrating computational-biology and AI tools into a single, reproducible workflow:

### Phase 0 — MAG Quality Control *(added)*
- Gate MAGs on completeness/contamination (CheckM2) before any modelling effort is spent — low-quality bins make bad models.

### Phase 1 — Automated Network Reconstruction
- Ingests FASTA files of MAGs (known strains and novel/unknown clusters); calls genes with Prodigal/Pyrodigal.
- **Defaults to CarveMe** (top-down, **BiGG** namespace — fast and consistent across MAGs, which is required for a shared community metabolite pool). **gapseq** is available as an alternative engine (ModelSEED namespace) but is *not mixed* with CarveMe models in one community.

### Phase 2 — Gap Filling
- Algorithmic **Linear-Programming gap filling** (cobra) so every draft model can mathematically produce biomass on a defined medium.
- Optional **MetaPathPredict** to predict probable KEGG modules in incomplete MAGs (an opt-in refinement; requires KEGG↔BiGG namespace mapping).

### Phase 3 — Predicting Kinetic Parameters *(implemented)*
- **$K_m$ + $V_{max}$ constrain substrate uptake** via Michaelis–Menten in the dynamic loop (this is the coupling between the intracellular LP and the extracellular ODEs).
- **$k_{cat}$ constrains intracellular reaction velocities** via enzyme-constrained (GECKO-lite) bounds, $|v_r| \le k_{cat,r}\cdot[E]$ — a *separate* layer applied once before the run (`muode.enzyme`).
- A **default `heuristic` predictor** (`muode.predict`) assigns literature $K_m$ for common substrates and a $k_{cat}$ around the genome-wide median (~13.7/s; Bar-Even 2011) — **dependency-free**, so the whole pipeline runs everywhere. Values are *placeholders*, not measurements.
- **Deep-learning predictors are opt-in** (`ml` extra): DLKcat for $k_{cat}$, Kroll et al. for $K_m$. They need a per-reaction *enzyme/substrate context* (sequence + SMILES), and `build_enzyme_context` provides the BiGG→(sequence, SMILES) mapping. **The pipeline runs end-to-end on the heuristic without them**; predictions *sharpen* the bounds.

### Phase 4 — Dynamic Community Simulation (dFBA)
- A **native Static-Optimization-Approach integrator** (Mahadevan et al., 2002) updates extracellular metabolite concentrations and per-species biomass at discrete time steps.
- Species interact through a shared extracellular pool (COMETS-style compartmentalised community). **MICOM** is used for what it is designed for — steady-state cross-feeding / abundance-constrained community snapshots — and can serve as a per-step solver; it is *not* itself the time integrator. **COMETS** is a planned alternative dynamic backend.
- Applies a defined "diet" (e.g., human gut Western diet) as the environmental boundary condition.

### Phase 5 — Perturbation Engine
- Dynamic knockouts (e.g. a broad-spectrum antibiotic constraining a target pathway's flux toward 0 at a given *efficacy*).
- **Secondary extinctions** caused by disrupted cross-feeding emerge naturally from the simulation.

## ⚙️ Installation

µODE's **core engine runs on a free, open-source solver (HiGHS)** — no commercial licence required. Gurobi or CPLEX are recommended (and drop-in) only for very large communities.

**Prerequisites:** Python 3.10+, Conda/Mamba. *(Optional: Gurobi/CPLEX academic licence for large-scale FBA.)*

```bash
# Clone
git clone https://github.com/cumbof/muODE.git
cd muODE

# Create and activate the environment (core engine + CLI + workflow driver)
mamba env create -f environment.yml
conda activate muode

# Install the muODE package
pip install -e .

# Verify
muode demo
pytest -q
```

External per-phase tools (CarveMe, CheckM2, gapseq) are provisioned automatically
by Snakemake from `workflow/envs/` when you run with `--use-conda`.

## 💻 Quick Start & CLI Usage

µODE exposes an intuitive CLI for each pipeline phase. Commands marked ✅ run
today; ⚙️ shell out to external reconstruction tools (driven at scale by the
Snakemake workflow).

**0. Try the engine immediately ✅**

```bash
muode demo --outdir results/demo
```

**1. Build Models from MAGs ⚙️**

```bash
muode build --mags ./data/raw_mags/ --outdir ./models/draft_gems/ --engine carveme
```

**2. Gap-Fill & (optionally) Predict Kinetics ⚙️**

```bash
# heuristic predictor (no extra deps) writes a {stem}.kinetics.json per model
muode refine --models ./models/draft_gems/ --outdir ./models/kinetic_gems/ --predict-kinetics
# (use --predictor dlkcat|km-ml with the `ml` extra for the deep-learning models)
```

Then feed the predicted kinetics into the simulation, optionally enabling the
enzyme-constraint ($k_{cat}$) layer:

```bash
muode simulate --community ./simulation_env/community.json \
               --kinetics ./models/kinetic_gems/ --enzyme-constraints --outdir ./results/
```

**3. Assemble Community & Diet ✅**

```bash
muode assemble --models ./models/kinetic_gems/ --abundance ./data/profile.tsv \
               --diet western_gut --outdir ./simulation_env/
# huge community? subsample by abundance, and read a MetaSBT profile directly:
muode assemble --models ./models/kinetic_gems/ --abundance ./data/metasbt_profile.tsv \
               --metasbt --coverage 0.95 --outdir ./simulation_env/
```

**4. Run the Dynamic Simulation ✅**

```bash
muode simulate --community ./simulation_env/community.json --time 24 --step 0.1 --outdir ./results/
```

**5. Validate Against a Known Community ✅**

```bash
muode validate --results ./results/ --expected examples/benchmarks/toy_cross_feeding.yaml
```

**6. Spatial (colony / biofilm) Simulation ✅**

```bash
# 2D reaction-diffusion dynamic FBA; no data needed (built-in cross-feeding demo)
muode spatial --nx 24 --time 12 --outdir ./results/spatial/
```

**Or run the whole pipeline with Snakemake (recommended for many MAGs):**

```bash
conda activate muode
snakemake --use-conda --cores 8 --configfile config/config.yaml
```

**Validate the whole pipeline locally first (no external tools, any architecture):**

```bash
# uses the dependency-free `stub` engine on two bundled toy MAGs; no --use-conda
snakemake --cores 4 --configfile config/config.demo.yaml
```

## 📊 Inputs and Outputs

**Inputs**
- **Genomes:** `.fasta`/`.fna` files containing your MAGs/strains.
- **Abundance Profile:** a 2-column TSV mapping MAG ids to relative abundance (e.g. a MetaSBT profile).
- **Diet/Media Profile:** a preset name, or a CSV defining nutrient concentrations (mmol/L) and optional influx.

**Outputs**
- **Metabolic Models:** standard `.xml` (SBML) / `.json` per species.
- **Time-course Data:** `biomass.csv`, `metabolites.csv`, `growth_rates.csv` tracking every species and extracellular metabolite at each step.
- **Interactions:** `cross_feeding.csv` (producer → metabolite → consumer).
- **Visualization:** stacked-area composition charts, metabolite trajectories, and cross-feeding network graphs.

## 🧮 Mathematical Framework

µODE relies on **dynamic Flux Balance Analysis (dFBA)** via the Static Optimization Approach.
Inside each cell, a pseudo-steady-state linear program maximizes the biomass objective $Z$:

$$ \text{Maximize } Z = c^T v $$
$$ \text{Subject to } S \cdot v = 0 $$
$$ v_{min} \leq v \leq v_{max} $$

where $S$ is the stoichiometric matrix and $v$ the flux vector.
Outside the cell, the environment evolves by ODEs. For biomass $X_i$ and metabolite $M_j$:

$$ \frac{dX_i}{dt} = (\mu_i - \delta - D)\, X_i $$
$$ \frac{dM_j}{dt} = \sum_{i} v_{j,i}\, X_i + \phi_j - D\, M_j $$

Here $\mu_i$ is the growth rate of species $i$ (the LP objective value), $v_{j,i}$ the exchange flux of metabolite $j$ by species $i$, $\delta$ an optional death rate, $D$ a chemostat dilution rate, and $\phi_j$ a nutrient influx. At each step the substrate uptake bound is set by Michaelis–Menten kinetics, $v^{up}_{j,i} = V_{max}\,M_j/(K_m + M_j)$, capped so a step cannot over-deplete a metabolite.

## 💊 Perturbation & Antibiotic Modeling

```bash
muode perturb --community ./simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" \
              --efficacy 0.95 \
              --outdir ./results_antibiotic/
```

This restricts the flux capacity of the targeted pathway by 95%, letting you watch the real-time cascade of primary deaths and secondary cross-feeding starvation. You can also `--knockout` explicit reactions or `--remove-species` entirely.

## 🗺️ Roadmap

See **[docs/EVALUATION.md](docs/EVALUATION.md)** for the full, justified plan.

- [x] **M0 — Core engine & scaffold:** native dFBA integrator, perturbation engine, kinetics/diet layer, cross-feeding inference, CLI, Snakemake workflow, tests.
- [~] **M1 — Reconstruction at scale (in progress):** ✅ end-to-end DAG fan-out per MAG with a `stub` engine that runs the *whole* pipeline locally; ✅ `reconstruction_summary.tsv` aggregation; ✅ QC-driven failure isolation (non-simulatable models dropped, not fatal); ✅ optional memote rule + universal-model gap-fill wiring. **Remaining:** validate CarveMe/Prodigal + CheckM2 on real MAG sets on a capable host; multi-threading tuning for large clusters (>500 MAGs).
- [~] **M2 — Kinetics refinement (in progress):** ✅ kinetic-parameter store with $k_{cat}$ + persistence; ✅ dependency-free heuristic predictor (default); ✅ GECKO-lite enzyme-constraint layer wired through CLI + workflow + dFBA; ✅ opt-in DLKcat/Kroll wrappers + BiGG→(sequence, SMILES) context. **Remaining:** bundle/validate real DLKcat & Kroll checkpoints on a GPU host; full protein-*pool* GECKO budget; ESM-2 embeddings.
- [~] **M3 — Validation & scale (in progress):** ✅ benchmark/validation framework (`muode validate`: relative-abundance MAE + Spearman, SCFA/metabolite error, cross-feeding edge F1, pass/fail vs. tolerances); ✅ abundance-aware subsampling (`top_n`/`min_abundance`/`coverage`) for huge communities; ✅ MetaSBT profile ingestion contract (`muode.metasbt`). **Remaining:** benchmark against real synthetic/gut datasets on a capable host; COMETS alternative dynamic backend.
- [~] **M4 — Reach (in progress):** ✅ spatiotemporal (PDE) 2D reaction-diffusion colony/biofilm engine (`muode.spatial`, `muode spatial`) — per-cell community FBA + metabolite diffusion, reproducing spatial cross-feeding gradients; ✅ spatial figures. **Remaining:** interactive GUI dashboard (deferred — not headless-testable); performance work for large GEMs on large grids; ESM-2 embeddings.

## 🙏 Acknowledgments

µODE stands on the shoulders of giants in the open-source systems-biology community, including the dynamic-FBA formulation of **Mahadevan et al. (2002)**, **COBRApy**, **CarveMe**, **gapseq**, **MICOM**, **COMETS**, **MetaPathPredict**, **DLKcat**, **CheckM2**, and **memote**. See [docs/EVALUATION.md](docs/EVALUATION.md) for citations.
