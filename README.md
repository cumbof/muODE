# µODE (muODE) 🦠📉

An automated, ODE-based simulation engine for modeling large-scale microbial community dynamics, AI-predicted kinetics, and metabolic cross-feeding directly from metagenome-assembled genomes.

µODE (pronounced micro-O-D-E) bridges the gap between compositional metagenomics (knowing *who* is there) and functional, predictive systems biology (knowing *what* they are doing, and *how* they will interact over time).

µODE automates the construction of Genome-Scale Metabolic Models (GEMs) from MAGs, uses deep learning to predict missing kinetic parameters, and simulates dynamic community evolution using **dynamic Flux Balance Analysis (dFBA)** — plus an optional **ecology layer** (pH, bile acids, sporulation, antibiotic pharmacokinetics, bacteriocins, oxygen, bacteriophage predation) and **multi-kingdom** support (bacteria, archaea, fungi/eukaryotes, viruses) that turn a metabolic core into a mechanistic model of community ecology.

## 📑 Table of Contents

- [Key Features](#-key-features)
- [What Works Today](#-what-works-today)
- [The µODE Pipeline Architecture](#-the-µode-pipeline-architecture)
- [Ecology & Multi-Kingdom Layers](#-ecology--multi-kingdom-layers)
- [Installation](#️-installation)
- [Quick Start & CLI Usage](#-quick-start--cli-usage)
- [Worked Examples](#-worked-examples)
- [Inputs and Outputs](#-inputs-and-outputs)
- [Mathematical Framework](#-mathematical-framework)
- [Perturbation & Antibiotic Modeling](#-perturbation--antibiotic-modeling)
- [Documentation](#-documentation)
- [Citation](#-citation)
- [License](#-license)

## 🚀 Key Features

- **Automated GEM Reconstruction:** high-throughput extraction of context-specific metabolic models directly from known and unknown MAGs (CarveMe/gapseq for prokaryotes; a MetaEuk→CarveFungi route for fungi; a dependency-free `stub` engine for end-to-end local testing).
- **AI-Driven Gap Filling:** deep learning to flag probable metabolic modules in fragmented metagenomic assemblies, seeding algorithmic LP gap filling.
- **Kinetic Parameter Prediction:** Graph Neural Networks and protein language models to predict species-specific $k_{cat}$ and $K_m$ values without experimental data (an *optional refinement layer*), feeding Michaelis–Menten uptake bounds and **GECKO-lite enzyme constraints with a shared protein-pool budget**.
- **Massive-Scale dynamic FBA:** solves community dynamics — cross-feeding, nutrient depletion, biomass growth — over time with a native Static-Optimization-Approach integrator.
- **Ecology layer:** pH/SCFA inhibition, bile-acid transformation, sporulation/germination, time-varying antibiotic PK/PD, and bacteriocin antagonism, all composable on top of the metabolic core.
- **Multi-kingdom:** fungi/protists simulate through the *same* engine as bacteria; oxygen-tolerance gating + scavenging couple the kingdoms; bacteriophages enter as Levin–Stewart infection ODEs (not faked as metabolic species).
- **Perturbation Sandbox:** simulate environmental stressors, targeted antibiotic knockouts, and probiotic interventions, and watch **secondary extinctions** cascade through the cross-feeding network.
- **Timed interventions:** transplant a whole community mid-run (FMT) or dose a probiotic, as an emergent competition outcome rather than a hard-coded result.
- **Spatial dFBA, validation, reporting, interoperability:** a 2D reaction–diffusion colony/biofilm engine; a benchmark framework that scores a run against known expectations; a one-file HTML report; and a COMETS export bridge for cross-validation.

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
- **Timed biomass injection** (`muode.inject`): transplants (FMT), probiotic
  doses and inoculation as mid-run state events.
- **Michaelis–Menten kinetics + diets**, cross-feeding inference, time-course CSVs
  and figures.
- **Kinetics refinement** (`muode.predict`): a dependency-free heuristic predictor
  for $K_m$/$k_{cat}$ (opt-in DLKcat/Kroll wrappers behind the `ml` extra) and a
  **GECKO-lite enzyme-constraint layer** (`muode.enzyme`) that caps intracellular
  reaction velocities from $k_{cat}$, including a **shared protein-pool budget**.
- **Ecology layer** (`muode.ecology`): pH/SCFA weak-acid inhibition, bile-acid
  transformation + secondary-bile-acid toxicity, sporulation/germination life
  cycle, one-compartment antibiotic PK/PD, and bacteriocin antagonism — composable
  modifiers that are a strict no-op when unused.
- **Multi-kingdom** (`muode.traits`, `muode.oxygen`, `muode.phage`): oxygen-
  tolerance gating + facultative O₂ scavenging, and bacteriophage predation as a
  coupled infection ODE; eukaryote metabolism runs on the unchanged engine.
- **Validation & scale:** a **benchmark framework** (`muode validate`) scoring a
  run against a known community (relative-abundance error, SCFA/metabolite error,
  cross-feeding-edge F1) and **abundance-aware subsampling** for huge communities.
- **Spatial dynamic FBA** (`muode.spatial`, `muode spatial`): a 2D reaction-
  diffusion colony/biofilm engine — per-cell community FBA with metabolite
  diffusion — reproducing **spatial cross-feeding gradients**.
- **Reporting & interoperability:** `muode report` bundles a run into one
  self-contained `report.html`; `muode export-comets` writes a COMETS layout +
  `cometspy` driver to cross-check the integrator against COMETS.
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

- A **test suite** (108 tests) asserting the biology (cross-feeding-supported
  growth; secondary extinction; the *E. coli* acetate switch; stub-MAG
  cross-feeding through SBML; pH/bile/spore/antibiotic/oxygen/phage mechanisms;
  protein-pool allocation; FMT engraftment).

## 🧬 The µODE Pipeline Architecture

µODE orchestrates a phased pipeline, integrating computational-biology and AI tools into a single, reproducible workflow:

### Phase 0 — MAG Quality Control
- Gate MAGs on completeness/contamination (CheckM2) before any modelling effort is spent — low-quality bins make bad models. CheckM2 stays in the pipeline but is gated by `run_checkm2` (off where it cannot run; all MAGs then pass).

### Phase 1 — Automated Network Reconstruction
- Ingests FASTA files of MAGs (known strains and novel/unknown clusters); calls genes with Prodigal/Pyrodigal (prokaryotes) or **MetaEuk** (eukaryotes).
- **Defaults to CarveMe** (top-down, **BiGG** namespace — fast and consistent across MAGs, which is required for a shared community metabolite pool). **gapseq** is available as an alternative engine (ModelSEED namespace) but is *not mixed* with CarveMe models in one community.
- **Domain-aware routing** (`muode.traits.reconstruction_route`, driven by a `traits.tsv` side file): bacteria/archaea → CarveMe; **fungi → MetaEuk + CarveFungi**; other eukaryotes → an experimental MetaEuk + eggNOG/ModelSEEDpy draft; **viruses → no GEM** (they become `PhageInfection` ecology layers). CarveMe and gapseq are prokaryote-only — gapseq even crashes on eukaryotic contigs — so eukaryotes are routed to eukaryote-aware tools, never faked.

### Phase 2 — Gap Filling
- Algorithmic **Linear-Programming gap filling** (cobra) so every draft model can mathematically produce biomass on a defined medium.
- Optional **MetaPathPredict** to predict probable KEGG modules in incomplete MAGs (an opt-in refinement; requires KEGG↔BiGG namespace mapping).

### Phase 3 — Predicting Kinetic Parameters
- **$K_m$ + $V_{max}$ constrain substrate uptake** via Michaelis–Menten in the dynamic loop (this is the coupling between the intracellular LP and the extracellular ODEs).
- **$k_{cat}$ constrains intracellular reaction velocities** via enzyme-constrained (GECKO-lite) bounds, $|v_r| \le k_{cat,r}\cdot[E]$, plus an optional **shared protein-pool budget** $\sum_r (MW_r/(k_{cat,r}\cdot 3600))\,|v_r| \le P$ that forces the cell to *allocate* one finite enzyme mass between pathways (the mechanism behind overflow metabolism).
- A **default `heuristic` predictor** (`muode.predict`) assigns literature $K_m$ for common substrates and a $k_{cat}$ around the genome-wide median (~13.7/s; Bar-Even 2011) — **dependency-free**, so the whole pipeline runs everywhere. Values are *placeholders*, not measurements.
- **Deep-learning predictors are opt-in** (`ml` extra): DLKcat for $k_{cat}$, Kroll et al. for $K_m$. **The pipeline runs end-to-end on the heuristic without them**; predictions *sharpen* the bounds.

### Phase 4 — Dynamic Community Simulation (dFBA)
- A **native Static-Optimization-Approach integrator** (Mahadevan et al., 2002) updates extracellular metabolite concentrations and per-species biomass at discrete time steps.
- Species interact through a shared extracellular pool (COMETS-style compartmentalised community). **MICOM** is used for what it is designed for — steady-state cross-feeding / abundance-constrained community snapshots — and can serve as a per-step solver; it is *not* itself the time integrator. **COMETS** is supported as an export target for independent cross-validation (`muode export-comets`).
- Applies a defined "diet" (e.g., human gut Western diet) as the environmental boundary condition.
- **Pluggable, batch-aware LP backends.** Each time step hands *all* of its FBA LPs (every species, and in
  `DynamicFBA.run_ensemble` every ensemble member) to a `SolverBackend` at once. The default backend reproduces the
  historical per-organism solves byte for byte. The [manylp](https://github.com/cumbof/manylp) backend
  (`make_backend("manylp-gpu")` / `"manylp-cpu"`) certifies whole batches against cached optimal bases, returns
  certified-unique fluxes (`flux_rule="pfba-unique"`), can persist its bases across runs (`atlas_dir=`), and can
  integrate alternative-optima envelopes (`policies=[...]`):

  ```python
  from muode import DynamicFBA
  from muode.backends import make_backend

  engine = DynamicFBA(t_end=48, dt=0.1)
  result = engine.run(community, diet, kinetics, backend=make_backend("manylp-gpu"))
  runs = engine.run_ensemble(community, diets, kinetics,
                             backend=make_backend("manylp-gpu", atlas_dir="~/.muode/atlas"))
  ```

### Phase 5 — Perturbation, Intervention & Ecology
- Dynamic knockouts (e.g. a broad-spectrum antibiotic constraining a target pathway's flux toward 0 at a given *efficacy*); **secondary extinctions** from disrupted cross-feeding emerge naturally.
- Timed **injections** (FMT, probiotic dose) introduce biomass mid-run.
- The **ecology layer** adds non-metabolic processes (pH, bile, spores, antibiotic PK, bacteriocins, oxygen, phages) — see below.

## 🌍 Ecology & Multi-Kingdom Layers

The dFBA core models **metabolism**. Real communities are also shaped by chemistry, pharmacology, antagonism, dormancy and predation. muODE adds these as composable, optional **ecology layers** (`muode.ecology.EcologyModel`) — an empty model is a strict no-op, numerically identical to the bare engine.

| Layer | Models | Module |
|-------|--------|--------|
| `WeakAcidInhibition` | pH dynamics + undissociated-SCFA toxicity (self-limiting fermentation; colonization resistance) | `muode.ph` |
| `BileAcidTransform` + `BileAcidInhibition` | BSH + 7α-dehydroxylation of the bile pool; secondary-bile-acid growth inhibition | `muode.bile` |
| `SporeForming` | vegetative↔spore life cycle; spores survive antibiotics, germinate on a bile germinant | `muode.lifecycle` |
| `Antibiotic` | one-compartment PK + Emax kill on vegetative cells (basis of recurrence) | `muode.antibiotic` |
| `Bacteriocin` | diffusible-toxin interference competition | `muode.antagonism` |
| `OxygenSensitivity` | O₂-tolerance growth gating + facultative O₂ scavenging (keeps the lumen anaerobic) | `muode.oxygen` |
| `PhageInfection` | Levin–Stewart lytic/temperate predation, coupled to host biomass | `muode.phage` |

**Multi-kingdom.** muODE's core is *paradigm-defined, not taxon-defined*: any GEM that satisfies the `OrganismModel` interface — bacterial, archaeal, fungal or protist — simulates with the same code. What differs between kingdoms is (a) reconstruction tooling (Phase 1 routing) and (b) the oxygen relationship (the `OxygenSensitivity` layer). Phages have no metabolism, so they are infection ODEs, not metabolic species. See **[docs/ecology.md](docs/ecology.md)** and **[docs/kingdoms.md](docs/kingdoms.md)**.

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

External per-phase tools (CarveMe, CheckM2, gapseq, CarveFungi) are provisioned
automatically by Snakemake from `workflow/envs/` when you run with `--use-conda`.

## 💻 Quick Start & CLI Usage

µODE exposes a CLI for each pipeline phase. Commands marked ✅ run today;
⚙️ call external reconstruction tools (driven at scale by the Snakemake
workflow). Each step below notes what it **consumes** and what it **produces**
so the chain is unambiguous.

**0. Try the engine immediately ✅**

No data needed — the built-in two-species cross-feeding demo runs with only
numpy/scipy:

```bash
muode demo --outdir results/demo
```

---

**1. Build Draft GEMs from MAGs ⚙️**

*Consumes:* MAG FASTA files (`.fna`/`.fasta`) you provide in `data/raw_mags/`.
*Produces:* one SBML model per MAG in `models/draft_gems/`.

```bash
muode build --mags ./data/raw_mags/ --outdir ./models/draft_gems/ --engine carveme
```

> For scale (dozens–thousands of MAGs), use the Snakemake workflow: it fans
> reconstruction out per MAG with optional CheckM2 QC gating and domain-aware
> routing (bacteria→CarveMe, fungi→CarveFungi, viruses→skipped).

---

**2. Gap-Fill & Predict Kinetics ⚙️**

*Consumes:* draft GEMs from step 1.
*Produces:* gap-filled GEMs + one `{stem}.kinetics.json` per model in `models/kinetic_gems/`.

```bash
# heuristic predictor (no extra deps) writes a {stem}.kinetics.json per model
muode refine --models ./models/draft_gems/ --outdir ./models/kinetic_gems/ --predict-kinetics
# deep-learning kinetics (needs the `ml` extra):
# muode refine ... --predictor dlkcat
```

---

**3. Assemble the Community ✅**

*Consumes:* gap-filled GEMs from step 2 + a user-supplied abundance TSV.
*Produces:* `simulation_env/community.json` — consumed by all downstream commands.

**Abundance TSV** is a simple 2-column file (`mag_id`, `rel_abundance`) that
you produce with whatever profiler you ran on the sample (MetaSBT, Bracken,
Kraken2, mOTUs, MetaPhlAn, or a manual estimate). muODE does not care how the
sample was profiled — the TSV is the interface. If you omit it, equal abundance
is assumed for every MAG.

```bash
muode assemble \
  --models    ./models/kinetic_gems/ \
  --abundance ./data/abundance.tsv \
  --diet      western_gut \
  --outdir    ./simulation_env/
```

For very large communities, subsample by abundance (pick one):

```bash
muode assemble ... --coverage 0.95       # fewest top MAGs covering 95% of abundance
muode assemble ... --max-species 100     # at most 100 MAGs
muode assemble ... --min-abundance 0.001 # drop anything below 0.1%
```

---

**4. Run the Dynamic Simulation ✅**

*Consumes:* `community.json` from step 3 (+ optionally the kinetics files from step 2).
*Produces:* `results/biomass.csv`, `metabolites.csv`, `growth_rates.csv`,
`cross_feeding.csv` and figures.

```bash
# basic simulation:
muode simulate --community ./simulation_env/community.json \
               --time 24 --step 0.1 --outdir ./results/

# with predicted kinetics + GECKO-lite enzyme constraints and a protein-pool budget:
muode simulate --community ./simulation_env/community.json \
               --kinetics ./models/kinetic_gems/ --enzyme-constraints --protein-pool 0.2 \
               --time 24 --step 0.1 --outdir ./results/

# transplant a donor community mid-run (FMT) at t = 24 h:
muode simulate --community ./recipient_env/community.json \
               --inject ./donor_env/community.json --inject-time 24 \
               --time 96 --step 0.1 --outdir ./results_fmt/
```

---

**5. Validate Against a Known Community ✅**

*Consumes:* results directory from step 4 + a benchmark expectation YAML.
A ready-made YAML for the built-in toy community is in `examples/benchmarks/`.

```bash
muode validate --results ./results/ \
               --expected examples/benchmarks/toy_cross_feeding.yaml
```

---

**6. Spatial (colony / biofilm) Simulation ✅**

Standalone — no previous steps required. Uses the same built-in cross-feeding
demo community on a 2D reaction-diffusion grid.

```bash
muode spatial --nx 24 --time 12 --outdir ./results/spatial/
```

---

**7. Report & COMETS export ✅**

```bash
# bundle a run into one self-contained, shareable HTML file:
muode report --results ./results/                 # writes ./results/report.html

# export the community to COMETS for an independent dynamic-FBA cross-check:
muode export-comets --community ./simulation_env/community.json --outdir ./comets_run/
```

---

**Run the whole pipeline with Snakemake (recommended for many MAGs):**

```bash
conda activate muode
snakemake --use-conda --cores 8 --configfile config/config.yaml
```

**Validate the full pipeline locally first (no external tools, any architecture):**

```bash
# dependency-free `stub` engine on two bundled toy MAGs — no --use-conda needed
snakemake --cores 4 --configfile config/config.demo.yaml
```

## 🔬 Worked Examples

Each example ships a `README.md`, an abundance TSV, a config, a `download_genomes.sh`,
and (where relevant) a dependency-light `mechanistic_demo.py` that runs **anywhere**.

| Example | What it shows |
|---------|---------------|
| [`examples/gut_western/`](examples/gut_western/) | A dense Western-diet gut community: carbon competition + SCFA cross-feeding. |
| [`examples/fmt_cdiff/`](examples/fmt_cdiff/) | **Fecal microbiota transplant** for recurrent *C. difficile*: control (recurrence) vs treatment (engraftment) differ only by a timed donor injection; the full ecology stack (bile acids, spores, antibiotic PK, pH) reproduces the recurrence-vs-cure contrast. |
| [`examples/multikingdom/`](examples/multikingdom/) | **Bacteria + fungus + phage:** dropping the fungus collapses the obligate anaerobe (oxygen), dropping the phage unleashes the pathobiont (predation) — both emergent. |
| [`examples/strain_competition/`](examples/strain_competition/) | **One species, three strains:** on glucose the best grower excludes the rest; a colicin producer overturns it (interference); a private-niche strain coexists — the winner depends on which competition dominates. |
| [`examples/benchmarks/`](examples/benchmarks/) | Benchmark-expectation YAMLs for `muode validate`. |

```bash
# run any mechanistic demo on any machine (toy models, no GEMs needed):
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/multikingdom/mechanistic_demo.py
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/strain_competition/mechanistic_demo.py
```

## 📊 Inputs and Outputs

**Inputs**
- **Genomes:** `.fasta`/`.fna` MAG files you provide (one per bin).
- **Abundance profile:** a 2-column TSV (`mag_id`, `rel_abundance`) you provide.
  Any upstream profiler works; muODE is agnostic to how the sample was profiled.
  Omit it to use equal abundance for all MAGs.
- **Traits (optional):** a side file (`mag_id`, `domain`, `oxygen`) carrying
  kingdom and oxygen relationship; keeps the abundance TSV a plain 2-column file.
- **Diet / media profile:** a preset name (e.g. `western_gut`) or a CSV of
  metabolite concentrations (mmol/L) with an optional influx column.

**Outputs**
- **Metabolic Models:** standard `.xml` (SBML) / `.json` per species.
- **Time-course Data:** `biomass.csv`, `metabolites.csv`, `growth_rates.csv`;
  with the ecology layer, also `environment.csv` (pH, drug, germination signal)
  and `spores.csv`.
- **Interactions:** `cross_feeding.csv` (producer → metabolite → consumer).
- **Visualization & report:** stacked-area composition charts, metabolite
  trajectories, cross-feeding network graphs, and a one-file `report.html`.

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

Here $\mu_i$ is the growth rate of species $i$ (the LP objective value), $v_{j,i}$ the exchange flux of metabolite $j$ by species $i$, $\delta$ an optional death rate, $D$ a chemostat dilution rate, and $\phi_j$ a nutrient influx. At each step the substrate uptake bound is set by Michaelis–Menten kinetics, $v^{up}_{j,i} = V_{max}\,M_j/(K_m + M_j)$, capped so a step cannot over-deplete a metabolite. Ecology layers compose by multiplying growth/uptake factors and summing death/metabolite rates (see [docs/ecology.md](docs/ecology.md)).

## 💊 Perturbation & Antibiotic Modeling

```bash
muode perturb --community ./simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" \
              --efficacy 0.95 \
              --outdir ./results_antibiotic/
```

This restricts the flux capacity of the targeted pathway by 95%, letting you watch the real-time cascade of primary deaths and secondary cross-feeding starvation. You can also `--knockout` explicit reactions or `--remove-species` entirely. For *time-varying* drug exposure with spore survival and recurrence, use the `Antibiotic` ecology layer (see [docs/ecology.md](docs/ecology.md) and the FMT example).

## 📚 Documentation

Detailed per-feature docs live in **[docs/](docs/)**:

| Document | Covers |
|----------|--------|
| [engine.md](docs/engine.md) | SOA integrator, solver abstraction, cross-feeding |
| [reconstruction.md](docs/reconstruction.md) | CarveMe / gapseq / stub, gene calling, namespaces, gap-fill, CheckM2 |
| [kinetics.md](docs/kinetics.md) | Km/Vmax uptake, kcat enzyme constraints, protein-pool budget, predictors |
| [assembly.md](docs/assembly.md) | Community manifest, abundance TSV, diet, subsampling |
| [perturbation.md](docs/perturbation.md) | Antibiotic / knockout / species-removal engine |
| [injection.md](docs/injection.md) | Timed biomass injection (FMT, probiotics) |
| [ecology.md](docs/ecology.md) | pH/SCFA, bile acids, spores, antibiotic PK, bacteriocins, oxygen, phage |
| [kingdoms.md](docs/kingdoms.md) | Fungi/protists, phages, archaea — what is and isn't modelled |
| [validation.md](docs/validation.md) | Benchmark framework, metrics, pass/fail report |
| [spatial.md](docs/spatial.md) | 2D reaction-diffusion colony/biofilm engine |
| [workflow.md](docs/workflow.md) | Snakemake pipeline, rules, config reference, SLURM |
| [visualization.md](docs/visualization.md) | Static figures + the self-contained HTML report |
| [LIMITATIONS.md](docs/LIMITATIONS.md) | Complete, honest inventory of what muODE can and cannot model |
| [EVALUATION.md](docs/EVALUATION.md) | Scientific assessment, design rationale, references |

## 📖 Citation

A manuscript describing µODE is **in preparation**. If you use µODE in your
research in the meantime, please cite this repository:

```bibtex
@software{muode,
  author = {Cumbo, Fabio},
  title  = {{µODE}: an ODE-based simulation engine for microbial community dynamics},
  url    = {https://github.com/cumbof/muODE},
  year   = {2026}
}
```

## 📄 License

µODE is released under the **MIT License**. See [LICENSE](LICENSE) for the full text.

© 2026 Fabio Cumbo.
