# Example: 12-species western gut community

This example builds and simulates a realistic human gut community representative
of a western dietary pattern (high simple carbohydrates, low dietary fibre).
All 12 genomes are public NCBI RefSeq reference assemblies, and the abundance
table was constructed from the literature to reflect the ecological fingerprint
of a western gut.

Run from the repository root unless noted otherwise. Requires CarveMe and
Prodigal (see [docs/reconstruction.md](../../docs/reconstruction.md));
on aarch64 or machines without CarveMe, substitute `engine: stub` — the
pipeline will run but the biology will not be meaningful.

---

## Community composition and ecological rationale

The 12 species cover the dominant phyla and the main metabolic guilds of the
human gut. The relative abundances reflect the published western gut pattern:
Bacteroides dominance (Firmicutes/Bacteroidetes ratio tilted toward Bacteroidetes),
reduced butyrate producers, low Akkermansia and Prevotella.

| MAG ID | Organism | Phylum | Functional role | Rel. abundance |
|--------|----------|--------|-----------------|---------------|
| `B_thetaiotaomicron_VPI5482` | *Bacteroides thetaiotaomicron* VPI-5482 | Bacteroidota | Primary carbohydrate degrader (PULs); secretes **acetate + succinate** | 0.22 |
| `B_fragilis_NCTC9343` | *Bacteroides fragilis* NCTC 9343 | Bacteroidota | Capsular polysaccharide; immunomodulatory; secretes **propionate + acetate** | 0.13 |
| `E_rectale_ATCC33656` | *Eubacterium rectale* ATCC 33656 | Bacillota | Major butyrate producer (acetyl-CoA pathway); **cross-feeds on acetate** | 0.12 |
| `Bl_obeum_A2162` | *Blautia obeum* A2-162 | Bacillota | Hydrogen-consuming **acetogen**; enables thermodynamically favourable fermentation by consuming H₂ | 0.08 |
| `R_intestinalis_L182` | *Roseburia intestinalis* L1-82 | Bacillota | Butyrate from starch/arabinoxylan; **cross-feeds on Bacteroides succinate** | 0.09 |
| `F_prausnitzii_A2165` | *Faecalibacterium prausnitzii* A2-165 | Bacillota | Anti-inflammatory butyrate producer; **marker species for gut health**; low in western diet | 0.08 |
| `R_bromii_L263` | *Ruminococcus bromii* L2-63 | Bacillota | **Keystone starch degrader**; releases glucose/maltose consumed by other species | 0.07 |
| `B_longum_NCC2705` | *Bifidobacterium longum* NCC2705 | Actinomycetota | Ferments oligosaccharides (GOS, FOS); secretes **acetate + lactate** cross-fed to butyrate producers | 0.07 |
| `P_copri_DSM18205` | *Prevotella copri* DSM 18205 | Bacteroidota | Hemicellulose and plant polysaccharide degrader; low in western diet | 0.04 |
| `A_muciniphila_BAA835` | *Akkermansia muciniphila* ATCC BAA-835 | Verrucomicrobiota | Mucin layer degrader; secretes **propionate + acetate**; abundance reduced by high-fat diet | 0.04 |
| `L_acidophilus_NCFM` | *Lactobacillus acidophilus* NCFM | Bacillota | Produces **lactate** from simple sugars; cross-feeds *Coprococcus* | 0.03 |
| `C_comes_ATCC27758` | *Coprococcus comes* ATCC 27758 | Bacillota | Converts **lactate → butyrate** (cross-feeds from *Bifidobacterium* and *Lactobacillus*) | 0.03 |

### Cross-feeding network (expected)

```
glucose/fructose (diet)
    │
    ├─▶ B. thetaiotaomicron ──acetate──▶ E. rectale ──butyrate──▶ host
    │         └──────────────succinate──▶ R. intestinalis ──butyrate──▶ host
    │
    ├─▶ R. bromii ──glucose/maltose──▶ E. rectale, R. intestinalis
    │
    ├─▶ B. longum ──acetate──▶ E. rectale ──butyrate──▶ host
    │       └────────lactate──▶ C. comes ──butyrate──▶ host
    │
    ├─▶ L. acidophilus ──lactate──▶ C. comes
    │
    └─▶ Bl. obeum (consumes H₂, enables other fermenters)

A. muciniphila: mucin ──▶ propionate + acetate
```

The diet for this example (`western_gut` preset) supplies glucose (10 mmol/L,
1 mmol/L/h influx), fructose (5 mmol/L, 0.5 mmol/L/h), and lactose (2 mmol/L)
as primary nutrients. For a production analysis, replace with a full VMH western
diet CSV (see [docs/assembly.md](../../docs/assembly.md)):
```yaml
diet: "path/to/vmh_western_diet.csv"
```

---

## Prerequisites

```bash
# muODE itself
conda activate muode

# CarveMe + Prodigal (for GEM reconstruction)
conda activate carveme    # or whichever env has CarveMe

# NCBI Datasets CLI (for genome download)
conda install -c conda-forge ncbi-datasets-cli
```

---

## Step 0: Download reference genomes

```bash
bash examples/gut_western/download_genomes.sh
```

This creates one `.fna` file per species in
`examples/gut_western/data/mags/`. The script downloads each genome from
NCBI RefSeq using `ncbi-datasets` and renames the FASTA to the MAG ID used
in `abundance.tsv`. Check the table at the top of `download_genomes.sh` if
any accession needs updating.

---

## Step 1: Reconstruct genome-scale metabolic models

Build one draft GEM per genome using CarveMe. This step requires
Prodigal (gene calling) and CarveMe (template-based reconstruction).

```bash
muode build \
  --mags examples/gut_western/data/mags \
  --engine carveme \
  --outdir examples/gut_western/models/draft
```

Expected output: 12 SBML files in `examples/gut_western/models/draft/`,
one per genome (e.g. `B_thetaiotaomicron_VPI5482.xml`). CarveMe may
print warnings for reactions not found in the template; this is normal.

---

## Step 2: Refine models (LP gap-fill + structural QC)

```bash
muode refine \
  --models examples/gut_western/models/draft \
  --outdir examples/gut_western/models/refined
```

Each model is gap-filled to ensure a biomass reaction is active. A QC
report (`*.qc.json`) and refined SBML are written for every model.
Models that still cannot grow after gap-filling are flagged in the report
and will be dropped during assembly (see `simulatable` column in
`reconstruction_summary.tsv`).

To also predict Km and kcat (heuristic predictor, no extra dependencies):

```bash
muode refine \
  --models examples/gut_western/models/draft \
  --outdir examples/gut_western/models/refined \
  --predict-kinetics \
  --predictor heuristic
```

---

## Step 3: Assemble the community

Combine the refined models with the abundance table and dietary conditions
into a community manifest (`community.json`).

```bash
muode assemble \
  --models examples/gut_western/models/refined \
  --abundance examples/gut_western/abundance.tsv \
  --diet western_gut \
  --total-biomass 0.01 \
  --outdir examples/gut_western/simulation_env
```

This writes `examples/gut_western/simulation_env/community.json`.
Each model is paired with its relative abundance from `abundance.tsv`; the
abundances are used to initialise species biomass proportionally from
`--total-biomass` (0.01 gDW/L total community).

The `--min-abundance 0.01` flag (in `config.yaml`) would drop anything below
1% relative abundance; at the CLI, add it with `--min-abundance 0.01`.

---

## Step 4: Simulate

Run the dynamic community ODE simulation for 48 hours of simulated time.

```bash
muode simulate \
  --community examples/gut_western/simulation_env/community.json \
  --time 48.0 \
  --step 0.1 \
  --outdir examples/gut_western/results
```

Outputs written to `examples/gut_western/results/`:

| File | Content |
|------|---------|
| `biomass.csv` | Biomass time-course per species (gDW/L × time) |
| `metabolites.csv` | Extracellular metabolite concentrations (mmol/L × time) |
| `cross_feeding.csv` | Detected cross-feeding interactions (mean flux product) |
| `meta.json` | Run metadata (solver, dt, t_end, diet, …) |
| `biomass.png` | Stacked-area chart of community composition over time |
| `metabolites.png` | Line chart of the 15 most dynamic metabolites |
| `cross_feeding.png` | Directed cross-feeding network (requires networkx) |

### What to look for

**`biomass.png`:**  
*B. thetaiotaomicron* and *B. fragilis* should dominate at steady state,
consistent with Bacteroides dominance in a western dietary pattern. Butyrate
producers (*E. rectale*, *R. intestinalis*) should bloom after a lag phase
once Bacteroides have seeded the acetate pool.

**`metabolites.png`:**  
Glucose and fructose decrease rapidly. Acetate rises first (primary
fermentation product), then partially depletes as butyrate producers consume
it. Butyrate and propionate accumulate as end products. The acetate
peak-and-decline is the signature cross-feeding signal: it appears, gets
cross-fed, then stabilises at a low pool.

**`cross_feeding.png`:**  
Expect directed edges from *B. thetaiotaomicron* → acetate → butyrate
producers, and from *B. longum* / *L. acidophilus* → lactate →
*C. comes*.

---

## Step 5: Validate against expected outcomes

```bash
muode validate \
  --results examples/gut_western/results \
  --expected examples/gut_western/benchmark.yaml
```

The benchmark checks that Bacteroides remain dominant (MAE ≤ 0.15),
that acetate and glucose show the expected dynamics (metabolite RMSE ≤
5 mmol/L), and that at least one of the expected cross-feeding edges is
detected (F1 ≥ 0.25). Tolerances are generous because exact outcomes
depend on GEM quality.

---

## Step 6: Perturbation — keystone species removal

*Ruminococcus bromii* is a keystone starch degrader: it breaks down
resistant starch into oligosaccharides that other species depend on.
Simulating its removal shows how losing a single ecological specialist
cascades through the cross-feeding network.

```bash
muode perturb \
  --community examples/gut_western/simulation_env/community.json \
  --remove-species R_bromii_L263 \
  --time 48.0 \
  --outdir examples/gut_western/results_perturb_Rbromii
```

Compare `results/biomass.png` with `results_perturb_Rbromii/biomass.png`.
Expected observation: species that directly consumed *R. bromii*-derived
products (e.g. *E. rectale*, *R. intestinalis*) show reduced final biomass
because the starch-derived glucose/maltose supply is gone. This is a
**secondary extinction** cascade — those species are not targeted by the
perturbation but starve indirectly.

A second meaningful perturbation: remove *Faecalibacterium prausnitzii*,
a marker species for gut dysbiosis and IBD:

```bash
muode perturb \
  --community examples/gut_western/simulation_env/community.json \
  --remove-species F_prausnitzii_A2165 \
  --time 48.0 \
  --outdir examples/gut_western/results_perturb_Fpr
```

Total community butyrate production should decrease in the perturbed run.

---

## Alternative: run the full pipeline with Snakemake

The single-config Snakemake run reconstructs, refines, assembles, simulates
and validates in one command (parallelises across all 12 MAGs):

```bash
conda activate muode
snakemake --use-conda --cores 8 --configfile examples/gut_western/config.yaml
```

The perturbation arm is disabled in the default `config.yaml`. To enable it,
uncomment the `perturbation:` block at the bottom of the config:

```yaml
perturbation:
  remove_species: "R_bromii_L263"
```

Results land in `examples/gut_western/results/`.

---

## Exploring results in a notebook

```python
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = "examples/gut_western/results"

bio = pd.read_csv(f"{OUTDIR}/biomass.csv", index_col=0)
met = pd.read_csv(f"{OUTDIR}/metabolites.csv", index_col=0)

# relative abundance over time (normalised)
rel = bio.div(bio.sum(axis=1), axis=0)
ax = rel.plot.area(figsize=(12, 5), title="Community composition over time")
ax.set_ylabel("relative abundance")
ax.set_xlabel("time (h)")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/relative_abundance.png", dpi=150)

# SCFA dynamics
scfa = met[["ac_e", "but_e", "ppa_e"]].rename(columns={
    "ac_e": "acetate", "but_e": "butyrate", "ppa_e": "propionate"
})
scfa.plot(figsize=(10, 4), title="Short-chain fatty acids over time")
plt.ylabel("concentration (mmol/L)")
plt.xlabel("time (h)")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/scfa.png", dpi=150)
```
