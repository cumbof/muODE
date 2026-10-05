# µODE — Snakemake workflow

The `workflow/` directory contains a Snakemake pipeline that runs µODE at scale:
fanning reconstruction and refinement across hundreds of MAGs in parallel (on a
cluster) and chaining them into the integrative simulation phases.

---

## Workflow vs CLI

| Situation | Use |
|-----------|-----|
| Dozens to thousands of MAGs | Workflow (per-MAG fan-out, SLURM, checkpoint) |
| A handful of MAGs or a local run | CLI (`muode build` / `muode refine` / …) |
| A dependency-free run on any machine | `snakemake --configfile config/config.demo.yaml` |

The CLI and the workflow are complementary: both call the same `muode` library
functions, and every Snakemake rule maps to one or more CLI commands.

---

## Running the workflow

```bash
conda activate muode

# full pipeline on a capable machine (CarveMe + CheckM2):
snakemake --use-conda --cores 8 --configfile config/config.yaml

# full pipeline with SLURM (edit profiles/slurm/config.yaml first):
snakemake --use-conda --profile profiles/slurm --configfile config/config.yaml

# dependency-free run — no CarveMe, no CheckM2, no bioinformatics stack:
snakemake --cores 4 --configfile config/config.demo.yaml
```

---

## Pipeline DAG

```
MAGs (*.fna)
    │
    ├─▶ [checkm2] ──▶ [select_mags checkpoint]
    │                         │
    │                         ▼ (one job per MAG)
    │              [reconstruct] ──▶ [refine]
    │                                    │
    │                     ┌──────────────┤
    │                     │              │
    │              [reconstruct_report]  │  (optional)
    │                     │         [memote_all]
    │                     ▼
    └──────────▶ [assemble] ──▶ [simulate] ──▶ [validate]
                                    │
                                    └──▶ [perturb]
```

The `select_mags` **checkpoint** determines which per-MAG jobs (reconstruct, refine)
are created at runtime, so only MAGs that survive QC are ever modeled.

---

## Rules

### Phase 0: `checkm2` + `select_mags`

- `checkm2`: runs CheckM2 on the MAG directory. Gated by `run_checkm2: true`.
- `select_mags` *(checkpoint)*: reads the CheckM2 report (or passes all MAGs if
  CheckM2 is disabled) and writes `{outdir}/qc/passing_mags.txt`.

### Phase 1: `reconstruct`

Per-MAG rule. Calls Prodigal + CarveMe (or the offline engine). Produces one
`{outdir}/models/draft/{mag}.xml` per MAG. The engine is chosen by `engine:` in the
config.

### Phase 2/3: `refine`

Per-MAG rule. Runs `muode refine` on each draft model: LP gap-fill, structural QC,
optional kinetics prediction. Produces:

- `{outdir}/models/refined/{mag}.xml` — the refined model
- `{outdir}/models/refined/{mag}.qc.json` — gap-fill + QC report
- `{outdir}/models/refined/{mag}.kinetics.json` — predicted kinetics (if `predict_kinetics: true`)

### Phase 3b: `reconstruct_report`

Aggregates all `*.qc.json` files into `{outdir}/qc/reconstruction_summary.tsv` after
all `refine` jobs complete.

### Phase 3c: `memote` + `memote_all` (optional)

Runs memote snapshot reports on all refined models if `run_memote: true`. Produces
per-model HTML reports and a summary file.

### Phase 4a: `assemble`

Calls `workflow/scripts/assemble.py`. Reads the reconstruction summary, optionally
an abundance TSV, applies subsampling, and writes `{outdir}/simulation/community.json`.

### Phase 4b: `simulate`

Calls `muode simulate`. Produces `{outdir}/simulation/biomass.csv` plus the full
output set. Depends on `kinetics` jobs if `predict_kinetics: true`.

### Phase 5: `perturb` (optional)

Calls `muode perturb`. Enabled by a `perturbation:` block in the config. Produces
`{outdir}/perturbation/biomass.csv`.

### Phase 6: `validate` (optional)

Calls `muode.validate`. Enabled by a `benchmark:` entry in the config. Produces
`{outdir}/validation/report.json` and fails the rule if benchmarks are not met.

---

## Config reference (`config/config.yaml`)

```yaml
# --- inputs / outputs -------------------------------------------------------
mags_dir: "data/raw_mags"     # directory of MAG FASTA files
mag_extension: "fna"          # file extension (without dot)
abundance: null               # 2-column TSV (mag_id, rel_abundance) or null
outdir: "results"

# --- Phase 0: MAG quality gate ----------------------------------------------
run_checkm2: false            # enable on a capable machine
min_completeness: 50.0        # % completeness threshold
max_contamination: 10.0       # % contamination threshold

# --- Phase 1: reconstruction ------------------------------------------------
engine: "carveme"             # carveme | gapseq | stub
universe: "bacteria"          # CarveMe template
gapfill_media: null           # gap-fill medium (e.g. "M9")

# --- Phase 2/3: refinement --------------------------------------------------
universal_model: null         # universal SBML for LP gap-filling
run_memote: false             # generate per-model memote HTML reports
predict_kinetics: false       # predict Km + kcat during refine
predictor: "heuristic"        # heuristic | dlkcat | km-ml
enzyme_constraints: false     # GECKO-lite kcat caps at simulate time
default_vmax: 10.0
default_km: 0.01

# --- Phase 4: simulation ----------------------------------------------------
diet: "western_gut"
total_biomass: 0.01
t_end: 24.0
dt: 0.1

# --- scale: abundance-aware subsampling ------------------------------------
max_species: null
min_abundance: null
abundance_coverage: null

# --- Phase 6: validation (optional) ----------------------------------------
benchmark: null               # path to benchmark expectation YAML

# --- compute ----------------------------------------------------------------
threads: 4

# --- Phase 5: perturbation (optional) ---------------------------------------
perturbation: null
# perturbation:
#   target_pathway: "Folate Biosynthesis"
#   efficacy: 0.95
#   remove_species: "Bin.42"
#   knockout: "FOLD,DHFR"
```

---

## Demo config (`config/config.demo.yaml`)

The demo config uses the offline `stub` engine and the bundled MAGs
(`data/demo_mags/`) so the entire DAG runs without CarveMe, CheckM2 or any conda env:

```bash
snakemake --cores 4 --configfile config/config.demo.yaml
```

Use it to confirm the workflow is wired correctly on a new machine before running
real data.

---

## SLURM profile

`profiles/slurm/config.yaml` configures Snakemake to submit each rule as a separate
SLURM job. Edit the partition, account and resource defaults before use:

```bash
snakemake --profile profiles/slurm --configfile config/config.yaml
```

Per-rule resource overrides (memory, CPUs) can be added to the `config.yaml` under
`cluster_config:` or via the `resources:` key within each rule.

---

## Conda environments

Each heavy external tool runs in its own conda environment (provisioned by
`--use-conda`):

| File | Tools |
|------|-------|
| `workflow/envs/carveme.yaml` | Prodigal, CarveMe, memote |
| `workflow/envs/checkm2.yaml` | CheckM2 + Diamond database |
| `workflow/envs/muode.yaml` | µODE + cobra + scipy |

The µODE-internal steps (refine, assemble, simulate, perturb, validate) run in the
active `muode` conda environment where you ran `pip install -e .`.
