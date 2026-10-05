# µODE — Genome-scale model reconstruction

Phase 1 of the µODE pipeline turns raw MAG FASTA files into genome-scale metabolic
models (GEMs) in SBML format. This document covers the reconstruction engines, the
namespace constraint, CheckM2 quality gating, and the LP gap-filling step.

---

## Overview

```
MAGs (*.fna)
    │
    ▼  [gene calling: Prodigal / Pyrodigal]
predicted proteins (*.faa)
    │
    ▼  [reconstruction: CarveMe or gapseq]
draft GEMs (*.xml, SBML)
    │
    ▼  [LP gap-filling: muode.gapfill]
refined GEMs (*.xml, SBML)
```

---

## Reconstruction engines

### CarveMe (default)

CarveMe performs **top-down** reconstruction: it starts from a universal template
model — a superset of all known bacterial reactions — and prunes it to fit the
genome's gene content, producing a context-specific, ready-to-use GEM in seconds
to minutes per genome.

- **Namespace: BiGG.** All metabolite and reaction identifiers follow the
  community-standard BiGG namespace. This is *required* for multi-species
  community modeling: every species must share one metabolite namespace so the
  shared extracellular pool works correctly.
- **Gene-calling prerequisite.** CarveMe expects a protein FASTA. µODE's Snakemake
  workflow calls **Prodigal** first to predict genes from the raw nucleotide
  FASTA.

CLI (single MAG or small set):

```bash
muode build --mags ./data/raw_mags/ --outdir ./models/draft_gems/ --engine carveme
```

Snakemake (scale — hundreds of MAGs in parallel):

```bash
snakemake --use-conda --cores 8 --configfile config/config.yaml
```

Config key: `engine: "carveme"`

### gapseq (alternative)

gapseq performs **bottom-up** reconstruction: it predicts metabolic-pathway
presence from a genome and builds a GEM from scratch. It is more sensitive to
incomplete genomes but much slower (hours per genome).

- **Namespace: ModelSEED** — *incompatible* with CarveMe models. Never mix CarveMe
  (BiGG) and gapseq (ModelSEED) models in a single µODE run; choose one engine per
  run.
- **When to use:** when a bottom-up approach is scientifically preferred (e.g.
  novel taxa with no close template). Requires the `gapseq` conda env.

Config key: `engine: "gapseq"`

### Offline engine (dependency-free, any architecture)

The built-in offline engine (`engine: "stub"`) is a dependency-free reconstruction
backend. It reads a tag in the MAG's first FASTA header (`muode-stub:<role>`) and
emits a compact, pre-built SBML model for that role (`producer`, `consumer`,
`both`). These are functional `LinprogOrganism` instances that run the full
pipeline (reconstruct → refine → QC → assemble → simulate) with no CarveMe,
gapseq, Prodigal or solver license — so the entire Snakemake DAG runs on any
machine, including one with no bioinformatics stack.

```bash
snakemake --cores 4 --configfile config/config.demo.yaml
```

Config key: `engine: "stub"`

### Eukaryote routes (fungi and protists)

CarveMe and gapseq are prokaryote-only, so µODE reconstructs eukaryotes through a
separate genome→proteins→GEM path, selected per-MAG by domain (see
[kingdoms.md](kingdoms.md)):

- **fungi** (`engine: "carvefungi"`) — MetaEuk gene calling, then CarveFungi, a
  fungal-specific reconstructor.
- **other eukaryotes** (`engine: "eukaryote_generic"`) — MetaEuk gene calling,
  then eggNOG-mapper orthology and a ModelSEEDpy draft. This route produces a
  draft model; prefer a curated template (e.g. Yeast8) when one is available.

Both need a MetaEuk protein reference database (`euk_ref_db`). The Snakemake
workflow routes fungus/eukaryote MAGs to these engines automatically when
`euk_ref_db` is set; a curated model supplied per-MAG in `eukaryote_models`
overrides automated reconstruction. Viruses (phages) have no GEM and are routed to
the ecology layer as `PhageInfection` models instead.

---

## Namespace constraint — why you cannot mix engines

A community simulation has one shared extracellular pool. When species A secretes
`glc__D_e` (BiGG for D-glucose) and species B tries to take up `cpd00027_e`
(ModelSEED for D-glucose), the engine sees them as *different* metabolites, no
cross-feeding occurs, and the simulation is meaningless. **Always use the same
engine for every MAG in a run.**

---

## LP gap-filling (`muode.gapfill`)

Draft GEMs frequently cannot produce biomass on a defined medium because of gaps
in the genome assembly or annotation. The `refine` step runs an LP gap-fill:

1. Identifies the reactions whose addition would restore a positive biomass flux.
2. Solves a minimal-cardinality LP to find the smallest set to add.
3. Writes the result to the refined model.

An optional **universal model** (SBML, any source) seeds the candidate reaction
set. Specify it with `--universal` (CLI) or `universal_model:` (config).

```bash
muode refine --models ./models/draft_gems/ --outdir ./models/refined_gems/
```

The `refine` command also runs fast structural QC on each model and writes a
`{stem}.qc.json` (see `muode qc`).

---

## CheckM2 quality gating (`muode.qc`, Phase 0)

Before spending compute on reconstruction, low-quality MAG bins can be filtered
out by CheckM2 completeness and contamination scores.

CheckM2 is **part of the pipeline** but **off by default** (`run_checkm2: false`)
because it requires a Diamond database and a bioinformatics stack not present in
every environment. Enable it on a capable machine:

```yaml
# config/config.yaml
run_checkm2: true
min_completeness: 50.0   # drop MAGs below 50% completeness
max_contamination: 10.0  # drop MAGs above 10% contamination
```

When CheckM2 is disabled, all MAGs pass the quality gate. The `select_mags`
checkpoint in the Snakemake workflow writes the passing list to
`{outdir}/qc/passing_mags.txt` regardless of whether CheckM2 ran.

---

## Reconstruction summary (`reconstruction_summary.tsv`)

After all per-MAG `refine` jobs finish, the `reconstruct_report` rule aggregates
every `{stem}.qc.json` into one TSV. Columns include:

| Column | Description |
|--------|-------------|
| `mag_id` | MAG stem name |
| `grows_now` | Can the model produce biomass after gap-fill? |
| `n_reactions_added` | Reactions added by LP gap-fill |
| `energy_generating_cycle` | True if an EGC artefact was detected |
| `simulatable` | Overall pass/fail used by `assemble` |

The `assemble` step reads the `simulatable` column and silently drops
non-simulatable models instead of aborting the whole run — the failure-isolation
mechanism for large-scale runs where a few bad MAGs are expected.

---

## References

- Machado et al. (2018) *Fast automated reconstruction of genome-scale metabolic
  models for microbial species and communities.* Nucleic Acids Res. 46(15):7542–7553.
- Zimmermann et al. (2021) *gapseq: informed prediction of bacterial metabolic
  pathways and reconstruction of accurate metabolic models.* Genome Biol. 22:81.
- Chklovski et al. (2023) *CheckM2: a rapid, scalable and accurate tool for
  assessing microbial genome quality using machine learning.* Nature Methods.
