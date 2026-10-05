# µODE — Community assembly

Assembly bundles refined GEMs, their relative abundances, a diet, and a total
biomass into a **community manifest** (`community.json`) that every downstream
command (`simulate`, `perturb`, `spatial`) consumes.

---

## What assembly does

1. Loads refined GEMs (SBML) from a directory.
2. Reads the reconstruction QC summary and **silently drops non-simulatable
   models** (the failure-isolation mechanism for large runs).
3. Optionally reads a user-supplied abundance TSV and applies **abundance-aware
   subsampling** to reduce very large communities.
4. Writes `community.json` with the resolved model paths, per-species abundances,
   diet choice, and total biomass.

---

## The abundance TSV

µODE accepts a simple **2-column TSV** (`mag_id`, `rel_abundance`):

```
bin1    0.45
bin2    0.30
bin3    0.15
bin4    0.10
```

Column order matters; the header line is optional. Values are interpreted as
relative abundances and are normalized automatically if they do not sum to 1.

**µODE is agnostic to how the sample was profiled.** Use MetaSBT, Bracken,
Kraken2, mOTUs, MetaPhlAn, manual estimates, or any combination — as long as you
produce this 2-column file. There is no built-in profiler; the TSV is the contract
boundary. If no abundance file is supplied, µODE assumes **equal relative
abundance** for every model that passed QC.

---

## Diet

A diet specifies the initial extracellular metabolite concentrations and influx
rates — the environmental boundary condition for the simulation.

### Built-in presets

| Preset | Description |
|--------|-------------|
| `western_gut` | A typical Western human gut diet (carbohydrates, proteins, fats) |
| `minimal_glucose` | Minimal medium, glucose as the sole carbon source |

```bash
muode assemble ... --diet western_gut
```

### Custom diet CSV

Supply a CSV with columns `metabolite`, `concentration` (mmol/L), and optionally
`influx` (mmol/L/h):

```csv
metabolite,concentration,influx
glc__D_e,10.0,0.5
o2_e,0.21,0.0
nh4_e,5.0,0.1
```

```bash
muode assemble ... --diet ./my_diet.csv
```

---

## Abundance-aware subsampling (`muode.subsample`)

For very large communities (hundreds of MAGs), simulating every species is
expensive. The `assemble` step can subsample by abundance before writing the
manifest:

```bash
# keep the top 100 most abundant MAGs
muode assemble ... --max-species 100

# drop any MAG below 0.1% relative abundance
muode assemble ... --min-abundance 0.001

# keep the fewest top MAGs that together account for 95% of total abundance
muode assemble ... --coverage 0.95
```

These options can be combined. Subsampling runs after QC filtering and before the
manifest is written; dropped MAGs are reported on stdout.

---

## The community manifest (`community.json`)

```json
{
  "models": [
    "/absolute/path/to/bin1.xml",
    "/absolute/path/to/bin2.xml"
  ],
  "abundances": {
    "bin1": 0.45,
    "bin2": 0.30
  },
  "diet": "western_gut",
  "total_biomass": 0.01
}
```

Paths are stored absolute so the manifest can be used from any working directory.
`total_biomass` (gDW/L) is the sum of all species' initial biomasses; individual
biomasses are proportional to relative abundance.

---

## CLI reference

```bash
muode assemble \
  --models    ./models/refined_gems/ \
  --abundance ./data/abundance.tsv \   # optional; equal weights if omitted
  --diet      western_gut \
  --outdir    ./simulation_env/ \
  --total-biomass 0.01 \
  --max-species 100 \                  # optional subsampling
  --min-abundance 0.001 \
  --coverage 0.95
```

Output: `simulation_env/community.json`

---

## Snakemake rule

The `assemble` rule in `workflow/Snakefile` is driven by `assemble.py` in
`workflow/scripts/`. It depends on:

- all refined GEMs that survived the `select_mags` checkpoint,
- `reconstruction_summary.tsv` for QC-driven failure isolation,
- the abundance TSV from `config["abundance"]` (if set).

Config keys relevant to assembly:

```yaml
abundance: null          # path to a 2-column abundance TSV, or null (equal weights)
diet: "western_gut"
total_biomass: 0.01
max_species: null
min_abundance: null
abundance_coverage: null
```
