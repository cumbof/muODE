# Example: phage therapy — a lytic phage reins in a pathobiont bloom

A fast-growing pathobiont out-competes a commensal for a shared carbon source; a **lytic
phage specific to the pathobiont** crashes its bloom and *releases* the commensal. The
predation, the phage amplification, and the commensal's competitive rebound are all
**emergent** — the only difference between the two arms is whether the phage is present.

The reusable mechanism — `muode.phage.PhageInfection`, a Levin–Stewart predator–prey
layer — lives in the package. A phage is **not** an FBA organism (no metabolism, no
biomass objective), so it is modelled where it belongs: adsorption removes host biomass,
a latent infected pool lyses into a burst of new virions, coupled to the metabolic world
only through host biomass.

## Two levels, same mechanism

### Toy (runs on any machine) — `phage_scenario.py`
Two dependency-light `LinprogOrganism` growers on one sugar; the phage targets the fast
one. Exercises the *mechanism* without a genome-scale model. Pinned by
`tests/test_phage.py::test_phage_therapy_scenario_controls_bloom` (host crashes to ~0,
commensal released, and the phage titer **amplifies** above its dose — so the crash is
predation, not coincidence).

```bash
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/phage_therapy/phage_scenario.py
```

### Genome-scale (real GEMs) — `genome_scale.py`
The same contrast on **real reconstructions**, so who out-blooms whom is decided by the
genomes' stoichiometry, not a dial:

| member | organism | role |
|--------|----------|------|
| `K_pneumoniae` | *Klebsiella pneumoniae* HS11286 | facultative pathobiont — the phage's host |
| `E_coli` | *E. coli* K-12 MG1655 | commensal competitor — **not** a host for this phage |
| *(phage)* | a lytic Klebsiella phage | `PhageInfection` layer, not a GEM |

Both are facultative gram-negatives → CarveMe (`-u gramneg`), one BiGG pool, competing on
`gut_glucose.csv` (a defined anaerobic glucose medium — verified to grow a real BiGG GEM).

```bash
bash examples/phage_therapy/download_genomes.sh          # host + commensal genomes
snakemake --use-conda --cores 8 --configfile examples/phage_therapy/config.yaml
python examples/phage_therapy/genome_scale.py --models examples/phage_therapy/results/refined
```

> **Status: runs end to end on real GEMs; both mechanisms emerge.** Full genome-scale run
> (2026-08-08, ts-02, 48 h; `genome_scale_result.json`) — final biomass, gDW/L:
>
> ```
>             K_pneumoniae   E_coli
> no phage          0.391    0.452
> + phage           0.000    0.800
> ```
>
> The phage **crashes the host completely** (K. pneumoniae 0.01 → 0 by t≈8 h) and
> **releases the commensal** (E. coli 0.45 → 0.80, +77%). It is genuine predation, not
> coincidence: the phage titer **amplifies above its dose** (0.5 → 0.62 at t8) on the host,
> then washes out once the host is gone (0.62 → 0.01 — fixed host range, no coexistence
> oscillation). The commensal rebound is a real competitive release — E. coli grows only
> *after* the host clears, inheriting the freed glucose. Both the interference (phage) and
> the resource-release (metabolic competition) arms come out of the real stoichiometry.
> Only the two CarveMe GEMs need a reconstruction host.

## Honest scope

- **The phage parameters are literature-scale, not a measured isolate.** Burst ~80,
  latent ~0.5 h, etc. set the *tempo* of the crash; the claim "phage present → host down,
  commensal up" is robust to their exact values. They are ASSUMED — cite a specific
  phage's kinetics if you pin the scenario to one (`PHAGE` in `genome_scale.py`).
- **No evolution.** Host-range shifts, receptor-loss resistance, and CRISPR immunity are
  out of paradigm and *not* faked — the host range is fixed, which is why the phage
  clears its host completely rather than settling into a coexistence oscillation. See
  `docs/LIMITATIONS.md`.
- **Well-mixed.** Spatial refuge (a biofilm the phage cannot penetrate) is not
  represented; that belongs in `muode.spatial`.
