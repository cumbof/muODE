# µODE — Visualization and reporting

µODE generates static figures with **matplotlib** after a simulation run. All
rendering is headless-safe (no display required — the `Agg` backend is used), so
figures are produced on remote clusters and HPC systems without a graphical
environment.

Figures are written as PNG files to the simulation output directory and are
generated automatically by `muode simulate`, `muode demo`, `muode perturb`, and
`muode spatial`.

---

## Standard simulation figures

Produced by every call to `muode simulate` (and `muode demo`, `muode perturb`) via
`muode.viz.save_all`:

### `biomass.png`

A **stacked-area chart** of community composition over time. The y-axis is total
community biomass (gDW/L), with each species' contribution shown as a filled area.
Species are sorted by final biomass (descending); the top 20 are shown by default
(configurable by calling `plot_biomass` directly with `top=N`).

Useful for identifying which species dominate at steady state, spotting extinctions,
and seeing competitive-exclusion dynamics.

### `metabolites.png`

A **line chart** of the most dynamic extracellular metabolites over time. The y-axis
is concentration (mmol/L). The top 15 metabolites by range (max − min over the run)
are shown, which picks up the metabolites that change most (substrates depleted,
products accumulated).

Useful for identifying substrate depletion, cross-feeding metabolite dynamics
(acetate, butyrate, succinate, …), and diagnosing starvation events.

### `cross_feeding.png`

A **directed network graph** of the cross-feeding interactions detected in the run.
Nodes are species; directed edges are annotated with the exchanged metabolite. Edge
width is proportional to interaction strength (mean flux product). The layout uses a
spring algorithm (NetworkX).

Requires `networkx` in addition to `matplotlib`; if networkx is not installed the
network figure is skipped.

---

## Spatial simulation figures

Produced by `muode spatial` via `muode.viz.save_spatial`:

### `spatial_total_biomass.png`

Line chart of **total biomass per species over time** (summed across all grid
cells) — the same information as `biomass.png` but for a spatially structured
community.

### `spatial_final.png`

A grid of **heatmaps** (or 1D profiles for single-row grids) showing the final state:

- one panel per species — the spatial distribution of that species' biomass
  (gDW/cell) at the last time point;
- one panel per metabolite with non-trivial spatial variation (range across the grid
  > 1e-9 at the final frame) — the concentration heatmap.

This is the key diagnostic for spatial cross-feeding: higher consumer biomass near
the producer (where the diffused metabolite concentration is highest) and a
metabolite gradient across the grid.

---

## Programmatic use

```python
from muode import viz

# after a standard simulation:
viz.save_all(result, outdir="results/")

# individual figures:
viz.plot_biomass(result, "results/biomass.png", top=30)
viz.plot_metabolites(result, "results/metabolites.png", top=10)
viz.plot_cross_feeding_network(result, "results/cross_feeding.png")

# after a spatial simulation:
viz.save_spatial(spatial_result, outdir="results/spatial/")
viz.plot_spatial_final(spatial_result, "results/spatial/final.png")
viz.plot_spatial_total_biomass(spatial_result, "results/spatial/biomass.png")
```

---

## Installing matplotlib

matplotlib is an optional dependency:

```bash
pip install "muode[viz]"
# or
pip install matplotlib networkx
```

Without it, µODE still runs all simulations; figure output is simply skipped (with a
warning in the CLI).

---

## Self-contained HTML report

`muode report` bundles a results directory into a single, portable `report.html` —
summary tables (final composition, fold-changes, extinctions, the cross-feeding edge
list, the most dynamic metabolites and the run parameters) with the PNG figures
embedded as base64. One file, no server, no extra dependencies (pure stdlib), so a
whole run can be shared or archived as a single artefact.

```bash
muode simulate --community simulation_env/community.json --outdir results/
muode report --results results/                 # writes results/report.html
muode report --results results/ --out run42.html
```

Programmatically:

```python
from muode.report import build_report
build_report("results/", "results/report.html")
```

---

## Exploring results interactively

Every run saves its full data as CSV and (for spatial runs) `spatial.npz`, so a
results directory can be loaded directly into a notebook for ad-hoc exploration:

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# standard simulation
bio = pd.read_csv("results/biomass.csv", index_col=0)
bio.plot.area(figsize=(10, 5), title="Community biomass over time")
plt.show()

# spatial simulation
data = np.load("results/spatial/spatial.npz", allow_pickle=True)
species = list(data["species"])
fig, axes = plt.subplots(1, len(species))
for ax, sp in zip(axes, species):
    ax.imshow(data[f"biomass_{sp}"][-1], origin="lower", cmap="viridis")
    ax.set_title(sp)
plt.tight_layout()
plt.show()
```
