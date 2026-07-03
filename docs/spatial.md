# µODE — Spatial Dynamic FBA (2D Reaction-Diffusion)

The spatial engine (`muode.spatial`) extends the community dFBA model to a 2D
grid, adding metabolite diffusion between grid cells. This enables simulation of
colony growth, biofilm formation, and spatially structured cross-feeding where
metabolite gradients form naturally.

---

## Physical model

### Grid

The simulation domain is a rectangular grid of `ny × nx` cells, each of size
`dx × dx` mm. Every grid cell contains its own biomass field (one value per
species) and metabolite concentrations. The cells are connected by diffusion.

### Per-cell FBA

At each time step, the engine solves a **separate FBA** for every grid cell:

1. For each cell $(y, x)$, set Michaelis–Menten uptake bounds from the local
   metabolite concentrations.
2. Solve each species' FBA in that cell (growth rate + exchange fluxes).
3. Integrate biomass and metabolites in that cell (Euler step).

This means the engine solves `ny × nx × n_species` FBAs per step, which can be
expensive for large grids.

### Metabolite diffusion

After each Euler step, metabolite concentrations diffuse across the grid using a
**finite-difference Laplacian** with **Neumann (no-flux) boundary conditions**
(zero gradient at the domain edges — equivalent to an impermeable wall):

$$\frac{\partial M_j}{\partial t} = D_j\, \nabla^2 M_j + \text{(FBA exchange fluxes)}$$

Discretised as:

$$M_j(y,x,t+\Delta t) = M_j(y,x,t) + D_j \frac{\Delta t}{\Delta x^2}
  \bigl[M_j(y+1,x) + M_j(y-1,x) + M_j(y,x+1) + M_j(y,x-1) - 4\,M_j(y,x)\bigr]$$

with reflected boundary conditions at the edges.

The diffusivity $D_j$ (mm²/h) is uniform by default (`default_diffusivity`);
per-metabolite values can be passed via the `SpatialDynamicFBA` constructor.

---

## Stability

The finite-difference diffusion scheme is conditionally stable: the time step
must satisfy the CFL condition $\Delta t \le \Delta x^2 / (4D)$. muODE warns
if this condition is violated and the step is set too large.

---

## The `SpatialDynamicFBA` class

```python
from muode.spatial import SpatialDynamicFBA

engine = SpatialDynamicFBA(
    nx=24,                    # grid width
    ny=1,                     # grid height (1 = 1D strip)
    dx=1.0,                   # cell size (mm)
    t_end=12.0,               # simulated time (h)
    dt=0.05,                  # Euler step (h)
    default_diffusivity=2.0,  # mm²/h for all metabolites
    n_jobs=1,                 # per-cell solver threads (-1 = all cores)
)
result = engine.run(community, diet, kinetics, inoculum=inoculum)
```

---

## Inoculum helpers

The inoculum specifies where each species starts in the grid.

```python
from muode.spatial import uniform_inoculum, halves_inoculum, point_inoculum

# all cells have the same starting biomass
inoculum = uniform_inoculum((ny, nx), species_ids, total_biomass=0.01)

# two species occupy opposite halves of the grid (demonstrates cross-feeding diffusion)
inoculum = halves_inoculum((ny, nx), id_left="producer", id_right="consumer",
                           amount=0.04, axis=1)

# a single cell is inoculated (colony growth from a point)
inoculum = point_inoculum((ny, nx), species_ids, cy=0, cx=12, amount=0.01)
```

The `inoculum` is a dict `{species_id: np.ndarray of shape (ny, nx)}`.

---

## `SpatialResult`

```python
result.biomass        # dict {species_id: np.ndarray of shape (n_frames, ny, nx)}
result.metabolites    # dict {metabolite_id: np.ndarray of shape (n_frames, ny, nx)}
result.total_biomass()  # pd.DataFrame (time × species), summed over grid
result.to_npz("results/spatial/")    # save compressed arrays
```

---

## CLI

```bash
# toy cross-feeding demo on a 1D strip (no data needed):
muode spatial --nx 24 --time 12 --outdir ./results/spatial/

# with real GEMs from a community manifest:
muode spatial --community simulation_env/community.json \
              --nx 30 --ny 30 --dx 0.5 \
              --time 24 --step 0.02 \
              --diffusivity 1.5 \
              --outdir ./results/spatial/
```

The CLI uses the **halves inoculum** for the toy demo (producer on the left
half, consumer on the right) to demonstrate that cross-feeding requires metabolite
diffusion across the spatial gradient.

---

## Spatial cross-feeding validation

The demo reproduces the expected spatial pattern: the consumer species grows
fastest **near the interface** with the producer (where the diffused metabolite
concentration is highest), not at the far edge. This gradient is a validation of
the diffusion model and confirms that spatial structure matters for cross-feeding
efficiency.

---

## Output files

| File | Content |
|------|---------|
| `spatial.npz` | Per-species biomass + per-metabolite field arrays (`n_frames × ny × nx`) |
| `spatial_total_biomass.csv` | Total biomass per species over time (summed over grid) |
| `spatial_total_biomass.png` | Time-course plot |
| `spatial_final.png` | Heatmaps of final biomass and metabolite fields |

---

## Performance notes

A `30 × 30` grid with 3 species and 20 metabolites over 24 h at `dt=0.05`
requires ~480 time steps × 900 cells × 3 FBA solves = ~1.3 million LP solves.
With HiGHS on a modern workstation each small LP takes ~0.1 ms, so this run
would take ~2 minutes. For larger GEMs or grids, the per-cell FBAs are the main
lever for speedup: set `n_jobs` (or `muode spatial -j -1`) to solve the populated
cells of each step across worker threads. The cells are independent, so the
fields are identical for any `n_jobs`; each worker holds private model copies, so
the gain scales with the grid and model size where the solver call dominates.
