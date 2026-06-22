# muODE examples

## Built-in cross-feeding demo (no data or solver licence needed)

```bash
muode demo --outdir results/demo
```

A two-species community:

- **A_glucose** ferments glucose and secretes acetate.
- **B_acetate** cannot use glucose; it grows *only* on the acetate A releases.

So B depends entirely on A through the shared metabolite pool — the canonical
cross-feeding motif muODE exists to capture. The run writes time-course CSVs, a
`cross_feeding.csv`, and figures, and prints the final community state. B blooms
*after* A has produced acetate.

### See a secondary extinction

```bash
muode demo --outdir results/demo_ko --remove A_glucose
```

Removing the feeder (A) starves the dependent species (B) even though nothing
touches B directly — B's biomass never grows. This is the perturbation engine's
core use case (an antibiotic that wipes one species cascading through the
cross-feeding network).

The demo is built from `muode.examples` using the dependency-light
`LinprogOrganism` backend, so it runs with only numpy/scipy.

## Real GEMs

To run on genome-scale models, point the pipeline at a directory of MAG FASTA
files and let the Snakemake workflow reconstruct, refine, assemble and simulate
(see the top-level README "Quick Start" and `config/config.yaml`). The same
engine that runs the toy demo runs the real community unchanged.
