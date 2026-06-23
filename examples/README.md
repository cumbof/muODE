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

## Whole-pipeline smoke test (the `stub` engine)

`data/mags/` holds two tiny toy MAGs (`fermenter.fna`, `consumer.fna`) tagged
with `muode-stub:<role>` in their FASTA headers. The dependency-free `stub`
reconstruction engine turns them into simulatable placeholder GEMs, so the
*entire* Snakemake DAG runs on any machine — including ones where CarveMe and
CheckM2 cannot:

```bash
conda activate muode
snakemake --cores 4 --configfile config/config.demo.yaml   # note: no --use-conda
```

This reconstructs, refines, QCs (writing `reconstruction_summary.tsv`),
assembles and simulates the two stubs into the same glucose→acetate cross-feeding
community as the toy demo — end to end through real SBML files. It validates the
pipeline plumbing before you swap `engine: carveme` on a capable host. The stub
models are **not** science.

## Real GEMs

To run on genome-scale models, point the pipeline at a directory of MAG FASTA
files and let the Snakemake workflow reconstruct, refine, assemble and simulate
(see the top-level README "Quick Start" and `config/config.yaml`). The same
engine that runs the toy demo runs the real community unchanged.
