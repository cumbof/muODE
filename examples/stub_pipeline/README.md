# Stub pipeline — the plumbing smoke test

This is the one example that is **deliberately not science**. It exists to prove the
*whole* Snakemake DAG — reconstruct → refine → QC → assemble → simulate → validate —
executes end to end on **any** machine, including ones where CarveMe, gapseq and
CheckM2 cannot run (no solver licence, no bioinformatics stack).

`data/mags/` holds two toy MAGs, `fermenter.fna` and `consumer.fna`, tagged
`[muode-stub:<role>]` in their FASTA headers. The dependency-free `stub`
reconstruction engine reads that tag and emits a placeholder GEM for each, giving the
same glucose → acetate cross-feeding pair as `muode demo` — but arriving there through
real SBML files and every rule in the workflow.

```bash
conda activate muode
snakemake --cores 4 --configfile config/config.demo.yaml   # note: no --use-conda
```

The run is checked against `../benchmarks/stub_pipeline.yaml`, which asserts the
intended biology (both members ≈ 0.5, acetate produced then cross-fed down). That
makes this a **regression test for the pipeline**, not a claim about biology.

> **The stub models are not science.** Two contigs of filler sequence carry no
> metabolism; the placeholder GEMs are hand-built to ferment and cross-feed. A green
> run here says the plumbing is sound and nothing more. For results that mean
> something, see `../fmt_cdiff/` (real gapseq GEMs) or swap `engine: carveme` on a
> capable host.
