# µODE documentation

Per-feature documentation for µODE. Start with the [README](../README.md) for a
quick orientation, then dive into any feature below.

| Document | What it covers |
|----------|---------------|
| [engine.md](engine.md) | Dynamic-FBA engine: the SOA integrator, solver abstraction, cross-feeding |
| [reconstruction.md](reconstruction.md) | CarveMe / gapseq / offline / eukaryote (CarveFungi, ModelSEEDpy) engines, gene calling, namespaces, LP gap-filling, CheckM2 |
| [kinetics.md](kinetics.md) | Km/Vmax uptake bounds, kcat enzyme constraints, protein-pool budget, heuristic/DLKcat/Kroll predictors |
| [assembly.md](assembly.md) | Community manifest, abundance TSV, diet, subsampling |
| [perturbation.md](perturbation.md) | Antibiotic / knockout / species-removal perturbation engine |
| [injection.md](injection.md) | Timed biomass injection: transplants (FMT), probiotic doses, inoculation |
| [ecology.md](ecology.md) | Ecology layer: pH/SCFA, bile acids, sporulation, antibiotic PK, bacteriocins, oxygen, phage |
| [kingdoms.md](kingdoms.md) | Multiple kingdoms: fungi/protists (eukaryotes), phages (viruses), archaea |
| [validation.md](validation.md) | Benchmark framework, metrics (MAE, Spearman, F1), pass/fail report |
| [spatial.md](spatial.md) | 2D reaction–diffusion colony/biofilm engine |
| [workflow.md](workflow.md) | Snakemake pipeline, all rules, config reference, SLURM |
| [visualization.md](visualization.md) | Static figures and the self-contained HTML report |
