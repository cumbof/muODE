# Genome-derived per-strain mu_max (gRodon)

Predicts each Clark strain's physiological maximum growth rate from its genome alone
(codon usage of ribosomal proteins) -- the data-free source for muODE's mu_max cap,
so the Clark mechanism can transfer to scenarios without measured OD/composition.

Pipeline (genomes = data/clark2021/genomes/*.fna):
1. `fetch_ref.py`         -- E. coli K-12 ribosomal proteins (UniProt) as the HE-gene reference.
2. `annotate_ribosomal.py` -- pyrodigal gene-call + pyhmmer(phmmer) vs the reference ->
   tag ribosomal genes in each genome's CDS (~54/genome).
3. `predict_growth.R`    -- gRodon2 `predictGrowth(mode="full")` on the tagged CDS ->
   minimal doubling time -> mu_max = ln(2)/d.

Result: `../grodon_mumax.json` (median ~0.235/h; range 0.086-0.668).
Dependencies: Python (pyrodigal, pyhmmer) for steps 1-2; R (r-base,
bioconductor-cordon/biostrings, gRodon2) for step 3.
