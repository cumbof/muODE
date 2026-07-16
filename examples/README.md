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

## Spatial (colony / biofilm) demo

```bash
muode spatial --nx 24 --time 12 --outdir results/spatial
```

Runs the 2D reaction-diffusion engine on the same two-species cross-feeding
community, but **spatially**: the glucose fermenter is inoculated in one half of
a strip and the acetate specialist in the other. The fermenter grows on the
ambient glucose and secretes acetate, which must *diffuse* across the grid before
the acetate specialist can use it — so the consumer grows fastest near the
interface and barely at the far edge (a spatial cross-feeding gradient). Writes
`spatial.npz`, `spatial_total_biomass.csv`, and field/time-course figures. The
engine reuses the same models, diet and kinetics as the well-mixed demo.

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

## Real GEMs — 12-species western gut community

`gut_western/` is a self-contained, end-to-end example using 12 public NCBI
RefSeq reference genomes representing the dominant phyla and metabolic guilds
of a human gut under a western dietary pattern:

| Taxon | Role | Rel. abundance |
|-------|------|---------------|
| *Bacteroides thetaiotaomicron* VPI-5482 | Primary carbohydrate degrader (PULs) | 0.22 |
| *Bacteroides fragilis* NCTC 9343 | Immunomodulatory Bacteroides | 0.13 |
| *Eubacterium rectale* ATCC 33656 | Butyrate producer (cross-feeds on acetate) | 0.12 |
| *Blautia obeum* A2-162 | Hydrogen-consuming acetogen | 0.08 |
| *Roseburia intestinalis* L1-82 | Butyrate from starch/arabinoxylan | 0.09 |
| *Faecalibacterium prausnitzii* A2-165 | Anti-inflammatory butyrate producer | 0.08 |
| *Ruminococcus bromii* L2-63 | Keystone starch degrader | 0.07 |
| *Bifidobacterium longum* NCC2705 | Ferments oligosaccharides → acetate + lactate | 0.07 |
| *Prevotella copri* DSM 18205 | Plant polysaccharide degrader (low in western diet) | 0.04 |
| *Akkermansia muciniphila* ATCC BAA-835 | Mucin degrader; reduced by high-fat diet | 0.04 |
| *Lactobacillus acidophilus* NCFM | Lactate producer from simple sugars | 0.03 |
| *Coprococcus comes* ATCC 27758 | Lactate → butyrate converter | 0.03 |

The abundance table in `gut_western/abundance.tsv` was built from the literature
to reflect the western gut microbiome signature (Bacteroides dominance,
reduced *F. prausnitzii* and *Akkermansia*, low *Prevotella*).

```bash
# Download 12 reference genomes from NCBI
bash examples/gut_western/download_genomes.sh

# Full pipeline (CarveMe + Prodigal required)
snakemake --use-conda --cores 8 --configfile examples/gut_western/config.yaml
```

For step-by-step CLI instructions, the species selection rationale,
cross-feeding network diagram, perturbation experiments (keystone species
removal), and guidance on interpreting the output figures, see
**[gut_western/README.md](gut_western/README.md)**.

## Fecal microbiota transplant (FMT) for recurrent *C. difficile*

`fmt_cdiff/` simulates an **FMT** as treatment for **recurrent *Clostridioides
difficile* infection** on **real genome-scale models**. The five members — the
pathogen and four donors — are reconstructed with **gapseq** (one ModelSEED
namespace, so they cross-feed with no id translation), and *no growth yield is a
dial*: each is the biomass stoichiometry that came out of the genome. The donor
community is transplanted into the recipient mid-simulation, and the experiment
**ablates** the bile and pH ecology layers to decompose *which* arm of
colonization resistance actually clears the pathogen:

```bash
# check readiness + the per-member yields that used to be dials (runs anywhere)
python examples/fmt_cdiff/scenario.py

# the mechanism-decomposition study (WORKSTATION: ~30 min/arm × 5 arms)
python examples/fmt_cdiff/run.py --outdir results/fmt
```

Whether the transplant clears *C. difficile* is an **emergent** outcome of the
shared-pool dynamic FBA, and the five arms (`no_fmt`, `fmt_full`, `fmt_no_bile`,
`fmt_no_ph`, `fmt_competition`) separate nutrient competition from SCFA
acidification and secondary bile acids. See
**[fmt_cdiff/README.md](fmt_cdiff/README.md)** for the community, the ablation
table, and how to read the result.

## Multi-kingdom community (bacteria + fungus + phage)

`multikingdom/` shows muODE handling the **non-bacterial** members of a gut
sample. A real metagenome is not only bacteria, and muODE's core is
*paradigm-defined, not taxon-defined*: a fungal GEM simulates with the same
engine as a bacterial one. The example demonstrates the two things that
genuinely differ between kingdoms:

| Member | Kingdom | Mechanism it adds |
|--------|---------|-------------------|
| *Bacteroides thetaiotaomicron* | bacteria (obligate anaerobe) | keystone fermenter |
| *Candida albicans* | **eukaryote** (fungus, facultative) | scavenges O₂ → keeps the niche anaerobic |
| *Klebsiella pneumoniae* | bacteria (facultative) | pathobiont; the phage host |
| *vB_Kpn* | **virus** (phage) | lytic predator (coupled infection ODE) |

```bash
# any machine — toy-model mechanistic demo (no GEMs, no CarveMe)
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/multikingdom/mechanistic_demo.py
```

It runs the community and two ablations: drop the **fungus** and the O₂ it
scavenged poisons the obligate-anaerobe keystone; drop the **phage** and the
*Klebsiella* pathobiont blooms unchecked. Both are emergent. The example is
explicit about reconstruction routing (CarveMe cannot build the fungus; the phage
has no GEM at all) and about what stays out of scope (phage/host *evolution*).
See **[multikingdom/README.md](multikingdom/README.md)**.

## Strain competition (one species, three strains)

`strain_competition/` asks **which strain of the same species survives** in a
shared environment. Three *E. coli* strains — differing only in accessory traits
— start at equal biomass on one glucose pool, and the winner depends on which
competition dominates:

| Strain | Differs by | |
|--------|------------|--|
| `Ecoli__glc_specialist` | best glucose grower (lowest R\*) | wins on resources |
| `Ecoli__colicinogenic` | secretes a **colicin** | wins by interference |
| `Ecoli__arabinose_user` | has a **private sugar** | coexists by niche |

```bash
# any machine — toy-model mechanistic demo (no GEMs, no CarveMe)
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/strain_competition/mechanistic_demo.py
```

Three emergent outcomes: on glucose alone the specialist **competitively
excludes** the rest (Gause / Tilman R\*); adding the colicin (`muode.Bacteriocin`)
lets the *weaker* producer **overturn** the specialist (interference beats
exploitation); giving one strain a private substrate yields **coexistence** (niche
differentiation). Its key prerequisite — and honest caveat — is that the inputs be
**strain-resolved** genomes: the `../metagenomics/` catalogue's dRep step (0.95
ANI) would merge the strains, so use reference isolates or strain-resolved
metagenomics. See **[strain_competition/README.md](strain_competition/README.md)**.
