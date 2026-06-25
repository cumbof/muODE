# µODE — Scientific Evaluation & Revised Plan

*Assessment of the original concept, the corrections applied, and the resulting
implementation plan.*

## TL;DR verdict

**The idea is reasonable and scientifically sound.** The vision — go from
metagenome-assembled genomes (MAGs) to a *dynamic, predictive* model of community
metabolism — sits squarely on the cutting edge of microbial-community systems
biology, and every phase maps onto real, peer-reviewed methods. It is ambitious
but not fantasy.

However, the original README contained **five technical inaccuracies / gaps** that
would have produced a pipeline that either doesn't run as described or is
scientifically misleading. None are fatal; all are now corrected in the plan and
(for the core) in code:

1. **MICOM is not a dynamic-FBA engine.** It computes a *steady-state* community
   solution, not a time course. The dynamics had to be supplied by a real dFBA
   integrator.
2. **kcat and Km play different roles** and were conflated. Only Km/Vmax enter the
   dynamic loop; kcat belongs to a separate enzyme-constrained layer.
3. **Namespace consistency was unaddressed.** Mixing CarveMe (BiGG) and gapseq
   (ModelSEED) models in one community breaks the shared metabolite pool.
4. **The hard Gurobi/CPLEX requirement** is an unnecessary barrier; a free solver
   default is now first-class.
5. **No quality control or validation** was specified — without it the output is
   an unfalsifiable simulation.

The rest of this document justifies each point and lays out the revised, phased
plan. A **working, tested core engine already exists** (see "What works today").

---

## Phase-by-phase assessment

### Phase 1 — Reconstruction (CarveMe / gapseq) ✅ sound, with a namespace caveat
Top-down (CarveMe) and bottom-up (gapseq) automated reconstruction are the right,
standard tools. The important subtlety the README missed:

- A **community** simulation requires every model to share one metabolite
  namespace, because species interact through a *common pool* of extracellular
  metabolites keyed by id. CarveMe emits **BiGG** ids; gapseq emits **ModelSEED**
  ids. **You cannot mix them** in one community without an explicit harmonisation
  layer (MetaNetX/ModelSEED cross-refs).
- **Decision:** default to **CarveMe** (fast — seconds to minutes per genome;
  consistent BiGG namespace; integrates with cobra/MICOM). Keep gapseq as an
  *alternative engine* selected per-run, never mixed. gapseq is too slow (hours
  per genome) to be the default at the ">500 MAGs" scale the README targets.
- **Added:** an explicit gene-calling step (Pyrodigal/Prodigal) before CarveMe.

### Phase 2 — Gap filling (MetaPathPredict + LP) ✅ sound, integration cost flagged
- Algorithmic **LP gap filling** (minimal reaction set to enable biomass on a
  medium) is standard and now *implemented* on top of cobra.
- **MetaPathPredict** (Geller-McGrath et al., 2024) genuinely predicts KEGG-module
  presence in incomplete MAGs. The honest caveat: it works in the **KEGG-module**
  namespace, while CarveMe models are **BiGG**. Turning a predicted module into
  reactions the LP can add requires a KEGG→BiGG/ModelSEED mapping. That mapping is
  the real engineering cost, so **AI gap filling is an opt-in refinement**, not a
  blocker for an end-to-end run.

### Phase 3 — Kinetic prediction (DLKcat / Km) ⚠️ needed disambiguation
This was the most conceptually muddled phase. FBA is constraint-based and does not
use Michaelis–Menten rate laws *internally*. So where do predicted kinetics go?

- **Km + Vmax → substrate uptake bounds.** In dynamic FBA the *upper bound* of each
  uptake reaction is recomputed every step from the extracellular concentration:
  `v = Vmax·M/(Km+M)`. **This** is where Km matters, and it is the classic coupling
  between the intracellular LP and the extracellular ODEs (Mahadevan 2002).
- **kcat → enzyme-constrained internal velocities.** `Vmax = kcat·[E]` bounds an
  *intracellular* reaction, the domain of enzyme-constrained models (GECKO,
  sMOMENT). This is a **separate, optional layer**, not part of the basic loop.
- **Decision:** the engine *requires* only Km/Vmax (with permissive defaults that
  recover substrate-unlimited dFBA). DLKcat/Km predictors are wrapped as an
  optional layer that *sharpens* these bounds. The pipeline runs end-to-end without
  ML; ML is a refinement. The README's `DLERKm` appears to be an approximate name —
  the relevant published Km predictor is Kroll et al. (2021).

### Phase 4 — Dynamic community simulation ❗ the key correction
> "Integrates individual GEMs … using the MICOM framework … Uses ODEs to
> dynamically update … at discrete time steps."

**MICOM does not do this.** MICOM (Diener et al., 2020) solves a *single*
regularised **steady-state** community FBA (a cooperative-tradeoff snapshot of who
grows on what). It is excellent for one time point but is **not a time
integrator.** The README's own math — `dX/dt = µX`, `dM/dt = Σ vX` — is the
**Static Optimization Approach (SOA)** to dynamic FBA (Mahadevan, Edwards & Doyle,
*Biophys. J.* 2002), which is a *different* algorithm from MICOM.

- **Decision:** muODE ships a **native SOA dynamic-FBA integrator** (implemented,
  tested) that loops: set Michaelis–Menten uptake bounds → solve each species' FBA
  → integrate the extracellular ODEs one Euler step → repeat. This *is* the
  dynamics. A MICOM/cobra community model can be plugged in as the per-step solver;
  **COMETS** (a purpose-built dynamic community-FBA engine) is a planned alternative
  backend. MICOM is retained for what it is good at: steady-state cross-feeding
  snapshots and abundance-constrained community growth.

### Phase 5 — Perturbation engine ✅ sound and implemented
Modelling an antibiotic as a capacity reduction (`Vmax ← (1−efficacy)·Vmax`) on a
target pathway, a reaction knockout, or a species removal is a defensible
first-order approach. Crucially, the interesting result — **secondary extinctions
from disrupted cross-feeding** — is then an *emergent* property of the dynamic run,
not something hard-coded. This is implemented and demonstrated in tests.

### Mathematical framework ✅ correct
The FBA LP and the dFBA ODEs are stated correctly. The implementation adds the
practical details the README omits: a CFL-style stability cap so a step cannot
over-deplete a metabolite, optional cell death and chemostat dilution, and
non-negativity clamping.

---

## Cross-cutting additions

- **Solver:** default to **HiGHS** (free, fast, bundled via scipy/optlang).
  Gurobi/CPLEX remain drop-in for very large communities. Removes the hard
  commercial-licence barrier.
- **MAG QC gate (Phase 0):** CheckM2 completeness/contamination filtering *before*
  modelling — low-quality bins make bad models.
- **Model QC:** fast structural checks (mass/charge balance, and an
  **energy-generating-cycle** test that catches the "ATP from nothing" artefact),
  plus optional **memote** reports.
- **Validation:** sanity targets (e.g. short-chain-fatty-acid production ranges for
  gut communities) so results are falsifiable — currently a roadmap item.
- **Reproducibility:** Snakemake DAG, per-rule conda envs, pinned config, a tested
  toy example, and deterministic runs.

## Deferred / re-scoped (kept in the roadmap, not the core path)

- gapseq as default → optional alternative engine.
- AI kinetic prediction → optional refinement, not a prerequisite.
- ESM-2 kinetics, a live interactive GUI dashboard, real GPU DLKcat/Kroll
  checkpoints, and an in-process COMETS backend → genuinely future / external.
  (Spatiotemporal PDE biofilms, the protein-pool budget, the COMETS *export*
  bridge and a static HTML report are now implemented.)

---

## Principal risks (honest)

1. **Computational scale.** A ">500 MAG" community is a very large repeated LP.
   Mitigations: abundance thresholding/subsampling, HiGHS/Gurobi, parallel
   per-species solves, and (future) a COMETS backend.
2. **Parameter & prediction uncertainty.** AI-predicted kcat/Km carry real error;
   treat them as priors, report sensitivity, never as ground truth.
3. **Compounding model error.** Draft GEMs from incomplete MAGs + gap filling +
   predicted kinetics compound. QC gates and validation are the antidote; outputs
   are hypotheses to test, not measurements.

---

## Revised roadmap (supersedes the README checklist)

**Milestone 0 — Engine & scaffold (DONE)**
- [x] Native SOA dynamic-FBA community integrator (cobra + linprog backends)
- [x] Perturbation engine (knockout / pathway attenuation / species removal)
- [x] Diet/medium + Michaelis–Menten kinetics layer with safe defaults
- [x] Cross-feeding inference, time-course + figure outputs
- [x] CLI (`build/refine/assemble/simulate/perturb/qc/demo`)
- [x] Snakemake workflow (per-MAG fan-out, QC checkpoint, SLURM profile)
- [x] Test suite (toy cross-feeding + real *E. coli* core dFBA)

**Milestone 1 — Reconstruction at scale (in progress)**
- [x] End-to-end per-MAG DAG fan-out; a dependency-free `stub` engine runs the
      *whole* pipeline (reconstruct → refine → QC → assemble → simulate) on any
      architecture, validated locally on bundled toy MAGs
- [x] `reconstruction_summary.tsv` aggregation across all MAGs (one QC table per run)
- [x] QC-driven failure isolation: non-simulatable models are dropped from the
      community rather than aborting the run
- [x] Optional memote rule; LP gap-fill universal-model wiring (`universal_model`)
- [x] CheckM2 kept in the pipeline but gated by `run_checkm2` (off where it can't run)
- [ ] Validate CarveMe + Prodigal + CheckM2 on a real MAG set on a capable host
- [ ] Multi-threading / resource tuning for very large clusters (>500 MAGs)

**Milestone 2 — Kinetics refinement (in progress)**
- [x] Kinetic-parameter store extended with per-reaction kcat + JSON persistence
- [x] Dependency-free `heuristic` predictor (default): literature Km for common
      substrates, kcat around the genome-wide median — runs everywhere
- [x] GECKO-lite enzyme-constraint layer (`muode.enzyme`) wired through the CLI
      (`simulate --enzyme-constraints`), the workflow, and the dFBA loop
- [x] Opt-in DLKcat / Kroll-Km wrappers (`ml` extra) + `build_enzyme_context`
      BiGG→(sequence, SMILES) namespace mapping
- [x] Shared protein-*pool* GECKO/sMOMENT budget (`apply_protein_pool_constraint`,
      `simulate --protein-pool`): one finite enzyme mass allocated across all
      kcat-constrained reactions
- [ ] Bundle/validate real DLKcat & Kroll checkpoints on a GPU host *(external)*
- [ ] Measured per-enzyme proteome allocation (per-enzyme MW + abundances); ESM-2
      embeddings *(needs proteomics / a GPU host)*

**Milestone 3 — Validation & scale (in progress)**
- [x] Benchmark/validation framework (`muode.validate`, `muode validate`):
      relative-abundance MAE + Spearman, metabolite/SCFA error, cross-feeding
      edge precision/recall/F1, pass-fail vs. explicit tolerances; wired as an
      optional workflow `validate` rule
- [x] Abundance-aware subsampling (`muode.subsample`, `Community.subsample`:
      top-N / min-abundance / cumulative-coverage) for huge communities
- [x] Upstream abundance contract simplified: muODE accepts any 2-column TSV
      (`mag_id, rel_abundance`) produced by whatever upstream profiler the user
      ran (MetaSBT, Bracken, Kraken2, custom scripts, manual). muODE is agnostic
      to the profiling method — the TSV is the interface.
- [x] COMETS **export bridge** (`muode export-comets`, `muode.comets`): writes a
      COMETS layout + params + `cometspy` driver so the same community can be run
      in the independent COMETS engine for cross-validation
- [ ] Benchmark against real synthetic/gut datasets on a capable host *(external)*
- [ ] COMETS as an *in-process* per-step backend (the solver abstraction is the
      integration point; deferred — heavy external Java/solver dependency)

**Milestone 4 — Reach (in progress)**
- [x] Spatiotemporal (PDE) 2D reaction-diffusion dynamic-FBA engine
      (`muode.spatial`): per-cell community FBA + metabolite diffusion under
      no-flux boundaries, with spatial-cross-feeding tests and figures
- [x] Self-contained static HTML report (`muode report`, `muode.report`): one
      portable `report.html` per run (composition, cross-feeding, metabolites,
      run metadata, embedded figures)
- [ ] Live, server-backed interactive dashboard (deferred — a separate web app,
      not headless-testable; the static report + CSV/`.npz` are its data)
- [ ] Performance for large GEMs on large grids; ESM-2 kinetics *(GPU host)*

---

## Key references

- Mahadevan, Edwards & Doyle (2002) *Dynamic FBA of diauxic growth in E. coli.*
  Biophys. J. — the SOA dynamic-FBA algorithm muODE implements.
- Machado et al. (2018) *Fast automated reconstruction (CarveMe).* NAR.
- Zimmermann et al. (2021) *gapseq.* Genome Biology.
- Geller-McGrath et al. (2024) *MetaPathPredict.* Nature Communications.
- Li et al. (2022) *DLKcat.* Nature Catalysis.
- Kroll et al. (2021) *Deep learning of Km.* Nature Communications.
- Bar-Even et al. (2011) *The moderately efficient enzyme.* Biochemistry — source of
  the genome-wide median kcat used by the heuristic predictor.
- Sánchez et al. (2017) *GECKO: enzyme-constrained models.* Mol. Syst. Biol.
- Diener, Gibbons & Resendis-Antonio (2020) *MICOM.* mSystems.
- Dukovski et al. (2021) *COMETS protocol.* Nature Protocols.
- Lieven et al. (2020) *memote.* Nature Biotechnology.
- Chklovski et al. (2023) *CheckM2.* Nature Methods.
