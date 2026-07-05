# rCDI / FMT ecological dynamics — a mechanistic muODE example

A worked, end-to-end demonstration of composing the **full muODE ecology stack**
to reproduce the mechanistic narrative of recurrent *Clostridioides difficile*
infection (rCDI) and its resolution by fecal microbiota transplant (FMT).

Every trajectory is a genuine emergent output of the muODE dynamic-flux-balance
engine and its composable ecology layers — antibiotic pharmacokinetics/dynamics,
bile-acid transformation, spore germination gating, SCFA/pH inhibition,
bacteriophage predation and nutrient competition. Nothing is scripted or drawn by
hand: change a guild's metabolism or an ecology constant and the ODE integration
— and the outcome — changes with it.

This is the *dynamics* companion to the genome-reconstruction pipeline in the
parent [`fmt_cdiff/`](../README.md) folder. Where the parent walks real genomes
through reconstruction → simulation, this folder runs a richer four-scenario
ecological analysis on transparent, guild-structured toy models so the mechanism
is visible (and runnable) on any machine, with no GEMs or solver licence needed.

---

## The four scenarios

| Script | What it runs | Emergent outcome |
|---|---|---|
| `antibiotic_relapse.py` | 10-day vancomycin course, **no** transplant | Commensal crash + drug PK/PD, bile shield lost, spores persist → germinate → **relapse** |
| `fmt_resolution.py` | Same course **+ a Day-12 FMT** | Donor engraftment, bile ratio reverts, *C. difficile* driven to **extinction** |
| `designed_consortium.py` | Naive bai-only core **vs.** a full 12-member consortium (into a depauperate lumen); plus the inferred cross-feeding network | Minimal core **fails**; rationally designed consortium **cures** |
| `failed_fmt_autopsy.py` | A successful FMT vs. two failure modes | Donor **metabolic insufficiency** and **bacteriophage predation** each abort engraftment |
| `common.py` | Thin adapter: re-exports the model from `muode.scenarios.rcdi` and adds the plotting style | *(imported by all scripts as `C`)* |
| `run_all.py` | Regenerate every figure + CSV in one command | — |

Each script writes a multi-panel PNG and the underlying CSV trajectories to
`./results/` (git-ignored).

---

## Run it

```bash
conda activate muode

# everything (~15 min; seven 42-day dynamic-FBA runs)
python run_all.py

# …or any single scenario
python antibiotic_relapse.py
python fmt_resolution.py
python designed_consortium.py
python failed_fmt_autopsy.py
```

The engine is deterministic, so reruns reproduce identical figures and CSVs.

---

## The model in brief

The model lives in the muODE package itself, at `muode/scenarios/rcdi.py`
(re-exported here through `common.py`); everything below is defined and
commented there.

**Members.** Each organism is a `LinprogOrganism` — a small linear program with
exchange reactions and one or more objective-weighted growth modes — grouped into
functional **guilds**: primary polysaccharide degraders (the trophic root),
`bai`-operon 7α-dehydroxylating effectors (build the secondary-bile shield),
amino-acid/Stickland scavengers (deplete the pathogen's proline/glycine),
butyrate/SCFA producers, a hydrogenotrophic sink, a fungal opportunist, and
*C. difficile* itself (a Stickland mode + a sugar mode).

**Ecology stack** (composable layers over the shared extracellular pool):
`WeakAcidInhibition` (SCFA→pH), `BileAcidTransform` (taurocholate→cholate→
deoxycholate, driven by guild biomass), `BileAcidInhibition` (graded suppression
of vegetative *C. difficile*), `SporeForming` (germination on primary bile,
blocked by secondary bile, re-sporulation under stress, slow washout),
`Antibiotic` (single-compartment vancomycin PK + Emax/EC50 PD), and — for the
autopsy — `PhageInfection` (host-biomass-coupled Levin–Stewart lysis).

**Environment.** A mucin-rich gut `Diet` with continuous influx and a
chemostat-style washout (`DILUTION`, colonic transit) that bounds every field to
a physiological steady state.

**Timeline.** 42 days: a 10-day course, an intervention at Day 12 (once the drug
has decayed below EC50), follow-up to Day 42. Scenarios differ *only* in whether
the resident community is present at t=0 and which members are delivered at
Day 12 — so outcomes are attributable.

The key emergent separations: relapse vs. FMT-cure and naive-vs-designed
consortium each differ by ~10⁴–10⁵× in residual pathogen burden at Day 42.

---

## Adapting this to real genomes (e.g. your own MAG cohort)

This example is deliberately synthetic so the mechanism is legible. To run the
**same ecological analysis on real assembled genomes** — say MAGs from an FMT
metagenomic study — the muODE engine and ecology stack stay *identical*; you
replace only the community and wire up the ecology roles. **A ready-to-edit
adapter, [`real_data_adapter.py`](real_data_adapter.py), does exactly this** —
it takes a directory of GEMs, an abundance TSV, and a role table, and builds the
very same ecology stack the toy scenarios use. The three steps below are what it
expects:

1. **Reconstruct GEMs.** Assemble + bin your metagenomes, then run the MAGs
   through the muODE reconstruction workflow (see the parent `fmt_cdiff/README.md`
   and `download_genomes.sh` for the pattern). This yields a `CobraOrganism` per
   MAG that the engine consumes exactly like the `LinprogOrganism`s here.

2. **Supply composition** via muODE's abundance contract — a two-column TSV of
   `(mag_id, relative_abundance)` (cf. `../donor_abundance.tsv`,
   `../recipient_abundance.tsv`). Profile your genomes across the samples to get
   these.

3. **Map genomes → ecology roles.** This is the real work and is **not**
   automatic. The ecology layers are keyed to functional roles — which MAG is the
   pathogen, which carry the `bai`/`bsh` bile operons, which are spore-formers,
   which are antibiotic-susceptible, and (for the phage layer) which host each
   phage infects. You record these in a small **role table**,
   [`roles.tsv`](roles.tsv) — one row per `mag_id` with 0/1 flags — derived from
   genome annotation and taxonomy (operon presence for `bai`/`bsh`, sporulation
   genes, resistance profile, and CRISPR-spacer / prophage matching for
   host–phage linkages). The adapter reads it and populates the layer id-sets.

Then run the whole thing on your data:

```bash
conda run -n muode python real_data_adapter.py \
    --models ./gems \                 # <mag_id>.xml GEMs from muode build/refine
    --abundance ../recipient_abundance.tsv \
    --roles roles.tsv \
    --donor-abundance ../donor_abundance.tsv \
    --arm fmt                         # relapse | fmt | phage
```

The **same string id** must appear in all three inputs for a member (GEM filename
stem = abundance id = role id); the adapter fails loudly on any mismatch. One id
space you *must* reconcile by hand: the carbon/amino-acid/SCFA metabolite ids in
the adapter's `CONFIG` block have to match your GEMs' exchange namespace (BiGG for
CarveMe, `cpd*****` for gapseq/ModelSEED) — otherwise members share no pool and no
cross-feeding (hence no colonisation resistance) can emerge.

In short: steps 1–2 are automated by the workflow; step 3 — the genome-to-role
mapping in `roles.tsv` — is the piece you provide from annotations. The dynamics,
once wired, are produced by the very same engine you see running here. Note that
the outcome is then a genuine *prediction*: the relapse-vs-cure separation appears
only if your assembled community actually carries the requisite functional
structure (e.g. a `bai`+ effector fed by primary degraders) — the shipped
`roles.tsv` deliberately lacks one, to show that you annotate by evidence.

> Note on antibiotic breadth: `muode/scenarios/rcdi.py` models the vancomycin course as broadly
> suppressive of the resident bacterial community (consistent with the profound,
> broad microbiome collapse reported clinically during vancomycin treatment of
> CDI), which is what leaves the depauperate post-antibiotic lumen the pathogen
> exploits. For a real cohort you would instead set susceptibility per MAG from
> its resistance profile.
