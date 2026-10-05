# µODE — Perturbation engine

The perturbation engine runs a simulation under a modified condition and shows how
the community responds. Typical uses include antibiotic treatment, gene knockouts,
and species removals.

> A perturbation changes reaction **bounds** at the start of the run. To introduce
> **biomass** at a chosen time instead (a transplant or probiotic dose), see
> [injection.md](injection.md); the two can be combined in one run. For
> *time-varying* drug exposure with spore survival and recurrence, use the
> `Antibiotic` ecology layer ([ecology.md](ecology.md)).

---

## How perturbations work

A perturbation is applied **at the start of the simulation** by modifying the flux
bounds of one or more reactions before the dFBA loop begins. The simulation then
runs normally, and the community response — growth-rate changes, metabolite
shifts, secondary extinctions — emerges from the dynamics.

This is a *first-order* model: the perturbation is instantaneous and its magnitude
is constant through the run. For a perturbation whose strength changes over time
(e.g. a drug that is gradually eliminated), use the `Antibiotic` ecology layer,
which adds one-compartment pharmacokinetics and an Emax kill model.

---

## Perturbation types

### Pathway attenuation (antibiotic model)

Reduces the flux capacity of every reaction in a named subsystem/pathway by a
fraction `efficacy` ($0 \le \text{efficacy} \le 1$):

$$v_{r,\max}' = (1 - \text{efficacy}) \cdot v_{r,\max} \quad \forall r \in \text{pathway}$$

An efficacy of 0 is no effect; 1 is a complete block (equivalent to a pathway
knockout). This is an interpretable model of an antibiotic that targets a specific
pathway (e.g. folate biosynthesis inhibited by trimethoprim).

```bash
muode perturb --community simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" \
              --efficacy 0.95 \
              --outdir results_antibiotic/
```

### Reaction knockout

Sets one or more named reaction fluxes to zero ($v_{r,\max} = v_{r,\min} = 0$).

```bash
muode perturb --community simulation_env/community.json \
              --knockout "DHFR,FOLD" \
              --outdir results_knockout/
```

Reaction ids must match the BiGG ids in the models.

### Species removal

Removes one or more species entirely before the simulation starts (initial biomass
set to zero and the organism removed from the LP). Useful for modeling a
species-level antibiotic effect or for identifying keystone species.

```bash
muode perturb --community simulation_env/community.json \
              --remove-species "Bacteroides_fragilis,Prevotella_copri" \
              --outdir results_removal/
```

---

## Secondary extinctions

The most informative result of a perturbation is often a **secondary extinction**:
a species the perturbation never touched goes extinct because a cross-feeding
partner it depended on was disrupted.

These emerge from the simulation — they are not hard-coded. When species A produces
metabolite M that species B needs, and you knock out A's pathway for M, species B
eventually starves. Use `result.extinct()` and `result.cross_feeding()` to
diagnose which species went extinct and why.

```bash
muode perturb --community simulation_env/community.json \
              --remove-species "Roseburia_intestinalis" \
              --outdir results/
```

---

## Combining perturbation types

Multiple types can be combined in one run; all modifications are applied before the
loop starts:

```bash
muode perturb --community simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" --efficacy 0.9 \
              --knockout "FOLD" \
              --remove-species "Prevotella_copri"
```

---

## Snakemake rule

The `perturb` rule is optional. Enable it with a `perturbation:` mapping in the
config:

```yaml
# config/config.yaml
perturbation:
  target_pathway: "Folate Biosynthesis"
  efficacy: 0.95
  # remove_species: "Bacteroides_fragilis"
  # knockout: "DHFR,FOLD"
```

Output: `{outdir}/perturbation/biomass.csv` (same format as the main simulation).

---

## Programmatic use

```python
from muode.dfba import DynamicFBA
from muode.perturb import Perturbation, PerturbationTarget

targets = [
    PerturbationTarget(subsystem="Folate Biosynthesis", efficacy=0.95),
    PerturbationTarget(reactions=["DHFR"], efficacy=1.0),
]
pert = Perturbation(name="antibiotic_combo", targets=targets)

engine = DynamicFBA(t_end=24.0, dt=0.1)
result = engine.run(community, diet, kinetics, perturbation=pert)
print(pert.describe())          # "Folate Biosynthesis @ 95% + DHFR knockout"
print(result.extinct())         # species that went extinct
```

---

## Output files

A perturbation run produces the same file set as a regular simulation, written to
`--outdir`:

| File | Content |
|------|---------|
| `biomass.csv` | Species biomass time course |
| `metabolites.csv` | Extracellular metabolite concentrations |
| `growth_rates.csv` | Per-species growth rates |
| `cross_feeding.csv` | Producer→consumer edges detected |
| `meta.json` | Run parameters |
| `biomass.png` | Stacked-area biomass chart |
| `metabolites.png` | Top metabolite trajectories |
| `cross_feeding.png` | Directed cross-feeding network |

Compare these with the unperturbed run to quantify the perturbation effect.
