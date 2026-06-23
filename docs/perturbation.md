# µODE — Perturbation Engine

The perturbation engine lets you run a simulation under a modified condition and
observe how the community responds. Common use cases include simulating
antibiotic treatment, gene knockouts, probiotic interventions, and species
removals.

---

## How perturbations work

A perturbation is applied **at the start of the simulation** by modifying the
flux bounds of one or more reactions before the dFBA loop begins. The simulation
then runs normally, and the community response (growth rate changes, metabolite
shifts, secondary extinctions) emerges from the dynamics.

This is a *first-order* model of a perturbation: it assumes the perturbation is
instantaneous and its magnitude is constant throughout the run. Time-varying
perturbations (e.g. an antibiotic that is gradually metabolised) would require
a custom hook into the loop — currently not implemented.

---

## Perturbation types

### Pathway attenuation (antibiotic model)

Reduces the flux capacity of every reaction in a named subsystem/pathway by a
fraction `efficacy` ($0 \le \text{efficacy} \le 1$):

$$v_{r,\max}' = (1 - \text{efficacy}) \cdot v_{r,\max} \quad \forall r \in \text{pathway}$$

An efficacy of 0 is no effect; an efficacy of 1 is a complete block (equivalent
to a pathway knockout). This is an interpretable model of an antibiotic that
targets a specific pathway (e.g. folate biosynthesis inhibited by trimethoprim).

```bash
muode perturb --community simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" \
              --efficacy 0.95 \
              --outdir results_antibiotic/
```

### Reaction knockout

Sets one or more named reaction fluxes to zero ($v_r = 0$, implemented as
$v_{r,\max} = 0$ and $v_{r,\min} = 0$).

```bash
muode perturb --community simulation_env/community.json \
              --knockout "DHFR,FOLD" \
              --outdir results_knockout/
```

Reaction ids must match the BiGG ids in the models.

### Species removal

Removes one or more species entirely from the community before the simulation
starts (sets initial biomass to zero and removes the organism from the LP).
Useful for modelling a species-level antibiotic effect or for studying which
species are keystones.

```bash
muode perturb --community simulation_env/community.json \
              --remove-species "Bacteroides_fragilis,Prevotella_copri" \
              --outdir results_removal/
```

---

## Secondary extinctions

The most scientifically interesting result of a perturbation is often a
**secondary extinction**: a species that was not directly targeted by the
perturbation goes extinct because a cross-feeding partner it depended on was
disrupted.

These emerge naturally from the simulation — they are not hard-coded. When
species A produces metabolite M that species B requires for growth, and you
knock out A's pathway for M, species B will eventually starve in the simulation.
Use `result.extinct()` and `result.cross_feeding()` to diagnose which species
went extinct and why.

```bash
muode perturb --community simulation_env/community.json \
              --remove-species "Roseburia_intestinalis" \
              --outdir results/
# then check the output
cat results/biomass.csv | head
```

---

## Combining perturbation types

Multiple perturbation types can be combined in one run:

```bash
muode perturb --community simulation_env/community.json \
              --target-pathway "Folate Biosynthesis" --efficacy 0.9 \
              --knockout "FOLD" \
              --remove-species "Prevotella_copri"
```

All three modifications are applied simultaneously before the loop starts.

---

## Snakemake rule

The `perturb` rule is optional. Enable it by providing a `perturbation:`
mapping in the config:

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
print(pert.describe())          # "Folate Biosynthesis @ 95 % + DHFR knockout"
print(result.extinct())         # species that went extinct
```

---

## Output files

The perturbation run produces the same file set as a regular simulation, written
to the `--outdir` directory:

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
