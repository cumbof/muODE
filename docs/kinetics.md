# µODE — Kinetic Parameter Prediction and Enzyme Constraints

Kinetics in muODE play two distinct roles that must not be conflated:

| Role | Parameter | Where it acts | Layer |
|------|-----------|--------------|-------|
| Substrate uptake | $K_m$, $V_{max}$ | Exchange reactions — bounds the rate at which a species can take up each substrate from the extracellular pool | **Required** (has safe defaults) |
| Enzyme constraints | $k_{cat}$ | Intracellular reactions — caps flux through enzyme-limited steps | **Optional** (GECKO-lite layer) |

---

## Why kinetics matter in dynamic FBA

Standard FBA is agnostic to extracellular metabolite concentrations: uptake
fluxes are bounded by a constant $v_{max}$ regardless of how much substrate is
available. In dFBA this is unrealistic — at low concentrations a species should
take up substrates more slowly.

muODE couples each uptake reaction to the extracellular pool via Michaelis–Menten:

$$v^{up}_{j,i} = \frac{V_{max}\, M_j}{K_m + M_j}$$

This single equation makes nutrient depletion, substrate switching (e.g. the
*E. coli* diauxic shift), and competitive exclusion emerge naturally from the
simulation without hard-coding any of those behaviours.

---

## The `KineticParameters` store (`muode.kinetics`)

All kinetic data lives in a typed `KineticParameters` object that can be built
programmatically or loaded from JSON:

```python
from muode.kinetics import KineticParameters

# start from defaults (Vmax=10.0, Km=0.01 mmol/gDW/h for all reactions)
kin = KineticParameters(default_vmax=10.0, default_km=0.01)

# override Km for a specific uptake reaction
kin.set_km("EX_glc__D_e", 0.005)

# store a predicted kcat for an intracellular reaction
kin.set_kcat("PGI", 312.0)  # /s

# persist and reload
kin.to_json("model.kinetics.json")
kin2 = KineticParameters.from_json("model.kinetics.json")

# merge predictions from multiple files
merged = KineticParameters()
for f in Path("models/").glob("*.kinetics.json"):
    merged.merge(KineticParameters.from_json(f))
```

---

## Predictors

### Heuristic predictor (default, no extra dependencies)

The built-in `HeuristicPredictor` assigns kinetic parameters from literature
without requiring any ML runtime or internet access:

- **$K_m$:** look-up table of experimentally measured values for the ~50 most
  common substrates in BiGG-namespace models (glucose, acetate, oxygen, etc.).
  Unknown substrates get the median value (~0.1 mmol/L).
- **$k_{cat}$:** the genome-wide median $k_{cat}$ from Bar-Even et al. (2011),
  ~13.7/s, which is a defensible prior when no enzyme-specific data is available.

These are **placeholders**, not measurements. They preserve the qualitative
dynamics (substrate limitation, depletion) without over-fitting to unknown
parameters. The pipeline runs end-to-end with only the heuristic predictor.

```bash
muode refine --models ./models/draft_gems/ --outdir ./models/refined_gems/ \
             --predict-kinetics --predictor heuristic
```

### DLKcat ($k_{cat}$ prediction, opt-in)

DLKcat (Li et al., 2022) predicts enzyme turnover numbers from protein sequence
+ substrate SMILES using a graph neural network. Requires the `ml` extra:

```bash
pip install "muode[ml]"
```

DLKcat needs an **enzyme context**: a mapping from each BiGG reaction id to the
catalysing enzyme's amino-acid sequence and the substrate's SMILES string.
`muode.predict.build_enzyme_context` provides this mapping by cross-referencing
the BiGG database.

```bash
muode refine --models ./models/draft_gems/ --outdir ./models/refined_gems/ \
             --predict-kinetics --predictor dlkcat
```

### Kroll $K_m$ predictor (opt-in)

Kroll et al. (2021) trained a deep-learning model for $K_m$ prediction from
protein sequence + metabolite fingerprint. Same `ml` extra requirement and same
enzyme-context mechanism as DLKcat.

```bash
muode refine ... --predictor km-ml
```

---

## GECKO-lite enzyme constraints (`muode.enzyme`)

Beyond Michaelis–Menten uptake bounds, muODE supports a second, optional layer
of kinetic constraints on *intracellular* reaction velocities. This is a
simplified version of the GECKO framework (Sánchez et al., 2017):

$$|v_r| \le k_{cat,r} \cdot [E_r]$$

where $[E_r]$ is the effective enzyme concentration (estimated from proteome
fractions or set to a default). These bounds are applied once before the dFBA
loop starts and survive each per-step reset.

Enzyme constraints tighten the solution space significantly and can alter which
metabolic strategies are optimal. They are only meaningful when $k_{cat}$ values
are available (from DLKcat or a curated source).

```bash
muode simulate --community simulation_env/community.json \
               --kinetics models/refined_gems/ \
               --enzyme-constraints
```

The `apply_enzyme_constraints` method of `CobraOrganism` accepts a
`KineticParameters` object and a per-organism scaling factor. It returns a
report dict with `n_constrained` (number of reactions bounded).

### Protein-pool budget (GECKO/sMOMENT, optional)

Per-reaction caps bound each reaction independently. The **shared protein-pool
budget** instead forces all enzyme-constrained reactions to draw on *one* finite
enzyme-mass budget, so the cell must *allocate* limited protein between competing
pathways:

$$\sum_r \frac{MW_r}{k_{cat,r}\cdot 3600}\,\lvert v_r\rvert \;\le\; P
\qquad (\text{g enzyme / gDW})$$

$MW_r$ is the enzyme molecular weight (kDa = g/mmol; a typical default unless
overridden) and $P$ is the proteome fraction available to these reactions
(e.g. 0.2–0.5 g/gDW). This single budget is what makes **overflow metabolism**
(acetate/ethanol secretion at high growth) and substrate hierarchies emerge from
enzyme economics rather than being hard-coded. muODE implements it exactly,
without reaction splitting, via an auxiliary usage variable $a_r\ge\lvert v_r\rvert$
per reaction and one pool constraint added to the solver (so it survives every
per-step bound reset).

```bash
muode simulate --community simulation_env/community.json \
               --kinetics models/refined_gems/ \
               --protein-pool 0.2
```

`apply_protein_pool_constraint(model, kinetics, organism_id, pool_budget)`
(in `muode.enzyme`) adds the constraint and returns `n_pooled` (reactions drawing
on the pool) and the `pool_budget`. It is still GECKO-*lite* — a single global
budget with default masses, not a measured per-enzyme proteome allocation —
so calibrate $P$, $MW_r$ and $k_{cat,r}$ against data before any quantitative
claim.

---

## Config reference

```yaml
# config/config.yaml
predict_kinetics: false      # write per-MAG {stem}.kinetics.json during refine
predictor: "heuristic"       # heuristic (default) | dlkcat | km-ml
enzyme_constraints: false    # apply GECKO-lite kcat caps at simulation time
default_vmax: 10.0           # fallback Vmax (mmol/gDW/h) for all uptake reactions
default_km: 0.01             # fallback Km (mmol/L) for all uptake reactions
```

---

## References

- Bar-Even et al. (2011) *The moderately efficient enzyme: evolutionary and
  physicochemical trends shaping enzyme parameters.* Biochemistry 50(21):4402–4410.
- Li et al. (2022) *Deep learning-based kcat prediction enables improved enzyme-
  constrained model reconstruction.* Nature Catalysis 5:662–672.
- Kroll et al. (2021) *Deep learning allows genome-scale prediction of Michaelis
  constants from structural features.* PLOS Biology 19(10):e3001402.
- Sánchez et al. (2017) *Improving the phenotype predictions of a yeast genome-
  scale metabolic model by incorporating enzymatic constraints.* Mol. Syst. Biol.
  13(8):935.
