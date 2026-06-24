# µODE — Timed Biomass Injection (transplants, probiotics, inoculation)

A [perturbation](perturbation.md) changes the *bounds* of organisms already in a
community (an antibiotic, a knockout, a removal). An **injection** is the
complementary primitive: it changes the *state* of the community at a chosen time
by introducing biomass mid-run. This is what lets muODE model a discrete
ecological **event** rather than a static "who is present":

- a **fecal microbiota transplant (FMT)** — a whole donor community delivered as
  a bolus into a dysbiotic recipient at time *t*;
- a **probiotic dose** — one or a few strains added at a chosen time;
- a staged **inoculation** — seeding a second species once the first has
  conditioned the medium.

See the [`fmt_cdiff`](../examples/fmt_cdiff/) example for a full, realistic
application (FMT for recurrent *C. difficile* infection).

---

## How it works

Mechanistically, nothing special happens at the injection. The introduced
organisms are ordinary members of the [`Community`](assembly.md): they simply
start at **zero biomass** and have it set to the inoculum amount when the event
fires. The dynamic-FBA loop already idles any species below its biomass floor
(`min_biomass`), so an injected organism contributes no flux until its bolus
arrives — then it competes for the shared metabolite pool like every other
species.

An injection fires at the **first integration step whose time is ≥ `time`**. The
biomass is added *before* that step is recorded and solved, so:

- the bolus is visible in the recorded time course as a vertical jump, and
- the injected organisms are metabolically active from that step onward.

Whether an injected community **engrafts** — and whether it displaces an
incumbent (e.g. donor commensals restoring colonization resistance against a
pathogen) — is an **emergent** result of the simulation, not something the
injection hard-codes.

---

## The `Injection` object

```python
from muode import Injection

# explicit per-organism biomass (gDW/L) added at t = 24 h
inj = Injection(time=24.0, biomass={"B_theta": 0.03, "E_rectale": 0.02}, name="FMT")

# or from relative abundances scaled to a total inoculum
inj = Injection.from_abundances(
    time=24.0,
    abundances={"B_theta": 0.6, "E_rectale": 0.4},
    total_biomass=0.05,
    name="FMT",
)
inj.total()      # 0.05  (gDW/L introduced)
inj.describe()   # "FMT@24h (+2 taxa, 0.05 gDW/L)"
```

| Field | Meaning |
|-------|---------|
| `time` | when the bolus is added (h) |
| `biomass` | mapping organism id → biomass added (gDW/L) |
| `name` | label used in `meta` and CLI output |

If an injected id is **already present** in the community, the amount is *added*
to its current biomass (a boost — e.g. a transplant reinforcing a strain the
recipient already carries at low abundance). If the id is new, the dormant
zero-biomass member is seeded.

---

## Composing two communities — `merge_for_injection`

The usual pattern is to assemble a **recipient** and a **donor** as two
independent communities, then transplant one into the other:

```python
from muode import DynamicFBA, merge_for_injection

merged, event = merge_for_injection(recipient, donor, time=24.0, name="FMT")
# `merged`  : recipient members at their initial biomass + donor members dormant (0)
# `event`   : an Injection that introduces donor.initial_biomass() at t=24 h

result = DynamicFBA(t_end=96.0, dt=0.1).run(merged, diet, kinetics, injections=[event])
```

Donor members whose id already exists in the recipient are **not** duplicated:
the recipient's model is kept and the injection boosts that species at `time`.

---

## Running an injection

`DynamicFBA.run()` takes an `injections=` list (alongside the optional
`perturbation=`), so a transplant and a perturbation can be combined in one run:

```python
result = engine.run(community, diet, kinetics,
                    perturbation=pert,        # optional bound changes
                    injections=[event])       # optional timed state events
result.meta["injections"]   # ["FMT@24h (+9 taxa, 0.05 gDW/L)"]
```

---

## CLI

Both `muode simulate` and `muode perturb` accept donor manifests to transplant
in:

```bash
# transplant a donor community.json into a recipient at t = 24 h
muode simulate --community recipient_env/community.json \
               --inject   donor_env/community.json \
               --inject-time 24 \
               --time 96 --outdir results_fmt
```

`--inject` is repeatable; pass one `--inject-time` to apply the same time to all,
or one per `--inject` to stagger them:

```bash
muode simulate --community recipient.json \
               --inject probiotic_A.json --inject-time 12 \
               --inject probiotic_B.json --inject-time 36 \
               --time 72 --outdir results
```

### Manifest-level injections

Alternatively a single `community.json` can fully describe the event with an
`injections` block (the referenced organisms must be among the manifest's
`models`):

```json
{
  "models": ["...recipient and donor .xml..."],
  "abundances": { "C_difficile_630": 0.5, "...recipient...": 0.0 },
  "diet": "western_gut",
  "total_biomass": 0.02,
  "injections": [
    { "name": "FMT", "time": 24.0, "total_biomass": 0.05,
      "abundances": { "B_thetaiotaomicron_VPI5482": 0.20, "...donor...": 0.05 } }
  ]
}
```

Each entry accepts either an explicit `biomass: {id: gDW/L}` map or
`total_biomass` + `abundances` (scaled like `Injection.from_abundances`).

---

## Output files

An injection run produces the same file set as any simulation (see
[visualization.md](visualization.md)): `biomass.csv`, `metabolites.csv`,
`growth_rates.csv`, `cross_feeding.csv`, `meta.json`, and the PNG figures. The
injection appears as a jump in `biomass.png` at `time`, and is listed under
`meta["injections"]`.

To quantify an effect, run the simulation **with and without** the injection and
compare the two `biomass.csv` files (the [`fmt_cdiff`](../examples/fmt_cdiff/)
README walks through a control-vs-treatment comparison).
