# µODE — Validation framework

The validation module (`muode.validate`) scores a completed simulation against a
set of prior expectations about the community. It turns a run from an uncheckable
computation into a testable hypothesis.

---

## Concepts

A **benchmark expectation** is a YAML file that specifies what a correctly simulated
community should look like. It can include:

- expected **relative abundances** at the end of the simulation;
- expected **metabolite concentrations** (or ranges) at steady state;
- expected **cross-feeding edges** (which species exchanges which metabolite with
  which other species);
- **tolerances** for each metric.

The validator computes quantitative metrics and produces a pass/fail report.

---

## Benchmark expectation format

```yaml
# examples/benchmarks/my_community.yaml
name: "Gut community benchmark"
description: "Expected behaviour of a 4-species gut community on a Western diet"

relative_abundances:
  Bacteroides_fragilis: 0.40
  Roseburia_intestinalis: 0.25
  Prevotella_copri: 0.20
  Escherichia_coli: 0.15
  tolerance: 0.10        # allow ±10 percentage points

metabolites:
  butyrate_e:
    min: 5.0             # mmol/L
    max: 20.0
  acetate_e:
    min: 1.0
    max: 50.0

cross_feeding:
  - producer: Roseburia_intestinalis
    consumer: Bacteroides_fragilis
    metabolite: butyrate_e
  - producer: Escherichia_coli
    consumer: Roseburia_intestinalis
    metabolite: acetate_e

tolerances:
  abundance_mae: 0.15    # mean absolute error threshold
  spearman_rho: 0.70     # minimum Spearman correlation
  metabolite_rmse: 5.0   # mmol/L
  cross_feeding_f1: 0.50 # minimum F1 on cross-feeding edge set
```

---

## Metrics computed

### Relative-abundance accuracy

- **MAE** (mean absolute error): average of |predicted − expected| across species.
- **Spearman ρ**: rank correlation between predicted and expected abundances.

Final-frame biomass values are normalized to relative abundances before comparison.
Species present in the expectation but absent in the simulation count as zero.

### Metabolite accuracy

For each expected metabolite, the final simulated concentration is compared to the
expected range. The metric is **RMSE** across all expected metabolites.

### Cross-feeding edge precision/recall/F1

The simulated cross-feeding network (from `result.cross_feeding()`) is compared to
the expected edge set. Precision and recall are computed at the (producer, consumer,
metabolite) triplet level, and their harmonic mean gives F1. A high F1 means the
simulation reproduces the expected interactions; a low F1 suggests missing
cross-feeding or spurious interactions.

---

## CLI

```bash
muode validate \
  --results  ./results/ \
  --expected examples/benchmarks/my_community.yaml \
  --report   ./results/validation/report.json
```

Exit code 0 if all tolerances are met; exit code 1 if any fail.

### Example output

```json
{
  "name": "Gut community benchmark",
  "passed": true,
  "abundance": {
    "mae": 0.08,
    "spearman_rho": 0.91,
    "passed": true
  },
  "metabolites": {
    "rmse": 2.3,
    "passed": true
  },
  "cross_feeding": {
    "precision": 0.83,
    "recall": 0.75,
    "f1": 0.79,
    "passed": true
  }
}
```

---

## Snakemake rule (optional)

Set a benchmark YAML in the config to add a validation step to the workflow:

```yaml
# config/config.yaml
benchmark: "examples/benchmarks/my_community.yaml"
```

The `validate` rule runs after `simulate` and writes `{outdir}/validation/report.json`.
The rule's exit code reflects pass/fail, so a failed validation fails the Snakemake
job.

---

## Programmatic use

```python
from muode.validate import BenchmarkExpectation, validate_outputs, compare

# load an expectation from a YAML file
expectation = BenchmarkExpectation.from_file("examples/benchmarks/toy_cross_feeding.yaml")

# run the full validation against a results directory
report = validate_outputs("results/", expectation)
print(report.to_dict())
print("PASSED" if report.passed else "FAILED")

# or validate a single SimulationResult object
report = compare(result, expectation)
```

---

## Writing your own benchmark YAML

Start from the bundled cross-feeding example:

```bash
cp examples/benchmarks/toy_cross_feeding.yaml examples/benchmarks/my_run.yaml
# edit the species names, expected abundances and metabolites to match your community
muode validate --results ./results/ --expected examples/benchmarks/my_run.yaml
```

You do not need to fill in all sections. The validator skips any section not present
in the YAML; a benchmark with only `cross_feeding:` expectations is valid and useful.
