#!/usr/bin/env python3
"""Run muODE against the Clark et al. 2021 measurements and report the score.

    # 1. genomes (26 strains, exact-strain matched)
    python examples/benchmarks/clark2021/fetch_genomes.py --outdir data/clark2021/genomes

    # 2. reconstruct GEMs *on DM38* -- the medium the cells were actually in
    snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \
      --config mags_dir=data/clark2021/genomes mag_extension=fna \
               diet=dm38 gapfill_media=diet outdir=results/clark2021

    # 3. score
    python examples/benchmarks/clark2021/run_benchmark.py \
      --models results/clark2021/refined --outdir results/clark2021/benchmark

Three tiers, reported separately on purpose
-------------------------------------------
**Tier 1 -- monoculture phenotype.**  Grow each strain alone on DM38 and ask which
fermentation products it secretes.  No community dynamics, no abundance
weighting, no kinetics worth arguing about: a failure here is a *reconstruction*
failure and nothing downstream can repair it.  A GEM of *B. thetaiotaomicron* that
secretes butyrate is simply wrong.  Get this right first.

**Tier 2 -- pairwise.**  Two strains: does the engine get the interaction?

**Tier 3 -- assembly.**  Up to 23 strains: composition and butyrate.

Reporting them together as one accuracy number would hide which layer is broken,
which is the entire diagnostic value of the dataset.

Honesty rails built in
----------------------
* Metabolites are scored as **net production** against the DM38 baseline, because
  the medium already contains 28.3 mM lactate (see `muode.benchmarks.clark2021`).
* Strains that did not grow in DM38 in vitro (FP) are excluded from tier 1 *by
  name*, and the exclusion is printed in the report rather than applied quietly.
* Communities whose members have no GEM are skipped and **counted**, so a
  benchmark run on half the panel cannot be mistaken for a good one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from muode.benchmarks import clark2021 as ck
from muode.community import Community
from muode.dfba import DynamicFBA
from muode.kinetics import KineticParameters
from muode.organism import CobraOrganism

#: The incubation the measurements came from: 48 h, sealed, anaerobic.
T_END = 48.0


def load_models(models_dir: Path) -> dict[str, Path]:
    """Map strain code -> GEM, by the filename stem the genomes were named with."""
    found = {}
    for path in sorted(models_dir.glob("*.xml")) + sorted(models_dir.glob("*.sbml")):
        code = path.stem.split(".")[0]
        if code in ck.STRAINS:
            found[code] = path
    return found


def simulate(codes, models, diet, kinetics, t_end=T_END, dt=0.1, n_jobs=1):
    """Run one community and return (net metabolite production, relative abundance)."""
    organisms = [CobraOrganism.from_file(str(models[c]), id=c) for c in codes]
    community = Community(organisms, abundances={c: 1.0 / len(codes) for c in codes},
                          total_biomass=0.01)
    result = DynamicFBA(t_end=t_end, dt=dt, n_jobs=n_jobs).run(community, diet, kinetics)

    final_mets = result.metabolites.iloc[-1].to_dict()
    net = {
        bigg: float(final_mets.get(bigg, 0.0)) - diet.initial_concentration(bigg)
        for bigg in ck.METABOLITES.values()
    }

    biomass = result.biomass.iloc[-1].to_dict()
    total = sum(max(0.0, v) for v in biomass.values())
    abundances = ({c: max(0.0, v) / total for c, v in biomass.items()} if total > 0
                  else {c: 0.0 for c in biomass})
    return net, abundances


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, required=True, help="directory of refined GEMs (<CODE>.xml)")
    ap.add_argument("--data", type=Path, default=None, help="Clark measurement table (downloaded if absent)")
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--max-richness", type=int, default=None, help="skip communities above this size")
    ap.add_argument("--tier", choices=["1", "2", "3", "all"], default="all")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = args.data or (args.outdir / "MasterDF.csv")
    if not data.exists():
        print(f"downloading the measurement table -> {data}")
        ck.fetch(data)

    diet = ck.dm38()
    kinetics = KineticParameters()
    df = ck.load(data)
    models = load_models(args.models)

    print(f"GEMs: {len(models)}/{len(ck.STRAINS)} strains -> {sorted(models)}")
    absent = sorted(set(ck.STRAINS) - set(models))
    if absent:
        print(f"MISSING GEMs for {absent} -- every community containing them is skipped")

    report: dict = {"n_models": len(models), "missing_models": absent, "t_end_h": T_END}

    # -- tier 1: monoculture secretion phenotype ----------------------------
    if args.tier in ("1", "all"):
        print("\n=== tier 1: monoculture phenotype on DM38 ===")
        truth = ck.monoculture_phenotypes(df)
        skip = ck.non_growers(df)
        predicted: dict = {}
        for code in sorted(models):
            net, _ = simulate([code], models, diet, kinetics, n_jobs=args.threads)
            predicted[code] = {b: v > 5.0 for b, v in net.items()}
            secretes = sorted(b for b, v in predicted[code].items() if v)
            print(f"  {code}  {ck.STRAINS[code].species:38} -> {secretes or '(nothing)'}")

        report["tier1"] = {
            bigg: ck.score_phenotypes(predicted, truth, bigg, exclude=skip)
            for bigg in ck.METABOLITES.values()
        }
        print(f"\n  excluded (did not grow in DM38 in vitro): {list(skip)}")
        for bigg, score in report["tier1"].items():
            print(f"  {bigg:10} accuracy={score['accuracy']:.2f}  F1={score['f1']:.2f}  "
                  f"(n={score['n_strains']})  wrong: {score['misclassified'] or 'none'}")

    # -- tiers 2/3: communities ---------------------------------------------
    if args.tier in ("2", "3", "all"):
        want = {"2": (2, 2), "3": (3, 99), "all": (2, 99)}[args.tier]
        observed = [
            o for o in ck.observations(df)
            if want[0] <= o.richness <= want[1]
            and all(c in models for c in o.species)
            and (args.max_richness is None or o.richness <= args.max_richness)
        ]
        print(f"\n=== tiers 2/3: {len(observed)} communities (richness "
              f"{min((o.richness for o in observed), default=0)}-"
              f"{max((o.richness for o in observed), default=0)}) ===")

        predicted_mets, abundance_error = {}, []
        for i, obs in enumerate(observed, 1):
            net, abundances = simulate(list(obs.species), models, diet, kinetics,
                                       n_jobs=args.threads)
            predicted_mets[obs.community] = net
            if obs.abundances:
                shared = [c for c in obs.species if c in obs.abundances]
                if shared:
                    abundance_error.append(sum(
                        abs(abundances.get(c, 0.0) - obs.abundances[c]) for c in shared
                    ) / len(shared))
            if i % 25 == 0 or i == len(observed):
                print(f"  {i}/{len(observed)}")

        report["community"] = {
            "n_communities": len(observed),
            "metabolites": {
                bigg: ck.score_metabolites(predicted_mets, observed, bigg, net=True)
                for bigg in ck.METABOLITES.values()
            },
            "abundance_mae": (sum(abundance_error) / len(abundance_error)
                              if abundance_error else None),
        }
        print()
        if not observed:
            print("  NO communities could be scored: every one contains a strain with no\n"
                  "  GEM.  This is not a result -- reconstruct the missing strains first.")
        for bigg, score in report["community"]["metabolites"].items():
            if "mae" not in score:                      # degenerate: too few pairs
                print(f"  {bigg:10} not scored ({score.get('error', 'no data')})")
                continue
            r = score.get("pearson_r")
            print(f"  {bigg:10} r={'n/a' if r is None else format(r, '.3f')}  "
                  f"MAE={score['mae']:.2f} mM  bias={score['bias']:+.2f}  (n={score['n']})")
        mae = report["community"]["abundance_mae"]
        if mae is not None:
            print(f"  relative-abundance MAE: {mae:.3f}")

    out = args.outdir / "benchmark_report.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
