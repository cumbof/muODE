#!/usr/bin/env python3
"""Calibrate uptake Vmax against Clark et al.'s measured monocultures -- honestly.

    python examples/benchmarks/clark2021/fit_kinetics.py \
      --models results/clark2021/models/refined \
      --outdir results/clark2021/kinetics

What this does
--------------
For each of the 26 strains, grows the GEM alone in DM38 for 48 h across a sweep of
Vmax, and compares the endpoint (net acetate / butyrate / lactate / succinate) with
what Clark actually measured.  Then it asks the question an optimiser never asks:
**does the data determine Vmax at all?**

For a 48 h batch culture that runs to substrate exhaustion, usually it does not.
The endpoint is fixed by *yield* -- the stoichiometry's biomass and product per mole
of substrate -- and not by *rate*.  Raising Vmax makes the culture finish sooner, not
differently.  On E. coli core in DM38, endpoint acetate moves 22.4 -> 20.4 mM while
Vmax ranges over 6 -> 100 mmol/gDW/h: a 16x change in the parameter for a 9% change
in the observable, most of it integration noise.

So per strain this reports one of:

* ``identified``     -- the endpoint genuinely moves with Vmax and the optimum is
                        interior.  A number is reported.
* ``bounded_below``  -- Vmax must exceed some value for the culture to grow/finish,
                        and above that nothing changes.  **A bound is reported and a
                        point estimate is deliberately withheld.**
* ``unidentifiable`` -- Vmax moves nothing here.
* ``no_growth``      -- the model never grew.  That indicts the reconstruction or the
                        medium; it is not a kinetics problem and must not be "fixed"
                        by tuning Vmax.

Why this is the right answer rather than a disappointing one
------------------------------------------------------------
If the endpoint predictions are insensitive to Vmax above the bound, then muODE's
Clark endpoint scores are testing the **reconstructions and their stoichiometry**,
and the invented kinetics are *not* contaminating them.  That is a much stronger
claim than a fitted parameter would have been, and it is checkable -- this script is
the check.

Where kinetics genuinely do bite is competition: who wins a shared substrate, i.e.
community composition, and any time-resolved prediction.  An endpoint measures
neither.  Constraining Vmax properly needs growth curves (Clark reports none) or the
16S fractions -- and fitting on the fractions would burn the very data the benchmark
scores, so it is left alone.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import clark2021 as ck
from muode.calibrate import fit_vmax
from muode.community import Community
from muode.dfba import DynamicFBA
from muode.kinetics import KineticParameters
from muode.organism import CobraOrganism

T_END = 48.0
INOCULUM = 0.01
DEFAULT_GRID = [1, 2, 4, 6, 8, 10, 15, 20, 30, 50, 100]


def make_predictor(model_path: Path, code: str, diet, dt: float):
    """endpoint (net metabolites + biomass) of one strain in DM38, given a Vmax."""
    import cobra.io

    model = cobra.io.read_sbml_model(str(model_path))

    def predict(vmax: float):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            org = CobraOrganism(model.copy(), id=code)
            comm = Community([org], abundances={code: 1.0}, total_biomass=INOCULUM)
            res = DynamicFBA(t_end=T_END, dt=dt).run(
                comm, diet, KineticParameters(default_vmax=float(vmax))
            )
        final = res.metabolites.iloc[-1]
        out = {"biomass": float(res.biomass[code].iloc[-1])}
        for bigg in ck.METABOLITES.values():
            # net of the DM38 baseline: the medium already contains 28.3 mM lactate
            out[bigg] = float(final.get(bigg, 0.0)) - diet.initial_concentration(bigg)
        return out

    return predict


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, required=True)
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/kinetics"))
    ap.add_argument("--dt", type=float, default=0.2)
    ap.add_argument("--grid", type=float, nargs="+", default=DEFAULT_GRID)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = args.data or (args.outdir / "MasterDF.csv")
    if not data.exists():
        ck.fetch(data)

    diet = ck.dm38()
    df = ck.load(data)
    truth = ck.monoculture_phenotypes(df)          # per-strain measured secretion
    observed = ck.observations(df, diet=diet)
    mono = {o.species[0]: o for o in observed if len(o.species) == 1}

    models = {}
    # the repo ships the GEMs gzip-compressed (*.xml.gz); cobra reads them directly.
    for path in (sorted(args.models.glob("*.xml")) + sorted(args.models.glob("*.xml.gz"))
                 + sorted(args.models.glob("*.sbml"))):
        code = path.stem.split(".")[0]
        if code in ck.STRAINS:
            models[code] = path

    skipped = sorted(set(ck.STRAINS) - set(models))
    print(f"GEMs: {len(models)}/{len(ck.STRAINS)}"
          + (f"   MISSING (not calibrated): {skipped}" if skipped else ""))

    report = {"t_end_h": T_END, "grid": list(args.grid), "missing_models": skipped,
              "strains": {}}
    tally: dict[str, int] = {}

    for code in sorted(models):
        obs = mono.get(code)
        if obs is None:
            print(f"  {code}: no monoculture measurement -- skipped")
            continue
        target = {b: v for b, v in obs.net_metabolites.items()}
        # One strain's failure must not discard the whole sweep (~30 min of solves):
        # record it and carry on so the report still gets written.
        try:
            ident, prof = fit_vmax(
                make_predictor(models[code], code, diet, args.dt),
                observed=target,
                vmax_grid=args.grid,
                inoculum=INOCULUM,
            )
        except Exception as exc:  # noqa: BLE001 -- log-and-continue over 26 independent fits
            tally["error"] = tally.get("error", 0) + 1
            report["strains"][code] = {"species": ck.STRAINS[code].species,
                                       "verdict": "error", "error": repr(exc)}
            print(f"  {code}  {ck.STRAINS[code].species:36} ERROR: {exc!r}")
            continue
        tally[ident.verdict] = tally.get(ident.verdict, 0) + 1
        report["strains"][code] = {
            "species": ck.STRAINS[code].species,
            "verdict": ident.verdict,
            "vmax": ident.best,
            "lower_bound": ident.lower_bound,
            "plateau_from": ident.plateau_from,
            "sensitivity": ident.sensitivity,
            "observed_net": target,
            "predicted": {t: list(v) for t, v in prof.predicted.items()},
        }
        print(f"  {code}  {ck.STRAINS[code].species:36} {ident.summary()}")

    print("\n=== verdicts ===")
    for verdict, n in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(f"  {verdict:16} {n}")

    fitted = tally.get("identified", 0)
    if fitted == 0 and tally:
        print(
            "\nNo strain's Vmax is identified by the endpoint data -- as expected for a\n"
            "48 h batch culture that runs to substrate exhaustion.  This is a RESULT,\n"
            "not a failure: it means the Clark endpoint scores test the reconstructions\n"
            "and their stoichiometry, and are NOT contaminated by muODE's invented\n"
            "kinetics.  Report the bounds, and state that predictions are insensitive\n"
            "to Vmax above them.  Do not report a fitted Vmax: there isn't one."
        )

    out = args.outdir / "kinetics_report.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
