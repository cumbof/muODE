#!/usr/bin/env python3
"""Run the genome-scale rCDI/FMT study -- the mechanism-decomposition experiment.

This is a WORKSTATION job.  Each arm integrates a 5-member genome-scale community over
~4 simulated days; on genome-scale LPs that is ~30 min per arm, and this runs five arms.
Nothing here needs a GPU or a cluster -- just time and one core (or -j for parallel arms).

WHAT IT ANSWERS
    Does the donor community resist colonization by C. difficile, and through which
    arm?  Every yield is gapseq stoichiometry, so the competition is decided by the
    genomes rather than by a dial; the arms then ablate one mechanism at a time:

        arm                what it isolates
        -----------------  ------------------------------------------------------------
        untreated          NOTHING done: does the pathogen colonize at all?  Read first
                           -- if it does not, no other arm means anything.
        abx_only           vancomycin alone: the standard of care, and the rCDI baseline
        fmt_only           FMT into an untreated gut: can donors displace an established
                           pathogen without the drug?
        fmt_full           drug + FMT: the treatment
        fmt_no_bile        the treatment WITHOUT the bile arm -> is bile load-bearing?
        fmt_no_ph          the treatment WITHOUT SCFA acidification
        fmt_competition    neither bile nor pH: pure nutrient competition

    Read `untreated` and `abx_only` before anything else.  An earlier version of this
    study dosed vancomycin in EVERY arm, so its "control" was really drug-without-FMT;
    the drug kills at ~4.5/h against a pathogen growing at 0.036/h, so it cleared the
    infection everywhere and all five arms reported CLEARED -- a result about the drug,
    read as a result about the community.

    Note what this roster can and cannot show.  The first genome-scale run found
    fmt_competition == abx_only and the obvious reading was that these donors compete
    badly -- but the cause was the DIET: its starch row failed to translate into the
    ModelSEED namespace, leaving a colonic medium with no fibre in it and every donor
    here is a fibre degrader.  With starch restored (media.py), B. theta overtakes the
    pathogen in monoculture (0.0416 vs 0.0359/h) and C. difficile gains nothing from it.
    Whether that advantage survives washout is what this run decides -- monoculture mu
    is not lambda.  If fmt_competition STILL looks no better than abx_only, that is the
    honest answer, not a bug.

RUN (on the workstation)
    python examples/fmt_cdiff/run.py --outdir results/fmt            # all seven arms
    python examples/fmt_cdiff/run.py --core --outdir results/fmt     # the decisive four
    # ~30 min per arm.  Start with --core.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import scenario as gs  # noqa: E402

#: arm name -> (fmt?, ablate)
#:
#: The first two arms are the ones the study was missing.  Every arm used to carry the
#: vancomycin course, so the "control" was really *drug, no FMT* -- and since the drug
#: kills at ~4.5/h against a pathogen growing at 0.036/h, it cleared the infection in
#: every arm before the community could matter.  A study of colonization resistance has
#: to contain an arm where nothing is done at all.
ARMS = {
    "untreated": (False, "abx"),          # no drug, no FMT: does the pathogen colonize?
    "abx_only": (False, ""),              # the standard of care -- and the rCDI baseline
    "fmt_only": (True, "abx"),            # FMT into an untreated gut
    "fmt_full": (True, ""),               # drug + FMT: the treatment
    "fmt_no_bile": (True, "bile"),
    "fmt_no_ph": (True, "ph"),
    "fmt_competition": (True, "bile ph"),
}

#: Arms worth running when you only want the question answered, in order of importance.
CORE_ARMS = ("untreated", "abx_only", "fmt_full", "fmt_competition")

#: Pathogen biomass (vegetative + spores) below which the infection is called cleared.
#: Shared with figures.py so the figure's threshold line and this verdict cannot drift.
#:
#: DEMOTED, deliberately.  This is an invented threshold, and crossed with a finite
#: horizon it produces a verdict about when we stopped looking rather than about the
#: biology: in the first working run every treated arm read CLEARED while the pathogen
#: was still growing exponentially at t_end.  At the measured rates abx_only crosses 0.1
#: at t~194 h and fmt_full at t~436 h -- same run, same biology, opposite verdicts,
#: chosen by the horizon.  Keep it as a coarse flag; read REBOUND_RATE for the result.
CLEARED_BELOW = 0.1

#: Hours at the end of the run over which the rebound rate is fitted.  Long enough to
#: average out the Euler wobble, short enough to sit well after the drug has cleared
#: (concentration ~0.02 by t=24 with a 2 h half-life).
REBOUND_WINDOW_H = 30.0


def rebound_rate(series, window_h: float = REBOUND_WINDOW_H) -> float:
    """Net specific growth rate (1/h) of the pathogen at the END of the run.

    THE outcome measure, and the reason it replaces ``cleared`` as the headline: it is
    free of both the threshold and the horizon.  A log-linear fit of the last
    ``window_h`` hours gives lambda = dln(X)/dt = mu - D - (losses), which is exactly the
    colonization-resistance question:

        lambda < 0   the community EXCLUDES the pathogen -- it washes out
        lambda ~ 0   it is held at a steady state
        lambda > 0   it is growing; the arm delays recurrence, it does not prevent it

    A final biomass answers "how far had it got when we stopped"; lambda answers "where
    is it going", which is the thing a clinician and a reviewer both actually want.
    """
    import numpy as np

    t = np.asarray(series.index, dtype=float)
    x = np.asarray(series.values, dtype=float)
    mask = (t >= t[-1] - window_h) & (x > 0)
    if mask.sum() < 3:
        return float("nan")
    # ln X ~ a + lambda t  -- the slope IS the net specific growth rate
    return float(np.polyfit(t[mask], np.log(x[mask]), 1)[0])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="results/fmt")
    ap.add_argument("--arms", nargs="+", choices=list(ARMS), default=list(ARMS))
    ap.add_argument("--core", action="store_true",
                    help=f"run only the decisive arms: {' '.join(CORE_ARMS)}")
    ap.add_argument("--t-end", type=float, default=120.0,
                    help="hours; default 120 = ~3 colonic transits (see build_scenario)")
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--dilution-rate", type=float, default=gs.DILUTION_RATE,
                    help="colonic washout 1/h; 0 restores the (unphysical) batch culture")
    args = ap.parse_args()
    arms = list(CORE_ARMS) if args.core else args.arms

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Report which members are actually driving this -- and refuse to run a CDI study
    # that is missing the pathogen, loudly, before spending hours on it.
    status = gs.readiness()
    print("community readiness:")
    for slug, st in status.items():
        print(f"  {'OK ' if st == 'ready' else '!! '} {slug:32} {st}")
    if status[gs.PATHOGEN] != "ready":
        print(f"\nABORT: pathogen {gs.PATHOGEN} is {status[gs.PATHOGEN]}.", file=sys.stderr)
        return 1

    summary = {"t_end": args.t_end, "dt": args.dt,
               "dilution_rate": args.dilution_rate, "arms": {}}
    for name in arms:
        fmt, ablate = ARMS[name]
        print(f"\n=== arm: {name}  (fmt={fmt}, ablate={ablate!r}) ===")
        t0 = time.time()
        res = gs.build_scenario(fmt=fmt, ablate=ablate, t_end=args.t_end,
                                dt=args.dt, dilution_rate=args.dilution_rate)
        secs = time.time() - t0

        res.biomass.to_csv(out / f"{name}_biomass.csv")
        res.metabolites.to_csv(out / f"{name}_metabolites.csv")
        if getattr(res, "spores", None) is not None:
            res.spores.to_csv(out / f"{name}_spores.csv")

        final_veg = float(res.biomass[gs.PATHOGEN].iloc[-1])
        min_veg = float(res.biomass[gs.PATHOGEN].min())

        # CLEARANCE MUST COUNT SPORES.  Sporulation is not clearance -- it is precisely
        # how C. difficile survives a drug course and comes back, so a metric that reads
        # only the vegetative pool scores dormancy as a cure and gets rCDI exactly
        # backwards.  The previous run declared CLEARED in all five arms while the
        # pathogen sat in a spore pool nothing ever looked at.
        final_spores = 0.0
        if getattr(res, "spores", None) is not None and gs.PATHOGEN in res.spores:
            final_spores = float(res.spores[gs.PATHOGEN].iloc[-1])
        final_total = final_veg + final_spores
        cleared = final_total < CLEARED_BELOW

        lam = rebound_rate(res.biomass[gs.PATHOGEN])

        latched = res.meta.get("spore_latched") or {}
        summary["arms"][name] = {
            "rebound_rate": lam,                    # THE outcome: 1/h, horizon-free
            "excluded": bool(lam < 0),              # lambda<0 = washout = resistance
            "pathogen_final": final_total,          # veg + spores: the reservoir
            "pathogen_final_vegetative": final_veg,
            "pathogen_final_spores": final_spores,
            "pathogen_min_vegetative": min_veg,
            "cleared": cleared,                     # coarse flag; see CLEARED_BELOW
            "spore_latched": latched,
            "runtime_s": round(secs, 1),
        }
        print(f"    pathogen: veg={final_veg:.4f}  spores={final_spores:.4f}  "
              f"total={final_total:.4f}  lambda={lam:+.4f}/h  "
              f"{'EXCLUDED' if lam < 0 else 'GROWING'}  ({secs/60:.1f} min)")
        if latched:
            print(f"    !! SPORULATION TRIGGER LATCHED for {list(latched)}: growth never "
                  f"reached mu_stress ({latched}). This arm cannot be read.",
                  file=sys.stderr)

    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # The headline table.  lambda is the result; the biomass columns are context.
    print("\n" + "=" * 78)
    print("MECHANISM DECOMPOSITION")
    print("=" * 78)
    print(f"  {'arm':20} {'total':>8} {'veg':>8} {'spores':>8} {'lambda/h':>9}  outcome")
    for name in arms:
        a = summary["arms"][name]
        print(f"  {name:20} {a['pathogen_final']:8.4f} "
              f"{a['pathogen_final_vegetative']:8.4f} {a['pathogen_final_spores']:8.4f} "
              f"{a['rebound_rate']:+9.4f}  "
              f"{'EXCLUDED' if a['excluded'] else 'GROWING'}")
    print("\n  lambda < 0: the community excludes the pathogen (washout) -- resistance.")
    print("  lambda > 0: it is still growing at t_end; the arm DELAYS recurrence.")
    print("  A high spore fraction is not a cure: it is the recurrence reservoir.")
    print(f"\nwrote {out}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
