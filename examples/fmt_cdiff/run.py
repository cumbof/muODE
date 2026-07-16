#!/usr/bin/env python3
"""Run the genome-scale rCDI/FMT study -- the mechanism-decomposition experiment.

This is a WORKSTATION job.  Each arm integrates a 5-member genome-scale community over
~4 simulated days; on genome-scale LPs that is ~30 min per arm, and this runs five arms.
Nothing here needs a GPU or a cluster -- just time and one core (or -j for parallel arms).

WHAT IT ANSWERS
    The toy scenario reached its "FMT clears the pathogen" conclusion through nutrient
    competition, while the bile mechanism it advertised was decorative and the yields
    that decided the competition were tuned in (see tests/test_provenance.py).  This
    rerun removes the tuning -- every yield is now gapseq stoichiometry -- and decomposes
    the clearance by ablating layers:

        arm                what it isolates
        -----------------  ------------------------------------------------------------
        no_fmt             recurrence baseline (spores germinate in the wiped gut)
        fmt_full           the claim: donor community clears the pathogen
        fmt_no_bile        clearance WITHOUT the bile arm  -> is bile load-bearing now?
        fmt_no_ph          clearance WITHOUT SCFA acidification
        fmt_competition    neither bile nor pH: pure nutrient competition

    If fmt_full clears but fmt_competition does not, the mechanism is real and muODE has
    decomposed colonisation resistance into its parts -- the actual paper.  If
    fmt_competition clears just as well, the bile story is still decorative and that,
    too, is a publishable (and honest) finding.

RUN (on the workstation)
    python examples/fmt_cdiff/run.py --outdir results/fmt
    # ~2-3 h total.  Add --arms fmt_full no_fmt to run a subset first.
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
ARMS = {
    "no_fmt": (False, ""),
    "fmt_full": (True, ""),
    "fmt_no_bile": (True, "bile"),
    "fmt_no_ph": (True, "ph"),
    "fmt_competition": (True, "bile ph"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="results/fmt")
    ap.add_argument("--arms", nargs="+", choices=list(ARMS), default=list(ARMS))
    ap.add_argument("--t-end", type=float, default=96.0)
    ap.add_argument("--dt", type=float, default=0.05)
    args = ap.parse_args()

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

    summary = {"t_end": args.t_end, "dt": args.dt, "arms": {}}
    for name in args.arms:
        fmt, ablate = ARMS[name]
        print(f"\n=== arm: {name}  (fmt={fmt}, ablate={ablate!r}) ===")
        t0 = time.time()
        res = gs.build_scenario(fmt=fmt, ablate=ablate, t_end=args.t_end, dt=args.dt)
        secs = time.time() - t0

        res.biomass.to_csv(out / f"{name}_biomass.csv")
        res.metabolites.to_csv(out / f"{name}_metabolites.csv")
        if getattr(res, "spores", None) is not None:
            res.spores.to_csv(out / f"{name}_spores.csv")

        final_path = float(res.biomass[gs.PATHOGEN].iloc[-1])
        min_path = float(res.biomass[gs.PATHOGEN].min())
        cleared = final_path < 0.1
        summary["arms"][name] = {
            "pathogen_final": final_path,
            "pathogen_min": min_path,
            "cleared": cleared,
            "runtime_s": round(secs, 1),
        }
        print(f"    pathogen: min={min_path:.4f}  final={final_path:.4f}  "
              f"{'CLEARED' if cleared else 'PERSISTS'}  ({secs/60:.1f} min)")

    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # The headline table: does clearance survive removing the advertised mechanism?
    print("\n" + "=" * 64)
    print("MECHANISM DECOMPOSITION (pathogen final biomass per arm)")
    print("=" * 64)
    for name in args.arms:
        a = summary["arms"][name]
        print(f"  {name:20} {a['pathogen_final']:8.4f}   "
              f"{'CLEARED' if a['cleared'] else 'PERSISTS'}")
    print(f"\nwrote {out}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
