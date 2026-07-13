#!/usr/bin/env python3
"""Can the models actually *eat* what a medium offers?

A diet row is a promise: "this nutrient is available."  The promise is empty
unless two things hold downstream:

1. **The CarveMe universe has an exchange for it.**  Gap-filling draws reactions
   from that universe.  If ``EX_<met>_e`` is not in it, no carved model can ever
   import the metabolite, and the diet row is a silent no-op -- it costs nothing,
   changes nothing, and looks like biology in the CSV.
2. **The carved models actually carry that exchange.**  The universe having it
   means a model *could*; only the models tell you whether they *do*.  A medium
   full of bile acids is worthless if none of the 89 MAGs got a bile-acid
   transporter from their genome annotation.

This script answers both, for any set of candidate ids.  Run it BEFORE adding a
metabolite to a diet -- that is the whole point.  It must run where CarveMe is
installed (i.e. NOT on aarch64).

    python examples/diets/check_medium_namespace.py \
        --models results_newdiet/models/draft \
        cholate dchac tchola gchola tdechola gdchola \
        pect xylan4 xylan8 bglc urea nh4 mqn7 mqn8 2dmmq8

With no ids given, it checks every metabolite in the named diet instead, which
is how you find rows that are already dead:

    python examples/diets/check_medium_namespace.py --diet western_gut \
        --models results_newdiet/models/draft
"""

from __future__ import annotations

import argparse
import glob
import gzip
import os
import re
import sys
from collections import Counter


def universe_exchanges() -> set[str]:
    """Every ``EX_*`` reaction id in the CarveMe universe, read straight from SBML.

    Parsed by regex rather than cobra: the universe is a large model and we only
    need reaction ids, so loading it into a solver-backed object would cost
    minutes to answer a question about strings.
    """
    try:
        import carveme
    except ImportError:
        sys.exit(
            "carveme is not importable -- run this on the workstation, not the "
            "aarch64 dev box (CarveMe has no aarch64 build)."
        )

    root = os.path.join(os.path.dirname(carveme.__file__), "data", "generated")
    candidates = sorted(glob.glob(os.path.join(root, "universe*.xml*")))
    if not candidates:
        sys.exit(f"no universe model found under {root} -- has CarveMe's layout changed?")

    # Prefer the bacterial universe (what `carve` uses by default); fall back to
    # whatever is there, and say which one we read.
    path = next((p for p in candidates if "bacteria" in p), candidates[0])
    print(f"universe: {path}")

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        text = fh.read()
    return set(re.findall(r'id="R_(EX_\w+)"', text))


def model_exchanges(model_dir: str) -> tuple[Counter, int]:
    """Count, across every SBML in ``model_dir``, how many models carry each EX_."""
    paths = sorted(glob.glob(os.path.join(model_dir, "*.xml")))
    if not paths:
        sys.exit(f"no *.xml models under {model_dir}")
    seen: Counter = Counter()
    for p in paths:
        with open(p) as fh:
            seen.update(set(re.findall(r'id="R_(EX_\w+)"', fh.read())))
    return seen, len(paths)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metabolites", nargs="*",
                    help="BiGG ids WITHOUT the _e suffix (e.g. dchac, pect, urea)")
    ap.add_argument("--diet", help="instead check every metabolite in this diet/preset")
    ap.add_argument("--models", help="directory of carved SBML models to scan too")
    args = ap.parse_args()

    if args.diet:
        from muode.diet import load_diet
        from muode.media import supplied_metabolites, to_compound
        diet = load_diet(args.diet)
        mets = sorted(to_compound(m) for m in supplied_metabolites(diet))
        print(f"diet {args.diet}: {len(mets)} supplied metabolites")
    elif args.metabolites:
        mets = args.metabolites
    else:
        ap.error("give some metabolite ids, or --diet")

    universe = universe_exchanges()
    print(f"          {len(universe)} exchange reactions\n")

    carried: Counter = Counter()
    n_models = 0
    if args.models:
        carried, n_models = model_exchanges(args.models)
        print(f"models:   {n_models} carved SBMLs in {args.models}\n")

    hdr = f"{'metabolite':<16} {'in universe':>11}"
    if n_models:
        hdr += f"{'models carrying it':>22}"
    print(hdr)
    print("-" * len(hdr))

    dead, unreachable = [], []
    for met in mets:
        ex = f"EX_{met}_e"
        in_uni = ex in universe
        row = f"{met:<16} {'yes' if in_uni else 'NO':>11}"
        if n_models:
            n = carried.get(ex, 0)
            row += f"{n:>13} / {n_models:<6}"
            if in_uni and n == 0:
                unreachable.append(met)
        if not in_uni:
            dead.append(met)
        print(row)

    print()
    if dead:
        print(f"NOT IN THE UNIVERSE ({len(dead)}) -- adding these to a diet does NOTHING:")
        print(f"  {', '.join(dead)}")
    if unreachable:
        print(f"IN THE UNIVERSE BUT NO MODEL CARRIES THE EXCHANGE ({len(unreachable)}):")
        print(f"  {', '.join(unreachable)}")
        print("  The medium *could* supply these; these genomes cannot take them up.")
        print("  Supplying them is honest but inert -- say so rather than claiming coverage.")
    if not dead and not unreachable:
        print("every metabolite is in the universe and reachable by at least one model.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
