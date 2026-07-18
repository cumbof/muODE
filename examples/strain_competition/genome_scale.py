#!/usr/bin/env python3
"""Genome-scale strain competition: three real E. coli strains, three modes.

The genome-scale counterpart of ``mechanistic_demo.py`` (the toy).  Same question --
which strain of one species wins a shared habitat, and how the answer flips with the mode
of competition -- but on REAL CarveMe reconstructions, so who wins on nutrients is an
OUTPUT of the genomes, not the toy's assigned label.

    # 1. genomes + GEMs (workstation) -- strain-resolved reference isolates
    bash examples/strain_competition/download_genomes.sh
    snakemake --snakefile workflow/Snakefile --use-conda --cores 8 \
              --configfile examples/strain_competition/config.yaml

    # 2. run the three modes on the refined GEMs
    python examples/strain_competition/genome_scale.py --models examples/strain_competition/results/refined

Three modes
-----------
  1. RESOURCE.   Bare metabolic competition on glucose: the emergent winner (lowest R*).
  2. INTERFERENCE.  A colicin (``Bacteriocin`` layer) from the producer strain suppresses
     the resource winner -- interference can overturn exploitation.
  3. NICHE.      Add L-arabinose.  IF a strain has a private catabolic route to it, it
     escapes exclusion and coexists.

Honest caveat (see README "What is / isn't modelled")
-----------------------------------------------------
On the toy, the traits are ASSIGNED to isolate each mode.  On real strains they come from
the GEMs, and the clean separation may NOT reproduce: closely related E. coli reconstruct
to near-identical GEMs (so resource competition can be degenerate), most E. coli catabolise
arabinose (so it is not a private niche), and the colicin operon is a specific accessory
gene -- here the PRODUCER is assigned as a scenario knob, not detected.  What this script
does is REPORT what the real GEMs actually do, which is the point of running it.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

EFFICIENT = "Ecoli__glc_specialist"    # expected best glucose competitor
COLICIN = "Ecoli__colicinogenic"       # assigned colicin producer (a scenario knob)
NICHE = "Ecoli__arabinose_user"        # expected arabinose user
STRAINS = (EFFICIENT, COLICIN, NICHE)


def load_models(models_dir: Path):
    from muode.organism import CobraOrganism

    found = {}
    for path in sorted(models_dir.glob("*.xml")) + sorted(models_dir.glob("*.xml.gz")):
        code = path.name.split(".")[0]
        if code in STRAINS:
            found[code] = CobraOrganism.from_file(str(path), id=code)
    return found


def _diet(base_csv: Path, with_arabinose: bool):
    from muode.diet import Diet

    d = Diet.from_csv(base_csv)
    if not with_arabinose:
        return d
    conc = dict(d.concentrations)
    conc["arab__L_e"] = 10.0                       # a second, potentially-private sugar
    influx = dict(d.influx); influx["arab__L_e"] = 1.0
    return Diet(conc, influx, dict(d.max_uptake), dict(d.source), name="defined+arabinose")


def _ecology(with_colicin: bool):
    from muode.antagonism import Bacteriocin
    from muode.ecology import EcologyModel

    if not with_colicin:
        return EcologyModel([])
    # narrow-spectrum colicin: producer immune, targets the resource winner
    return EcologyModel([Bacteriocin(name="colicin", producers={COLICIN},
                                     targets={EFFICIENT}, metabolite="colicin_e",
                                     production=0.9, decay=0.15, ki=0.04)])


def scenario(models, base_csv, with_colicin: bool, with_arabinose: bool,
             t_end: float = 48.0, dt: float = 0.05):
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.kinetics import KineticParameters

    organisms = [models[c] for c in STRAINS if c in models]
    community = Community(organisms, {c: 1.0 / len(organisms) for c in models},
                          total_biomass=0.03)
    return DynamicFBA(t_end=t_end, dt=dt).run(
        community, _diet(base_csv, with_arabinose), KineticParameters(),
        ecology=_ecology(with_colicin))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, required=True,
                    help="directory of refined GEMs (<strain>.xml[.gz])")
    ap.add_argument("--diet", type=Path,
                    default=Path(__file__).resolve().parent / "defined_medium.csv")
    ap.add_argument("--outdir", type=Path, default=Path("results/strain_competition_gs"))
    ap.add_argument("--t-end", type=float, default=48.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    models = load_models(args.models)
    missing = [c for c in STRAINS if c not in models]
    if missing:
        print(f"ABORT: missing GEM(s) for {missing} in {args.models}", file=sys.stderr)
        return 1

    modes = [("1. resource", False, False), ("2. + colicin", True, False),
             ("3. + arabinose niche", True, True)]
    print(f"{'mode':22s}{'specialist':>12s}{'colicin':>9s}{'niche':>9s}")
    print("-" * 52)
    for label, col, ara in modes:
        res = scenario(models, args.diet, col, ara, t_end=args.t_end)
        fb = res.final_biomass()
        print(f"{label:22s}{fb.get(EFFICIENT,0):12.3f}{fb.get(COLICIN,0):9.3f}"
              f"{fb.get(NICHE,0):9.3f}")
        res.to_csv(args.outdir / label.split(".")[0].strip())
    print(f"\nwrote time-courses to {args.outdir}/  (read against README's honest caveat)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
