#!/usr/bin/env python3
"""Genome-scale phage therapy: a real host + a real commensal, one lytic phage.

The genome-scale counterpart of ``phage_scenario.py`` (which is the toy).  A facultative
pathobiont (*Klebsiella pneumoniae*) out-blooms a commensal (*E. coli*) for the shared
carbon; a lytic phage specific to the pathobiont crashes the bloom and RELEASES the
commensal.  Unlike the toy, the yields are the CarveMe reconstructions' stoichiometry --
who out-blooms whom is decided by the genomes, not by a dial -- and the only difference
between the two arms is whether the phage is present.

    # 1. genomes + GEMs (workstation)
    bash examples/phage_therapy/download_genomes.sh
    snakemake --snakefile workflow/Snakefile --use-conda --cores 8 \
              --configfile examples/phage_therapy/config.yaml

    # 2. run the therapy vs no-therapy contrast on the refined GEMs
    python examples/phage_therapy/genome_scale.py --models examples/phage_therapy/results/refined

The phage
---------
The phage is not a GEM -- it is a Levin-Stewart ``PhageInfection`` layer.  Its parameters
are literature-scale for a lytic Klebsiella phage (burst ~80 virions, latent ~0.5 h),
NOT a specific measured isolate: they set the tempo of the crash, and the falsifiable
claim ("phage present -> host down, commensal up") is robust to their exact values.
Treat them as ASSUMED, and cite a specific phage's kinetics if you pin this to one.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HOST = "K_pneumoniae"       # the pathobiont the phage targets
COMMENSAL = "E_coli"        # competes for carbon; not a host for this phage

#: PhageInfection parameters -- literature-scale for a lytic Klebsiella phage, ASSUMED.
PHAGE = dict(name="vB_Kpn-like", adsorption_rate=10.0, burst_size=80.0,
             latent_period=0.5, decay_rate=0.1, initial_titer=0.5)


def load_models(models_dir: Path):
    from muode.organism import CobraOrganism

    found = {}
    for path in sorted(models_dir.glob("*.xml")) + sorted(models_dir.glob("*.xml.gz")):
        code = path.name.split(".")[0]
        if code in (HOST, COMMENSAL):
            found[code] = CobraOrganism.from_file(str(path), id=code)
    return found


def scenario(models, therapy: bool, diet, t_end: float = 48.0, dt: float = 0.05):
    """One arm: host + commensal on ``diet``, with the phage present iff ``therapy``."""
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.ecology import EcologyModel
    from muode.kinetics import KineticParameters
    from muode.phage import PhageInfection

    organisms = [models[c] for c in (HOST, COMMENSAL) if c in models]
    community = Community(organisms, {c: 1.0 / len(organisms) for c in models},
                          total_biomass=0.02)
    layers = [PhageInfection(host=HOST, **PHAGE)] if therapy else []
    return DynamicFBA(t_end=t_end, dt=dt).run(
        community, diet, KineticParameters(), ecology=EcologyModel(layers))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, required=True,
                    help="directory of refined GEMs (<CODE>.xml[.gz])")
    ap.add_argument("--diet", type=Path,
                    default=Path(__file__).resolve().parent / "gut_glucose.csv")
    ap.add_argument("--outdir", type=Path, default=Path("results/phage_therapy"))
    ap.add_argument("--t-end", type=float, default=48.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    from muode.diet import Diet

    models = load_models(args.models)
    missing = [c for c in (HOST, COMMENSAL) if c not in models]
    if missing:
        print(f"ABORT: missing GEM(s) for {missing} in {args.models}", file=sys.stderr)
        return 1
    diet = Diet.from_csv(args.diet)

    no = scenario(models, therapy=False, diet=diet, t_end=args.t_end)
    yes = scenario(models, therapy=True, diet=diet, t_end=args.t_end)
    nb, yb = no.final_biomass(), yes.final_biomass()

    print(f"{'':12s}{HOST:>14s}{COMMENSAL:>10s}")
    print(f"{'no phage':12s}{nb.get(HOST, 0):14.4f}{nb.get(COMMENSAL, 0):10.4f}")
    print(f"{'+ phage':12s}{yb.get(HOST, 0):14.4f}{yb.get(COMMENSAL, 0):10.4f}")
    host_knockdown = yb.get(HOST, 0) < nb.get(HOST, 0)
    commensal_released = yb.get(COMMENSAL, 0) > nb.get(COMMENSAL, 0)
    print(f"\nphage knocks the host down: {host_knockdown}; releases the commensal: "
          f"{commensal_released}")
    no.to_csv(args.outdir / "no_phage")
    yes.to_csv(args.outdir / "with_phage")
    print(f"wrote time-courses to {args.outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
