#!/usr/bin/env python
"""Strain competition: one species, three strains, one shared environment.

Runs anywhere (numpy/scipy only). It asks a single ecological question — *which
strain of the same species survives when they share a habitat and a carbon
source?* — and shows that the answer is not fixed: it depends on which
interaction dominates. All three strains are the *same species* (toy
`Escherichia coli`-like models); they differ only in a handful of accessory
traits, exactly the kind of strain-level variation that distinguishes clones of
one species.

Three scenarios, each an emergent outcome rather than a scripted one:

  1. RESOURCE COMPETITION.  On a single limiting sugar the strains obey Gause's
     law / Tilman's R*: the strain that grows best on glucose (highest yield +
     affinity, i.e. the lowest break-even concentration R*) competitively
     EXCLUDES the others. One winner.

  2. INTERFERENCE COMPETITION.  Add a colicin: a *metabolically weaker* producer
     strain secretes a toxin that suppresses the otherwise-dominant strain
     (muode.Bacteriocin). The producer now wins DESPITE growing slower — chemical
     warfare overturns the resource-competition winner.

  3. NICHE DIFFERENTIATION.  Give one strain a private accessory catabolic
     pathway (a sugar the others cannot touch). It persists on its own resource
     while the rivals fight over glucose -> COEXISTENCE instead of exclusion.

Everything is a dependency-light LinprogOrganism so the *mechanisms* — not a
particular genome-scale model — are exercised. The genome-scale version would
swap these for real strain GEMs reconstructed from strain-resolved genomes
(CarveMe per strain; see README.md), and the colicin/immunity and accessory
sugar operon would come from each strain's accessory genome.

Run with:  PYTHONPATH=<repo> python examples/strain_competition/mechanistic_demo.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from muode import (
    Bacteriocin,
    Community,
    Diet,
    DynamicFBA,
    EcologyModel,
    KineticParameters,
    LinprogOrganism,
)

# Three strains of ONE species. Ids keep the muODE "<species>__<strain>"-style
# convention so they read as strains, not distinct species.
EFFICIENT = "Ecoli__glc_specialist"   # best glucose competitor (lowest R*)
COLICIN = "Ecoli__colicinogenic"      # weaker grower, but secretes a colicin
NICHE = "Ecoli__arabinose_user"       # mediocre on glucose, has a private sugar

GLC = "glc_e"    # the shared, limiting carbon source
ARA = "ara_e"    # a private accessory substrate only NICHE can use
COLICIN_POOL = "colicin_e"


def _strain(id: str, substrates: dict[str, float]) -> LinprogOrganism:
    """A toy grower that can grow on one or more substrates.

    ``substrates`` maps a metabolite id to the biomass yield on it. Each
    substrate gets an uptake exchange (dynamically bounded by Michaelis-Menten
    in the engine) and a growth reaction consuming it.
    """
    reactions = []
    exchanges = {}
    objective = {}
    for met, yld in substrates.items():
        ex = f"EX_{met}"
        reactions.append((ex, {met: -1.0}, -1000.0, 1000.0))
        exchanges[met] = ex
        grow = f"GROW_{met}"
        reactions.append((grow, {met: -1.0}, 0.0, 1000.0))
        objective[grow] = yld
    return LinprogOrganism(id=id, reactions=reactions, objective=objective,
                           exchanges=exchanges)


def build_community() -> Community:
    organisms = [
        _strain(EFFICIENT, {GLC: 0.42}),                 # highest glucose yield
        _strain(COLICIN, {GLC: 0.30}),                   # lower glucose yield
        _strain(NICHE, {GLC: 0.26, ARA: 0.38}),          # weak on glc, strong on ara
    ]
    # Equal starting inocula: nobody is handed the win up front.
    abundances = {EFFICIENT: 1 / 3, COLICIN: 1 / 3, NICHE: 1 / 3}
    return Community(organisms, abundances, total_biomass=0.03)


def kinetics() -> KineticParameters:
    kp = KineticParameters(metabolite_defaults={GLC: (10.0, 0.5), ARA: (10.0, 0.5)})
    # Reinforce the R* ordering on glucose with uptake affinity (Km): the
    # specialist scavenges glucose to a lower residual concentration.
    kp.set(EFFICIENT, GLC, vmax=10.0, km=0.20)
    kp.set(COLICIN, GLC, vmax=10.0, km=0.80)
    kp.set(NICHE, GLC, vmax=10.0, km=0.80)
    return kp


def diet(with_arabinose: bool) -> Diet:
    conc = {GLC: 10.0}
    influx = {GLC: 1.0}
    if with_arabinose:
        conc[ARA] = 5.0
        influx[ARA] = 0.5
    return Diet(concentrations=conc, influx=influx, name="shared_habitat")


def ecology(with_colicin: bool) -> EcologyModel:
    if not with_colicin:
        return EcologyModel([])
    # A narrow-spectrum colicin: the producer is immune, and here it targets the
    # specific rival that would otherwise dominate (the glucose specialist).
    return EcologyModel([
        Bacteriocin(name="colicin", producers={COLICIN}, targets={EFFICIENT},
                    metabolite=COLICIN_POOL, production=0.9, decay=0.15, ki=0.04),
    ])


def run(with_colicin: bool, with_arabinose: bool):
    comm = build_community()
    return DynamicFBA(t_end=48.0, dt=0.01).run(
        comm, diet(with_arabinose), kinetics(), ecology=ecology(with_colicin))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=Path("results/strain_competition"))
    args = ap.parse_args()

    scenarios = [
        ("1. resource competition (glucose only)", False, False),
        ("2. + colicin (interference)",            True,  False),
        ("3. + private niche (glucose+arabinose)", True,  True),
    ]
    print(f"{'scenario':42s}  {'specialist':>10s} {'colicin':>8s} {'niche':>8s}")
    print("-" * 74)
    results = {}
    for label, colicin, arabinose in scenarios:
        res = run(colicin, arabinose)
        results[label] = res
        fb = res.final_biomass()
        print(f"{label:42s}  {fb.get(EFFICIENT, 0.0):10.3f} "
              f"{fb.get(COLICIN, 0.0):8.3f} {fb.get(NICHE, 0.0):8.3f}")

    print("\nReading (all start at equal biomass on the same glucose pool):")
    print("  1. On one limiting sugar the best glucose competitor excludes the")
    print("     others (competitive exclusion / lowest R* wins).")
    print("  2. The colicin producer suppresses the specialist and takes over,")
    print("     though it is the weaker grower -- interference beats exploitation.")
    print("  3. The niche strain lives on its private sugar and coexists; niche")
    print("     differentiation is the escape from competitive exclusion.")

    for label, res in results.items():
        tag = label.split(".", 1)[0].strip()
        res.to_csv(args.outdir / f"scenario_{tag}")
    print(f"\nWrote time-courses to {args.outdir}/")


if __name__ == "__main__":
    main()
