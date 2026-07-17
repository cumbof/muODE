#!/usr/bin/env python
"""Multi-kingdom gut snippet: bacteria + a fungus + a phage, on toy models.

Runs anywhere (numpy/scipy only) and shows the two mechanisms that make a sample
*multi-kingdom* rather than bacteria-only behave differently, both as emergent
outcomes rather than scripted results:

  1. EUKARYOTE / oxygen.  A facultative gut yeast (Saccharomyces) scavenges the
     oxygen diffusing in from the mucosa.  That draws O2 down and *protects* the
     obligate-anaerobe keystone (Bacteroides-like), which is otherwise poisoned.
     -> compare the community with vs without the fungus.

  2. VIRUS / phage.  A lytic phage specific to a facultative pathobiont
     (Klebsiella-like) crashes its bloom and releases the carbon back to the
     commensals.
     -> compare the community with vs without the phage.

Everything is a dependency-light LinprogOrganism so the *mechanisms* — not a
particular genome-scale model — are what is exercised.  The genome-scale version
swaps these for real GEMs: Bacteroides/Klebsiella via CarveMe, and the fungus as
curated Yeast8 (S. cerevisiae) with its exchanges harmonized to BiGG so it shares
the pool — CarveMe cannot build a eukaryote.  The phage is not a GEM at all but a
PhageInfection layer.  See README.md.

Run with:  PYTHONPATH=<repo> python examples/multikingdom/mechanistic_demo.py
"""

from __future__ import annotations

from muode import (
    Community,
    Diet,
    Domain,
    DynamicFBA,
    EcologyModel,
    KineticParameters,
    LinprogOrganism,
    MicrobeTraits,
    OxygenSensitivity,
    OxygenTolerance,
    PhageInfection,
)

BACTEROIDES = "B_thetaiotaomicron"   # obligate-anaerobe keystone (bacteria)
YEAST = "S_cerevisiae"             # facultative oxygen scavenger (eukaryote); Yeast8
KLEBSIELLA = "K_pneumoniae"          # facultative pathobiont, the phage host (bacteria)


def _grower(id: str, yld: float, product: str | None = None) -> LinprogOrganism:
    stoich = {"glc_e": -1.0}
    reactions = [("EX_glc_e", {"glc_e": -1.0}, -1000.0, 1000.0)]
    exchanges = {"glc_e": "EX_glc_e"}
    if product is not None:
        stoich[product] = 1.0
        reactions.append((f"EX_{product}", {product: -1.0}, 0.0, 1000.0))
        exchanges[product] = f"EX_{product}"
    reactions.append(("GROW", stoich, 0.0, 1000.0))
    return LinprogOrganism(id=id, reactions=reactions, objective={"GROW": yld},
                           exchanges=exchanges)


def build_community() -> Community:
    organisms = [
        _grower(BACTEROIDES, yld=0.18, product="ac_e"),   # slow, ferments to acetate
        _grower(YEAST, yld=0.12),                        # facultative, modest grower
        _grower(KLEBSIELLA, yld=0.32),                     # fast pathobiont bloom
    ]
    traits = {
        BACTEROIDES: MicrobeTraits(domain=Domain.BACTERIA,
                                   oxygen=OxygenTolerance.OBLIGATE_ANAEROBE),
        YEAST: MicrobeTraits(domain=Domain.EUKARYOTE,
                               oxygen=OxygenTolerance.FACULTATIVE,
                               cell_mass_pg=40.0, notes="S. cerevisiae; curated Yeast8"),
        KLEBSIELLA: MicrobeTraits(domain=Domain.BACTERIA,
                                  oxygen=OxygenTolerance.FACULTATIVE),
    }
    return Community(organisms,
                     {BACTEROIDES: 0.45, YEAST: 0.15, KLEBSIELLA: 0.40},
                     total_biomass=0.03, traits=traits)


def diet() -> Diet:
    # carbon plus a steady trickle of oxygen from the mucosa
    return Diet(concentrations={"glc_e": 8.0, "o2_e": 0.3},
                influx={"glc_e": 1.2, "o2_e": 0.15}, name="mucosal_gut")


def kinetics() -> KineticParameters:
    return KineticParameters(metabolite_defaults={"glc_e": (10.0, 0.5)})


def ecology(community: Community, with_fungus: bool, with_phage: bool) -> EcologyModel:
    layers = []
    # oxygen layer built straight from the declared traits; if we ablate the
    # fungus its scavenging simply disappears and O2 stays high.
    traits = dict(community.traits)
    if not with_fungus:
        traits.pop(YEAST, None)
    layers.append(OxygenSensitivity.from_traits(
        traits, o2_initial=0.3, o2_influx=0.15, consumption=18.0))
    if with_phage:
        layers.append(PhageInfection(host=KLEBSIELLA, name="vB_Kpn",
                                     adsorption_rate=11.0, burst_size=60.0,
                                     latent_period=0.4, decay_rate=0.1,
                                     initial_titer=0.5))
    return EcologyModel(layers)


def run(with_fungus: bool, with_phage: bool):
    comm = build_community()
    if not with_fungus:
        comm = Community([o for o in comm.organisms if o.id != YEAST],
                         {k: v for k, v in comm.abundances.items() if k != YEAST},
                         total_biomass=comm.total_biomass,
                         traits={k: v for k, v in comm.traits.items() if k != YEAST})
    eco = ecology(build_community(), with_fungus, with_phage)
    return DynamicFBA(t_end=36.0, dt=0.02).run(comm, diet(), kinetics(), ecology=eco)


def main() -> None:
    scenarios = [
        ("full community (fungus + phage)", True, True),
        ("no fungus (O2 not scavenged)", False, True),
        ("no phage (Klebsiella unchecked)", True, False),
    ]
    print(f"{'scenario':36s}  {'Bacteroides':>11s} {'yeast':>8s} "
          f"{'Klebsiella':>10s} {'O2_final':>8s}")
    print("-" * 80)
    for label, fungus, phage in scenarios:
        res = run(fungus, phage)
        fb = res.final_biomass()
        o2 = res.environment["oxygen"].iloc[-1]
        bact = fb.get(BACTEROIDES, 0.0)
        cand = fb.get(YEAST, float("nan"))
        kleb = fb.get(KLEBSIELLA, 0.0)
        print(f"{label:36s}  {bact:11.3f} {cand:8.3f} {kleb:10.3f} {o2:8.3f}")

    print("\nReading:")
    print("  * Removing the fungus leaves O2 high -> the obligate-anaerobe")
    print("    Bacteroides is suppressed (cross-kingdom protection lost).")
    print("  * Removing the phage lets the Klebsiella pathobiont bloom unchecked.")


if __name__ == "__main__":
    main()
