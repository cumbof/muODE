#!/usr/bin/env python3
"""Genome-scale multi-kingdom gut community: real GEMs decide the cross-kingdom outcome.

The same two mechanisms as mechanistic_demo.py (the toy), but on REAL reconstructions so
the result rests on stoichiometry, not toy yields:
  * B. thetaiotaomicron + K. pneumoniae -- CarveMe (BiGG), committed under models/;
  * S. cerevisiae -- curated Yeast8 (SysBioChalmers), exchanges harmonized to BiGG so it
    shares the o2_e / glc__D_e pool (CarveMe cannot build a eukaryote). See
    harmonize_fungal_gem.py.
  * the phage is a PhageInfection layer, not a GEM.

Two emergent mechanisms, each isolated by an ablation:
  1. the facultative FUNGUS scavenges mucosal O2 -> draws it down -> protects the
     obligate-anaerobe keystone (remove the fungus and O2 stays high, anaerobe suppressed);
  2. the lytic PHAGE crashes the Klebsiella pathobiont bloom (remove it and Klebsiella booms).

Workstation job (3 genome-scale LPs x 720 steps x 3 arms). Run on ts-02.
"""
import warnings, logging, sys
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path

import cobra
from muode import (Community, Diet, Domain, DynamicFBA, EcologyModel, KineticParameters,
                   MicrobeTraits, OxygenSensitivity, OxygenTolerance, PhageInfection)
from muode.organism import CobraOrganism

HERE = Path(__file__).resolve().parent
BACTEROIDES = "B_thetaiotaomicron"
YEAST = "S_cerevisiae"
KLEBSIELLA = "K_pneumoniae"

MODELS = {
    BACTEROIDES: HERE / "models" / "B_thetaiotaomicron.xml.gz",
    KLEBSIELLA:  HERE / "models" / "K_pneumoniae.xml.gz",
    YEAST:       HERE / "data" / "eukaryote_models" / "S_cerevisiae_bigg.xml.gz",
}
TRAITS = {
    BACTEROIDES: MicrobeTraits(domain=Domain.BACTERIA, oxygen=OxygenTolerance.OBLIGATE_ANAEROBE),
    YEAST:       MicrobeTraits(domain=Domain.EUKARYOTE, oxygen=OxygenTolerance.FACULTATIVE,
                               cell_mass_pg=40.0, notes="S. cerevisiae; curated Yeast8"),
    KLEBSIELLA:  MicrobeTraits(domain=Domain.BACTERIA, oxygen=OxygenTolerance.FACULTATIVE),
}
ABUND = {BACTEROIDES: 0.45, YEAST: 0.15, KLEBSIELLA: 0.40}


def _org(code):
    return CobraOrganism(cobra.io.read_sbml_model(str(MODELS[code])), id=code)


def build_community(members):
    return Community([_org(c) for c in members],
                     {c: ABUND[c] for c in members}, total_biomass=0.03,
                     traits={c: TRAITS[c] for c in members})


def diet():
    return Diet.from_csv(HERE / "mucosal_aerobic.csv", name="mucosal_aerobic")


def kinetics():
    return KineticParameters(metabolite_defaults={"glc__D_e": (10.0, 0.5)})


# O2 handling, GROUNDED and no-double-count (see o2_ground.py sweep):
#   ki_o2    -- B. theta is a moderately aerotolerant anaerobe (cytochrome bd oxidase,
#               nanaerobic respiration); inhibited at trace O2, anchored at ~5 uM.
#   consume=False -- the real GEMs already exchange o2_e via FBA, so O2 draw-down is left
#               entirely to the fungus's actual stoichiometric respiration. Adding the
#               layer's phenomenological consumption on top would DOUBLE-COUNT (the
#               oxygen.py docstring prescribes consume=False for real GEMs). The sweep
#               confirms it: consume=True gives an inflated 1.73x protection; the honest
#               consume=False gives 1.25x -- modest, real, and apt for B. theta's O2
#               tolerance. Result is influx-independent (0.15..4.0 identical under the
#               layer's scavenging), so o2_influx is just the mucosal-diffusion baseline.
O2_KI = 0.005
O2_INFLUX = 1.0


def ecology(members, with_phage, ki_o2=O2_KI, o2_influx=O2_INFLUX, consumption=18.0,
            consume=False):
    traits = {c: TRAITS[c] for c in members}
    layers = [OxygenSensitivity.from_traits(traits, o2_initial=0.3, o2_influx=o2_influx,
                                            ki_o2=ki_o2, consumption=consumption,
                                            consume=consume)]
    if with_phage:
        layers.append(PhageInfection(host=KLEBSIELLA, name="vB_Kpn", adsorption_rate=11.0,
                                     burst_size=60.0, latent_period=0.4, decay_rate=0.1,
                                     initial_titer=0.5))
    return EcologyModel(layers)


def run(with_fungus, with_phage, **eco_kw):
    members = [c for c in (BACTEROIDES, YEAST, KLEBSIELLA) if with_fungus or c != YEAST]
    return DynamicFBA(t_end=36.0, dt=0.05).run(
        build_community(members), diet(), kinetics(),
        ecology=ecology(members, with_phage, **eco_kw))


def main():
    print(f"{'scenario':36s}  {'Bacteroides':>11s} {'yeast':>8s} "
          f"{'Klebsiella':>10s} {'O2_final':>8s}", flush=True)
    print("-" * 80, flush=True)
    for label, fungus, phage in [
            ("full community (fungus + phage)", True, True),
            ("no fungus (O2 not scavenged)", False, True),
            ("no phage (Klebsiella unchecked)", True, False)]:
        res = run(fungus, phage)
        fb = res.final_biomass()
        o2 = res.environment["oxygen"].iloc[-1] if (
            res.environment is not None and "oxygen" in res.environment) else float("nan")
        print(f"{label:36s}  {fb.get(BACTEROIDES, 0.0):11.3f} "
              f"{fb.get(YEAST, float('nan')):8.3f} {fb.get(KLEBSIELLA, 0.0):10.3f} "
              f"{o2:8.3f}", flush=True)
    print("\nReading: fungus scavenges O2 -> protects the anaerobe; phage checks the pathobiont.")


if __name__ == "__main__":
    main()
