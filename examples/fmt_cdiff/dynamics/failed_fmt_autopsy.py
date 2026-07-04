#!/usr/bin/env python
"""
failed_fmt_autopsy.py
=============================================================================
Scenario 4 -- identifying the mechanistic bottlenecks of failed FMTs.

This scenario performs an in-silico "algorithmic autopsy" of transplants that
fail to clear the pathogen, and isolates two distinct mechanistic bottlenecks.
This script reproduces both by comparing a SUCCESSFUL FMT against each failure
mode, holding everything else fixed. Each comparison changes exactly one thing,
so the failure is attributable.

Bottleneck 1 -- Donor metabolic insufficiency (hidden cross-feeding dependency).
  A donor that carries the bai operon can still fail if it lacks the upstream
  primary polysaccharide degraders (Bacteroides spp.) that liberate the
  monomeric sugars the bai guild needs. Starved of carbon in the competitive
  lumen, the 7alpha-dehydroxylating taxa collapse before they can transform the
  bile pool -- so the secondary-bile shield is never restored and the pathogen
  persists. We compare the full donor against the same donor with the primary
  degraders removed.

Bottleneck 2 -- Bacteriophage predation of the engrafting community.
  A recipient-derived, strictly lytic phage cocktail specific to the engrafting
  bai-operon guild (the 7alpha-dehydroxylating Clostridia that build the
  secondary-bile shield) amplifies on the incoming host biomass and lyses it
  within days of the bolus. With the effector guild decimated, the shield never
  reforms and the pathogen is not cleared -- the transplant fails despite a
  metabolically complete donor. muODE models this with the Levin-Stewart
  infection ODE (muode.phage.PhageInfection) coupled to host biomass -- a class
  of dynamics that pure metabolic FBA cannot represent. We compare the full donor
  against the same donor delivered into a phage-carrying recipient.

Both scenarios inject into the depauperate Day-12 post-antibiotic lumen.

Outputs (./results/): failed_fmt_autopsy.png + CSV trajectories.

Run:  conda run -n muode python failed_fmt_autopsy.py
"""

from __future__ import annotations

import warnings

import numpy as np

from muode.phage import PhageInfection

import common as C

# The primary polysaccharide degraders whose absence is bottleneck 1.
# (common.PRIMARY_DEGRADERS)
DONOR_WITHOUT_DEGRADERS = [m for m in C.CONSORTIUM_12 if m not in C.PRIMARY_DEGRADERS]

# The engrafting bai-operon guild the recipient phage cocktail preys on
# (bottleneck 2): the 7alpha-dehydroxylating effectors that build the shield.
PHAGE_HOSTS = [s.member_id for s in C.CONSORTIUM_SPECS if s.bai]


def run_full_donor():
    """The successful reference FMT (full 12-member donor)."""
    return C.run_scenario(commensal_level=0.0, inject_members=C.CONSORTIUM_12)


def run_degrader_deficient():
    """Bottleneck 1: donor lacking the primary polysaccharide degraders."""
    return C.run_scenario(commensal_level=0.0, inject_members=DONOR_WITHOUT_DEGRADERS)


def run_phage_predation():
    """Bottleneck 2: full donor delivered into a phage-carrying recipient.

    A strictly lytic phage cocktail (one PhageInfection layer per bai effector)
    is present from t=0 as a resident virome. Each phage amplifies once its host
    is seeded by the FMT at Day 12, lysing the engrafting effector guild before
    it can rebuild the secondary-bile shield.
    """
    phages = [
        PhageInfection(
            host=host, name=f"phi-{host.split('_')[0][:4]}",
            adsorption_rate=14.0,     # (gDW/L)^-1 day^-1: aggressive predation
            burst_size=60.0,          # virions released per unit host biomass lysed
            latent_period=0.4,        # days from infection to lysis
            decay_rate=0.15,          # free-phage decay (1/day) -- persists to Day 12
            initial_titer=0.6,        # resident virome present before the FMT
            lysogeny_fraction=0.0,    # strictly lytic
        )
        for host in PHAGE_HOSTS
    ]
    return C.run_scenario(commensal_level=0.0, inject_members=C.CONSORTIUM_12,
                          extra_layers=phages)


def make_figure(full, deficient, phage_res):
    plt = C._mpl()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

    def mark_fmt(ax):
        ax.axvline(C.FMT_DAY, color=C.PALETTE["muted"], lw=1.2, ls=":")

    t = full.biomass.index.values

    # ---- Panel A: donor insufficiency -------------------------------------
    C.shade_dosing(axA)
    axA.plot(t, full.biomass[C.CDIFF], color=C.PALETTE["commensal"], lw=2,
             label="C. difficile -- full donor (cured)")
    axA.plot(t, deficient.biomass[C.CDIFF], color=C.PALETTE["pathogen"], lw=2,
             label="C. difficile -- donor lacks degraders (fails)")
    axA.plot(t, C.biomass_sum(deficient, [s.member_id for s in C.CONSORTIUM_SPECS if s.bai]),
             color=C.PALETTE["secondary"], lw=1.6, ls="-.",
             label="bai guild (degrader-deficient donor)")
    axA.set_xlabel("time (days)")
    axA.set_ylabel("biomass (gDW/L)")
    axA.set_title("A  Bottleneck 1: donor metabolic insufficiency", loc="left")
    mark_fmt(axA)
    axA.legend(loc="upper left", fontsize=8)

    # ---- Panel B: phage predation -----------------------------------------
    C.shade_dosing(axB)
    axB.plot(t, full.biomass[C.CDIFF], color=C.PALETTE["commensal"], lw=2,
             label="C. difficile -- no phage (cured)")
    axB.plot(t, phage_res.biomass[C.CDIFF], color=C.PALETTE["pathogen"], lw=2,
             label="C. difficile -- recipient phage (fails)")
    axB.plot(t, C.biomass_sum(phage_res, PHAGE_HOSTS), color=C.PALETTE["secondary"],
             lw=1.6, ls="-.", label="engrafting bai guild (lysed)")
    axB.set_xlabel("time (days)")
    axB.set_ylabel("biomass (gDW/L)")
    axB.set_title("B  Bottleneck 2: bacteriophage predation", loc="left")
    mark_fmt(axB)
    # overlay the total free-phage titre on a secondary axis (a population count,
    # a legitimately different quantity/scale from biomass -- summed over the
    # cocktail's per-host phage columns in the environment frame)
    if phage_res.environment is not None:
        titre_cols = [c for c in phage_res.environment.columns
                      if c.startswith("phage[")]
        if titre_cols:
            ax2 = axB.twinx()
            ax2.grid(False)
            ax2.plot(t, phage_res.environment[titre_cols].sum(axis=1),
                     color=C.PALETTE["drug"], lw=1.4, ls=":", label="phage titre")
            ax2.set_ylabel("free phage titre (a.u.)", color=C.PALETTE["drug"])
            ax2.tick_params(axis="y", labelcolor=C.PALETTE["drug"])
    axB.legend(loc="upper left", fontsize=8)

    fig.suptitle("Failed-FMT autopsy -- two mechanistic bottlenecks to engraftment",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def main():
    warnings.filterwarnings("ignore")
    outdir = C.output_dir()
    print("[Autopsy] running the successful reference FMT ...")
    full, _ = run_full_donor()
    print("[Autopsy] bottleneck 1: degrader-deficient donor ...")
    deficient, _ = run_degrader_deficient()
    print("[Autopsy] bottleneck 2: recipient bacteriophage ...")
    phage_res, _ = run_phage_predation()

    fig = make_figure(full, deficient, phage_res)
    png = outdir / "failed_fmt_autopsy.png"
    fig.savefig(png, bbox_inches="tight")
    for name, res in [("full", full), ("degrader_deficient", deficient),
                      ("phage", phage_res)]:
        res.biomass.to_csv(outdir / f"autopsy_{name}_biomass.csv", index_label="day")

    def cd(res):
        return res.biomass[C.CDIFF].iloc[-1]
    print(f"  C. difficile end -- full donor (cure)      : {cd(full):.3g} gDW/L")
    print(f"  C. difficile end -- degrader-deficient     : {cd(deficient):.3g} gDW/L (FAILS)")
    print(f"  C. difficile end -- recipient phage        : {cd(phage_res):.3g} gDW/L (FAILS)")
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
