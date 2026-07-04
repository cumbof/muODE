#!/usr/bin/env python
"""
fmt_resolution.py
=============================================================================
Scenario 2 -- fecal microbiota transplant (FMT) engraftment and resolution of
recurrent C. difficile infection.

Runs the CURE arm: the identical 10-day vancomycin course as the relapse
scenario, but at Day 12 (once the luminal drug has decayed below EC50) a bolus
of the healthy donor community is injected. The ONLY difference from the relapse
arm is that timed transplant -- everything downstream (engraftment, shield
restoration, pathogen clearance) is emergent.

Three panels:
  A  Engraftment. Biomass trajectories after the Day-12 FMT: the donor commensals
     (and specifically the bai-operon guild) expand and aggressively outcompete
     the germinating C. difficile for shared mucosal carbohydrates and amino
     acids. The FMT time is marked.
  B  Bile-acid shield restoration. As donor biomass expands, cumulative
     7alpha-dehydroxylation crosses threshold and the primary-to-secondary bile
     ratio reverts toward the healthy baseline -- the reverse of the relapse arm.
  C  Pathogen clearance. C. difficile vegetative (solid) and spore (dashed)
     biomass. In stark contrast to the relapse arm, the germinating pathogen
     meets a dual constraint (Stickland-substrate starvation + secondary bile
     inhibition) and is driven toward extinction.

Outputs (./results/): fmt_resolution.png + CSV trajectories.

Run:  conda run -n muode python fmt_resolution.py
"""

from __future__ import annotations

import warnings

import numpy as np

import common as C

# The bai-operon guild whose expansion restores the shield (highlighted in A).
BAI_GUILD = [s.member_id for s in C.CONSORTIUM_SPECS if s.bai]


def run():
    """Execute the FMT scenario: 10-day vancomycin + healthy-donor bolus at Day 12."""
    return C.run_scenario(commensal_level=C.COMMENSAL_LEVEL,
                          inject_members=C.CONSORTIUM_12, with_antibiotic=True)


def make_figure(result, antibiotic):
    plt = C._mpl()
    t = result.biomass.index.values

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.2))

    def mark_fmt(ax):
        ax.axvline(C.FMT_DAY, color=C.PALETTE["muted"], lw=1.2, ls=":")
        ax.annotate("FMT", xy=(C.FMT_DAY, ax.get_ylim()[1]),
                    xytext=(C.FMT_DAY + 0.6, 0.92), textcoords=("data", "axes fraction"),
                    fontsize=8, color=C.PALETTE["muted"])

    # ---- Panel A: engraftment (donor expansion vs pathogen) ---------------
    C.shade_dosing(axA)
    donor = C.biomass_sum(result, C.CONSORTIUM_12)
    bai = C.biomass_sum(result, BAI_GUILD)
    axA.plot(t, donor, color=C.PALETTE["commensal"], lw=2, label="Donor commensals (total)")
    axA.plot(t, bai, color=C.PALETTE["secondary"], lw=2, ls="-.",
             label="bai-operon guild")
    axA.plot(t, result.biomass[C.CDIFF], color=C.PALETTE["pathogen"], lw=2,
             label="C. difficile (vegetative)")
    axA.set_xlabel("time (days)")
    axA.set_ylabel("biomass (gDW/L)")
    axA.set_title("A  Donor engraftment", loc="left")
    mark_fmt(axA)
    axA.legend(loc="upper left", fontsize=8)

    # ---- Panel B: bile-acid shield restoration ----------------------------
    C.shade_dosing(axB)
    axB.plot(t, C.primary_bile(result), color=C.PALETTE["primary"], lw=2, ls="-",
             label="Primary (taurocholate + cholate)")
    axB.plot(t, C.secondary_bile(result), color=C.PALETTE["secondary"], lw=2, ls="--",
             label="Secondary (deoxycholate + lithocholate)")
    axB.set_xlabel("time (days)")
    axB.set_ylabel("bile-acid concentration (mmol/L)")
    axB.set_title("B  Bile-acid shield restored", loc="left")
    mark_fmt(axB)
    axB.legend(loc="upper left", fontsize=8)

    # ---- Panel C: pathogen clearance --------------------------------------
    C.shade_dosing(axC)
    veg = result.biomass[C.CDIFF]
    spore = result.spores[C.CDIFF] if result.spores is not None else veg * 0.0
    axC.plot(t, veg, color=C.PALETTE["pathogen"], lw=2, ls="-",
             label="C. difficile vegetative")
    axC.plot(t, spore, color=C.PALETTE["pathogen"], lw=2, ls="--",
             label="C. difficile spores")
    axC.set_xlabel("time (days)")
    axC.set_ylabel("biomass (gDW/L)")
    axC.set_title("C  Pathogen driven to extinction", loc="left")
    mark_fmt(axC)
    axC.legend(loc="upper right", fontsize=8)

    fig.suptitle("FMT engraftment restores the bile-acid shield and "
                 "clears C. difficile", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def main():
    warnings.filterwarnings("ignore")
    outdir = C.output_dir()
    print("[fmt] running the FMT scenario (10-day vancomycin + Day-12 donor bolus) ...")
    result, antibiotic = run()

    fig = make_figure(result, antibiotic)
    png = outdir / "fmt_resolution.png"
    fig.savefig(png, bbox_inches="tight")
    result.biomass.to_csv(outdir / "fmt_resolution_biomass.csv", index_label="day")
    result.metabolites.to_csv(outdir / "fmt_resolution_metabolites.csv", index_label="day")

    veg = result.biomass[C.CDIFF]
    print(f"  C. difficile vegetative biomass: start={veg.iloc[0]:.3g}, "
          f"end={veg.iloc[-1]:.3g} gDW/L  (cure => end -> 0)")
    print(f"  secondary bile acids at end: {C.secondary_bile(result).iloc[-1]:.3g} "
          "mmol/L  (shield RESTORED)")
    print(f"  donor biomass at end: {C.biomass_sum(result, C.CONSORTIUM_12).iloc[-1]:.3g} "
          "gDW/L  (engrafted)")
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
