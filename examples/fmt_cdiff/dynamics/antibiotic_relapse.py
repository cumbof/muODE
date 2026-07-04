#!/usr/bin/env python
"""
antibiotic_relapse.py
=============================================================================
Scenario 1 -- antibiotic-induced community collapse and spore-driven
Clostridioides difficile relapse.

Runs the RELAPSE arm: the pre-intervention rCDI microbiome is exposed to a
10-day vancomycin course and then left to recover on its own (NO fecal
transplant). Every trajectory below is produced by the muODE dynamic-FBA engine
and its composable ecology layer -- nothing is scripted.

Three panels:
  A  Pharmacokinetics + commensal clearance. The single-compartment vancomycin
     concentration (right axis, first-order decay) drives a saturating death
     rate on the vegetative biomass of the susceptible obligate anaerobes
     (Ruminococcaceae / Lachnospiraceae-like Firmicutes), which crash (left
     axis). The dosing window is shaded.
  B  Bile-acid pool. Primary bile acids (taurocholate + cholate, solid) versus
     secondary bile acids (deoxycholate + lithocholate, dashed). Loss of the
     commensal 7alpha-dehydroxylation (bai) capacity inverts the primary-to-
     secondary ratio.
  C  Pathogen life cycle. C. difficile vegetative (solid) and dormant spore
     (dashed) biomass. Spores are insulated from the drug and persist through
     the course; once the drug clears in a primary-bile-rich, secondary-bile-
     poor lumen they germinate and drive a secondary bloom -- the relapse.

Outputs (written to ./results/):
  antibiotic_relapse.png            -- the three-panel figure
  antibiotic_relapse_biomass.csv    -- full biomass trajectories
  antibiotic_relapse_metabolites.csv-- full extracellular-pool trajectories
  antibiotic_relapse_pathogen.csv   -- C. difficile vegetative + spore + drug

Run:
  conda run -n muode python antibiotic_relapse.py
"""

from __future__ import annotations

import warnings

import numpy as np

import common as C


def run():
    """Execute the relapse scenario (10-day vancomycin, no FMT)."""
    # commensal_level > 0  -> the functional pre-antibiotic patient microbiome.
    # inject_members = None -> no transplant: the community must recover alone.
    return C.run_scenario(commensal_level=C.COMMENSAL_LEVEL, inject_members=None,
                          with_antibiotic=True)


def make_figure(result, antibiotic):
    """Render the three-panel relapse figure from a completed simulation."""
    plt = C._mpl()
    t = result.biomass.index.values                    # time axis (days)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.2))

    # ---- Panel A: drug PK + susceptible-commensal crash -------------------
    # NOTE ON THE DUAL AXIS: mixing two y-scales is normally discouraged, but a
    # drug concentration overlaid on its pharmacodynamic effect is the standard
    # PK/PD convention chosen for this panel. The
    # axes are colour-keyed and labelled to keep the reading unambiguous.
    C.shade_dosing(axA)                                # shade the 10-day course
    susceptible = C.biomass_sum(result, C.SUSCEPTIBLE_COMMENSALS)
    axA.plot(t, susceptible, color=C.PALETTE["commensal"], lw=2,
             label="Susceptible obligate anaerobes")
    axA.set_xlabel("time (days)")
    axA.set_ylabel("commensal biomass (gDW/L)", color=C.PALETTE["commensal"])
    axA.tick_params(axis="y", labelcolor=C.PALETTE["commensal"])
    axA.set_title("A  Vancomycin PK/PD & commensal crash", loc="left")

    axA2 = axA.twinx()
    axA2.grid(False)
    axA2.plot(t, C.drug_curve(antibiotic, t), color=C.PALETTE["drug"], lw=1.8,
              ls="--", label="Luminal drug")
    axA2.set_ylabel("drug concentration (a.u.)", color=C.PALETTE["drug"])
    axA2.tick_params(axis="y", labelcolor=C.PALETTE["drug"])
    # one combined legend for the twinned axes
    lines = axA.get_lines() + axA2.get_lines()
    axA.legend(lines, [ln.get_label() for ln in lines], loc="upper right", fontsize=8)

    # ---- Panel B: bile-acid pool inversion --------------------------------
    C.shade_dosing(axB)
    axB.plot(t, C.primary_bile(result), color=C.PALETTE["primary"], lw=2,
             ls="-", label="Primary (taurocholate + cholate)")
    axB.plot(t, C.secondary_bile(result), color=C.PALETTE["secondary"], lw=2,
             ls="--", label="Secondary (deoxycholate + lithocholate)")
    axB.set_xlabel("time (days)")
    axB.set_ylabel("bile-acid concentration (mmol/L)")
    axB.set_title("B  Bile-acid pool: shield lost", loc="left")
    axB.legend(loc="upper left", fontsize=8)

    # ---- Panel C: pathogen vegetative vs spore ----------------------------
    C.shade_dosing(axC)
    veg = result.biomass[C.CDIFF]
    spore = result.spores[C.CDIFF] if result.spores is not None else veg * 0.0
    axC.plot(t, veg, color=C.PALETTE["pathogen"], lw=2, ls="-",
             label="C. difficile vegetative")
    axC.plot(t, spore, color=C.PALETTE["pathogen"], lw=2, ls="--",
             label="C. difficile spores")
    axC.set_xlabel("time (days)")
    axC.set_ylabel("biomass (gDW/L)")
    axC.set_title("C  Spore persistence -> relapse", loc="left")
    axC.legend(loc="upper left", fontsize=8)

    fig.suptitle("Antibiotic-induced collapse and spore-driven "
                 "C. difficile relapse (no FMT)", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def save_tables(result, antibiotic, outdir):
    """Persist the underlying trajectories as CSVs for full reproducibility."""
    import pandas as pd

    result.biomass.to_csv(outdir / "antibiotic_relapse_biomass.csv", index_label="day")
    result.metabolites.to_csv(outdir / "antibiotic_relapse_metabolites.csv", index_label="day")
    t = result.biomass.index.values
    pathogen = pd.DataFrame({
        "vegetative": result.biomass[C.CDIFF].values,
        "spores": (result.spores[C.CDIFF].values
                   if result.spores is not None else np.zeros_like(t)),
        "drug": C.drug_curve(antibiotic, t),
    }, index=t)
    pathogen.to_csv(outdir / "antibiotic_relapse_pathogen.csv", index_label="day")


def main():
    warnings.filterwarnings("ignore")     # silence benign depletion warnings
    outdir = C.output_dir()
    print("[relapse] running the relapse scenario (10-day vancomycin, no FMT) ...")
    result, antibiotic = run()

    fig = make_figure(result, antibiotic)
    png = outdir / "antibiotic_relapse.png"
    fig.savefig(png, bbox_inches="tight")
    save_tables(result, antibiotic, outdir)

    # concise, quantitative summary of the emergent outcome
    veg = result.biomass[C.CDIFF]
    print(f"  C. difficile vegetative biomass: start={veg.iloc[0]:.3g}, "
          f"end={veg.iloc[-1]:.3g} gDW/L  (relapse => end >> start)")
    print(f"  secondary bile acids at end: {C.secondary_bile(result).iloc[-1]:.3g} "
          "mmol/L  (shield NOT restored)")
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
