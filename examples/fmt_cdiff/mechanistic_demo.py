#!/usr/bin/env python
"""Mechanistic recurrent-CDI / FMT demo -- runs anywhere (no GEMs needed).

This drives the full muODE ecology stack (antibiotic pharmacokinetics, spore
survival, bile-acid transformation + germination gating, SCFA/pH inhibition and
nutrient competition) on dependency-light toy models, so the *mechanisms* behind
FMT for recurrent C. difficile can be seen on any machine -- including the ones
where CarveMe cannot run.  The genome-scale version is the rest of this folder
(see README.md).

It runs two arms that differ only by the transplant:

  * RECURRENCE -- antibiotics clear vegetative C. difficile, but the spore
    reservoir survives and germinates in the wiped gut -> the infection returns.
  * FMT        -- a donor community is transplanted after the antibiotic; it
    competes for carbon, acidifies via SCFA, and restores 7α-dehydroxylation
    (cholate -> deoxycholate), which blocks germination and inhibits the
    pathogen -> the infection clears.

Usage:
    python examples/fmt_cdiff/mechanistic_demo.py [--outdir DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from scenario import BAI, CDIFF, COMPETITOR, cdi_scenario  # noqa: E402


def _summary(tag: str, result) -> None:
    cd = result.biomass[CDIFF]
    sp = result.spores[CDIFF]
    env = result.environment
    print(f"\n=== {tag} ===")
    print(f"  C. difficile vegetative : start {cd.iloc[0]:.3f}  "
          f"nadir {cd.min():.3f}  final {cd.iloc[-1]:.3f} gDW/L")
    print(f"  C. difficile spores     : start {sp.iloc[0]:.3f}  final {sp.iloc[-1]:.3f} gDW/L")
    print(f"  donor competitor / bai  : {result.biomass[COMPETITOR].iloc[-1]:.3f} / "
          f"{result.biomass[BAI].iloc[-1]:.3f} gDW/L")
    print(f"  final pH                : {env['pH'].iloc[-1]:.2f}")
    print(f"  final deoxycholate      : {result.metabolites['dca_e'].iloc[-1]:.2f} mmol/L")
    print(f"  final germination signal: {env['germination_signal'].iloc[-1]:.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=Path("results/fmt_mechanistic"))
    ap.add_argument("--t-end", type=float, default=96.0)
    args = ap.parse_args()

    recurrence = cdi_scenario(fmt=False, t_end=args.t_end)
    fmt = cdi_scenario(fmt=True, t_end=args.t_end)

    _summary("RECURRENCE (antibiotic only)", recurrence)
    _summary("FMT (antibiotic + transplant)", fmt)

    cd_rec = recurrence.biomass[CDIFF].iloc[-1]
    cd_fmt = fmt.biomass[CDIFF].iloc[-1]
    print(f"\nC. difficile final burden:  recurrence {cd_rec:.3f}  vs  FMT {cd_fmt:.3f} gDW/L")
    print("-> the only difference between the arms is the timed donor injection.")

    for tag, res in [("recurrence", recurrence), ("fmt", fmt)]:
        outdir = args.outdir / tag
        res.to_csv(outdir)
        try:
            from muode import viz
            viz.save_all(res, outdir)
        except ImportError:
            pass
    print(f"\nWrote time-courses (and figures if matplotlib is installed) to {args.outdir}/")


if __name__ == "__main__":
    main()
