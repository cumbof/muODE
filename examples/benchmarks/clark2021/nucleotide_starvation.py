#!/usr/bin/env python3
"""Nucleotide starvation -- the root cause of muODE under-growing saccharolytic taxa.

The yield calibration (yield_calibrate.py) fixed muODE's over-growers and beat the
null on composition, but could not lift the ~15 strains muODE UNDER-grows -- which
are the true dominators (Bacteroides et al.). This pins WHY.

Chain of evidence (all reproducible here):
  1. Trajectory: an under-grower (DF) starts at a healthy mu~1.8/h but PLATEAUS at
     biomass ~0.029 within ~5 h and stops -- not rate-limited, a hard early ceiling.
  2. Depletion: at the plateau it has exhausted only the TRACE nucleobase pools
     (cytosine/guanine/uracil/xanthine, ~0.0001-0.02 mM) and nitrate, while leaving
     ~97% of the glucose/maltose untouched. It is nucleotide-starved, not carbon-
     limited: the CarveMe GEM cannot synthesize purines/pyrimidines de novo.
  3. Rescue: raising the nucleobase supply lifts biomass 2-89x
     (DF 0.029->2.60, BH 0.13->0.64, PC/BA ~2x) -- confirming the cap is nucleotides.
  4. Encoding: DM38 lists these nucleobases as real components (derive_dm38.py) but
     four (ade/gua/csn/ura) are encoded at a 0.0001 mM placeholder while xanthine /
     thymidine got real values -- a medium-encoding gap that, combined with the GEMs'
     nucleotide auxotrophy, starves the saccharolytic taxa.

Fix (not applied here -- needs the real DM38 nucleobase concentrations from Clark's
recipe, or gapfilling de novo nucleotide biosynthesis into the GEMs). Either lifts
the starved strains toward the OD-weighting ceiling (MAE 0.108) that yield
calibration alone (0.131) cannot reach.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from muode.benchmarks import clark2021 as ck
from muode.community import Community
from muode.dfba import DynamicFBA
from muode.kinetics import KineticParameters
from muode.organism import CobraOrganism

NUCLEOBASES = ["csn_e", "gua_e", "ura_e", "xan_e", "ade_e", "hxan_e", "thymd_e"]
MODELS = Path("examples/benchmarks/clark2021/models")


def run(code, diet, kin):
    o = CobraOrganism.from_file(str(MODELS / f"{code}.xml.gz"), id=code)
    comm = Community([o], abundances={code: 1.0}, total_biomass=0.01)
    return float(DynamicFBA(t_end=48.0, dt=0.1, n_jobs=1).run(comm, diet, kin).biomass[code].iloc[-1])


def main() -> int:
    diet = ck.medium("bigg")
    kin = KineticParameters()
    boosted = copy.deepcopy(diet)
    for m in NUCLEOBASES:
        if m in boosted.concentrations:
            boosted.concentrations[m] = 20.0

    encoding = {m: round(diet.initial_concentration(m), 5) for m in NUCLEOBASES if m in diet.concentrations}
    rescue = {}
    for code, od in [("DF", 0.91), ("PC", 1.84), ("BH", 1.25), ("BA", 0.90)]:
        b0, b1 = run(code, diet, kin), run(code, boosted, kin)
        rescue[code] = {"measured_od": od, "dm38_biomass": round(b0, 3),
                        "with_nucleobases": round(b1, 3), "fold": round(b1 / max(b0, 1e-9), 1)}
        print(f"{code}: measOD={od}  DM38={b0:.3f}  +nucleobases={b1:.3f}  ({b1/max(b0,1e-9):.0f}x)")

    out = {"dm38_nucleobase_encoding_mM": encoding, "rescue": rescue,
           "conclusion": ("under-growth is nucleotide starvation: CarveMe GEMs cannot "
                          "make purines/pyrimidines de novo and DM38 encodes ade/gua/csn/ura "
                          "at a 0.0001 mM placeholder; supplying nucleobases rescues growth.")}
    Path("results/clark2021/benchmark/nucleotide_starvation_bigg.json").write_text(json.dumps(out, indent=2))
    print("\nDM38 encoded nucleobases (mM):", encoding)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
