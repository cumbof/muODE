#!/usr/bin/env python3
"""Hybrid-reconstruction monoculture panel -- the fast decisive test.

The reconstruction comparison (reconstruction_comparison_bigg_vs_gapseq.json)
showed CarveMe/BiGG and gapseq/ModelSEED have COMPLEMENTARY taxonomic strengths on
anaerobic DM38: gapseq grows the Bacteroidetes (the community dominators BiGG
starves via respiratory/nitrate dependence), BiGG grows several Firmicutes/Bifido.

This builds the per-strain BEST-reconstruction panel -- gapseq for the 9
Bacteroidetes (genus Prevotella/Parabacteroides/Phocaeicola/Bacteroides), BiGG for
the other 17 -- using ONLY monoculture data. For each strain it reads the
ATPM-calibrated monoculture biomass (biomass at the fitted ATP-maintenance that
matches the strain's own measured OD; under-growers stay at ATPM=0 and their
residual is reported honestly). It then asks the one question that decides whether
the 2-3h community gate re-run is worth it:

  does the hybrid panel's monoculture biomass track measured OD across the board,
  and does OD-weighting with the hybrid yields beat the BiGG-only panel toward the
  0.108 measured-OD ceiling?

No community sims here -- pure per-strain yields + the same abundance-MAE proxy
used in yield_calibrate.py --proxy, so the numbers are directly comparable.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import clark2021 as ck

OUT = Path("results/clark2021/benchmark")
BACT_GENERA = ("Prevotella", "Parabacteroides", "Phocaeicola", "Bacteroides")


def bacteroidetes() -> set[str]:
    return {c for c, s in ck.STRAINS.items() if s.species.startswith(BACT_GENERA)}


def calibrated_biomass(cal_path: Path) -> dict[str, float]:
    """Biomass at the fitted ATPM, interpolated from the calibration curve
    (identical logic to yield_calibrate.py --proxy)."""
    out: dict[str, float] = {}
    for line in cal_path.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        curve, a = r["curve"], r["atpm"]
        if len(curve) == 1:
            b = curve[0][1]
        else:
            b = curve[-1][1]
            for (a1, b1), (a2, b2) in zip(curve, curve[1:]):
                if a1 <= a <= a2:
                    b = b1 + (b2 - b1) * ((a - a1) / (a2 - a1) if a2 != a1 else 0)
                    break
        out[r["code"]] = max(0.0, b)
    return out


def main() -> int:
    bigg = calibrated_biomass(OUT / "atpm_calibration_bigg.jsonl")
    gseq = calibrated_biomass(OUT / "atpm_calibration_modelseed.jsonl")
    bact = bacteroidetes()

    # per-strain best reconstruction
    hybrid, source = {}, {}
    for c in set(bigg) | set(gseq):
        if c in bact and c in gseq:
            hybrid[c], source[c] = gseq[c], "gapseq"
        elif c in bigg:
            hybrid[c], source[c] = bigg[c], "bigg"
        elif c in gseq:
            hybrid[c], source[c] = gseq[c], "gapseq"

    obs = ck.observations(ck.load(str(OUT / "MasterDF.csv")),
                          diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg"))
    od = {o.species[0]: max(0.0, o.od) for o in obs if o.richness == 1}

    def mae(w):
        e = []
        for o in obs:
            if o.richness < 2 or not o.abundances:
                continue
            sh = [c for c in o.species if c in o.abundances and c in w]
            if not sh:
                continue
            tot = sum(w[c] for c in sh) or 1.0
            p = {c: w[c] / tot for c in sh}
            e.append(sum(abs(p[c] - o.abundances[c]) for c in sh) / len(sh))
        return round(float(np.mean(e)), 3), len(e)

    codes = sorted(c for c in hybrid if c in od)
    sp_bigg = spearmanr([bigg[c] for c in codes if c in bigg],
                        [od[c] for c in codes if c in bigg])[0]
    sp_hyb = spearmanr([hybrid[c] for c in codes], [od[c] for c in codes])[0]

    m_null, n = mae({c: 1.0 for c in od})
    rep = {
        "n_strains": len(codes),
        "panel_source": {c: source[c] for c in codes},
        "mono_biomass_vs_OD_spearman": {
            "bigg_only": round(float(sp_bigg), 3),
            "hybrid": round(float(sp_hyb), 3)},
        "abundance_MAE": {
            "uniform_null": m_null,
            "bigg_only_calibrated": mae(bigg)[0],
            "hybrid_calibrated": mae(hybrid)[0],
            "measured_OD_ceiling": mae(od)[0],
            "muODE_dFBA_reference": 0.210,
            "n_communities": n},
    }
    (OUT / "hybrid_panel_report.json").write_text(json.dumps(rep, indent=2))

    # per-strain table: does each strain now track its OD?
    print("code src     OD   calib_bio   under?")
    for c in sorted(codes, key=lambda x: -od[x]):
        u = "UNDER" if hybrid[c] < 0.5 * od[c] else ""
        print(f"{c:>3}  {source[c]:<6} {od[c]:>5.2f}  {hybrid[c]:>8.3f}   {u}")
    print("\n" + json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
