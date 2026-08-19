#!/usr/bin/env python3
"""Yield-based competitive ability -- the breakthrough diagnostic.

After per-species pH failed to reproduce Clark composition (the failure is rooted
in muODE's growth-rate RANKING, not the ecology layer), this asks whether the
competitive currency is monoculture YIELD (carrying capacity / OD), and whether
muODE simply mispredicts it.

Findings (data only, no sims):
  * measured monoculture OD predicts real dominance well (spearman ~0.58); muODE
    monoculture biomass barely does (~0.26) and barely tracks measured OD (~0.19).
  * higher measured OD predicts the pair winner 79% of the time; muODE biomass 52%.
  * abundance predicted as OD_i / sum(OD) scores MAE ~0.108 across 549 communities,
    vs uniform-null 0.164 and muODE dFBA 0.210 -- i.e. a monoculture-yield model
    with ZERO interaction terms beats both, using less information than gLV.

Conclusion: Clark composition is largely set by monoculture yields; muODE fails
because its FBA yields are wrong, not because interactions are missing. The fix is
to calibrate per-strain monoculture yield to the measured OD (independent of the
community outcomes) -- strictly less fitting than gLV, which also fits the pairwise
interaction coefficients.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, pearsonr

import clark2021 as ck

OUT = Path("results/clark2021/benchmark")


def main() -> int:
    obs = ck.observations(ck.load(str(OUT / "MasterDF.csv")),
                          diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg"))
    mono = {o.species[0]: o for o in obs if o.richness == 1}
    od = {c: max(0.0, o.od) for c, o in mono.items()}
    muode_bio = {json.loads(l)["code"]: json.loads(l)["mono"]["biomass"]
                 for l in (OUT / "mono_secretion_bigg.jsonl").read_text().splitlines() if l.strip()}

    # -- (1) dominance predictors -------------------------------------------------
    dom, cnt = {}, {}
    for o in obs:
        for c in o.species:
            if c in o.abundances:
                dom[c] = dom.get(c, 0) + o.abundances[c]; cnt[c] = cnt.get(c, 0) + 1
    dom = {c: dom[c] / cnt[c] for c in dom}
    codes = [c for c in muode_bio if c in od and c in dom]
    mb = np.array([muode_bio[c] for c in codes]); odv = np.array([od[c] for c in codes])
    dm = np.array([dom[c] for c in codes])

    # -- (2) pair-winner prediction ----------------------------------------------
    recs = [json.loads(l) for l in (OUT / "pairwise_interactions_bigg.jsonl").read_text().splitlines()]
    byp = {}
    for r in recs:
        byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r
    od_hit, mu_hit, npair = 0, 0, 0
    for (a, b), d in byp.items():
        if a not in d or b not in d or a not in od or b not in od:
            continue
        mp = d[a]["meas_pair"] + d[b]["meas_pair"]
        if mp <= 0:
            continue
        mw = a if d[a]["meas_pair"] > d[b]["meas_pair"] else b
        od_hit += (a if od[a] > od[b] else b) == mw
        mu_hit += (a if muode_bio[a] > muode_bio[b] else b) == mw
        npair += 1

    # -- (3) composition MAE: OD-weighted vs null vs muODE ------------------------
    def mae(weight):
        errs = []
        for o in obs:
            if o.richness < 2 or not o.abundances:
                continue
            sh = [c for c in o.species if c in o.abundances and c in od]
            if not sh:
                continue
            w = weight(sh); tot = sum(w.values()) or 1.0
            p = {c: w[c] / tot for c in sh}
            errs.append(sum(abs(p[c] - o.abundances[c]) for c in sh) / len(sh))
        return float(np.mean(errs)), len(errs)

    mae_uniform, n = mae(lambda sh: {c: 1.0 for c in sh})
    mae_od, _ = mae(lambda sh: {c: od[c] for c in sh})

    report = {
        "dominance_vs_muODE_biomass_spearman": round(float(spearmanr(dm, mb)[0]), 3),
        "dominance_vs_measured_OD_spearman": round(float(spearmanr(dm, odv)[0]), 3),
        "muODE_biomass_vs_measured_OD_spearman": round(float(spearmanr(mb, odv)[0]), 3),
        "pair_winner_accuracy": {
            "by_measured_OD": round(od_hit / npair, 3),
            "by_muODE_biomass": round(mu_hit / npair, 3), "n_pairs": npair},
        "composition_abundance_MAE": {
            "uniform_null": round(mae_uniform, 3),
            "measured_OD_weighted": round(mae_od, 3),
            "muODE_dFBA_reference": 0.210, "n_communities": n},
    }
    (OUT / "yield_competition_bigg.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
