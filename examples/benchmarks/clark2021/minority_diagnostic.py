#!/usr/bin/env python3
"""Why does the minority get over-excluded? -- the no-sim decomposition.

The hybrid dynamic gate lifted minority-share median 0.061 -> 0.133, but measured is
0.262, so muODE still over-excludes the loser ~2x. Before picking a mechanism to fix
it, this localizes the gap using data already on disk (no new sims): for every
scorable pair compare three minority shares --

  * OD-proportional:  min(OD_i,OD_j)/(OD_i+OD_j)   -- what the calibrated YIELDS alone
                      predict (static, zero dynamics)
  * muODE dynamic:    min(pa,pb)/(pa+pb)            -- the dFBA pair outcome
  * measured:         min(mpa,mpb)/(mpa+mpb)        -- reality

Two hypotheses, distinguished by where the gap sits:

  H-dynamics : muODE_dynamic << OD_proportional ~ measured
      -> the dFBA COMPETITION adds exclusion beyond the yields; reality is roughly
         yield-proportional. Fix = moderate the winner-take-all dynamics
         (enzyme-pool rate/yield trade-off is the right lever).

  H-yields   : muODE_dynamic ~ OD_proportional << measured
      -> the YIELDS themselves predict too much exclusion; reality is FLATTER than
         yield-proportion (an anti-exclusion / niche-partitioning mechanism reality
         has that yields don't). Enzyme-pool alone will NOT reach 0.26; need
         partitioning/density-dependence.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import clark2021 as ck

OUT = Path("results/clark2021/benchmark")


def main() -> int:
    obs = ck.observations(ck.load(str(OUT / "MasterDF.csv")),
                          diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg"))
    od = {o.species[0]: max(0.0, o.od) for o in obs if o.richness == 1}
    byp = {}
    for l in (OUT / "pairwise_interactions_bigg.jsonl").read_text().splitlines():
        r = json.loads(l)
        byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r
    done = {}
    for l in (OUT / "hybrid_gate_biomass.jsonl").read_text().splitlines():
        if l.strip():
            r = json.loads(l); done["-".join(sorted(r["codes"]))] = r["biomass"]

    rows = []
    for (a, b), d in byp.items():
        kp = "-".join(sorted((a, b)))
        if kp not in done or a not in done[kp] or b not in done[kp]:
            continue
        if a not in d or b not in d or a not in od or b not in od:
            continue
        pa, pb = done[kp][a], done[kp][b]
        mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
        oa, ob = od[a], od[b]
        if pa + pb <= 0 or mpa + mpb <= 0 or oa + ob <= 0:
            continue
        rows.append({
            "pair": kp,
            "muode_min": min(pa, pb) / (pa + pb),
            "od_min": min(oa, ob) / (oa + ob),
            "meas_min": min(mpa, mpb) / (mpa + mpb),
            # does muODE's dynamic winner match the OD (yield) winner?
            "muode_winner_is_od_winner": (pa > pb) == (oa > ob),
        })

    mu = np.array([r["muode_min"] for r in rows])
    odm = np.array([r["od_min"] for r in rows])
    ms = np.array([r["meas_min"] for r in rows])
    rep = {
        "n_pairs": len(rows),
        "minority_share_median": {
            "muode_dynamic": round(float(np.median(mu)), 3),
            "od_proportional": round(float(np.median(odm)), 3),
            "measured": round(float(np.median(ms)), 3)},
        "minority_share_mean": {
            "muode_dynamic": round(float(mu.mean()), 3),
            "od_proportional": round(float(odm.mean()), 3),
            "measured": round(float(ms.mean()), 3)},
        "gap_decomposition": {
            "yields_add_exclusion (od_prop below measured)": round(float(np.median(ms) - np.median(odm)), 3),
            "dynamics_add_exclusion (muode below od_prop)": round(float(np.median(odm) - np.median(mu)), 3)},
        "muode_dynamic_winner_matches_yield_winner_frac": round(float(np.mean(
            [r["muode_winner_is_od_winner"] for r in rows])), 3),
    }
    verdict = ("H-dynamics: dFBA competition over-excludes beyond yields; reality ~ "
               "yield-proportional -> enzyme-pool is the right lever"
               if (np.median(odm) - np.median(mu)) > (np.median(ms) - np.median(odm))
               else "H-yields: yields themselves over-exclude; reality flatter than "
               "yield-proportion -> need partitioning/density-dependence, not just enzyme-pool")
    rep["verdict"] = verdict
    (OUT / "minority_diagnostic.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
