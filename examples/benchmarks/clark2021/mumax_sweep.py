#!/usr/bin/env python3
"""Pin the physiological mu_max ceiling on a representative pair subset.

The full mu_max=0.35 gate showed the cap works but OVER-corrects (running minority
median ~0.43 vs measured 0.262) -- 0.35/h is too tight. Rather than run sequential
2.5 h full gates guessing values, this sweeps mu_max on the same representative
subset and reports two metrics that must move together:

  * minority-share median   -> should approach the measured 0.262 (too high = the
    cap over-equalizes; too low = winner-take-all not yet relaxed)
  * pair-winner accuracy     -> must STAY high (~0.73 baseline, 0.79 OD ceiling):
    an over-tight cap that pushes every pair to 50/50 lifts minority but destroys
    the winner prediction. The right mu_max maximizes coexistence WITHOUT losing
    directionality.

The subset deliberately mixes over-excluded (coexisting) pairs and reality-skewed
(dominated) pairs so the sweep is penalized for wrongly equalizing a true winner.
Pick the mu_max that best matches BOTH, then run one full gate at that value.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import logging
logging.disable(logging.CRITICAL)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck                # noqa: E402
import hybrid_gate as HG                                    # noqa: E402

OUT = Path("results/clark2021/benchmark")


def scorable_pairs():
    """Measured pairs with both strains in the panel and measured abundances."""
    byp = {}
    for l in (OUT / "pairwise_interactions_bigg.jsonl").read_text().splitlines():
        r = json.loads(l)
        byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r
    paths, _, _ = HG.build_panel(Path("examples/benchmarks/clark2021/models"),
                                 Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    out = {}
    for (a, b), d in byp.items():
        if a in paths and b in paths and a in d and b in d:
            mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
            if mpa + mpb > 0:
                out[(a, b)] = min(mpa, mpb) / (mpa + mpb)   # measured minority
    return out, byp


def representative_subset(meas_min, n):
    """Even spread across the measured-minority range (dominated -> coexisting),
    so the subset spans reality-skewed and over-excluded pairs alike."""
    ordered = sorted(meas_min, key=lambda p: meas_min[p])
    if n >= len(ordered):
        return ordered
    idx = [round(i * (len(ordered) - 1) / (n - 1)) for i in range(n)]
    return [ordered[i] for i in sorted(set(idx))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mumaxes", default="0.5,0.7,1.0")
    ap.add_argument("--n", type=int, default=18)
    args = ap.parse_args()
    mumaxes = [float(x) for x in args.mumaxes.split(",")]

    meas_min, byp = scorable_pairs()
    subset = representative_subset(meas_min, args.n)
    paths, source, atpm = HG.build_panel(Path("examples/benchmarks/clark2021/models"),
                                         Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    HG._init(paths, atpm)          # sets _G; we override mumax per sweep value below

    print(f"subset: {len(subset)} pairs (measured minority "
          f"{meas_min[subset[0]]:.2f}..{meas_min[subset[-1]]:.2f})")
    print(f"measured minority median (subset) = "
          f"{statistics.median([meas_min[p] for p in subset]):.3f}\n")
    print(f"{'mu_max':>7} {'minority_med':>12} {'pair_winner':>12}   (measured min 0.262)")

    results = {}
    for mm in mumaxes:
        HG._G["mumax"] = mm
        mins, winhit, npair = [], 0, 0
        rows = {}
        for (a, b) in subset:
            _, bio = HG._gate_one([a, b])
            pa, pb = bio.get(a, 0.0), bio.get(b, 0.0)
            tot = pa + pb
            if tot <= 0:
                continue
            mn = min(pa, pb) / tot
            mins.append(mn)
            d = byp[(a, b)]
            mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
            mw = a if mpa > mpb else b
            pw = a if pa > pb else b
            winhit += (pw == mw); npair += 1
            rows["-".join((a, b))] = {"pred_min": round(mn, 3),
                                      "meas_min": round(meas_min[(a, b)], 3)}
        med = statistics.median(mins) if mins else float("nan")
        acc = winhit / npair if npair else float("nan")
        results[mm] = {"minority_median": round(med, 3), "pair_winner": round(acc, 3),
                       "n": npair, "pairs": rows}
        print(f"{mm:>7} {med:>12.3f} {acc:>12.3f}", flush=True)

    (OUT / "mumax_sweep.json").write_text(json.dumps(
        {"measured_minority_median_subset":
         round(statistics.median([meas_min[p] for p in subset]), 3),
         "baseline_pair_winner": 0.733, "OD_ceiling_winner": 0.791,
         "results": results}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
