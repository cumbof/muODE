#!/usr/bin/env python3
"""Full 549-community DYNAMIC composition MAE -- the flagship end-to-end test.

Everything so far scored monocultures + pairs, or a STATIC OD-weighting proxy
(abundance proportional to calibrated yield, MAE 0.118). This is the real thing:
run the actual dynamic dFBA on every measured Clark community (richness >= 2) with
the full mechanism -- per-strain best reconstruction (harmonized gapseq Bacteroidetes
+ BiGG), ATPM yield calibration, physiological mu_max cap, pFBA for determinism, one
shared BiGG DM38 pool -- and score predicted composition (biomass_i / sum biomass)
against measured relative abundances.

Scored against the reference points established across the project:
    uniform null          0.164
    muODE raw dFBA        0.210   (worse than null -- where we started)
    static hybrid proxy   0.118
    measured-OD ceiling   0.108
The question: does the mechanism, run dynamically at full community richness, beat
the null and approach the static proxy / OD ceiling?

Reuses hybrid_gate's _init/_gate_one (mumax + pfba via _G). Simulates each unique
species set once (communities repeat across replicate observations), checkpoints,
then scores the MAE over all qualifying observations exactly as the proxy did.
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck                # noqa: E402
import hybrid_gate as HG                                    # noqa: E402


def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bigg-models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--harm-models", type=Path, default=Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--mumax", type=float, default=0.6)
    ap.add_argument("--pfba", action="store_true", default=True)
    ap.add_argument("--no-pfba", dest="pfba", action="store_false")
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    paths, source, atpm = HG.build_panel(args.bigg_models, args.harm_models)
    workers = args.workers or HG.default_workers()

    obs = ck.observations(ck.load(str(args.outdir / "MasterDF.csv")),
                          diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg"))
    # qualifying observations: richness>=2, measured abundances, all strains in panel
    qual = [o for o in obs if o.richness >= 2 and o.abundances
            and all(s in paths for s in o.species)]
    # unique species sets to simulate (replicates share a set)
    unique = sorted({tuple(sorted(o.species)) for o in qual})
    print(f"community MAE: {len(qual)} scorable observations, {len(unique)} unique "
          f"species sets | mumax={args.mumax} pfba={args.pfba} workers={workers}", flush=True)

    ckpt = args.outdir / f"community_biomass_mumax{args.mumax:g}{'_pfba' if args.pfba else ''}.jsonl"
    done: dict = {}
    if ckpt.exists():
        for l in ckpt.read_text().splitlines():
            if l.strip():
                r = json.loads(l); done["-".join(sorted(r["codes"]))] = r["biomass"]
    # Longest-processing-time-first scheduling: start the heaviest (highest-richness)
    # communities first so an R23 monster runs packed among the busy first wave
    # instead of alone on one core at the end (which would dominate wall time).
    todo = sorted((u for u in unique if "-".join(sorted(u)) not in done),
                  key=len, reverse=True)
    print(f"  unique sets: done={len(done)} todo={len(todo)} "
          f"(heaviest-first: {[len(u) for u in todo[:5]]}...)", flush=True)

    if todo:
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=HG._init,
                initargs=(paths, atpm, None, args.mumax, None, args.pfba)) as ex:
            futs = {ex.submit(HG._gate_one, list(u)): u for u in todo}
            for k, fut in enumerate(as_completed(futs), 1):
                codes, bio = fut.result()
                fh.write(json.dumps({"codes": sorted(codes), "biomass": bio}) + "\n"); fh.flush()
                done["-".join(sorted(codes))] = bio
                if k % 25 == 0 or k == len(todo):
                    print(f"    {k}/{len(todo)}", flush=True)

    # ---- score composition MAE over all qualifying observations -----------------
    errs = []
    for o in qual:
        key = "-".join(sorted(o.species))
        bio = done.get(key)
        if not bio:
            continue
        sh = [c for c in o.species if c in o.abundances]
        tot = sum(max(0.0, bio.get(c, 0.0)) for c in sh)
        if tot <= 0:
            # community failed to grow: predict uniform (least-committal)
            pred = {c: 1.0 / len(sh) for c in sh}
        else:
            pred = {c: max(0.0, bio.get(c, 0.0)) / tot for c in sh}
        errs.append(sum(abs(pred[c] - o.abundances[c]) for c in sh) / len(sh))

    mae = float(np.mean(errs)) if errs else float("nan")
    rep = {
        "n_observations_scored": len(errs),
        "n_unique_communities": len(unique),
        "dynamic_community_MAE": round(mae, 4),
        "reference": {"uniform_null": 0.164, "muODE_raw_dFBA": 0.210,
                      "static_hybrid_proxy": 0.118, "measured_OD_ceiling": 0.108},
        "config": {"mumax": args.mumax, "pfba": args.pfba},
    }
    (args.outdir / f"community_mae_report_mumax{args.mumax:g}{'_pfba' if args.pfba else ''}.json"
     ).write_text(json.dumps(rep, indent=2))
    print("\n=== DYNAMIC COMMUNITY COMPOSITION MAE ===")
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
