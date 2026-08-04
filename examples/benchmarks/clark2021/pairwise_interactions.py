#!/usr/bin/env python3
"""Pairwise-interaction validation -- the go/no-go gate for an FBA-parameterized
ecological (gLV) layer.

Context
-------
The Clark 2021 community compositions are reproduced in silico only by *fitted*
ecological models: the paper's own generalized Lotka-Volterra (gLV) learns its
interaction coefficients ``a_ij`` by regression from monoculture + pairwise
co-culture data, and the leading community-FBA tool (MICOM) takes relative
abundance as an *input*, not an output.  A pure de-novo dFBA (muODE's Tier-2/3
run) does not beat the uniform null at composition -- an open problem for the
whole mechanistic-FBA class, not a muODE bug.

The proposed fix is to let muODE's dFBA *derive* the gLV interaction coefficients
mechanistically and then run the gLV forward for composition.  Before building
that layer, this script asks the one falsifiable question that gates it:

    Do muODE's mechanistic pairwise interactions agree with the *measured* ones?

If yes, the ecology layer is licensed (mechanistic interactions can stand in for
the fitted ones).  If no, we have found exactly where the mechanism breaks -- in
~1 h of small (mono + 2-species) sims, not a 45 h community run.

What it computes
----------------
For every strain that grew in monoculture, and every measured pair, the
interaction of j on i as a log growth-ratio:

    interaction[i<-j] = log( abundance_i(in pair i,j) / abundance_i(alone) )

both for muODE (absolute final biomass, ``result.biomass.iloc[-1]``, which the
scored pipeline normalizes away) and for the data (absolute abundance = OD x 16S
fraction).  Sign: <0 = j suppresses i (competition), >0 = j facilitates i
(cross-feeding), ~0 = neutral.  A ``+log(2)`` initial-condition correction sets
the neutral point at 0 (both experiment and sim inoculate each pair member at
half its monoculture inoculum, assuming equal-fraction inoculation); this
constant does not affect the correlation, only the zero-point / sign stats.

Reuses ``simulate``'s exact organism/community/dFBA construction so the numbers
match the benchmark; checkpointed + resumable like ``tier23_parallel.py``.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck              # noqa: E402
from muode.community import Community                      # noqa: E402
from muode.dfba import DynamicFBA                          # noqa: E402
from muode.kinetics import KineticParameters              # noqa: E402
from muode.organism import CobraOrganism                  # noqa: E402
from run_benchmark import T_END, load_models              # noqa: E402

#: total inoculum (matches ``simulate``); split equally across pair members.
TOTAL_BIOMASS = 0.01
#: floors below which a monoculture "did not grow" and its ratio is meaningless.
BIO_FLOOR = 1e-4          # absolute biomass (sim); >~1% of the 0.01 inoculum
OD_FLOOR = 0.05           # measured monoculture abundance (OD x fraction)

_G: dict = {}


def _init(namespace: str, model_paths: dict[str, str]) -> None:
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["mets"] = ck.metabolite_map(namespace)
    _G["kin"] = KineticParameters()
    _G["models"] = {c: Path(p) for c, p in model_paths.items()}


def sim_biomass(codes, models, diet, kinetics, t_end=T_END, dt=0.1, n_jobs=1):
    """Run one community; return absolute final biomass per code.

    Identical construction to ``run_benchmark.simulate`` (same inoculum, same
    equal starting fractions), but returns the *absolute* endpoint biomass the
    scored pipeline discards -- the quantity a growth-ratio interaction needs.
    """
    organisms = [CobraOrganism.from_file(str(models[c]), id=c) for c in codes]
    community = Community(organisms, abundances={c: 1.0 / len(codes) for c in codes},
                          total_biomass=TOTAL_BIOMASS)
    result = DynamicFBA(t_end=t_end, dt=dt, n_jobs=n_jobs).run(community, diet, kinetics)
    return {c: float(v) for c, v in result.biomass.iloc[-1].to_dict().items()}


def _run_one(codes: list[str]) -> tuple[list[str], dict]:
    return codes, sim_biomass(codes, _G["models"], _G["diet"], _G["kin"])


def default_workers() -> int:
    cores = os.cpu_count() or 2
    avail_mb = 2000.0
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable"):
                avail_mb = int(line.split()[1]) / 1024.0
                break
    except Exception:
        pass
    by_mem = int((avail_mb - 600.0) / 500.0)   # mono/pair peak ~0.5 GB
    return max(1, min(cores - 1, by_mem))


def load_checkpoint(path: Path) -> dict[str, dict]:
    """Read {'-'.join(sorted(codes)) -> biomass dict}; tolerates a torn last line."""
    done: dict[str, dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[rec["key"]] = rec["biomass"]
    return done


def _key(codes) -> str:
    return "-".join(sorted(codes))


def build_units(observed, models):
    """Monocultures + measured pairs whose members all have a GEM.

    Uses the *measured* pairs (135 of the 325 possible) so every muODE prediction
    has a data counterpart to be scored against.
    """
    have = set(models)
    monos = sorted({o.species[0] for o in observed
                    if o.richness == 1 and o.species[0] in have})
    pairs = sorted({tuple(sorted(o.species)) for o in observed
                    if o.richness == 2 and all(s in have for s in o.species)})
    units = [[c] for c in monos] + [list(p) for p in pairs]
    return monos, pairs, units


def measured_abundance(observed):
    """{community key -> {code -> absolute abundance = OD x 16S fraction}}."""
    out = {}
    for o in observed:
        if o.od != o.od:                                   # NaN OD
            continue
        out[_key(o.species)] = {c: o.od * o.abundances[c]
                                for c in o.species if c in o.abundances}
    return out


def interactions(monos, pairs, sim_bio, meas_abund):
    """Build paired (measured, predicted) interaction records for every direction."""
    log2 = math.log(2.0)
    recs = []
    for (a, b) in pairs:
        pk = _key((a, b))
        sim_pair = sim_bio.get(pk)
        meas_pair = meas_abund.get(pk)
        if sim_pair is None or not meas_pair:
            continue
        for i, j in ((a, b), (b, a)):                      # j acts on i
            sim_mono = sim_bio.get(_key([i]))
            meas_mono = meas_abund.get(_key([i]))
            if not sim_mono or not meas_mono:
                continue
            si_alone = sim_mono.get(i, 0.0)
            si_pair = sim_pair.get(i, 0.0)
            mi_alone = meas_mono.get(i, 0.0)
            mi_pair = meas_pair.get(i, 0.0)
            # require growth alone in BOTH sim and data -- else the ratio is noise
            if si_alone < BIO_FLOOR or mi_alone < OD_FLOOR:
                continue
            pred = math.log(max(si_pair, 1e-9) / si_alone) + log2
            meas = math.log(max(mi_pair, 1e-9) / mi_alone) + log2
            recs.append({"focal": i, "partner": j,
                         "predicted": pred, "measured": meas,
                         "sim_alone": si_alone, "sim_pair": si_pair,
                         "meas_alone": mi_alone, "meas_pair": mi_pair})
    return recs


def correlate(recs):
    import numpy as np
    from scipy.stats import pearsonr, spearmanr
    p = np.array([r["predicted"] for r in recs])
    m = np.array([r["measured"] for r in recs])
    sign_agree = float(np.mean(np.sign(p) == np.sign(m))) if len(p) else float("nan")
    out = {"n": int(len(p)), "sign_agreement": sign_agree,
           "frac_measured_negative": float(np.mean(m < 0)) if len(m) else float("nan"),
           "frac_predicted_negative": float(np.mean(p < 0)) if len(p) else float("nan")}
    if len(p) >= 3:
        pr, pp = pearsonr(p, m)
        sr, sp = spearmanr(p, m)
        out.update({"pearson_r": float(pr), "pearson_p": float(pp),
                    "spearman_r": float(sr), "spearman_p": float(sp)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--data", type=Path, default=Path("results/clark2021/benchmark/MasterDF.csv"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--workers", type=int, default=None, help="default: memory-aware")
    ap.add_argument("--limit", type=int, default=None, help="cap units run (smoke test)")
    ap.add_argument("--analyze-only", action="store_true",
                    help="skip sims; just (re)build interactions from the checkpoint")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    diet = ck.medium(args.namespace)
    mets = ck.metabolite_map(args.namespace)
    df = ck.load(args.data)
    observed = ck.observations(df, diet=diet, metabolites=mets)
    models = load_models(args.models)
    model_paths = {c: str(p) for c, p in models.items()}

    monos, pairs, units = build_units(observed, models)
    print(f"namespace={args.namespace}  monocultures={len(monos)}  "
          f"measured pairs={len(pairs)}  units={len(units)}", flush=True)

    ckpt = args.outdir / f"pairwise_biomass_{args.namespace}.jsonl"
    done = load_checkpoint(ckpt)
    todo = [u for u in units if _key(u) not in done]
    if args.limit is not None:
        todo = todo[:args.limit]

    if not args.analyze_only and todo:
        workers = args.workers or default_workers()
        print(f"done={len(done)}  todo={len(todo)}  workers={workers}", flush=True)
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(args.namespace, model_paths)) as ex:
            futs = {ex.submit(_run_one, list(u)): u for u in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                codes, biomass = fut.result()
                fh.write(json.dumps({"key": _key(codes), "codes": sorted(codes),
                                     "biomass": biomass}) + "\n")
                fh.flush()
                done[_key(codes)] = biomass
                if n % 10 == 0 or n == len(todo):
                    print(f"  {n}/{len(todo)}", flush=True)

    # -- analysis ---------------------------------------------------------------
    meas_abund = measured_abundance(observed)
    recs = interactions(monos, pairs, done, meas_abund)
    stats = correlate(recs)

    report = {"namespace": args.namespace, "t_end_h": T_END,
              "n_monocultures": len(monos), "n_measured_pairs": len(pairs),
              "n_units_simulated": len(done), "interaction_stats": stats}
    (args.outdir / f"pairwise_report_{args.namespace}.json").write_text(
        json.dumps(report, indent=2, default=str))
    with open(args.outdir / f"pairwise_interactions_{args.namespace}.jsonl", "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")

    print("\n=== pairwise interaction validation ===", flush=True)
    for k, v in stats.items():
        print(f"  {k:24} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
