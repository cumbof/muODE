#!/usr/bin/env python3
"""Community-parallel driver for the Clark 2021 Tiers 2/3 benchmark.

The stock ``run_benchmark.py`` scores every community in one process. Each community
is independent, peaks at <1 GB RSS, and takes minutes, so the run parallelises
cleanly *across communities as separate processes* -- which is GIL-immune (unlike
``DynamicFBA``'s in-community ``n_jobs`` threads, whose speedup depends on whether
GLPK releases the GIL) and memory-feasible on this box.

Reliability choices for a multi-hour run on a modest device:
  * **Checkpointed + resumable.** Every finished community is appended to a JSONL
    immediately; a re-run skips what is already there, so a kill/restart loses at
    most the in-flight communities.
  * **Biggest first.** Communities are scheduled largest-richness-first: longest-
    processing-time scheduling minimises makespan, and any memory blow-up on the
    23-member community happens in the first minutes, not 40 h in.
  * **Memory-aware worker count.** Defaults to what MemAvailable can hold at
    ~0.8 GB/worker, capped at cores-1.

Scoring reuses ``simulate`` and ``ck.score_metabolites`` unchanged, so the numbers
are identical to a serial ``run_benchmark.py`` run -- this only changes wall-time.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import clark2021 as ck            # noqa: E402
from muode.kinetics import KineticParameters             # noqa: E402
from run_benchmark import T_END, load_models, simulate   # noqa: E402

#: Per-process context, built once in each worker (avoids re-reading the diet/map
#: for every community). Populated by ``_init``.
_G: dict = {}


def _init(namespace: str, model_paths: dict[str, str]) -> None:
    """ProcessPoolExecutor initializer: silence GLPK, build the per-process context."""
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["mets"] = ck.metabolite_map(namespace)
    _G["kin"] = KineticParameters()
    _G["models"] = {c: Path(p) for c, p in model_paths.items()}


def _run_one(codes: list[str]) -> tuple[dict, dict]:
    """Worker task: one community -> (net metabolite production, relative abundance)."""
    return simulate(codes, _G["models"], _G["diet"], _G["kin"], _G["mets"], n_jobs=1)


def default_workers() -> int:
    """How many communities to run at once, bounded by free RAM (~0.8 GB each)."""
    cores = os.cpu_count() or 2
    avail_mb = 2000.0
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable"):
                avail_mb = int(line.split()[1]) / 1024.0
                break
    except Exception:
        pass
    by_mem = int((avail_mb - 600.0) / 800.0)
    return max(1, min(cores - 1, by_mem))


def load_checkpoint(path: Path) -> dict[str, dict]:
    """Read the JSONL checkpoint into {community: record}; tolerates a torn last line."""
    done: dict[str, dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue                                  # torn write from a kill; skip
            done[rec["community"]] = rec
    return done


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--data", type=Path, default=Path("results/clark2021/benchmark/MasterDF.csv"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--tier", choices=["2", "3", "all"], default="all")
    ap.add_argument("--max-richness", type=int, default=None)
    ap.add_argument("--workers", type=int, default=None, help="default: memory-aware")
    ap.add_argument("--limit", type=int, default=None, help="cap communities run (testing)")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    diet = ck.medium(args.namespace)
    mets = ck.metabolite_map(args.namespace)
    df = ck.load(args.data)
    models = load_models(args.models)
    model_paths = {c: str(p) for c, p in models.items()}

    want = {"2": (2, 2), "3": (3, 99), "all": (2, 99)}[args.tier]
    observed = [
        o for o in ck.observations(df, diet=diet, metabolites=mets)
        if want[0] <= o.richness <= want[1]
        and all(c in models for c in o.species)
        and (args.max_richness is None or o.richness <= args.max_richness)
    ]
    observed.sort(key=lambda o: o.richness, reverse=True)   # biggest first

    ckpt = args.outdir / f"tier23_partial_{args.namespace}_{args.tier}.jsonl"
    done = load_checkpoint(ckpt)
    todo = [o for o in observed if o.community not in done]
    if args.limit is not None:
        todo = todo[:args.limit]

    workers = args.workers or default_workers()
    print(f"namespace={args.namespace} tier={args.tier}  communities={len(observed)}  "
          f"done={len(done)}  todo={len(todo)}  workers={workers}", flush=True)

    if todo:
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(args.namespace, model_paths)) as ex:
            futs = {ex.submit(_run_one, list(o.species)): o for o in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                o = futs[fut]
                net, ab = fut.result()
                rec = {"community": o.community, "species": list(o.species),
                       "richness": o.richness, "net": net, "abundances": ab}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                done[o.community] = rec
                if n % 10 == 0 or n == len(todo):
                    print(f"  {n}/{len(todo)}  (last richness {o.richness})", flush=True)

    # -- score exactly like run_benchmark.py, over whatever is finished ----------
    predicted_mets = {c: rec["net"] for c, rec in done.items()}
    scored_obs = [o for o in observed if o.community in predicted_mets]
    obs_by_comm = {o.community: o for o in scored_obs}

    abundance_error = []
    for c in predicted_mets:
        o = obs_by_comm.get(c)
        if o is None or not o.abundances:
            continue
        shared = [s for s in o.species if s in o.abundances]
        if shared:
            ab = done[c]["abundances"]
            abundance_error.append(
                sum(abs(ab.get(s, 0.0) - o.abundances[s]) for s in shared) / len(shared))

    report = {
        "namespace": args.namespace, "tier": args.tier, "t_end_h": T_END,
        "community": {
            "n_communities": len(scored_obs),
            "metabolites": {
                exch: ck.score_metabolites(predicted_mets, scored_obs, exch, net=True)
                for exch in mets.values()
            },
            "abundance_mae": (sum(abundance_error) / len(abundance_error)
                              if abundance_error else None),
        },
    }
    out = args.outdir / f"tier23_report_{args.namespace}_{args.tier}.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {out}  ({len(scored_obs)}/{len(observed)} communities scored)", flush=True)
    for exch, sc in report["community"]["metabolites"].items():
        if "mae" not in sc:
            print(f"  {exch:10} not scored ({sc.get('error', 'no data')})")
            continue
        r = sc.get("pearson_r")
        print(f"  {exch:10} r={'n/a' if r is None else format(r, '.3f')}  "
              f"MAE={sc['mae']:.2f} mM  bias={sc['bias']:+.2f}  (n={sc['n']})")
    mae = report["community"]["abundance_mae"]
    if mae is not None:
        print(f"  relative-abundance MAE: {mae:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
