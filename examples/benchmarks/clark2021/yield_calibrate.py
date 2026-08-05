#!/usr/bin/env python3
"""Yield calibration -- fit per-strain ATP maintenance to measured monoculture OD.

The breakthrough (yield_competition.py): Clark composition is largely set by
monoculture YIELD (OD-weighted abundance beats the null and muODE's dFBA), and
muODE fails only because its FBA yields are wrong. The existing calibration tunes
Vmax (rate), which is unidentifiable; the biomass YIELD was never calibrated, and
every GEM ships with ATPM lower bound = 0 (no maintenance -> maximal, unrealistic
yield).

This fits the one mechanistic knob that sets yield -- non-growth ATP maintenance
(ATPM) -- per strain, so the monoculture dFBA biomass matches the measured OD.
Uses ONLY monoculture data (independent of the community/pairwise outcomes we
predict); strictly less fitting than gLV, which also fits pairwise coefficients.
ATPM can only LOWER yield, so it corrects over-growers (ER, RI, BU...) exactly;
under-growers (PC, BL: a substrate-access limit, not maintenance) stay at ATPM=0
and their residual is reported honestly, not hidden.

Two modes:
  --calibrate : fit ATPM per strain -> atpm_calibration_bigg.json
  --gate      : re-run 25 monos + 126 pairs with calibrated ATPM, score vs measured
                (pair-winner accuracy toward the 79% OD ceiling; interaction corr;
                 minority share) and vs the baseline dFBA.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import cobra

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck              # noqa: E402
from muode.community import Community                      # noqa: E402
from muode.dfba import DynamicFBA                          # noqa: E402
from muode.kinetics import KineticParameters              # noqa: E402
from muode.organism import CobraOrganism                  # noqa: E402
from run_benchmark import T_END, load_models              # noqa: E402

ATPM_GRID = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
_G: dict = {}


def _init(namespace, model_paths, target_od, atpm=None):
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["kin"] = KineticParameters()
    _G["paths"] = model_paths
    _G["target"] = target_od
    _G["atpm"] = atpm or {}


def _organism(code, atpm_value):
    model = cobra.io.read_sbml_model(str(_G["paths"][code]))
    if "ATPM" in [r.id for r in model.reactions]:
        model.reactions.ATPM.lower_bound = float(atpm_value)
    return CobraOrganism(model, id=code)


def _mono_biomass(org, codes=None):
    code = org.id
    comm = Community([org], abundances={code: 1.0}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, _G["diet"], _G["kin"])
    return float(res.biomass.iloc[-1].to_dict().get(code, 0.0))


def _calibrate_one(code):
    """Adaptive: run ATPM=0 first (1 dFBA for under-growers), then raise ATPM only
    until the target OD is bracketed, then interpolate. Biomass is monotone-down
    in ATPM, so an early stop at the first point below target is exact enough."""
    model = cobra.io.read_sbml_model(str(_G["paths"][code]))
    has_atpm = "ATPM" in [r.id for r in model.reactions]
    target = _G["target"].get(code, 0.0)

    def biomass(a):
        if has_atpm:
            model.reactions.ATPM.lower_bound = a
        org = CobraOrganism(model, id=code)
        comm = Community([org], abundances={code: 1.0}, total_biomass=0.01)
        res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, _G["diet"], _G["kin"])
        return float(res.biomass.iloc[-1].to_dict().get(code, 0.0))

    curve = [(0.0, biomass(0.0))]
    b0 = curve[0][1]
    fitted = 0.0
    if has_atpm and b0 > target:                       # over-grower: sweep up to bracket
        fitted = ATPM_GRID[-1]                          # cap default
        for a in ATPM_GRID[1:]:
            b = biomass(a)
            curve.append((a, b))
            if b <= target:
                a1, b1 = curve[-2]
                fitted = a1 + (a - a1) * (b1 - target) / (b1 - b) if b1 != b else a
                break
    return code, {"atpm": round(fitted, 3), "target_od": round(target, 3),
                  "biomass_at_atpm0": round(b0, 3), "curve": curve, "has_atpm": has_atpm}


def _gate_one(codes):
    orgs = [_organism(c, _G["atpm"].get(c, 0.0)) for c in codes]
    comm = Community(orgs, abundances={c: 1.0 / len(codes) for c in codes}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, _G["diet"], _G["kin"])
    return codes, {c: float(v) for c, v in res.biomass.iloc[-1].to_dict().items()}


def default_workers():
    cores = os.cpu_count() or 2
    avail = 2000.0
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable"):
                avail = int(line.split()[1]) / 1024.0
                break
    except Exception:
        pass
    return max(1, min(cores - 1, int((avail - 600.0) / 500.0)))


def measured_od(namespace):
    obs = ck.observations(ck.load("results/clark2021/benchmark/MasterDF.csv"),
                          diet=ck.medium(namespace), metabolites=ck.metabolite_map(namespace))
    return {o.species[0]: max(0.0, o.od) for o in obs if o.richness == 1}


def score_gate(done, namespace, outdir):
    import math
    import numpy as np
    from scipy.stats import spearmanr, pearsonr
    base = [json.loads(l) for l in (outdir / f"pairwise_interactions_{namespace}.jsonl").read_text().splitlines()]
    byp = {}
    for r in base:
        byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r
    meas, basep, newp = [], [], []
    for r in base:
        i, j = r["focal"], r["partner"]
        ka, kp = "-".join([i]), "-".join(sorted((i, j)))
        if ka not in done or kp not in done:
            continue
        ba, bp = done[ka]["biomass"].get(i, 0.0), done[kp]["biomass"].get(i, 0.0)
        if ba <= 0:
            continue
        meas.append(r["measured"]); basep.append(r["predicted"])
        newp.append(math.log(max(bp, 1e-9) / ba) + math.log(2))
    win_new, win_base, mino_new, mino_meas = 0, 0, [], []
    n = 0
    for (a, b), d in byp.items():
        kp = "-".join(sorted((a, b)))
        if kp not in done or a not in done[kp]["biomass"] or b not in done[kp]["biomass"]:
            continue
        if a not in d or b not in d:
            continue
        pa, pb = done[kp]["biomass"][a], done[kp]["biomass"][b]
        mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
        if mpa + mpb <= 0 or pa + pb <= 0:
            continue
        mw = a if mpa > mpb else b
        win_new += (a if pa > pb else b) == mw
        ba0 = done.get("-".join([a]), {}).get("biomass", {}).get(a, 0)
        bb0 = done.get("-".join([b]), {}).get("biomass", {}).get(b, 0)
        win_base += (a if ba0 > bb0 else b) == mw   # not used; kept for parity
        mino_new.append(min(pa, pb) / (pa + pb)); mino_meas.append(min(mpa, mpb) / (mpa + mpb))
        n += 1
    meas, basep, newp = np.array(meas), np.array(basep), np.array(newp)
    rep = {
        "namespace": namespace, "n_interactions": int(len(meas)), "n_pairs": n,
        "interaction_corr_with_measured": {
            "baseline_dFBA": round(float(spearmanr(basep, meas)[0]), 3),
            "yield_calibrated": round(float(spearmanr(newp, meas)[0]), 3)},
        "pair_winner_accuracy": {
            "baseline_dFBA": 0.523, "yield_calibrated": round(win_new / n, 3),
            "OD_ceiling": 0.791},
        "minority_share_median": {
            "baseline_dFBA": 0.061, "yield_calibrated": round(float(np.median(mino_new)), 3),
            "measured": round(float(np.median(mino_meas)), 3)},
    }
    (outdir / f"yield_gate_report_{namespace}.json").write_text(json.dumps(rep, indent=2))
    print("\n=== YIELD-CALIBRATED GATE ===")
    print(json.dumps(rep, indent=2))
    return rep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--proxy", action="store_true",
                    help="static yield-allocation proxy from the calibration (no sims)")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    models = load_models(args.models)
    model_paths = {c: str(p) for c, p in models.items()}
    target = measured_od(args.namespace)
    workers = args.workers or default_workers()
    cal_path = args.outdir / f"atpm_calibration_{args.namespace}.jsonl"

    if args.calibrate:
        done = set()
        if cal_path.exists():
            done = {json.loads(l)["code"] for l in cal_path.read_text().splitlines() if l.strip()}
        todo = [c for c in sorted(models) if c not in done]
        print(f"calibrate: {len(todo)} strains, workers={workers}", flush=True)
        with open(cal_path, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(args.namespace, model_paths, target)) as ex:
            futs = {ex.submit(_calibrate_one, c): c for c in todo}
            for k, fut in enumerate(as_completed(futs), 1):
                code, rec = fut.result()
                fh.write(json.dumps({"code": code, **rec}) + "\n"); fh.flush()
                print(f"  {k}/{len(todo)} {code}: OD={rec['target_od']} "
                      f"bio@0={rec['biomass_at_atpm0']} -> ATPM={rec['atpm']}", flush=True)

    if args.proxy:
        import numpy as np
        from scipy.stats import spearmanr
        cal = {}
        for l in cal_path.read_text().splitlines():
            if not l.strip():
                continue
            r = json.loads(l); curve = r["curve"]; a = r["atpm"]
            if len(curve) == 1:
                b = curve[0][1]
            else:
                b = curve[-1][1]
                for (a1, b1), (a2, b2) in zip(curve, curve[1:]):
                    if a1 <= a <= a2:
                        b = b1 + (b2 - b1) * ((a - a1) / (a2 - a1) if a2 != a1 else 0); break
            cal[r["code"]] = max(0.0, b)
        raw = {json.loads(l)["code"]: json.loads(l)["mono"]["biomass"]
               for l in (args.outdir / f"mono_secretion_{args.namespace}.jsonl").read_text().splitlines() if l.strip()}
        obs = ck.observations(ck.load("results/clark2021/benchmark/MasterDF.csv"),
                              diet=ck.medium(args.namespace), metabolites=ck.metabolite_map(args.namespace))
        od = {o.species[0]: max(0., o.od) for o in obs if o.richness == 1}

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

        codes = [c for c in cal if c in od and c in raw]
        m_null, n = mae({c: 1.0 for c in od})
        rep = {
            "mono_biomass_vs_OD_spearman": {
                "raw": round(float(spearmanr([raw[c] for c in codes], [od[c] for c in codes])[0]), 3),
                "calibrated": round(float(spearmanr([cal[c] for c in codes], [od[c] for c in codes])[0]), 3)},
            "abundance_MAE": {
                "uniform_null": m_null, "raw_muODE_biomass": mae(raw)[0],
                "ATPM_calibrated": mae(cal)[0], "measured_OD_ceiling": mae(od)[0],
                "muODE_dFBA_reference": 0.210, "n_communities": n},
        }
        (args.outdir / f"yield_proxy_report_{args.namespace}.json").write_text(json.dumps(rep, indent=2))
        print(json.dumps(rep, indent=2))
        return 0

    if args.gate:
        atpm = {json.loads(l)["code"]: json.loads(l)["atpm"]
                for l in cal_path.read_text().splitlines() if l.strip()}
        obs = ck.observations(ck.load("results/clark2021/benchmark/MasterDF.csv"),
                              diet=ck.medium(args.namespace), metabolites=ck.metabolite_map(args.namespace))
        monos = sorted({o.species[0] for o in obs if o.richness == 1 and o.species[0] in models})
        pairs = sorted({tuple(sorted(o.species)) for o in obs
                        if o.richness == 2 and all(s in models for s in o.species)})
        units = [[c] for c in monos] + [list(p) for p in pairs]
        ckpt = args.outdir / f"yield_biomass_{args.namespace}.jsonl"
        done = {}
        if ckpt.exists():
            for l in ckpt.read_text().splitlines():
                if l.strip():
                    r = json.loads(l); done["-".join(sorted(r["codes"]))] = r
        todo = [u for u in units if "-".join(sorted(u)) not in done]
        print(f"gate: units={len(units)} done={len(done)} todo={len(todo)} workers={workers}", flush=True)
        if todo:
            with open(ckpt, "a") as fh, ProcessPoolExecutor(
                    max_workers=workers, initializer=_init,
                    initargs=(args.namespace, model_paths, target, atpm)) as ex:
                futs = {ex.submit(_gate_one, list(u)): u for u in todo}
                for k, fut in enumerate(as_completed(futs), 1):
                    codes, bio = fut.result()
                    rec = {"codes": sorted(codes), "biomass": bio}
                    fh.write(json.dumps(rec) + "\n"); fh.flush()
                    done["-".join(sorted(codes))] = rec
                    if k % 10 == 0 or k == len(todo):
                        print(f"  {k}/{len(todo)}", flush=True)
        score_gate(done, args.namespace, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
