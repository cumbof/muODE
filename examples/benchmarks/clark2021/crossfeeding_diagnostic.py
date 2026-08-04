#!/usr/bin/env python3
"""Cross-feeding diagnostic -- WHY muODE over-predicts competitive exclusion.

The pairwise gate (``pairwise_interactions.py``) showed muODE's dFBA pairwise
interactions are ~uncorrelated with the measured ones and *over-produce
competition*: 39% of pairs called competitive vs 26% real, where 74% of real
pairs coexist.  Two hypotheses for the over-competition:

  (A) starved cross-feeding -- the BiGG GEMs under-secrete fermentation products
      (the Tier-1 alternate-optima limitation), so the facilitation that lets
      real pairs coexist never flows, leaving only competition for shared carbon.
      => fixable in per-organism FBA flux routing.
  (B) genuine hard resource competition -- pairs really do fight over the same
      substrates and muODE captures that; reality softens it by niche/temporal
      partitioning muODE lacks.  => needs ecological mechanism (diauxie), not
      secretion.

This script tests them.  It (1) re-runs the 25 monocultures capturing the FULL
net secretion profile + the strain's importable exchange set, then (2) builds a
directional cross-feeding potential  xfeed(i<-j) = sum of what j secretes that i
can take up, and asks whether it explains the measured facilitation and the
muODE residual (predicted - measured).  It also quantifies under-secretion
directly: muODE monoculture acid production vs the MEASURED values in MasterDF.

Read:
  * if measured interaction rises with xfeed potential but the muODE residual
    falls (muODE increasingly *too competitive* where cross-feeding should
    happen) AND muODE under-secretes the acids -> hypothesis (A), the bottleneck
    is FBA flux routing / secretion, not the ecology layer.
  * if xfeed potential is unrelated to measured coexistence -> hypothesis (B),
    competition is real and the fix is niche/diauxie in the ecology layer.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck              # noqa: E402
from muode.community import Community                      # noqa: E402
from muode.dfba import DynamicFBA                          # noqa: E402
from muode.kinetics import KineticParameters              # noqa: E402
from muode.organism import CobraOrganism                  # noqa: E402
from run_benchmark import T_END, load_models              # noqa: E402

TOTAL_BIOMASS = 0.01
SECR_FLOOR = 0.3          # mM net production to count as "secreted"
#: canonical cross-feedable fermentation products (BiGG ids)
FERMENT = ["ac_e", "lac__L_e", "lac__D_e", "succ_e", "for_e", "etoh_e",
           "ppa_e", "but_e", "h2_e", "co2_e", "pyr_e", "mal__L_e", "fum_e",
           "akg_e", "ala__L_e", "acald_e"]

_G: dict = {}


def _init(namespace: str, model_paths: dict[str, str]) -> None:
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["kin"] = KineticParameters()
    _G["models"] = {c: Path(p) for c, p in model_paths.items()}


def _run_mono(code: str) -> tuple[str, dict]:
    """One monoculture -> biomass, full net secretion, importable exchange set."""
    diet = _G["diet"]
    org = CobraOrganism.from_file(str(_G["models"][code]), id=code)
    importable = sorted({list(r.metabolites)[0].id for r in org.model.exchanges})
    community = Community([org], abundances={code: 1.0}, total_biomass=TOTAL_BIOMASS)
    result = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(community, diet, _G["kin"])
    final = result.metabolites.iloc[-1].to_dict()
    secr = {m: float(final[m]) - diet.initial_concentration(m)
            for m in final
            if float(final[m]) - diet.initial_concentration(m) > SECR_FLOOR}
    biomass = float(result.biomass.iloc[-1].to_dict().get(code, 0.0))
    return code, {"biomass": biomass, "secr": secr, "importable": importable}


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
    return max(1, min(cores - 1, int((avail_mb - 600.0) / 500.0)))


def load_ckpt(path: Path) -> dict:
    done = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[rec["code"]] = rec["mono"]
    return done


def analyze(mono: dict, namespace: str, outdir: Path):
    import numpy as np
    from scipy.stats import spearmanr, pearsonr

    diet = ck.medium(namespace)
    mets_map = ck.metabolite_map(namespace)                # measured col -> exch id
    df = ck.load("results/clark2021/benchmark/MasterDF.csv")
    observed = ck.observations(df, diet=diet, metabolites=mets_map)
    mono_obs = {o.species[0]: o for o in observed if o.richness == 1}

    secr = {c: m["secr"] for c, m in mono.items()}
    imp = {c: set(m["importable"]) for c, m in mono.items()}

    def xfeed(i, j):
        """what partner j secretes that focal i can take up (mM), excl. medium."""
        return sum(amt for m, amt in secr.get(j, {}).items()
                   if m in imp.get(i, set()) and m not in diet.concentrations)

    def xfeed_ferment(i, j):
        return sum(amt for m, amt in secr.get(j, {}).items()
                   if m in imp.get(i, set()) and m in FERMENT)

    recs = [json.loads(l) for l in
            open(outdir / f"pairwise_interactions_{namespace}.jsonl")]
    pot, potf, meas, pred, resid = [], [], [], [], []
    for r in recs:
        i, j = r["focal"], r["partner"]
        if i not in secr or j not in secr:
            continue
        pot.append(xfeed(i, j)); potf.append(xfeed_ferment(i, j))
        meas.append(r["measured"]); pred.append(r["predicted"])
        resid.append(r["predicted"] - r["measured"])
    pot, potf = np.array(pot), np.array(potf)
    meas, pred, resid = np.array(meas), np.array(pred), np.array(resid)

    def corr(x, y):
        if len(x) < 3:
            return None
        sr, sp = spearmanr(x, y); pr, pp = pearsonr(x, y)
        return {"spearman": round(float(sr), 3), "spearman_p": round(float(sp), 4),
                "pearson": round(float(pr), 3), "n": int(len(x))}

    # -- decisive cross-check: does REAL cross-feeding explain REAL coexistence? --
    # Potential built from the MEASURED monoculture acid secretion (MasterDF), not
    # muODE's broken secretion. donor->acceptor = j secretes acid net>5 AND i
    # actually consumes it in monoculture (net<-2): genuine cross-feeding, from data.
    meas_secr = {c: {ex: o.net_metabolites.get(ex, 0.0) for ex in mets_map.values()}
                 for c, o in mono_obs.items()}

    def meas_xfeed(i, j):
        return sum(v for ex, v in meas_secr.get(j, {}).items()
                   if v > 5.0 and ex in imp.get(i, set()))

    def donor_acceptor(i, j):
        return sum(v for ex, v in meas_secr.get(j, {}).items()
                   if v > 5.0 and meas_secr.get(i, {}).get(ex, 0.0) < -2.0)

    mxf, dac, mmeas = [], [], []
    for r in recs:
        i, j = r["focal"], r["partner"]
        if i not in meas_secr or j not in meas_secr:
            continue
        mxf.append(meas_xfeed(i, j)); dac.append(donor_acceptor(i, j))
        mmeas.append(r["measured"])
    mxf, dac, mmeas = np.array(mxf), np.array(dac), np.array(mmeas)
    measured_crossfeed = {
        "measured_coexistence_vs_measured_xfeed_potential": corr(mxf, mmeas),
        "measured_coexistence_vs_donor_acceptor": corr(dac, mmeas),
        "n_pairs_with_genuine_crossfeed": int((dac > 0).sum()),
        "note": ("H(A) as the COMPOSITION-failure driver is refuted if real "
                 "coexistence does not rise with real cross-feeding: then the "
                 "system is competition-dominated and fixing muODE's secretion "
                 "would not fix composition."),
    }

    # under-secretion: muODE monoculture acids vs measured
    undersec = {}
    for c, o in mono_obs.items():
        if c not in secr:
            continue
        row = {}
        for col, exch in mets_map.items():
            meas_net = o.net_metabolites.get(exch)
            muode_net = secr[c].get(exch, 0.0)
            row[col] = {"muode": round(muode_net, 2),
                        "measured": round(meas_net, 2) if meas_net is not None else None}
        undersec[c] = row

    report = {
        "namespace": namespace, "n_monocultures": len(mono),
        "hypothesis_A_crossfeeding": {
            "measured_vs_xfeed_potential": corr(pot, meas),
            "measured_vs_xfeed_ferment": corr(potf, meas),
            "predicted_vs_xfeed_potential": corr(pot, pred),
            "residual(pred-meas)_vs_xfeed_potential": corr(pot, resid),
            "note": ("H(A) supported if measured rises with potential (+spearman) "
                     "while residual falls (-spearman): muODE too-competitive "
                     "exactly where cross-feeding should occur."),
        },
        "measured_crossfeeding_crosscheck": measured_crossfeed,
        "undersecretion_muode_vs_measured_mM": undersec,
        "secretion_budget": {c: round(sum(v for v in s.values()), 1)
                             for c, s in secr.items()},
    }
    (outdir / f"crossfeeding_report_{namespace}.json").write_text(
        json.dumps(report, indent=2, default=str))
    print("\n=== cross-feeding diagnostic ===")
    print(json.dumps(report["hypothesis_A_crossfeeding"], indent=2))
    print("\n--- decisive cross-check (REAL secretion) ---")
    print(json.dumps(measured_crossfeed, indent=2))
    print("\nunder-secretion (muODE vs measured, net mM) -- acetate/butyrate/lactate/succinate:")
    for c in sorted(undersec):
        cells = "  ".join(f"{col[:3]} {d['muode']}/{d['measured']}"
                          for col, d in undersec[c].items())
        print(f"  {c}: {cells}")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--analyze-only", action="store_true")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    models = load_models(args.models)
    df = ck.load("results/clark2021/benchmark/MasterDF.csv")
    observed = ck.observations(df, diet=ck.medium(args.namespace),
                               metabolites=ck.metabolite_map(args.namespace))
    monos = sorted({o.species[0] for o in observed
                    if o.richness == 1 and o.species[0] in models})
    print(f"monocultures to profile: {len(monos)}", flush=True)

    ckpt = args.outdir / f"mono_secretion_{args.namespace}.jsonl"
    done = load_ckpt(ckpt)
    todo = [c for c in monos if c not in done]

    if not args.analyze_only and todo:
        workers = args.workers or default_workers()
        model_paths = {c: str(p) for c, p in models.items()}
        print(f"done={len(done)} todo={len(todo)} workers={workers}", flush=True)
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(args.namespace, model_paths)) as ex:
            futs = {ex.submit(_run_mono, c): c for c in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                code, rec = fut.result()
                fh.write(json.dumps({"code": code, "mono": rec}) + "\n")
                fh.flush()
                done[code] = rec
                print(f"  {n}/{len(todo)}  {code}  biomass={rec['biomass']:.2f}  "
                      f"secreted {len(rec['secr'])} mets", flush=True)

    analyze(done, args.namespace, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
