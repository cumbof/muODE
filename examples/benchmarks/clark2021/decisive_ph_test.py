#!/usr/bin/env python3
"""Decisive test: does pH negative-feedback break muODE's winner-take-all?

The go/no-go confirmed acidification is a real, growth-rank-independent driver of
the MEASURED pairwise interactions (partner acidification vs measured interaction
spearman -0.335; R^2 0.003 -> 0.084 added over growth rank). muODE has the exact
mechanism (muode.ph.WeakAcidInhibition: SCFA -> pH -> undissociated-acid growth
inhibition) but it is INERT because the FBA secretes ~0 acid.

This is the controlled experiment: give muODE realistic acids and turn the pH
layer on, then re-run the pairwise gate. Acids are produced by a source layer
coupled to each species' GROWTH at that species' MEASURED monoculture yield
(mM acid per unit biomass), so acid accumulates realistically *without* relying
on the broken FBA secretion -- and without hand-tuning to the pairwise outcome.
WeakAcidInhibition then acidifies the medium and inhibits growth; because the
dominant grower is also the biggest acidifier, it caps itself and (via the shared
substrate pool) frees resources for the slower partner.

Baseline (no ecology, commit 1912d24): faster monoculture-grower wins 94% of
pairs, minority member's share median 0.06. The pH mechanism WORKS for
composition iff those move toward the measured 52% / 0.26.

NOT a fit: acid yields come from monoculture data; Ki is the literature default,
uniform. If uniform Ki already helps, per-species acid tolerance (from
independent phylogeny/literature, not Clark) would help more -- a later step.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))
import clark2021 as ck              # noqa: E402
from muode.community import Community                      # noqa: E402
from muode.dfba import DynamicFBA                          # noqa: E402
from muode.ecology import EcologyModel                    # noqa: E402
from muode.kinetics import KineticParameters              # noqa: E402
from muode.organism import CobraOrganism                  # noqa: E402
from muode.ph import WeakAcidInhibition                    # noqa: E402
from muode.secretion import GrowthCoupledSecretion        # noqa: E402
from run_benchmark import T_END, load_models              # noqa: E402

DT = 0.1
ACID_EXCH = ("ac_e", "lac__L_e", "lac__D_e", "succ_e", "but_e", "ppa_e")

#: Per-species acid tolerance (undissociated-acid Ki, mmol/L) assigned in 3 tiers
#: PURELY from independent physiology -- Duncan et al. 2009 Environ. Microbiol.
#: (growth of colonic bacteria at pH 5.5/6.2/6.7): all Bacteroidetes grow poorly
#: at pH 5.5 (sensitive); butyrate-producing Firmicutes + Bifidobacterium dominate
#: at pH 5.5 (tolerant). Set once by taxon; NOT tuned to any Clark pairwise outcome.
KI_SENSITIVE, KI_INTERMEDIATE, KI_TOLERANT = 3.0, 10.0, 30.0
ACID_TOLERANCE = {
    # Bacteroidetes -- acid-sensitive (Duncan: Bacteroides poor at pH 5.5)
    "BT": KI_SENSITIVE, "BU": KI_SENSITIVE, "BV": KI_SENSITIVE, "BC": KI_SENSITIVE,
    "BF": KI_SENSITIVE, "BO": KI_SENSITIVE, "BY": KI_SENSITIVE, "PJ": KI_SENSITIVE,
    "PC": KI_SENSITIVE,
    "DP": KI_SENSITIVE,   # Desulfovibrio piger -- neutralophilic sulfate reducer
    # butyrate-producing Firmicutes + Bifidobacterium -- acid-tolerant (Duncan)
    "ER": KI_TOLERANT, "RI": KI_TOLERANT, "FP": KI_TOLERANT, "AC": KI_TOLERANT,
    "CC": KI_TOLERANT, "BL": KI_TOLERANT, "BA": KI_TOLERANT, "BP": KI_TOLERANT,
    "CA": KI_TOLERANT,    # Collinsella aerofaciens -- acid-tolerant Actinobacterium
    # remaining Clostridia/Lachnospiraceae/Actinobacteria -- intermediate
    "DL": KI_INTERMEDIATE, "DF": KI_INTERMEDIATE, "BH": KI_INTERMEDIATE,
    "CG": KI_INTERMEDIATE, "CH": KI_INTERMEDIATE, "EL": KI_INTERMEDIATE,
    "HB": KI_INTERMEDIATE,
}


def measured_yields(namespace: str, mono_biomass: Mapping[str, float]):
    """{species -> {acid_exch -> mM per unit biomass}} from measured monoculture."""
    diet = ck.medium(namespace)
    mmap = ck.metabolite_map(namespace)
    df = ck.load("results/clark2021/benchmark/MasterDF.csv")
    obs = ck.observations(df, diet=diet, metabolites=mmap)
    mono = {o.species[0]: o for o in obs if o.richness == 1}
    acid_exch = [e for e in ACID_EXCH if e in mmap.values()]
    out = {}
    for c, o in mono.items():
        b = mono_biomass.get(c, 0.0)
        if b <= 0:
            continue
        out[c] = {e: max(0.0, o.net_metabolites.get(e, 0.0)) / b
                  for e in acid_exch if o.net_metabolites.get(e, 0.0) > 0.5}
    return out


_G: dict = {}


def _init(namespace, model_paths, yields, per_species_ki):
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["kin"] = KineticParameters()
    _G["models"] = {c: Path(p) for c, p in model_paths.items()}
    _G["yields"] = yields
    _G["per_species_ki"] = per_species_ki


def _ecology(codes):
    y = {c: _G["yields"].get(c, {}) for c in codes}
    if _G.get("per_species_ki"):
        ki = {c: ACID_TOLERANCE[c] for c in codes if c in ACID_TOLERANCE}
        acid_layer = WeakAcidInhibition(ki=ki)
    else:
        acid_layer = WeakAcidInhibition()
    return EcologyModel([GrowthCoupledSecretion(yields=y, dt=DT), acid_layer])


def sim_biomass_ph(codes):
    organisms = [CobraOrganism.from_file(str(_G["models"][c]), id=c) for c in codes]
    community = Community(organisms, abundances={c: 1.0 / len(codes) for c in codes},
                          total_biomass=0.01)
    result = DynamicFBA(t_end=T_END, dt=DT, n_jobs=1).run(
        community, _G["diet"], _G["kin"], ecology=_ecology(codes))
    bio = {c: float(v) for c, v in result.biomass.iloc[-1].to_dict().items()}
    pH = None
    env = result.environment
    if env is not None and "pH" in getattr(env, "columns", []):
        pH = float(env["pH"].iloc[-1])
    return codes, {"biomass": bio, "final_pH": pH}


def _run_one(codes):
    return sim_biomass_ph(codes)


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


def analyze_full(done, namespace, outdir, per_species_ki):
    """Aggregate verdict: does the pH layer improve the pairwise gate vs baseline?"""
    import math
    import numpy as np
    from scipy.stats import spearmanr, pearsonr

    base = [json.loads(l) for l in
            open(outdir / f"pairwise_interactions_{namespace}.jsonl")]
    meas, basep, newp = [], [], []
    fast_wins_new, fast_wins_meas, minority_new, minority_meas = [], [], [], []
    by_pair = {}
    for r in base:
        by_pair.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r
    for r in base:
        i, j = r["focal"], r["partner"]
        ka, kpair = "-".join([i]), "-".join(sorted((i, j)))
        if ka not in done or kpair not in done:
            continue
        bi_alone = done[ka]["biomass"].get(i, 0.0)
        bi_pair = done[kpair]["biomass"].get(i, 0.0)
        if bi_alone <= 0:
            continue
        meas.append(r["measured"]); basep.append(r["predicted"])
        newp.append(math.log(max(bi_pair, 1e-9) / bi_alone) + math.log(2))
    for (a, b), d in by_pair.items():
        kp = "-".join(sorted((a, b)))
        if kp not in done or a not in done[kp]["biomass"] or b not in done[kp]["biomass"]:
            continue
        pa, pb = done[kp]["biomass"][a], done[kp]["biomass"][b]
        if a not in d or b not in d:
            continue
        mono_a = done.get("-".join([a]), {}).get("biomass", {}).get(a, 0)
        mono_b = done.get("-".join([b]), {}).get("biomass", {}).get(b, 0)
        fast = a if mono_a > mono_b else b
        if pa + pb > 0:
            new_win = a if pa > pb else b
            fast_wins_new.append(fast == new_win)
            minority_new.append(min(pa, pb) / (pa + pb))
        mp_a, mp_b = d[a]["meas_pair"], d[b]["meas_pair"]
        if mp_a + mp_b > 0:
            meas_win = a if mp_a > mp_b else b
            fast_wins_meas.append(fast == meas_win)
            minority_meas.append(min(mp_a, mp_b) / (mp_a + mp_b))
    meas, basep, newp = np.array(meas), np.array(basep), np.array(newp)

    def cc(x, y):
        sr, _ = spearmanr(x, y); pr, _ = pearsonr(x, y)
        return {"spearman": round(float(sr), 3), "pearson": round(float(pr), 3)}

    report = {
        "namespace": namespace, "per_species_ki": per_species_ki, "n_interactions": int(len(meas)),
        "correlation_with_measured": {
            "baseline_dFBA": cc(basep, meas), "pH_muODE": cc(newp, meas)},
        "faster_grower_wins_pair": {
            "baseline": 0.942, "pH_muODE": round(float(np.mean(fast_wins_new)), 3),
            "measured": round(float(np.mean(fast_wins_meas)), 3)},
        "minority_share_median": {
            "baseline": 0.061, "pH_muODE": round(float(np.median(minority_new)), 3),
            "measured": round(float(np.median(minority_meas)), 3)},
    }
    tag = "phki" if per_species_ki else "phuniform"
    (outdir / f"ph_gate_report_{namespace}_{tag}.json").write_text(json.dumps(report, indent=2))
    print("\n=== pH GATE aggregate verdict ===")
    print(json.dumps(report, indent=2))
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--smoke", action="store_true", help="only AC, BT, AC-BT")
    ap.add_argument("--pairs", type=str, default=None,
                    help="custom smoke: comma-sep pairs e.g. BT-BL,BU-BA (runs their monos too)")
    ap.add_argument("--per-species-ki", action="store_true",
                    help="use Duncan-2009 per-taxon acid tolerance (else uniform Ki)")
    ap.add_argument("--tag", type=str, default="", help="checkpoint suffix")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    models = load_models(args.models)
    mono_bio = {json.loads(l)["code"]: json.loads(l)["mono"]["biomass"]
                for l in (args.outdir / f"mono_secretion_{args.namespace}.jsonl").read_text().splitlines() if l.strip()}
    yields = measured_yields(args.namespace, mono_bio)

    if args.pairs:
        pairset = [tuple(p.split("-")) for p in args.pairs.split(",")]
        monoset = sorted({s for p in pairset for s in p})
        units = [[c] for c in monoset] + [list(p) for p in pairset]
    elif args.smoke:
        units = [["AC"], ["BT"], ["AC", "BT"]]
    else:
        obs = ck.observations(df=ck.load("results/clark2021/benchmark/MasterDF.csv"),
                              diet=ck.medium(args.namespace),
                              metabolites=ck.metabolite_map(args.namespace))
        monos = sorted({o.species[0] for o in obs if o.richness == 1 and o.species[0] in models})
        pairs = sorted({tuple(sorted(o.species)) for o in obs
                        if o.richness == 2 and all(s in models for s in o.species)})
        units = [[c] for c in monos] + [list(p) for p in pairs]

    suffix = args.tag or ("_smoke" if args.smoke else ("_pairs" if args.pairs else ""))
    ckpt = args.outdir / f"ph_biomass_{args.namespace}{suffix}.jsonl"
    done = {}
    if ckpt.exists():
        for l in ckpt.read_text().splitlines():
            if l.strip():
                r = json.loads(l); done["-".join(sorted(r["codes"]))] = r
    todo = [u for u in units if "-".join(sorted(u)) not in done]

    workers = args.workers or default_workers()
    model_paths = {c: str(p) for c, p in models.items()}
    print(f"units={len(units)} done={len(done)} todo={len(todo)} workers={workers}", flush=True)
    if todo:
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(args.namespace, model_paths, yields, args.per_species_ki)) as ex:
            futs = {ex.submit(_run_one, list(u)): u for u in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                codes, rec = fut.result()
                fh.write(json.dumps({"codes": sorted(codes), **rec}) + "\n"); fh.flush()
                done["-".join(sorted(codes))] = {"codes": sorted(codes), **rec}
                print(f"  {n}/{len(todo)} {'-'.join(sorted(codes))} "
                      f"bio={ {k: round(v,3) for k,v in rec['biomass'].items()} } pH={rec.get('final_pH')}",
                      flush=True)

    if args.pairs:
        import math
        base = {}
        for l in open(args.outdir / f"pairwise_interactions_{args.namespace}.jsonl"):
            r = json.loads(l); base[(r["focal"], r["partner"])] = (r["measured"], r["predicted"])
        print(f"\n=== per-species-Ki={args.per_species_ki} : interaction(i<-j) ===")
        print("  pair        measured  baseline-muODE  pH-muODE   final_pH")
        for p in pairset:
            i, j = p
            pk = "-".join(sorted(p))
            if pk not in done or "-".join([i]) not in done or "-".join([j]) not in done:
                continue
            bi_alone = done["-".join([i])]["biomass"].get(i, 0)
            bi_pair = done[pk]["biomass"].get(i, 0)
            newint = (math.log(max(bi_pair, 1e-9) / bi_alone) + math.log(2)) if bi_alone > 0 else float("nan")
            meas, basepred = base.get((i, j), (float("nan"), float("nan")))
            print(f"  {i}<-{j:5}  {meas:+7.2f}   {basepred:+7.2f}       {newint:+7.2f}    {done[pk].get('final_pH')}")
        return 0

    if args.smoke:
        a = done["AC"]["biomass"]["AC"]; b = done["BT"]["biomass"]["BT"]
        pab = done["AC-BT"]["biomass"]
        print("\n=== AC-BT smoke (baseline no-pH: AC 1.36 / BT 0.10, BT excluded) ===")
        print(f"  mono AC={a:.3f}  BT={b:.3f}")
        print(f"  pair  AC={pab['AC']:.3f}  BT={pab['BT']:.3f}  final_pH={done['AC-BT'].get('final_pH')}")
        minority = min(pab.values()) / sum(pab.values()) if sum(pab.values()) > 0 else 0
        print(f"  minority share = {minority:.3f}  (baseline 0.06, measured ~0.26)")
        return 0

    analyze_full(done, args.namespace, args.outdir, args.per_species_ki)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
