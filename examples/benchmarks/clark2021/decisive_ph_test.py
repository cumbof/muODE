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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck              # noqa: E402
from muode.community import Community                      # noqa: E402
from muode.dfba import DynamicFBA                          # noqa: E402
from muode.ecology import EcologyLayer, EcologyModel       # noqa: E402
from muode.kinetics import KineticParameters              # noqa: E402
from muode.organism import CobraOrganism                  # noqa: E402
from muode.ph import WeakAcidInhibition                    # noqa: E402
from run_benchmark import T_END, load_models              # noqa: E402

DT = 0.1
ACID_EXCH = ("ac_e", "lac__L_e", "lac__D_e", "succ_e", "but_e", "ppa_e")


@dataclass
class GrowthCoupledAcid(EcologyLayer):
    """Produce each species' acids in proportion to its growth increment.

    ``yields[species][acid]`` = mM acid per unit biomass, calibrated so a
    monoculture reaches its measured endpoint acid. Tracks per-species biomass
    between steps and releases yield * dX each step (a controlled stand-in for
    the fermentation the FBA fails to route to secretion).
    """

    name: str = "growth_coupled_acid"
    yields: Dict[str, Dict[str, float]] = field(default_factory=dict)
    dt: float = DT
    _last: Dict[str, float] = field(default_factory=dict)

    def reset(self, community) -> None:
        self._last = {o.id: 0.0 for o in community.organisms}

    def metabolite_rates(self, t, M, X) -> Dict[str, float]:
        rates: Dict[str, float] = {}
        for sp, y in self.yields.items():
            if sp not in X:
                continue
            dX = X[sp] - self._last.get(sp, 0.0)
            self._last[sp] = X[sp]
            if dX <= 0:
                continue
            for acid, k in y.items():
                rates[acid] = rates.get(acid, 0.0) + k * dX / self.dt   # *dt in engine
        return rates


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


def _init(namespace, model_paths, yields):
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    _G["diet"] = ck.medium(namespace)
    _G["kin"] = KineticParameters()
    _G["models"] = {c: Path(p) for c, p in model_paths.items()}
    _G["yields"] = yields


def _ecology(codes):
    y = {c: _G["yields"].get(c, {}) for c in codes}
    return EcologyModel([GrowthCoupledAcid(yields=y, dt=DT), WeakAcidInhibition()])


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--namespace", choices=["bigg", "modelseed"], default="bigg")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--smoke", action="store_true", help="only AC, BT, AC-BT")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    models = load_models(args.models)
    mono_bio = {json.loads(l)["code"]: json.loads(l)["mono"]["biomass"]
                for l in (args.outdir / f"mono_secretion_{args.namespace}.jsonl").read_text().splitlines() if l.strip()}
    yields = measured_yields(args.namespace, mono_bio)

    if args.smoke:
        units = [["AC"], ["BT"], ["AC", "BT"]]
    else:
        obs = ck.observations(df=ck.load("results/clark2021/benchmark/MasterDF.csv"),
                              diet=ck.medium(args.namespace),
                              metabolites=ck.metabolite_map(args.namespace))
        monos = sorted({o.species[0] for o in obs if o.richness == 1 and o.species[0] in models})
        pairs = sorted({tuple(sorted(o.species)) for o in obs
                        if o.richness == 2 and all(s in models for s in o.species)})
        units = [[c] for c in monos] + [list(p) for p in pairs]

    ckpt = args.outdir / f"ph_biomass_{args.namespace}{'_smoke' if args.smoke else ''}.jsonl"
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
                initargs=(args.namespace, model_paths, yields)) as ex:
            futs = {ex.submit(_run_one, list(u)): u for u in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                codes, rec = fut.result()
                fh.write(json.dumps({"codes": sorted(codes), **rec}) + "\n"); fh.flush()
                done["-".join(sorted(codes))] = {"codes": sorted(codes), **rec}
                print(f"  {n}/{len(todo)} {'-'.join(sorted(codes))} "
                      f"bio={ {k: round(v,3) for k,v in rec['biomass'].items()} } pH={rec.get('final_pH')}",
                      flush=True)

    if args.smoke:
        a = done["AC"]["biomass"]["AC"]; b = done["BT"]["biomass"]["BT"]
        pab = done["AC-BT"]["biomass"]
        print("\n=== AC-BT smoke (baseline no-pH: AC 1.36 / BT 0.10, BT excluded) ===")
        print(f"  mono AC={a:.3f}  BT={b:.3f}")
        print(f"  pair  AC={pab['AC']:.3f}  BT={pab['BT']:.3f}  final_pH={done['AC-BT'].get('final_pH')}")
        minority = min(pab.values()) / sum(pab.values()) if sum(pab.values()) > 0 else 0
        print(f"  minority share = {minority:.3f}  (baseline 0.06, measured ~0.26)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
