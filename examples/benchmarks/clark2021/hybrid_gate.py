#!/usr/bin/env python3
"""Dynamic community gate for the hybrid gapseq/BiGG panel -- the real proof.

The hybrid MONOCULTURE panel (hybrid_panel_eval.py) beat BiGG-only on composition,
but that was a static OD-weighting proxy. This runs the actual dynamic dFBA on the
harmonized hybrid panel -- gapseq-derived Bacteroidetes (namespace-translated to
BiGG by harmonize.py, invariance-verified) + native BiGG Firmicutes, every strain
at its ATP-maintenance-calibrated yield, all sharing ONE BiGG DM38 pool -- and asks
the question the whole arc has been driving at:

  now that the yields are RIGHT, does muODE's dynamic competition reproduce
  coexistence, or is the per-step max-biomass FBA still structurally too
  winner-take-all (measured minority share 0.26 vs muODE's 0.06)?

Stage 1 (this script): 25 monos + 126 measured pairs. Scores pair-winner accuracy
(toward the 0.79 OD ceiling), minority-share median (toward 0.26), and interaction
correlation with measured -- the cheap, decisive test before the full 549-community
assembly run. Reuses the scoring in yield_calibrate.score_gate; measured data comes
from pairwise_interactions_bigg.jsonl (namespace-independent, from MasterDF).
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
import clark2021 as ck                # noqa: E402
from muode.community import Community                        # noqa: E402
from muode.dfba import DynamicFBA                            # noqa: E402
from muode.kinetics import KineticParameters                # noqa: E402
from muode.organism import CobraOrganism, MuMaxCapped  # noqa: E402
from run_benchmark import T_END                             # noqa: E402
from yield_calibrate import _maint_rxn, score_gate, measured_od, default_workers  # noqa: E402

BACT_GENERA = ("Prevotella", "Parabacteroides", "Phocaeicola", "Bacteroides")
_G: dict = {}


def _init(paths, atpm, capacity=None, mumax=None, lp_timeout=None, pfba=False, mumax_map=None):
    try:
        import swiglpk
        swiglpk.glp_term_out(swiglpk.GLP_OFF)
    except Exception:
        pass
    import logging
    logging.disable(logging.CRITICAL)
    _G["diet"] = ck.medium("bigg")          # everything speaks BiGG post-harmonization
    _G["kin"] = KineticParameters()
    _G["paths"] = paths
    _G["atpm"] = atpm
    _G["capacity"] = capacity or {}         # per-strain carrying capacity (= measured OD)
    _G["mumax"] = mumax                      # uniform physiological max growth-rate ceiling (1/h)
    _G["mumax_map"] = mumax_map or {}        # per-strain mu_max (e.g. genome-derived, gRodon)
    _G["lp_timeout"] = lp_timeout            # per-LP wall-clock cap (s); guards degenerate LPs
    _G["pfba"] = pfba                        # parsimonious FBA -> deterministic fluxes


def _organism(code):
    model = cobra.io.read_sbml_model(str(_G["paths"][code]))
    rid = _maint_rxn(model)
    if rid:
        model.reactions.get_by_id(rid).lower_bound = float(_G["atpm"].get(code, 0.0))
    org = CobraOrganism(model, id=code, parsimonious=_G.get("pfba", False))
    # physiological max-growth-rate cap, applied POST-solve (see MuMaxCapped) so a
    # fast grower cannot run away and monopolize shared substrate -- without the
    # degenerate-LP blowup that a biomass-ub cap causes. Per-strain (genome-derived)
    # value overrides the uniform one when a mumax_map is supplied.
    mm = _G.get("mumax_map", {}).get(code, _G.get("mumax"))
    if mm is not None:
        org = MuMaxCapped(org, mm)
    return org


def _ecology(codes):
    """Density-dependent self-limitation keyed to measured OD, if capacities given."""
    if not _G.get("capacity"):
        return None
    from muode.density import LogisticCarryingCapacity
    from muode.ecology import EcologyModel
    cap = {c: _G["capacity"][c] for c in codes if c in _G["capacity"]}
    if not cap:
        return None
    return EcologyModel([LogisticCarryingCapacity(capacity=cap)])


def _gate_one(codes):
    orgs = [_organism(c) for c in codes]
    comm = Community(orgs, abundances={c: 1.0 / len(codes) for c in codes}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(
        comm, _G["diet"], _G["kin"], ecology=_ecology(codes))
    return codes, {c: float(v) for c, v in res.biomass.iloc[-1].to_dict().items()}


def build_panel(bigg_dir: Path, harm_dir: Path):
    """Per-strain best-reconstruction path + merged calibrated ATPM."""
    bact = {c for c, s in ck.STRAINS.items() if s.species.startswith(BACT_GENERA)}
    paths, source = {}, {}
    for c in ck.STRAINS:
        if c in bact and (harm_dir / f"{c}.xml.gz").exists():
            paths[c], source[c] = harm_dir / f"{c}.xml.gz", "gapseq"
        elif (bigg_dir / f"{c}.xml.gz").exists():
            paths[c], source[c] = bigg_dir / f"{c}.xml.gz", "bigg"
    outdir = Path("results/clark2021/benchmark")
    gcal = {json.loads(l)["code"]: json.loads(l)["atpm"]
            for l in (outdir / "atpm_calibration_modelseed.jsonl").read_text().splitlines() if l.strip()}
    bcal = {json.loads(l)["code"]: json.loads(l)["atpm"]
            for l in (outdir / "atpm_calibration_bigg.jsonl").read_text().splitlines() if l.strip()}
    atpm = {c: (gcal.get(c, 0.0) if source[c] == "gapseq" else bcal.get(c, 0.0)) for c in paths}
    return paths, source, atpm


def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bigg-models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--harm-models", type=Path, default=Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    ap.add_argument("--outdir", type=Path, default=Path("results/clark2021/benchmark"))
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--logistic", action="store_true",
                    help="add density-dependent self-limitation, capacity = measured OD")
    ap.add_argument("--mumax", type=float, default=None,
                    help="uniform physiological max growth-rate ceiling (1/h) on the biomass reaction")
    ap.add_argument("--mumax-file", type=Path, default=None,
                    help="JSON {code: mumax} per-strain ceilings (e.g. genome-derived gRodon)")
    ap.add_argument("--mumax-scale", type=float, default=1.0,
                    help="multiply every per-strain mumax by this (rescale absolute level)")
    ap.add_argument("--lp-timeout", type=float, default=None,
                    help="per-LP wall-clock cap (s); use with --mumax to guard degenerate LPs")
    ap.add_argument("--pfba", action="store_true",
                    help="parsimonious FBA per step -> deterministic fluxes across platforms")
    ap.add_argument("--only", default=None, help="comma-separated pair subset e.g. BO-ER,AC-BT (debug)")
    ap.add_argument("--uniform-atpm", type=float, default=None,
                    help="override per-strain calibrated ATPM with one literature value "
                         "(mmol/gDW/h) -> DATA-FREE yields (no OD calibration)")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    paths, source, atpm = build_panel(args.bigg_models, args.harm_models)
    if args.uniform_atpm is not None:                    # data-free: drop the OD calibration
        atpm = {c: args.uniform_atpm for c in paths}
    workers = args.workers or default_workers()

    obs = ck.observations(ck.load(str(args.outdir / "MasterDF.csv")),
                          diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg"))
    monos = sorted({o.species[0] for o in obs if o.richness == 1 and o.species[0] in paths})
    pairs = sorted({tuple(sorted(o.species)) for o in obs
                    if o.richness == 2 and all(s in paths for s in o.species)})
    units = [[c] for c in monos] + [list(p) for p in pairs]

    # carrying capacity K_i = measured monoculture OD (only when --logistic)
    capacity = None
    tag = ""
    if args.logistic:
        capacity = {o.species[0]: max(0.0, o.od) for o in obs
                    if o.richness == 1 and o.species[0] in paths}
        tag = "_logistic"
        print(f"logistic layer ON: {len(capacity)} OD carrying capacities", flush=True)
    mumax_map = {}
    if args.mumax_file:
        raw = json.loads(args.mumax_file.read_text())
        mumax_map = {c: float(v) * args.mumax_scale for c, v in raw.items()}
        tag = f"_mumaxfile{('x%g' % args.mumax_scale) if args.mumax_scale != 1 else ''}"
        import statistics as _st
        print(f"per-strain mu_max ON ({args.mumax_file.name}, x{args.mumax_scale}): "
              f"{len(mumax_map)} strains, median {_st.median(mumax_map.values()):.3f}/h", flush=True)
    elif args.mumax is not None:
        tag = f"_mumax{args.mumax:g}"
        print(f"mu_max cap ON: {args.mumax}/h ceiling on the biomass reaction", flush=True)
    if args.uniform_atpm is not None:
        tag += f"_uatpm{args.uniform_atpm:g}"
        print(f"DATA-FREE: uniform ATPM {args.uniform_atpm} for all strains "
              "(no OD calibration)", flush=True)
    if args.pfba:
        tag += "_pfba"
        print("pFBA ON: parsimonious flux (deterministic across platforms)", flush=True)

    if args.only:                             # debug: restrict to named pairs/monos
        want = {tuple(sorted(p.split("-"))) for p in args.only.split(",")}
        units = [u for u in units if tuple(sorted(u)) in want]
        tag += "_only"

    ckpt = args.outdir / f"hybrid_gate_biomass{tag}.jsonl"
    done = {}
    if ckpt.exists():
        for l in ckpt.read_text().splitlines():
            if l.strip():
                r = json.loads(l); done["-".join(sorted(r["codes"]))] = r
    todo = [u for u in units if "-".join(sorted(u)) not in done]
    print(f"hybrid gate: gapseq={sum(v=='gapseq' for v in source.values())} "
          f"bigg={sum(v=='bigg' for v in source.values())} | "
          f"units={len(units)} done={len(done)} todo={len(todo)} workers={workers}", flush=True)

    if todo:
        with open(ckpt, "a") as fh, ProcessPoolExecutor(
                max_workers=workers, initializer=_init,
                initargs=(paths, atpm, capacity, args.mumax, args.lp_timeout, args.pfba, mumax_map)) as ex:
            futs = {ex.submit(_gate_one, list(u)): u for u in todo}
            for k, fut in enumerate(as_completed(futs), 1):
                codes, bio = fut.result()
                rec = {"codes": sorted(codes), "biomass": bio}
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                done["-".join(sorted(codes))] = rec
                if k % 10 == 0 or k == len(todo):
                    print(f"  {k}/{len(todo)}", flush=True)

    if args.only:                             # debug subset: just report endpoints, no scoring
        for u in units:
            k = "-".join(sorted(u))
            print(f"  {k}: {done.get(k, {}).get('biomass')}", flush=True)
        return 0
    # score with the bigg measured-interaction reference (namespace-independent data)
    rep = score_gate(done, "bigg", args.outdir)
    (args.outdir / f"hybrid_gate_report{tag}.json").write_text(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
