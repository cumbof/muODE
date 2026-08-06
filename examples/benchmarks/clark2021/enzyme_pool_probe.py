#!/usr/bin/env python3
"""Does a proteome-budget (rate-yield trade-off) lift the over-excluded minority?

The minority diagnostic localized the coexistence gap to the dFBA COMPETITION
DYNAMICS (not yields): reality is ~yield-proportional (minority 0.35 OD-prop, 0.26
measured) but the per-step max-biomass dFBA crushes it to 0.13. Enzyme-pool is the
indicated lever -- a proteome budget imposes a rate-yield trade-off that caps the
fastest grower's advantage, the classic coexistence mechanism.

The package's enzyme-pool (muode.enzyme.apply_protein_pool_constraint) is present
but UNPARAMETERIZED (0 kcats), so this is the directional proxy I flagged: uniform
kcat + MW reduce the proteome budget to a per-organism TOTAL-INTERNAL-FLUX cap,

    sum_r |v_r|  <=  f * Sigma_ref ,   Sigma_ref = pFBA min-flux at max growth,

with one tightness knob f (f=1 ~ baseline; f<1 forces the cell below its greedy
flux optimum -> lower max growth -> compressed growth-rate advantage). It sweeps f
on the pure-dynamics over-excluded pairs (both strains grow, no under-growers) and
asks the single question that decides whether to run the full gate: does minority
share rise monotonically toward the measured ~0.26 as f tightens?
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import cobra
from cobra.flux_analysis import pfba

warnings.filterwarnings("ignore")
import logging
logging.disable(logging.CRITICAL)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck                # noqa: E402
from muode.community import Community                        # noqa: E402
from muode.dfba import DynamicFBA                            # noqa: E402
from muode.enzyme import _is_internal, apply_protein_pool_constraint  # noqa: E402
from muode.kinetics import KineticParameters                # noqa: E402
from muode.organism import CobraOrganism                    # noqa: E402
from run_benchmark import T_END                             # noqa: E402
from hybrid_gate import build_panel, _maint_rxn             # noqa: E402

SUBSET = ["CA-ER", "AC-CA", "AC-BP", "ER-PC", "BO-ER", "AC-BT"]
_REF: dict = {}


def sigma_ref(model) -> float:
    """Total abs internal flux at (parsimonious) max growth -- the flux budget base."""
    sol = pfba(model)
    return float(sum(abs(sol.fluxes[r.id]) for r in model.reactions if _is_internal(r)))


def make_org(code, paths, atpm, f):
    model = cobra.io.read_sbml_model(str(paths[code]))
    rid = _maint_rxn(model)
    if rid:
        model.reactions.get_by_id(rid).lower_bound = float(atpm.get(code, 0.0))
    if f is not None and f < 1.0 + 1e-9:
        ref = _REF.get(code)
        if ref is None:
            ref = _REF[code] = sigma_ref(model)
        # uniform kcat=1/s + MW=3600 kDa => pool coeff = 1 => constraint is sum|v| <= budget
        kin = KineticParameters()
        for r in model.reactions:
            if _is_internal(r):
                kin.set_kcat(code, r.id, 1.0)
        apply_protein_pool_constraint(model, kin, code, pool_budget=f * ref, default_mw=3600.0)
    return CobraOrganism(model, id=code)


def pair_biomass(codes, paths, atpm, f, kin_run):
    orgs = [make_org(c, paths, atpm, f) for c in codes]
    comm = Community(orgs, abundances={c: 1.0 / len(codes) for c in codes}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, ck.medium("bigg"), kin_run)
    bio = res.biomass.iloc[-1].to_dict()
    return {c: float(bio.get(c, 0.0)) for c in codes}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fs", default="1.0,0.7,0.5,0.3")
    ap.add_argument("--subset", default=",".join(SUBSET))
    args = ap.parse_args()
    fs = [float(x) for x in args.fs.split(",")]
    pairs = [p.split("-") for p in args.subset.split(",")]

    paths, source, atpm = build_panel(Path("examples/benchmarks/clark2021/models"),
                                      Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    kin_run = KineticParameters()

    # measured minority per pair for reference
    byp = {}
    for l in (Path("results/clark2021/benchmark/pairwise_interactions_bigg.jsonl")
              ).read_text().splitlines():
        r = json.loads(l); byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r

    out = {}
    print(f"{'pair':<8} {'meas_min':>8} " + " ".join(f"f={f:<4}" for f in fs))
    for a, b in pairs:
        d = byp[tuple(sorted((a, b)))]
        mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
        meas_min = min(mpa, mpb) / (mpa + mpb)
        row = {"measured_min": round(meas_min, 3), "by_f": {}}
        cells = []
        for f in fs:
            bio = pair_biomass([a, b], paths, atpm, f, kin_run)
            tot = sum(bio.values()) or 1.0
            mn = min(bio.values()) / tot
            row["by_f"][f] = {"biomass": {k: round(v, 3) for k, v in bio.items()},
                              "minority_share": round(mn, 3)}
            cells.append(f"{mn:<6.3f}")
        out["-".join(sorted((a, b)))] = row
        print(f"{a}-{b:<6} {meas_min:>8.3f} " + " ".join(cells), flush=True)

    Path("results/clark2021/benchmark/enzyme_pool_probe.json").write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
