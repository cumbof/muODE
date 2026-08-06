#!/usr/bin/env python3
"""Does moderating the winner's growth-rate advantage restore coexistence?

Mechanism-agnostic, cheap test that de-risks the whole minority-survival question.
The diagnostic localized the over-exclusion to the competition DYNAMICS; enzyme-pool
turned out blocked/mismatched for this substrate-limited regime. Underneath both
candidates (enzyme-pool and substrate-allocation) is one question: if we compress
the growth-rate distribution so the fastest grower loses its runaway edge, does the
minority survive?

This caps every organism's per-step growth at an ABSOLUTE ceiling mu_cap (biomass
reaction upper bound; folded into the baseline so it survives the dFBA's per-step
reset). A cap that binds the winner but not the slower minority compresses the
advantage without a uniform slowdown (which would cancel in relative share). It
sweeps mu_cap on the 6 pure-dynamics over-excluded pairs and asks whether minority
share climbs toward the measured ~0.26.

  climbs  -> rate-moderation is the lever; an enzyme/ec implementation is worth the
             cost, and this cap is a fast stand-in to tune against.
  flat    -> the winner-take-all is NOT about rate; it's substrate ALLOCATION
             (dfba.py independent-greedy draw) -> pivot there.

No LP bloat (one bound per model), so it is fast and definitely binding.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import cobra

warnings.filterwarnings("ignore")
import logging
logging.disable(logging.CRITICAL)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from muode.benchmarks import clark2021 as ck                # noqa: E402
from muode.community import Community                        # noqa: E402
from muode.dfba import DynamicFBA                            # noqa: E402
from muode.kinetics import KineticParameters                # noqa: E402
from muode.organism import CobraOrganism                    # noqa: E402
from run_benchmark import T_END                             # noqa: E402
from hybrid_gate import build_panel, _maint_rxn             # noqa: E402

SUBSET = ["CA-ER", "AC-CA", "AC-BP", "ER-PC", "BO-ER", "AC-BT"]


def _biomass_rxns(model):
    return [r for r in model.reactions if getattr(r, "objective_coefficient", 0.0) != 0.0]


def make_org(code, paths, atpm, mu_cap):
    model = cobra.io.read_sbml_model(str(paths[code]))
    rid = _maint_rxn(model)
    if rid:
        model.reactions.get_by_id(rid).lower_bound = float(atpm.get(code, 0.0))
    if mu_cap is not None:
        for r in _biomass_rxns(model):
            if r.upper_bound > mu_cap:
                r.upper_bound = float(mu_cap)   # snapshotted into _base_bounds by CobraOrganism
    return CobraOrganism(model, id=code)


def pair_biomass(codes, paths, atpm, mu_cap, kin):
    orgs = [make_org(c, paths, atpm, mu_cap) for c in codes]
    comm = Community(orgs, abundances={c: 1.0 / len(codes) for c in codes}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, ck.medium("bigg"), kin)
    bio = res.biomass.iloc[-1].to_dict()
    return {c: float(bio.get(c, 0.0)) for c in codes}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--caps", default="none,0.5,0.3,0.2,0.1",
                    help="mu_cap sweep; 'none' = uncapped control")
    ap.add_argument("--subset", default=",".join(SUBSET))
    args = ap.parse_args()
    caps = [None if c == "none" else float(c) for c in args.caps.split(",")]
    pairs = [p.split("-") for p in args.subset.split(",")]

    paths, source, atpm = build_panel(Path("examples/benchmarks/clark2021/models"),
                                      Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    kin = KineticParameters()
    byp = {}
    for l in (Path("results/clark2021/benchmark/pairwise_interactions_bigg.jsonl")
              ).read_text().splitlines():
        r = json.loads(l); byp.setdefault(tuple(sorted((r["focal"], r["partner"]))), {})[r["focal"]] = r

    out = {}
    hdr = " ".join(f"cap={('none' if c is None else c)!s:<5}" for c in caps)
    print(f"{'pair':<8} {'meas':>5}  {hdr}")
    for a, b in pairs:
        d = byp[tuple(sorted((a, b)))]
        mpa, mpb = d[a]["meas_pair"], d[b]["meas_pair"]
        meas_min = min(mpa, mpb) / (mpa + mpb)
        row = {"measured_min": round(meas_min, 3), "by_cap": {}}
        cells = []
        for c in caps:
            bio = pair_biomass([a, b], paths, atpm, c, kin)
            tot = sum(bio.values()) or 1.0
            mn = min(bio.values()) / tot
            key = "none" if c is None else c
            row["by_cap"][str(key)] = {"biomass": {k: round(v, 3) for k, v in bio.items()},
                                       "minority_share": round(mn, 3)}
            cells.append(f"{mn:<9.3f}")
        out["-".join(sorted((a, b)))] = row
        print(f"{a}-{b:<6} {meas_min:>5.2f}  " + " ".join(cells), flush=True)

    Path("results/clark2021/benchmark/growth_cap_probe.json").write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
