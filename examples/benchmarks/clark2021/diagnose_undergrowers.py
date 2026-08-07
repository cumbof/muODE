#!/usr/bin/env python3
"""Why do the 6 strains under-grow in BOTH reconstructions? -- the triage diagnostic.

BL, BA, BH, DF, CH, DP under-grow on anaerobic DM38 in both BiGG and gapseq (they
are the residual gap after the hybrid panel: dynamic MAE 0.114 vs ceiling 0.108).
They span distinct hard metabolisms -- Bifidobacteria (bifid shunt), an acetogen
(BH, Wood-Ljungdahl), a sulfate reducer (DP), Clostridia (DF/CH). This classifies
each block before we try to fix it.

For each strain x reconstruction it runs FBA (single solves, fast) and asks:
  * g_dm38  : max growth on DM38 (medium exchanges open, others closed)
  * g_open  : max growth with EVERY importable metabolite opened
    -> g_open >> g_dm38  => AUXOTROPHY: the network works given nutrients DM38 lacks
       (then rank which single extra nutrient rescues growth most)
    -> g_open ~ g_dm38 (both low) => NETWORK GAP: an internal pathway is broken/missing
       (a reaction to add, like DF's anaerobic NAD-DHOD), not a medium problem
"""
from __future__ import annotations

import argparse
import json
import logging
import warnings
from pathlib import Path

import cobra

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

from muode.benchmarks import clark2021 as ck                       # noqa: E402

UNDERGROWERS = ["BL", "BA", "BH", "DF", "CH", "DP"]
UPTAKE = 10.0     # uniform uptake bound for the diagnostic (mmol/gDW/h)


def _set_medium(model, open_ids, bound=UPTAKE):
    for rxn in model.exchanges:
        met = list(rxn.metabolites)[0].id
        rxn.lower_bound = -bound if (open_ids is None or met in open_ids) else 0.0


def diagnose(model, diet_ids):
    _set_medium(model, diet_ids)
    g_dm38 = float(model.slim_optimize() or 0.0)
    _set_medium(model, None)                 # open everything importable
    g_open = float(model.slim_optimize() or 0.0)
    rescuers = []
    # if opening the medium helps a lot, find which single extra nutrient rescues most
    if g_open > max(1e-4, 1.5 * g_dm38):
        closed = [list(r.metabolites)[0].id for r in model.exchanges
                  if list(r.metabolites)[0].id not in diet_ids]
        for met in closed:
            _set_medium(model, set(diet_ids) | {met})
            g = float(model.slim_optimize() or 0.0)
            if g > g_dm38 + 0.01:
                rescuers.append((round(g, 3), met))
        rescuers.sort(reverse=True)
    return g_dm38, g_open, rescuers[:8]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigg", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--gapseq", type=Path, default=Path("examples/benchmarks/clark2021/gems_gapseq"))
    args = ap.parse_args()

    od = {o.species[0]: max(0.0, o.od) for o in ck.observations(
        ck.load("results/clark2021/benchmark/MasterDF.csv"),
        diet=ck.medium("bigg"), metabolites=ck.metabolite_map("bigg")) if o.richness == 1}
    bigg_diet = set(ck.medium("bigg").metabolites())
    ms_diet = set(ck.medium("modelseed").metabolites())

    out = {}
    for code in UNDERGROWERS:
        rec = {"species": ck.STRAINS[code].species, "measured_OD": round(od.get(code, 0.0), 3)}
        for ns, mdir, diet in [("bigg", args.bigg, bigg_diet), ("gapseq", args.gapseq, ms_diet)]:
            p = mdir / f"{code}.xml.gz"
            if not p.exists():
                rec[ns] = {"error": "model missing"}
                continue
            m = cobra.io.read_sbml_model(str(p))
            g_dm38, g_open, rescuers = diagnose(m, diet)
            verdict = ("NETWORK-GAP" if g_open <= max(1e-4, 1.5 * g_dm38)
                       else "AUXOTROPHY")
            rec[ns] = {"g_dm38": round(g_dm38, 4), "g_open": round(g_open, 4),
                       "verdict": verdict, "top_rescuers": rescuers}
        out[code] = rec
        b = rec.get("bigg", {})
        print(f"{code} ({rec['species']}, OD {rec['measured_OD']}): "
              f"BiGG g_dm38={b.get('g_dm38')} g_open={b.get('g_open')} -> {b.get('verdict')}",
              flush=True)
        if b.get("top_rescuers"):
            print(f"    rescuers: {b['top_rescuers']}", flush=True)

    Path("results/clark2021/benchmark/undergrowers_diagnostic.json").write_text(json.dumps(out, indent=2))
    print("\nwrote undergrowers_diagnostic.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
