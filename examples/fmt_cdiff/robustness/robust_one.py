#!/usr/bin/env python3
"""Run ONE robustness config of the fmt_cdiff scenario and dump a small JSON.

One config per process (the sweep_vmax_bai.py pattern), so a shell `xargs -P` fans the
whole design out across cores with no multiprocessing/solver-state hazards.

A config is a JSON file with any of:
  label (required), fmt(bool), ablate(str), donors(list of member slugs to SEED via FMT),
  dilution_rate, vmax_bai, vmax_bsh, fmt_time, fmt_biomass, t_end, dt.

Defaults reproduce the PUBLISHED decomposition operating point (library vmax_bai, so the
baseline lambda/C.diff match fmt_decomposition_result.json: lam~0.0087, C.diff~0.022 with
bile; lam~0.0119, C.diff~0.039 without). `donors` lets us drop the keystone: seeding every
donor EXCEPT C_scindens_ATCC35704 is the Buffie-2015 bad case -- an FMT that engrafts but
does not restore secondary-bile-acid colonisation resistance.
"""
import sys, json, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path
HERE = Path(__file__).resolve().parent
EX = HERE.parent
sys.path.insert(0, str(EX))
import numpy as np
import scenario as gs
from muode.dfba import DynamicFBA
from muode.inject import Injection
from muode.kinetics import KineticParameters
from muode.lifecycle import SporeForming


def rebound(series, window_h=30.0):
    t = np.asarray(series.index, float); x = np.asarray(series.values, float)
    m = (t >= t[-1] - window_h) & (x > 0)
    return float(np.polyfit(t[m], np.log(x[m]), 1)[0]) if m.sum() >= 2 else float("nan")


def run_cfg(cfg):
    donors = cfg.get("donors", list(gs.DONORS))
    ablate = cfg.get("ablate", "")
    D = cfg.get("dilution_rate", gs.DILUTION_RATE)
    vmax_bai = cfg.get("vmax_bai", None)        # None -> library default (published point)
    vmax_bsh = cfg.get("vmax_bsh", None)
    fmt_time = cfg.get("fmt_time", 12.0)
    fmt_biomass = cfg.get("fmt_biomass", 0.05)
    t_end = cfg.get("t_end", 120.0); dt = cfg.get("dt", 0.05)
    fmt = cfg.get("fmt", True)

    community = gs.build_community(
        abundances={gs.PATHOGEN: 1.0, **{d: 0.0 for d in gs.DONORS}})
    injections = None
    if fmt and donors:
        injections = [Injection.from_abundances(
            fmt_time, {d: 1.0 / len(donors) for d in donors},
            total_biomass=fmt_biomass, name="FMT")]
    ecology = gs.cdi_ecology(ablate, vmax_bai=vmax_bai, vmax_bsh=vmax_bsh)
    result = DynamicFBA(t_end=t_end, dt=dt, dilution_rate=D).run(
        community, gs.cdi_diet(dilution_rate=D), KineticParameters(),
        injections=injections, ecology=ecology)
    for layer in ecology.layers:
        if isinstance(layer, SporeForming):
            result.meta["spore_latched"] = layer.latched()

    B = result.biomass; M = result.metabolites
    veg = float(B[gs.PATHOGEN].iloc[-1])
    spores = 0.0
    if getattr(result, "spores", None) is not None and gs.PATHOGEN in result.spores:
        spores = float(result.spores[gs.PATHOGEN].iloc[-1])
    return {
        "label": cfg["label"],
        "lam": round(rebound(B[gs.PATHOGEN]), 5),
        "cdiff_veg": round(veg, 5),
        "cdiff_spores": round(spores, 5),
        "cdiff_total": round(veg + spores, 5),
        "dca_final": round(float(M["dca_e"].iloc[-1]), 4) if "dca_e" in M else None,
        "ca_final": round(float(M["ca_e"].iloc[-1]), 4) if "ca_e" in M else None,
        "scindens_seeded": ("C_scindens_ATCC35704" in donors) and fmt,
        "donors": donors, "ablate": ablate, "dilution_rate": D,
        "vmax_bai": vmax_bai, "fmt_biomass": fmt_biomass, "fmt_time": fmt_time,
        "latched": result.meta.get("spore_latched") or {},
    }


if __name__ == "__main__":
    cfg = json.loads(Path(sys.argv[1]).read_text())
    out = run_cfg(cfg)
    od = HERE / "results"; od.mkdir(parents=True, exist_ok=True)
    (od / f"{cfg['label']}.json").write_text(json.dumps(out))
    print(json.dumps(out), flush=True)
