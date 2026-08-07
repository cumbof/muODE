#!/usr/bin/env python3
"""Anchor vmax_bai to the measured deoxycholate concentration (not fit-to-clearance).

vmax_bai (community 7a-dehydroxylation rate) is a library INVENTED parameter -- there
is no published whole-community bai turnover rate -- and muode/bile.py demands that any
claim resting on it cite a sweep. C. scindens is a genuine low-abundance keystone
(0.0008/h on the published western_gut diet), so at the default vmax_bai=1.0 the model
makes only 0.093 mM deoxycholate while 3.9 mM cholate sits unconverted, well below the
inhibitory Ki (0.5 mM) and the in-vivo secondary-bile-acid level (~0.3-0.7 mM DCA in a
healthy colon).

This sweeps vmax_bai on the fmt_full arm and records, per value, the FINAL deoxycholate
concentration and the rebound lambda (log-linear slope of the last 30 h; <0 = washout).
The grounding is: pick the vmax_bai that reproduces the MEASURED in-vivo DCA level, then
report whether the bile arm achieves durable clearance THERE -- anchored to the observed
metabolite, not tuned to the outcome. One vmax_bai per process (run in parallel).
"""
import sys, json, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import scenario as gs


def rebound(series, window_h=30.0):
    t = np.asarray(series.index, float); x = np.asarray(series.values, float)
    m = (t >= t[-1] - window_h) & (x > 0)
    return float(np.polyfit(t[m], np.log(x[m]), 1)[0]) if m.sum() >= 2 else float("nan")


def main():
    vb = float(sys.argv[1])
    res = gs.build_scenario(fmt=True, ablate="", vmax_bai=vb)
    M = res.metabolites
    out = {
        "vmax_bai": vb,
        "dca_final_mM": round(float(M["dca_e"].iloc[-1]), 4) if "dca_e" in M else None,
        "ca_final_mM": round(float(M["ca_e"].iloc[-1]), 4) if "ca_e" in M else None,
        "cdiff_final": round(float(res.biomass[gs.PATHOGEN].iloc[-1]), 4),
        "rebound_lambda": round(rebound(res.biomass[gs.PATHOGEN]), 4),
    }
    out["excluded"] = out["rebound_lambda"] < 0
    Path("results/fmt").mkdir(parents=True, exist_ok=True)
    Path(f"results/fmt/sweep_vmaxbai_{vb:g}.json").write_text(json.dumps(out))
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
