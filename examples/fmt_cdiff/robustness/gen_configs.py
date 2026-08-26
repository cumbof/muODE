#!/usr/bin/env python3
"""Generate the robustness config set: specificity (leave-one-donor-out / Buffie bad case)
and sensitivity (perturb the invented/derived knobs; keep BOTH bile-on and bile-off at each
point so the SIGN of the bile benefit can be checked, not just the absolute lambda)."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CFG = HERE / "configs"; CFG.mkdir(exist_ok=True)

DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
          "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]
KEYSTONE = "C_scindens_ATCC35704"

configs = []

# --- reference points (reproduce the published decomposition operating point) ----------
configs.append({"label": "ref_bile_on", "ablate": ""})
configs.append({"label": "ref_bile_off", "ablate": "bile"})

# --- specificity / Buffie bad case: leave each donor out of the FMT bolus --------------
# bile layer stays ON; a donor absent from the bolus never establishes, so leaving out
# C. scindens removes the 7a-dehydroxylase in the biologically real way.
for d in DONORS:
    configs.append({"label": f"loo_{d}", "donors": [x for x in DONORS if x != d]})

# --- sensitivity: colonic washout (derived 0.021-0.042/h; 0.025 is 40 h transit) -------
for D in [0.0125, 0.02, 0.025, 0.033, 0.05]:
    configs.append({"label": f"dil_{D:g}_on",  "dilution_rate": D, "ablate": ""})
    configs.append({"label": f"dil_{D:g}_off", "dilution_rate": D, "ablate": "bile"})

# --- sensitivity: FMT dose (bolus biomass) --------------------------------------------
for fb in [0.02, 0.035, 0.05, 0.08]:
    configs.append({"label": f"dose_{fb:g}_on",  "fmt_biomass": fb, "ablate": ""})
    configs.append({"label": f"dose_{fb:g}_off", "fmt_biomass": fb, "ablate": "bile"})

# --- sensitivity: FMT timing (when the bolus goes in) ---------------------------------
for ft in [6.0, 12.0, 24.0, 36.0]:
    configs.append({"label": f"time_{ft:g}_on",  "fmt_time": ft, "ablate": ""})
    configs.append({"label": f"time_{ft:g}_off", "fmt_time": ft, "ablate": "bile"})

for c in configs:
    (CFG / f"{c['label']}.json").write_text(json.dumps(c))
print(f"wrote {len(configs)} configs to {CFG}")
for c in configs:
    print(" ", c["label"])
