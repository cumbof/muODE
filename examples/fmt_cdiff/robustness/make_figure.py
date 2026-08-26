#!/usr/bin/env python3
"""Two-panel robustness figure for the FMT / C. difficile keystone analysis.
Panel A: leave-one-donor-out specificity -- DCA + rebound lambda; only C. scindens removal
         abolishes deoxycholate (the Buffie-2015 keystone bad case).
Panel B: operating envelope -- bile benefit (d_lambda = lambda_off - lambda_on) across
         washout / dose / timing sweeps, shaded where the benefit holds.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
R = {json.loads(p.read_text())["label"]: json.loads(p.read_text()) for p in RES.glob("*.json")}

DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
          "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]
SHORT = {"R_intestinalis_L182": "R. intestinalis", "F_prausnitzii_A2165": "F. prausnitzii",
         "B_thetaiotaomicron_VPI5482": "B. theta", "C_scindens_ATCC35704": "C. scindens"}

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})

# ---- Panel A: leave-one-out specificity -------------------------------------
labels = ["all 4\ndonors"] + [f"− {SHORT[d]}" for d in DONORS]
keys = ["ref_bile_on"] + [f"loo_{d}" for d in DONORS]
dca = [R[k]["dca_final"] for k in keys]
lam = [R[k]["lam"] for k in keys]
x = np.arange(len(labels))
# keystone (C. scindens removal) highlighted
colors = ["#4c78a8"] + ["#9ecae1", "#9ecae1", "#9ecae1", "#d62728"]

barsA = axA.bar(x, dca, color=colors, width=0.62, edgecolor="black", linewidth=0.5)
axA.set_ylabel("secondary bile acid\n(deoxycholate, mM)")
axA.set_xticks(x)
axA.set_xticklabels(labels, fontsize=8.5)
axA.set_title("A  Keystone specificity (leave-one-donor-out)", loc="left", fontweight="bold", fontsize=11)
axA.axhline(0, color="black", lw=0.6)
for xi, v in zip(x, dca):
    axA.text(xi, v + 0.004, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
# overlay lambda (rebound rate) on twin axis
axA2 = axA.twinx()
axA2.spines["top"].set_visible(False)
axA2.plot(x, lam, "o-", color="#333333", lw=1.3, ms=5, label="C. diff rebound λ")
axA2.set_ylabel("C. difficile rebound λ (1/h)")
axA2.set_ylim(0, max(lam) * 1.35)
axA.annotate("only C. scindens removal\nabolishes DCA → worst rebound",
             xy=(4, dca[4]), xytext=(2.1, max(dca) * 0.7),
             fontsize=8, color="#d62728",
             arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.1))
axA2.legend(loc="upper left", fontsize=8, frameon=False)

# ---- Panel B: operating envelope --------------------------------------------
def sweep(name, vals):
    out = []
    for v in vals:
        on = R.get(f"{name}_{v:g}_on"); off = R.get(f"{name}_{v:g}_off")
        if on and off:
            out.append((v, off["lam"] - on["lam"]))
    return out

sweeps = [
    ("washout (1/h)", "dil", [0.0125, 0.02, 0.025, 0.033, 0.05], "#4c78a8"),
    ("FMT dose (gDW/L)", "dose", [0.02, 0.035, 0.05, 0.08], "#59a14f"),
    ("FMT timing (h)", "time", [6.0, 12.0, 24.0, 36.0], "#b07aa1"),
]
# normalise each x to its center value's index for a shared axis
for title, name, vals, col in sweeps:
    rows = sweep(name, vals)
    xs = [f"{v:g}" for v, _ in rows]
    dl = [d for _, d in rows]
    axB.plot(range(len(xs)), dl, "o-", color=col, lw=1.4, ms=5, label=title)

axB.axhline(0, color="black", lw=0.8)
axB.axhspan(0, axB.get_ylim()[1] if axB.get_ylim()[1] > 0 else 0.01, color="#e8f5e9", zorder=0)
axB.set_ylabel("bile benefit  Δλ = λ$_{off}$ − λ$_{on}$")
axB.set_xlabel("parameter value (low → high within each sweep)")
axB.set_title("B  Operating envelope (washout / dose / timing)", loc="left", fontweight="bold", fontsize=11)
axB.text(0.02, 0.94, "benefit region (Δλ > 0)", transform=axB.transAxes,
         fontsize=8, color="#2e7d32", va="top")
axB.legend(fontsize=8, frameon=False, loc="lower left")
axB.set_xticks(range(4))

fig.tight_layout()
out = HERE / "fig_fmt_robustness.png"
fig.savefig(out, bbox_inches="tight")
fig.savefig(HERE / "fig_fmt_robustness.pdf", bbox_inches="tight")
print("wrote", out)
