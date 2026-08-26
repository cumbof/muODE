#!/usr/bin/env python3
"""Two SEPARATE robustness figures for the FMT / C. difficile keystone analysis.

  fig_fmt_specificity.{png,pdf}  -- MAIN TEXT: leave-one-donor-out. Only C. scindens
      removal abolishes deoxycholate and gives the worst C. diff rebound (Buffie-2015
      keystone bad case; the 'not by chance' specificity control).
  fig_fmt_envelope.{png,pdf}     -- SUPPLEMENTARY: bile benefit d_lambda across
      washout / dose / timing sweeps; robust in the physiological range, fades at extremes.
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

# ============================ FIGURE 1 (MAIN): specificity ====================
labels = ["all 4\ndonors"] + [f"− {SHORT[d]}" for d in DONORS]
keys = ["ref_bile_on"] + [f"loo_{d}" for d in DONORS]
dca = [R[k]["dca_final"] for k in keys]
lam = [R[k]["lam"] for k in keys]
x = np.arange(len(labels))
colors = ["#4c78a8"] + ["#9ecae1", "#9ecae1", "#9ecae1", "#d62728"]

fig1, ax = plt.subplots(figsize=(6.0, 4.4))
ax.bar(x, dca, color=colors, width=0.62, edgecolor="black", linewidth=0.5)
ax.set_ylabel("secondary bile acid\n(deoxycholate, mM)")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
ax.axhline(0, color="black", lw=0.6)
for xi, v in zip(x, dca):
    ax.text(xi, v + 0.004, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
ax2 = ax.twinx(); ax2.spines["top"].set_visible(False)
ax2.plot(x, lam, "o-", color="#333333", lw=1.3, ms=5, label="C. difficile rebound λ")
ax2.set_ylabel("C. difficile rebound λ (1/h)")
ax2.set_ylim(0, max(lam) * 1.35)
ax.annotate("only C. scindens removal\nabolishes DCA → worst rebound",
            xy=(4, dca[4]), xytext=(1.85, max(dca) * 0.72),
            fontsize=8, color="#d62728",
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.1))
ax2.legend(loc="upper left", fontsize=8, frameon=False)
fig1.tight_layout()
fig1.savefig(HERE / "fig_fmt_specificity.png", bbox_inches="tight")
fig1.savefig(HERE / "fig_fmt_specificity.pdf", bbox_inches="tight")
print("wrote fig_fmt_specificity.png/pdf")

# ============================ FIGURE 2 (SUPP): envelope =======================
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
fig2, axB = plt.subplots(figsize=(6.4, 4.4))
maxn = max(len(v) for _, _, v, _ in sweeps)
for title, name, vals, col in sweeps:
    rows = sweep(name, vals)
    xs = [f"{v:g}" for v, _ in rows]
    dl = [d for _, d in rows]
    axB.plot(range(len(xs)), dl, "o-", color=col, lw=1.4, ms=5, label=title)
    for i, (lbl, d) in enumerate(zip(xs, dl)):
        axB.annotate(lbl, (i, d), textcoords="offset points", xytext=(0, -11 if d < 0 else 7),
                     ha="center", fontsize=6.5, color=col)
axB.axhline(0, color="black", lw=0.8)
top = axB.get_ylim()[1]
axB.axhspan(0, max(top, 0.001), color="#e8f5e9", zorder=0)
axB.set_ylabel("bile benefit  Δλ = λ$_{off}$ − λ$_{on}$  (1/h)")
axB.set_xlabel("sweep step (low → high)")
axB.text(0.02, 0.95, "benefit region (Δλ > 0)", transform=axB.transAxes,
         fontsize=8, color="#2e7d32", va="top")
axB.legend(fontsize=8, frameon=False, loc="lower left")
axB.set_xticks(range(maxn))
fig2.tight_layout()
fig2.savefig(HERE / "fig_fmt_envelope.png", bbox_inches="tight")
fig2.savefig(HERE / "fig_fmt_envelope.pdf", bbox_inches="tight")
print("wrote fig_fmt_envelope.png/pdf")
