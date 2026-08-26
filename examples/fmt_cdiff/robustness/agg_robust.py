#!/usr/bin/env python3
"""Aggregate the robustness result JSONs into (1) a specificity / Buffie-bad-case table
and (2) a sensitivity summary that reports whether the SIGN of the bile benefit
(lambda_off > lambda_on, and C.diff_off > C.diff_on) holds across every perturbation."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results"

R = {}
for p in sorted(RES.glob("*.json")):
    d = json.loads(p.read_text())
    R[d["label"]] = d


def show(label):
    d = R.get(label)
    if not d:
        return f"{label:22} MISSING"
    return (f"{label:26} lam={d['lam']:+.4f}  C.diff={d['cdiff_total']:.4f}  "
            f"DCA={d['dca_final']}  scindens_seeded={d['scindens_seeded']}")


print("=" * 78)
print("SPECIFICITY / BUFFIE BAD CASE  (leave-one-donor-out; bile layer ON throughout)")
print("=" * 78)
print(show("ref_bile_on"))
print(show("ref_bile_off"))
DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
          "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]
for d in DONORS:
    print(show(f"loo_{d}"))
on = R.get("ref_bile_on"); off = R.get("ref_bile_off")
loo_k = R.get("loo_C_scindens_ATCC35704")
if on and off and loo_k:
    print(f"\n  interpretation: removing C. scindens drives lambda "
          f"{on['lam']:+.4f} -> {loo_k['lam']:+.4f} (target = bile-off {off['lam']:+.4f}); "
          f"removing any other donor should stay near bile-on {on['lam']:+.4f}.")

print("\n" + "=" * 78)
print("SENSITIVITY  (bile benefit = lambda_off - lambda_on > 0, C.diff_off - C.diff_on > 0)")
print("=" * 78)


def sweep(name, vals, fmt):
    rows = []
    for v in vals:
        onk = R.get(f"{name}_{fmt(v)}_on"); offk = R.get(f"{name}_{fmt(v)}_off")
        if not onk or not offk:
            continue
        dlam = offk["lam"] - onk["lam"]
        dcd = offk["cdiff_total"] - onk["cdiff_total"]
        rows.append((v, onk["lam"], offk["lam"], dlam, dcd,
                     onk["cdiff_total"], offk["cdiff_total"]))
    return rows


def report(title, rows):
    print(f"\n{title}")
    print(f"  {'value':>8} {'lam_on':>8} {'lam_off':>8} {'d_lam':>8} "
          f"{'Cd_on':>7} {'Cd_off':>7} {'d_Cd':>7}  bile_helps")
    ok = True
    for v, lon, loff, dl, dc, con, coff in rows:
        helps = (dl > 0) and (dc > 0)
        ok = ok and helps
        print(f"  {v:>8g} {lon:>+8.4f} {loff:>+8.4f} {dl:>+8.4f} "
              f"{con:>7.4f} {coff:>7.4f} {dc:>+7.4f}  {'YES' if helps else 'NO'}")
    print(f"  -> bile benefit sign STABLE across this sweep: {ok}")
    return ok


g = lambda v: f"{v:g}"
allok = True
allok &= report("colonic washout (dilution_rate, 1/h)",
                sweep("dil", [0.0125, 0.02, 0.025, 0.033, 0.05], g))
allok &= report("FMT dose (bolus biomass, gDW/L)",
                sweep("dose", [0.02, 0.035, 0.05, 0.08], g))
allok &= report("FMT timing (bolus time, h)",
                sweep("time", [6.0, 12.0, 24.0, 36.0], g))
print(f"\n==> bile benefit sign stable across ALL {sum(1 for _ in R)} perturbations: {allok}")

# machine-readable roll-up for the figure / manuscript
roll = {"specificity": {k: R.get(k) for k in
                        ["ref_bile_on", "ref_bile_off"] + [f"loo_{d}" for d in DONORS]},
        "sign_stable_all": allok}
(HERE / "robustness_summary.json").write_text(json.dumps(roll, indent=2))
print(f"\nwrote {HERE/'robustness_summary.json'}")
