"""Validation figure 4: the engine's gradients (what HMC runs on) against finite differences (check V5).

a  every V5 comparison (results/validate_engine.json, 2026-10-06): relative error
   |fd - autodiff| / max(|fd|, 1), one row per quantity, log axis, against the 1e-5 pass line.
b  the one entry over the line (tree 16, d/dGam_T at Gam_T = 0): finite-difference error against step size h,
   from results/v5_ladder.json (14_probe_grad_ladder.py, 2026-10-13, reproduces the recorded values exactly).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# palette: reference categorical slots 1-2 (validated all-pairs, light mode), as in src/07
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
TOL = 1e-5

V5 = json.load(open(ROOT / "results" / "validate_engine.json"))["detail"]["V5"]["rows"]
LAD = json.load(open(ROOT / "results" / "v5_ladder.json"))

groups = [("x", "node times\n(240 coordinates)"), ("Lam_T", "editing rate Λ$_T$"),
          ("Gam_T", "silencing exposure Γ$_T$"), ("b", "birth rate b"), ("d", "death rate δ"),
          ("beta", "tape dropout offset β"), ("alpha", "cell dropout offset α")]

fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.35, 1]})

# ---- a: every comparison ----
rng = np.random.default_rng(0)
for gi, (key, lab) in enumerate(groups):
    y = len(groups) - 1 - gi
    rows = [r for r in V5 if r["param"] == key]
    e = np.array([r["err"] for r in rows])
    out = e > TOL
    jit = rng.uniform(-0.22, 0.22, e.size)
    axa.scatter(e[~out], y + jit[~out], s=22, color=BLUE, edgecolors=SURF, lw=0.6, zorder=3)
    if out.any():
        axa.scatter(e[out], y + jit[out] * 0, s=70, color=ORANGE, marker="D", edgecolors=SURF, lw=1, zorder=4)
        axa.annotate("at Γ$_T$ = 0 exactly\n(a clone with no silencing):\nsee panel b",
                     (e[out][0], y), xytext=(1.6e-5, y - 1.0), fontsize=8.5, color=INK,
                     ha="left", va="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    if gi < len(groups) - 1:
        axa.axhline(y - 0.5, color=GRID, lw=0.8, zorder=0)
axa.axvline(TOL, color=INK, lw=1.2, ls="--", zorder=1)
axa.text(TOL * 1.3, len(groups) - 0.55, "pass line 10⁻⁵", fontsize=9, va="top")
axa.set_yticks(range(len(groups))); axa.set_yticklabels([g[1] for g in groups][::-1])
axa.set_ylim(-0.6, len(groups) - 0.4)
axa.set_xscale("log"); axa.set_xlim(1e-12, 1e-2)
axa.tick_params(axis="y", length=0); axa.spines["left"].set_visible(False)
axa.grid(axis="x", color=GRID, lw=0.8); axa.set_axisbelow(True)
axa.set_xlabel("relative error, autodiff vs finite difference (log scale)")
axa.set_title(f"a  {sum(r['err'] <= TOL for r in V5)} of {len(V5)} gradient entries agree", loc="left")
axa.text(0, -0.14, "20 simulated clones (25–64 cells, 30 tapes); constant, silencing and time-varying settings.\n"
         "Finite differences at fixed step h = 10⁻⁵ are themselves only accurate to ~10⁻⁸, which sets the floor here.",
         transform=axa.transAxes, fontsize=8.5, color=INK2, va="top")

# ---- b: the outlier's step ladder ----
lad = [r for r in LAD["ladder"] if r["rel_err_vs_autodiff"] is not None]
h = np.array([r["h"] for r in lad]); err = np.array([r["rel_err_vs_autodiff"] for r in lad])
axb.plot(h, err, color=ORANGE, lw=2, zorder=2)
axb.scatter(h, err, s=60, color=ORANGE, marker="D", edgecolors=SURF, lw=1, zorder=3)
k = int(np.argmin(np.abs(h - 1e-5)))
axb.annotate("step used in the check\n(h = 10⁻⁵)", (h[k], err[k]), xytext=(h[k] / 40, err[k] * 6),
             fontsize=8.5, ha="left", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axb.axhline(TOL, color=INK, lw=1.2, ls="--", zorder=1)
axb.text(h.min(), TOL * 1.6, "pass line 10⁻⁵", fontsize=9, ha="right")
axb.annotate("rounding error\ntakes over", (h[-1], err[-1]), xytext=(h[-1] * 3, err[-1] * 25), fontsize=8.5,
             ha="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axb.set_xscale("log"); axb.set_yscale("log"); axb.invert_xaxis()
axb.grid(color=GRID, lw=0.8); axb.set_axisbelow(True)
axb.set_xlabel("finite-difference step h (smaller →)")
axb.set_ylabel("relative error vs autodiff")
axb.set_title("b  The outlier: the reference converges onto autodiff", loc="left")
axb.text(0, -0.14, f"d(log-posterior)/dΓ$_T$ at Γ$_T$ = 0; autodiff = {LAD['autodiff']:,.2f} nats per unit Γ$_T$.\n"
         "Error falls ~100× per 10× smaller step (the coarse reference, not the gradient).",
         transform=axb.transAxes, fontsize=8.5, color=INK2, va="top")

fig.tight_layout()
f = OUT / "fig_val4_gradients.png"
fig.savefig(f, dpi=200, bbox_inches="tight")
print(f"wrote {f}")
for r in LAD["ladder"]:
    print(f"  h={r['h']:.0e}  rel err {r['rel_err_vs_autodiff']}")
