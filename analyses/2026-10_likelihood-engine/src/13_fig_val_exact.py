"""Validation figure 2: the exact checks of the engine's likelihood on one error axis.

Replots results/validate_engine.json (job 17612197, 2026-10-06); computes nothing new. Rows (Justin, 2026-10-13):
V1 literal enumeration, V1 full-state-space pruning, V2 normalisation, V3 JAX vs reference twin. (V4 kernel vs
expm and V8a tree prior dropped from the figure; still PASS, README Session 22/25.) One mark = the WORST error
over all of that check's comparisons (V1 saved only its worst). Error is relative for log-likelihoods
(V1, V3); for V2 it is |sum - 1|. V1-V2 ran with technical dropout, heritable silencing, depth-0 censoring, a
time-varying editing rate and two composition epochs ON at random parameter values; V3 per SETTINGS.
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

# palette: reference categorical slot 1 (validated, light mode), as in src/07
BLUE = "#2a78d6"
INK, INK2, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#ecebe7"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
EPS = np.finfo(float).eps

E = json.load(open(ROOT / "results" / "validate_engine.json"))["detail"]

# (label, what it was compared on, worst error, tolerance)
rows = [
    ("Brute force: add up every possible\nhistory of the ancestors' tapes",
     f"{E['V1']['vs_enumeration']['configs']} tapes, trees of 3–4 cells (up to 41³ = 68,921 histories)",
     E["V1"]["vs_enumeration"]["worst_rel"], 1e-10),
    ("Same calculation without the prefix\nshortcut: every tape string a state",
     f"{E['V1']['vs_full_space']['configs']} tapes, trees of 3–6 cells",
     E["V1"]["vs_full_space"]["worst_rel"], 1e-10),
    ("Probabilities of every possible\nreadout add up to exactly 1",
     f"{E['V2']['configs']} fixed trees + parameters, every readout listed",
     E["V2"]["worst_abs"], 1e-10),
    ("Fast engine (JAX) vs\nslow reference code",
     f"{E['V3']['tape_comparisons']:,} tapes, simulated clones of 25–210 cells",
     E["V3"]["worst_rel"], 1e-10),
]

fig, ax = plt.subplots(figsize=(10.5, 4.3))
ax.axvspan(EPS, 1e-14, color=BAND, zorder=0)
ax.text(np.sqrt(EPS * 1e-14), len(rows) - 0.15, "rounding level", ha="center", va="top",
        fontsize=8.5, color=INK2)
for i, (lab, what, worst, tol) in enumerate(rows):
    y = len(rows) - 1 - i
    ax.plot([tol, tol], [y + 0.02, y + 0.36], color=INK, lw=1.4, ls=(0, (2, 1.5)), zorder=2)
    ax.scatter([worst], [y], s=70, color=BLUE, edgecolors=SURF, lw=1.2, zorder=3)
    ax.text(worst * 3, y + 0.05, f"worst {worst:.0e}", va="bottom", fontsize=9, color=INK)
    ax.text(worst * 3, y - 0.05, what, va="top", fontsize=8, color=INK2)
    ax.text(tol * 1.6, y + 0.2, f"pass line {tol:.0e}", va="center", fontsize=8, color=INK2)
    if i < len(rows) - 1:
        ax.axhline(y - 0.5, color=GRID, lw=0.8, zorder=0)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows][::-1])
ax.set_ylim(-0.6, len(rows) - 0.1)
ax.set_xscale("log"); ax.set_xlim(1e-17, 1e-5)
ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_xlabel("worst disagreement in each check (log scale)")
ax.set_title("The likelihood matches brute-force calculation to rounding error", loc="left")
ax.text(0, -0.22, "Top three rows: technical dropout, heritable silencing, censoring of blank tapes, a time-varying "
        "editing rate and a change of symbol\ncomposition all ON, at random parameter values. JAX row: dropout and "
        "censoring in all 90 simulated clones, silencing in 67, time-varying\nrate and composition in 22. "
        "Error is relative for log-likelihoods; for the sum-to-1 row it is the difference from 1.",
        transform=ax.transAxes, fontsize=8.5, color=INK2, va="top")

fig.tight_layout()
f = OUT / "fig_val2_exact.png"
fig.savefig(f, dpi=200, bbox_inches="tight")
print(f"wrote {f}")
for lab, what, worst, tol in rows:
    print(f"  {lab.splitlines()[0]:42s} worst {worst:.2e}  tol {tol:.0e}  margin {tol / worst:.1e}x")
