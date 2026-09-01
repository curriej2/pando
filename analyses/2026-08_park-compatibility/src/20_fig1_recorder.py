#!/usr/bin/env python3
r"""FIGURE 1 -- what the recorder can actually record.  Two panels.

(a) q accumulated symbol by symbol, against the uniform-xi counterfactual. The
    reader sees where q comes from, that a few symbols supply most of it, and that
    uneven usage is what inflates it.
(b) q across all five experimental arms.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Initial", "Subclone", "Mouse1", "Mouse2", "Mouse3"]
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
REF = "Mouse1"

xi_all = json.loads((RES / "xi_vectors.json").read_text())
xi = np.sort(np.array(list(xi_all[REF]["xi"].values())))[::-1]
xi = xi / xi.sum()
M = xi.size
q_obs = float((xi ** 2).sum())
q_uni = 1.0 / M

cum_obs = np.cumsum(xi ** 2)
cum_uni = np.arange(1, M + 1) / M ** 2          # k * (1/M)^2
k_half = int(np.searchsorted(cum_obs, q_obs / 2)) + 1

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 13.5})

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.2), gridspec_kw={"width_ratios": [1.35, 1]})
r = np.arange(1, M + 1)

axA.plot(r, cum_uni, color=AXIS, lw=2.4)
axA.plot(r, cum_obs, color=BLUE, lw=2.6)
axA.axvline(k_half, color=ORANGE, lw=1.5, ls=(0, (5, 3)))
axA.plot([k_half], [q_obs / 2], "o", ms=9, color=ORANGE,
         markeredgecolor=SURF, markeredgewidth=1.8, zorder=5)
axA.text(M * 0.97, q_obs, f"observed  $q$ = {q_obs:.4f}", color=BLUE,
         ha="right", va="bottom", fontsize=12, fontweight="bold")
axA.text(M * 0.97, q_uni, f"uniform $\\xi$   $q$ = {q_uni:.4f}", color=INK2,
         ha="right", va="bottom", fontsize=12)
axA.text(k_half + 3, q_obs / 2, f"{k_half} symbols give\nhalf of $q$", color=ORANGE,
         fontsize=11.5, va="center")
axA.set_xlim(0, M + 1); axA.set_ylim(0, q_obs * 1.16)
axA.set_xlabel("insertion symbols, ranked by frequency")
axA.set_ylabel("$q=\\sum_i \\xi_i^{\\,2}$   accumulated")
axA.set_title("a   Uneven use of the alphabet nearly doubles $q$", loc="left", pad=10)

qs = []
for a in ARMS:
    v = np.array(list(xi_all[a]["xi"].values())); qs.append(float(((v / v.sum()) ** 2).sum()))
axB.bar([LABEL[a] for a in ARMS], qs, color=BLUE, width=0.62, edgecolor="none")
for i, v in enumerate(qs):
    axB.text(i, v + 0.0004, f"{v:.4f}", ha="center", fontsize=11.5, color=INK)
axB.set_ylim(0, max(qs) * 1.22)
axB.set_ylabel("$q=\\sum_i \\xi_i^{\\,2}$")
axB.set_title("b   $q$ is stable across all five arms", loc="left", pad=10)
axB.tick_params(axis="x", labelsize=11)

fig.tight_layout()
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "fig1_recorder.png", dpi=200, facecolor=SURF)
print(f"q_obs={q_obs:.6f}  q_uniform={q_uni:.6f}  ratio={q_obs/q_uni:.2f}x  "
      f"M={M}  half of q from top {k_half}")
