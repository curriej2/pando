#!/usr/bin/env python3
r"""
FIGURE 3, PANEL a -- tape dropout is not a coin flip on entries.

NULL: every (cell, tape) instance is recovered independently with one common
probability p, so a cell's recovered-tape count is R_c ~ Bin(k, p), k = 166,
Var(R_c) = k p (1-p).  Plotted as the orange curve; it is a pmf over integers, so
it is directly comparable to the observed fraction-of-cells-per-integer bars with
no rescaling.

READ-OFF: the two horizontal bars beneath the distributions are +/- 1 sd of each,
drawn to the same x-scale.  Their ratio squared is the variance inflation factor

    VIF = Var_obs(R_c) / [k p (1-p)] = 1 + (k-1) rho

with rho the correlation between two dropout indicators in the same cell.
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import binom
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARM, THR = "Mouse1", 20

z = np.load(RES / f"dropout_matrix_{ARM}.npz", allow_pickle=False)
Y = z["recovered"]
n, k = Y.shape
R = Y.sum(1)
p = float(Y.mean())
obs = np.bincount(R, minlength=k + 1)[: k + 1] / n
x = np.arange(k + 1)
null = binom.pmf(x, k, p)

sd_obs = float(R.std(ddof=1))
sd_null = float(np.sqrt(k * p * (1 - p)))
vif = sd_obs**2 / sd_null**2
rho = (vif - 1) / (k - 1)
sigma = np.sqrt(rho * p * (1 - p))
print(f"{ARM}: n={n:,} p={p:.4f} mean R_c={R.mean():.1f} sd={sd_obs:.1f} "
      f"null sd={sd_null:.2f} VIF={vif:.1f} rho={rho:.4f} sigma={sigma:.3f} "
      f"below100={100*(R<100).mean():.1f}%")

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 14})

fig, ax = plt.subplots(figsize=(9.0, 5.6))
top = max(obs.max(), null.max()) * 1.02
band = -0.085 * top                      # y of the sd read-off bars
ax.set_ylim(band * 1.9, top)

ax.bar(x, obs, width=1.0, color=BLUE, linewidth=0, zorder=2, label="observed")
ax.fill_between(x, 0, null, color=ORANGE, alpha=0.15, zorder=1)
ax.plot(x, null, color=ORANGE, lw=2.6, zorder=3,
        label="if dropout were a coin flip")

# the +/- 1 sd read-off, both to the same x-scale
for y0, c, s, lab in ((band, ORANGE, sd_null, f"null   sd {sd_null:.1f}"),
                      (band * 1.62, BLUE, sd_obs, f"observed   sd {sd_obs:.0f}")):
    m = R.mean() if c == BLUE else k * p
    ax.plot([m - s, m + s], [y0, y0], color=c, lw=3.4, solid_capstyle="butt", zorder=4)
    ax.plot([m - s, m - s, np.nan, m + s, m + s], [y0 - 0.012 * top, y0 + 0.012 * top,
            np.nan, y0 - 0.012 * top, y0 + 0.012 * top], color=c, lw=1.6, zorder=4)
    ax.text(m + s + 3, y0, lab, color=c, fontsize=11.5, va="center", ha="left")

ax.axhline(0, color=AXIS, lw=0.9, zorder=1)
ax.axvline(THR, color=INK2, lw=1.3, ls=(0, (4, 3)), zorder=2)
ax.text(THR + 2.5, top * 0.5, f"cells admitted\nat $\\geq${THR} tapes",
        color=INK2, fontsize=10.5, va="center")

ax.text(0.985, 0.97,
        "$R_c \\sim \\mathrm{Bin}(k,p)$,   $\\mathrm{Var}=k\\,p(1-p)$\n"
        f"$k=166$,  $p={p:.3f}$\n\n"
        "$\\mathrm{VIF}=\\dfrac{\\mathrm{Var_{obs}}(R_c)}{k\\,p(1-p)}"
        f"={vif:.0f}$",
        transform=ax.transAxes, ha="right", va="top", fontsize=12.5, color=INK,
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=AXIS, lw=0.8))
ax.legend(frameon=False, fontsize=11.5, loc="upper left", bbox_to_anchor=(0.02, 0.99))

ax.set_xlim(-2, k + 2)
ax.set_xlabel("tapes recovered per cell   $R_c$   (of 166)")
ax.set_ylabel("fraction of cells")
ax.set_yticks([t for t in ax.get_yticks() if t >= 0])
ax.set_title("a   Tape dropout is not a coin flip on entries", loc="left", pad=12)

fig.tight_layout()
fig.savefig(FIG / "fig3a_null.png", dpi=200, facecolor=SURF)
print("wrote figures/fig3a_null.png")
