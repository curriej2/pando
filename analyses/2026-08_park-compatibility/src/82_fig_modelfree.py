#!/usr/bin/env python3
r"""
FIG C, model-free.  Two stacked panels sharing the x axis, from `81`.

x = the NUMBER of a cell's 20 neighbours that are missing a given tape, 0 to 20.
No model, no residuals, no stratification, no training mask -- so the denominator
is exactly 20 and x is a genuine count.  The A/B tape split stays (relatedness from
half A, dropout from half B), because without it a pair that both lack tapes 1-50
would look related BECAUSE their dropout matches.

⚠⚠ WHY TWO PANELS.  The obvious single panel -- rate against j -- DOES NOT SHOW THE
EFFECT, and shows it BACKWARDS.  At matched j the arbitrary clone-mates are the
slightly better predictor in all five arms (median -0.6 to -13.3 percentage points),
and that is correct, not a bug: if 15 of 20 RANDOM cells lack a tape the clone-wide
rate must be high, whereas 15 of 20 RELATIVES can be a local cluster in a clone
where the tape is otherwise fine.  So at matched j the random sample implies a worse
tape.
⇒ the effect lives entirely in the DISTRIBUTION of j.  Relatives cluster, so they
reach unanimity far more often, and unanimity is where the rate is decisive.
The two panels make that decomposition explicit:

  top     how often each count occurs        -- these differ, and that IS the effect
  bottom  the dropout rate at each count     -- these coincide

reading in one sentence: a given level of neighbour agreement means the same thing
whoever the neighbours are, but your relatives reach agreement far more often.

Usage: 82_fig_modelfree.py [arm]        (default Subclone)
"""
from pathlib import Path
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARM = sys.argv[1] if len(sys.argv) > 1 else "Subclone"
NICE = {"Initial": "Pre-TX", "Mouse1": "Mouse 1", "Mouse2": "Mouse 2",
        "Mouse3": "Mouse 3", "Subclone": "Subclone"}
COL = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
       "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
NEUT, NEUT2 = "#cfcec6", "#8f8e86"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

d = json.loads((RES / f"modelfree_{ARM}_mincl100_k20.json").read_text())
K = d["k"]
c = {t: np.array(d["counts"][t]).sum(0) for t in ("rel", "rnd")}
h = {t: np.array(d["hits"][t]).sum(0) for t in ("rel", "rnd")}
share = {t: 100 * c[t] / c[t].sum() for t in c}
rate = {t: 100 * h[t] / np.maximum(c[t], 1) for t in c}
una = {t: share[t][0] + share[t][K] for t in c}
xs = np.arange(K + 1)
SER = (("rel", f"the {K} nearest relatives", COL[ARM]),
       ("rnd", f"{K} arbitrary cells from the same clone", NEUT2))

fig, (a0, a1) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.16))
for t, lab, col in SER:
    a0.plot(xs, share[t], color=col, lw=2.0, marker="o", ms=4.0, label=lab, zorder=3)
a0.set_ylabel("share of all cell-tape pairs (%)")
a0.set_title("how often that many neighbours are missing the tape",
             fontsize=10.5, pad=8, loc="left", color=INK)
a0.legend(frameon=False, fontsize=9.0, loc="upper center", handlelength=1.6)
for t_, col_, tx, ty in (("rel", COL[ARM], K - 2.4, share["rel"][K] + 6.0),
                         ("rnd", NEUT2, K - 4.6, 7.2)):
    a0.annotate(f"{share[t_][K]:.1f}%", xy=(K, share[t_][K]), xytext=(tx, ty),
                fontsize=8.4, color=col_, ha="right", va="center",
                arrowprops=dict(arrowstyle="-", color=col_, lw=0.9))
a0.text(0.25, 0.58, f"all {K} agree it is MISSING:\n"
        f"{share['rel'][K]/share['rnd'][K]:.1f}\u00d7 as often for relatives\n\n"
        f"unanimous either way: {una['rel']:.0f}% vs {una['rnd']:.0f}%",
        transform=a0.transAxes, ha="left", va="top", fontsize=8.4, color=INK2,
        linespacing=1.5)

for t, lab, col in SER:
    a1.plot(xs, rate[t], color=col, lw=2.0, marker="o", ms=4.0, zorder=3)
a1.set_ylabel("observed dropout rate (%)")
a1.set_xlabel(f"number of those {K} neighbours missing the tape")
a1.set_title("what that many neighbours tells you  —  the two curves coincide",
             fontsize=10.5, pad=8, loc="left", color=INK)
a1.set_xticks(np.arange(0, K + 1, 2))
a1.set_ylim(-6, 106)
a1.text(0.02, 0.95, "a given level of agreement means the same thing\n"
        "whoever the neighbours are", transform=a1.transAxes, ha="left", va="top",
        fontsize=8.4, color=INK2)
fig.suptitle(f"{NICE[ARM]} · {d['n_cells']:,} cells in {d['n_clones']} clones · "
             f"relatedness from one half of the tapes, dropout from the other",
             fontsize=8.6, color=INK2, y=0.955)
fig.savefig(FIG / f"fig6c_modelfree_{ARM}.png"); plt.close(fig)
print(f"wrote figures/fig6c_modelfree_{ARM}.png")
print(f"  unanimous: relatives {una['rel']:.1f}%  arbitrary {una['rnd']:.1f}%")
print(f"  j=20 share: relatives {share['rel'][K]:.1f}%  arbitrary {share['rnd'][K]:.1f}%"
      f"  ({share['rel'][K]/share['rnd'][K]:.1f}x)")
