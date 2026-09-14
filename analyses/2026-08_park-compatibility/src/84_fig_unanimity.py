#!/usr/bin/env python3
r"""
FIG C -- per-tape unanimity, two panels.  From `83`.

  LEFT   how often is a tape missing in ALL k cells of a set
  RIGHT  how often is it PRESENT in all k -- the discriminator

x = capture-matched random sets, y = sets of close relatives, log-log, one point
per tape.  THE DIAGONAL IS THE NULL, drawn from data: under no lineage structure
the two set types are exchangeable and every tape sits on the line.

WHY TWO PANELS.  Shared capture would cluster BOTH unanimity types -- uniformly
well-captured cells agree on presence, uniformly badly-captured cells agree on
absence.  One-directional loss should lift MISSING and leave PRESENT nearly alone,
since presence is merely the default state.

⭐ COMMON-TAPE CONVENTION (Justin, 2026-09-14).  A tape enters the MISSING panel
only if all-k-missing happened at least once in both set types, which needs a high
missing rate; the PRESENT panel needs the opposite.  The two tests therefore
exclude OPPOSITE ENDS of the per-tape rate range, and the panels were being
computed on different tape populations -- Pre-TX 96 vs 152.  All statistics are now
computed on the tapes eligible for BOTH, so the two medians are comparable.  It
also strengthens the result, because the common set drops near-dead tapes where
both set types reach unanimity automatically and the fold is pinned near 1.

⚠ Tapes where the matched control NEVER reached unanimity have no x to plot on a
log axis.  They are drawn at the left-edge floor 0.5/nsets.  That column is a
CENSORING CONVENTION, not a measurement -- it is stated in the caption.
⚠ Tapes where the RELATIVES never reached unanimity but the control did lie against
the effect and also cannot be plotted; they are counted in the caption.
⚠ Anchors are sampled and sets overlap, so no confidence interval is drawn.

Usage: 84_fig_unanimity.py [arm] [k]        (default Subclone 5)
"""
from pathlib import Path
import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARM = sys.argv[1] if len(sys.argv) > 1 else "Subclone"
KUSE = sys.argv[2] if len(sys.argv) > 2 else "5"
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

d = json.loads((RES / f"unanimity_{ARM}_mincl100.json").read_text())
ns = np.array(d["nset"][KUSE], float)
ok0 = ns > 0
FLOOR = 0.5 / ns.max()
P = {}
for key in ("miss", "pres"):
    pr = np.where(ok0, np.array(d[key][KUSE]["rel"], float) / np.maximum(ns, 1), 0.0)
    pm = np.where(ok0, np.array(d[key][KUSE]["mat"], float) / np.maximum(ns, 1), 0.0)
    P[key] = (pr, pm)
BOTH = ok0.copy()
for key in ("miss", "pres"):
    BOTH &= (P[key][0] > 0) & (P[key][1] > 0)

fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.9), sharex=True, sharey=True)
lim = (FLOOR / 2.2, 1.35)
MED, NOTE = {}, {}
for ax, key, ttl in ((axs[0], "miss", f"missing in ALL {KUSE} cells"),
                     (axs[1], "pres", f"present in ALL {KUSE} cells")):
    pr, pm = P[key]
    ax.plot(lim, lim, color=NEUT2, lw=1.1, zorder=2)
    ax.scatter(pm[BOTH], pr[BOTH], s=26, color=COL[ARM], alpha=0.8, zorder=4,
               edgecolor=SURFACE, linewidth=0.5)
    cens = ok0 & (pr > 0) & (pm == 0)                  # control never unanimous
    # ⭐ OPEN markers, restored 2026-09-14 on Justin's second thought.  This column
    # sits at a CENSORING FLOOR, not at a measured x, and with filled markers a
    # reader takes it for a measurement.  The open style is the only visual signal
    # that these tapes are 'below the resolution of the control', which is also
    # where the effect is strongest -- 49 tapes on Pre-TX at k = 5.
    if cens.sum():
        ax.scatter(np.full(int(cens.sum()), FLOOR), pr[cens], s=30, facecolor="none",
                   edgecolor=COL[ARM], linewidth=1.0, zorder=4)
    MED[key] = float(np.median(pr[BOTH] / pm[BOTH]))
    NOTE[key] = (int((pr[BOTH] > pm[BOTH]).sum()), int(BOTH.sum()),
                 int(cens.sum()), int((ok0 & (pr == 0) & (pm > 0)).sum()))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
    ax.set_xlabel("capture-matched random sets")
    ax.set_title(ttl, fontsize=10.5, pad=8, color=INK)
axs[0].set_ylabel("sets of close relatives")
cm = d["cap_mean"][KUSE]
am, bm, ac, ab = NOTE["miss"]
_, _, pc, pb = NOTE["pres"]
cap = (
    f"{NICE[ARM]} · one point per tape · $k$ = {KUSE} · the grey diagonal is the null: with no "
    f"lineage structure the two set types are exchangeable and every tape sits on it.\n"
    f"Median fold change, relatives over capture-matched random: "
    f"{MED['miss']:.2f}× for MISSING against {MED['pres']:.2f}× for PRESENT, "
    f"a ratio of {MED['miss']/MED['pres']:.1f}×. {am} of {bm} tapes lie above the line on the left.\n"
    f"Statistics use the {bm} tapes measurable in BOTH panels — a tape needs a high missing rate to "
    f"reach all-{KUSE}-missing and a low one to reach all-{KUSE}-present, so the two tests otherwise "
    f"exclude opposite ends of the range.\n"
    f"\u26a0 The open circles at the left edge are a censoring floor, not a measurement: "
    f"{ac} tapes (left) and {pc} (right) where the matched control never reached "
    f"unanimity, plotted at 0.5/n. "
    f"{ab} and {pb} tapes lie below the line off the log axis.\n"
    f"Capture is matched cell by cell on rank: {cm['rel']:.1f} vs {cm['mat']:.1f} A-tapes recovered "
    f"per set, against {cm['rnd']:.1f} for unmatched random sets."
)
fig.text(0.5, -0.02, cap, ha="center", va="top", fontsize=8.0, color=INK2,
         linespacing=1.65, wrap=True)
fig.savefig(FIG / f"fig6c_unanimity_{ARM}_k{KUSE}.png"); plt.close(fig)
print(f"wrote figures/fig6c_unanimity_{ARM}_k{KUSE}.png  "
      f"miss {MED['miss']:.2f} pres {MED['pres']:.2f} ratio {MED['miss']/MED['pres']:.2f} "
      f"n_common {bm} censored {ac}/{pc}")
