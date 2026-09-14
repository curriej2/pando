#!/usr/bin/env python3
r"""
FIG C v3 -- per-tape unanimity, two panels.  From `83`.

  LEFT   how often is a tape missing in ALL k cells of a set
  RIGHT  how often is it PRESENT in all k -- the discriminator

x = capture-matched random sets, y = sets of close relatives, log-log, one point
per tape.  THE DIAGONAL IS THE NULL, and it is drawn from data: under no lineage
structure the two set types are exchangeable and every tape sits on the line.

WHY TWO PANELS.  Shared capture quality would cluster BOTH unanimity types --
uniformly well-captured cells agree on presence, uniformly badly-captured cells
agree on absence.  One-directional loss should lift MISSING and leave PRESENT
alone, because presence is merely the default state.  So the left panel rising off
the diagonal while the right panel stays on it is evidence the effect is not
capture, obtained without any model.

⚠ Tapes with zero matched-random unanimity cannot be placed on a log axis.  They
are drawn as OPEN markers on the left edge at 0.5/nsets and counted in the caption:
they are tapes where the effect is strongest, so dropping them silently would
understate it.

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
FOLD = {}
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.9), sharex=True, sharey=True)
for ax, key, ttl in ((axs[0], "miss", f"missing in ALL {KUSE} cells"),
                     (axs[1], "pres", f"present in ALL {KUSE} cells")):
    pr = np.where(ok0, np.array(d[key][KUSE]["rel"], float) / np.maximum(ns, 1), np.nan)
    pm = np.where(ok0, np.array(d[key][KUSE]["mat"], float) / np.maximum(ns, 1), np.nan)
    live = ok0 & (pr > 0)
    on = live & (pm > 0)                       # both measurable
    off = live & (pm == 0)                     # matched-random never unanimous
    lim = (FLOOR / 2.2, 1.35)
    ax.plot(lim, lim, color=NEUT2, lw=1.1, zorder=2)
    ax.scatter(pm[on], pr[on], s=26, color=COL[ARM], alpha=0.8, zorder=4,
               edgecolor=SURFACE, linewidth=0.5)
    if off.sum():
        ax.scatter(np.full(int(off.sum()), FLOOR), pr[off], s=30, facecolor="none",
                   edgecolor=COL[ARM], linewidth=1.0, zorder=4)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("capture-matched random sets")
    ax.set_title(ttl, fontsize=10.5, pad=8, color=INK)
    fold = pr[on] / pm[on]
    med = float(np.median(fold))
    above = int((pr[on] > pm[on]).sum())
    FOLD[key] = med
    # ⚠ "present stays ON the diagonal" was too strong: every present tape is
    # above the line too, just by far less.  Draw the median fold as a line
    # PARALLEL to the diagonal so the magnitude is visible rather than only
    # annotated, and let the two medians carry the discriminator between them.
    xx = np.array(lim)
    ax.plot(xx, med * xx, color=COL[ARM] if key == "miss" else INK2, lw=1.1,
            ls=(0, (5, 3)), zorder=3)
    ax.text(0.03, 0.97, f"median {med:.2f}x the matched random rate\n"
            f"{above} of {int(on.sum())} tapes above the line"
            + (f"\n{int(off.sum())} more off-scale (open)" if off.sum() else ""),
            transform=ax.transAxes, ha="left", va="top", fontsize=8.6,
            color=COL[ARM] if key == "miss" else INK2, linespacing=1.5)
axs[0].set_ylabel("sets of close relatives")
axs[1].text(0.97, 0.03, f"missing clusters {FOLD['miss']/FOLD['pres']:.1f}x more"
            f"\nthan presence does", transform=axs[1].transAxes, ha="right",
            va="bottom", fontsize=8.8, color=INK, linespacing=1.5)
cm = d["cap_mean"][KUSE]
fig.suptitle(f"{NICE[ARM]} · one point per tape · the diagonal is the null\n"
             f"capture matched cell by cell: {cm['rel']:.1f} vs {cm['mat']:.1f} A-tapes "
             f"recovered per set (unmatched random would be {cm['rnd']:.1f})",
             fontsize=8.8, color=INK2, y=1.03)
fig.savefig(FIG / f"fig6c_unanimity_{ARM}_k{KUSE}.png"); plt.close(fig)
print(f"wrote figures/fig6c_unanimity_{ARM}_k{KUSE}.png")
