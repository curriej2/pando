#!/usr/bin/env python3
r"""
COMPANION TO FIG C -- what the matching equalised, and what it did not.

Three distributions over SETS (not cells), from `83`:
  left    set MEAN capture                    -- equal by construction
  middle  WITHIN-set sd of capture            -- equal too, and THIS is what drives
                                                 unanimity, since a uniformly badly
                                                 captured set agrees easily
  right   WITHIN-set mean pairwise relatedness -- ⭐ THE MANIPULATION, and the only
                                                 thing that differs

⇒ left to right: the two set types are indistinguishable in capture, on both its
level and its within-set spread, and differ only in how related their members are.
That is the whole design in one figure.

⚠ Only TWO series are drawn (Justin, 2026-09-14).  A third, plain unmatched random,
was distracting: the reader needs to see the two groups are even, not how far apart
they would have been.  The unmatched numbers are kept in the caption, because they
are the evidence that the matching was necessary rather than decorative.

⚠ Capture is measured on the A tapes only -- the same half that defines relatedness
and disjoint from the B tapes where unanimity is scored -- so the matching variable
never touches the outcome.

Usage: 85_fig_setdiagnostics.py [arm] [k]      (default Subclone 5)
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
NEUT2 = "#8f8e86"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

d = json.loads((RES / f"unanimity_{ARM}_mincl100.json").read_text())
B = d["bins"]
SERS = (("rel", f"sets of {KUSE} close relatives", COL[ARM], 2.2, "-"),
        ("mat", "capture-matched random sets", "#111111", 1.5, (0, (4, 2))))
PAN = (("hist_cap", "cap", "mean capture of the set\n(A-tapes recovered per cell)"),
       ("hist_sd", "sd", "spread of capture WITHIN the set\n(sd of A-tapes recovered)"),
       ("hist_rel", "rel", "mean pairwise relatedness WITHIN the set\n(shared prefix depth per A-tape)"))
fig, axs = plt.subplots(1, 3, figsize=(11.4, 3.7))
M = {}
for ax, (key, bk, xlab) in zip(axs, PAN):
    e = np.array(B[bk]); ctr = 0.5 * (e[1:] + e[:-1])
    for s, lab, col, lw, ls in SERS:
        h = np.array(d[key][KUSE][s], float)
        if h.sum() == 0:
            continue
        ax.plot(ctr, 100 * h / h.sum(), color=col, lw=lw, ls=ls,
                label=lab if key == "hist_cap" else None, zorder=3)
    M[bk] = {s: float(np.sum(ctr * np.array(d[key][KUSE][s], float)) /
                      max(np.sum(d[key][KUSE][s]), 1)) for s in ("rel", "mat", "rnd")}
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel("% of sets" if key == "hist_cap" else "")
    w = np.maximum(np.array(d[key][KUSE]["mat"], int) + np.array(d[key][KUSE]["rel"], int), 0)
    lo, hi = np.percentile(np.repeat(ctr, w), [0.2, 99.8])
    ax.set_xlim(max(e[0], lo - 0.06 * (hi - lo)), min(e[-1], hi + 0.06 * (hi - lo)))
axs[0].legend(frameon=False, fontsize=8.6, loc="upper left", handlelength=1.8)
cap = (
    f"{NICE[ARM]} · $k$ = {KUSE} · distributions over SETS, not cells. The matched control replaces "
    f"every cell of a related set by one drawn within a narrow window of its own capture rank.\n"
    f"Left and middle: the two set types are equal in capture — mean {M['cap']['rel']:.1f} vs "
    f"{M['cap']['mat']:.1f} A-tapes recovered, within-set spread {M['sd']['rel']:.2f} vs "
    f"{M['sd']['mat']:.2f}. Right: they differ in relatedness, {M['rel']['rel']:.2f} vs "
    f"{M['rel']['mat']:.2f}.\n"
    f"The middle panel matters most: what drives all-{KUSE} unanimity is within-set HOMOGENEITY of "
    f"capture, since a uniformly badly-captured set agrees easily, so matching the mean alone would "
    f"not have controlled it.\n"
    f"Matching was necessary rather than decorative — unmatched random sets would average "
    f"{M['cap']['rnd']:.1f} A-tapes with spread {M['sd']['rnd']:.2f}. "
    f"Capture is measured on the A tapes only, disjoint from the B tapes where unanimity is scored."
)
fig.text(0.5, -0.06, cap, ha="center", va="top", fontsize=8.0, color=INK2, linespacing=1.65)
fig.tight_layout()
fig.savefig(FIG / f"fig6d_setdiag_{ARM}_k{KUSE}.png"); plt.close(fig)
print(f"wrote figures/fig6d_setdiag_{ARM}_k{KUSE}.png")
for bk in ("cap", "sd", "rel"):
    print(f"  {bk:>4}: rel {M[bk]['rel']:.2f}  mat {M[bk]['mat']:.2f}  (unmatched {M[bk]['rnd']:.2f})")
