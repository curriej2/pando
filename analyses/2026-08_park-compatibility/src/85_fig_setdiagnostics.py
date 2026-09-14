#!/usr/bin/env python3
r"""
COMPANION TO FIG C -- what the matching equalised, and what it did not.

Three panels, all distributions over SETS (not cells), from `83`:
  left    set MEAN capture          -- matched by construction, so rel and mat
                                       coincide; plain random sits apart
  middle  WITHIN-set sd of capture  -- the quantity that actually drives unanimity,
                                       since a uniformly badly-captured set reaches
                                       agreement easily.  Matching equalises this too.
  right   WITHIN-set mean pairwise relatedness -- ⭐ THE MANIPULATION.  rel is far
                                       to the right of mat and rnd, which coincide.

⇒ read left to right: the two set types are indistinguishable in capture, on both
its level and its within-set spread, and differ only in how related their members
are.  That is the whole design in one figure.

⚠ Capture is measured on the A tapes only -- the same half used for relatedness and
disjoint from the B tapes where unanimity is scored, so the matching variable never
touches the outcome.

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
NEUT, NEUT2 = "#cfcec6", "#8f8e86"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

d = json.loads((RES / f"unanimity_{ARM}_mincl100.json").read_text())
B = d["bins"]
SERS = (("rel", f"sets of {KUSE} close relatives", COL[ARM], 2.2, "-"),
        ("mat", "capture-matched random sets", "#111111", 1.5, (0, (4, 2))),
        ("rnd", "plain random sets", NEUT2, 1.5, "-"))
PAN = (("hist_cap", "cap", "mean capture of the set\n(A-tapes recovered per cell)",
        "matched by construction"),
       ("hist_sd", "sd", "spread of capture WITHIN the set\n(sd of A-tapes recovered)",
        "matched too — and this is what drives unanimity"),
       ("hist_rel", "rel", "mean pairwise relatedness WITHIN the set\n(shared prefix depth per A-tape)",
        "the manipulation: this is the only thing that differs"))
fig, axs = plt.subplots(1, 3, figsize=(11.4, 3.9))
for ax, (key, bk, xlab, note) in zip(axs, PAN):
    e = np.array(B[bk]); ctr = 0.5 * (e[1:] + e[:-1])
    for s, lab, col, lw, ls in SERS:
        h = np.array(d[key][KUSE][s], float)
        if h.sum() == 0:
            continue
        ax.plot(ctr, 100 * h / h.sum(), color=col, lw=lw, ls=ls,
                label=lab if key == "hist_cap" else None, zorder=3)
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel("% of sets" if key == "hist_cap" else "")
    ax.set_title(note, fontsize=9.0, color=INK2, pad=7, loc="left")
    lo, hi = np.percentile(np.repeat(ctr, np.maximum(
        np.array(d[key][KUSE]["rnd"], int) + np.array(d[key][KUSE]["rel"], int), 0)), [0.2, 99.8])
    ax.set_xlim(max(e[0], lo - 0.06 * (hi - lo)), min(e[-1], hi + 0.06 * (hi - lo)))
axs[0].legend(frameon=False, fontsize=8.4, loc="upper left", handlelength=1.8)
fig.suptitle(f"{NICE[ARM]} · k = {KUSE} · capture measured on the A tapes only, "
             f"disjoint from the B tapes where unanimity is scored",
             fontsize=8.8, color=INK2, y=1.04)
fig.tight_layout()
fig.savefig(FIG / f"fig6d_setdiag_{ARM}_k{KUSE}.png"); plt.close(fig)
print(f"wrote figures/fig6d_setdiag_{ARM}_k{KUSE}.png")
for key, bk, _, _ in PAN:
    e = np.array(B[bk]); ctr = 0.5 * (e[1:] + e[:-1])
    ms = []
    for s, *_ in SERS:
        h = np.array(d[key][KUSE][s], float)
        ms.append(np.sum(ctr * h) / max(h.sum(), 1))
    print(f"  {bk:>4} mean: rel {ms[0]:.2f}  mat {ms[1]:.2f}  rnd {ms[2]:.2f}")
