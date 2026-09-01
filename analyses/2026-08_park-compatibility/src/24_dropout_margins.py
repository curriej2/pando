#!/usr/bin/env python3
r"""
================================================================================
 Is dropout a coin flip?  --  the two margins of the (cell x tape) matrix
================================================================================

NULL.  "Dropout is random at the entry level": every (cell, tape) instance is
recovered independently with one common probability p.  A cell carries k = 166
tapes, so

    R_c = sum_t Y_ct,   Y_ct ~iid Bern(p)   =>   R_c ~ Bin(k, p),
    Var(R_c) = k p (1-p),      p_hat = (1/nk) sum_ct Y_ct.

MEASURING THE DEPARTURE.  Give each cell its own propensity pi_c with mean p and
variance sigma^2, R_c | pi_c ~ Bin(k, pi_c).  Law of total variance:

    Var(R_c) = E[k pi(1-pi)]        + Var(k pi)
             = k (p - p^2 - sigma^2) + k^2 sigma^2
             = k p (1-p) + k (k-1) sigma^2

so the variance inflation factor is the classical design effect

    VIF_cell = Var(R_c) / [k p (1-p)] = 1 + (k-1) rho_cell,
    rho_cell = sigma^2_cell / [p (1-p)]          (invert: rho = (VIF-1)/(k-1))

rho is the INTRACLASS CORRELATION of two dropout indicators in the same cell:
the excess probability that tape t' is missing given tape t is missing there.

The tape margin is the same algebra with the roles swapped -- R_t over n cells,
VIF_tape = 1 + (n-1) rho_tape.  ONLY THE rho's ARE COMPARABLE ACROSS MARGINS: the
VIFs differ by the (k-1) vs (n-1) multiplier, which has nothing to do with the
strength of the effect.

CAVEATS carried with the numbers, not fixed here:
  * under a genuine two-way structure p_ct = sigma(alpha_c + beta_t) each margin's
    rho is contaminated by the other at second order; the clean split is the
    two-way variance-components fit (later script).
  * rho_cell > 0 is EXPECTED from per-cell sequencing depth alone.  It establishes
    that dropout is non-ignorable, not that it is biological.
  * the paper's cell filter (>=100 tapes Initial/Subclone, >=20 mice) is a
    selection on R_c itself, so every statistic is reported both unfiltered (the
    assay) and filtered (what an inference would actually see).

Outputs: results/dropout_margins.json, figures/fig3b_draft.png
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Initial", "Subclone", "Mouse1", "Mouse2", "Mouse3"]
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
FOCAL = "Mouse1"


def margins(Y):
    """Both margins of a boolean cell x tape recovery matrix."""
    n, k = Y.shape
    p = float(Y.mean())
    out = {"n_cells": n, "n_tapes": k, "p_recovered": p}
    for name, R, m in (("cell", Y.sum(1).astype(float), k),
                       ("tape", Y.sum(0).astype(float), n)):
        var, null = float(R.var(ddof=1)), m * p * (1 - p)
        vif = var / null if null > 0 else float("nan")
        rho = (vif - 1) / (m - 1)
        out[name] = {
            "items": int(R.size), "per_item": int(m),
            "mean": float(R.mean()), "sd": float(np.sqrt(var)),
            "null_sd": float(np.sqrt(null)), "vif": vif, "rho": rho,
            "sigma": float(np.sqrt(max(rho, 0) * p * (1 - p))),
            "quantiles": {str(q): float(np.quantile(R, q))
                          for q in (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)},
            "min": float(R.min()), "max": float(R.max())}
    return out


D = {}
for a in (sys.argv[1:] or ARMS):
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y = z["recovered"]
    keep = Y.sum(1) >= THR[a]
    D[a] = {"all": margins(Y), "filtered": margins(Y[keep]),
            "threshold": THR[a], "dropped_by_filter": int((~keep).sum()),
            "hist": np.bincount(Y.sum(1), minlength=Y.shape[1] + 1).tolist()}
    for tag in ("all", "filtered"):
        m = D[a][tag]
        print(f"{a:9s} {tag:9s} n={m['n_cells']:6,d}  p={m['p_recovered']:.4f}  "
              f"R_c sd {m['cell']['sd']:6.1f} vs null {m['cell']['null_sd']:5.2f}  "
              f"VIF {m['cell']['vif']:8.1f}  rho_cell {m['cell']['rho']:.4f}  |  "
              f"rho_tape {m['tape']['rho']:.4f}", flush=True)
(RES / "dropout_margins.json").write_text(json.dumps(D, indent=1))

# ------------------------------------------------------------------ draft panel
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE, AQUA = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 13.5})
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.0))

# left: focal arm, observed vs the Bernoulli null
h = np.array(D[FOCAL]["hist"], dtype=float)
k = len(h) - 1
x = np.arange(k + 1)
m = D[FOCAL]["all"]
p = m["p_recovered"]
axL.bar(x, h / h.sum(), width=1.0, color=BLUE, linewidth=0)
from math import lgamma
lb = (np.array([lgamma(k + 1) - lgamma(i + 1) - lgamma(k - i + 1) for i in x])
      + x * np.log(p) + (k - x) * np.log1p(-p))
axL.plot(x, np.exp(lb), color=ORANGE, lw=2.4)
axL.axvline(THR[FOCAL], color=INK2, lw=1.4, ls=(0, (4, 3)))
axL.text(THR[FOCAL] + 2, axL.get_ylim()[1] * 0.94, f"QC cut\n$\\geq${THR[FOCAL]}",
         color=INK2, fontsize=11, va="top")
axL.text(0.97, 0.94,
         f"observed sd {m['cell']['sd']:.0f} tapes\n"
         f"null sd {m['cell']['null_sd']:.1f} tapes\n"
         f"VIF {m['cell']['vif']:.0f}   $\\rho$ = {m['cell']['rho']:.3f}",
         transform=axL.transAxes, ha="right", va="top", fontsize=12, color=INK,
         bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))
axL.set_xlabel("tapes recovered per cell   $R_c$   (of 166)")
axL.set_ylabel("fraction of cells")
axL.set_title(f"a   {LABEL[FOCAL]}: observed vs a coin flip on entries",
              loc="left", pad=10)

# right: every arm
for a, col in zip(ARMS, [BLUE, AQUA, ORANGE, "#8a6fd4", "#c02f4b"]):
    h = np.array(D[a]["hist"], dtype=float)
    axR.plot(np.arange(len(h)), h / h.sum(), color=col, lw=2.0, label=LABEL[a])
axR.legend(frameon=False, fontsize=11)
axR.set_xlabel("tapes recovered per cell   $R_c$")
axR.set_ylabel("fraction of cells")
axR.set_title("b   The same shape in every arm", loc="left", pad=10)

fig.tight_layout()
fig.savefig(FIG / "fig3b_draft.png", dpi=200, facecolor=SURF)
print("wrote figures/fig3b_draft.png")
