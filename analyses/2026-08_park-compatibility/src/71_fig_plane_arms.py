#!/usr/bin/env python3
r"""
The detection plane, all five arms.

x = pi_hat = k/m, how complete the loss is.
y = -log10(p), how surprising it is.  p = max(normal, hypergeometric) from the
    exact finite-population permutation moments.

⚑ One facet per ARM, because the Bonferroni cut is EFALSE / N_test and N_test
differs per arm -- so the horizontal line in each facet IS that arm's decision
rule, the same principle that made the clade-size faceting work.  A line drawn
across pooled arms would be wrong for every one of them.

⚠ Supersedes fig5a/fig5b, which put the MARGIN over a permutation maximum on y.
That statistic was retired with the margin route; do not mix the two figures.

Only TESTABLE combos are drawn (sigma^2 > 0 and m < n).  The untestable share is
printed per facet because it varies 53-100% across arms and is clone-size
structure, not noise -- it caps how much any arm can show.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Initial", "Mouse1", "Mouse2", "Mouse3", "Subclone"]
NICE = {"Initial": "Pre-TX", "Mouse1": "Mouse 1", "Mouse2": "Mouse 2",
        "Mouse3": "Mouse 3", "Subclone": "Subclone"}
YMAX = float(sys.argv[sys.argv.index("--ymax") + 1]) if "--ymax" in sys.argv else 60.0

SURFACE, INK, INK2, GRID, NULL = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#eb6834"
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seqblue", RAMP)
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

Z = {a: np.load(RES / f"exact_{a}_d6_plane.npz", allow_pickle=False) for a in ARMS}
import json
J = {a: json.load(open(RES / f"exact_{a}_d6.json")) for a in ARMS}
PI = Z[ARMS[0]]["pi_edges"]; NLP = Z[ARMS[0]]["nlp_edges"]
keep = NLP[:-1] < YMAX
ymesh = np.append(NLP[:-1][keep], YMAX)
vmax = max(Z[a]["H"].sum(0)[:, keep].max() for a in ARMS)

fig, axes = plt.subplots(1, len(ARMS), figsize=(2.55 * len(ARMS) + 1.4, 4.2),
                         sharey=True, sharex=True)
X, Y = np.meshgrid(PI, ymesh, indexing="ij")
for ax, a in zip(axes, ARMS):
    H = Z[a]["H"].sum(0)[:, keep].astype(float)
    thr = float(Z[a]["nlp_threshold"])
    pc = ax.pcolormesh(X, Y, np.ma.masked_where(H <= 0, H), cmap=CMAP,
                       norm=LogNorm(vmin=1, vmax=vmax), shading="flat", rasterized=True)
    ax.axhline(thr, color=NULL, lw=2.0, zorder=4)
    ax.annotate(f"{thr:.1f}", xy=(0.03, thr), xytext=(0, 4), textcoords="offset points",
                color=NULL, fontsize=8.5, fontweight="bold")
    d = J[a]
    ax.set_title(f"{NICE[a]}\n{d['n_events']:,} events · {100*d['frac_missing_in_events']:.1f}% "
                 f"of missing", fontsize=9.5, color=INK, loc="center")
    ax.annotate(f"{100*d['n_testable']/d['n_valid']:.0f}% testable",
                xy=(0.97, 0.97), xycoords="axes fraction", ha="right", va="top",
                fontsize=8, color=INK2)
    ax.set_xlim(0, 1); ax.set_ylim(0, YMAX); ax.set_xticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=8.5)
axes[0].set_ylabel(r"$-\log_{10} p$")
fig.supxlabel(r"$\hat\pi = k/m$   (fraction of the clade missing the tape)",
              fontsize=10, y=0.01)
fig.suptitle("Heritable tape loss, five arms — complete, and far past the threshold",
             fontsize=11.5, fontweight="bold", x=0.055, ha="left", y=1.10)
fig.text(0.055, 1.045, "each combo is one (clade, tape); orange = that arm's Bonferroni cut "
         "at one expected false positive", fontsize=8.5, color=INK2, ha="left")
cb = fig.colorbar(pc, ax=list(axes), pad=0.012, fraction=0.022)
cb.set_label("combos per bin", color=INK2, fontsize=9)
cb.outline.set_edgecolor(GRID); cb.ax.tick_params(colors=INK2, labelsize=8)
fig.savefig(FIG / "fig4d_plane_arms.png"); plt.close(fig)
for a in ARMS:
    H = Z[a]["H"].sum(0); thr = float(Z[a]["nlp_threshold"])
    above = H[:, NLP[:-1] >= thr].sum()
    hi = H[:, NLP[:-1] >= YMAX].sum()
    print(f"  {NICE[a]:>9}: threshold -log10(p) = {thr:5.2f}; {above:>9,.0f} combos above; "
          f"{hi:>8,.0f} beyond the y-limit ({100*hi/max(above,1):.1f}% of those above)")
print("wrote figures/fig4d_plane_arms.png")
