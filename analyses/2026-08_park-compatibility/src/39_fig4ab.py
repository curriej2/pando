#!/usr/bin/env python3
r"""
FIGURE 4, PANELS a and b -- related cells lose the same tapes.

(a) IS IT HERITABLE, MEASURED AGAINST SOMETHING WE KNOW IS?
    rho depends on a feature's marginal frequency -- Var(pi) <= p(1-p), so a
    feature present in 97% of cells has 8x less room to vary between clones than
    one at 40%.  A control at an unmatched frequency is therefore meaningless.
    So the control is a CURVE: the same statistic on "tape z has reached site L",
    L = 2..6, whose marginals sweep 0.97 down to 0.35 and overlap missingness.
    Missingness is read against its own arm's curve at its own frequency.

    ⚑ Plotted on the EQUAL-CLONE-WEIGHTED estimator, not the pooled one.  The
    pooled rho is a ratio of sums, so clones enter weighted by n_C(n_C-1): the
    largest clone holds 98.9% of Mouse2's weight and 84.0% of Mouse1's, and a
    huge clone spans so much of the tree that within-clone concordance is diluted.
    Pooling therefore HID the signal (Mouse2 +0.009 pooled vs +0.246 equal-weight).

(b) IS IT ONE CLONE, OR ALL OF THEM?
    The distribution of per-clone excess, one point per clone with >= 5 cells,
    against the control measured identically.  This is the panel that answers
    "which clone is contributing" -- which the pooled statistic cannot.

Excess is always observed minus the mean of a within-sample permutation null that
preserves each clone's size and sample composition.

Output: figures/fig4a_heritable.png, figures/fig4b_perclone.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]
SERIES = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
          "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, ORANGE = "#e1e0d9", "#c3c2b7", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 14})

D = {a: json.loads((RES / f"pool_robustness_{a}.json").read_text()) for a in ARMS}

# ----------------------------------------------------------------- panel a
fig, ax = plt.subplots(figsize=(11.2, 5.6))
ratio = {}; ymax = 0.0
for a in ARMS:
    fam = D[a]["families"]
    ctrl = [(fam[k]["marginal"], fam[k]["equal_weight_excess"])
            for k in fam if k.startswith("depth")]
    ctrl.sort()
    cx, cy = np.array([c[0] for c in ctrl]), np.array([c[1] for c in ctrl])
    mx, my = fam["missing"]["marginal"], fam["missing"]["equal_weight_excess"]
    ax.plot(cx, cy, "-o", color=SERIES[a], lw=1.6, ms=5.5, alpha=0.85,
            markerfacecolor=SURF, markeredgecolor=SERIES[a], zorder=3)
    ax.plot([mx], [my], "D", ms=13, color=SERIES[a], markeredgecolor=SURF,
            markeredgewidth=1.8, zorder=5)
    j = int(np.argmin(np.abs(cx - mx)))          # control at the nearest marginal
    ratio[a] = my / cy[j] if cy[j] > 0 else np.nan
    ax.plot([mx, cx[j]], [my, cy[j]], color=SERIES[a], lw=1.0, ls=(0, (2, 2)),
            alpha=0.6, zorder=2)
    ymax = max(ymax, cy.max(), my)

ax.axhline(0, color=ORANGE, lw=2.0, zorder=4)
ax.text(0.995, 0.012, "permutation null", color=ORANGE, fontsize=11,
        ha="right", va="bottom", transform=ax.get_yaxis_transform())
from matplotlib.lines import Line2D
h = [Line2D([], [], color=SERIES[a], lw=1.6, marker="o", ms=5.5,
            markerfacecolor=SURF, label=f"{LABEL[a]}      {100*ratio[a]:.0f}%")
     for a in ARMS]
h += [Line2D([], [], color=INK2, lw=0, marker="D", ms=11, label="  tape missing"),
      Line2D([], [], color=INK2, lw=1.4, marker="o", ms=5, markerfacecolor=SURF,
             label="  control: tape reached site $L$")]
leg = ax.legend(handles=h, frameon=False, fontsize=10.5, loc="upper left",
                bbox_to_anchor=(1.015, 1.0), handlelength=1.8,
                labelspacing=0.42, borderpad=0.6,
                title="missingness, as a % of the\ncontrol at matched frequency",
                title_fontsize=10.5, alignment="left")
leg.get_title().set_color(INK2)
for t in leg.get_texts():
    t.set_color(INK)
ax.set_xlim(0, 1.02); ax.set_ylim(-0.04, ymax * 1.06)
ax.set_xlabel("marginal frequency of the feature")
ax.set_ylabel("excess within-clone correlation\n(observed $-$ permutation null)")
ax.set_title("a   Related cells lose the same tapes", loc="left", pad=30)
ax.text(0, 1.015, "equal weight per clone; $\\rho$ depends on frequency, so the "
        "control is a curve", transform=ax.transAxes, ha="left", va="bottom",
        fontsize=11.5, color=INK2)
fig.tight_layout(); fig.savefig(FIG / "fig4a_heritable.png", dpi=200, facecolor=SURF)
print("wrote fig4a_heritable.png  ratios: " +
      "  ".join(f"{LABEL[a]} {100*ratio[a]:.0f}%" for a in ARMS))

# ----------------------------------------------------------------- panel b
fig, ax = plt.subplots(figsize=(10.6, 5.8))
rng = np.random.default_rng(0)
allv = np.concatenate([np.array(D[a]["families"][f]["clone_excess"])
                       for a in ARMS for f in ("missing", "depth>=6")])
lo, hi = allv.min(), allv.max()
ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.20 * (hi - lo))   # headroom for labels
ytxt = ax.get_ylim()[1] - 0.01 * (hi - lo)
for i, a in enumerate(ARMS):
    for j, (fam, off, fill) in enumerate(
            [("missing", -0.17, True), ("depth>=6", 0.17, False)]):
        e = np.array(D[a]["families"][fam]["clone_excess"])
        sz = np.array(D[a]["families"][fam]["clone_sizes"])
        x = i + off + rng.normal(0, 0.045, e.size)
        ax.scatter(x, e, s=np.clip(3 + sz / 6, 3, 90), alpha=0.30,
                   color=SERIES[a] if fill else INK2, edgecolors="none", zorder=2)
        ax.plot([i + off - 0.13, i + off + 0.13], [np.median(e)] * 2,
                color=INK if fill else INK2, lw=2.6, zorder=4)
        # each arm's largest clone -- the one that dominates the pooled estimate
        b_ = int(np.argmax(sz))
        ax.plot([x[b_]], [e[b_]], "o", ms=13, markerfacecolor="none",
                markeredgecolor=INK, markeredgewidth=1.6, zorder=6)
    n = D[a]["families"]["missing"]["n_clones_tested"]
    pos = 100 * D[a]["families"]["missing"]["frac_clones_positive"]
    ax.text(i, ytxt, f"{n:,} clones\n{pos:.0f}% positive", ha="center",
            va="top", fontsize=10, color=INK2)
ax.axhline(0, color=ORANGE, lw=2.0, zorder=3)
ax.text(0.995, 0.01, "permutation null", color=ORANGE, fontsize=11, ha="right",
        va="bottom", transform=ax.get_yaxis_transform())
ax.set_xticks(range(len(ARMS)))
ax.set_xticklabels([LABEL[a] for a in ARMS], fontsize=11.5)
ax.set_xlim(-0.55, len(ARMS) - 0.45)
ax.set_ylabel("per-clone excess correlation")
ax.set_title("b   …in essentially every clone, not one big one", loc="left", pad=30)
ax.text(0, 1.015, "one clone per point ($\\geq$5 cells), sized by cell count  ·  "
        "colour = tape missing, grey = control  ·  $\\bigcirc$ = largest clone",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=11, color=INK2)
fig.tight_layout(); fig.savefig(FIG / "fig4b_perclone.png", dpi=200, facecolor=SURF)
print("wrote fig4b_perclone.png")
