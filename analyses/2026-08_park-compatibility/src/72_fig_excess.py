#!/usr/bin/env python3
r"""
Fig 4d -- how far past chance the whole catalogue sits, five arms.

THE IDEA.  Under H0 the p-values of the testable combos are uniform on (0,1), so
the expected number with -log10(p) >= t is exactly

    E_0(t) = N_test * 10^(-t)

a CLOSED FORM.  On a log y-axis that is a straight line of slope -1, fixed by one
number per arm (N_test).  Plotting the observed survival count against it makes
the null explicit instead of hiding it inside a p-value -- which is the standing
guidance in this project ("lead with the comparison, not the p-value").

⚑ Two facts that make the panel self-explaining:
   * the Bonferroni cut at one expected false positive is p <= 1/N_test, i.e.
     t = log10(N_test) -- EXACTLY where each arm's null line crosses y = 1.  So
     the threshold is not an extra annotation, it is a readable feature of the
     figure.
   * the vertical gap between an arm's observed curve and its own null line, read
     at that crossing, IS the excess factor quoted in the text.

WHAT IS PLOTTED, per arm:
   observed  #{testable combos with -log10(p) >= t}, from the binned plane
   null      N_test * 10^(-t)
   marker    at t = log10(N_test), on the observed curve

⚠ These are COMBOS, not events.  One inherited loss is detected through many
overlapping (anchor, depth) clades -- the collapse is 25-121x -- so this panel
answers "is there a population that chance cannot produce?", NOT "how many
distinct losses are there?".  The event counts live in exact_{arm}_d6.json.
⚠ Only testable combos enter (sigma^2 > 0 and m < n); the testable share is
53-100% across arms and is clone-size structure, so N_test differs accordingly.

Palette: five-arm categorical slots, with the null NEUTRAL -- the house rule,
because orange fails the normal-vision floor against magenta and yellow.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Subclone", "Mouse2", "Mouse1", "Mouse3", "Initial"]
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

fig, ax = plt.subplots(figsize=(7.4, 5.2))
TMAX = 60.0
tgrid = np.linspace(0, TMAX, 400)
rows = []
for a in ARMS:
    z = np.load(RES / f"exact_{a}_d6_plane.npz", allow_pickle=False)
    H = z["H"].sum(axis=(0, 1))                 # counts per -log10(p) bin
    left = z["nlp_edges"][:-1]                  # bin left edges; last bin is [60, inf)
    surv = H[::-1].cumsum()[::-1].astype(float)  # #{combos with -log10(p) >= left}
    NT = int(z["n_testable"]); thr = float(np.log10(NT))
    # ⚠ read the FIRST bin at or above the threshold, not the bin containing it:
    # bins are 0.5 wide in t, and the containing bin adds combos below the cut
    # (9% on Subclone).  This matches the counts the catalogue reports.
    obs_at = float(surv[np.searchsorted(left, thr, side="left")])
    ax.plot(tgrid, NT * 10.0 ** (-tgrid), color=NEUT, lw=1.4, zorder=1)
    ax.step(left, np.where(surv > 0, surv, np.nan), where="post",
            color=COL[a], lw=2.4, zorder=3)
    ax.plot([thr], [obs_at], "o", ms=8, mfc=SURFACE, mec=COL[a], mew=2.2, zorder=5)
    rows.append((a, NT, thr, obs_at))
ax.axhline(1.0, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
ax.annotate("one expected false positive", xy=(1.0, 1.35), fontsize=8.5, color=INK2)
ax.annotate("expected under $H_0$\n(uniform $p$: slope $-1$)", xy=(9.6, 1.5e-2),
            fontsize=8.5, color=NEUT2, ha="left")
lo, hi = min(r[2] for r in rows), max(r[2] for r in rows)
ax.axvspan(lo, hi, color=NEUT, alpha=0.35, zorder=0)
ax.annotate("Bonferroni cuts\n$t=\\log_{10}N_{\\rm test}$", xy=((lo + hi) / 2, 4e7),
            ha="center", fontsize=8.5, color=INK2)
ax.set_yscale("log"); ax.set_xlim(0, TMAX); ax.set_ylim(3e-3, 3e8)
ax.set_xlabel(r"$t$  =  $-\log_{10} p$")
ax.set_ylabel(r"combos with $-\log_{10} p \geq t$")
ax.set_title("What chance would produce, and what is there", pad=22,
             fontsize=11.5, fontweight="bold", loc="left")
ax.annotate("marker = each arm's Bonferroni cut, where its own null line crosses 1",
            xy=(0, 1.012), xycoords="axes fraction", fontsize=8.5, color=INK2)
h = [plt.Line2D([], [], color=COL[a], lw=2.4) for a in ARMS]
lab = [f"{NICE[a]}  ({r[3]:,.0f} vs 1)" for a, r in zip(ARMS, rows)]
h.append(plt.Line2D([], [], color=NEUT, lw=1.4)); lab.append("null, one line per arm")
ax.legend(h, lab, frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2,
          handlelength=1.8, bbox_to_anchor=(0.995, 0.985))
ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
fig.savefig(FIG / "fig4d_excess.png"); plt.close(fig)
print(f"{'arm':>9}{'N_test':>13}{'threshold t':>13}{'observed':>12}{'expected':>10}{'excess':>13}")
for a, NT, thr, obs in rows:
    print(f"{NICE[a]:>9}{NT:>13,}{thr:>13.2f}{obs:>12,.0f}{1.0:>10.2f}{obs:>12,.0f}x")
print("wrote figures/fig4d_excess.png")
