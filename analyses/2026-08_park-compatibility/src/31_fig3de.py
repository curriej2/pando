#!/usr/bin/env python3
r"""
FIGURE 3, PANELS d and e -- is dropout INFORMATIVE about what the tape recorded?

These two are one argument and must be read together, so they share a y-axis:
mean edit depth, in sites, of the tapes we DID recover.  If the axis along which
data go missing is also an axis along which the recorder read less, then the tapes
we fail to see are not a random sample of the tapes that exist, and every
time-like quantity downstream (lambda, branch lengths, node ages) inherits the
bias.

(d) PER TAPE.  x = fraction of cells recovering tape z; y = its mean depth.
    Rising => a tape that is hard to see is also a tape that recorded less.  This
    is row A9's prediction: pegRNA and TargetBC are co-integrated, so a closed or
    silenced locus is simultaneously unreadable and unedited -- one latent cause,
    two symptoms.
    CONFOUND: a low-recovery tape is seen mostly in HIGH-capture cells.  The open
    circles recompute y using only cells with R_c >= 140, where every tape is
    scored in a comparable cell population.

(e) PER CELL.  x = R_c; y = mean depth over a FIXED reference set of the 30
    easiest tapes, identical for every cell -- without that control a low-R_c cell
    would be scored on its easy tapes only, which is a different question.
    Flat => the cell axis carries no information about editing, so it can be
    marginalised rather than modelled.

Both panels are Mouse 1.  rho values are checked against 26_dropout_depth.py.

CEILING CAVEAT, stated on the panel: Mouse1 tapes average 4.9 of 6 sites, so
depth has little room left to vary; panel e's flatness is partly attenuation.
Pre-TX, at ~3 of 6, is the arm with dynamic range and gives rho = +0.11.

Output: figures/fig3de_informative.png, results/fig3de.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARM, BLUE = "Mouse1", "#2a78d6"
JUNK, HI_CELL, TOP_TAPES = 2, 140, 30

z = np.load(RES / f"dropout_matrix_{ARM}.npz", allow_pickle=False)
Y, dep, trm = z["recovered"], z["depth"].astype(float), z["term"]
n, k = Y.shape
ok = Y & (trm != JUNK)                       # recovered AND depth not censored
Rc = Y.sum(1)
rate = Y.mean(0)

# ---- (d) per tape
depth_t = dep.sum(0, where=ok) / np.maximum(ok.sum(0), 1)
hi = Rc >= HI_CELL
okh = ok[hi]
depth_hi = np.where(okh.sum(0) > 20, dep[hi].sum(0, where=okh) / np.maximum(okh.sum(0), 1), np.nan)
rho_d = spearmanr(rate, depth_t).statistic
rho_dc = spearmanr(rate[np.isfinite(depth_hi)], depth_hi[np.isfinite(depth_hi)]).statistic
o = np.argsort(rate); q = max(k // 10, 1)
lo_d, hi_d = depth_t[o[:q]].mean(), depth_t[o[-q:]].mean()

# ---- (e) per cell, on the fixed reference set
ref = np.argsort(rate)[::-1][:TOP_TAPES]
okr = ok[:, ref]
nref = okr.sum(1)
depth_c = np.where(nref >= 5, dep[:, ref].sum(1, where=okr) / np.maximum(nref, 1), np.nan)
m = np.isfinite(depth_c)
rho_e = spearmanr(Rc[m], depth_c[m]).statistic
edges = np.arange(20, 168, 15)
bx, by = [], []
for a_, b_ in zip(edges[:-1], edges[1:]):
    sel = m & (Rc >= a_) & (Rc < b_)
    if sel.sum() >= 30:
        bx.append(Rc[sel].mean()); by.append(np.nanmean(depth_c[sel]))

# The decile trend shows the coupling is a THRESHOLD, not a gradient: depth
# climbs steeply up to beta ~0.3 and is flat above it.  So the informative part
# of dropout is confined to a removable minority of tapes -- quantify that.
CUT = 0.30
sel = rate >= CUT
rho_above = spearmanr(rate[sel], depth_t[sel]).statistic
print(f"    threshold: {(~sel).sum()} tapes below beta={CUT} "
      f"({100*(~sel).mean():.0f}%), mean depth {depth_t[~sel].mean():.2f}; "
      f"{sel.sum()} above, mean depth {depth_t[sel].mean():.2f}, "
      f"rho among them {rho_above:+.3f}")

out = {"arm": ARM, "threshold": {"cut": CUT, "n_below": int((~sel).sum()),
                                 "depth_below": float(depth_t[~sel].mean()),
                                 "depth_above": float(depth_t[sel].mean()),
                                 "rho_above_cut": float(rho_above)}, "per_tape": {"rho": float(rho_d), "rho_controlled": float(rho_dc),
                                "decile_lo": float(lo_d), "decile_hi": float(hi_d)},
       "per_cell": {"rho_controlled": float(rho_e), "n_ref_tapes": TOP_TAPES,
                    "binned_x": bx, "binned_y": by}}
(RES / "fig3de.json").write_text(json.dumps(out, indent=1))
print(f"(d) per tape: rho {rho_d:+.3f}  controlled {rho_dc:+.3f}  "
      f"deciles {lo_d:.2f} -> {hi_d:.2f} sites")
print(f"(e) per cell: rho {rho_e:+.3f} on the top-{TOP_TAPES} reference tapes; "
      f"binned depth {by[0]:.3f} -> {by[-1]:.3f} sites")

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 14})
fig, (axD, axE) = plt.subplots(1, 2, figsize=(13.4, 5.4), sharey=True)
ylo, yhi = 1.4, 5.6

axD.plot(rate, depth_t, "o", ms=6.5, color=BLUE, markeredgecolor=SURF,
         markeredgewidth=0.8, alpha=0.9, zorder=3)
axD.plot(rate[np.isfinite(depth_hi)], depth_hi[np.isfinite(depth_hi)], "o", ms=6.5,
         markerfacecolor="none", markeredgecolor=INK2, markeredgewidth=1.0,
         alpha=0.55, zorder=2)
dx, dy = [], []
for i in range(10):
    idx = o[i * k // 10:(i + 1) * k // 10]
    dx.append(rate[idx].mean()); dy.append(depth_t[idx].mean())
axD.plot(dx, dy, "-", color="#eb6834", lw=2.4, zorder=4)
axD.plot(dx[1:-1], dy[1:-1], "o", ms=8, color="#eb6834", markeredgecolor=SURF,
         markeredgewidth=1.5, zorder=5)
axD.plot([dx[0], dx[-1]], [dy[0], dy[-1]], "s", ms=12, color="#eb6834",
         markeredgecolor=SURF, markeredgewidth=1.6, zorder=6)
axD.text(0.40, 2.85, f"decile means: {lo_d:.2f} sites in the worst\n"
         f"decile, $\\approx${depth_t[sel].mean():.1f} everywhere above "
         f"$\\hat\\beta_t\\approx{CUT:.1f}$",
         color="#eb6834", fontsize=11.5, ha="left", va="top")
axD.text(0.975, 0.05, f"Spearman $\\rho={rho_d:+.2f}$   "
         f"({rho_dc:+.2f} on high-capture cells)\n"
         f"but only ${rho_above:+.2f}$ among the {100*sel.mean():.0f}% of tapes "
         f"above $\\hat\\beta_t={CUT:.1f}$",
         transform=axD.transAxes, ha="right", va="bottom", fontsize=11.5, color=INK,
         bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))
axD.set_xlim(-0.03, 1.0); axD.set_ylim(ylo, yhi)
axD.set_xlabel("fraction of cells recovering the tape   $\\hat\\beta_t$")
axD.set_ylabel("mean edit depth of recovered tapes\n(sites of 6)")
axD.set_title("d   A tape we cannot see is a tape that recorded less",
              loc="left", pad=12)

sub = np.where(m)[0]
sub = sub[np.random.default_rng(0).permutation(sub.size)[:4000]]
axE.plot(Rc[sub] + np.random.default_rng(1).normal(0, 0.8, sub.size), depth_c[sub],
         "o", ms=2.6, color=BLUE, alpha=0.12, markeredgecolor="none", zorder=2)
axE.plot(bx, by, "-o", color="#eb6834", lw=2.4, ms=9, markeredgecolor=SURF,
         markeredgewidth=1.6, zorder=4)
axE.text(0.975, 0.05, f"Spearman $\\rho={rho_e:+.2f}$\n"
         f"{by[0]:.2f} $\\rightarrow$ {by[-1]:.2f} sites across the range",
         transform=axE.transAxes, ha="right", va="bottom", fontsize=11.5, color=INK,
         bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))
axE.text(0.03, 0.95, "depth scored on one fixed set of\n30 tapes, the same for every cell",
         transform=axE.transAxes, ha="left", va="top", fontsize=11, color=INK2)
axE.set_xlim(15, 167)
axE.set_xlabel("tapes recovered per cell   $R_c$")
axE.set_title("e   …but a cell we see badly recorded just as much", loc="left", pad=12)

fig.tight_layout()
fig.savefig(FIG / "fig3de_informative.png", dpi=200, facecolor=SURF)
print("wrote figures/fig3de_informative.png")
