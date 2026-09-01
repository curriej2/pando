#!/usr/bin/env python3
r"""
FIGURE 3, PANEL b -- the tape axis is large, and it is a fixed property of the tape.

MARGIN.  Panel a summed the rows of the (cell x tape) matrix; this sums the
columns:  R_t = sum_c Y_ct,  beta_t = R_t / n.  Under the same coin-flip null all
166 tapes share one p, so R_t ~ Bin(n, p) and Var(beta_t) = p(1-p)/n -- for Mouse1
an sd of 0.0044, i.e. every tape within ~1.3 points of 60.6%.

WHY NOT COMPARE VIFs.  VIF = 1 + (m-1) rho, and m is 166 for a cell but n = 12,232
for a tape, so the tape VIF is inflated by the shape of the experiment.  Only the
rho's compare: rho_cell 0.131 vs rho_tape 0.250.

RELIABILITY.  Measure the same 166 tapes in two independent libraries,
X = beta + e1, Y = beta + e2 with independent noise.  Then

    Cov(X, Y) = Var(beta) = sigma_beta^2         (noise cancels -- unbiased)
    corr(X, Y) = sigma_beta^2 / sqrt((sigma_beta^2 + v1)(sigma_beta^2 + v2))

so r is the fraction of the observed between-tape spread that is real.  Under the
null sigma_beta = 0 and r = 0.  We report the measured r against the r predicted
from counting noise alone (v_i = p_i(1-p_i)/n_i); the shortfall is the
library-specific component, which we convert back to an sd in probability units.

Output: figures/fig3b_tape_axis.png, results/fig3b_tape_axis.json
"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARM = "Mouse1"
OTHER = ["Mouse2", "Mouse3", "Initial", "Subclone"]
# Categorical slots 1,3,4,5,7 of the dataviz reference palette, in this order.
# Orange (slot 2) is skipped because it is reserved for the coin-flip null in
# panel a, and it fails the normal-vision floor against magenta (dE 12.9) and
# yellow (13.7); red (slot 8) fails it worse (7.1).  With orange unavailable to
# the series, panel b's null is drawn as a NEUTRAL band instead -- here it is
# reference furniture, not a competing series, so orange never means two things.
# Validated (adjacent pairlist, the documented one for line charts, light mode):
# CVD dE 9.1, normal-vision dE 19.6, lightness and chroma pass.  Aqua/yellow/
# magenta sit below 3:1 on the light surface -> relief rule, met by the legend
# labels plus the per-arm table in README.md.
SERIES = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
          "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
NULLFILL, NULLLINE = "#cfcec6", "#8f8e86"
REP = ("Initial", "Initial_1", "Initial_2")      # true replicate libraries

rates, P, Nc = {}, {}, {}
for a in [ARM] + OTHER:
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y = z["recovered"]
    rates[a] = Y.mean(0); P[a] = float(Y.mean()); Nc[a] = Y.shape[0]

n, p, b = Nc[ARM], P[ARM], rates[ARM]
sd_null = np.sqrt(p * (1 - p) / n)
sd_obs = float(b.std(ddof=1))
vif = sd_obs**2 / sd_null**2
rho = (vif - 1) / (n - 1)
sigma_beta_analytic = float(np.sqrt(max(sd_obs**2 - sd_null**2, 0.0)))
out = {"arm": ARM, "n_cells": n, "p": p, "sd_null_rate": float(sd_null),
       "sigma_beta_analytic": sigma_beta_analytic,
       "sd_obs_rate": sd_obs, "vif_tape": float(vif), "rho_tape": float(rho),
       "rate_min": float(b.min()), "rate_max": float(b.max()),
       "decile_lo": float(np.quantile(b, 0.1)), "decile_hi": float(np.quantile(b, 0.9))}
print(f"{ARM}: n={n:,}  p={p:.4f}  per-tape rate sd {sd_obs:.4f} vs null {sd_null:.4f} "
      f"({sd_obs/sd_null:.0f}x)  VIF_tape {vif:,.0f}  rho_tape {rho:.4f}")
print(f"  rates span {b.min():.3f}..{b.max():.3f}, deciles "
      f"{np.quantile(b,0.1):.3f}/{np.quantile(b,0.9):.3f}")
print(f"  noise-corrected between-tape sd  sqrt({sd_obs:.4f}^2 - {sd_null:.4f}^2) "
      f"= {sigma_beta_analytic:.4f}  (noise is {100*sd_null**2/sd_obs**2:.3f}% of the variance)")

# ---- reliability from two independent libraries of the same population
za = np.load(RES / f"dropout_matrix_{REP[0]}.npz", allow_pickle=False)
Ya, sm = za["recovered"], za["sample"]
names = [str(s) for s in za["sample_names"]]
m1, m2 = sm == names.index(REP[1]), sm == names.index(REP[2])
X, Yv = Ya[m1].mean(0), Ya[m2].mean(0)
n1, n2 = int(m1.sum()), int(m2.sum())
p1, p2 = float(Ya[m1].mean()), float(Ya[m2].mean())
r = float(np.corrcoef(X, Yv)[0, 1])
cov = float(np.cov(X, Yv, ddof=1)[0, 1])            # unbiased for sigma_beta^2
v1, v2 = float(X.var(ddof=1)) - cov, float(Yv.var(ddof=1)) - cov   # total noise
c1, c2 = p1 * (1 - p1) / n1, p2 * (1 - p2) / n2                    # counting only
r_pred = cov / np.sqrt((cov + c1) * (cov + c2))
excess = [max(v1 - c1, 0.0), max(v2 - c2, 0.0)]
# The two libraries differ in overall depth (p1 vs p2), which is a global property
# of the library, not a disagreement about WHICH tapes are bad.  Regressing one on
# the other absorbs that shift/scale; the residual sd is the tape-specific
# disagreement, to be compared with the counting noise it must contain.
slope, icpt = np.polyfit(X, Yv, 1)
resid = Yv - (icpt + slope * X)
sd_resid = float(resid.std(ddof=2))
sd_resid_expected = float(np.sqrt(slope**2 * c1 + c2))
out["reliability"] = {
    "arm": REP[0], "libraries": [REP[1], REP[2]], "n": [n1, n2], "p": [p1, p2],
    "r_measured": r, "r_predicted_counting_only": float(r_pred),
    "sigma_beta": float(np.sqrt(cov)),
    "sd_counting": [float(np.sqrt(c1)), float(np.sqrt(c2))],
    "sd_library_excess": [float(np.sqrt(excess[0])), float(np.sqrt(excess[1]))],
    "slope": float(slope), "sd_residual": sd_resid,
    "sd_residual_expected_counting": sd_resid_expected}
print(f"\nreliability, {REP[1]} vs {REP[2]} (n={n1:,}/{n2:,}):")
print(f"  r measured {r:.4f}   predicted from counting noise alone {r_pred:.4f}")
print(f"  sigma_beta (real between-tape sd) {np.sqrt(cov):.4f}")
print(f"  counting-noise sd {np.sqrt(c1):.4f}/{np.sqrt(c2):.4f}   "
      f"library-specific sd {np.sqrt(excess[0]):.4f}/{np.sqrt(excess[1]):.4f}")
print(f"  after regressing out the global library depth (slope {slope:.3f}): "
      f"residual sd {sd_resid:.4f}, of which counting noise alone would give "
      f"{sd_resid_expected:.4f}")

cross = {a: float(np.corrcoef(rates[ARM], rates[a])[0, 1]) for a in OTHER}
out["cross_arm_r"] = cross
print("\ncross-arm r vs " + ARM + ": " + "  ".join(f"{a} {v:.3f}" for a, v in cross.items()))
(RES / "fig3b_tape_axis.json").write_text(json.dumps(out, indent=1))

# ------------------------------------------------------------------- the panel
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE, AQUA = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 14})
fig, ax = plt.subplots(figsize=(9.0, 5.6))

rank = np.arange(1, len(b) + 1)
# Each arm is sorted on its OWN rates: the curves show that the shape is
# universal, not that rank 40 is the same tape everywhere.  The claim that it is
# the same tapes rests on the cross-arm correlations printed above (0.76-0.96),
# not on this plot -- plotting the other arms in Mouse1's order was tried and is
# unreadable.
ax.fill_between(rank, p - 3 * sd_null, p + 3 * sd_null, color=NULLFILL,
                lw=0, zorder=2)
ax.plot(rank, np.full_like(rank, p, dtype=float), color=NULLLINE, lw=1.4, zorder=2)
for a in [ARM] + OTHER:
    ax.plot(rank, np.sort(rates[a]), color=SERIES[a], lw=2.2, zorder=4,
            label=LABEL[a])

ax.text(136, 0.475, "coin-flip null (Mouse 1)\nevery tape within $\\pm$1.3 points",
        color=INK2, fontsize=11, ha="center", va="top")
ax.annotate("", xy=(136, p - 3 * sd_null - 0.004), xytext=(136, 0.495),
            arrowprops=dict(arrowstyle="-|>", color=NULLLINE, lw=1.5))
# legend text stays in ink; the coloured handle beside it carries identity
leg = ax.legend(frameon=False, fontsize=11, loc="upper left",
                bbox_to_anchor=(0.30, 0.34), handlelength=1.7, borderaxespad=0,
                labelspacing=0.45)
for t in leg.get_texts():
    t.set_color(INK)

ax.text(0.985, 0.035,
        "$\\hat\\beta_t = R_t/n$,   null $\\mathrm{Var}=p(1-p)/n$\n"
        f"sd  {sd_obs:.3f}  observed   vs   {sd_null:.4f}  null\n"
        f"$\\rho_{{tape}}={rho:.2f}$   (vs  $\\rho_{{cell}}=0.13$)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11.5, color=INK,
        linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=AXIS, lw=0.8))

ax.set_xlim(0, len(b) + 1); ax.set_ylim(0, 1.0)
ax.set_xlabel("the 166 tapes, ranked by recovery within each arm")
ax.set_ylabel("fraction of cells recovering the tape   $\\hat\\beta_t$")
# title claims only what the panel shows -- reproducibility is now text, not plot
ax.set_title("b   Which tape it is matters more than which cell", loc="left", pad=12)

fig.tight_layout()
fig.savefig(FIG / "fig3b_tape_axis.png", dpi=200, facecolor=SURF)
print("\nwrote figures/fig3b_tape_axis.png")
