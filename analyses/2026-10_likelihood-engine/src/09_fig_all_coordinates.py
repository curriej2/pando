"""Teaching figure for Betancourt §1.2-1.4, "intuition 3": being near the mode means being near the
centre in EVERY coordinate at once.

a  P(all D coordinates within +-c SD of the mode) = p_c^D for a standard Gaussian, c = 1, 2, 3,
   with Monte Carlo points where the probability is measurable.
b  One draw at D = 400, coordinate by coordinate: nothing unusual in any single coordinate.
c  The same draw's distance from the mode as coordinates are added: r(k) = sqrt(sum_{i<=k} q_i^2)
   grows like sqrt(k) -- many ordinary offsets add up to a large distance.
Notes: sciphy_notes.md §G.4.11.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, MUTED, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#a3a29c", "#e4e3df", "#fcfcfb", "#efeeea"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
rng = np.random.default_rng(13102026)
DMAX, DDRAW = 1000, 400

fig = plt.figure(figsize=(16.5, 5.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.45, 1.05], wspace=0.36)

# ---------------------------------------------------------------- a: p_c^D
ax = fig.add_subplot(gs[0])
D = np.arange(1, DMAX + 1)
ax.axvspan(200, 800, color=BAND, zorder=0)
ax.text(515, 1e-77, "our fixed-topology\nposterior: D ≈ 200–800", ha="left", va="bottom", fontsize=8.5,
        color=INK2)
ax.axhline(1e-6, color=MUTED, lw=0.8, ls=":")
ax.text(DMAX, 1.6e-6, "one in a million", ha="right", va="bottom", fontsize=8, color=INK2)
d1pct = {}
for c, col in zip([1, 2, 3], C):
    p = stats.norm.cdf(c) - stats.norm.cdf(-c)
    y = p ** D
    ax.plot(D, y, color=col)
    d1pct[c] = np.log(0.01) / np.log(p)
    yl = y[-1] if y[-1] > 1e-80 else None
    lab = f"±{c} SD  (each coordinate: {p:.1%})"
    if c == 1:
        ax.text(170, 1e-27, lab, color=col, fontsize=9, fontweight="bold", va="bottom")
    elif c == 2:
        ax.text(170, 1e-20, lab, color=col, fontsize=9, fontweight="bold", va="bottom")
    else:
        ax.text(DMAX * 0.98, y[-1] * 6, lab, color=col, fontsize=9, fontweight="bold", ha="right", va="bottom")
# Monte Carlo check for +-1 SD where measurable
NMC = 2_000_000
Dmc = [1, 2, 5, 10, 15, 20, 25]
z = np.abs(rng.standard_normal((NMC, max(Dmc)))) < 1
mc = [np.mean(z[:, :d].all(axis=1)) for d in Dmc]
ax.plot(Dmc, mc, "o", color=INK, ms=5, zorder=5, label=f"Monte Carlo, {NMC:,} draws (±1 SD)")
ax.set_yscale("log")
ax.set_ylim(1e-80, 3)
ax.set_xlim(0, DMAX)
ax.set_xlabel("D = number of parameters")
ax.set_ylabel("P(every coordinate within the box)")
ax.set_title("a  Near the mode in every coordinate", loc="left", pad=30)
ax.text(0, 1.02, "probability that one draw falls inside a box of\n±c SD around the mode, in all D coordinates", transform=ax.transAxes,
        fontsize=8.5, color=INK2, va="bottom")
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), frameon=False, fontsize=8)
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------------------------------------------------------------- b: one draw, coordinate by coordinate
q = rng.standard_normal(DDRAW)
out = np.abs(q) > 1
axb = fig.add_subplot(gs[1])
axb.axhspan(-1, 1, color=BAND, zorder=0)
axb.axhline(0, color=INK2, lw=1)
idx = np.arange(1, DDRAW + 1)
axb.plot(idx[~out], q[~out], "o", color=MUTED, ms=3.2, mew=0)
axb.plot(idx[out], q[out], "o", color=INK, ms=3.2, mew=0)
axb.text(DDRAW + 4, 0, "the mode\n(every\ncoordinate 0)", fontsize=8, color=INK2, va="center")
axb.text(DDRAW + 4, 1, "±1 SD", fontsize=8, color=INK2, va="center")
axb.set_xlim(0, DDRAW + 60)
axb.set_ylim(-4, 4)
axb.set_xlabel("coordinate i (one draw, D = 400)")
axb.set_ylabel("$q_i$  (SDs from the mode)")
axb.set_title("b  One typical draw (D = 400), coordinate by coordinate", loc="left", pad=30)
axb.text(0, 1.02, f"{out.sum()} of {DDRAW} coordinates outside ±1 SD (expected {0.3173*DDRAW:.0f} ± "
         f"{np.sqrt(DDRAW*0.3173*0.6827):.0f})\nlargest |q_i| = {np.abs(q).max():.2f}: none extreme, yet outside even the ±3 SD box",
         transform=axb.transAxes, fontsize=8.5, color=INK2, va="bottom")
axb.grid(axis="y", color=GRID, lw=0.8); axb.set_axisbelow(True)
axb.spines["bottom"].set_position(("data", -4))

# ---------------------------------------------------------------- c: distance accumulates
axc = fig.add_subplot(gs[2])
k = np.arange(1, DDRAW + 1)
lo, hi = stats.chi.ppf(0.05, k), stats.chi.ppf(0.95, k)
axc.fill_between(k, lo, hi, color=BAND, lw=0)
axc.plot(k, np.sqrt(k), color=INK2, lw=1.2, ls="--")
r = np.sqrt(np.cumsum(q ** 2))
axc.plot(k, r, color=INK, lw=1.8)
axc.text(DDRAW * 0.97, np.sqrt(DDRAW) + 1.3, "this draw", ha="right", fontsize=9, color=INK, fontweight="bold")
axc.text(DDRAW * 0.55, np.sqrt(DDRAW * 0.55) - 3.2, "√k  (shaded: 90% of draws)", fontsize=8.5,
         color=INK2, ha="left")
axc.set_xlim(0, DDRAW); axc.set_ylim(0, 24)
axc.set_xlabel("k = coordinates included so far")
axc.set_ylabel("distance from the mode, r(k)  (SDs)")
axc.set_title("c  The same draw: small offsets add up", loc="left", pad=30)
axc.text(0, 1.02, f"r(k) = √(q₁² + … + q_k²) = distance from the mode\nusing the first k coordinates; at k = {DDRAW}: {r[-1]:.1f} SDs",
         transform=axc.transAxes, fontsize=8.5, color=INK2, va="bottom")
axc.grid(color=GRID, lw=0.8); axc.set_axisbelow(True)

fig.suptitle("Why no draw lands near the mode: it would have to be near the centre in every coordinate at once",
             x=0.06, y=1.04, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.975, "Standard Gaussian, D independent coordinates.  Reading rule: in a, a curve near the "
         "bottom means 'essentially never'; b and c are the same single draw.", fontsize=9, color=INK2)
fig.savefig(OUT / "fig_volume_3_allcoords.png", dpi=200, bbox_inches="tight")
plt.close(fig)

for c in (1, 2, 3):
    p = stats.norm.cdf(c) - stats.norm.cdf(-c)
    print(f"±{c} SD: per-coordinate {p:.4f};  P(all) at D=10 {p**10:.3g}, 100 {p**100:.3g}, 200 {p**200:.3g}, "
          f"400 {p**400:.3g}, 800 {p**800:.3g};  drops below 1% at D = {d1pct[c]:.0f}")
print("MC (±1 SD):", dict(zip(Dmc, np.round(mc, 5))), " exact:", {d: round(0.682689**d, 5) for d in Dmc})
print(f"draw: {out.sum()} outside ±1 SD, max |q| {np.abs(q).max():.2f}, r(400) {r[-1]:.2f}")
print("wrote", OUT / "fig_volume_3_allcoords.png")
