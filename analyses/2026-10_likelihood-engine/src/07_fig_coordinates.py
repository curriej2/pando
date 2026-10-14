"""Teaching figures: which posterior summaries depend on the coordinates ("the ruler").

Fig A  lambda ~ Exp(1), summarised in lambda and in u = log(lambda).
Fig B  an illustrative weakly-informed late node (height h ~ Beta(1.5, 6), NOT a fitted
       posterior), summarised in three rulers: calendar time (the engine's point estimate),
       the logit/ratio coordinate x = logit(h) (where the sampler runs, parent at h = 1), and
       recorder depth e = E(day) under Park's in-vitro clock E(t) = A(1 - exp(-t/tau)).

Every summary is computed on a fine grid in its own coordinate (density there = density in h
times |dh/dcoord|), then mapped back to the reporting scale. A Monte Carlo cross-check is printed.
Notes: sciphy_notes.md §G.4.10.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats, special
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)

# palette: reference categorical slots 1-3 (validated all-pairs, light mode)
C = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
MASS = 0.90


def summaries(grid, dens, to_report):
    """Mode, HPD, median, equal-tailed interval, mean-of-coordinate on a grid; mapped to report scale."""
    w = dens * np.gradient(grid)
    w = w / w.sum()
    cdf = np.cumsum(w)
    q = lambda p: grid[np.searchsorted(cdf, p)]
    order = np.argsort(dens)[::-1]                       # HPD: highest-density cells first
    keep = order[: np.searchsorted(np.cumsum(w[order]), MASS) + 1]
    return {
        "mode": to_report(grid[np.argmax(dens)]),
        "hpd": sorted([to_report(grid[keep].min()), to_report(grid[keep].max())]),
        "median": to_report(q(0.5)),
        "et": sorted([to_report(q((1 - MASS) / 2)), to_report(q((1 + MASS) / 2))]),
        "coord_mean": to_report(np.sum(w * grid)),
    }


def strip(ax, rows, rulers, xlab, xlim, note):
    """One row per summary; one dodged mark per ruler; right-hand verdict text."""
    n = len(rulers)
    for i, (key, label) in enumerate(rows):
        y = len(rows) - 1 - i
        vals = [r["s"][key] for r in rulers]
        for j, (r, v) in enumerate(zip(rulers, vals)):
            yy = y + (j - (n - 1) / 2) * 0.22
            if isinstance(v, (list, tuple)):
                ax.plot(v, [yy, yy], color=r["c"], lw=2.5, solid_capstyle="round")
                ax.plot(v, [yy, yy], "|", color=r["c"], ms=9, mew=2)
            else:
                ax.plot(v, yy, "o", color=r["c"], ms=8, mec=SURF, mew=1.5)
        flat = np.ravel(vals)
        same = np.ptp(flat.reshape(n, -1), axis=0).max() < 0.01 * (xlim[1] - xlim[0])
        ax.text(1.01, y, "same in every ruler" if same else "depends on the ruler",
                transform=ax.get_yaxis_transform(), va="center", fontsize=9,
                color=INK2 if same else INK, fontweight="normal" if same else "bold")
        ax.axhline(y - 0.5, color=GRID, lw=0.8, zorder=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for _, lab in rows][::-1])
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlab)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(0, 1.02, note, transform=ax.transAxes, fontsize=9, color=INK2, va="bottom")


def density_panel(ax, x, d, s_native, color, title, xlab, label):
    ax.plot(x, d, color=color)
    lo, hi = s_native["hpd"]
    m = (x >= lo) & (x <= hi)
    ax.fill_between(x[m], 0, d[m], color=color, alpha=0.18, lw=0)
    ax.axvline(s_native["mode"], color=color, lw=1.2, ls="--")
    ax.text(s_native["mode"], d.max() * 1.04, " mode", color=INK2, fontsize=9, va="bottom")
    ax.set_title(title, loc="left", pad=20)
    ax.text(0, 1.03, label, transform=ax.transAxes, fontsize=9, color=INK2, va="bottom")
    ax.set_xlabel(xlab)
    ax.set_ylabel("density (per unit of this axis)")
    ax.set_ylim(0, d.max() * 1.18)
    ax.grid(axis="y", color=GRID, lw=0.8)


KEY = ("In panels a-c: shaded = 90% highest-density region and dashed line = mode, each computed in "
       "that panel's own coordinates.")


ROWS = [("mode", "most probable value (mode / MAP)"),
        ("hpd", f"{int(MASS*100)}% highest-density interval"),
        ("median", "median"),
        ("et", f"{int(MASS*100)}% equal-tailed interval"),
        ("mean", "mean of the reported quantity"),
        ("coord_mean", "average the coordinate, then convert")]

# ---------------------------------------------------------------- Fig A: lambda ~ Exp(1)
lam = np.linspace(1e-9, 12, 400_001)
u = np.linspace(-14, np.log(12), 400_001)
d_lam = np.exp(-lam)
d_u = np.exp(u - np.exp(u))                      # e^{-lambda} * |dlambda/du|
sA = [summaries(lam, d_lam, lambda v: v), summaries(u, d_u, np.exp)]
for s in sA:
    s["mean"] = 1.0                               # E[lambda] -- invariant by construction
sA_native_u = summaries(u, d_u, lambda v: v)

fig = plt.figure(figsize=(11, 7.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.55, wspace=0.25)
density_panel(fig.add_subplot(gs[0, 0]), lam[lam < 5], d_lam[lam < 5], sA[0], C[0],
              "a  Posterior density in λ coordinates", "λ", "ruler: λ")
density_panel(fig.add_subplot(gs[0, 1]), u[u > -5], d_u[u > -5], sA_native_u, C[1],
              "b  The same posterior in u = log λ coordinates", "u = log λ", "ruler: log λ")
axc = fig.add_subplot(gs[1, :])
strip(axc, ROWS, [{"s": sA[0], "c": C[0]}, {"s": sA[1], "c": C[1]}], "λ (every summary converted back to λ)",
      (-0.1, 4.2), "c  Each summary computed in each ruler, then reported on the λ scale.\n"
      "Reading rule: where the two colours coincide, the summary is a property of the posterior; "
      "where they split, it is a property of the coordinates.")
axc.plot([], [], "o", color=C[0], label="computed in λ")
axc.plot([], [], "o", color=C[1], label="computed in log λ")
axc.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, frameon=False, fontsize=9)
fig.suptitle("One posterior, two rulers: λ ~ Exponential(1)", x=0.06, y=1.0, ha="left", fontsize=13,
             fontweight="bold")
fig.text(0.06, 0.955, KEY, fontsize=9, color=INK2)
fig.savefig(OUT / "fig_coord_A_exp1.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig B: a late node, three rulers
a_, b_ = 1.5, 6.0                                 # illustrative posterior on height h (fraction of 45 d)
T_EXP, A, TAU = 45.0, 4.9553, 10.7089            # Park clock constants (§H.9.5)
day_of_h = lambda h: T_EXP * (1 - h)
E = lambda day: A * (1 - np.exp(-day / TAU))
Einv = lambda e: -TAU * np.log(1 - e / A)

h = np.linspace(1e-7, 1 - 1e-7, 400_001)
d_h = stats.beta.pdf(h, a_, b_)
x = np.linspace(-14, 8, 400_001)
d_x = stats.beta.pdf(special.expit(x), a_, b_) * special.expit(x) * (1 - special.expit(x))
e = np.linspace(E(1e-6), E(T_EXP) - 1e-9, 400_001)
day_e = Einv(e)
d_e = stats.beta.pdf(1 - day_e / T_EXP, a_, b_) * (1 / T_EXP) * TAU / (A - e)   # |dh/dday|*|dday/de|

sB = [summaries(h, d_h, day_of_h),
      summaries(x, d_x, lambda v: day_of_h(special.expit(v))),
      summaries(e, d_e, Einv)]
mean_day = day_of_h(a_ / (a_ + b_))
for s in sB:
    s["mean"] = mean_day

fig = plt.figure(figsize=(12, 8))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.2], hspace=0.6, wspace=0.3)
nat = [summaries(h, d_h, lambda v: v), summaries(x, d_x, lambda v: v), summaries(e, d_e, lambda v: v)]
density_panel(fig.add_subplot(gs[0, 0]), h[h < 0.8], d_h[h < 0.8], nat[0], C[0],
              "a  Calendar height h", "h (fraction of 45 d before harvest)", "ruler: calendar time (engine point estimate)")
density_panel(fig.add_subplot(gs[0, 1]), x[(x > -6) & (x < 3)], d_x[(x > -6) & (x < 3)], nat[1], C[1],
              "b  Ratio coordinate x = logit h", "x", "ruler: ratio coordinate (where HMC runs)")
me = e > E(25)
density_panel(fig.add_subplot(gs[0, 2]), e[me], d_e[me], nat[2], C[2],
              "c  Recorder depth (Park clock)", "E(day), edits/tape", "ruler: recorder time (engine option F)")
axd = fig.add_subplot(gs[1, :])
strip(axd, ROWS, [{"s": s, "c": c} for s, c in zip(sB, C)], "dissemination-node date (day of the 45-day experiment)",
      (20, 45.5), "d  Each summary computed in each ruler, reported as a calendar day.  Reading rule: where the colours coincide, the summary is a\n"
      "property of the posterior; where they split, a property of the coordinates.  Illustrative posterior h ~ Beta(1.5, 6) — not fitted to any data.")
for c, lab in zip(C, ["computed in calendar time", "computed in ratio coordinate", "computed in recorder time"]):
    axd.plot([], [], "o", color=c, label=lab)
axd.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncol=3, frameon=False, fontsize=9)
fig.suptitle("A weakly informed late node: where is its 'most probable' date?", x=0.06, y=1.0, ha="left",
             fontsize=13, fontweight="bold")
fig.text(0.06, 0.955, KEY, fontsize=9, color=INK2)
fig.savefig(OUT / "fig_coord_B_node.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- printed record + MC check
rng = np.random.default_rng(13102026)
hs = rng.beta(a_, b_, 2_000_000)
ds = day_of_h(hs)
print("Fig A (lambda scale):")
for name, s in zip(["lambda", "log lambda"], sA):
    print(f"  {name:10s} mode {s['mode']:.3f}  HPD {np.round(s['hpd'],3)}  median {s['median']:.4f}"
          f"  ET {np.round(s['et'],3)}  coord-mean->lambda {s['coord_mean']:.4f}")
print("Fig B (calendar day):")
for name, s in zip(["calendar", "ratio x", "recorder e"], sB):
    print(f"  {name:10s} mode {s['mode']:.2f}  HPD {np.round(s['hpd'],2)}  median {s['median']:.2f}"
          f"  ET {np.round(s['et'],2)}  coord-mean->day {s['coord_mean']:.2f}")
print(f"  MC (2e6 draws): median day {np.median(ds):.2f}  ET {np.round(np.quantile(ds,[.05,.95]),2)}"
      f"  mean day {ds.mean():.2f}  Einv(mean E(day)) {Einv(E(ds).mean()):.2f}"
      f"  expit(mean logit h)->day {day_of_h(special.expit(special.logit(hs).mean())):.2f}")
print("wrote", OUT / "fig_coord_A_exp1.png", OUT / "fig_coord_B_node.png")
