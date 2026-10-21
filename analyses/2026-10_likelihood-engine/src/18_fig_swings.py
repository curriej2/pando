"""Teaching figure (notes §G.4.21): the swing picture of HMC, and what it buys over a random walk.
1-D standard Gaussian (V = q^2/2, K = p^2/2, exact flow) for a-c; a 100-D standard Gaussian for d.

a  Swing size = pull-back + push: amplitude R = sqrt(q0^2 + p0^2) from the starting position q0 and the
   drawn momentum p0. Same push from further out -> bigger swing; big push from the middle -> big swing.
b  The push does not depend on position: momentum draws when the chain sits near the mode vs far out are
   the same N(0,1) (b1); the resulting swing sizes differ only because of the pull-back (b2).
c  No starting points are supplied: a chain started at q0 = 4 falls into the posterior on its own.
d  The advantage over a random walk: one coordinate of a 100-D Gaussian traced against computational
   cost (evaluations), HMC vs random-walk Metropolis, both started at stationarity.
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
rng = np.random.default_rng(20102026)


def subtitle(a, txt):
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.5, color=INK2, va="bottom")


fig = plt.figure(figsize=(16, 10.4))
gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 1, 1], hspace=1.25, wspace=1.0)

# ---------------------------------------------------------------- a: pull-back + push
ax = fig.add_subplot(gs[0:2, 0:3])
t = np.linspace(0, 2 * np.pi, 400)
cases = [(0.2, 0.5, "small pull-back, small push"), (2.0, 0.5, "big pull-back, same small push"),
         (0.2, 2.0, "small pull-back, big push")]
for (q0, p0, lab), col in zip(cases, C):
    q = q0 * np.cos(t) + p0 * np.sin(t)
    Rr = np.hypot(q0, p0)
    ax.plot(t, q, color=col, label=f"{lab}:  q₀ = {q0}, p₀ = {p0}  →  swing size R = {Rr:.2f}")
    ax.plot(0, q0, "o", color=col, ms=7, mec=SURF, mew=1.5, zorder=5)
ax.axhline(0, color=INK2, lw=0.8, ls=":")
ax.text(6.35, 0.05, "mode", fontsize=8, color=INK2, va="bottom", ha="right")
ax.set_xlim(0, 2 * np.pi); ax.set_ylim(-2.6, 3.2)
ax.set_xlabel("fictitious time t  (one full swing shown)")
ax.set_ylabel("q")
ax.set_title("a  Swing size = how far you start out  +  how hard you are pushed", loc="left", pad=28)
subtitle(ax, "● = where the swing starts (the current q);  budget H = q₀²/2 + p₀²/2, so swing size R = √(2H) = √(q₀² + p₀²)")
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------------------------------------------------------------- b: push independent of position
n = 400_000
qs = rng.standard_normal(n)                              # where the chain sits (stationary)
ps = rng.standard_normal(n)                              # the push -- drawn WITHOUT looking at q
Rs = np.hypot(qs, ps)
near, far = np.abs(qs) < 0.5, np.abs(qs) > 1.5
axb1 = fig.add_subplot(gs[0, 3:6])
bins = np.linspace(-4, 4, 61)
axb1.hist(ps[near], bins=bins, density=True, histtype="step", color=C[0], lw=2, label="chain near the mode (|q| < 0.5)")
axb1.hist(ps[far], bins=bins, density=True, histtype="step", color=C[1], lw=2, label="chain far out (|q| > 1.5)")
axb1.plot(bins, stats.norm.pdf(bins), color=INK, lw=1, ls="--", label="N(0, 1)")
axb1.set_xlabel("momentum drawn (the push)"); axb1.set_ylabel("density")
axb1.set_title("b  The push does NOT depend on where you are…", loc="left", pad=28)
subtitle(axb1, f"SD of pushes: near the mode {ps[near].std():.3f}, far out {ps[far].std():.3f}")
axb1.legend(frameon=False, fontsize=7.5, loc="upper right")
axb1.grid(axis="y", color=GRID, lw=0.8); axb1.set_axisbelow(True)
axb2 = fig.add_subplot(gs[1, 3:6])
rb = np.linspace(0, 4.5, 46)
axb2.hist(Rs[near], bins=rb, density=True, histtype="step", color=C[0], lw=2, label="chain near the mode")
axb2.hist(Rs[far], bins=rb, density=True, histtype="step", color=C[1], lw=2, label="chain far out")
axb2.set_xlabel("swing size R"); axb2.set_ylabel("density")
axb2.set_title("…but the swing does, through the pull-back", loc="left", pad=28)
subtitle(axb2, f"median swing: near the mode {np.median(Rs[near]):.2f}, far out {np.median(Rs[far]):.2f}  "
         "(a swing always reaches at least as far as its start)")
axb2.legend(frameon=False, fontsize=7.5, loc="upper right")
axb2.grid(axis="y", color=GRID, lw=0.8); axb2.set_axisbelow(True)

# ---------------------------------------------------------------- c: chain started far out falls in
axc = fig.add_subplot(gs[2:4, 0:3])
for rep, col in zip(range(3), C):
    q = 4.0; path_t = [0.0]; path_q = [q]; kept_t = [0.0]; kept_q = [q]; t0 = 0.0
    for i in range(12):
        p = rng.standard_normal(); T = rng.uniform(0.5, 2.5)
        tt = np.linspace(0, T, 60)
        seg = q * np.cos(tt) + p * np.sin(tt)
        path_t += list(t0 + tt[1:]); path_q += list(seg[1:])
        q = seg[-1]; t0 += T
        kept_t.append(t0); kept_q.append(q)
    axc.plot(path_t, path_q, color=col, lw=1.2, alpha=0.8)
    axc.plot(kept_t, kept_q, "s", color=col, ms=4.5, zorder=5)
axc.axhspan(-1.96, 1.96, color=BAND, zorder=0)
axc.text(0.3, -2.25, "shaded: where 95% of the posterior lies", fontsize=8.5, color=INK2)
axc.plot(0, 4, "o", color=INK, ms=7, zorder=6)
axc.annotate("all three chains start here,\nfar outside the posterior", (0, 4), (2.2, 3.9), fontsize=8.5,
             arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axc.set_ylim(-4.6, 4.6)
axc.set_xlabel("fictitious time, transitions laid end to end")
axc.set_ylabel("q")
axc.set_title("c  Nobody supplies good starting points: the chain falls in by itself", loc="left", pad=34)
subtitle(axc, "3 chains from q₀ = 4; ■ = kept samples. The first swings are huge (big pull-back) and carry the\n"
         "chain inside within one or two transitions; from then on the swings are posterior-sized")
axc.grid(color=GRID, lw=0.8); axc.set_axisbelow(True)

# ---------------------------------------------------------------- d: HMC vs RWM at equal cost, D = 100
D, BUDGET = 100, 3000
x = rng.standard_normal(D)                               # stationary start
lp = -0.5 * x @ x; s = 2.38 / np.sqrt(D); rwm_trace = []
for i in range(BUDGET):                                  # one likelihood evaluation per iteration
    xp = x + s * rng.standard_normal(D); lpp = -0.5 * xp @ xp
    if np.log(rng.random()) < lpp - lp:
        x, lp = xp, lpp
    rwm_trace.append(x[0])
x = rng.standard_normal(D); eps = D ** -0.25; used = 0; hmc_cost = [0]; hmc_trace = [x[0]]
while used < BUDGET:                                     # one gradient evaluation per leapfrog step
    p = rng.standard_normal(D); L = int(np.ceil(rng.uniform(0.8, 1.6) / eps))
    x0, H0 = x.copy(), 0.5 * x @ x + 0.5 * p @ p
    p = p - 0.5 * eps * x
    for l in range(L):
        x = x + eps * p
        if l < L - 1:
            p = p - eps * x
    p = p - 0.5 * eps * x
    if np.log(rng.random()) >= H0 - (0.5 * x @ x + 0.5 * p @ p):
        x = x0
    used += L; hmc_cost.append(used); hmc_trace.append(x[0])
ht = np.array(hmc_trace)
r_rwm = np.corrcoef(rwm_trace[:-1], rwm_trace[1:])[0, 1]; r_hmc = np.corrcoef(ht[:-1], ht[1:])[0, 1]
axd1 = fig.add_subplot(gs[2, 3:6])
axd1.axhspan(-1.96, 1.96, color=BAND, zorder=0)
axd1.plot(np.arange(1, BUDGET + 1), rwm_trace, color=C[1], lw=1.4)
axd1.set_xlim(0, BUDGET); axd1.set_ylim(-5, 5)
axd1.set_ylabel("q₁")
axd1.set_title("d  The payoff, at equal cost: random walk (top) vs HMC (bottom)", loc="left", pad=28)
subtitle(axd1, f"random-walk Metropolis, 1 evaluation per step: consecutive samples correlate {r_rwm:.3f} — it creeps")
axd1.grid(color=GRID, lw=0.8); axd1.set_axisbelow(True)
axd2 = fig.add_subplot(gs[3, 3:6], sharex=axd1)
axd2.axhspan(-1.96, 1.96, color=BAND, zorder=0)
axd2.plot(hmc_cost, ht, "-o", color=C[0], lw=0.8, ms=2.5)
axd2.set_ylim(-5, 5)
axd2.set_xlabel("computational cost  (evaluations of the posterior or its gradient)")
axd2.set_ylabel("q₁")
subtitle(axd2, f"HMC, ≈{np.mean(np.diff(hmc_cost)):.0f} gradients per sample: consecutive samples correlate {r_hmc:.3f} "
         "— it jumps across")
axd2.grid(color=GRID, lw=0.8); axd2.set_axisbelow(True)
axd2.text(BUDGET * 0.99, -4.6, "q₁ = one coordinate of a 100-D posterior; shaded = 95% of it", ha="right",
          fontsize=8, color=INK2)

fig.suptitle("The swing picture: where the swing starts, how hard it is pushed, and why that beats a random walk",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.965, "a–c: 1-D standard Gaussian, exact swings.  d: 100-D standard Gaussian, leapfrog steps.",
         fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_swings.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"b: push SD near {ps[near].std():.4f} far {ps[far].std():.4f}; median R near {np.median(Rs[near]):.3f} "
      f"far {np.median(Rs[far]):.3f}; min R far {Rs[far].min():.3f}")
print(f"d: RWM range of q1 over {BUDGET} evals: [{min(rwm_trace):.2f}, {max(rwm_trace):.2f}], "
      f"HMC samples {len(hmc_trace)-1}, range [{min(hmc_trace):.2f}, {max(hmc_trace):.2f}]")
print(f"d: HMC trace index of max {ht.argmax()} value {ht.max():.2f}; first value {ht[0]:.2f}; SD {ht.std():.3f}; RWM SD {np.std(rwm_trace):.3f}, lag-1 corr RWM {np.corrcoef(rwm_trace[:-1], rwm_trace[1:])[0,1]:.3f} HMC {np.corrcoef(ht[:-1], ht[1:])[0,1]:.3f}")
print("wrote", OUT / "fig_hmc_swings.png")
