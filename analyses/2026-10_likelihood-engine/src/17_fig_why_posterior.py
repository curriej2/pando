"""Teaching figure (notes §G.4.20): why kept HMC positions are distributed as the posterior even though the
particle moves FASTEST through the mode. 1-D standard Gaussian, K = p^2/2, exact flow (orbits are circles
of radius R = sqrt(q^2 + p^2) in the (q, p) plane; q swings between -R and +R).

a  One orbit: share of time spent near each q (arcsine law, U-shaped -- slow at the edges, fast through the
   mode), for three orbit sizes.
b  How big the orbits are: R from a stationary start q ~ N(0,1) and a fresh p ~ N(0,1) (Rayleigh) --
   mostly small.
c  Stack all orbits, weighted by how often each size is drawn: the U-shaped pieces sum to the bell curve.
d  The same fact as head-count x dwell time: the fraction of orbits that reach q is exp(-q^2/2), and among
   those the average time spent at q is the same for every q.
e  What HMC keeps: one endpoint per transition (random flow time), histogram vs target.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import integrate, stats
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#a3a29c", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
rng = np.random.default_rng(19102026)
N = 2_000_000
q0 = rng.standard_normal(N); p0 = rng.standard_normal(N)
R = np.hypot(q0, p0)                                     # orbit radius set by (current q, drawn p)
theta = rng.uniform(0, 2 * np.pi, N)                     # a uniformly random moment along the orbit
q_on_orbit = R * np.cos(theta)
xg = np.linspace(-3.5, 3.5, 701)

fig = plt.figure(figsize=(16, 10.2))
gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 0.9, 0.9], hspace=1.1, wspace=0.95)


def subtitle(a, txt):
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.5, color=INK2, va="bottom")


# ---------------------------------------------------------------- a: one orbit's dwell time
ax = fig.add_subplot(gs[0:2, 0:2])
for r, col in zip([0.6, 1.2, 2.0], C):
    x = np.linspace(-r * 0.999, r * 0.999, 600)
    ax.plot(x, 1 / (np.pi * np.sqrt(r * r - x * x)), color=col, label=f"orbit size R = {r}")
ax.annotate("fast through the mode:\nlittle time here", (0, 1 / (np.pi * 2.0)), (-1.9, 0.95), fontsize=8.5,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.annotate("slow at the turning\npoints: lots of time", (1.95, 0.9), (0.55, 1.05), fontsize=8.5,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.set_xlim(-2.4, 2.4); ax.set_ylim(0, 1.7)
ax.legend(frameon=False, fontsize=8, loc="upper center")
ax.set_xlabel("q"); ax.set_ylabel("share of one orbit's time, per unit q")
ax.set_title("a  ONE orbit lingers at its edges", loc="left", pad=28)
subtitle(ax, "an orbit of size R swings between −R and +R;\nU-shaped curve = where it lingers")
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------------------------------------------------------------- b: orbit sizes
axb = fig.add_subplot(gs[0:2, 2:4])
rr = np.linspace(0, 4.5, 400)
axb.hist(R, bins=90, range=(0, 4.5), density=True, color=SEQ[2], alpha=0.8, label="2,000,000 momentum draws")
axb.plot(rr, rr * np.exp(-rr**2 / 2), color=INK, lw=1.8, label="exact (Rayleigh)")
axb.axvline(np.median(R), color=INK2, ls="--", lw=1)
axb.text(np.median(R) + 0.08, 0.50, f"median orbit\nsize {np.median(R):.2f}", fontsize=8.5, color=INK)
axb.set_xlabel("orbit size R  (how far q can swing from the mode)")
axb.set_ylabel("how often an orbit this size is drawn")
axb.set_title("b  …but most orbits are SMALL", loc="left", pad=28)
subtitle(axb, f"budget = cost of current q + energy of drawn p;\n{np.mean(R < 1):.0%} of orbits never swing past q = ±1")
axb.legend(frameon=False, fontsize=8, loc="upper right")
axb.grid(axis="y", color=GRID, lw=0.8); axb.set_axisbelow(True)

# ---------------------------------------------------------------- c: stacked contributions
axc = fig.add_subplot(gs[0:2, 4:6])
bands = [0, 0.5, 1.0, 1.5, 2.0, 2.5, np.inf]
edges = np.linspace(-3.5, 3.5, 71)
w = np.diff(edges)[0]
bottom = np.zeros(len(edges) - 1)
for i in range(len(bands) - 1):
    m = (R >= bands[i]) & (R < bands[i + 1])
    h, _ = np.histogram(q_on_orbit[m], bins=edges)
    h = h / (N * w)
    lab = f"R {bands[i]:.1f}–{bands[i+1]:.1f}" if np.isfinite(bands[i + 1]) else f"R ≥ {bands[i]:.1f}"
    axc.bar(edges[:-1], h, width=w, bottom=bottom, align="edge", color=SEQ[i], edgecolor=SURF, lw=0.3, label=lab)
    bottom += h
axc.plot(xg, stats.norm.pdf(xg), color=INK, lw=2, label="target posterior")
axc.set_xlim(-3.5, 3.5)
axc.set_xlabel("q"); axc.set_ylabel("density")
axc.set_title("c  Stack them: the U's add up to the bell", loc="left", pad=28)
subtitle(axc, "each colour = orbits of one size, at a random moment;\nsmall orbits (light) pile their edges near the mode")
axc.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=1)
axc.grid(axis="y", color=GRID, lw=0.8); axc.set_axisbelow(True)

# ---------------------------------------------------------------- d: head-count x dwell
xs = np.linspace(0, 3.2, 33)
reach_mc = np.array([np.mean(R >= x) for x in xs])
dwell = []
for x in xs:
    val, _ = integrate.quad(lambda r: r * np.exp(-r * r / 2) / (np.pi * np.sqrt(r * r - x * x)), x, np.inf, limit=200)
    dwell.append(val / np.exp(-x * x / 2))
dwell = np.array(dwell)
axd1 = fig.add_subplot(gs[2, 0:3])
axd1.plot(xs, np.exp(-xs**2 / 2), color=C[0])
axd1.plot(xs, reach_mc, "o", color=INK, ms=3.5)
axd1.set_ylabel("fraction of orbits\nthat reach q")
axd1.set_title("d  Head-count × dwell time  =  posterior", loc="left", pad=28)
subtitle(axd1, "top: every orbit passes the mode; only big ones reach the tails (dots = Monte Carlo)\n"
         "bottom: among orbits that reach q, the average time spent there is the SAME at every q")
axd1.set_xlim(0, 3.2); axd1.grid(color=GRID, lw=0.8); axd1.set_axisbelow(True)
axd2 = fig.add_subplot(gs[3, 0:3], sharex=axd1)
axd2.plot(xs, dwell, color=C[1])
axd2.set_ylim(0, 0.6)
axd2.set_ylabel("average time at q,\namong those orbits")
axd2.set_xlabel("distance of q from the mode")
axd2.text(1.6, 0.47, f"flat at {dwell.mean():.3f} = 1/√(2π)", fontsize=8.5, color=INK)
axd2.grid(color=GRID, lw=0.8); axd2.set_axisbelow(True)

# ---------------------------------------------------------------- e: what HMC keeps
axe = fig.add_subplot(gs[2:4, 3:6])
n_tr = 50_000
q = 2.5; kept = np.empty(n_tr)
for i in range(n_tr):
    p = rng.standard_normal(); T = rng.uniform(0.5, 2.5)
    q = q * np.cos(T) + p * np.sin(T)                    # exact flow, endpoint only
    kept[i] = q
axe.hist(kept, bins=70, range=(-3.5, 3.5), density=True, color=C[0], alpha=0.55,
         label=f"{n_tr:,} kept endpoints (one per momentum draw)")
axe.plot(xg, stats.norm.pdf(xg), color=INK, lw=2, label="target posterior")
axe.set_xlabel("q"); axe.set_ylabel("density")
axe.set_title("e  What HMC keeps: ONE endpoint per trajectory", loc="left", pad=28)
subtitle(axe, f"flow time drawn between 0.5 and 2.5;  mean {kept.mean():+.3f} (truth 0),  SD {kept.std():.3f} (truth 1)")
axe.legend(frameon=False, fontsize=8, loc="upper right")
axe.grid(axis="y", color=GRID, lw=0.8); axe.set_axisbelow(True)

fig.suptitle("Why the kept positions look like the posterior, even though the particle is fastest at the mode",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.965, "1-D standard Gaussian, exact orbits.  Read a → b → c (one orbit, how many of each size, "
         "all of them stacked), then d (the same fact as head-count × dwell time), then e (the sampler's output).",
         fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_why_posterior.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"b: P(R<1) {np.mean(R<1):.4f} (exact {1-np.exp(-0.5):.4f}); median R {np.median(R):.3f}")
print(f"c: stacked total vs N(0,1): max abs diff {np.max(np.abs(bottom - stats.norm.pdf(edges[:-1]+w/2))):.4f}")
print(f"d: dwell among reaching orbits: min {dwell.min():.5f} max {dwell.max():.5f} (1/sqrt(2pi) {1/np.sqrt(2*np.pi):.5f})")
print(f"e: kept mean {kept.mean():.4f} SD {kept.std():.4f}")
print("wrote", OUT / "fig_hmc_why_posterior.png")
