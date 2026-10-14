"""Teaching figure for Betancourt §2.3: a pinched typical set (Neal's funnel) defeats a chain silently.

Target: v ~ N(0, 3^2), x_i | v ~ N(0, e^v), i = 1..9  (D = 10). Exact draws are available, so truth is
known: P(v < -3) = Phi(-1) = 0.159, E[v] = 0.
Sampler: random-walk Metropolis, isotropic step s = 0.5 (acceptance ~0.37 -- looks well tuned), 4 chains
started from EXACT draws of the target.
Contrast: the same sampler, same s, on the non-centred parameterisation (v, z), x_i = e^{v/2} z_i,
z_i ~ N(0, 1) -- the identical distribution with the funnel straightened out.
Panels: a geometry + where the centred chains went; b one chain's v trace; c running estimate of
P(v < -3), centred, 4 chains; d the same, non-centred. Split R-hat printed per parameterisation.
Notes: sciphy_notes.md §G.4.14.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
CH = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # palette slots 1-4 (adjacent-pair use)
INK, INK2, MUTED, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#a3a29c", "#e4e3df", "#fcfcfb", "#efeeea"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
})
DX, S, N, NCH, WARM, SEED = 9, 0.5, 300_000, 4, 30_000, 1
TRUE_P = stats.norm.cdf(-1.0)


def lp_centred(v, x):
    return -v**2 / 18 - 0.5 * DX * v - 0.5 * np.sum(x**2, axis=-1) * np.exp(-v)


def lp_noncentred(v, z):
    return -v**2 / 18 - 0.5 * np.sum(z**2, axis=-1)


def rwm(logp, v0, y0, seed):
    rng = np.random.default_rng(seed)
    v, y = v0.copy(), y0.copy()
    lp = logp(v, y)
    V = np.empty((N, NCH)); Y1 = np.empty((N, NCH)); acc = 0.0
    for t in range(N):
        vp = v + S * rng.standard_normal(NCH)
        yp = y + S * rng.standard_normal((NCH, DX))
        lpp = logp(vp, yp)
        a = np.log(rng.random(NCH)) < lpp - lp
        v = np.where(a, vp, v); y = np.where(a[:, None], yp, y); lp = np.where(a, lpp, lp)
        V[t] = v; Y1[t] = y[:, 0]; acc += a.mean()
    return V, Y1, acc / N


def split_rhat(chains):                                   # chains: (n, m)
    n = chains.shape[0] // 2
    seqs = np.concatenate([chains[:n], chains[n:2 * n]], axis=1)
    W = seqs.var(axis=0, ddof=1).mean()
    B = n * seqs.mean(axis=0).var(ddof=1)
    return np.sqrt(((n - 1) / n * W + B / n) / W)


rng0 = np.random.default_rng(SEED)
v0 = rng0.normal(0, 3, NCH)                               # starts are exact draws of the target
z0 = rng0.standard_normal((NCH, DX))
Vc, X1c, acc_c = rwm(lp_centred, v0, z0 * np.exp(v0 / 2)[:, None], SEED + 1)
Vn, _, acc_n = rwm(lp_noncentred, v0, z0, SEED + 1)
it = np.arange(1, N + 1)
run = lambda V: np.vstack([np.full((WARM, NCH), np.nan),
                           np.cumsum(V[WARM:] < -3, axis=0) / np.arange(1, N - WARM + 1)[:, None]])
rh_c, rh_n = split_rhat(Vc[WARM:]), split_rhat(Vn[WARM:])
rh_ci, rh_ni = split_rhat((Vc[WARM:] < -3).astype(float)), split_rhat((Vn[WARM:] < -3).astype(float))

fig = plt.figure(figsize=(15, 9.2))
gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.3, width_ratios=[1, 1.35])

# a -- geometry
ax = fig.add_subplot(gs[0, 0])
ex = np.random.default_rng(7)
ve = ex.normal(0, 3, 6000); xe = ex.standard_normal(6000) * np.exp(ve / 2)
ax.scatter(xe, ve, s=3, color=MUTED, alpha=0.5, lw=0, label="exact draws of the target")
thin = slice(WARM, None, 150)
ax.scatter(X1c[thin].ravel(), Vc[thin].ravel(), s=3, color=CH[0], alpha=0.5, lw=0,
           label="where the 4 random-walk chains went")
vmin_c = Vc[WARM:].min()
ax.axhline(vmin_c, color=CH[0], lw=1, ls="--")
ax.text(-9.6, vmin_c - 0.35, f"deepest any chain reached: v = {vmin_c:.1f}", color=INK, fontsize=8.5, va="top")
ax.text(9.6, -7.6, f"true probability below that depth: {stats.norm.cdf(vmin_c/3):.1%}", color=INK2,
        fontsize=8.5, ha="right")
ax.set_xlim(-10, 10); ax.set_ylim(-9, 9)
ax.set_xlabel("$x_1$  (one of the 9 lower-level parameters)")
ax.set_ylabel("v  (log variance of the $x_i$)")
ax.set_title("a  The funnel: wide mouth above, narrow neck below", loc="left", pad=16)
ax.text(0, 1.02, "v ~ N(0, 3²),  xᵢ | v ~ N(0, eᵛ),  D = 10", transform=ax.transAxes, fontsize=8.5,
        color=INK2, va="bottom")
ax.legend(loc="upper left", frameon=False, fontsize=8, markerscale=3)
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)

# b -- one chain's trace
axb = fig.add_subplot(gs[0, 1])
k = int(np.argmin(Vc[WARM:].min(axis=0)))                 # the chain that went deepest
axb.axhspan(-9, -3, color=BAND, zorder=0)
axb.plot(it[::20], Vc[::20, k], color=CH[k], lw=0.8)
axb.text(N * 0.995, -8.6, "shaded: the neck (v < −3), 15.9% of the true probability", ha="right", va="bottom",
         fontsize=8.5, color=INK2)
axb.set_ylim(-9, 9); axb.set_xlim(0, N)
axb.set_xlabel("iteration"); axb.set_ylabel("v")
axb.set_title(f"b  One chain's trace (chain {k+1}): long sticky visits to the neck's edge", loc="left", pad=16)
axb.text(0, 1.02, f"random-walk Metropolis, step 0.5, acceptance {acc_c:.0%} — tuning looks healthy",
         transform=axb.transAxes, fontsize=8.5, color=INK2, va="bottom")
axb.grid(color=GRID, lw=0.8); axb.set_axisbelow(True)


def running_panel(a, R, title, sub, rh, labels=True):
    ends = R[-1].copy()
    order = np.argsort(ends)
    for a_, b_ in zip(order[:-1], order[1:]):                # keep end labels >= 0.018 apart
        if ends[b_] - ends[a_] < 0.018:
            ends[b_] = ends[a_] + 0.018
    for j in range(NCH):
        a.plot(it[::50], R[::50, j], color=CH[j], lw=1.6)
        if labels:
            a.text(N * 1.005, ends[j], f"chain {j+1}", color=INK, fontsize=8, va="center")
    a.axhline(TRUE_P, color=INK, lw=1.2, ls="--")
    a.text(N * 0.12, TRUE_P + 0.006, f"truth {TRUE_P:.3f}", fontsize=8.5, color=INK, va="bottom")
    a.axvspan(0, WARM, color=BAND, zorder=0)
    a.text(WARM / 2, 0.375, "warm-up\n(discarded)", ha="center", va="top", fontsize=8, color=INK2)
    a.set_xlim(0, N); a.set_ylim(0, 0.4)
    a.set_xlabel("iteration")
    a.set_ylabel("running estimate of P(v < −3),\nafter warm-up")
    a.set_title(title, loc="left", pad=16)
    a.text(0, 1.02, sub + f";  split R̂: v {rh[0]:.2f},  1[v < −3] {rh[1]:.2f}",
           transform=a.transAxes, fontsize=8.5, color=INK2, va="bottom")
    a.grid(color=GRID, lw=0.8); a.set_axisbelow(True)


running_panel(fig.add_subplot(gs[1, 0]), run(Vc), "c  Funnel as written: estimates oscillate and disagree",
              "4 chains started from exact draws", (rh_c, rh_ci))
running_panel(fig.add_subplot(gs[1, 1]), run(Vn),
              "d  Same distribution, non-centred (xᵢ = e^(v/2) zᵢ): same sampler, same step",
              f"acceptance {acc_n:.0%}; all 4 chains overlap", (rh_n, rh_ni), labels=False)

fig.suptitle("A pinched typical set: the chain cannot enter the neck, so its estimates are biased and unstable",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.957, "Neal's funnel (exact answer known).  Reading rule: in c and d, a healthy run converges onto "
         "the dashed truth with all chains together; R̂ near 1.00 means the chains agree.", fontsize=9, color=INK2)
fig.savefig(OUT / "fig_funnel_pathology.png", dpi=200, bbox_inches="tight")
plt.close(fig)

post = slice(WARM, None)
print(f"centred:     acc {acc_c:.2f}  P(v<-3) per chain {np.round((Vc[post] < -3).mean(0), 3)}  "
      f"mean v {np.round(Vc[post].mean(0), 2)}  min v {Vc[post].min():.2f}  split-Rhat v {rh_c:.3f} ind {rh_ci:.3f}")
print(f"non-centred: acc {acc_n:.2f}  P(v<-3) per chain {np.round((Vn[post] < -3).mean(0), 3)}  "
      f"mean v {np.round(Vn[post].mean(0), 2)}  min v {Vn[post].min():.2f}  split-Rhat v {rh_n:.3f} ind {rh_ni:.3f}")
print(f"truth: P(v<-3) {TRUE_P:.4f}, E[v] 0;  P(v < deepest centred v) {stats.norm.cdf(vmin_c/3):.4f}")
print("starts v0:", np.round(v0, 2))
print("wrote", OUT / "fig_funnel_pathology.png")
