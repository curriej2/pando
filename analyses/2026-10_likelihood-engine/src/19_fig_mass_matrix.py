"""Teaching figure for Betancourt §4.2.1 (notes §G.4.27): what the mass matrix M does and how it is obtained.

Row 1, an elongated valley (posterior SDs 10 and 1): one HMC trajectory with the SAME number of leapfrog
  steps (a) with no mass matrix (M = I) and (b) with a diagonal M (M^-1 = posterior variances); (c) the
  trajectory of (b) redrawn in rescaled coordinates q / SD, where the valley is round.
Row 2, a tilted ridge (unit SDs, correlation 0.98): (d) diagonal M (cannot rotate), (e) dense M
  (M^-1 = posterior covariance), (f) the trajectory of (e) in whitened coordinates, where the ridge is round.
Row 3, how M is obtained: warm-up on an 8-D posterior whose coordinates have widths 0.25 ... 6 (think: early,
  well-recorded node dates ... late, weakly recorded ones). Start with M = I; after each warm-up window,
  set M^-1 = (regularised) variance of that window's draws; freeze M for sampling.
  (g) trace of the widest coordinate through the windows, (h) what the adapted M reads out: estimated
  posterior SD per coordinate vs the truth.
Step size in every case = 0.9 / (fastest swing frequency under that M), so each run is as fast as stability
allows; leapfrog with accept/reject throughout.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.3, color=INK2, va="bottom")


def eps_for(Prec, Minv):
    return 0.9 / np.sqrt(np.max(np.linalg.eigvals(Minv @ Prec).real))


def trajectory(Prec, Minv, q0, p0, L):
    eps = eps_for(Prec, Minv)
    q, p = q0.copy(), p0.copy(); path = [q.copy()]
    p = p - 0.5 * eps * Prec @ q
    for l in range(L):
        q = q + eps * Minv @ p
        if l < L - 1:
            p = p - eps * Prec @ q
        path.append(q.copy())
    return np.array(path), eps


def ellipse(Sig, k=2.0, n=300):
    t = np.linspace(0, 2 * np.pi, n)
    Lc = np.linalg.cholesky(Sig)
    return (Lc @ np.vstack([np.cos(t), np.sin(t)]) * k).T


def valley(ax, Sig, lim):
    for k, a in [(1, 0.9), (2, 0.6), (3, 0.35)]:
        e = ellipse(Sig, k)
        ax.plot(e[:, 0], e[:, 1], color=MUTED, lw=1.0, alpha=a + 0.1)
    ax.set_xlim(-lim[0], lim[0]); ax.set_ylim(-lim[1], lim[1])
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(True)


def draw_path(ax, path, col):
    ax.plot(path[:, 0], path[:, 1], "-o", color=col, ms=3, lw=1.6)
    ax.plot(*path[0], "o", color=INK, ms=7, zorder=5)
    ax.plot(*path[-1], "s", color=INK, ms=7, zorder=5)


L_STEPS = 20
fig = plt.figure(figsize=(16, 15))
gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.05], hspace=0.62, wspace=1.0)

# ---------------------------------------------------------------- row 1: elongated valley
S1 = np.diag([10.0**2, 1.0**2]); P1 = np.linalg.inv(S1)
q0 = np.array([-12.0, 0.6])
z = rng.standard_normal(2)
axa = fig.add_subplot(gs[0, 0:2])
valley(axa, S1, (32, 3.6))
pa, ea = trajectory(P1, np.eye(2), q0, z, L_STEPS)               # p ~ N(0, I)
draw_path(axa, pa, C[1])
axa.set_title("a  No mass matrix (M = I)", loc="left", pad=28)
subtitle(axa, f"{L_STEPS} steps of size {ea:.2f} (set by the NARROW direction):\nzig-zags across, barely moves along")
axa.set_xlabel("q₁  (wide direction, SD 10)"); axa.set_ylabel("q₂  (narrow, SD 1)")

axb = fig.add_subplot(gs[0, 2:4])
valley(axb, S1, (32, 3.6))
Mb_inv = S1.copy()                                                  # M^-1 = posterior variances
pb, eb = trajectory(P1, Mb_inv, q0, np.linalg.cholesky(np.linalg.inv(Mb_inv)) @ z, L_STEPS)  # p ~ N(0, M)
draw_path(axb, pb, C[0])
axb.set_title("b  Diagonal M:  M⁻¹ = posterior variances", loc="left", pad=28)
subtitle(axb, f"same {L_STEPS} steps, size {eb:.2f}: the wide coordinate is 'light' and moves fast,\n"
         "the narrow one 'heavy' and slow — it sweeps across the valley")
axb.set_xlabel("q₁  (wide)"); axb.set_ylabel("q₂  (narrow)")

axc = fig.add_subplot(gs[0, 4:6])
sd1 = np.sqrt(np.diag(S1))
valley(axc, np.eye(2), (3.6, 3.6))
draw_path(axc, pb / sd1, C[0])
axc.set_aspect("equal")
axc.set_title("c  …which is the same as rescaling", loc="left", pad=28)
subtitle(axc, "the run in (b) redrawn in units of each coordinate's SD:\nthe valley is round, the swing is a plain arc")
axc.set_xlabel("q₁ / 10"); axc.set_ylabel("q₂ / 1")

# ---------------------------------------------------------------- row 2: tilted ridge
r = 0.98
S2 = np.array([[1.0, r], [r, 1.0]]); P2 = np.linalg.inv(S2)
q0r = np.array([-1.8, -1.6])
axd = fig.add_subplot(gs[1, 0:2])
valley(axd, S2, (3.4, 3.4)); axd.set_aspect("equal")
Md_inv = np.diag(np.diag(S2))                                       # diagonal = variances = 1 each
pd_, ed = trajectory(P2, Md_inv, q0r, z, L_STEPS)
draw_path(axd, pd_, C[1])
axd.set_title("d  Tilted ridge, diagonal M", loc="left", pad=28)
subtitle(axd, f"both coordinates have SD 1, so diagonal M = identity:\nstep {ed:.2f} set by the thin direction; zig-zags, stays put")
axd.set_xlabel("q₁"); axd.set_ylabel("q₂")

axe = fig.add_subplot(gs[1, 2:4])
valley(axe, S2, (3.4, 3.4)); axe.set_aspect("equal")
pe, ee = trajectory(P2, S2, q0r, np.linalg.cholesky(P2) @ z, L_STEPS)
draw_path(axe, pe, C[0])
axe.set_title("e  Tilted ridge, dense M:  M⁻¹ = covariance", loc="left", pad=28)
subtitle(axe, f"same {L_STEPS} steps, size {ee:.2f}: M also ROTATES, so the motion\nlines up with the ridge and travels along it")
axe.set_xlabel("q₁"); axe.set_ylabel("q₂")

axf = fig.add_subplot(gs[1, 4:6])
Lc2 = np.linalg.cholesky(S2)
valley(axf, np.eye(2), (3.6, 3.6)); axf.set_aspect("equal")
draw_path(axf, np.linalg.solve(Lc2, pe.T).T, C[0])
axf.set_title("f  …which is rescaling + rotating", loc="left", pad=28)
subtitle(axf, "the run in (e) in 'whitened' coordinates (ridge straightened\nand squashed to unit width): round again")
axf.set_xlabel("whitened coordinate 1"); axf.set_ylabel("whitened coordinate 2")

# ---------------------------------------------------------------- row 3: how M is obtained (warm-up)
W = np.array([0.25, 0.4, 0.6, 1.0, 1.5, 2.5, 4.0, 6.0]); D = len(W)
P3 = np.diag(1 / W**2)
windows = [25, 50, 100, 200, 400]                                   # warm-up windows (iterations)
N_SAMPLE = 800
Minv = np.eye(D)
q = rng.standard_normal(D) * W
trace, bounds, eps_hist = [], [0], []
L3_RANGE = (4, 16)                                                  # jittered trajectory length (avoids returning to the start)


def hmc_step(q, Minv):
    eps = eps_for(P3, Minv)
    Mch = np.linalg.cholesky(np.linalg.inv(Minv))
    p = Mch @ rng.standard_normal(D)
    H0 = 0.5 * q @ P3 @ q + 0.5 * p @ Minv @ p; qo = q.copy()
    p = p - 0.5 * eps * P3 @ q
    L3 = rng.integers(*L3_RANGE)
    for l in range(L3):
        q = q + eps * Minv @ p
        if l < L3 - 1:
            p = p - eps * P3 @ q
    p = p - 0.5 * eps * P3 @ q
    if np.log(rng.random()) >= H0 - (0.5 * q @ P3 @ q + 0.5 * p @ Minv @ p):
        q = qo
    return q, eps


for wlen in windows:
    draws = []
    for i in range(wlen):
        q, e = hmc_step(q, Minv); draws.append(q.copy()); trace.append(q.copy())
    eps_hist.append(e)
    draws = np.array(draws); n = len(draws)
    var = draws.var(axis=0, ddof=1)
    Minv = np.diag((n / (n + 5)) * var + 1e-3 * (5 / (n + 5)))     # Stan-style regularised update
    bounds.append(bounds[-1] + wlen)
samples = []
for i in range(N_SAMPLE):
    q, e = hmc_step(q, Minv); trace.append(q.copy()); samples.append(q.copy())
eps_hist.append(e)
trace = np.array(trace); samples = np.array(samples)

axg = fig.add_subplot(gs[2, 0:4])
it = np.arange(1, len(trace) + 1)
axg.axhspan(-2 * W[-1], 2 * W[-1], color=BAND, zorder=0)
axg.plot(it, trace[:, -1], color=C[0], lw=1.2)
for k, b in enumerate(bounds[1:]):
    axg.axvline(b, color=INK2, lw=0.9, ls="--")
    lab = "M updated" if k < len(windows) - 1 else "M frozen"
    axg.text(b * 1.04, 17.2, lab, fontsize=8, color=INK2, rotation=90, va="top")
for k in range(len(windows)):
    axg.text(np.sqrt(max(bounds[k], 1) * bounds[k + 1]), -20.5, f"w{k+1}", ha="center", fontsize=8, color=INK2)
axg.text(np.sqrt(bounds[-1] * len(trace)), -20.5, "sampling (M fixed)", ha="center", fontsize=8.5, color=INK)
axg.set_ylim(-22, 18); axg.set_xscale("log"); axg.set_xlim(1, len(trace))
axg.set_xlabel("iteration (log scale; 4–15 leapfrog steps each)    w1–w5 = warm-up windows")
axg.set_ylabel("widest coordinate (true SD 6)")
axg.set_title("g  How M is learned: warm-up in windows", loc="left", pad=28)
subtitle(axg, "8-D posterior, widths 0.25 … 6. Start with M = I: the wide coordinate creeps. After each window,\n"
         f"M⁻¹ ← variance of that window's draws. Step size: {eps_hist[0]:.2f} at start → {eps_hist[-1]:.2f} once M is learned.  "
         "Shaded = ±2 true SD.")
axg.grid(color=GRID, lw=0.8); axg.set_axisbelow(True)

axh = fig.add_subplot(gs[2, 4:6])
est_sd = np.sqrt(np.diag(Minv))
xs = np.arange(D)
axh.bar(xs - 0.18, W, width=0.36, color=MUTED, label="true posterior SD")
axh.bar(xs + 0.18, est_sd, width=0.36, color=C[0], label="read off the adapted M  (√ of M⁻¹ diagonal)")
axh.set_yscale("log")
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
axh.yaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4, 8])); axh.yaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2", "4", "8"])); axh.yaxis.set_minor_locator(NullLocator())
axh.set_xticks(xs)
axh.set_xticklabels([f"{i+1}" for i in range(D)])
axh.set_xlabel("coordinate  (1 = 'early, well-recorded node' … 8 = 'late, weak')")
axh.set_ylabel("posterior SD (log scale)")
axh.set_title("h  What M tells you", loc="left", pad=28)
subtitle(axh, "after warm-up, M⁻¹ holds each coordinate's estimated\nposterior variance — which unknowns the data pin down")
axh.legend(frameon=False, fontsize=7.8, loc="upper left")
axh.grid(axis="y", color=GRID, lw=0.6); axh.set_axisbelow(True)

fig.suptitle("The mass matrix M: tune every coordinate's swing to the same rhythm — what it does, and how it is learned",
             x=0.06, y=0.995, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.968, "Grey ellipses = posterior contours.  ● = start of the trajectory, ■ = end.  Every trajectory in a row "
         "uses the same number of steps (same cost); the step size is as large as stability allows under that M.",
         fontsize=9, color=INK2)
fig.subplots_adjust(top=0.92)
fig.savefig(OUT / "fig_hmc_mass_matrix.png", dpi=200, bbox_inches="tight")
plt.close(fig)

span = lambda pth, v: np.ptp(pth @ v)
print(f"row1: eps identity {ea:.3f}, diag {eb:.3f}; distance along wide axis: identity {np.ptp(pa[:,0]):.1f}, diag {np.ptp(pb[:,0]):.1f}")
u = np.array([1, 1]) / np.sqrt(2)
print(f"row2: eps diag {ed:.3f}, dense {ee:.3f}; travel along ridge: diag {span(pd_, u):.2f}, dense {span(pe, u):.2f} (ridge SD {np.sqrt(1+r):.2f})")
print(f"row3: eps per window {np.round(eps_hist, 3)}; est SD {np.round(est_sd, 2)} vs true {W}")
print(f"sampling-phase SD of widest coord {samples[:, -1].std():.2f} (true 6)")
print("wrote", OUT / "fig_hmc_mass_matrix.png")
