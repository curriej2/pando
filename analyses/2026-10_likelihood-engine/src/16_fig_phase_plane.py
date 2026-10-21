"""Teaching figure for Betancourt §3.2-3.3 (notes §G.4.19): Hamilton's equations as "walk along the contours
of H in phase space", seen in the (q, p) plane and projected back onto q.

a  Gaussian posterior, V = q^2/2, K = p^2/2: contours of H (circles) and the rotated-gradient velocity field
   (dq/dt, dp/dt) = (dH/dp, -dH/dq) = (p, -q).
b  One trajectory taken in small steps: kick p <- p - dt dV/dq (vertical), drift q <- q + dt p (horizontal),
   dt = 0.3, from (q, p) = (1.2, -0.8) -- the slow-motion table of §G.4.18 -- over the exact contour.
c  A lopsided posterior: q = log(lambda), lambda ~ Exp(1), so pi(q) ~ exp(q - e^q), V = e^q - q (the Fig A
   example of §G.4.10). Contours are no longer circles; one trajectory at budget H = 2.62 follows one.
d  The projection onto q over (fictitious) time for that posterior: q(t) along successive trajectories, each
   started with a freshly drawn momentum (new budget) -- the actual HMC chain (leapfrog, small step,
   Metropolis-corrected).
e  Where the chain's kept positions land (many transitions) against the target density.
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
rng = np.random.default_rng(19102026)

# lopsided potential
V2 = lambda q: np.exp(q) - q
dV2 = lambda q: np.exp(q) - 1.0


def leapfrog(q, p, dV, eps, L):
    qs, ps = [q], [p]
    p = p - 0.5 * eps * dV(q)
    for i in range(L):
        q = q + eps * p
        if i < L - 1:
            p = p - eps * dV(q)
        qs.append(q); ps.append(p - (0.5 * eps * dV(q) if i == L - 1 else 0.0))
    p = p - 0.5 * eps * dV(q)
    return q, p, np.array(qs), np.array(ps)


fig = plt.figure(figsize=(16, 10.5))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 0.85], hspace=0.42, wspace=0.9)

# ---------------------------------------------------------------- a: Gaussian phase plane
ax = fig.add_subplot(gs[0, 0:2])
qq, pp = np.meshgrid(np.linspace(-3, 3, 400), np.linspace(-3, 3, 400))
H = 0.5 * (qq**2 + pp**2)
levels = [0.5, 1.04, 2.0, 3.5]
ax.contour(qq, pp, H, levels=levels, colors=MUTED, linewidths=1.2)
gq, gp = np.meshgrid(np.linspace(-2.75, 2.75, 12), np.linspace(-2.75, 2.75, 12))
ax.quiver(gq, gp, gp, -gq, color=INK2, angles="xy", scale=38, width=0.004, alpha=0.8)
ax.axvline(0, color=INK2, lw=0.8, ls=":")
ax.plot(0, 0, "o", color=INK, ms=4)
ax.annotate("top of a circle: at the mode,\nmoving right, fastest", (0, np.sqrt(2)), (-2.9, 2.55),
            fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.annotate("right edge: turning point\n(p = 0, pushed back)", (np.sqrt(2), 0), (1.0, -2.75),
            fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_aspect("equal")
ax.set_xlabel("q  (the parameter)"); ax.set_ylabel("p  (its momentum)")
ax.set_title("a  Gaussian posterior: walk along the circles", loc="left", pad=30)
ax.text(0, 1.02, "grey circles = contours of H = q²/2 + p²/2\narrows = (p, −q): the gradient of H turned 90° (length = speed)",
        transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
ax.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- b: kick / drift staircase
axb = fig.add_subplot(gs[0, 2:4])
th = np.linspace(0, 2 * np.pi, 400)
r0 = np.sqrt(2 * 1.04)
axb.plot(r0 * np.cos(th), r0 * np.sin(th), color=MUTED, lw=1.2, ls="--")
q, p, dt = 1.2, -0.8, 0.3
pts = [(q, p)]
for s in range(10):
    p_new = p - dt * q                                   # kick: vertical
    axb.annotate("", (q, p_new), (q, p), arrowprops=dict(arrowstyle="-|>", color=C[1], lw=1.8,
                 shrinkA=0, shrinkB=0, mutation_scale=9))
    q_new = q + dt * p_new                               # drift: horizontal
    axb.annotate("", (q_new, p_new), (q, p_new), arrowprops=dict(arrowstyle="-|>", color=C[0], lw=1.8,
                 shrinkA=0, shrinkB=0, mutation_scale=9))
    q, p = q_new, p_new
    pts.append((q, p))
pts = np.array(pts)
axb.plot(pts[:, 0], pts[:, 1], "o", color=INK, ms=4, zorder=5)
for s, (dx, dy) in {0: (0.12, -0.05), 3: (0.05, -0.32), 8: (-0.55, 0.12)}.items():
    axb.text(pts[s, 0] + dx, pts[s, 1] + dy, f"step {s}", fontsize=8.5, color=INK)
axb.plot([], [], color=C[1], lw=2, label="kick: p ← p − δt·∂V/∂q  (vertical)")
axb.plot([], [], color=C[0], lw=2, label="drift: q ← q + δt·p  (horizontal)")
axb.plot([], [], color=MUTED, lw=1.2, ls="--", label="exact contour, H = 1.04")
axb.legend(loc="upper left", frameon=False, fontsize=8)
axb.axvline(0, color=INK2, lw=0.8, ls=":")
axb.set_xlim(-2.2, 2.2); axb.set_ylim(-2.2, 2.2); axb.set_aspect("equal")
axb.set_xlabel("q"); axb.set_ylabel("p")
axb.set_title("b  The same walk in small steps", loc="left", pad=30)
axb.text(0, 1.02, "start (q, p) = (1.2, −0.8), δt = 0.3\neach kick + drift = one step round the circle",
         transform=axb.transAxes, fontsize=8, color=INK2, va="bottom")
axb.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- c: lopsided phase plane
axc = fig.add_subplot(gs[0, 4:6])
qq, pp = np.meshgrid(np.linspace(-4.5, 2.2, 500), np.linspace(-3, 3, 400))
Hc = V2(qq) + 0.5 * pp**2
axc.contour(qq, pp, Hc, levels=[1.5, 2.0, 2.62, 3.5, 4.5], colors=MUTED, linewidths=1.2)
gq, gp = np.meshgrid(np.linspace(-4.2, 1.8, 12), np.linspace(-2.7, 2.7, 11))
fq, fp = gp, -dV2(gq)
nrm = np.maximum(np.hypot(fq, fp), 1e-9)
axc.quiver(gq, gp, fq / nrm, fp / nrm, color=INK2, angles="xy", scale=30, width=0.004, alpha=0.8)
H0 = 2.62
qT, pT, qs, ps = leapfrog(0.0, np.sqrt(2 * (H0 - 1.0)), dV2, 0.005, 2200)
axc.plot(qs, ps, color=C[2], lw=2.2)
lo = qs.min(); hi = qs.max()
axc.axvline(0, color=INK2, lw=0.8, ls=":")
axc.annotate(f"gentle side: weak force,\nlong slow swing out to q = {lo:.2f}", (lo, 0), (-4.4, -2.85),
             fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axc.annotate(f"steep side: strong force,\nsharp turn at q = {hi:.2f}", (hi, 0), (0.15, 2.35),
             fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
axc.set_xlim(-4.5, 2.2); axc.set_ylim(-3, 3)
axc.set_xlabel("q = log λ"); axc.set_ylabel("p")
axc.set_title("c  A lopsided posterior: the contours bend", loc="left", pad=30)
axc.text(0, 1.02, "V = e^q − q (λ ~ Exp(1) on the log scale)\ngreen = one trajectory at budget H = 2.62; arrows: direction only", transform=axc.transAxes, fontsize=8, color=INK2, va="bottom")
axc.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- d/e: the chain on the lopsided posterior
EPS = 0.02


def hmc_chain(n, q0, record_path=False):
    q = q0; out = np.empty(n); paths = []; acc = 0
    for i in range(n):
        p0 = rng.standard_normal()
        T = rng.uniform(1.0, 3.0)
        L = int(np.ceil(T / EPS))
        q1, p1, qs_, _ = leapfrog(q, p0, dV2, EPS, L)
        dH = (V2(q1) + 0.5 * p1**2) - (V2(q) + 0.5 * p0**2)
        ok = np.log(rng.random()) < -dH
        if record_path:
            paths.append((qs_, ok))
        if ok:
            q = q1; acc += 1
        out[i] = q
    return out, paths, acc / n


_, paths, _ = hmc_chain(9, -2.0, record_path=True)
axd = fig.add_subplot(gs[1, 0:4])
t0 = 0.0
for j, (qs_, ok) in enumerate(paths):
    tt = t0 + EPS * np.arange(len(qs_))
    axd.plot(tt, qs_, color=C[2], lw=1.8)
    axd.plot(tt[0], qs_[0], "o", color=INK, ms=5, zorder=5)
    axd.plot(tt[-1], qs_[-1], "s", color=INK, ms=5, zorder=5)
    t0 = tt[-1]
    axd.axvline(t0, color=MUTED, lw=0.8, ls=":")
axd.axhline(0, color=INK2, lw=0.8, ls=":")
axd.text(0.15, 0.25, "mode", fontsize=8, color=INK2)
axd.set_xlabel("fictitious time t, trajectories laid end to end")
axd.set_ylabel("q")
axd.set_title("d  Projected onto q: the chain is a sequence of swings", loc="left", pad=18)
axd.text(0, 1.02, "● = momentum redrawn (new budget, new contour);  ■ = end of trajectory, momentum thrown away, "
         "q kept;  flow times drawn between 1 and 3", transform=axd.transAxes, fontsize=8, color=INK2,
         va="bottom")
axd.grid(color=GRID, lw=0.8); axd.set_axisbelow(True)

chain, _, acc = hmc_chain(20000, 0.0)
axe = fig.add_subplot(gs[1, 4:6])
grid = np.linspace(-6, 2.5, 400)
axe.hist(chain[500:], bins=70, range=(-6, 2.5), density=True, color=C[0], alpha=0.55,
         label=f"kept q from {len(chain)-500:,} transitions")
axe.plot(grid, np.exp(grid - np.exp(grid)), color=INK, lw=2, label="target π(q) ∝ exp(q − e^q)")
axe.set_xlabel("q"); axe.set_ylabel("density")
axe.set_title("e  …and the kept positions fill in the posterior", loc="left", pad=18)
axe.text(0, 1.02, f"mean {chain[500:].mean():.3f} (truth −0.577),  SD {chain[500:].std():.3f} (truth 1.283)",
         transform=axe.transAxes, fontsize=8, color=INK2, va="bottom")
axe.legend(loc="upper left", frameon=False, fontsize=8)
axe.grid(axis="y", color=GRID, lw=0.8); axe.set_axisbelow(True)

fig.suptitle("Hamilton's equations = walk along the contours of H in (parameter, momentum) space; "
             "project onto the parameter to get the sampler", x=0.06, y=1.0, ha="left", fontsize=13,
             fontweight="bold")
fig.text(0.06, 0.965, "Read top row left to right (the rule, the rule in steps, the rule on a lopsided posterior), "
         "then the bottom row (what the parameter does, and where it ends up).", fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_phase_plane.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"c: budget H={H0}, turning points q = {lo:.3f} (left) and {hi:.3f} (right); "
      f"energy drift along trajectory {np.ptp(V2(qs) + 0.5*ps**2):.2e}")
print(f"e: acceptance {acc:.3f}; mean {chain[500:].mean():.4f} (truth {-np.euler_gamma:.4f}); "
      f"SD {chain[500:].std():.4f} (truth {np.pi/np.sqrt(6):.4f})")
print("wrote", OUT / "fig_hmc_phase_plane.png")
