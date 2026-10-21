"""Teaching figure for Betancourt §5.1 (notes §G.4.29): where the step size enters the leapfrog integrator,
why leapfrog does not drift, and what a divergence is. 1-D posteriors, unit mass.

a  Anatomy of ONE leapfrog step (spring V = q^2/2, eps = 0.6, start (1.2, -0.8)): half kick (vertical,
   length eps/2 x force), drift (horizontal, length eps x velocity), half kick at the new position.
b  The same stretch of time T = 1.8 covered with eps = 0.6 (3 steps) and eps = 0.2 (9 steps): the step size
   is how finely the swing is chopped; number of steps = T / eps; each step = one gradient.
c  Three full swings: Euler (both updates from the old values) spirals outward; leapfrog stays on its loop.
d  Energy along those runs: Euler climbs; leapfrog wobbles in a band that shrinks ~ eps^2.
e  The stability limit on a spring of frequency 1: eps = 1.9 stays bounded, eps = 2.1 explodes.
f  Why divergences happen in practice: a posterior that is gentle on the right (frequency 1) and stiff on
   the left (frequency 5, limit eps < 0.4). With eps = 0.5, fine on the right, the swing explodes on
   entering the stiff side.
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


def subtitle(a, txt):
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.3, color=INK2, va="bottom")


def circle(ax, r, **kw):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(t), r * np.sin(t), **kw)


def leapfrog_sub(q, p, eps, L, force):
    """Return every sub-step point: start, after half kick, after drift, after half kick, ..."""
    pts = [(q, p)]
    for _ in range(L):
        p = p + 0.5 * eps * force(q); pts.append((q, p))
        q = q + eps * p;              pts.append((q, p))
        p = p + 0.5 * eps * force(q); pts.append((q, p))
    return np.array(pts)


def arrow(ax, a, b, col, lw=2.0):
    ax.annotate("", b, a, arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, shrinkA=0, shrinkB=0,
                                          mutation_scale=11))


spring = lambda q: -q                                           # force = d log pi / dq for a standard normal
fig = plt.figure(figsize=(16, 10.8))
gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.32, width_ratios=[1.35, 1, 1])

# ---------------------------------------------------------------- a: one step, annotated
ax = fig.add_subplot(gs[0, 0])
q0, p0, eps = 1.2, -0.8, 0.6
r0 = np.hypot(q0, p0)
circle(ax, r0, color=MUTED, lw=1.2, ls="--")
pts = leapfrog_sub(q0, p0, eps, 1, spring)
arrow(ax, pts[0], pts[1], C[1]); arrow(ax, pts[1], pts[2], C[0]); arrow(ax, pts[2], pts[3], C[1])
ax.plot(*pts[0], "o", color=INK, ms=7, zorder=5); ax.plot(*pts[3], "s", color=INK, ms=7, zorder=5)
qe, pe = q0 * np.cos(eps) + p0 * np.sin(eps), p0 * np.cos(eps) - q0 * np.sin(eps)
ax.plot(qe, pe, "o", mfc="none", mec=INK, ms=9, mew=1.5, zorder=6)
f0, f1 = spring(q0), spring(pts[2][0])
ax.text(pts[0][0] + 0.03, (pts[0][1] + pts[1][1]) / 2, f"① half kick (change speed only)\nΔp = (ε/2) × force(q)\n   = 0.3 × ({f0:+.2f}) = {0.5*eps*f0:+.2f}",
        fontsize=8.5, color=C[1], va="center")
ax.text((pts[1][0] + pts[2][0]) / 2, pts[1][1] + 0.05, f"② drift (change position only)\nΔq = ε × velocity = 0.6 × ({pts[1][1]:+.2f}) = {eps*pts[1][1]:+.2f}",
        fontsize=8.5, color=C[0], ha="center", va="bottom")
ax.text(pts[2][0] + 0.04, pts[2][1] - 0.05, f"③ half kick, force at the NEW q\nΔp = 0.3 × ({f1:+.2f}) = {0.5*eps*f1:+.2f}",
        fontsize=8.5, color=C[1], ha="left", va="top")
ax.annotate("exact position\nafter time ε", (qe, pe), (0.3, -1.52), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.set_xlim(-0.1, 1.85); ax.set_ylim(-1.6, -0.6)
ax.set_xlabel("q  (position)"); ax.set_ylabel("p  (momentum)")
ax.set_title("a  ONE leapfrog step, ε = 0.6", loc="left", pad=28)
subtitle(ax, "● start  ■ end of the step  ○ where the EXACT motion would be after time ε\n"
         "dashed = the exact energy contour;  ε sets every arrow's length")
ax.grid(color=GRID, lw=0.6); ax.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- b: same time, two step sizes
axb = fig.add_subplot(gs[0, 1])
T = 1.8
circle(axb, r0, color=MUTED, lw=1.2, ls="--")
tt = np.linspace(0, T, 200)
axb.plot(q0 * np.cos(tt) + p0 * np.sin(tt), p0 * np.cos(tt) - q0 * np.sin(tt), color=INK, lw=3, alpha=0.25,
         label=f"exact motion for time T = {T}")
for e_, col, lab in [(0.6, C[1], "ε = 0.6  →  3 steps (3 gradients)"), (0.2, C[0], "ε = 0.2  →  9 steps (9 gradients)")]:
    P = leapfrog_sub(q0, p0, e_, int(round(T / e_)), spring)
    axb.plot(P[:, 0], P[:, 1], color=col, lw=1.6, label=lab)
    axb.plot(P[::3, 0], P[::3, 1], "o", color=col, ms=4)
axb.plot(q0, p0, "o", color=INK, ms=7, zorder=5)
axb.set_xlim(-1.75, 1.75); axb.set_ylim(-1.75, 1.75); axb.set_aspect("equal")
axb.set_xlabel("q"); axb.set_ylabel("p")
axb.set_title("b  Same stretch of time, two step sizes", loc="left", pad=28)
subtitle(axb, "number of steps = T / ε.  Small ε: finer staircase, closer to the\nexact arc, more gradients.  Large ε: cheaper, cruder")
axb.legend(frameon=False, fontsize=7.8, loc="upper right")
axb.grid(color=GRID, lw=0.6); axb.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- c/d: Euler vs leapfrog over 3 swings
def euler(q, p, eps, n, force):
    out = [(q, p)]
    for _ in range(n):
        q, p = q + eps * p, p + eps * force(q); out.append((q, p))
    return np.array(out)


eps = 0.3; n3 = int(round(3 * 2 * np.pi / eps))
E = euler(1.0, 0.0, eps, n3, spring)
Lf = leapfrog_sub(1.0, 0.0, eps, n3, spring)[::3]
axc = fig.add_subplot(gs[0, 2])
circle(axc, 1.0, color=MUTED, lw=1.2, ls="--")
axc.plot(E[:, 0], E[:, 1], color=C[1], lw=1.4, label="Euler (q and p updated together)")
axc.plot(Lf[:, 0], Lf[:, 1], color=C[0], lw=2.2, label="leapfrog (half kick, drift, half kick)")
axc.plot(1.0, 0.0, "o", color=INK, ms=7, zorder=5)
axc.set_xlim(-16, 16); axc.set_ylim(-16, 16); axc.set_aspect("equal")
axc.set_xlabel("q"); axc.set_ylabel("p")
axc.set_title("c  Three full swings, ε = 0.3", loc="left", pad=28)
subtitle(axc, "same start (●), same step size, same number of steps:\nEuler's errors add up and it spirals away; leapfrog stays on its loop")
axc.legend(frameon=False, fontsize=7.8, loc="upper right")
axc.grid(color=GRID, lw=0.6); axc.spines[["top", "right"]].set_visible(True)

axd = fig.add_subplot(gs[1, 0])
H = lambda P: 0.5 * P[:, 0]**2 + 0.5 * P[:, 1]**2
tE = eps * np.arange(len(E))
axd.plot(tE, H(E), color=C[1], lw=1.6, label="Euler, ε = 0.3")
for e_, col, ls in [(0.3, C[0], "-"), (0.1, C[2], "-")]:
    nn = int(round(3 * 2 * np.pi / e_))
    P = leapfrog_sub(1.0, 0.0, e_, nn, spring)[::3]
    axd.plot(e_ * np.arange(len(P)), H(P), color=col, lw=1.8, ls=ls, label=f"leapfrog, ε = {e_}")
axd.axhline(0.5, color=INK, lw=1, ls="--")
axd.set_yscale("log")
axd.set_xlabel("time (three full swings = 18.8)"); axd.set_ylabel("energy H  (log scale; true value 0.5)")
axd.set_title("d  Energy along the way", loc="left", pad=28)
subtitle(axd, "Euler: 0.5 → 114.  Leapfrog: wobbles in a band (worst error 0.011 at\nε = 0.3, 0.0012 at ε = 0.1 — 3× smaller step, ≈9× smaller error)")
axd.legend(frameon=False, fontsize=7.8, loc="upper left")
axd.grid(color=GRID, lw=0.6)
ins = axd.inset_axes([0.42, 0.42, 0.55, 0.3])
for e_, col in [(0.3, C[0]), (0.1, C[2])]:
    nn = int(round(3 * 2 * np.pi / e_))
    P = leapfrog_sub(1.0, 0.0, e_, nn, spring)[::3]
    ins.plot(e_ * np.arange(len(P)), H(P), color=col, lw=1.4)
ins.axhline(0.5, color=INK, lw=0.8, ls="--")
ins.set_ylim(0.485, 0.505); ins.tick_params(labelsize=7)
ins.set_title("zoom: leapfrog only", fontsize=8, color=INK2)

# ---------------------------------------------------------------- e: stability limit
axe = fig.add_subplot(gs[1, 1])
for e_, col in [(1.9, C[0]), (2.1, C[1])]:
    P = leapfrog_sub(1.0, 0.0, e_, 30, spring)[::3]
    err = np.abs(H(P) - 0.5) + 1e-16
    axe.plot(np.arange(len(P)), err, "-o", color=col, ms=3.5, lw=1.6,
             label=f"ε = {e_}  ({'stable' if e_ < 2 else 'past the limit: diverges'})")
axe.set_yscale("log"); axe.set_ylim(1e-3, 1e14)
axe.set_xlabel("leapfrog step"); axe.set_ylabel("|energy error|  (log scale)")
axe.set_title("e  The hard limit: ε × (swing frequency) < 2", loc="left", pad=28)
subtitle(axe, "spring with frequency 1: just under the limit the error stays bounded;\njust over it, every step flings the seat further out")
axe.legend(frameon=False, fontsize=7.8, loc="center right")
axe.grid(color=GRID, lw=0.6)

# ---------------------------------------------------------------- f: divergence where the posterior gets stiff
k_left = 25.0                                                     # stiffness on the left: frequency 5
force_mixed = lambda q: -q if q >= 0 else -k_left * q
axf = fig.add_subplot(gs[1, 2])
eps = 0.5
P = leapfrog_sub(2.0, 0.0, eps, 14, force_mixed)
Pk = P[::3]
Hm = lambda q, p: 0.5 * p**2 + (0.5 * q**2 if q >= 0 else 0.5 * k_left * q**2)
blow = next((i for i in range(len(Pk)) if Hm(*Pk[i]) > 50), len(Pk))
show = Pk[:min(blow + 1, len(Pk))]
qq, pp = np.meshgrid(np.linspace(-2.5, 2.5, 400), np.linspace(-6, 6, 400))
Hg = 0.5 * pp**2 + np.where(qq >= 0, 0.5 * qq**2, 0.5 * k_left * qq**2)
axf.contour(qq, pp, Hg, levels=[2.0], colors=MUTED, linewidths=1.2, linestyles="--")
axf.axvspan(-2.5, 0, color=BAND, zorder=0)
axf.plot(show[:, 0], show[:, 1], "-o", color=C[1], ms=4, lw=1.6)
axf.plot(2.0, 0.0, "o", color=INK, ms=7, zorder=5)
axf.text(-2.4, 5.4, "stiff side\n(freq. 5: needs ε < 0.4)", fontsize=8.3, color=INK, va="top")
axf.text(0.2, 5.4, "gentle side\n(freq. 1: ε = 0.5 fine)", fontsize=8.3, color=INK, va="top")
axf.annotate("flung out:\na DIVERGENCE", (show[-1][0], np.clip(show[-1][1], -5.8, 5.8)), (0.6, -4.8), fontsize=8.5,
             color=C[1], fontweight="bold", arrowprops=dict(arrowstyle="-", color=C[1], lw=0.8), annotation_clip=False)
axf.set_xlim(-2.5, 2.5); axf.set_ylim(-6, 6)
axf.set_xlabel("q"); axf.set_ylabel("p")
axf.set_title("f  Why divergences happen: stiffness varies", loc="left", pad=28)
subtitle(axf, "one step size, tuned for the bulk, ε = 0.5: fine on the gentle\nside, past the limit once the swing enters the stiff side")
axf.grid(color=GRID, lw=0.6); axf.spines[["top", "right"]].set_visible(True)

fig.suptitle("The leapfrog integrator: where the step size ε goes, why it does not drift, and what a divergence is",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.965, "Each step:  ① p ← p + (ε/2)·force(q)   ② q ← q + ε·p/mass   ③ p ← p + (ε/2)·force(new q)."
         "   force = gradient of the log posterior (for us: the engine).  1-D posteriors, mass 1.", fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_leapfrog.png", dpi=200, bbox_inches="tight")
plt.close(fig)

a = leapfrog_sub(q0, p0, 0.6, 1, spring)
print("a: start", a[0], "after half kick", np.round(a[1], 3), "after drift", np.round(a[2], 3), "end", np.round(a[3], 3),
      "exact after eps", (round(qe, 3), round(pe, 3)))
print(f"f: energy at steps {[round(Hm(*x),2) for x in Pk[:min(blow+2,len(Pk))]]}")
print("wrote", OUT / "fig_hmc_leapfrog.png")
