"""Teaching figure (notes §G.4.33-34): leapfrog follows a slightly wider, shallower "shadow" landscape, and the
exp(-H) pick weights undo that. 1-D standard normal posterior (V = q^2/2), unit mass, deliberately coarse step
eps = 0.9 so the effect is visible. For this posterior leapfrog conserves exactly
    H_shadow = p^2/2 + (1 - eps^2/4) q^2/2,
so its swing is the exact swing of the wider valley V_shadow = (1 - eps^2/4) q^2/2.

a  Phase plane: the true energy contour through the start vs the leapfrog points (one long run), which lie
   exactly on the shadow ellipse; points coloured by their exp(-H) weight.
b  The two landscapes: true V vs shadow V, with the energy budget of this trajectory -> where each swing turns.
c  The pick weight exp(-H) along the trajectory as a function of position, against exp(-eps^2 q^2 / 8).
d  The q-distributions: what a uniform pick samples (the shadow posterior), the weight, and the result of
   weighting (the true posterior) -- long runs of both rules.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#a3a29c", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
rng = np.random.default_rng(22102026)


def subtitle(a, txt):
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.4, color=INK2, va="bottom")


EPS = 0.9
s2 = 1 - EPS**2 / 4                                        # shadow stiffness factor
H = lambda q, p: 0.5 * q * q + 0.5 * p * p
Hs = lambda q, p: 0.5 * p * p + 0.5 * s2 * q * q


def lf(q, p, e):
    p = p - 0.5 * e * q; q = q + e * p; p = p - 0.5 * e * q
    return q, p


q0, p0 = 1.2, -0.8
H0 = H(q0, p0); Hs0 = Hs(q0, p0)
pts = [(q0, p0)]; q, p = q0, p0
for i in range(60):
    q, p = lf(q, p, EPS); pts.append((q, p))
pts = np.array(pts)
Hpts = H(pts[:, 0], pts[:, 1])
w = np.exp(-(Hpts - H0))                                   # weight relative to the start

fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.25)

# ---------------------------------------------------------------- a: phase plane
ax = fig.add_subplot(gs[0, 0])
t = np.linspace(0, 2 * np.pi, 400)
rt = np.sqrt(2 * H0)
ax.plot(rt * np.cos(t), rt * np.sin(t), color=INK, lw=2, label=f"TRUE swing (energy contour, H = {H0:.2f})")
ax.plot(np.sqrt(2 * Hs0 / s2) * np.cos(t), np.sqrt(2 * Hs0) * np.sin(t), color=C[1], lw=1.6, ls="--",
        label="shadow swing (wider valley) — leapfrog follows this exactly")
cmap = LinearSegmentedColormap.from_list("seq", SEQ)
sc = ax.scatter(pts[:, 0], pts[:, 1], c=w, cmap=cmap, s=42, edgecolor=INK, lw=0.4, zorder=5,
                label=f"leapfrog points (ε = {EPS}), coloured by weight e^(−H)")
ax.plot(q0, p0, "o", mfc="none", mec=INK, ms=12, mew=1.8, zorder=6)
ax.annotate("start", (q0, p0), (1.45, -1.6), fontsize=8.5, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
qmax_true, qmax_sh = rt, np.sqrt(2 * Hs0 / s2)
for sgn in (1, -1):
    ax.annotate("", (sgn * qmax_sh, 0), (sgn * qmax_true, 0),
                arrowprops=dict(arrowstyle="-|>", color=C[1], lw=1.5, mutation_scale=10))
ax.text(qmax_sh + 0.05, 0.12, f"overshoot\n{qmax_true:.2f} → {qmax_sh:.2f}", fontsize=8, color=C[1])
cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("weight e^(−H) relative to the start", fontsize=8.5)
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.0, 2.0); ax.set_aspect("equal")
ax.set_xlabel("q  (position)"); ax.set_ylabel("p  (momentum)")
ax.set_title("a  Leapfrog swings on a slightly WIDER loop", loc="left", pad=28)
subtitle(ax, "points far out in q (beyond the true loop) are pale = low weight;\n"
         "points near the middle are heavier — the weight depends on position")
ax.legend(frameon=False, fontsize=7.6, loc="upper left")
ax.grid(color=GRID, lw=0.6); ax.spines[["top", "right"]].set_visible(True)

# ---------------------------------------------------------------- b: the two landscapes
axb = fig.add_subplot(gs[0, 1])
xq = np.linspace(-2.4, 2.4, 400)
axb.plot(xq, 0.5 * xq**2, color=INK, lw=2.2, label="TRUE landscape  V = q²/2  (posterior turned upside down)")
axb.plot(xq, 0.5 * s2 * xq**2, color=C[1], lw=2, ls="--", label=f"shadow landscape  (1 − ε²/4)·q²/2  — wider, shallower")
axb.axhline(H0, color=MUTED, lw=1.2)
axb.text(-2.35, H0 + 0.05, f"energy budget of this trajectory  H = {H0:.2f}", fontsize=8.3, color=INK2)
for qm, col in [(qmax_true, INK), (qmax_sh, C[1])]:
    for sgn in (1, -1):
        axb.plot([sgn * qm] * 2, [0, 0.5 * (1 if col == INK else s2) * qm**2], color=col, lw=1, ls=":")
axb.text(0, -0.17, f"true swing turns at ±{qmax_true:.2f};   leapfrog's (shadow) swing turns at ±{qmax_sh:.2f}",
         fontsize=8.3, ha="center", color=INK)
axb.set_ylim(-0.3, 2.6); axb.set_xlim(-2.4, 2.4)
axb.set_xlabel("q"); axb.set_ylabel("height V  (nats)")
axb.set_title("b  …because it rolls on a slightly shallower valley", loc="left", pad=28)
subtitle(axb, f"leapfrog with ε = {EPS} behaves EXACTLY like the true swing on the dashed valley:\n"
         "a shallower valley lets the same energy carry the seat further out")
axb.legend(frameon=False, fontsize=7.8, loc="upper center")
axb.grid(color=GRID, lw=0.6); axb.set_axisbelow(True)

# ---------------------------------------------------------------- c: weight vs position
axc = fig.add_subplot(gs[1, 0])
axc.scatter(pts[:, 0], w, c=w, cmap=cmap, s=40, edgecolor=INK, lw=0.4, zorder=5, label="leapfrog points: e^(−H) relative to start")
axc.plot(xq, np.exp(-EPS**2 * (xq**2 - q0**2) / 8), color=C[1], lw=1.8,
         label="e^(−ε²(q² − q₀²)/8)  =  true density ÷ shadow density (rescaled)")
axc.axvline(0, color=INK2, lw=0.8, ls=":")
axc.set_xlim(-2.3, 2.3)
axc.set_xlabel("position q of the leapfrog point"); axc.set_ylabel("pick weight e^(−H)  (start = 1)")
axc.set_title("c  The weight is a function of POSITION", loc="left", pad=28)
subtitle(axc, "every point sits on the curve: the weight is the ratio (true posterior density) ÷ (leapfrog's\n"
         "shadow density) at that q — an importance weight that pulls the over-wide swing back in")
axc.legend(frameon=False, fontsize=7.8, loc="lower center")
axc.grid(color=GRID, lw=0.6); axc.set_axisbelow(True)

# ---------------------------------------------------------------- d: the q-distributions
def run(rule, n_it=80000, L=4):
    q = 0.0; out = np.empty(n_it)
    for it in range(n_it):
        p = rng.standard_normal(); k = rng.integers(0, L + 1)
        Q = [q]; P = [p]; qq, pp = q, p
        for j in range(k):
            qq, pp = lf(qq, pp, -EPS); Q.insert(0, qq); P.insert(0, pp)
        qq, pp = q, p
        for j in range(L - k):
            qq, pp = lf(qq, pp, EPS); Q.append(qq); P.append(pp)
        Q = np.array(Q); P = np.array(P)
        if rule == "uniform":
            j = rng.integers(0, len(Q))
        else:
            ww = np.exp(-H(Q, P)); ww /= ww.sum(); j = rng.choice(len(Q), p=ww)
        q = Q[j]; out[it] = q
    return out[1000:]


qu = run("uniform"); qw = run("weighted")
axd = fig.add_subplot(gs[1, 1])
xg = np.linspace(-4.5, 4.5, 400); bins = np.linspace(-4.5, 4.5, 61)
sd_sh = 1 / np.sqrt(s2)
axd.plot(xg, stats.norm.pdf(xg), color=INK, lw=2.4, label="TRUE posterior (SD 1)")
axd.plot(xg, stats.norm.pdf(xg, scale=sd_sh), color=C[1], lw=1.8, ls="--", label=f"shadow posterior (SD {sd_sh:.3f})")
for arr, col, lab in [(qu, C[1], f"uniform pick, long run: SD {qu.std():.3f}"), (qw, C[0], f"e^(−H) pick, long run: SD {qw.std():.3f}")]:
    h, e = np.histogram(arr, bins=bins, density=True)
    axd.step(e[:-1], h, where="post", color=col, lw=1.4, alpha=0.9, label=lab)
axd.set_xlabel("q"); axd.set_ylabel("density")
axd.set_title("d  Uniform pick samples the WRONG valley; e^(−H) fixes it", loc="left", pad=28)
subtitle(axd, "uniform picks are photos of the shadow swing, so they follow the shadow posterior;\n"
         "weighting by e^(−H) = by true ÷ shadow density converts them back to the true posterior")
axd.legend(frameon=False, fontsize=7.8, loc="upper right")
axd.grid(axis="y", color=GRID, lw=0.6); axd.set_axisbelow(True)

fig.suptitle("What the e^(−H) weight corrects: leapfrog swings on a slightly wider, shallower valley than the true one",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.965, f"1-D standard normal posterior, deliberately coarse step ε = {EPS}.  For this posterior leapfrog "
         "conserves EXACTLY the shadow energy p²/2 + (1 − ε²/4)·q²/2, so the energy error at a point depends only on its position.",
         fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_shadow.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"shadow conserved along 60 steps: max |H_shadow - H_shadow0| = {np.max(np.abs(Hs(pts[:,0],pts[:,1]) - Hs0)):.2e}")
print(f"turning points: true {qmax_true:.3f}, shadow {qmax_sh:.3f}; max |q| of leapfrog points {np.abs(pts[:,0]).max():.3f}")
print(f"weight vs formula max abs diff {np.max(np.abs(w - np.exp(-EPS**2*(pts[:,0]**2-q0**2)/8))):.2e}")
print(f"d: uniform SD {qu.std():.4f} (shadow {sd_sh:.4f}); weighted SD {qw.std():.4f} (truth 1)")
print("wrote", OUT / "fig_hmc_shadow.png")
