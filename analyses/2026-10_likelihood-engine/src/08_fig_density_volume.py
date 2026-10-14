"""Teaching figures for Betancourt §1.2-1.4: probability = density x volume.

Fig 1  A 2-D standard Gaussian cut into rings of width 0.5 around the mode: density per unit area
       (falls outward), ring area (grows outward), probability in the ring = their product (peaks in a
       middle ring). Monte Carlo draws are counted per ring as a check.
Fig 2  The same product in D dimensions: at distance r from the mode, probability is proportional to
       e^{-r^2/2} (density) x r^{D-1} (size of the shell), which peaks at r* = sqrt(D-1).
Notes: sciphy_notes.md §G.4.11.
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

C = ["#2a78d6", "#eb6834", "#1baf7a"]                     # palette slots 1-3
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]  # blue 100->700
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
rng = np.random.default_rng(13102026)

# ---------------------------------------------------------------- Fig 1: rings in 2-D
W, NR = 0.5, 7
edges = np.arange(NR + 1) * W
r_in, r_out = edges[:-1], edges[1:]
p_ring = np.exp(-r_in**2 / 2) - np.exp(-r_out**2 / 2)        # exact P(r_in < R < r_out), 2-D
area = np.pi * (r_out**2 - r_in**2)
dens = p_ring / area                                          # average density per unit area
N_DRAW = 4000
pts = rng.standard_normal((N_DRAW, 2))
big = rng.standard_normal((2_000_000, 2))
frac_mc = np.histogram(np.hypot(*big.T), bins=edges)[0] / len(big)

fig = plt.figure(figsize=(12, 6.2))
gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.05], hspace=0.65, wspace=0.28)
ax = fig.add_subplot(gs[:, 0])
g = np.linspace(-3.6, 3.6, 600)
X, Y = np.meshgrid(g, g)
Z = np.exp(-(X**2 + Y**2) / 2) / (2 * np.pi)
cmap = LinearSegmentedColormap.from_list("seq", [SURF] + SEQ)
im = ax.imshow(Z, extent=[g[0], g[-1], g[0], g[-1]], origin="lower", cmap=cmap)
th = np.linspace(0, 2 * np.pi, 400)
for r in edges[1:]:
    ax.plot(r * np.cos(th), r * np.sin(th), color=INK2, lw=0.8)
ax.scatter(*pts.T, s=2, color=INK, alpha=0.35, lw=0)
ax.plot(0, 0, "+", color=SURF, ms=12, mew=2.5)
ax.text(0.12, 0.12, "mode", color=SURF, fontsize=9, fontweight="bold", va="bottom")
for i in range(NR):
    rm = (r_in[i] + r_out[i]) / 2
    ax.text(-rm * np.cos(np.pi / 4), -rm * np.sin(np.pi / 4), f"{i+1}", color=SURF if i < 3 else INK,
            fontsize=8, ha="center", va="center", fontweight="bold")
ax.set_aspect("equal")
ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.6, 3.6)
ax.set_xlabel("$q_1$ (SDs from the mode)"); ax.set_ylabel("$q_2$ (SDs from the mode)")
ax.set_title("a  2-D standard Gaussian, cut into rings of width 0.5", loc="left", pad=20)
ax.text(0, 1.015, f"shading = density (darker = higher);  dots = {N_DRAW:,} draws",
        transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
ax.spines[["top", "right"]].set_visible(True)

ring_lab = [f"{i+1}" for i in range(NR)]
xb = np.arange(NR)
panels = [
    (dens, C[0], "b  Density in the ring (probability per unit area)", "per unit area",
     "highest at the centre, falls outward"),
    (area, C[1], "c  Area of the ring", "area (SD²)", "grows outward: ring k has area ∝ (2k − 1)"),
    (p_ring, C[2], "d  Probability in the ring = density × area", "probability",
     "peaks in rings 2-3 (r = 0.5-1.5), not at the mode; ring 1 holds 12%"),
]
for row, (vals, col, title, ylab, note) in enumerate(panels):
    a = fig.add_subplot(gs[row, 1])
    a.bar(xb, vals, color=col, width=0.78)
    if row == 2:
        a.plot(xb, frac_mc, "o", color=INK, ms=4, label="fraction of 2,000,000 draws")
        a.legend(frameon=False, fontsize=8, loc="upper right")
        for x_, v in zip(xb, vals):
            a.text(x_, v, f"{v:.0%}" if v >= 0.01 else f"{v:.1%}", ha="center", va="bottom",
                   fontsize=8, color=INK2)
    a.set_title(title, loc="left", fontsize=10, pad=15)
    a.text(0, 1.03, note, transform=a.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK2)
    a.set_ylabel(ylab, fontsize=9)
    a.set_xticks(xb)
    a.set_xticklabels([f"{i+1}\n{r_in[i]:.1f}–{r_out[i]:.1f}" for i in range(NR)] if row == 2 else ring_lab,
                      fontsize=8)
    a.grid(axis="y", color=GRID, lw=0.8); a.set_axisbelow(True)
    a.set_ylim(0, vals.max() * 1.3)
    if row == 2:
        a.set_xlabel("ring (number, and distance from the mode in SDs)")
fig.suptitle("Where is the probability?  Density is highest at the mode, but the mode's neighbourhood is small",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.savefig(OUT / "fig_volume_1_rings2d.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig 2: shells in D dimensions
DS = [1, 2, 10, 100]
r = np.linspace(0, 13, 4001)
fig, axs = plt.subplots(1, 4, figsize=(14, 4.4), sharey=True)
for a, D in zip(axs, DS):
    d_ = np.exp(-r**2 / 2)
    s_ = r ** (D - 1) if D > 1 else np.ones_like(r)
    prod = stats.chi.pdf(r, D)                               # exact density of the distance R
    for y, col, lab in [(d_ / d_.max(), C[0], "density"),
                        (s_ / s_.max(), C[1], "shell size"),
                        (prod / prod.max(), C[2], "probability\n(product)")]:
        a.plot(r, y, color=col, lw=2.5 if col == C[2] else 2)
    a.fill_between(r, 0, prod / prod.max(), color=C[2], alpha=0.15, lw=0)
    rs = np.sqrt(D - 1)
    a.axvline(rs, color=INK2, lw=1, ls="--")
    a.set_title(f"D = {D}", loc="left")
    right = D >= 100
    a.text(rs - 0.25 if right else rs + 0.25, 1.06,
           f"peak at r = √(D−1) = {rs:.1f}" if D > 1 else "peak at the mode (r = 0)",
           fontsize=8.5, color=INK, va="bottom", ha="right" if right else "left")
    a.set_xlabel("r = distance from the mode (SDs)")
    a.set_xlim(0, 13); a.set_ylim(0, 1.2)
    a.grid(color=GRID, lw=0.8); a.set_axisbelow(True)
    mc = np.linalg.norm(rng.standard_normal((200_000, D)), axis=1)
    a.set_title(f"draws: median r = {np.median(mc):.2f}", loc="right", fontsize=8.5, color=INK2,
                fontweight="normal")
axs[0].set_ylabel("each curve scaled to its own maximum on this axis")
# direct labels on the first panel where all three are separated
axs[0].text(2.2, 0.2, "density\n∝ e^(−r²/2)", color=C[0], fontsize=9)
axs[1].text(9.0, 0.78, "shell size\n∝ r^(D−1)", color=C[1], fontsize=9)
axs[2].text(4.3, 0.70, "probability\n= product", color=C[2], fontsize=9, fontweight="bold")
handles = [plt.Line2D([], [], color=c, lw=2.5) for c in C]
fig.legend(handles, ["density at distance r  (e^(−r²/2))", "size of the shell at distance r  (∝ r^(D−1))",
                     "probability at distance r  (density × shell size)"],
           loc="lower left", bbox_to_anchor=(0.06, -0.08), ncol=3, frameon=False, fontsize=9)
fig.suptitle("The same product in D dimensions: more directions push the probability away from the mode",
             x=0.06, y=1.02, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.945, "Standard Gaussian in D dimensions.  Reading rule: the shaded curve is where draws actually "
         "land; the dashed line is its peak.  At D = 1 it sits on the mode; by D = 100 it is ~10 SDs away.",
         fontsize=9, color=INK2)
fig.subplots_adjust(top=0.8, wspace=0.12)
fig.savefig(OUT / "fig_volume_2_shellsD.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("Fig 1 rings (width 0.5): exact P per ring", np.round(p_ring, 4))
print("                         MC fraction      ", np.round(frac_mc, 4))
print("                         density/unit area", np.round(dens, 4), " area", np.round(area, 3))
for D in DS:
    mc = np.linalg.norm(rng.standard_normal((200_000, D)), axis=1)
    h, e = np.histogram(mc, bins=200)
    print(f"Fig 2 D={D:3d}: r* = {np.sqrt(D-1):.2f}  MC peak of R ~ {0.5*(e[h.argmax()]+e[h.argmax()+1]):.2f}"
          f"  MC median {np.median(mc):.2f}")
print("wrote", OUT / "fig_volume_1_rings2d.png", OUT / "fig_volume_2_shellsD.png")
