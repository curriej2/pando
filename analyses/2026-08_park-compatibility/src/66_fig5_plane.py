#!/usr/bin/env python3
r"""
Fig 5 -- the detection plane, and what a silencing event is.

Three STANDALONE panels (house convention: no composed grid).

 a  the plane        : every scorable (clade, tape) combo as log-density in
                       (pi_hat, score margin); the orange outline is everywhere
                       ANY of the calibration permutations reached.  Density
                       outside it cannot be produced by relabelling cells.
 b  the margin curve : observed vs null counts against the margin, log y.  The
                       two agree below ~0.5 and separate by 10^4 above 2.
 c  the payoff       : pi_hat of the called events against the rate the null
                       model predicted for those same cells.

Palette: dataviz skill reference instance.  Blue = observed/data, ORANGE = the
null (house convention, figure-wide).  Surface #fcfcfb, ink #0b0b0b/#52514e.
"""
import gzip, csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
arm = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"
D = int(sys.argv[2]) if len(sys.argv) > 2 else 6

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
NULL = "#eb6834"                                   # orange = the null, figure-wide
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seqblue", BLUE_RAMP)
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

z = np.load(RES / f"percombo_score_{arm}_d{D}_plane.npz", allow_pickle=False)
PI, MG = z["pi_edges"], z["margin_edges"]
ND = int(z["n_cal_draws"])
HoS = z["H_obs"].astype(float)                     # [stratum, pi, margin]
Ho, Hc = z["H_obs"].sum(0).astype(float), z["H_cal"].sum(0) / ND
# drop the +-1e30 under/overflow bins; nothing sits above margin 8 anyway
inner = slice(1, len(MG) - 1)
MGf = MG[inner.start:inner.stop + 1]
Ho, Hc, HoS = Ho[:, inner], Hc[:, inner], HoS[:, :, inner]
YLO, YHI = -2.0, 8.0        # the bulk below -2 is null-like and is shown in panel b
keep = (MGf[:-1] >= YLO) & (MGf[:-1] < YHI)
Ho, Hc, HoS, MGk = Ho[:, keep], Hc[:, keep], HoS[:, :, keep], MGf[:-1][keep]
ymesh = np.append(MGk, MGk[-1] + (MGk[1] - MGk[0]))

ev = list(csv.DictReader(gzip.open(RES / f"percombo_score_{arm}_d{D}.tsv.gz", "rt"),
                         delimiter="\t"))
e_pi = np.array([float(r["pi_hat"]) for r in ev])
e_mg = np.array([float(r["margin"]) for r in ev])
e_eb = np.array([float(r["expected_rate"]) for r in ev])
e_st = np.array([r["size_stratum"] for r in ev])

# ---------------------------------------------------------------- panel a
# ⚠⚠ The first version drew a single "null ceiling" POOLED over clade sizes. That
# line is not the decision rule -- the threshold is applied per clade-size stratum
# against each combo's OWN null -- and 58 of 127 called events fell below it while
# every one cleared its own stratum's cut.  Root cause: the decision variable
# (clade size) was on neither axis, so NO curve in that plane can be the boundary.
# Faceting puts clade size on the facet, and then the one line per panel IS the rule.
LAB = [str(x) for x in z["size_labels"]]
THRS = np.asarray(z["thresholds"], float)
present = [i for i in range(len(LAB)) if HoS[i].sum() > 0]
vmax = max(HoS[i].max() for i in present)
fig, axes = plt.subplots(1, len(present), figsize=(2.35 * len(present) + 1.5, 3.9),
                         sharey=True, sharex=True)
X, Y = np.meshgrid(PI, ymesh, indexing="ij")
for ax, i in zip(np.atleast_1d(axes), present):
    M = np.ma.masked_where(HoS[i] <= 0, HoS[i])
    pc = ax.pcolormesh(X, Y, M, cmap=CMAP, norm=LogNorm(vmin=1, vmax=vmax),
                       shading="flat", rasterized=True)
    if np.isfinite(THRS[i]):
        ax.axhline(THRS[i], color=NULL, lw=2.0, zorder=4)
        ax.annotate(f"{THRS[i]:.2f}", xy=(0.03, THRS[i]), xytext=(0.03, THRS[i] + 0.35),
                    color=NULL, fontsize=8.5, fontweight="bold")
    m = e_st == LAB[i]
    ax.scatter(e_pi[m], e_mg[m], s=15, facecolor="none", edgecolor=INK,
               linewidths=0.9, zorder=5)
    ax.set_title(f"{LAB[i]} cells   ({int(m.sum())})", fontsize=9.5)
    ax.set_xlim(0, 1); ax.set_ylim(YLO, YHI)
    ax.set_xticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=8.5)
axes0 = np.atleast_1d(axes)
axes0[0].set_ylabel(r"margin  $z_{\rm obs} - \max_b z_b$")
fig.supxlabel(r"$\hat\pi = k/m$   (fraction of the clade missing this tape)",
              fontsize=10, y=0.02)
fig.suptitle("A silencing event is a margin above the bar set for its clade size",
             fontsize=11, fontweight="bold", x=0.055, ha="left", y=0.99)
cb = fig.colorbar(pc, ax=axes0.tolist(), pad=0.012, fraction=0.028)
cb.set_label("combos per bin", color=INK2, fontsize=9)
cb.outline.set_edgecolor(GRID); cb.ax.tick_params(colors=INK2, labelsize=8)
axes0[0].legend([plt.Line2D([], [], color=NULL, lw=2.0),
                 plt.Line2D([], [], marker="o", ls="none", mfc="none", mec=INK, ms=5)],
                ["threshold: expected false $\\leq$ 2", "called event"],
                frameon=False, loc="upper left", fontsize=8, labelcolor=INK2,
                handlelength=1.5, borderpad=0.1)
fig.savefig(FIG / f"fig5a_plane_{arm}.png"); plt.close(fig)

# ---------------------------------------------------------------- panel b
fig, ax = plt.subplots(figsize=(6.0, 4.2))
# survival curves -- exactly the two counts the FDR is a ratio of
Ho_f, Hc_f = z["H_obs"].sum(0)[:, inner], z["H_cal"].sum(0)[:, inner] / ND
MGa = MG[inner.start:inner.stop]
o1 = Ho_f.sum(0)[::-1].cumsum()[::-1]
c1 = Hc_f.sum(0)[::-1].cumsum()[::-1]
ax.step(MGa, np.where(o1 > 0, o1, np.nan), where="post", color="#2a78d6", lw=2.0)
ax.step(MGa, np.where(c1 > 0, c1, np.nan), where="post", color=NULL, lw=2.0)
ax.set_yscale("log")
lastnull = MGa[c1 > 0].max()
ax.axvline(lastnull, color=NULL, lw=0.9, ls=(0, (2, 2)))
ax.annotate(f"no permutation exceeds {lastnull:.1f}", xy=(lastnull, 2.2),
            xytext=(lastnull - 0.4, 2.2), color=NULL, fontsize=8.5,
            va="center", ha="right")
ax.set_xlabel(r"margin  $z_{\rm obs} - \max_b z_b$")
ax.set_ylabel("combos with margin $\\geq$ x")
ax.set_title("Observed and null separate above a margin of ~0.5", pad=10)
ax.set_xlim(-8, 8); ax.set_ylim(0.3, 2e7)
ax.grid(axis="y", color=GRID, lw=0.7)
ax.set_axisbelow(True)
ax.legend([plt.Line2D([], [], color="#2a78d6", lw=2.0),
           plt.Line2D([], [], color=NULL, lw=2.0)],
          ["observed", "permuted (per draw)"], frameon=False, loc="lower left",
          fontsize=9, labelcolor=INK2, handlelength=1.6)
fig.savefig(FIG / f"fig5b_margin_{arm}.png"); plt.close(fig)

# ---------------------------------------------------------------- panel c
fig, ax = plt.subplots(figsize=(6.0, 4.2))
bins = np.linspace(0, 1, 26)
ax.hist(e_eb, bins=bins, color=NULL, alpha=0.85, label="predicted by the null model  "
        r"($\bar e$)")
ax.hist(e_pi, bins=bins, histtype="step", color="#2a78d6", lw=2.2,
        label=r"observed  ($\hat\pi$)")
ax.set_xlabel(r"fraction of the clade missing the tape")
ax.set_ylabel("called events")
ax.set_title(f"What was lost, against what the model predicted — {arm}, {len(ev)} events", pad=10)
ax.set_xlim(0, 1)
ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper center", fontsize=9, labelcolor=INK2)
ax.annotate(f"median $\\hat\\pi$ = {np.median(e_pi):.3f}\nmedian $\\bar e$ = {np.median(e_eb):.3f}",
            xy=(0.03, 0.72), xycoords="axes fraction", fontsize=9, color=INK2, va="top")
fig.savefig(FIG / f"fig5c_completeness_{arm}.png"); plt.close(fig)
print(f"wrote fig5a/b/c for {arm} (d{D}); {len(ev)} events, "
      f"max observed margin {e_mg.max():.2f}, null reach {lastnull:.2f}")
