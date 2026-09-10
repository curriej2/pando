#!/usr/bin/env python3
r"""
Fig 4c -- the worked example: what a loss looks like, how the test works, and
what it declines to call.  Mouse2 clone 76 throughout (the designated crux: near
zero CLONE-level excess, yet 98.9% of the arm's weight, so anything visible
inside it cannot be a clone or batch effect).

 a  the data      : cells x tapes for one clone, rows ordered by the ANCHOR's
                    depth-5 prefix so the clade is contiguous.  Tape 40 shows a
                    solid 38-cell block against a background where it is 93%
                    present.  Tape 102 shows an IDENTICAL block -- and is not
                    called.
 b  the test      : P(K >= k) for a clade of 38 drawn without replacement from
                    the block population, for both tapes.  Same observed k = 38,
                    p = 3e-46 vs p = 0.049, because the two curves differ.
 c  the orthogonal: the co-integrated symbol.  Tape 102's real loss sits on a
                    1,304-cell clade in the same clone, where the symbol its
                    pegRNA writes (AAGCGGA) is depleted in EVERY OTHER tape.

⚠ The 38-cell clade and the 1,304-cell clade share 32 of 38 cells but neither
contains the other -- they are defined by different anchors, so their prefix
partitions cross.  That is "a prefix only partially defines a clade", visible.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
D = np.load("/data1/choij10/justin/tmp/example_data.npz", allow_pickle=False)
DEP = json.load(open("/data1/choij10/justin/tmp/example_depletion.json"))["ev"]

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
BLUE, DARK, NULL = "#2a78d6", "#184f95", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

M = D["M"]; cols = list(D["cols"]); nc = int(D["n_clade"]); ns = int(D["n_samp"])
TT, TN, AN = int(D["tape_tested"]), int(D["tape_near"]), int(D["anchor"])
npop, m = int(D["npop"]), nc

# ------------------------------------------------------------------ panel a
# ⚠ v1 put the three key tapes at the far right among 36 background columns, where
# they were invisible, and stacked their labels on top of each other.  Split into
# two axes: the key tapes WIDE on the left, the background thin on the right.
KEY = [AN, TN, TT]
kj = [cols.index(t) for t in KEY]
bj = [j for j in range(len(cols)) if j not in kj]
fig, (axk, axb) = plt.subplots(1, 2, figsize=(8.6, 5.4), sharey=True,
                               gridspec_kw=dict(width_ratios=[1.2, 2.0], wspace=0.06))
cm = ListedColormap(["#eeeeea", DARK])
axk.imshow(M[:, kj], aspect="auto", interpolation="nearest", cmap=cm, vmin=0, vmax=1)
axb.imshow(M[:, bj], aspect="auto", interpolation="nearest", cmap=cm, vmin=0, vmax=1)
for ax in (axk, axb):
    ax.axhline(nc - 0.5, color=INK, lw=1.6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_edgecolor(GRID)
for i, (t, lab, col, bold) in enumerate((
        (AN, f"anchor {AN}\n(defines\nthe clade)", INK2, False),
        (TN, f"tape {TN}\nnot\ncalled", INK2, False),
        (TT, f"tape {TT}\n\nCALLED", BLUE, True))):
    axk.add_patch(plt.Rectangle((i - 0.5, -0.5), 1, M.shape[0], fill=False,
                                edgecolor=col, lw=2.2, zorder=5))
    axk.annotate(lab, xy=(i, -0.5), xytext=(0, 8), textcoords="offset points",
                 ha="center", va="bottom", fontsize=8.0, color=col,
                 fontweight="bold" if bold else "normal", annotation_clip=False)
axb.annotate(f"{len(bj)} other tapes, for background", xy=(len(bj) / 2, -0.5),
             xytext=(0, 8), textcoords="offset points", ha="center", va="bottom",
             fontsize=8.5, color=INK2, annotation_clip=False)
axk.annotate("the clade\n38 cells", xy=(-0.75, nc / 2), ha="right", va="center",
             fontsize=9, color=INK, fontweight="bold", annotation_clip=False)
axk.annotate(f"rest of the clone\n{ns} of {npop - nc}\nsampled", xy=(-0.75, nc + ns / 2),
             ha="right", va="center", fontsize=9, color=INK2, annotation_clip=False)
fig.suptitle("A tape lost on one clade — and a block that looks identical but is not",
             fontsize=11, fontweight="bold", x=0.13, ha="left", y=1.085)
fig.text(0.13, 1.022, "dark = tape missing in that cell    ·    Mouse2 clone 76, "
         "rows ordered by the anchor's depth-5 prefix", fontsize=8.5, color=INK2, ha="left")
fig.savefig(FIG / "fig4c_a_example_raster.png"); plt.close(fig)

# ------------------------------------------------------------------ panel b
# ⚠ v1 stacked the two curve annotations on the same point and over the title.
fig, ax = plt.subplots(figsize=(6.6, 4.6))
ks = np.arange(0, m + 1)
hh, ll = [], []
for tape, kp, col, lw in ((TT, int(D["kp40"]), BLUE, 2.4), (TN, int(D["kpN"]), INK2, 2.0)):
    sf = stats.hypergeom.sf(ks - 1, npop, kp, m)
    ax.step(ks, sf, where="post", color=col, lw=lw)
    hh.append(plt.Line2D([], [], color=col, lw=lw))
    ll.append(f"tape {tape} — missing in {kp/npop:.0%} of the {npop} cells")
ax.axhline(float(D["bonf"]), color=NULL, lw=1.6, ls=(0, (5, 3)))
ax.annotate(f"threshold  p = {float(D['bonf']):.1e}", xy=(0.6, float(D["bonf"]) * 3),
            fontsize=8.5, color=NULL, va="bottom")
ax.axvline(m, color=INK, lw=0.9, ls=(0, (2, 2)))
ax.annotate("observed\nk = 38", xy=(m - 1.0, 1e-10), ha="right", va="center",
            fontsize=8.5, color=INK)
for tape, p_, col, dy in ((TT, float(D["p40"]), BLUE, -14), (TN, float(D["pN"]), INK2, 12)):
    ax.plot([m], [p_], "o", ms=9, mfc="none", mec=col, mew=2.2)
    ax.annotate(f"p = {p_:.1e}", xy=(m, p_), xytext=(-10, dy), textcoords="offset points",
                fontsize=9.5, color=col, ha="right", fontweight="bold")
ax.set_yscale("log"); ax.set_ylim(1e-50, 20); ax.set_xlim(0, m + 1)
ax.set_xlabel("k = cells in the clade missing the tape")
ax.set_ylabel(r"$P(K \geq k)$ under the permutation null")
ax.set_title("Same clade, same k = 38, opposite verdicts", pad=8)
ax.legend(hh, ll, frameon=False, loc="lower left", fontsize=8.5, labelcolor=INK2,
          handlelength=1.8, bbox_to_anchor=(0.02, 0.02))
ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
fig.savefig(FIG / "fig4c_b_example_test.png"); plt.close(fig)

# ------------------------------------------------------------------ panel c
# ⚠ v1 coloured only AAGCGGA's inside dot blue while the legend implied all of them,
# and used the reserved null-orange for a comparison group.  Neutral = rest of clone,
# blue = inside; AAGCGGA carries emphasis through the connector and the label only.
OUT = "#8f8e86"
fig, ax = plt.subplots(figsize=(6.6, 4.6))
r = DEP["ranked"][:8]
fi = np.array([x["n_in"] / DEP["W_in"] for x in r]) * 100
fo = np.array([x["n_out"] / DEP["W_out"] for x in r]) * 100
lab = [x["symbol"] for x in r]
y = np.arange(len(r))[::-1]
for i in range(len(r)):
    top = lab[i] == "AAGCGGA"
    ax.plot([fo[i], fi[i]], [y[i], y[i]], color=BLUE if top else GRID,
            lw=3.0 if top else 1.8, zorder=3 if top else 1, solid_capstyle="round")
    ax.plot([fo[i]], [y[i]], "o", ms=7, color=OUT, zorder=4)
    ax.plot([fi[i]], [y[i]], "o", ms=7, color=BLUE, zorder=4)
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=9)
for t, l in zip(ax.get_yticklabels(), lab):
    if l == "AAGCGGA":
        t.set_fontweight("bold"); t.set_color(BLUE)
ax.set_xlabel("share of independent post-MRCA writes (%)")
ax.set_title("The co-integrated symbol vanishes too", pad=8)
ax.legend([plt.Line2D([], [], marker="o", ls="none", color=OUT),
           plt.Line2D([], [], marker="o", ls="none", color=BLUE)],
          ["rest of the clone", "inside the clade that lost tape 102"],
          frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK2)
ax.annotate(f"$\\Lambda$ = {DEP['top_lambda']:.1f} nats against a permutation null\n"
            f"whose best of 200 draws reached {DEP['null_max_lambda_max']:.1f}",
            xy=(0.97, 0.26), xycoords="axes fraction", ha="right", fontsize=8.5, color=INK)
ax.grid(axis="x", color=GRID, lw=0.7); ax.set_axisbelow(True)
ax.set_xlim(left=0)
fig.savefig(FIG / "fig4c_c_example_symbol.png"); plt.close(fig)
print(f"wrote fig4c_a/b/c; clade {m} of population {npop}; "
      f"tape {TT} p={float(D['p40']):.2e}, tape {TN} p={float(D['pN']):.2e}")
