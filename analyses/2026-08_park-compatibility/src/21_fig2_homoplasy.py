#!/usr/bin/env python3
r"""FIGURE 2 -- homoplasy is rare and concentrated.

Two versions, selected by argv[1]:
  simple  (default) the count-only estimator: E[s|m] = sum_i [1-(1-xi_i)^m], inverted
  mle               the Poissonised MLE, which also uses WHICH symbols were seen

The simple version is the one to present; the MLE is held in reserve.
"""
import gzip, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]
MODE = (sys.argv[1] if len(sys.argv) > 1 else "simple").lower()
MLE = MODE == "mle"

xi_all = json.loads((RES / "xi_vectors.json").read_text())
q = float(np.mean([((np.array(list(xi_all[a]["xi"].values())) /
                     np.array(list(xi_all[a]["xi"].values())).sum()) ** 2).sum() for a in ARMS]))
mstar = np.sqrt(2 / q)
xi = np.array(list(xi_all["Mouse1"]["xi"].values())); xi = xi / xi.sum()

rows = []
for a in ARMS:
    f = RES / (f"m_mle_{a}.tsv.gz" if MLE else f"m_estimates_{a}.tsv.gz")
    ic = (3, 4, 7, 9) if MLE else (4, 5, 6, 7)      # n_next, s, m, rec
    with gzip.open(f, "rt") as fh:
        next(fh)
        for line in fh:
            p = line.split("\t")
            rows.append((int(p[ic[0]]), int(p[ic[1]]), float(p[ic[2]]), float(p[ic[3]])))
A = np.array(rows); nn, s_obs, m, rec = A[:, 0], A[:, 1], A[:, 2], A[:, 3]

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE, AQUA = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 13.5})
fig, ax = plt.subplots(2, 2, figsize=(13.6, 9.2))
axA, axB, axC, axD = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]

# ---------------- (a) how m-hat is obtained
if not MLE:
    grid = np.unique(np.round(np.logspace(0, 4, 900)).astype(int)).astype(float)
    Es = (1.0 - np.power(1.0 - xi[None, :], grid[:, None])).sum(axis=1)
    s_ex = 44.0; m_ex = float(np.interp(s_ex, Es, grid))
    axA.plot(grid, Es, color=BLUE, lw=2.6)
    axA.plot([grid[0], m_ex], [s_ex, s_ex], color=ORANGE, lw=1.6, ls=(0, (5, 3)))
    axA.plot([m_ex, m_ex], [0, s_ex], color=ORANGE, lw=1.6, ls=(0, (5, 3)))
    axA.plot([m_ex], [s_ex], "o", ms=10, color=ORANGE, markeredgecolor=SURF,
             markeredgewidth=2, zorder=5)
    axA.text(1.35, s_ex + 3, "count $s$", color=ORANGE, fontsize=12)
    axA.text(m_ex * 1.3, 4, "read off $\\hat m$", color=ORANGE, fontsize=12)
    axA.text(0.035, 0.965, r"$\mathbb{E}[s\,|\,m]=\sum_i\left(1-(1-\xi_i)^m\right)$",
             transform=axA.transAxes, ha="left", va="top", fontsize=13, color=INK,
             bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))
    axA.set_xscale("log"); axA.set_xlim(1, 1e4); axA.set_ylim(0, xi.size * 1.08)
    axA.set_xlabel("$m$   independent write events  (not observable)")
    axA.set_ylabel("$s$   distinct symbols  (countable)")
    axA.set_title("a   Counting symbols recovers the number of events", loc="left", pad=10)
else:
    mm = np.logspace(0, 3.2, 400)
    order = np.argsort(xi)[::-1]
    for lab, idx, col in (("common symbols", order[:30], BLUE),
                          ("rare symbols", order[-30:], AQUA)):
        x = xi[idx]; W = x.sum()
        ll = np.array([np.sum(np.log1p(-np.exp(-t * x))) - t * (1 - W) for t in mm])
        ll = ll - ll.max()
        axA.plot(mm, ll, color=col, lw=2.6, label=f"$s$=30, {lab}")
        axA.plot([mm[ll.argmax()]], [0], "o", ms=10, color=col,
                 markeredgecolor=SURF, markeredgewidth=2, zorder=5)
    axA.legend(frameon=False, fontsize=11, loc="lower right")
    axA.set_xscale("log"); axA.set_ylim(-14, 3)
    axA.set_xlabel("$m$")
    axA.set_ylabel("log-likelihood  (relative)")
    axA.set_title("a   The MLE also uses WHICH symbols were seen", loc="left", pad=10)
    axA.text(0.03, 0.05,
             r"$\log L(m)=\sum_{i\in A}\log\left(1-e^{-m\xi_i}\right)-m\,(1-W_A)$",
             transform=axA.transAxes, fontsize=12.5, color=INK,
             bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))

# ---------------- (b) distribution against the threshold, LINEAR y
bins = np.logspace(0, np.log10(max(m.max(), 12)) + 0.05, 34)
axB.hist(m, bins=bins, color=BLUE, edgecolor=SURF, linewidth=0.4)
axB.axvline(mstar, color=ORANGE, lw=2.2, ls=(0, (5, 3)))
axB.set_xscale("log")
axB.set_xlabel("$\\hat m$   inferred write events per node")
axB.set_ylabel("prefix nodes")
axB.set_title("b   Almost all nodes fall below the collision threshold", loc="left", pad=10)
axB.text(mstar * 1.45, axB.get_ylim()[1] * 0.93,
         f"{100*(m<mstar).mean():.0f}% of nodes\nbelow $m^*$", color=ORANGE,
         fontsize=12.5, va="top")
axB.text(0.985, 0.62,
         "expected collisions among $m$ events\n"
         r"$=\binom{m}{2}q\;\;\Rightarrow\;\;$" + "set to 1:\n"
         r"$m^*=\sqrt{2/q}=$" + f"{mstar:.0f}",
         transform=axB.transAxes, ha="right", va="top", fontsize=12, color=INK,
         bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=AXIS, lw=0.8))

# ---------------- (c) recurrences vs clade size
edges = [3, 4, 6, 11, 21, 51, 101, 501, 10**9]
cx, cy = [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (nn >= lo) & (nn < hi)
    if sel.sum() < 20: continue
    cx.append(np.median(nn[sel])); cy.append(rec[sel].mean())
axC.plot(cx, cy, color=BLUE, lw=2.6, marker="o", ms=9,
         markerfacecolor=BLUE, markeredgecolor=SURF, markeredgewidth=2)
axC.set_xscale("log"); axC.set_yscale("log")
axC.set_xlabel("clade size   (cells carrying the prefix)")
axC.set_ylabel("mean recurrences per node\n$\\langle\\,\\hat m - s\\,\\rangle$")
axC.set_title("c   Homoplasy is a large-clade phenomenon", loc="left", pad=10)

# ---------------- (d) concentration
o = np.sort(rec)[::-1]; cum = np.cumsum(o) / o.sum()
fr = 100 * np.arange(1, o.size + 1) / o.size
axD.plot(fr, 100 * cum, color=BLUE, lw=2.6)
axD.plot([0, 100], [0, 100], color=AXIS, lw=1.8, ls=(0, (5, 3)))
k = max(int(o.size * 0.01) - 1, 0)
axD.plot([1], [100 * cum[k]], "o", ms=10, color=ORANGE,
         markeredgecolor=SURF, markeredgewidth=2, zorder=5)
axD.text(6, 100 * cum[k] - 4, f"top 1% of nodes\ncarry {100*cum[k]:.0f}%",
         color=ORANGE, fontsize=12, va="top")
axD.set_xlim(0, 100); axD.set_ylim(0, 102)
axD.set_xlabel("% of nodes, ranked by recurrences")
axD.set_ylabel("% of all homoplasy")
axD.set_title("d   …and it is concentrated in a few nodes", loc="left", pad=10)

fig.tight_layout()
out = FIG / (f"fig2_homoplasy_{'mle' if MLE else 'simple'}.png")
fig.savefig(out, dpi=200, facecolor=SURF)
print(f"[{MODE}] nodes {len(A):,}  q={q:.5f}  m*={mstar:.2f}  "
      f"below={100*(m<mstar).mean():.2f}%  top1%={100*cum[k]:.1f}%  -> {out.name}")
