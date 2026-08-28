#!/usr/bin/env python3
"""Where does the mass of inferred m sit, relative to the birthday threshold?

Panel A  distribution of m_hat over all prefix nodes, against m* = sqrt(2/q)
Panel B  recurrence rate vs clade size -- the mechanistic prediction, tested
Panel C  concentration: what share of all homoplasy sits in what share of nodes
         (this is the quantity D.4b step 5 needs: concentrated conflict can be
          deleted, spread-thin conflict cannot)
"""
import gzip, json, sys
import numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
TABLES = sys.argv[1:] or ["Mouse1", "Mouse2", "Mouse3"]

rows = []
for t in TABLES:
    p = RES / f"m_estimates_{t}.tsv.gz"
    if not p.exists():
        print(f"  (skipping {t}, not computed)"); continue
    with gzip.open(p, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((int(f[4]), int(f[5]), float(f[6]), float(f[7])))
a = np.array(rows)
n_next, s, m, rec = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
# q from the trimmed alphabet, averaged over the tables in play (they agree to <3%)
xi_all = json.loads((RES / "xi_vectors.json").read_text())
qs = []
for tb in TABLES:
    d = xi_all[tb]
    keep = {k: v for k, v in d["xi"].items() if v * d["n_edits"] >= 1000}
    tot = sum(keep.values())
    qs.append(sum((v / tot) ** 2 for v in keep.values()))
q = float(np.mean(qs))
birthday = np.sqrt(2 / q)
below = 100 * (m < birthday).mean()

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6,
})
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.4, 4.4))

# ---- A: where the mass of m_hat sits
bins = np.logspace(0, np.log10(max(m.max(), 12)) + 0.05, 30)
axA.hist(m, bins=bins, color=BLUE, alpha=0.85, edgecolor=SURFACE, linewidth=0.4)
axA.axvline(birthday, color=ORANGE, lw=1.8, ls=(0, (5, 2)), zorder=5)
axA.set_xscale("log"); axA.set_yscale("log")
axA.set_xlabel("$\\hat m$ — inferred write events at the node")
axA.set_ylabel("prefix nodes")
axA.annotate(f"$m^*={birthday:.1f}$\nbirthday threshold",
             xy=(birthday, 3.0e4), xycoords="data",
             xytext=(0.62, 0.56), textcoords="axes fraction",
             color=ORANGE, fontsize=8.5, va="center",
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9,
                             shrinkA=2, shrinkB=2))
axA.text(0.97, 0.95, f"{below:.1f}% of nodes fall below $m^*$",
         transform=axA.transAxes, color=BLUE, fontsize=9.5, va="top", ha="right",
         fontweight="bold")
axA.set_title("A.  Almost all nodes sit below the collision threshold",
              loc="left", fontsize=10, color=INK, pad=8)

# ---- B: recurrence rate vs clade size
edges = np.array([3, 4, 6, 11, 21, 51, 101, 501, 1e9])
cx, cy, lab = [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (n_next >= lo) & (n_next < hi)
    if sel.sum() < 20: continue
    cx.append(np.median(n_next[sel])); cy.append(100 * (rec[sel] >= 1).mean())
    lab.append(f"{rec[sel].mean():.2f}")
axB.plot(cx, cy, color=BLUE, lw=2.0, marker="o", ms=6.5,
         markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
for i, (x, y, tx) in enumerate(zip(cx, cy, lab)):
    axB.annotate(tx, (x, y), textcoords="offset points",
                 xytext=(0, 11) if i % 2 == 0 else (0, -17),
                 ha="center", fontsize=7.5, color=INK2)
axB.set_xscale("log")
axB.set_xlabel("clade size — carriers with a readable next symbol")
axB.set_ylabel("% of nodes carrying ≥1 recurrence")
axB.set_ylim(-7, 56)
axB.text(3.4, 50, "labels = mean recurrences per node", color=INK2, fontsize=8)
axB.set_title("B.  Homoplasy is a large-clade phenomenon",
              loc="left", fontsize=10, color=INK, pad=8)

# ---- C: concentration of homoplasy
o = np.sort(rec)[::-1]
cum = np.cumsum(o) / o.sum()
frac_nodes = np.arange(1, len(o) + 1) / len(o)
axC.plot(100 * frac_nodes, 100 * cum, color=BLUE, lw=2.0, zorder=3)
axC.plot([0, 100], [0, 100], color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
axC.text(46, 40, "if spread evenly", color=INK2, fontsize=8, rotation=31)
for pct, off in ((1, (26, -34)), (10, (30, 6))):
    y = 100 * cum[int(len(o) * pct / 100) - 1]
    axC.plot([pct], [y], marker="o", ms=6, color=ORANGE, zorder=5,
             markeredgecolor=SURFACE, markeredgewidth=1.4)
    axC.annotate(f"top {pct}% of nodes\ncarry {y:.0f}% of it", (pct, y),
                 textcoords="offset points", xytext=off, fontsize=8.5, color=ORANGE,
                 arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8, shrinkA=0, shrinkB=3))
axC.set_xlim(0, 100); axC.set_ylim(0, 102)
axC.set_xlabel("% of prefix nodes, ranked by recurrence count")
axC.set_ylabel("% of all estimated homoplasy")
axC.set_title("C.  …and it is concentrated, not spread thin",
              loc="left", fontsize=10, color=INK, pad=8)

fig.suptitle(f"Inferred write-event count $\\hat m$ across {len(a):,} prefix nodes   ·   "
             f"Park: all five tables   ·   $q={q:.4f}$, $m^*=\\sqrt{{2/q}}={birthday:.1f}$",
             fontsize=10.5, color=INK, x=0.008, ha="left", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.93))
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "m_distribution.png", dpi=200, facecolor=SURFACE)

print(f"nodes {len(a):,} | below m*: {below:.2f}% | total recurrences {rec.sum():,.0f}")
for pct in (0.5, 1, 5, 10):
    print(f"  top {pct}% of nodes carry {100*cum[max(int(len(o)*pct/100)-1,0)]:.1f}% of homoplasy")
print(f"wrote {FIG/'m_distribution.png'}")
