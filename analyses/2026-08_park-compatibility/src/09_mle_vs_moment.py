#!/usr/bin/env python3
r"""Compare the Poissonised MLE against the count-only estimator, and use the
comparison as a MODEL CHECK.

The two estimators are both unbiased under the assumed model (script 08's Monte
Carlo confirms this), so where they DISAGREE on real data, the model is wrong,
not one of the estimators.

The diagnostic is the observed mass W_A = sum_{i in A} xi_i.

    *** A TRAP, WHICH THE FIRST VERSION OF THIS SCRIPT FELL INTO ***

The obvious null is E[W_A | m] = sum_i xi_i * (1 - (1-xi_i)^m), evaluated at the
m implied by s.  That is WRONG, because it conditions on m while the data are
selected on s.  For fixed m the realised s is random, and a node showing an
unusually LARGE s got there by hitting unusually many RARE symbols -- which drags
W_A down.  So E[W_A | s] < E[W_A | m(s)] even when the model is perfectly true,
and the naive comparison manufactures a spurious deficit at large s.

The correct null conditions on s, exactly as the data do.  We get it by simulation:
draw m over a wide grid, record (s, W_A) for each replicate, then bin by the
REALISED s.  That yields the model's distribution of W_A given s, with the
selection effect included, and the observed W_A can be compared against it.
"""
import gzip, json
import numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
TABLES = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]

rows = []
for t in TABLES:
    with gzip.open(RES / f"m_mle_{t}.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((int(f[4]), float(f[5]), float(f[6]), float(f[7])))
A = np.array(rows)
s, W, mom, mle = A[:, 0], A[:, 1], A[:, 2], A[:, 3]

# reference xi (Mouse1, trimmed) for the predicted-mass curve
d = json.loads((RES / "xi_vectors.json").read_text())["Mouse1"]
xi = np.array([v for k, v in d["xi"].items() if v * d["n_edits"] >= 1000]); xi /= xi.sum()
# ---- correct null: simulate, then condition on the REALISED s
rng = np.random.default_rng(0)
sim_s, sim_W = [], []
for m_true in np.unique(np.round(np.logspace(0, 3.6, 90)).astype(int)):
    reps = 400 if m_true < 300 else 150
    for _ in range(reps):
        idx = np.unique(rng.choice(len(xi), size=int(m_true), p=xi))
        sim_s.append(len(idx)); sim_W.append(xi[idx].sum())
sim_s, sim_W = np.array(sim_s), np.array(sim_W)
print(f"null simulated: {len(sim_s):,} replicates, s spans {sim_s.min()}-{sim_s.max()}")

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.6, 4.5))

bins = [(1,1),(2,3),(4,6),(7,12),(13,25),(26,50),(51,80),(81,120)]
cx, med, lo, hi, obsW, predW, predLo, predHi = [], [], [], [], [], [], [], []
for a, b in bins:
    sel = (s >= a) & (s <= b)
    if sel.sum() < 20: continue
    r = mle[sel] / np.maximum(mom[sel], 1e-9)
    cx.append(np.median(s[sel])); med.append(np.median(r))
    lo.append(np.percentile(r, 25)); hi.append(np.percentile(r, 75))
    obsW.append(np.median(W[sel]))
    nsel = (sim_s >= a) & (sim_s <= b)                 # condition on s, as the data do
    predW.append(float(np.median(sim_W[nsel])) if nsel.sum() >= 30 else np.nan)
    predLo.append(float(np.percentile(sim_W[nsel], 5)) if nsel.sum() >= 30 else np.nan)
    predHi.append(float(np.percentile(sim_W[nsel], 95)) if nsel.sum() >= 30 else np.nan)

axA.axhline(1.0, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
axA.fill_between(cx, lo, hi, color=BLUE, alpha=0.18, lw=0, zorder=2)
axA.plot(cx, med, color=BLUE, lw=2.0, marker="o", ms=6.5, markerfacecolor=BLUE,
         markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
axA.set_xscale("log")
axA.set_xlabel("$s$ — distinct symbols at the node")
axA.set_ylabel(r"$\hat m_{\rm MLE}\,/\,\hat m_{\rm count}$")
axA.text(1.05, 0.885, "the two estimators agree\nwherever $s$ is small", color=INK2, fontsize=8.5)
axA.annotate("at large $s$ the MLE\nreads 20% fewer events", xy=(cx[-1], med[-1]),
             xytext=(0.42, 0.24), textcoords="axes fraction", color=ORANGE, fontsize=8.5,
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9, shrinkB=4))
axA.set_title("A.  Identical where $s$ is small, divergent where it is large",
              loc="left", fontsize=10, color=INK, pad=8)

axB.fill_between(cx, predLo, predHi, color=AXIS, alpha=0.30, lw=0, zorder=2)
axB.plot(cx, predW, color=AXIS, lw=2.0, marker="s", ms=5.5, markerfacecolor=SURFACE,
         markeredgecolor=AXIS, markeredgewidth=1.6, zorder=3)
axB.plot(cx, obsW, color=BLUE, lw=2.0, marker="o", ms=6.5, markerfacecolor=BLUE,
         markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=4)
axB.text(1.15, 0.80, "OBSERVED", color=BLUE,
         fontsize=8.5, fontweight="bold")
axB.text(1.15, 0.62, "NULL: iid draws, conditioned\non the same $s$ (90% band)", color=INK2, fontsize=8.5)
axB.set_xscale("log"); axB.set_ylim(0, 1.05)
axB.set_xlabel("$s$ — distinct symbols at the node")
axB.set_ylabel("$W_A$ — total $\\xi$-mass of the observed set")
axB.set_title("B.  Observed mass sits inside the null — no misspecification detected",
              loc="left", fontsize=10, color=INK, pad=8)

fig.suptitle("Poissonised MLE vs count-only estimator, 1,567,321 prefix nodes   ·   "
             "they part company only near saturation, where $\\hat m$ is acutely sensitive to the unseen mass",
             fontsize=10.5, color=INK, x=0.008, ha="left", y=0.985)
fig.tight_layout(rect=(0, 0, 1, 0.92))
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "mle_vs_moment.png", dpi=200, facecolor=SURFACE)

print(f"\n{'s':>8} {'nodes':>10} {'obs W_A':>9} {'null med':>9} {'null 5-95%':>16} {'inside?':>8} {'m ratio':>9}")
for (a, b), x, ow, pw, mr in zip([x for x in bins if ((s>=x[0])&(s<=x[1])).sum()>=20],
                                 cx, obsW, predW, med):
    sel = (s >= a) & (s <= b)
    i = cx.index(x); band = f"[{predLo[i]:.3f}, {predHi[i]:.3f}]"
    ok = "yes" if predLo[i] <= ow <= predHi[i] else "NO"
    print(f"{f'{a}-{b}':>8} {sel.sum():>10,} {ow:>9.3f} {pw:>9.3f} {band:>16} {ok:>8} {mr:>9.3f}")
print(f"\nwrote {FIG/'mle_vs_moment.png'}")
