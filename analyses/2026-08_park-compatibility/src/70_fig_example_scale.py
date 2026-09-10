#!/usr/bin/env python3
r"""
Panel a at ANY scale -- the raster of one clone with a clade bracketed and a few
tapes singled out.  Generalises 69's panel a so both examples come from one place.

  small scale : --anchor 122 --depth 5 --clade allmiss:40 --keys 122,102,40
                a 38-cell clade; tape 40 missing in 7% of the population, so the
                block is VISIBLE.
  large scale : --anchor 125 --depth 3 --clade largest   --keys 125,2,102
                a 1,304-cell clade; tape 102 is 99.4% missing inside against
                88.3% outside, so the block is NOT visible -- and the test still
                resolves it.  Tape 2 is missing in all 3,141 cells: sigma^2 = 0,
                so it is UNTESTABLE, the third state a tape can be in.

⚑ Why no clean block exists at the large scale: the clade is 42% of its block
population, so a complete loss forces the population rate above 0.42.  A large
clade IS most of its own background.  That is a property of the geometry, not of
this dataset, and it is the reason the eye stops being a useful detector as the
scale grows.

Usage: 70_fig_example_scale.py <arm> --clone C --anchor A --depth D
         --clade largest|allmiss:T --keys a,b,c --tag NAME [--maxrows 150]
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
CL, AN, DEP = arg("--clone", 76, int), arg("--anchor", 122, int), arg("--depth", 5, int)
PICK, TAG = arg("--clade", "largest"), arg("--tag", "example")
KEYS = [int(x) for x in arg("--keys", f"{AN}").split(",")]
MAXR, NBG = arg("--maxrows", 150, int), arg("--nbg", 36, int)
BONF = arg("--bonf", 1.715e-07, float)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
K = Y.shape[1]
_, g = np.unique(clone, return_inverse=True)
miss = ~Y
codes = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"]
cd = codes[:, AN, DEP - 1]
pop = np.flatnonzero((cd >= 0) & (g == CL))
inclone = np.flatnonzero(g == CL)
n_exc = inclone.size - pop.size
anchor_missing_clone = int(miss[inclone, AN].sum())
vals, inv = np.unique(cd[pop], return_inverse=True)
sizes = np.bincount(inv)
if PICK == "largest":
    ci = int(np.argmax(sizes))
else:
    t0 = int(PICK.split(":")[1])
    ci = [i for i in range(vals.size) if miss[pop[inv == i], t0].all()
          and sizes[i] >= 4][int(0)]
clade, rest = pop[inv == ci], pop[inv != ci]
m, npop = clade.size, pop.size
print(f"{arm} clone {CL}, anchor {AN} depth {DEP}: population {npop}, clade {m}, rest {rest.size}")
info = {}
for t in KEYS:
    k = int(miss[clade, t].sum()); kp = int(miss[pop, t].sum())
    out = (kp - k) / max(npop - m, 1)
    p = float(stats.hypergeom.sf(k - 1, npop, kp, m))
    testable = 0 < kp < npop
    info[t] = dict(k=k, kp=kp, pi=k / m, rate=kp / npop, out=out, p=p, testable=testable)
    print(f"   tape {t:>4}: k={k}/{m} (pi={k/m:.3f})  outside {out:.3f}  "
          f"population {kp/npop:.3f}  p={p:.2e}  {'testable' if testable else 'UNTESTABLE'}")

rng = np.random.default_rng(0)
def take(a, n):
    return a if a.size <= n else np.sort(rng.choice(a, size=n, replace=False))
rc, rr = take(clade, MAXR), take(rest, MAXR)
rr = rr[np.argsort(cd[rr])]
rowidx = np.concatenate([rc, rr])
bgpool = [t for t in range(K) if t not in KEYS]
bg = sorted(rng.choice(bgpool, size=NBG, replace=False))
M = miss[np.ix_(rowidx, KEYS + bg)]

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
BLUE, DARK = "#2a78d6", "#184f95"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})
nk = len(KEYS)
fig, (axk, axb) = plt.subplots(1, 2, figsize=(3.0 + 0.55 * nk + 5.0, 5.4), sharey=True,
                               gridspec_kw=dict(width_ratios=[0.42 * nk, 2.0], wspace=0.06))
cm = ListedColormap(["#eeeeea", DARK])
axk.imshow(M[:, :nk], aspect="auto", interpolation="nearest", cmap=cm, vmin=0, vmax=1)
axb.imshow(M[:, nk:], aspect="auto", interpolation="nearest", cmap=cm, vmin=0, vmax=1)
for ax in (axk, axb):
    ax.axhline(rc.size - 0.5, color=INK, lw=1.6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_edgecolor(GRID)
for i, t in enumerate(KEYS):
    d = info[t]
    if t == AN:
        lab, col, bold = f"anchor {t}\n(defines\nthe clade)", INK2, False
    elif not d["testable"]:
        lab, col, bold = f"tape {t}\nmissing in\nall cells\nUNTESTABLE", INK2, False
    elif d["p"] <= BONF:
        lab, col, bold = f"tape {t}\n\nCALLED", BLUE, True
    else:
        lab, col, bold = f"tape {t}\nnot\ncalled", INK2, False
    axk.add_patch(plt.Rectangle((i - 0.5, -0.5), 1, M.shape[0], fill=False,
                                edgecolor=col, lw=2.2, zorder=5))
    axk.annotate(lab, xy=(i, -0.5), xytext=(0, 8), textcoords="offset points",
                 ha="center", va="bottom", fontsize=8.0, color=col,
                 fontweight="bold" if bold else "normal", annotation_clip=False)
    if t == AN:
        foot = "missing\n0% inside\n0% outside"
    elif not d["testable"]:
        foot = f"missing\n{d['pi']:.0%} inside\n{d['out']:.0%} outside"
    else:
        foot = (f"missing\n{d['pi']:.1%} inside\n{d['out']:.1%} outside\n"
                f"p = {d['p']:.0e}")
    axk.annotate(foot, xy=(i, M.shape[0] - 0.5), xytext=(0, -9),
                 textcoords="offset points", ha="center", va="top", fontsize=8.0,
                 color=col, fontweight="bold" if bold else "normal",
                 annotation_clip=False)
axb.annotate(f"{NBG} other tapes, for background", xy=(NBG / 2, -0.5), xytext=(0, 8),
             textcoords="offset points", ha="center", va="bottom", fontsize=8.5,
             color=INK2, annotation_clip=False)
axk.annotate(f"the clade\n{m:,} cells" + (f"\n{rc.size} shown" if rc.size < m else "\nall shown"),
             xy=(-0.75, rc.size / 2), ha="right", va="center", fontsize=9,
             color=INK, fontweight="bold", annotation_clip=False)
axk.annotate("rest of the population\n" +
             (f"{rr.size} of {rest.size:,} shown" if rr.size < rest.size
              else f"{rest.size:,} cells"),
             xy=(-0.75, rc.size + rr.size / 2), ha="right", va="center", fontsize=9,
             color=INK2, annotation_clip=False)
sub = (f"Called when p ≤ {BONF:.1e} — one expected false positive across the whole arm. "
       f"p is the exact hypergeometric tail; the catalogue reports the more conservative "
       f"max(normal, hypergeometric).\n"
       f"⚠ The anchor is recovered in all {npop:,} cells shown BY CONSTRUCTION — a cell needs "
       f"it to have a prefix. Across the whole clone it is itself missing in "
       f"{anchor_missing_clone:,} of {inclone.size:,} ({anchor_missing_clone/inclone.size:.0%}), "
       f"which is why {n_exc:,} cells are excluded here.")
fig.suptitle(arg("--title", "A tape lost on one clade").replace("_", " "), fontsize=11, fontweight="bold",
             x=0.11, ha="left", y=1.085)
fig.text(0.11, 1.022, f"dark = tape missing    ·    {arm} clone {CL}, rows ordered by the "
         f"anchor's depth-{DEP} prefix", fontsize=8.5, color=INK2, ha="left")
fig.text(0.11, -0.105, sub, fontsize=8.5, color=INK, ha="left")
fig.savefig(FIG / f"fig4c_a_{TAG}.png"); plt.close(fig)
print(f"wrote figures/fig4c_a_{TAG}.png")
