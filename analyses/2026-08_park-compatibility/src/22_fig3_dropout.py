#!/usr/bin/env python3
r"""FIGURE 3 -- dropout, decomposed.

(a) the string "None" in the delivered tables means two different things: a tape
    that was never recovered, and a site inside a recovered tape that simply has
    not been edited yet. Pooling them gives the usual "~50% missing"; separating
    them shows most of it is one specific failure -- whole-tape dropout -- while
    the rest is informative.
(b) what that distinction buys: the fraction of cells whose membership in a
    depth-L character is actually KNOWN. Treating every None as missing makes this
    collapse with depth; recognising unedited sites as informative keeps it flat.
(c) whole-tape dropout per arm, against the cell-inclusion threshold each arm was
    filtered at -- part of the difference between arms is that QC choice.
"""
import csv, json, re
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Initial", "Subclone", "Mouse1", "Mouse2", "Mouse3"]
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_COUNT, N = 1000, 6

cache = RES / "dropout_decomposition.json"
if cache.exists():
    D = json.loads(cache.read_text())
else:
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    D = {}
    for a in ARMS:
        keep = {k for k, v in xi_all[a]["xi"].items() if v * xi_all[a]["n_edits"] >= MIN_COUNT}
        n_absent = n_tape = 0
        edits = unedited = junkent = 0
        det_corr = np.zeros(N + 1); det_naive = np.zeros(N + 1)
        with open(DATA / f"{a}_EditTable_filtered.csv", newline="") as fh:
            r = csv.reader(fh); hdr = next(r)[1:]
            blocks = defaultdict(list)
            for i, c in enumerate(hdr):
                tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
            blocks = {t: [i for _, i in sorted(v)] for t, v in blocks.items()}
            for row in r:
                v = row[1:]
                for idxs in blocks.values():
                    n_tape += 1
                    L, term = 0, None
                    for i in idxs:
                        if v[i] in keep: L += 1
                        else: term = "none" if v[i] == "None" else "junk"; break
                    if term is None: term = "complete"
                    if L == 0 and term == "none":
                        n_absent += 1; continue
                    edits += L
                    if term == "none": unedited += (N - L)
                    elif term == "junk": junkent += 1; unedited += (N - L - 1)
                    for k in range(1, N + 1):
                        if L >= k: det_naive[k] += 1
                        if L >= k or term in ("none", "complete"): det_corr[k] += 1
        tot = n_tape * N
        D[a] = dict(tape_absent=n_absent * N / tot, edits=edits / tot,
                    unedited=unedited / tot, junk=junkent / tot,
                    absent_frac=n_absent / n_tape,
                    det_corr=list(det_corr[1:] / n_tape), det_naive=list(det_naive[1:] / n_tape))
        print(f"  {a}: tape-absent {100*n_absent/n_tape:.1f}% of tapes | "
              f"entries edits {100*edits/tot:.1f}% unedited {100*unedited/tot:.1f}%", flush=True)
    cache.write_text(json.dumps(D, indent=1))

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE, AQUA = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 13.5})
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.4, 5.4))

# (a) the decomposition
x = np.arange(len(ARMS))
ed = np.array([100 * D[a]["edits"] for a in ARMS])
un = np.array([100 * D[a]["unedited"] for a in ARMS])
ab = np.array([100 * D[a]["tape_absent"] for a in ARMS])
axA.bar(x, ed, color=BLUE, width=0.62, label="edits")
axA.bar(x, un, bottom=ed, color=AQUA, width=0.62, label="unedited site (informative)")
axA.bar(x, ab, bottom=ed + un, color=ORANGE, width=0.62, label="tape not recovered")
axA.set_xticks(x); axA.set_xticklabels([LABEL[a] for a in ARMS], fontsize=11)
axA.set_ylabel("% of cell × tape × site entries")
axA.set_ylim(0, 100)
axA.legend(frameon=False, fontsize=11, loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=1)
axA.set_title("a   'None' means two different things", loc="left", pad=10)

# (b) determined fraction vs depth
lv = np.arange(1, N + 1)
for a in ARMS:
    axB.plot(lv, 100 * np.array(D[a]["det_naive"]), color=AXIS, lw=2.0)
    axB.plot(lv, 100 * np.array(D[a]["det_corr"]), color=BLUE, lw=2.4)
axB.text(N, 100 * D["Mouse1"]["det_corr"][-1] + 3, "unedited counted as informative",
         color=BLUE, ha="right", fontsize=11.5, fontweight="bold")
axB.text(N, 100 * D["Mouse1"]["det_naive"][-1] - 5, "every None counted as missing",
         color=INK2, ha="right", fontsize=11.5)
axB.set_xlabel("prefix depth $L$"); axB.set_ylabel("% of cells whose membership is known")
axB.set_ylim(0, 100); axB.set_xticks(lv)
axB.set_title("b   Recovering the distinction keeps the data usable", loc="left", pad=10)

# (c) dropout per arm vs the QC threshold
absf = [100 * D[a]["absent_frac"] for a in ARMS]
axC.bar(x, absf, color=ORANGE, width=0.62)
for i, v in enumerate(absf):
    axC.text(i, v + 1.2, f"{v:.0f}%", ha="center", fontsize=12, color=INK)
for i, a in enumerate(ARMS):
    axC.text(i, 3, f"≥{THR[a]}\ntapes", ha="center", fontsize=10.5, color="white")
axC.set_xticks(x); axC.set_xticklabels([LABEL[a] for a in ARMS], fontsize=11)
axC.set_ylabel("% of cell × tape instances not recovered")
axC.set_ylim(0, max(absf) * 1.35)
axC.set_title("c   Dropout differs by arm — partly by QC choice", loc="left", pad=10)

fig.tight_layout()
fig.savefig(FIG / "fig3_dropout.png", dpi=200, facecolor=SURF)
print("wrote fig3_dropout.png")
