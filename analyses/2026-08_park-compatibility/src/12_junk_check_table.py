#!/usr/bin/env python3
r"""Clone-restriction test for junk values, on any table, against a calibrated null.

Script 11 established on Mouse1 that non-NNNNGGA values behave like heritable
lineage events rather than readout artifacts.  This repeats the decisive test on
another table.

The null MATTERS and differs by table.  Top-clone share is inflated whenever a few
clones dominate: Subclone has 15 clones with one holding ~28% of cells, so random
cells already score ~0.28.  The null is therefore computed by drawing random cell
sets of the same size from the same table, and the statistic reported is the
EXCESS over that null.

Confound checked: a per-sequencing-run artifact would be SAMPLE-restricted but not
CLONE-restricted.  Reporting both, plus their ratio, separates the two.
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = Path(__file__).resolve().parents[1] / "results"
MIN_COUNT, MIN_CARRIERS = 1000, 10
TBL = sys.argv[1] if len(sys.argv) > 1 else "Subclone"

d = json.loads((RES / "xi_vectors.json").read_text())[TBL]
keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}

clone, samp = {}, {}
with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if r["ClonalBC"] and r["ClonalBC"] != "None":
            clone[r["CellID"]] = r["ClonalBC"]; samp[r["CellID"]] = r["Sample"]

jt, rt, cells = defaultdict(list), defaultdict(list), []
with open(DATA / f"{TBL}_EditTable_filtered.csv", newline="") as fh:
    r = csv.reader(fh); hdr = next(r)[1:]
    blocks = defaultdict(list)
    for i, c in enumerate(hdr):
        tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
    blocks = {tp: [i for _, i in sorted(v)] for tp, v in blocks.items()}
    for row in r:
        cl = clone.get(row[0])
        if cl is None: continue
        cells.append((cl, samp[row[0]])); vals = row[1:]
        for tp, idxs in blocks.items():
            v = [vals[i] for i in idxs]
            for k, x in enumerate(v):
                if x in keep: rt[(tp, k+1, x)].append((cl, samp[row[0]])); continue
                if x != "None": jt[(tp, k+1, x)].append((cl, samp[row[0]]))
                break

def share(v, i):
    c = Counter(x[i] for x in v); return c.most_common(1)[0][1] / len(v)

J = [v for v in jt.values() if len(v) >= MIN_CARRIERS]
R = [v for v in rt.values() if len(v) >= MIN_CARRIERS]
rng = np.random.default_rng(0)
def null_curve(i=0, reps=200):
    """Median top-share for random cell sets, on a log grid of sizes; interpolated."""
    grid = np.unique(np.round(np.logspace(1, np.log10(min(len(cells), 20000)), 18)).astype(int))
    med = []
    for n in grid:
        idx = rng.integers(0, len(cells), size=(reps, int(n)))
        med.append(float(np.median([share([cells[j] for j in row], i) for row in idx])))
    return grid, np.array(med)

print(f"{TBL}: {len(cells):,} barcoded cells, {len(set(c for c,_ in cells))} clones, "
      f"{len(set(s for _,s in cells))} samples")
print(f"  junk triples >={MIN_CARRIERS} carriers: {len(J):,} | real: {len(R):,}\n")

if not J:
    print("  no junk triples pass the carrier threshold — test not informative here"); sys.exit()

g, m = null_curve(0)
nullof = lambda v: float(np.interp(len(v), g, m))
jc = np.array([share(v,0) for v in J]); jn = np.array([nullof(v) for v in J])
rc = np.array([share(v,0) for v in R]); rn = np.array([nullof(v) for v in R])
js = np.array([share(v,1) for v in J]); rs = np.array([share(v,1) for v in R])

print(f"{'':>8} {'n':>8} {'top-clone':>10} {'null':>8} {'EXCESS':>8} {'top-sample':>11} {'clone/sample':>13}")
for lab, c, nl, s in [("junk", jc, jn, js), ("real", rc, rn, rs)]:
    print(f"{lab:>8} {len(c):>8,} {np.median(c):>10.3f} {np.median(nl):>8.3f} "
          f"{np.median(c-nl):>+8.3f} {np.median(s):>11.3f} {np.median(c/s):>13.3f}")
print(f"\n  junk excess over null is {'GREATER' if np.median(jc-jn) > np.median(rc-rn) else 'SMALLER'} "
      f"than real symbols' ({np.median(jc-jn):+.3f} vs {np.median(rc-rn):+.3f})")
print(f"  => junk behaves {'like a heritable lineage event' if np.median(jc-jn) > 0.5*np.median(rc-rn) else 'like an artifact'}")
