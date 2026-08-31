#!/usr/bin/env python3
r"""C = how many of the n-1 internal nodes the compatible skeleton fixes (D.4b).

D.4b is blunt that the compatibility FRACTION is only a means: "a high compatibility
fraction with a low C would be a null result -- lots of agreeing characters that all
describe the same few clades. Measure C or the check has not answered the question."

The maximal compatible set is max-clique in general. But a rigorous LOWER BOUND on C
needs no search: characters with conflict degree 0 conflict with nothing, so they are
pairwise compatible by construction, and their clades form a laminar family. Counting
its distinct non-trivial clades gives C_lower.

Reported per clone against n-1, plus the redundancy (characters per distinct clade),
which is exactly the "lots of characters describing the same few clades" failure D.4b
warns about.

Degrees are the missing-excluded convention (script 14 saves that one).
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parents[1] / "results"
MIN_CARRIERS = 3
TBL = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"

z = np.load(RES / f"characters_{TBL}.npz", allow_pickle=True)
g = np.load(RES / f"conflict_degrees_{TBL}.npz")
sizes = z["clone_sizes"]
deg_by_clone = defaultdict(dict)
for cl, ch, d in zip(g["clone"], g["char"], g["degree"]):
    deg_by_clone[int(cl)][int(ch)] = int(d)

rows = []
for ci, degs in sorted(deg_by_clone.items()):
    n = int(sizes[ci])
    sel = np.flatnonzero((z["ch_clone"] == ci) & (z["ch_len"] >= MIN_CARRIERS))
    if sel.size == 0 or n < 3: continue
    free = [k for k in range(sel.size) if degs.get(k, 0) == 0]
    clades = defaultdict(int)
    for k in free:
        i = sel[k]
        car = z["carriers"][z["ch_start"][i]:z["ch_start"][i] + z["ch_len"][i]]
        if 2 <= car.size < n:                       # non-trivial: not a singleton, not everything
            clades[frozenset(car.tolist())] += 1
    C = len(clades)
    rows.append(dict(clone=ci, cells=n, characters=int(sel.size),
                     conflict_free=len(free), C=C,
                     frac=C / (n - 1) if n > 1 else float("nan"),
                     redundancy=(sum(clades.values()) / C) if C else float("nan"),
                     clades=[sorted(s) for s in clades] if n <= 40 else None))

# laminarity check on a sample of clones -- must hold by construction
bad = 0
for r in rows[:40]:
    if not r["clades"]: continue
    S = [set(c) for c in r["clades"]]
    for a in range(len(S)):
        for b in range(a + 1, len(S)):
            x, y = S[a], S[b]
            if (x & y) and not (x <= y or y <= x): bad += 1
print(f"{TBL}: laminarity of the conflict-free clades — "
      f"{'PASS (all nested or disjoint)' if bad == 0 else f'FAIL {bad}'}\n")

tot_n = sum(r["cells"] for r in rows); tot_C = sum(r["C"] for r in rows)
print(f"{'cells':>9} {'clones':>7} {'chars':>10} {'conflict-free':>14} {'C':>8} {'C/(n-1)':>9} {'redundancy':>11}")
for lo, hi in [(3, 4), (5, 10), (11, 20), (21, 50), (51, 400)]:
    grp = [r for r in rows if lo <= r["cells"] <= hi]
    if not grp: continue
    nn = sum(r["cells"] - 1 for r in grp); CC = sum(r["C"] for r in grp)
    red = np.mean([r["redundancy"] for r in grp if r["C"]])
    print(f"{f'{lo}-{hi}':>9} {len(grp):>7} {sum(r['characters'] for r in grp):>10,} "
          f"{sum(r['conflict_free'] for r in grp):>14,} {CC:>8,} {CC/nn:>9.3f} {red:>11.1f}")
nn = sum(r["cells"] - 1 for r in rows)
print(f"\n  TOTAL: {len(rows)} clones, {tot_n:,} cells, C = {tot_C:,} of {nn:,} internal nodes "
      f"= {tot_C/nn:.3f}")
print(f"  conflict-free characters: {sum(r['conflict_free'] for r in rows):,} of "
      f"{sum(r['characters'] for r in rows):,} "
      f"({100*sum(r['conflict_free'] for r in rows)/sum(r['characters'] for r in rows):.1f}%)")
(RES / f"skeleton_C_{TBL}.json").write_text(json.dumps(
    [{k: v for k, v in r.items() if k != "clades"} for r in rows], indent=1))
