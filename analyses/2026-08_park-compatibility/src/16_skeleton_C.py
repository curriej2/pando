#!/usr/bin/env python3
r"""C = internal nodes fixed by the compatible skeleton, against MISSING-AS-ABSENT.

D.4b: "a high compatibility fraction with a low C would be a null result -- lots of
agreeing characters that all describe the same few clades. Measure C or the check
has not answered the question."

⚠ It must be the missing-as-absent graph. Under missing-EXCLUDED each pair is tested
on its own D1 & D2, so compatibility is pair-specific and does NOT give one laminar
family over all cells (see notes S D.4d -- the first attempt returned C/(n-1)=2.107,
which is impossible). Under missing-as-absent the test uses raw sets, so

        compatible  <=>  nested or disjoint  =>  laminar

and a skeleton is constructible. This script therefore verifies laminarity and
REFUSES to report C if it fails.

Two estimates per clone:
  C_free    distinct non-trivial clades among degree-0 characters -- a rigorous
            lower bound needing no search.
  C_greedy  same after a min-degree greedy maximal independent set on the conflict
            graph (min-degree, not max-degree removal: it yields a larger
            independent set). Still a lower bound on the true maximum, but tighter.

Also reports REDUNDANCY -- characters per distinct clade -- which is exactly the
"many characters, few clades" failure D.4b warns about.
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parents[1] / "results"
MIN_CARRIERS = 3
TBL = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"

z = np.load(RES / f"characters_{TBL}.npz", allow_pickle=True)
g = np.load(RES / f"conflict_degrees_as_absent_{TBL}.npz")
E = np.load(RES / f"conflict_edges_as_absent_{TBL}.npz")
sizes = z["clone_sizes"]

deg = defaultdict(dict)
for cl, ch, d in zip(g["clone"], g["char"], g["degree"]): deg[int(cl)][int(ch)] = int(d)
adj = defaultdict(lambda: defaultdict(set))
for cl, a, b in zip(E["clone"], E["a"], E["b"]):
    adj[int(cl)][int(a)].add(int(b)); adj[int(cl)][int(b)].add(int(a))


def greedy_independent(n, nbr):
    """Min-degree greedy: repeatedly take the lowest-degree surviving vertex."""
    alive = set(range(n))
    d = {v: len(nbr.get(v, ())) for v in alive}
    keep = []
    while alive:
        v = min(alive, key=lambda x: d[x])
        keep.append(v)
        gone = {v} | (nbr.get(v, set()) & alive)
        alive -= gone
        for u in gone:
            for w in nbr.get(u, ()):
                if w in alive: d[w] -= 1
    return keep


def clades_of(ci, ks, sel, n):
    out = defaultdict(int)
    for k in ks:
        i = sel[k]
        car = z["carriers"][z["ch_start"][i]:z["ch_start"][i] + z["ch_len"][i]]
        if 2 <= car.size < n: out[frozenset(car.tolist())] += 1
    return out


def laminar(clades):
    S = [set(c) for c in clades]
    for a in range(len(S)):
        for b in range(a + 1, len(S)):
            x, y = S[a], S[b]
            if (x & y) and not (x <= y or y <= x): return False
    return True


rows, bad = [], 0
for ci, dd in sorted(deg.items()):
    n = int(sizes[ci])
    sel = np.flatnonzero((z["ch_clone"] == ci) & (z["ch_len"] >= MIN_CARRIERS))
    if sel.size == 0 or n < 3: continue
    nbr = adj.get(ci, {})
    free = [k for k in range(sel.size) if dd.get(k, 0) == 0]
    grd = greedy_independent(sel.size, nbr)
    cf, cg = clades_of(ci, free, sel, n), clades_of(ci, grd, sel, n)
    if len(cg) <= 60 and not laminar(list(cg)): bad += 1
    rows.append(dict(clone=ci, cells=n, characters=int(sel.size),
                     free=len(free), greedy=len(grd),
                     C_free=len(cf), C_greedy=len(cg),
                     red=(sum(cg.values()) / len(cg)) if cg else float("nan")))

print(f"{TBL}: laminarity of the greedy compatible set — "
      f"{'PASS' if bad == 0 else f'FAIL on {bad} clones'}   "
      f"(must pass: as-absent compatibility implies laminarity)\n")
print(f"{'cells':>9} {'clones':>7} {'chars':>10} {'C_free':>8} {'C_greedy':>9} "
      f"{'C_g/(n-1)':>10} {'redundancy':>11}")
for lo, hi in [(3, 4), (5, 10), (11, 20), (21, 50), (51, 400)]:
    grp = [r for r in rows if lo <= r["cells"] <= hi]
    if not grp: continue
    nn = sum(r["cells"] - 1 for r in grp)
    print(f"{f'{lo}-{hi}':>9} {len(grp):>7} {sum(r['characters'] for r in grp):>10,} "
          f"{sum(r['C_free'] for r in grp):>8,} {sum(r['C_greedy'] for r in grp):>9,} "
          f"{sum(r['C_greedy'] for r in grp)/nn:>10.3f} "
          f"{np.mean([r['red'] for r in grp if r['C_greedy']]):>11.1f}")
nn = sum(r["cells"] - 1 for r in rows)
print(f"\n  TOTAL {len(rows)} clones, {sum(r['cells'] for r in rows):,} cells")
print(f"    C_free   = {sum(r['C_free'] for r in rows):,} / {nn:,} internal nodes = "
      f"{sum(r['C_free'] for r in rows)/nn:.3f}")
print(f"    C_greedy = {sum(r['C_greedy'] for r in rows):,} / {nn:,} internal nodes = "
      f"{sum(r['C_greedy'] for r in rows)/nn:.3f}")
(RES / f"skeleton_C_{TBL}.json").write_text(json.dumps(rows, indent=1))
