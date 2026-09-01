#!/usr/bin/env python3
r"""
================================================================================
 Near-linear compatible skeleton: deduplicate clades, then laminar insertion
================================================================================

Scripts 13/14 compute the full pairwise compatibility matrix, which is O(C^2) in
characters and is why the large clones OOM. But the SKELETON does not need that
matrix -- the pairwise fraction and conflict graph are DIAGNOSTICS. Building a
maximal compatible set is near-linear.

PHASE 1 -- DEDUPLICATE.  Many characters assert the identical clade (measured
redundancy 7-106x). They are trivially compatible and describe the same internal
node, so testing them against each other is wasted work. Group by clade and keep,
per distinct clade: multiplicity, and the set of TAPES asserting it.

  ⚑ The within-tape / cross-tape split of that redundancy is recorded, because it
  decides whether per-node support means anything. Two characters on the SAME tape
  (a prefix and its extension, with no branching between) are not independent
  evidence. Only distinct tapes are.

PHASE 2 -- INSERT, largest clade first.  Sorting by size descending pins the
nesting direction: every already-inserted clade K has |K| >= |S|, so "nested" can
only mean S subset of K, never the reverse.

Maintain owner[cell] = the SMALLEST inserted clade containing that cell (ROOT if
none). Then the entire insertion test is

        S is insertable  <=>  all cells of S share one owner

  Why: if all of S's cells have owner O then S subset of O. Any inserted K disjoint
  from O is disjoint from S; any K containing O contains S; and a child of O
  meeting S would have claimed those cells, contradicting the single owner. So S
  nests-or-disjoins with everything already placed. Two distinct owners means S
  straddles a boundary -- partial overlap, reject.

  Worked example, cells 1-8:
      A={1..6}  owners all ROOT      -> ACCEPT, owner[1..6]=A
      B={1,2,3} owners all A         -> ACCEPT, owner[1,2,3]=B
      C={3,4}   owner[3]=B, [4]=A    -> REJECT  (straddles B)
      D={4,5}   owners both A        -> ACCEPT

Cost: O(sum of clade sizes) for both phases, plus O(m log m) to sort. No conflict
graph, no C^2 anywhere; memory is O(cells) plus the clades themselves.

The result is a MULTIFURCATING tree. Nodes with >2 children are polytomies -- an
explicit statement that those lineages diverged and the characters cannot order
them, which is what we want to marginalise over rather than resolve arbitrarily.

VALIDATION: on Mouse3 this must land near script 16's C_greedy = 1,003, computed
by a completely different route (explicit conflict graph + min-degree greedy).
================================================================================
"""
import json, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parents[1] / "results"
MIN_CARRIERS = 3
ROOT = -1
TBL = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"
RESTARTS = int(sys.argv[2]) if len(sys.argv) > 2 else 200


def skeleton(clades, n, seed=None):
    """One greedy pass. Size-descending is REQUIRED (it pins the nesting direction);
    the order within each size class is free, which is what the restarts explore."""
    if seed is None:
        key = lambda i: (-clades[i][0], -clades[i][1])
        order = sorted(range(len(clades)), key=key)
    else:
        rng = np.random.default_rng(seed)
        jitter = rng.random(len(clades))
        order = sorted(range(len(clades)), key=lambda i: (-clades[i][0], jitter[i]))
    owner = np.full(n, ROOT, np.int32)
    accepted, parent = [], []
    for i in order:
        cells = clades[i][2]
        o = owner[cells]
        if o.size and (o == o[0]).all():
            nid = len(accepted)
            accepted.append(i); parent.append(int(o[0]))
            owner[cells] = nid
    return accepted, parent, owner


def best_skeleton(clades, n, restarts=0):
    """Deterministic support-ordered pass, then `restarts` randomised ones; keep the best."""
    acc, par, own = skeleton(clades, n)
    for s in range(restarts):
        a, p, o = skeleton(clades, n, seed=s)
        if len(a) > len(acc): acc, par, own = a, p, o
    return acc, par, own


def run(tbl):
    z = np.load(RES / f"characters_{tbl}.npz", allow_pickle=True)
    sizes = z["clone_sizes"]
    ch_clone, ch_tape, ch_start, ch_len = (z["ch_clone"], z["ch_tape"],
                                           z["ch_start"], z["ch_len"])
    car = z["carriers"]
    keep = np.flatnonzero(ch_len >= MIN_CARRIERS)
    by_clone = defaultdict(list)
    for i in keep:
        by_clone[int(ch_clone[i])].append(i)

    rows, t0 = [], time.time()
    tot_ch = tot_dedup = 0
    within_tape_extra = cross_tape_support = 0
    supp_hist = defaultdict(int)
    for ci, idx in sorted(by_clone.items()):
        n = int(sizes[ci])
        if n < 3: continue
        # ---- phase 1: deduplicate
        groups = defaultdict(lambda: [0, set(), None])       # key -> [mult, tapes, cells]
        for i in idx:
            cells = car[ch_start[i]:ch_start[i] + ch_len[i]]
            k = cells.tobytes()
            g = groups[k]
            g[0] += 1; g[1].add(int(ch_tape[i]))
            if g[2] is None: g[2] = cells
        tot_ch += len(idx); tot_dedup += len(groups)
        for mult, tapes, cells in groups.values():
            within_tape_extra += mult - len(tapes)
            cross_tape_support += len(tapes)
            supp_hist[min(len(tapes), 20)] += 1
        # ---- phase 2: insert
        clades = [(g[2].size, len(g[1]), g[2]) for g in groups.values()
                  if 2 <= g[2].size < n]
        acc, par, owner = best_skeleton(clades, n, RESTARTS)
        # polytomy structure: children = accepted nodes + cells still owned directly
        nkids = defaultdict(int)
        for p in par: nkids[p] += 1
        for c in range(n): nkids[int(owner[c])] += 1
        poly = [k for node, k in nkids.items() if node != ROOT and k > 2]
        rows.append(dict(clone=ci, cells=n, characters=len(idx),
                         distinct=len(groups), candidates=len(clades),
                         C=len(acc), frac=len(acc) / (n - 1),
                         polytomies=len(poly),
                         max_polytomy=int(max(poly)) if poly else 0,
                         mean_support=float(np.mean([clades[i][1] for i in acc])) if acc else 0.0))

    el = time.time() - t0
    nn = sum(r["cells"] - 1 for r in rows); CC = sum(r["C"] for r in rows)
    print(f"{tbl}: {len(rows)} clones, {RESTARTS} restarts, {sum(r['cells'] for r in rows):,} cells, {el:.1f}s\n")
    print(f"  characters {tot_ch:,} -> distinct clades {tot_dedup:,} "
          f"({tot_ch/max(tot_dedup,1):.1f}x redundancy)")
    print(f"  redundancy split: {cross_tape_support:,} clade-tape assertions, "
          f"{within_tape_extra:,} extra within-tape repeats "
          f"({100*within_tape_extra/max(tot_ch,1):.1f}% of characters)")
    print(f"  support (distinct tapes per clade): " +
          ", ".join(f"{k}{'+' if k==20 else ''}:{v:,}" for k,v in sorted(supp_hist.items())[:8]))
    print(f"\n  C = {CC:,} / {nn:,} internal nodes = {CC/nn:.3f}")
    print(f"  polytomies: {sum(r['polytomies'] for r in rows):,}, "
          f"largest {max(r['max_polytomy'] for r in rows):,}")
    print(f"\n  {'cells':>9} {'clones':>7} {'chars':>10} {'distinct':>9} {'C':>8} {'C/(n-1)':>9} {'supp':>6}")
    for lo, hi in [(3,4),(5,10),(11,20),(21,50),(51,400),(401,20000)]:
        g = [r for r in rows if lo <= r["cells"] <= hi]
        if not g: continue
        m = sum(r["cells"]-1 for r in g)
        print(f"  {f'{lo}-{hi}':>9} {len(g):>7} {sum(r['characters'] for r in g):>10,} "
              f"{sum(r['distinct'] for r in g):>9,} {sum(r['C'] for r in g):>8,} "
              f"{sum(r['C'] for r in g)/m:>9.3f} {np.mean([r['mean_support'] for r in g]):>6.1f}")
    (RES / f"skeleton_linear_{tbl}.json").write_text(json.dumps(rows, indent=1))
    return CC


if __name__ == "__main__":
    C = run(TBL)
    if TBL == "Mouse3":
        ref = json.loads((RES / "skeleton_C_Mouse3.json").read_text())
        rc = sum(r["C_greedy"] for r in ref)
        print(f"\n  CROSS-CHECK vs script 16 (conflict graph + min-degree greedy): "
              f"{C:,} vs {rc:,}  ({100*(C-rc)/rc:+.1f}%)")
